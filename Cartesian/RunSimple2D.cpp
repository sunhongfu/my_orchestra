// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <MDArray/Utils.h>

#include <Dicom/MR/Image.h>
#include <Dicom/Network.h>

#include <Orchestra/Cartesian2D/LxControlSource.h>

#include <Orchestra/Cartesian2D/KSpaceTransformer.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/SliceOrientation.h>
#include <Orchestra/Common/SliceCorners.h>
#include <Orchestra/Common/ReconTrace.h>

#include <Orchestra/Control/ProcessingControl.h>
#include <Orchestra/Control/CommonTypes.h>

#include <Orchestra/Core/SumOfSquares.h>
#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CommandLine.h"
#include "Rehearsal.h"

// Include this to avoid having to type fully qualified names
using namespace GERecon;
using namespace MDArray;

/**
 * Function that does a simple 2D Recon by creating a ProcessingControl object from a Pfile and
 * looping through each slice/echo/channel to reconstruct an image from the Pfile raw data.
 * 
 * Note: this will attempt to reconstruct all acquisitions associated with the specified Pfile.
 *
 * Limitations: 
 *  No parallel imaging (Asset or Arc)
 *  No retrospective phase correction
 *
 * @author Matt Bingen
 */
void GERecon::RunSimple2D(const Legacy::PfilePointer& pfile)
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("Cartesian2DRehearsal");

    // Create DICOM series to save images into
    const Legacy::DicomSeries dicomSeries(pfile);
    
    // Get DICOM network, series number, and series description (if specified) from the command line.
    // If not specified, the optionals and pointer will be empty and no attempt will be made to insert
    // the values or store the images.
    const GEDicom::NetworkPointer dicomNetwork = CommandLine::DicomNetwork();
    const boost::optional<int> seriesNumber = CommandLine::SeriesNumber();
    const boost::optional<std::string> seriesDescription = CommandLine::SeriesDescription();

    // Create a ProcessingControl from the Legacy Pfile. Note that the return type is a 
    // "ProcessingControlPointer". This type is defined in CommonTypes.h and should be thought of
    // as a pointer to a ProcessingControl object. However, it's a little more than a standard 
    // pointer that would be declared like "ProcessingControl* processingControl". These pointers
    // are referred to as "smart pointers", specifically in this case it is of type 
    // boost::shared_ptr. These pointers do automatic reference counting and are typically referred
    // to as memory managing objects. These pointers still point to an underlying object and the
    // "->" operator is used just like a normal "dumb" pointer to access the members of the pointed
    // to object. In this case a ProcessingControl object. Smart pointers should be used whenever
    // possible, and with this there should almost never be calls to 'malloc', 'free', or delete.
    // Developers should even strive to avoid 'new' and prefer boost::make_shared<>() instead!
    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfile->CreateOrchestraProcessingControl<Cartesian2D::LxControlSource>();

    // Pull info directly out of ProcessingControl using the Value function. Parameters are held
    // in ProcessingControl objects as a Name/Value pair. The name is held as a std::string and the
    // value can be of any type, which is determined when adding the parameter. Note that 
    // ProcessingControl can contain more than just the basic types, it can hold user defined
    // types/classes
    //
    // The support for variable types is accomplished with the C++ feature known as templates. Both
    // functions and classes can be templated and in this case the member function of 
    // ProcessingControl, Value(), is templated with a variable type.
    //
    // The signature of the function is as follows:
    //
    // template<class T>
    // const T& Value(const std::string& paramName) const;
    //
    // Don't get caught up with all the "const" in the function, they could be removed and the idea
    // would be the same. But for completeness, the first "const" means the return type is to be
    // constant, the second means the input string paramName will remain constant throughout the
    // function, and the last "const" means that the function does not change any data members
    // within the class in which it's defined.
    // 
    // The "template<class T>" in front of the declaration means that the function contains a
    // variable (templated) type (either as a return type and/or a parameter). Here the templated
    // parameter is just the return type. So, for a given string value we can get any number of 
    // return types. Of course the function needs to know what type to return when calling the 
    // function and since it cannot determine the type automatically, it must be specified. To 
    // specify the type or templated parameter, the "<type>" syntax is appended after the function
    // name and before the parameter list. Example ValueStrict<myFavoriteType>("SomeKey"); This may 
    // look a little funny at first, but it becomes more obvious once it becomes familiar and the 
    // general concept is understood.
    //
    // A few caveats: 
    // What if the string name doesn't exist in ProcessingControl?
    // -An exception will be thrown and a nice message will be displayed.
    // What if the templated parameter is incorrect?
    // -Another exception! You may get lucky with an implicit conversion.
    // 
    // This is a good example of how to get used to the template syntax. And since the names in
    // ProcessingControl are known, along with the types, trouble should be avoided.

    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");

    // Pulled directly from Pfile API and not ProcessingControl since the pfile class knows exactly
    // how many slices and acquisitions it has loaded and can produce raw data for.
    const int numPhases = pfile->PhaseCount();
    const int numSlices = pfile->SliceCount();
    const int numEchoes = pfile->EchoCount();
    const int numChannels = pfile->ChannelCount();

    // Create a transformer object to complete the K-Space to image-space transformation.
    // This will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformer(*processingControl);

    // Create channel combiner object that will do the channel combining work in channel loop.
    SumOfSquares channelCombiner(channelWeights);

    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient); 

    // Storage for transformed image data, thus the image sizes.
    ComplexFloatMatrix imageData(imageXRes, imageYRes);

    for(int currentPhase = 0; currentPhase < numPhases; ++currentPhase)
    {
        for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
        {
            for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
            {
                // Zero out channel combiner buffer for the next set of channels.
                channelCombiner.Reset();

                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    trace.ConsoleMsg("Processing Slice[%d of %d], Echo[%d of %d], Channel[%d of %d]", 
                        currentSlice+1, numSlices, currentEcho+1, numEchoes, currentChannel+1, numChannels);

                    // Extract a single channel of data from the Pfile
                    const ComplexFloatMatrix kSpace = pfile->KSpaceData<float>(currentSlice, currentEcho, currentChannel, currentPhase);

                    // Transform to image space. Data will be zipped from acquired size.
                    transformer.Apply(imageData, kSpace);

                    // Accumulate channel data in channel combiner.
                    channelCombiner.Accumulate(imageData, currentChannel);
                }

                // Get the combined channel image from the ChannelCombiner
                const ComplexFloatMatrix combinedImage = channelCombiner.GetCombinedImage();

                // Create storage for final image represented as a float matrix
                FloatMatrix magnitudeImage(combinedImage.shape());

                // Convert complex data to specified image type.
                MDArray::ComplexToReal(magnitudeImage, combinedImage, MDArray::MagnitudeData);

                // Get information for current slice
                const SliceOrientation& sliceOrientation = pfile->Orientation(currentSlice);
                const SliceCorners& prescribedCorners = pfile->Corners(currentSlice);
                const SliceCorners& acquiredCorners = pfile->AcquiredCorners(currentSlice);

                // Perform Gradwarp
                gradwarp.Run(magnitudeImage, currentSlice, prescribedCorners, acquiredCorners);

                // Rotate/Transpose image accordingly
                FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

                // Clip the image
                Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

                // Create storage for final image represented as a short matrix - form expected for DICOM.
                ShortMatrix finalImage(rotatedImage.shape());

                // Cast final image from float to short
                finalImage = MDArray::cast<short>(rotatedImage);

                // Create DICOM image
                const int imageNumber = ImageNumber2D(currentPhase, currentSlice, currentEcho, pfile);
                const ImageCorners imageCorners(prescribedCorners, sliceOrientation);
                const GEDicom::MR::ImagePointer dicom = dicomSeries.NewImage(finalImage, imageNumber, imageCorners);

                // Add custom annotation if specified on command line.
                // Note, boost::optional<> types only gets inserted if set.
                dicom->Insert<GEDicom::LongString>(0x0008, 0x103E, seriesDescription); // Series description, showed in image browser
                dicom->Insert<GEDicom::LongString>(0x0008, 0x1090, seriesDescription); // Manufacturer's model name, showed in annotation
                dicom->Insert<GEDicom::IntegerString>(0x0020, 0x0011, seriesNumber);

                // Save DICOM to file and also store it if network is active
                std::ostringstream strm;
                strm << "Image" << imageNumber << ".dcm";
                dicom->Save(strm.str());
                dicom->Store(dicomNetwork);
            }
        }
    }
}

int GERecon::ImageNumber2D(const int phase, const int slice, const int echo, const Legacy::PfilePointer& pfile)
{
    const int slicesPerPhase = pfile->SliceCount() * pfile->EchoCount();
    const int imageNumber = phase * slicesPerPhase + slice * pfile->EchoCount() + echo;

    return imageNumber;
}
