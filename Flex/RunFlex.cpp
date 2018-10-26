// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <boost/thread.hpp>

#include <MDArray/Utils.h>

#include <Dicom/MR/Image.h>
#include <Dicom/Network.h>

#include <Orchestra/Arc/Arc.h>
#include <Orchestra/Arc/SamplingPattern.h>

#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian3D/ZTransformer.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/SliceOrientation.h>
#include <Orchestra/Common/SliceCorners.h>
#include <Orchestra/Common/ReconTrace.h>

#include <Orchestra/Control/ProcessingControl.h>
#include <Orchestra/Control/CommonTypes.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Flex/Flex.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CommandLine.h"
#include "Rehearsal.h"

// Include this to avoid having to type fully qualified names
using namespace GERecon;
using namespace MDArray;

/**
 * Function that does a simple 3D ARC Recon by creating a ProcessingControl object from a Pfile and
 * looping through each slice/echo/channel to reconstruct an image from the Pfile raw data.
 *
 * After that, water-fat separation is performed using the complex in-phase and out-of-phase echo images.
 * 
 * This function is more memory intensive than the RunFlex() function, because it takes 5D array as input.
 * But it is simpler to use (only one function call to do the 2-pt Dixon water-fat separation)
 * 
 * @author Kang Wang
 */
void GERecon::RunFlexSimple()
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("FlexRehearsal");

    // Read Pfile from command line
    const boost::filesystem::path pfilePath = CommandLine::PfilePath();
    const Legacy::PfilePointer pfile = Legacy::Pfile::Create(pfilePath, Legacy::Pfile::AllAvailableAcquisitions, AnonymizationPolicy(AnonymizationPolicy::None));

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
    const Control::ConstProcessingControlPointer processingControl = Flex::CreateProcessingControlWithFlexParamsPopulated(*pfile);

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

    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int acqZRes = processingControl->Value<int>("AcquiredZRes");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");

    // Pulled directly from Pfile API and not ProcessingControl since the pfile class knows exactly
    // how many slices and acquisitions it has loaded and can produce raw data for.
    const int numSlices = pfile->SliceCount();
    const int numEchoes = pfile->EchoCount();
    const int numChannels = pfile->ChannelCount();

    const bool arcEnabled = processingControl->Value<bool>("ArcScan");
    
    const Range all = Range::all();

    // Storage for raw K-space data to be Arc processed. A zero-padded
    // and an acquired reference are created in the case of slice zipping.
    // In this case the range in Z actually acquired is dependent upon the
    // zip factor.
    Range acquiredZRange = all;

    if(numSlices > acqZRes)
    {
        const int offset = (numSlices - acqZRes) / 2;
        const int zStart = offset;
        const int zEnd = zStart + acqZRes - 1;

        acquiredZRange.setRange(zStart, zEnd);
    }

    // Create a transformer object to complete the K-Space to image-space transformation.
    // This will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformer(*processingControl);

    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);   

    // Storage for transformed image data, thus the image sizes.
    // If memory limited, try the non-simple RunFlex() function
    Array<std::complex<float>,5> imageDataAll(imageXRes, imageYRes, numSlices, numEchoes, numChannels);

    // Assume single phase for this release
    const int currentPhase = 0;
    for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
    {
        Array<std::complex<float>,4> paddedKSpace(acqXRes, acqYRes, numSlices, numChannels);
        Array<std::complex<float>,4> acqKSpace = paddedKSpace(all, all, acquiredZRange, all);
        paddedKSpace = 0;

        if(!pfile->IsZEncoded())
        {
            // Load all channel data (not z-encoded) for this echo
            for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
            {
                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    paddedKSpace(all, all, currentSlice, currentChannel) = pfile->KSpaceData<float>(Legacy::Pfile::PassSlicePair(currentPhase, currentSlice), currentEcho, currentChannel);
                }
            }
        }
        else
        {
            // Load all channel data (z-encoded) for this echo
            for(int currentSlice = 0; currentSlice < acqZRes; ++currentSlice)
            {
                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    acqKSpace(all, all, currentSlice, currentChannel) = pfile->KSpaceData<float>(Legacy::Pfile::PassSlicePair(currentPhase, currentSlice), currentEcho, currentChannel);
                }
            }

            // Perform Arc and z-FFT if Arc is enabled
            if (arcEnabled)
            {
                // Determine if the kacq file was specified. If not (or its
                // not valid) use the best-guess solution.
                const boost::optional<std::string> kacq = CommandLine::KacqFile();
                const bool distribute = true;
            
                trace.ConsoleMsg("Running Arc for Echo[%d of %d]", currentEcho+1, numEchoes);

                if(kacq)
                {
                    const Arc::SamplingPattern samplingPattern(*kacq);
                    GERecon::Arc::Process(acqKSpace, samplingPattern, distribute);
                }
                else
                {
                    trace.ConsoleMsg("No valid kacq file specified. Using best-guess processing.");

                    const int minTEPeak = processingControl->Value<int>("MinTEPeak");
                    GERecon::Arc::Process(acqKSpace, minTEPeak, distribute);
                }
            }

            // Z-Transform and Filter.
            Cartesian3D::ZTransformer zTransformer(*processingControl);
            trace.ConsoleMsg("Running Z-Transform and Filter for Echo[%d of %d]", currentEcho+1, numEchoes);

            for(int channel = 0; channel < numChannels; ++channel)
            {
                // Filter and transform the K-Space
                const ComplexFloatCube volume = acqKSpace(all, all, all, channel);
                ComplexFloatCube paddedVolume = paddedKSpace(all, all, all, channel);

                zTransformer.Apply(paddedVolume, volume);
            }
        }

        // Perform in-plane FFT
        for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
        {
            for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
            {
                // Extract a single channel of data from the Pfile
                const ComplexFloatMatrix kSpace = paddedKSpace(all, all, currentSlice, currentChannel);

                // Transform to image space. Data will be zipped from acquired size.
                ComplexFloatMatrix imageData = imageDataAll(all, all, currentSlice, currentEcho, currentChannel);
                transformer.Apply(imageData, kSpace);
            }
        }
    }

    // Run water-fat separation
    const int numSpecies = 2;
    Array<float, 4> imageDataWaterFat(imageXRes, imageYRes, numSlices, numSpecies);
    FloatCube imageDataWater = imageDataWaterFat(all, all, all, 0);
    FloatCube imageDataFat = imageDataWaterFat(all, all, all, 1);
    trace.ConsoleMsg("Running Flex water-fat separation");
    Flex::GenerateWaterFat(imageDataWater, imageDataFat, imageDataAll, *pfile);

    GenerateWaterFatDicom(imageDataWaterFat, pfile, gradwarp);
}

/**
 * Function that does a simple 3D ARC Recon by creating a ProcessingControl object from a Pfile.
 * 
 * It first generates field map with a slice-loop. Then it processes the field map.
 *
 * After that, water-fat separatation is performed with a slice loop,
 * using the processed field map and raw k-space data.
 *
 * @author Kang Wang
 */
void GERecon::RunFlex()
{
    GERecon::Trace trace("FlexRehearsal");

    // Read Pfile from command line
    const boost::filesystem::path pfilePath = CommandLine::PfilePath();
    const Legacy::PfilePointer pfile = Legacy::Pfile::Create(pfilePath, Legacy::Pfile::AllAvailableAcquisitions, AnonymizationPolicy(AnonymizationPolicy::None));

    const Control::ConstProcessingControlPointer processingControl = Flex::CreateProcessingControlWithFlexParamsPopulated(*pfile);

    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int acqZRes = processingControl->Value<int>("AcquiredZRes");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");

    // Pulled directly from Pfile API and not ProcessingControl since the pfile class knows exactly
    // how many slices and acquisitions it has loaded and can produce raw data for.
    const int numSlices = pfile->SliceCount();
    const int numEchoes = pfile->EchoCount();
    const int numChannels = pfile->ChannelCount();

    const bool arcEnabled = processingControl->Value<bool>("ArcScan");

    const Range all = Range::all();

    // Storage for raw K-space data to be Arc processed. A zero-padded
    // and an acquired reference are created in the case of slice zipping.
    // In this case the range in Z actually acquired is dependent upon the
    // zip factor.
    Range acquiredZRange = all;

    if(numSlices > acqZRes)
    {
        const int offset = (numSlices - acqZRes) / 2;
        const int zStart = offset;
        const int zEnd = zStart + acqZRes - 1;

        acquiredZRange.setRange(zStart, zEnd);
    }
    // This 5D array should be a few times smaller than the big 5D array in the RunFlexSimple() function,
    // because it is not zipped in-plane
    Array<std::complex<float>,5> paddedKSpace(acqXRes, acqYRes, numSlices, numEchoes, numChannels);
    Array<std::complex<float>,5> acqKSpace = paddedKSpace(all, all, acquiredZRange, all, all);
    paddedKSpace = 0;

    // Create a transformer object to complete the K-Space to image-space transformation.
    // This will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformer(*processingControl);

    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);   

    // Assume single phase for this release
    const int currentPhase = 0;
    for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
    {
        // Load all channel data (not z-encoded) for this echo
        if(!pfile->IsZEncoded())
        {
            for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
            {
                for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
                {
                    paddedKSpace(all, all, currentSlice, currentEcho, currentChannel) = pfile->KSpaceData<float>(Legacy::Pfile::PassSlicePair(currentPhase, currentSlice), currentEcho, currentChannel);
                }
            }
        }
        else
        {
            // Load all channel data (z-encoded) for this echo
            for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
            {
                for(int currentSlice = 0; currentSlice < acqZRes; ++currentSlice)
                {
                    acqKSpace(all, all, currentSlice, currentEcho, currentChannel) = pfile->KSpaceData<float>(Legacy::Pfile::PassSlicePair(currentPhase, currentSlice), currentEcho, currentChannel);
                }
            }

            // Perform Arc and z-FFT if Arc is enabled
            if (arcEnabled)
            {
                // Determine if the kacq file was specified. If not (or its
                // not valid) use the best-guess solution.
                const boost::optional<std::string> kacq = CommandLine::KacqFile();
                const bool distribute = true;
                trace.ConsoleMsg("Running Arc for Echo[%d of %d]", currentEcho+1, numEchoes);
                Array<std::complex<float>,4> imageSpaceOneEchoAllSliceChannel = acqKSpace(all, all, all, currentEcho, all);

                if(kacq)
                {
                    const Arc::SamplingPattern samplingPattern(*kacq);
                    GERecon::Arc::Process(imageSpaceOneEchoAllSliceChannel, samplingPattern, distribute);
                }
                else
                {
                    trace.ConsoleMsg("No valid kacq file specified. Using best-guess processing.");

                    const int minTEPeak = processingControl->Value<int>("MinTEPeak");
                    GERecon::Arc::Process(imageSpaceOneEchoAllSliceChannel, minTEPeak, distribute);
                }
            }
            // Z-Transform and Filter.
            Cartesian3D::ZTransformer zTransformer(*processingControl);
            trace.ConsoleMsg("Running Z-Transform and Filter for Echo[%d of %d]", currentEcho+1, numEchoes);

            for(int channel = 0; channel < numChannels; ++channel)
            {
                // Filter and transform the K-Space
                const ComplexFloatCube volume = acqKSpace(all, all, all, currentEcho, channel);
                ComplexFloatCube paddedVolume = paddedKSpace(all, all, all, currentEcho, channel);

                zTransformer.Apply(paddedVolume, volume);
            }
        }
    }

    // Generate field map, with a slice loop
    trace.ConsoleMsg("Generating field map");
    ComplexFloatCube fieldMap(imageXRes, imageYRes, numSlices);
    Array<std::complex<float>,4> imageSpaceOneSliceAllEchoChannel(imageXRes, imageYRes, numEchoes, numChannels);
    for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
    {
        for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
        {
            for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
            {
                // Extract a single channel of data from the Pfile
                const ComplexFloatMatrix kSpace = paddedKSpace(all, all, currentSlice, currentEcho, currentChannel);

                // Transform to image space. Data will be zipped from acquired size.
                ComplexFloatMatrix imageData = imageSpaceOneSliceAllEchoChannel(all, all, currentEcho, currentChannel);
                transformer.Apply(imageData, kSpace);
            }
        }

        ComplexFloatMatrix fieldMapOneSlice = fieldMap(all, all, currentSlice);
        Flex::GenerateFieldMap(fieldMapOneSlice, imageSpaceOneSliceAllEchoChannel, *processingControl);
    }

    // Process field map
    trace.ConsoleMsg("Processing field map");
    Flex::ProcessFieldMap(fieldMap, *processingControl);

    // Generate species using the processed field map, with a slice loop
    trace.ConsoleMsg("Separating species");
    const int numSpecies = 2;
    Array<float, 4> imageDataWaterFat(imageXRes, imageYRes, numSlices, numSpecies);
    FloatCube imageDataWater = imageDataWaterFat(all, all, all, 0);
    FloatCube imageDataFat = imageDataWaterFat(all, all, all, 1);
    for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
    {
        for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
        {
            for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
            {
                // Extract a single channel of data from the Pfile
                const ComplexFloatMatrix kSpace = paddedKSpace(all, all, currentSlice, currentEcho, currentChannel);

                // Transform to image space. Data will be zipped from acquired size.
                ComplexFloatMatrix imageData = imageSpaceOneSliceAllEchoChannel(all, all, currentEcho, currentChannel);
                transformer.Apply(imageData, kSpace);
            }
        }

        const ComplexFloatMatrix fieldMapOneSlice = fieldMap(all, all, currentSlice);
        FloatMatrix waterImage = imageDataWater(all, all, currentSlice);
        FloatMatrix fatImage = imageDataFat(all, all, currentSlice);
        Flex::Separate(waterImage, fatImage, fieldMapOneSlice, imageSpaceOneSliceAllEchoChannel, *processingControl);
    }

    trace.ConsoleMsg("Identify water and fat");
    Flex::IdentifyWaterFat(imageDataWater, imageDataFat, *processingControl);

    trace.ConsoleMsg("Generate water and fat DICOM");
    GenerateWaterFatDicom(imageDataWaterFat, pfile, gradwarp);
}

void GERecon::GenerateWaterFatDicom(const Array<float, 4>& imageDataWaterFat, const Legacy::PfilePointer& pfile, GradwarpPlugin& gradwarp)
{
    const Range all = Range::all();
    // Create DICOM series to save images into
    const Legacy::DicomSeries dicomSeries(pfile);

    // Get DICOM network, series number, and series description (if specified) from the command line.
    // If not specified, the optionals and pointer will be empty and no attempt will be made to insert
    // the values or store the images.
    const GEDicom::NetworkPointer dicomNetwork = CommandLine::DicomNetwork();

    const Control::ConstProcessingControlPointer processingControl = Flex::CreateProcessingControlWithFlexParamsPopulated(*pfile);
    const int originalSeriesNumber = processingControl->Value<int>("SeriesNumber");

    int baseSeriesNumber = originalSeriesNumber;
    int specifiedSeriesNumber = -1;
    if( CommandLine::SeriesNumber(specifiedSeriesNumber) )
    {
        baseSeriesNumber = specifiedSeriesNumber;
    }

    const boost::optional<std::string> seriesDescription = CommandLine::SeriesDescription();

    const int numSlices = imageDataWaterFat.extent(thirdDim);
    const int numSpecies = imageDataWaterFat.extent(fourthDim);

    for(int currentEcho = 0; currentEcho < numSpecies; ++currentEcho)
    {
        const int seriesNumberSpecies = (currentEcho == 1) ? baseSeriesNumber * 100 : baseSeriesNumber;
        for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
        {
            FloatMatrix magnitudeImage = imageDataWaterFat(all, all, currentSlice, currentEcho);

            // Get information for current slice
            const SliceOrientation& sliceOrientation = pfile->Orientation(currentSlice);
            const SliceCorners& sliceCorners = pfile->Corners(currentSlice);
            const SliceCorners& acquiredCorners = pfile->AcquiredCorners(currentSlice);

            // Perform Gradwarp
            gradwarp.Run(magnitudeImage, currentSlice, sliceCorners, acquiredCorners);

            // Rotate/Transpose image accordingly
            FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

            // Clip the image
            Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

            // Create storage for final image represented as a short matrix - form expected for DICOM.
            ShortMatrix finalImage(rotatedImage.shape());

            // Cast final image from float to short
            finalImage = MDArray::cast<short>(rotatedImage);

            // Create DICOM image
            const int imageNumber = currentSlice;
            const ImageCorners imageCorners(sliceCorners, sliceOrientation);
            const GEDicom::MR::ImagePointer dicom = dicomSeries.NewImage(finalImage, imageNumber, imageCorners);

            // Add custom annotation if specified on command line.
            // Note, boost::optional<> types only gets inserted if set.
            dicom->Insert<GEDicom::LongString>(0x0008, 0x103E, seriesDescription); // Series description, showed in image browser
            dicom->Insert<GEDicom::LongString>(0x0008, 0x1090, seriesDescription); // Manufacturer's model name, showed in annotation
            dicom->Insert<GEDicom::IntegerString>(0x0020, 0x0011, seriesNumberSpecies);

            // Save DICOM to file and also store it if network is active
            std::ostringstream strm;
            strm << "Image_S" << seriesNumberSpecies << "_I" << imageNumber << ".dcm";
            dicom->Save(strm.str());
            dicom->Store(dicomNetwork);
        }
    }
}
