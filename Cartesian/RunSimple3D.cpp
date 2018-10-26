// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.


#include <boost/scoped_ptr.hpp>

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

#include <Orchestra/Core/SumOfSquares.h>
#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>
#include <Orchestra/Legacy/SliceEntry.h>

#include "CommandLine.h"
#include "Rehearsal.h"

// Include this to avoid having to type fully qualified names
using namespace GERecon;
using namespace MDArray;

/**
 * Function that does a simple 3D Recon by creating a ProcessingControl object from a Pfile and
 * looping through each slice/echo/channel to reconstruct an image from the Pfile raw data.
 * 
 * Note: this will attempt to reconstruct all acquisitions associated with the specified Pfile.
 *
 * Limitations: 
 *  No Asset support
 *
 * @author Jason Darby
 */
void GERecon::RunSimple3D(const Legacy::PfilePointer& pfile)
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("Cartesian3DRehearsal");

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
    const Control::ProcessingControlPointer processingControl = pfile->CreateOrchestraProcessingControl();

    // Get basic scan information
    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int acqZRes = processingControl->Value<int>("AcquiredZRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const bool isHalfNex = processingControl->Value<bool>("HalfNex");

    const int numPasses = pfile->PassCount();
    const int numReconstructedSlices = pfile->ReconstructedSlicesPerAcq();
    const int numEchoes = pfile->EchoCount();
    const int numChannels = pfile->ChannelCount();

    const Range all = Range::all();

    boost::scoped_ptr<GERecon::Arc::SamplingPattern> samplingPattern;
    if(pfile->IsArc())
    {
        const int patternID = processingControl->Value<int>("ArcSamplingPatternID");
        const int kYPeak = processingControl->Value<int>("ArcKYPeak");
        samplingPattern.reset(new Arc::SamplingPattern(patternID, kYPeak));
    }

    // Storage for raw K-space data to be Arc processed. A zero-padded
    // and an acquired reference are created in the case of slice zipping.
    // In this case the range in Z actually acquired is dependent upon the
    // zip factor.
    Range acquiredZRange = all;

    if(numReconstructedSlices > acqZRes)
    {
        const int offset = (numReconstructedSlices - acqZRes) / 2;
        const int zStart = offset;
        const int zEnd = zStart + acqZRes - 1;

        acquiredZRange.setRange(zStart, zEnd);
    }

    MDArray::Array<std::complex<float>,4> paddedKSpace(acqXRes, acqYRes, numReconstructedSlices, numChannels);
    MDArray::Array<std::complex<float>,4> kSpace = paddedKSpace(all, all, acquiredZRange, all);
    paddedKSpace = 0;

    // Create transformers
    Cartesian2D::KSpaceTransformer transformer(*processingControl);
    Cartesian3D::ZTransformer zTransformer(*processingControl);

    // Create channel combiner object that will do the channel combining work in channel loop.
    const FloatVector& noiseValues = processingControl->ValueStrict<FloatVector>("ChannelWeights");
    SumOfSquares channelCombiner(noiseValues);

    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);

    for(int currentPass = 0; currentPass < numPasses; ++currentPass)
    {
        for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
        {
            for(int currentSlice = 0; currentSlice < acqZRes; ++currentSlice)
            {
                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    kSpace(all, all, currentSlice, currentChannel) = pfile->KSpaceData<float>(Legacy::Pfile::PassSlicePair(currentPass, currentSlice), currentEcho, currentChannel);
                }
            }

            // perform 3D Transform if Z-Encoded
            if(pfile->IsZEncoded())
            {
                // perform ARC processing if enabled
                if(pfile->IsArc())
                {
                    GERecon::Arc::Process(kSpace, *samplingPattern);
                }

                for(int channel = 0; channel < kSpace.extent(fourthDim); ++channel)
                {
                    // Filter and transform the K-Space
                    const ComplexFloatCube volume = kSpace(all, all, all, channel);
                    ComplexFloatCube paddedVolume = paddedKSpace(all, all, all, channel);
        
                    zTransformer.Apply(paddedVolume, volume);
                }
            }
        
            // Storage for centered data that will be reconstructed - thus the recon sizes.
            ComplexFloatMatrix centeredChannelData(reconXRes, reconYRes);

            for(int currentSlice = 0; currentSlice < numReconstructedSlices; ++currentSlice)
            {
                // Zero channel combiner
                channelCombiner.Reset();

                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    trace.ConsoleMsg("Processing Slice[%d of %d], Echo[%d of %d], Channel[%d of %d]", 
                                        currentSlice+1, numReconstructedSlices, currentEcho+1, numEchoes, currentChannel+1, numChannels);

                    // Extract a single channel of data from the pfile
                    const ComplexFloatMatrix channelData = paddedKSpace(all, all, currentSlice, currentChannel);

                    // Apply transformation
                    transformer.Apply(centeredChannelData, channelData);

                    // Accumulate Channel data in channel combiner.
                    channelCombiner.Accumulate(centeredChannelData, currentChannel);
                }

                // Get the combined channel image from the ChannelCombiner
                const ComplexFloatMatrix accumulatedChannels = channelCombiner.GetCombinedImage();

                // Create storage for final image represented as a float matrix
                FloatMatrix magnitudeImage(accumulatedChannels.shape());

                // Convert complex data to specified image type.
                MDArray::ComplexToReal(magnitudeImage, accumulatedChannels, MDArray::MagnitudeData);

                // Get physical slice number and its orientation and corner points for current acquired slice
                const SliceCorners& sliceCorners = pfile->Corners(currentPass, currentSlice);
                const SliceOrientation& sliceOrientation = pfile->Orientation(currentPass, currentSlice);
                const unsigned int sliceNumber = pfile->Info(currentPass, currentSlice)->LegacyRealSliceNumber();

                // Perform Gradwarp
                gradwarp.Run(magnitudeImage, sliceCorners, sliceNumber);

                // Rotate/Transpose image accordingly
                FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

                // Scale the image if necessary
                if(isHalfNex)
                {
                    //Image is scaled if Homodyne is enabled
                    const float imageScale = 256.0f / (reconYRes * reconXRes);
                    rotatedImage *= imageScale;
                }

                // Clip the image
                Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

                // Create storage for final image represented as a short matrix - what is sent to host.
                ShortMatrix finalImage(rotatedImage.shape());

                // Cast final image from float to short
                finalImage = MDArray::cast<short>(rotatedImage);

                // Save Image
                const int imageNumber = ImageNumber3D(currentPass, sliceNumber, currentEcho, pfile);

                const ImageCorners imageCorners(sliceCorners, sliceOrientation);
                const std::string fileName = GenerateFileName(pfile->RunNumber(), 0, imageNumber);
                dicomSeries.SaveImage(fileName, finalImage, imageNumber, imageCorners);
            }
        }
    }
}

// Helper function to determine image number
int GERecon::ImageNumber3D(const int pass, const int slice, const int echo, const Legacy::PfilePointer& pfile)
{
    // Need to map the legacy "pass" number to a phase number
    const int numPassesPerPhase = pfile->PassCount() / pfile->PhaseCount();
    const int phase = pass / numPassesPerPhase;

    const int numSlicesPerPass = pfile->ReconstructedSlicesPerAcq();
    const int numEchoes = pfile->EchoCount();
    const int slicesPerPhase = numSlicesPerPass * numPassesPerPhase * numEchoes;
    
    const int imageNumber = phase * slicesPerPhase + slice * numEchoes + echo;

    return imageNumber;
}
