// Copyright 2018 General Electric Company.  All rights reserved.

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Hdf5/Snap.h>

#include <Orchestra/Control/CommonTypes.h>
#include <Orchestra/Control/ProcessingControl.h>

// Algorithms
#include <Orchestra/Core/PhaseCorrectAndCombine.h>
#include <Orchestra/Core/SumOfSquares.h>
#include <Orchestra/Core/Windows.h>

#include <Orchestra/Spiral/LxControlSource.h>
#include <Orchestra/Spiral/SpiralWorker.h>

#include "SpiralUtil.h"

using namespace GERecon;
using namespace GERecon::Spiral;
using namespace GERecon::Legacy;
using namespace MDArray;

void GERecon::SpiralUtil::SpiralRecon(MDArray::FloatCube& imageData, // parasoft-suppress  METRICS-22 "High line count required for algorithm"
                                      const MDArray::ComplexFloat4D& rawData,
                                      const Control::ProcessingControlPointer& processingControl,
                                      const std::string& dataDirectory,
                                      const float scalingFactor)
{
    // Create an HDF5 file to save debug data to. This debug file will be enabled if the Snap Key
    // specified in this constructor is included on the command line for this rehearsal application:
    //      --GEHdf5.Snap.Key SpiralReconRehearsalDebug
    const GEHdf5::Snap debugFile("SpiralReconRehearsal.h5", "SpiralReconRehearsalDebug");
    debugFile.Write(rawData, "rawData");

    // Create a trace object to log messages to the console for this reeharsal
    GERecon::Trace trace("SpiralReconRehearsal");
    trace.ConsoleMsg("Initiating Spiral recon");

    // Acquisition and Recon Parameters.
    const int kissoffs = processingControl->Value<int>("Kissoffs");
    const int numAcquisitions = processingControl->Value<int>("NumAcquisitions");
    const int numTotalSlices = processingControl->Value<int>("NumSlices");
    const int numSlices = (numTotalSlices + (numAcquisitions-1)) / numAcquisitions - 2 * kissoffs;
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int acqXRes = processingControl->Value<int>("NumPointsPerArm");
    const int acqYRes = processingControl->Value<int>("NumArms");
    const float fermiRadius = processingControl->Value<float>("FermiRadius");
    const float fermiWidth = processingControl->Value<float>("FermiWidth");
    const float fermiEccentricity = processingControl->Value<float>("FermiEccentricity");
    const bool chopX = true; // Always chop in x-direction
    const bool chopY = !processingControl->Value<bool>("ChoppedData"); // Note inverse logic
    const int tediff = processingControl->Value<int>("B0MapArmsTEDifference");
    const int gridFovFactor = (processingControl->Value<bool>("IsGridSizeDoubled")) ? 2 : 1;
    const int gridWeight1 = (processingControl->Value<bool>("IsSampleDensityCorrectionEnabled")) ? 1 : 0;
    const int gridWeight2 = (processingControl->Value<bool>("IsApplyRolloffEnabled")) ? 1 : 0;
    const int functionWidth = processingControl->Value<int>("KaiserBesselFunctionWidth");
    const int shapeFactor = processingControl->Value<int>("KaiserBesselShapeFactor");
    const bool isApplyRolloffEnabled = (gridWeight2 != 0);
    const bool isB0MapEnabled = processingControl->Value<bool>("IsB0MapEnabled");
    const int imageSize = processingControl->Value<int>("TransformXRes");
    const float maxGradient = processingControl->Value<float>("MaxGradientAmplitude");
    const float maxSlewRate = processingControl->Value<float>("MaxGradientSlewRate");
    const float scanFOV = processingControl->Value<float>("ScanFOV");
    float displayFOV = processingControl->Value<float>("DisplayFOV");
    displayFOV = (displayFOV > 0.0f) ? displayFOV : scanFOV;
    const float sampleTime = processingControl->Value<float>("A2DSampleTime") / 1000;    // Convert from usecs to msecs
    const int numInterleaves = (processingControl->Value<bool>("HaveTwoInterleavesPerExcitation")) ? 2 : 1;
    const int numPointsBeforeEchoCenter = processingControl->Value<int>("NumPointsBeforeEchoCenter");
    const float transitionRadiusA = processingControl->Value<float>("GriddingRadiusA");
    const float transitionRadiusB = processingControl->Value<float>("GriddingRadiusB");
    const float densityFactor = processingControl->Value<float>("DensityFactor");
    const bool needToCreateTrajectory = processingControl->Value<bool>("ShouldReconCreateTrajectory");
    const bool useComplexCombine = processingControl->Value<bool>("UseComplexCombine");
    const bool createRealImages = processingControl->Value<bool>("CreateRealImages");

    // Warning message
    if( numChannels != rawData.extent(thirdDim) )
    {
        trace.ConsoleMsg("Raw data points number per arm: expected %d and received %d", acqXRes, rawData.extent(firstDim));
    }
    if( numSlices != rawData.extent(fourthDim) )
    {
        trace.ConsoleMsg("Raw data arm number: expected %d and received %d", acqYRes, rawData.extent(secondDim));
    }
    if( numChannels != rawData.extent(thirdDim) )
    {
        trace.ConsoleMsg("Raw data channel number: expected %d and received %d", numChannels, rawData.extent(thirdDim));
    }
    if( numSlices != rawData.extent(fourthDim) )
    {
        trace.ConsoleMsg("Raw data slice number: expected %d and received %d", numSlices, rawData.extent(fourthDim));
    }

    // Create a SpiralWorker and generate trajectories of data arms
    Spiral::SpiralWorker spiral;
    FloatCube trajectoryTable( 3, acqXRes, acqYRes );

    if ( needToCreateTrajectory )
    {
        trace.ConsoleMsg("Creating Spiral trajectory");
        spiral.CreateTrajectory(
            trajectoryTable,
            imageSize,
            maxGradient,
            maxSlewRate,
            scanFOV,
            sampleTime,
            isB0MapEnabled,
            numInterleaves,
            numPointsBeforeEchoCenter,
            transitionRadiusA,
            transitionRadiusB,
            densityFactor );
    }
    else
    {
        trace.ConsoleMsg("Loading Spiral trajectory");
        const std::string trajFile = dataDirectory + "/kspace";
        spiral.ReadTrajectoryFromDisk( trajectoryTable, trajFile );
    }
    debugFile.Write(trajectoryTable, "trajectoryTable");

    // Create an interpolation function for spiral regridding
    const int gridSize = imageSize * gridFovFactor;
    const int interFunctionLength = 256;
    FloatVector interpFunction( interFunctionLength );
    FloatVector interpFunctionInverse( gridSize );
    spiral.CreateInterpolationFunctionAndInverse( interpFunction, interpFunctionInverse, static_cast<float>(shapeFactor), functionWidth );

    // Calculate scaling parameters
    const float rationalScalingFactor = processingControl->Value<float>("RationalScalingFactor");
    const float lxScalingFudgeFactor = 256.0F;
    const float rawScalingFactor = scalingFactor * lxScalingFudgeFactor * rationalScalingFactor / (gridSize * gridSize);

    float densityScaling = 1.0f;
    float coordScaling = 1.0f;
    spiral.CalculateDisplayDensityAndCoordFactors( densityScaling, coordScaling, gridFovFactor, functionWidth, shapeFactor, displayFOV, scanFOV );

    // Create a Fermi filter to apply to grid
    FloatMatrix fermiFilter(gridSize, gridSize);
    Windows::Fermi(fermiFilter, fermiRadius * gridFovFactor, fermiWidth * gridFovFactor, fermiEccentricity ); // initialize the window

    // Create phase shift filter.
    // Create half pixel shift scaled by reconstructed size, not acquisition size.
    const bool gridChopY = true;
    ComplexFloatMatrix phaseShiftFilter(fermiFilter.shape());
    Windows::ChoppingWithShift(phaseShiftFilter, chopX, gridChopY, gridSize, gridSize);

    // Create a single K-Space filter composed of the Fermi window, phase shift filter, and rational scaling factor.
    // Since they are just multiplications, do it once instead of 3 different times in the channel loop.
    ComplexFloatMatrix kSpaceFilter(fermiFilter.shape());
    kSpaceFilter = fermiFilter * phaseShiftFilter; // This is point wise multiplication - MATLAB equivalent ".*"
    debugFile.Write(kSpaceFilter, "kSpaceFilter");

    // Algorithm plugin for Channel Accumulation
    boost::scoped_ptr<ChannelCombiner> channelCombiner;
    if ( useComplexCombine )
    {
        trace.ConsoleMsg("Using PhaseCorrectAndCombine");
        channelCombiner.reset(new PhaseCorrectAndCombine(processingControl->ValueStrict<FloatVector>("ChannelWeights"), *processingControl));
    }
    else
    {
        trace.ConsoleMsg("Using SumOfSquares");
        channelCombiner.reset(new SumOfSquares(processingControl->ValueStrict<FloatVector>("ChannelWeights")));
    }
    
    // Image type    
    MDArray::DataType imageType;
    if(createRealImages)
    {
        imageType = MDArray::RealData;
    }
    else
    {
        imageType = MDArray::MagnitudeData;
    }

    std::stringstream sliceNameString;
    std::stringstream channelNameString;
    for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
    {
        sliceNameString.str("");
        sliceNameString << "slice_" << std::setw(3) << std::setfill('0') << currentSlice;

        // Zero out channel combiner buffer for the next set of channels.
        channelCombiner->Reset();

        for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
        {
            channelNameString.str("");
            channelNameString << "channel_" << std::setw(3) << std::setfill('0') << currentChannel;

            // Extract an echo for a single channel of data from the pfile
            const ComplexFloatMatrix allChannelData = rawData(Range::all(), Range::all(), currentChannel, currentSlice);
            const bool isSampleDensityCorrectionEnabled = (gridWeight1 == 1);
            float frequencyOffset = 0.0f;      // Output frequency offset (Hz/sample)
            float xSlopeOfFit = 0.0f;         // Output x and y slope of fit (Hz-cm/sample)
            float ySlopeOfFit = 0.0f;

            const int numB0Arms = ( isB0MapEnabled ) ? 2 : 0;
            if ( isB0MapEnabled )
            {
                trace.ConsoleMsg("Calculating B0 map correction coefficients: channel %d in slice %d", currentChannel, currentSlice);

                // When B0 map correction is enabled, the first two arms are used to calculate correction values.
                // Create references to the first arm of data and to the first arm in the trajectory table,
                // along with a complex matrix to receive the reconstructed (and cropped) image.
                // Go through the reconstruction process up to creation of the complex image.
                ComplexFloatMatrix firstB0Image(imageSize, imageSize);
                const ComplexFloatMatrix firstB0Data = allChannelData( Range::all(), Range(0,0) );
                const FloatCube firstB0Traj = trajectoryTable( Range::all(), Range::all(), Range(0,0) );

                // Use a dummy demodulation vector and displacement Matrix when computing B0 map corrections
                const ComplexFloatVector demodVectorDummy;
                const FloatMatrix displacementVectorDummy;

                spiral.CreateSpiralChannelImage(
                    firstB0Image,                // Output is complex channel image
                    firstB0Data,
                    firstB0Traj,
                    rawScalingFactor,
                    densityScaling,
                    coordScaling,
                    gridSize,
                    interpFunction,             // Input        256floats
                    displacementVectorDummy,
                    demodVectorDummy,
                    functionWidth,
                    kSpaceFilter,
                    isSampleDensityCorrectionEnabled,
                    SpiralWorker::NoArmChop);         // No chopping for B0 map arms

                ComplexFloatMatrix secondB0Image(imageSize, imageSize);
                const ComplexFloatMatrix secondB0Data = allChannelData( Range::all(), Range(1,1) );
                const FloatCube secondB0Traj = trajectoryTable( Range::all(), Range::all(), Range(1,1) );
                spiral.CreateSpiralChannelImage(
                    secondB0Image,                // Output is complex channel image
                    secondB0Data,
                    secondB0Traj,
                    rawScalingFactor,
                    densityScaling,
                    coordScaling,
                    gridSize,
                    interpFunction,             // Input        256floats
                    displacementVectorDummy,
                    demodVectorDummy,
                    functionWidth,
                    kSpaceFilter,
                    isSampleDensityCorrectionEnabled,
                    SpiralWorker::NoArmChop);         // No chopping for B0 map arms

                // Use the two images to determine B0 map correction values for this channel
                FloatMatrix b0Image(imageSize, imageSize);
                spiral.CalculateB0MapCorrections(b0Image, frequencyOffset, xSlopeOfFit, ySlopeOfFit, firstB0Image, secondB0Image, tediff, sampleTime, gridSize);

                debugFile.Write(b0Image, sliceNameString.str() + "/" + channelNameString.str() + "/b0Image");
            }

            // Create demodulation vector. If we have a frequency offset then fill it in 
            ComplexFloatVector demodVector;
            if ( abs(frequencyOffset) > 0.0f )
            {
                demodVector.resize( acqXRes );
                spiral.CreateDemodulationVector( demodVector, frequencyOffset );
            }

            // Create a matrix containing kspace displacements to add to each trajectory arm
            FloatMatrix displacementVector;
            if ( (abs(xSlopeOfFit) > 0.0f) || (abs(ySlopeOfFit) > 0.0f) )
            {
                displacementVector.resize( 3, acqXRes );
                spiral.CreateKSpaceDisplacementVector( displacementVector, xSlopeOfFit, ySlopeOfFit );
            }

            // Storage to receive cropped center of image data
            ComplexFloatMatrix centralChannelData(imageSize, imageSize);

            // Create a reference to this echo's data arms, skipping any B0 map arms
            const ComplexFloatMatrix channelData = allChannelData( Range::all(), Range(numB0Arms,allChannelData.extent(secondDim)-1) );
            const FloatCube dataTraj = trajectoryTable( Range::all(), Range::all(), Range(numB0Arms,trajectoryTable.extent(thirdDim)-1) );
            const SpiralWorker::ChopType armChop = (chopY) ? SpiralWorker::OddArmChop : SpiralWorker::NoArmChop;

            spiral.CreateSpiralChannelImage(
                centralChannelData,        // Output is complex channel image
                channelData,
                dataTraj,
                rawScalingFactor,
                densityScaling,
                coordScaling,
                gridSize,
                interpFunction,             // Input        256floats
                displacementVector,         // Input        3floats x points
                demodVector,
                functionWidth,
                kSpaceFilter,
                isSampleDensityCorrectionEnabled,
                armChop);

            debugFile.Write(channelData, sliceNameString.str() + "/" + channelNameString.str() + "/data");
            debugFile.Write(dataTraj, sliceNameString.str() + "/" + channelNameString.str() + "/traj");
            debugFile.Write(centralChannelData, sliceNameString.str() + "/" + channelNameString.str() + "/image");

            // Apply rolloff filter to individual channel image before channel combine
            if ( isApplyRolloffEnabled )
            {
                spiral.ApplyRolloffFilter(centralChannelData, interpFunctionInverse);
            }
            debugFile.Write(centralChannelData, sliceNameString.str() + "/" + channelNameString.str() + "/imageFiltered");

            // Accumulate Channel data in channel combiner.
            channelCombiner->Accumulate(centralChannelData, currentChannel);
        }   // end for currentChannel

        // Get the combined channel image from the ChannelCombiner
        ComplexFloatMatrix accumulatedChannels = channelCombiner->GetCombinedImage();
        debugFile.Write(accumulatedChannels, sliceNameString.str() + "/combinedComplexImage");

        // Convert complex data to specified image type.
        FloatMatrix floatImage(accumulatedChannels.shape());
        MDArray::ComplexToReal(floatImage, accumulatedChannels, imageType);
        debugFile.Write(floatImage, sliceNameString.str() + "/finalImage");

        trace.ConsoleMsg("Spiral reconstuction done in slice %d", currentSlice);

        // Arrange output image data
        imageData(Range::all(), Range::all(), currentSlice) = floatImage;
    }
    debugFile.Write(imageData, "imageData");
}
