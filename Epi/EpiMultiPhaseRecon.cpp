// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <iostream>
#include <string>

#include <boost/filesystem.hpp>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Arc/Arc.h>
#include <Orchestra/Arc/Calibration.h>
#include <Orchestra/Arc/CalibrationManager.h>

#include <Orchestra/Asset/Calibration.h>
#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>

#include <Orchestra/Calibration/Common/RawFile.h>
#include <Orchestra/Calibration/Common/RawFileInfo.h>

#include <Orchestra/Cartesian2D/Homodyne.h>
#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian2D/RampSamplingKernel.h>
#include <Orchestra/Cartesian2D/RampSamplingPlugin.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Common/SliceOrientation.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Epi/LxControlSource.h>
#include <Orchestra/Epi/MultibandWorker.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/RowFlipPlugin.h>
#include <Orchestra/Epi/SelfNavDynamicPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/VariableRampSamplingKernel.h>
#include <Orchestra/Epi/VariableRampSamplingPlugin.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "EpiMultiPhaseRecon.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Legacy;
using namespace GERecon::Epi;

bool GERecon::Epi::EpiMultiPhaseRecon(const GERecon::Legacy::PfilePointer& pfile) // parasoft-suppress  METRICS-22 "Procedural rehearsal code handles many cases in a single int main pipeline"
{
    // Setup a pfile object. Also, extract the Pfile directory from the full Pfile path passed
    // into this function. The Pfile directory will be used as a location to save DICOM images.
    // It will also be used to determine if a vrgf.dat or vrgf_kernels.dat file can be found in
    // the Pfile directory
    const boost::filesystem::path pFileDirectory = pfile->File().parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "CppRehearsalDicoms";

    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfile->CreateOrchestraProcessingControl<Epi::LxControlSource>();    
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int acquiredXRes = processingControl->Value<int>("AcquiredXRes");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const int extraFramesBottom = processingControl->Value<int>("ExtraFramesBottom");
    const int totalNumReferenceViews = extraFramesTop + extraFramesBottom;
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const int kissoffViews = processingControl->Value<int>("TransformKissoffViews") + processingControl->Value<int>("FinalImageKissoffViews");
    const bool assetPhaseCorrectionOptimization = processingControl->Value<bool>("AssetPhaseCorrectionOptimizationEnabled");
    const int numPhases = pfile->PhaseCount();
    const int numSlices = pfile->SliceCount();
    const bool referenceViewsFlipped = pfile->EpiReferenceViewsFlipped();
    const int numPcCoefficients = processingControl->Value<int>("NumPcCoefficients");
    
    // Initialize Asset (if needed)
    boost::shared_ptr<Asset::Worker> asset;
    boost::shared_ptr<Asset::CalibrationFile> assetCalibration;
    const bool isAssetScan = processingControl->Value<bool>("Asset");
    if(isAssetScan)
    {
        // If this is an asset scan then an asset calibration file must exist next to the Pfile and have a name of: AssetCalibration.h5
        const unsigned int examNumber = processingControl->Value<unsigned int>("ExamNumber");
        const unsigned int coilID = processingControl->Value<unsigned int>("CoilConfigUID");
        assetCalibration = boost::make_shared<Asset::CalibrationFile>(examNumber, coilID, Asset::Calibration::Regular, *processingControl);
        if( !assetCalibration->UseForReadingIfContentsMatch(GERecon::Path::InputAppData()/"AssetCalibration.h5", GERecon::Path::InputAppData()/"AssetCalibration.h5") )
        {
            throw GERecon::Exception(__SOURCE__, "Invalid Asset Calibration file!");
        }

        assetCalibration->ReadyForReading();
        asset.reset(new Asset::Worker(assetCalibration, *processingControl));
    }

    // Multiband calibration (if needed)
    const bool multibandEnabled = processingControl->Value<bool>("MultibandEnabled");
    boost::shared_ptr<GERecon::Calibration::RawFile> calibrationFile;

    if(multibandEnabled)
    {
        const GERecon::Calibration::AcquisitionType calAcquisitionType = GERecon::Calibration::ThreeD;
        calibrationFile = boost::make_shared<GERecon::Calibration::RawFile>(*processingControl,calAcquisitionType);
        // If this is a Multiband scan then a Raw calibration file must exist next to the Pfile and have a name of: RawCalibration.h5 
        if( !calibrationFile->UseForReadingIfContentsMatch(GERecon::Path::InputAppData()/"RawCalibration.h5", GERecon::Path::InputAppData()/"RawCalibration.h5") )
        {
            throw GERecon::Exception(__SOURCE__, "Invalid Raw Calibration file!");
        }
        calibrationFile->ReadyForReading();
    }
    // Set up variables for Multiband processing
    std::vector<int> unaliasedGeometricSliceNumbers;
    MDArray::Array<std::complex<float>, 4> multibandUnaliasedSliceData;
    const int numAcquiredSlice = multibandEnabled ? processingControl->ValueStrict<int>("MultibandNumAcquiredSlices"): numSlices;

    ComplexFloatCube homodyneHighPassFilteredData;
    ComplexFloatCube homodyneLowPassFilteredData;
    ComplexFloatCube aliasedData;
    ComplexFloatMatrix unaliasedHighPassData;
    ComplexFloatMatrix unaliasedLowPassData;
    if (isAssetScan)
    {
        if(homodyneEnabled)
        {
            homodyneHighPassFilteredData.resize(reconXRes, reconYRes, numChannels);
            homodyneLowPassFilteredData.resize(reconXRes, reconYRes, numChannels);
            unaliasedHighPassData.resize(reconXRes, reconYRes);
            unaliasedLowPassData.resize(reconXRes, reconYRes);
        }
        else
        {
            aliasedData.resize(reconXRes, reconYRes, numChannels);
        }
    }

    // Match product scaling - for partial ky scans apply the following scalar
    float finalImageScaler = 1.0f;
    if(homodyneEnabled)
    {
        const int xRes = processingControl->Value<int>("TransformXRes");
        const int yRes = processingControl->Value<int>("TransformYRes");
        finalImageScaler = 256.0f / (xRes * yRes);
    }

    // Create the Reference Data object for Static Phase Correction
    boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin> nnpcPluginPtr;
    boost::shared_ptr<Epi::SelfNavDynamicPhaseCorrectionPlugin> sndpcPluginPtr;
    boost::shared_ptr<PhaseCorrectionReferenceFile> phaseCorrectionReferenceFile;
    phaseCorrectionReferenceFile = boost::make_shared<GERecon::Epi::PhaseCorrectionReferenceFile>(*processingControl, PhaseCorrectionReferenceFile::ReadMode);
    
    // The following MDArrays will be used to propagate data from phase to phase for dynamic phase
    // correction. Additional descriptions follow.
    // The baselineImages and baselineReferenceViews arrays holds x,ky data for
    // all slices, all channels of the baseline phase (phase 0). This array is 
    // filled on phase 0 for each slice.
    MDArray::Array<std::complex<float>, 4> baselineImages;
    MDArray::Array<std::complex<float>, 4> baselineReferenceViews;

    // The referenceViewConstantPhase array holds the constant phase component along
    // the readout direction for each reference view. This array is automatically filled
    // by the dynamic phase correction plugin. Thus, the responsibility for this rehearsal 
    // pipeline is to allocate space for this array and pass it to the dynamic
    // phase correction plugin.
    FloatCube referenceViewConstantPhase;

    // The dynamicCoefficientDeltas array holds coefficient deltas applied on top
    // of the static coefficients computed from the reference scan. These coefficients
    // are propagated from phase to phase and automatically updated by the dynamic
    // phase correction plugin. Thus, similar to the referenceViewConstantPhase, the only
    // responsibility for this rehearsal pipeline is to allocate space for these coefficients
    // and pass this vector to the dynamic phase correction plugin.
    DynamicCoefficientsType dynamicCoefficientDeltas(numSlices);

    if(totalNumReferenceViews > 0)
    {
        // If reference views were acquired, then dynamic phase correction will be applied.
        // Initialize the dynamic phase correction plugin and allocate space for the
        // data that will be passed from phase to phase for dynamic phase correction
        // processing. These matrices are described above.
        sndpcPluginPtr = boost::make_shared<Epi::SelfNavDynamicPhaseCorrectionPlugin>(*processingControl);

        baselineImages.resize(acquiredXRes, acquiredYRes, numChannels, numSlices);
        baselineReferenceViews.resize(acquiredXRes, totalNumReferenceViews, numChannels, numSlices);
        referenceViewConstantPhase.resize(totalNumReferenceViews, numChannels, numSlices);
        referenceViewConstantPhase = 0.0f;

        for(int slice = 0; slice < numSlices; ++slice)
        {
            // Allocate space for dynamic coefficient deltas
            // and initialize to zero
            dynamicCoefficientDeltas[slice].resize(numChannels);
            for(int channel = 0; channel < numChannels; ++channel)
            {
                dynamicCoefficientDeltas[slice][channel] = 0.0f;
            }
        }
    }
    else
    {
        nnpcPluginPtr = boost::make_shared<Epi::NearestNeighborStaticPhaseCorrectionPlugin>(*processingControl);
    }
   
    // Create matrices needed for intermediate reconstruction data
    ComplexFloatMatrix transformedData(reconXRes, reconYRes);
    ComplexFloatMatrix accumulatedChannels(reconXRes, reconYRes);  

    // Setup gradwarp
    GradwarpPlugin gwPlugin(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);

    // Setup transformer
    Cartesian2D::KSpaceTransformer transformer(*processingControl);
    
    // Determine if vrgf.dat or vrgf_kernels.dat is accessible in the Pfile directory.
    // If neither file can be found then the ramp sampling plugin is left uninitialized
    // and ramp sampling is not performed
    boost::shared_ptr<GERecon::RampSamplingPlugin> rampSamplingPlugin;
    const boost::filesystem::path vrgfDotDatPath = pFileDirectory / "vrgf.dat";
    const boost::filesystem::path vrgfKernelsDotDatPath = pFileDirectory / "vrgf_kernels.dat";
    if(boost::filesystem::exists(vrgfKernelsDotDatPath))
    {
        boost::shared_ptr<GERecon::Epi::VariableRampSamplingKernel> variableRampSamplingKernel = boost::make_shared<GERecon::Epi::VariableRampSamplingKernel>(*processingControl);
        rampSamplingPlugin = boost::make_shared<GERecon::Epi::VariableRampSamplingPlugin>(variableRampSamplingKernel, *processingControl);
    }
    else if(boost::filesystem::exists(vrgfDotDatPath))
    {
        boost::shared_ptr<GERecon::RampSamplingKernel> rampSamplingKernel = boost::make_shared<GERecon::RampSamplingKernel>(*processingControl);
        rampSamplingPlugin = boost::make_shared<GERecon::RampSamplingPlugin>(rampSamplingKernel, *processingControl);
    }
    
    // Create a sum of squares channel combiner object.  We will insert channel data into this channel combiner
    // as the data is processed.  After processing all channel data for a particular slice we will retrieve
    // the combined channel data from this channel object (the sum of squares is done for us under the hood).
    SumOfSquares channelCombiner(processingControl->ValueStrict<FloatVector>("ChannelWeights")); 
    const SliceInfoTable sliceTable = processingControl->ValueStrict<SliceInfoTable>("SliceTable");

    // Allocate space to hold input data
    ComplexFloatCube rawImageData(acquiredXRes, acquiredYRes, numChannels);
    ComplexFloatCube rawReferenceData(acquiredXRes, totalNumReferenceViews, numChannels);
    FloatMatrix constantCoefficients(acquiredYRes, numChannels);
    FloatMatrix linearCoefficients(acquiredYRes, numChannels);
    FloatCube pcCoefficients(acquiredYRes, numChannels, numPcCoefficients);

    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(pfile);

    // Initialize image number to zero
    int imageNumber = 0;
    
    for(int phase = 0; phase < numPhases; ++phase)
    {
        for(int logicalSlice = 0; logicalSlice < numAcquiredSlice; ++logicalSlice)
        {
            // Determine if rowflip.param can be found in the Pfile directory
            const boost::filesystem::path rowFlipParamFile = pFileDirectory / "rowflip.param";

            const int geometricSliceIndex = multibandEnabled ? sliceTable.GeometricSliceNumber(logicalSlice) : logicalSlice;
            // Dynamic phase correction works across all channels; thus, pull
            // all channels for the current slice out of the Pfile and
            // perform phase correction
            for(int channel = 0; channel < numChannels; ++channel)
            {
                const int echo = 0;
                const Pfile::PassSlicePair sliceInfo(phase,logicalSlice);
                const ComplexFloatMatrix rawData = multibandEnabled ? pfile->KSpaceData<float>(sliceInfo, echo, channel) : pfile->KSpaceData<float>(geometricSliceIndex, echo, channel, phase);
                rawImageData(Range::all(), Range::all(), channel) = rawData(Range::all(), Range(extraFramesTop, extraFramesTop + acquiredYRes - 1));

                if(extraFramesTop > 0)
                {
                    rawReferenceData(Range::all(), Range::all(), channel) = rawData(Range::all(), Range(fromStart, extraFramesTop - 1));
                }
                else if(extraFramesBottom > 0)
                {
                    rawReferenceData(Range::all(), Range::all(), channel) = rawData(Range::all(), Range(acquiredYRes, extraFramesBottom + acquiredYRes - 1));
                }

                const int channelIndex = (assetPhaseCorrectionOptimization ? phaseCorrectionReferenceFile->MaximumSnrChannel(geometricSliceIndex) : channel);
                constantCoefficients(Range::all(), channel) = phaseCorrectionReferenceFile->ConstantCoefficients(geometricSliceIndex, channelIndex);
                linearCoefficients(Range::all(), channel) = phaseCorrectionReferenceFile->LinearCoefficients(geometricSliceIndex, channelIndex);
            }
            pcCoefficients = phaseCorrectionReferenceFile->PcCoefficients(geometricSliceIndex);

            if(totalNumReferenceViews > 0)
            {
                if(!referenceViewsFlipped && boost::filesystem::exists(rowFlipParamFile))
                {
                    const boost::shared_ptr<GERecon::Epi::RowFlipParameters> rowFlipParams = boost::make_shared<GERecon::Epi::RowFlipParameters>(*processingControl, rowFlipParamFile);
                    const boost::shared_ptr<GERecon::Epi::RowFlipPlugin> rowFlipPlugin = boost::make_shared<GERecon::Epi::RowFlipPlugin>(rowFlipParams, *processingControl);
                    for(int i = 0; i < numChannels; ++i)
                    {
                        ComplexFloatMatrix currentChannel = rawReferenceData(Range::all(), Range::all(), i);
                        rowFlipPlugin->ApplyReferenceDataRowFlip(currentChannel);
                    }
                }

                // Chop to center transforms
                rawImageData(Range(fromStart,toEnd,2), Range::all(), Range::all()) *= -1.0f;
                rawReferenceData(Range(fromStart,toEnd,2), Range::all(), Range::all()) *= -1.0f;

                Fourier::Ifft(rawImageData, MDArray::firstDim);
                Fourier::Ifft(rawReferenceData, MDArray::firstDim);

                rawImageData /= static_cast<float>(rawImageData.extent(MDArray::firstDim));
                rawReferenceData /= static_cast<float>(rawReferenceData.extent(MDArray::firstDim));

                if(phase == 0)
                {
                    // If this is the first phase, store the x,ky data from this phase
                    // as the baseline data. All future phases will compare image and
                    // reference data from the current phase with this baseline data
                    // to compute updated phase correction coefficients.
                    baselineImages(Range::all(), Range::all(), Range::all(), geometricSliceIndex) = rawImageData;
                    baselineReferenceViews(Range::all(), Range::all(), Range::all(), geometricSliceIndex) = rawReferenceData;
                }

                // If reference views were acquired, use them to apply dynamic phase correction
                // On output, the linearCoefficients and constantCoefficients arrays will contain
                // the updated coefficients that were applied to the rawImageData. Note that if
                // assetPhaseCorrectionOptimization is enabled then the same set of coefficients
                // are applied to all channels of a given slice. In this case, only the first 
                // channel index of the linearCoefficients and constantCoefficients arrays contains
                // updated coefficients.
                // The referenceViewConstantPhase and dynamicCoefficientDeltas values are 
                // propagated from phase to phase for each slice. These values are automatically
                // updated by the dynamic phase correction plugin and do not need to be modified here.
                ComplexFloatCube baselineImageData = baselineImages(Range::all(), Range::all(), Range::all(), geometricSliceIndex);
                ComplexFloatCube baselineReferenceData = baselineReferenceViews(Range::all(), Range::all(), Range::all(), geometricSliceIndex);
                FloatMatrix referenceViewConstantPhaseCurrentSlice = referenceViewConstantPhase(Range::all(), Range::all(), geometricSliceIndex);
                const int shotIndex = 0; //Only support single shot
                sndpcPluginPtr->ApplyPhaseCorrection(referenceViewConstantPhaseCurrentSlice, dynamicCoefficientDeltas[geometricSliceIndex], rawImageData, rawReferenceData, baselineImageData, baselineReferenceData, pcCoefficients, shotIndex);

                // Transform phase corrected data back to kx,ky for apodization filtering and 2D kSpace transform
                Fourier::Fft(rawImageData, MDArray::firstDim);
                rawImageData(Range(fromStart, toEnd, 2), Range::all(), Range::all()) *= -1.0f;
            }
            else
            {
                // If reference views were not acquired, then apply static phase correction using
                // coefficients from the reference scan
                for(int channel = 0; channel < numChannels; ++channel)
                {
                    const FloatVector linearCoefficientsCurrentChannel = linearCoefficients(Range::all(), channel);
                    const FloatVector constantCoefficientsCurrentChannel = constantCoefficients(Range::all(), channel);
                    ComplexFloatMatrix dataToCorrect = rawImageData(Range::all(), Range::all(), channel);
                    nnpcPluginPtr->ApplyPhaseCorrection(dataToCorrect, linearCoefficientsCurrentChannel, constantCoefficientsCurrentChannel);
                }
            }

            std::cout << "Current Phase : " << phase << " Current Logical Slice: " << logicalSlice << std::endl;
            int geometricSliceLoopcount = 0;

            if(multibandEnabled) //for Multiband scans
            {
                // Create pointers for Multiband Calibration and Multiband Worker
                GERecon::Arc::CalibrationPtr calPtr;
                GERecon::Epi::MultibandCalibrationProcessor multibandCalibrationProcessor(*processingControl);
                Epi::MultibandWorker multibandWorker(*processingControl);

                // Call the Multiband Calibration processor 
                const bool flipDataInPhaseEncodeDirection = processingControl->ValueStrict<bool>("PhaseFlip");
                calPtr = CalibrationProcessorMultiband(calibrationFile,logicalSlice,multibandCalibrationProcessor, flipDataInPhaseEncodeDirection);    

                // Apply ramp sampling interpolation and unpack in plane accelerated data
                Range all = Range::all();
                const int resampledXRes = (rampSamplingPlugin.get() == NULL) ? rawImageData.extent(thirdDim) : rampSamplingPlugin->InterpolatedXRes();
                const int inplaneUnpackedYRes = multibandWorker.InPlaneUnpackedYRes(); // yRes after putting in zero for unacquired readouts for in plane accelerated scans
                const int numChannels = rawImageData.extent(thirdDim);

                ComplexFloatCube resampledKSpaceToUnalias(resampledXRes, inplaneUnpackedYRes, numChannels);
                resampledKSpaceToUnalias = 0.0f;
                const int inPlaneAcceleration = multibandWorker.InPlaneAccelerationFactor();
                if(rampSamplingPlugin.get() != NULL)
                {
                    for(int channel = 0; channel < rawImageData.extent(thirdDim); ++channel)
                    {
                        // Ramp sampling is enabled, do the interpolation here. Only interpolate the acquired lines.
                        const ComplexFloatMatrix kSpaceToInterpolate = rawImageData(all, all, channel);
                        resampledKSpaceToUnalias(all, Range(0,toEnd,inPlaneAcceleration), channel) = rampSamplingPlugin->ApplyRampSampling(kSpaceToInterpolate);
                    }
                }
                else
                {
                    // No need to do ramp sampling just unpack data to put zeros in place for data acquired with in plane acceleration
                    resampledKSpaceToUnalias(all, Range(0,toEnd,inPlaneAcceleration), all) = rawImageData;
                }

                // Perform Multiband Unaliasing   
                multibandWorker.UnaliasSliceData(multibandUnaliasedSliceData, unaliasedGeometricSliceNumbers, 
                    resampledKSpaceToUnalias, calPtr, logicalSlice);

                geometricSliceLoopcount = boost::numeric_cast<int>(unaliasedGeometricSliceNumbers.size());

            }
            else // for Non-Multiband scans
            {
                geometricSliceLoopcount = 1;
            }

            for(int unaliasedSliceIndex = 0; unaliasedSliceIndex < geometricSliceLoopcount; ++unaliasedSliceIndex)
            {
                const int totalgeometricslices = sliceTable.GeometricSliceLocations();
                const int geometricSliceNumber = multibandEnabled ? unaliasedGeometricSliceNumbers[unaliasedSliceIndex]: logicalSlice;

                if(geometricSliceNumber < totalgeometricslices) // Discarding extended slices for Multiband scans
                {
                    const SliceCorners sliceCornersForAsset = sliceTable.AcquiredSliceCorners(geometricSliceNumber);

                    if(!isAssetScan)
                    {
                        // Use a sum of squares channel combiner
                        channelCombiner.Reset();
                    }

                    for(int channel = 0; channel < numChannels; ++channel)
                    {
                        const ComplexFloatMatrix kSpaceToTransform = multibandEnabled ? multibandUnaliasedSliceData(Range::all(), Range::all(), unaliasedSliceIndex, channel): (rampSamplingPlugin ? rampSamplingPlugin->ApplyRampSampling(rawImageData(Range::all(), Range::all(), channel)) : rawImageData(Range::all(), Range::all(), channel));

                        if(isAssetScan && !homodyneEnabled) // parasoft-suppress  METRICS-23 "rehearsal code handles many cases in a single int main pipeline"
                        {
                            ComplexFloatMatrix aliasedDataSlice = aliasedData(Range::all(), Range::all(), channel);
                            transformer.Apply(aliasedDataSlice, kSpaceToTransform);
                        }
                        else if(homodyneEnabled && isAssetScan) // parasoft-suppress  METRICS-23 "Procedural rehearsal code handles many cases in a single int main pipeline"
                        {
                            ComplexFloatMatrix highPassImage = homodyneHighPassFilteredData(Range::all(), Range::all(), channel);
                            ComplexFloatMatrix lowPassImage = homodyneLowPassFilteredData(Range::all(), Range::all(), channel);
                            transformer.Apply(highPassImage, lowPassImage, kSpaceToTransform);
                        }
                        else if(!isAssetScan) // parasoft-suppress  METRICS-23 "Procedural, for-loop rehearsal code"
                        {
                            // Handles the homodyne only, zerofilling, or full kx,ky cases
                            transformer.Apply(transformedData, kSpaceToTransform);
                            channelCombiner.Accumulate(transformedData, channel);
                        }

                    } // End Channel Loop

                    if (isAssetScan && homodyneEnabled)
                    {
                        asset->Unalias(unaliasedHighPassData, unaliasedLowPassData, homodyneHighPassFilteredData, homodyneLowPassFilteredData, geometricSliceNumber, sliceCornersForAsset);
                        accumulatedChannels = unaliasedHighPassData;
                        Cartesian2D::Homodyne::ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
                    }
                    else if(isAssetScan && !homodyneEnabled) // parasoft-suppress  METRICS-23 "rehearsal code handles many cases in a single int main pipeline"
                    {   
                        asset->Unalias(accumulatedChannels, aliasedData, geometricSliceNumber, sliceCornersForAsset);
                    }
                    else if(!isAssetScan)
                    {
                        accumulatedChannels = channelCombiner.GetCombinedImage();
                    }

                    // Create a magnitude image
                    FloatMatrix finalImage(accumulatedChannels.shape());
                    MDArray::ComplexToReal(finalImage, accumulatedChannels, MDArray::MagnitudeData );

                    if(kissoffViews > 0)
                    {
                        finalImage(Range::all(), Range(fromStart,kissoffViews-1)) = 0.0f;
                        finalImage(Range::all(), Range(finalImage.extent(secondDim)-kissoffViews, finalImage.extent(secondDim)-1)) = 0.0f;
                    }

                    const SliceCorners sliceCorners = sliceTable.SliceCorners(geometricSliceNumber);
                    gwPlugin.Run(finalImage, sliceCorners, geometricSliceNumber);

                    finalImage *= finalImageScaler;

                    const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(geometricSliceNumber);
                    FloatMatrix rotatedImage = RotateTranspose::Apply<float>(finalImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

                    Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

                    ShortMatrix shortImage(rotatedImage.shape());
                    shortImage = cast<short>(rotatedImage);

                    imageNumber = phase * totalgeometricslices + geometricSliceNumber;

                    std::stringstream dicomFileName;
                    dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(5) << std::setfill('0') << imageNumber << ".dcm";
                    const ImageCorners imageCorners(sliceCorners, sliceOrientation);
                    dicomSeries.SaveImage(dicomFileName.str(), shortImage, imageNumber, imageCorners);

                }
            } // End Geometric Slice Loop
        } // End Logical  Loop
    } // End Phase Loop
    
    return true;
}

GERecon::Arc::CalibrationPtr GERecon::Epi::CalibrationProcessorMultiband(const boost::shared_ptr<GERecon::Calibration::RawFile>& calibrationfile,int slice, GERecon::Epi::MultibandCalibrationProcessor& multibandCalibrationProcessor, const bool flipDataInPhaseEncodeDirection)
{
    // Prepare a single slice of calibration data and its associated parameters so that we can call Calibration() later
    const float calScanCenter = calibrationfile->Info().ScanCenter();
    const MDArray::Array<std::complex<float>, 4> imageSpaceToCalibrateWith = multibandCalibrationProcessor.CalibrationImageSpace(calibrationfile->ImageSpaceData(GERecon::Calibration::SurfaceImageSpaceAllPassSets), 
                                                                                                                              calibrationfile->Info().FirstSliceCorners(),
                                                                                                                              calibrationfile->Info().LastSliceCorners(),
                                                                                                                              calScanCenter,
                                                                                                                              slice,
                                                                                                                              flipDataInPhaseEncodeDirection);
    
    const MDArray::Array<std::complex<float>,4>& kSpaceToCalibrateWith = multibandCalibrationProcessor.CalibrationKSpace(imageSpaceToCalibrateWith);

    const Arc::EncodesVector calAcquiredLocations = multibandCalibrationProcessor.CalAcquiredLocations();
    const Arc::EncodesVector undersampledDataAcquiredLocations = multibandCalibrationProcessor.UndersampledDataAcquiredLocations();
    const MDArray::TinyVector<int, 4>& undersampledDataDimensionSizes = multibandCalibrationProcessor.UndersampledDataDimensionSizes();
    const int yAcceleration = multibandCalibrationProcessor.OverallYAcceleration();
    const int zAcceleration = 1;
    const int xPeakPosition = multibandCalibrationProcessor.XPeakPosition();
    const int yPeakPosition = multibandCalibrationProcessor.YPeakPosition();
    const int zPeakPosition = 0;
    const int calibrationIndex = slice;
    const bool distribute = false;
    const int leftBlanks = multibandCalibrationProcessor.XCalPointsToSkipLeft();
    const int rightBlanks = multibandCalibrationProcessor.XCalPointsToSkipRight();
    const int unacceleratedDimKernelSize = multibandCalibrationProcessor.XKernelSize();
    const int acceleratedDimKernelSize = multibandCalibrationProcessor.YKernelSize();
    const int maxCalXSize = kSpaceToCalibrateWith.extent(MDArray::firstDim);
    const int maxCalYSize = kSpaceToCalibrateWith.extent(MDArray::secondDim);
    const int maxCalZSize = kSpaceToCalibrateWith.extent(MDArray::thirdDim);

    // The following function call will run ARC calibration and store the calibration in the calibration manager
    GERecon::Arc::CalibrationPtr calPtr = Arc::CalibrationManager::Instance()->Calibration(kSpaceToCalibrateWith,
                                                                                            calAcquiredLocations, 
                                                                                            undersampledDataAcquiredLocations,
                                                                                            undersampledDataDimensionSizes,
                                                                                            yAcceleration, 
                                                                                            zAcceleration,
                                                                                            xPeakPosition,
                                                                                            yPeakPosition,
                                                                                            zPeakPosition,
                                                                                            calibrationIndex,
                                                                                            distribute,
                                                                                            leftBlanks, 
                                                                                            rightBlanks,
                                                                                            unacceleratedDimKernelSize,
                                                                                            acceleratedDimKernelSize,
                                                                                            maxCalXSize,
                                                                                            maxCalYSize,
                                                                                            maxCalZSize);
    return calPtr;
}
