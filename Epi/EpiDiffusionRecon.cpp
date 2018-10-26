// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <iostream>
#include <string>

#include <boost/filesystem.hpp>

#include <Dicom/MR/Image.h>
#include <Dicom/MR/ImageModule.h>
#include <Dicom/MR/PrivateParameterModule.h>
#include <Dicom/Entity.h>
#include <Dicom/Instance.h>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Asset/Calibration.h>
#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>

#include <Orchestra/Cartesian2D/Homodyne.h>
#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian2D/Pocs.h>
#include <Orchestra/Cartesian2D/RampSamplingKernel.h>
#include <Orchestra/Cartesian2D/RampSamplingPlugin.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ImageType.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Common/SliceOrientation.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Epi/Diffusion/ComplexImageNex.h>
#include <Orchestra/Epi/Diffusion/HoecReconCorrection.h>
#include <Orchestra/Epi/Diffusion/TensorVectors.h>

#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/LxControlSource.h>
#include <Orchestra/Epi/VariableRampSamplingKernel.h>
#include <Orchestra/Epi/VariableRampSamplingPlugin.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>
#include <Orchestra/Legacy/PrivateDicomUtils.h>
#include <Orchestra/Legacy/LegacyImageDB.h>

#include "EpiDiffusionRecon.h"
#include "EpiReferenceScanRecon.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Legacy;
using namespace GERecon::Epi;
using namespace GERecon::Epi::Diffusion;

bool GERecon::Epi::EpiDiffusionRecon(const GERecon::Legacy::PfilePointer& pfile) // parasoft-suppress  METRICS-22 "Procedural rehearsal code handles many cases in a single int main pipeline"
{
    // Setup a pfile object. Also, extract the Pfile directory from the full Pfile path passed
    // into this function. The Pfile directory will be used as a location to save DICOM images.
    // It will also be used to determine if a vrgf.dat or vrgf_kernels.dat file can be found in
    // the Pfile directory
    const boost::filesystem::path pFileDirectory = pfile->File().parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "CppRehearsalDicoms"; // parasoft-suppress OPT-20 "Procedural, DICOM directory determined later"

    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfile->CreateOrchestraProcessingControl<Epi::LxControlSource>();
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const bool halfNexRecon = processingControl->ValueStrict<bool>("HalfNex") && !processingControl->ValueStrict<bool>("ZeroFillPartialKy");
    const int kissoffViews = processingControl->Value<int>("TransformKissoffViews") + processingControl->Value<int>("FinalImageKissoffViews");
    const bool assetPhaseCorrectionOptimization = processingControl->Value<bool>("AssetPhaseCorrectionOptimizationEnabled");
    const int numDiffusionDirections = processingControl->Value<int>("NumDiffusionDirections");
    const int numT2Images = processingControl->Value<int>("NumT2ImagesPerLocation");
    const int numBValues = processingControl->Value<int>("NumBValues");
    const int t2Nex = processingControl->Value<int>("DiffusionT2NumberOfNexes");
    const std::vector<int> diffusionNexTable = processingControl->Value<std::vector<int> >("DiffusionNexTable"); // parasoft-suppress OPT-20 "Procedural, just for rehearsal"
    const LxControlSource::NexType nexType = processingControl->Value<LxControlSource::NexType>("DiffusionNexType");
    const int vrgfInterpolatedXRes = processingControl->Value<int>("VrgfInterpolatedXRes");
    const int seriesNumber = processingControl->Value<int>("SeriesNumber");
    const int centerToEdge = processingControl->Value<int>("HomodyneCenter");
    const int transitionWidth = processingControl->Value<int>("HomodyneTransitionWidth");
    const int windowWidth = processingControl->Value<int>("HomodyneWindowWidth");
    const int iterations = processingControl->Value<int>("HomodyneIterations");
    const int assetAcceleration = processingControl->Value<int>("PhaseAssetAcceleration");
    const int numPhases = pfile->PhaseCount();
    const int numSlices = pfile->SliceCount();

    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(seriesNumber*100+3, pfile); // parasoft-suppress OPT-20 "Procedural, DICOM series number determined in process"

    // Initialize vectors of B-Values and diffusion direction annotation values
    // This information is not stored in the Pfile and is thus set to default
    // values here. If users would like to change these diffusion-specific
    // DICOM fields, they can alter the values held in these vectors. These
    // values are passed to the FillDiffusionDicomFields function below.
    std::vector<int> bValues(numBValues);
    for(size_t bValueIndex = 0; bValueIndex < bValues.size(); ++bValueIndex)
    {
        bValues[bValueIndex] = 0;
    }

    // Each diffusion direction is acquired as its own temporal phase in a scan. The diffusion
    // direction information is not stored in the Pfile header. Thus, for the diffusion direction
    // annotation to be correct, a user must specify the diffusion direction annotation for
    // each temporal phase of the scan. By default, the diffusion direction annotation field
    // (0x0043,0x1030) is set to the 'Dir1' annotation. To change this, a user can alter
    // the values held in this diffusionImageTypesPerPhase vector.
    std::vector<GERecon::DiffusionImageType> diffusionImageTypesPerPhase(numBValues * numDiffusionDirections + numT2Images);
    for(size_t phaseIndex = 0; phaseIndex < diffusionImageTypesPerPhase.size(); ++phaseIndex)
    {
        // Possible values are:
        //      DiffusionRightLeftImage
        //      DiffusionAnteriorPosteriorImage
        //      DiffusionSuperiorInferiorImage
        //      DiffusionT2Image
        //      DiffusionCombinedImage
        //      DiffusionDtiImage
        //      DiffusionDirection1Image
        //      DiffusionDirection2Image
        //      DiffusionDirection3Image
        //      DiffusionDirection4Image

        diffusionImageTypesPerPhase[phaseIndex] = GERecon::DiffusionDirection1Image;
    }

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

    ComplexFloatCube halfNexHighPassData;  //All-pass data for POCS
    ComplexFloatCube halfNexLowPassData;
    ComplexFloatCube aliasedData;
    ComplexFloatMatrix unaliasedHighPassData;
    ComplexFloatMatrix unaliasedLowPassData;
    const bool pocsEnabled = processingControl->ValueStrict<bool>("Pocs");
    boost::shared_ptr<GERecon::Cartesian2D::Pocs> pocs;
    if (isAssetScan)
    {
        // Scans with Partial Fourier reconstruction enabled (homodyne or pocs)
        if(halfNexRecon)
        {
            halfNexHighPassData.resize(reconXRes, reconYRes, numChannels);
            halfNexLowPassData.resize(reconXRes, reconYRes, numChannels);
            unaliasedHighPassData.resize(reconXRes, reconYRes);
            unaliasedLowPassData.resize(reconXRes, reconYRes);
        }
        // Full nex or zero-filled scans
        else
        {
            aliasedData.resize(reconXRes, reconYRes, numChannels);
        }
    }

    // Match product scaling - for partial ky scans apply the following scalar
    float finalImageScaler = 1.0f;
    if(halfNexRecon)
    {
        finalImageScaler = 256.0f / (reconXRes * reconYRes);
    }

    // Initialize POCS Partial Fourier reconstruction with ASSET (if needed)
    if(pocsEnabled && isAssetScan)
    {
        const int yResAfterAsset = acquiredYRes * assetAcceleration;
        const bool generateHighPassFilter = false;  // Pocs doesn't use HPF
        pocs.reset(new Cartesian2D::Pocs(Cartesian2D::PartialFourier::HalfNex, centerToEdge, transitionWidth, windowWidth, reconYRes, iterations));  
        pocs->Setup(reconXRes, yResAfterAsset, imageXRes, imageYRes, generateHighPassFilter);
    }
   
    // Create matrices needed for intermediate reconstruction data
    ComplexFloatMatrix transformedData(reconXRes, reconYRes);
    ComplexFloatMatrix accumulatedChannels(reconXRes, reconYRes);  

    // Setup transformer
    Cartesian2D::KSpaceTransformer transformer(*processingControl);

    // Setup gradwarp
    GradwarpPlugin gwPlugin(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);

    // Setup HOEC Correction
    HoecReconCorrection realtimeFieldAdjustment(*processingControl);

    // Setup complex image nex
    ComplexImageNex complexImageNex(*processingControl);

    // If tensor.dat exists in the Pfile directory and the number of diffusions
    // directions is greater than or equal to 6 then initialize a TensorVectors
    // object to read the tensor vector components from tensor.dat.
    const int tensorFileNumber = processingControl->Value<int>("TensorFileNumber");
    boost::shared_ptr<TensorVectors> tensorVectorsPtr;
    if(numDiffusionDirections >= 6)
    {
        tensorVectorsPtr = boost::make_shared<TensorVectors>(numDiffusionDirections, tensorFileNumber);
    }
    
    // Determine if vrgf.dat or vrgf_kernels.dat is accessible in the Pfile directory.
    // If neither file can be found then the ramp sampling plugin is left uninitialized
    // and ramp sampling is not performed
    boost::shared_ptr<GERecon::RampSamplingPlugin> rampSamplingPlugin;
    const boost::filesystem::path vrgfDotDatPath = pFileDirectory / "vrgf.dat";
    if(boost::filesystem::exists(vrgfDotDatPath))
    {
        boost::shared_ptr<GERecon::RampSamplingKernel> rampSamplingKernel = boost::make_shared<GERecon::RampSamplingKernel>(*processingControl);
        rampSamplingPlugin = boost::make_shared<GERecon::RampSamplingPlugin>(rampSamplingKernel, *processingControl);
    }
    
    // Create a sum of squares channel combiner object.  We will insert channel data into this channel combiner
    // as the data is processed.  After processing all channel data for a particular slice we will retrieve
    // the combined channel data from this channel object (the sum of squares is done for us under the hood).
    SumOfSquares channelCombiner(processingControl->ValueStrict<FloatVector>("ChannelWeights")); 
    const SliceInfoTable sliceTable = processingControl->ValueStrict<SliceInfoTable>("SliceTable");
    
    // Create the Reference Data object for Static Phase Correction
    // If this is an integrated reference scan diffusion scan then run
    // the reference scan pipeline to generate a ref.h5 file in the Pfile
    // directory prior to initializing the phase correction reference file
    boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin> nnpcPluginPtr;
    boost::shared_ptr<PhaseCorrectionReferenceFile> phaseCorrectionReferenceFile;
    FloatVector linearCoefficients(acquiredYRes);
    FloatVector constantCoefficients(acquiredYRes);
    nnpcPluginPtr = boost::make_shared<Epi::NearestNeighborStaticPhaseCorrectionPlugin>(*processingControl);
    int firstImageAcquisitionPhase = 0;
    if(processingControl->Value<bool>("IntegratedReferenceScan"))
    {
        EpiReferenceScanRecon(pfile);

        // First phase was a reference phase, set firstImagePhase
        // to 1 so we start with the first image phase when 
        // reconstructing images below
        firstImageAcquisitionPhase = 1;
    }
    phaseCorrectionReferenceFile = boost::make_shared<GERecon::Epi::PhaseCorrectionReferenceFile>(*processingControl, PhaseCorrectionReferenceFile::ReadMode);

    int imageNumber = 0;
    int currentBValueIndex = -1;
    int currentDirectionIndex = -1;
    int currentNexCount = 0;
    int complexNexCount = 0;
    int magnitudeNexCount = 0;
    
    for(int phase = firstImageAcquisitionPhase; phase < numPhases; ++phase)
    {
        // Determine T2 or BValue/Direction Indices
        // The following code is for product EPI diffusion sequences.
        // Product sequences acquire the diffusion directions in the
        // following order (for X number of BValues and Y number directions
        // per BValue):
        //    Reference Phase (if integrated ref scan)
        //    T2 Phase(s)
        //     B0, Dir0
        //     B0, Dir1
        //     ...
        //     BX, Dir0
        //     BX, DirY
        // The BValue and direction indices may be used to generate combined
        // images for each b
        const bool isT2Phase = phase < (numT2Images + firstImageAcquisitionPhase);
        if(isT2Phase)
        {
            currentBValueIndex = -1;
            currentDirectionIndex = -1;
            currentNexCount = t2Nex > 0 ? t2Nex: 1;
        }
        else            
        {
            currentBValueIndex = (phase - firstImageAcquisitionPhase - numT2Images) / numDiffusionDirections;
            currentDirectionIndex = (phase - firstImageAcquisitionPhase - numT2Images) % numDiffusionDirections;
            currentNexCount = diffusionNexTable[currentBValueIndex] > 0 ? diffusionNexTable[currentBValueIndex] : 1;
        }

        // If complex image nex'ing is enabled, accumulate nex'd data on
        // a channel by channel basis. The output will be an 
        // accumulated kSpace matrix for each channel. This accumulated
        // channel data is used by the remaining recon steps as if 
        // this were a single nex scan.
        ComplexFloatCube complexNexChannelData;
        if(nexType == GERecon::Epi::LxControlSource::ComplexImageSpace)
        {
            complexNexCount = currentNexCount;
            magnitudeNexCount = 1;
            complexNexChannelData.resize(vrgfInterpolatedXRes, acquiredYRes, numChannels);
        }
        else
        {
            complexNexCount = 1;
            magnitudeNexCount = currentNexCount;
        }

        for(int slice = 0; slice < numSlices; ++slice)
        {
            const SliceCorners sliceCornersForAsset = sliceTable.AcquiredSliceCorners(slice);
            std::cout << "Current Slice: " << slice << std::endl;
            
            if(nexType == GERecon::Epi::LxControlSource::ComplexImageSpace)
            {
                for(int channel = 0; channel < numChannels; ++channel)
                {
                    complexImageNex.Reset();
                    for(int echo = 0; echo < complexNexCount; ++echo)
                    {
                        // Retrieve input kSpace                        
                        ComplexFloatMatrix rawData = pfile->KSpaceData<float>(slice, echo, channel, phase);
                        
                        // Apply phase correction
                        const int channelIndex = (assetPhaseCorrectionOptimization ? phaseCorrectionReferenceFile->MaximumSnrChannel(slice) : channel);
                        linearCoefficients = phaseCorrectionReferenceFile->LinearCoefficients(slice, channelIndex);
                        constantCoefficients = phaseCorrectionReferenceFile->ConstantCoefficients(slice, channelIndex);
                        nnpcPluginPtr->ApplyPhaseCorrection(rawData, linearCoefficients, constantCoefficients);

                        // Interpolate ramp sampled data
                        const ComplexFloatMatrix kSpaceToTransform = (rampSamplingPlugin ? rampSamplingPlugin->ApplyRampSampling(rawData) : rawData);
                        
                        complexImageNex.ProcessAndAccumulate(kSpaceToTransform);
                    }
                    
                    complexNexChannelData(Range::all(), Range::all(), channel) = complexImageNex.CombinedKSpace();
                }
            }

            FloatMatrix magnitudeImage(imageXRes, imageYRes);
            magnitudeImage = 0.0f;
            for(int echo = 0; echo < magnitudeNexCount; ++echo)
            {
                if(!isAssetScan)
                {
                    // Use a sum of squares channel combiner
                    channelCombiner.Reset();
                }

                for(int channel = 0; channel < numChannels; ++channel)
                {
                    ComplexFloatMatrix kSpaceToTransform;
                    if(nexType == GERecon::Epi::LxControlSource::ComplexImageSpace)
                    {
                        kSpaceToTransform.reference(complexNexChannelData(Range::all(), Range::all(), channel));
                    }
                    else
                    {
                        ComplexFloatMatrix rawData = pfile->KSpaceData<float>(slice, echo, channel, phase);

                        const int channelIndex = (assetPhaseCorrectionOptimization ? phaseCorrectionReferenceFile->MaximumSnrChannel(slice) : channel);
                        linearCoefficients = phaseCorrectionReferenceFile->LinearCoefficients(slice, channelIndex);
                        constantCoefficients = phaseCorrectionReferenceFile->ConstantCoefficients(slice, channelIndex);
                        nnpcPluginPtr->ApplyPhaseCorrection(rawData, linearCoefficients, constantCoefficients);

                        const ComplexFloatMatrix interpolatedKSpace = (rampSamplingPlugin ? rampSamplingPlugin->ApplyRampSampling(rawData) : rawData);
                        kSpaceToTransform.reference(interpolatedKSpace);
                    }

                    if(isAssetScan && !halfNexRecon)
                    {
                        ComplexFloatMatrix aliasedDataSlice = aliasedData(Range::all(), Range::all(), channel);                   
                        transformer.Apply(aliasedDataSlice, kSpaceToTransform);                    
                    }
                    else if(halfNexRecon && isAssetScan) // parasoft-suppress  METRICS-23 "Procedural rehearsal code handles many cases in a single int main pipeline"
                    {
                        ComplexFloatMatrix highPassImage = halfNexHighPassData(Range::all(), Range::all(), channel);
                        ComplexFloatMatrix lowPassImage = halfNexLowPassData(Range::all(), Range::all(), channel);
                        transformer.Apply(highPassImage, lowPassImage, kSpaceToTransform);
                    }
                    else if(!isAssetScan) // parasoft-suppress  METRICS-23 "Procedural, for-loop rehearsal code"
                    {
                        // Handles the homodyne-only, pocs-only, zerofilling, or full kx,ky cases
                        transformer.Apply(transformedData, kSpaceToTransform);
                        channelCombiner.Accumulate(transformedData, channel);
                    }
                } // End Channel Loop

                if(isAssetScan && halfNexRecon)
                {
                    asset->Unalias(unaliasedHighPassData, unaliasedLowPassData, halfNexHighPassData, halfNexLowPassData, slice, sliceCornersForAsset);
                    accumulatedChannels = unaliasedHighPassData;
                    if(pocsEnabled)
                    {
                        // Unchop input data to avoid need for FftShift
                        MDArray::Chop(accumulatedChannels);
                        MDArray::Chop(unaliasedLowPassData);

                        // Apply Pocs to unaliased all-pass data, using phase information from low-pass data
                        pocs->ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
                    }
                    else //homodyneEnabled
                    {
                        Cartesian2D::Homodyne::ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
                    }
                }
                else if(isAssetScan && !halfNexRecon)
                {   
                    asset->Unalias(accumulatedChannels, aliasedData, slice, sliceCornersForAsset);
                }
                else if(!isAssetScan) // parasoft-suppress METRICS-23 "Procedural, if-condition in rehearsal code"
                {
                    accumulatedChannels = channelCombiner.GetCombinedImage();
                }

                magnitudeImage += abs(accumulatedChannels);
            }

            // Assure nex count is not zero, and then scale
            if(magnitudeNexCount != 0)
            {
                magnitudeImage /= static_cast<float>(magnitudeNexCount);
            }

            if(kissoffViews > 0)
            {
                magnitudeImage(Range::all(), Range(0,kissoffViews-1)) = 0.0f;
                magnitudeImage(Range::all(), Range(magnitudeImage.extent(secondDim)-kissoffViews, magnitudeImage.extent(secondDim)-1)) = 0.0f;
            }

            realtimeFieldAdjustment.Apply(magnitudeImage, phase, slice);

            const SliceCorners sliceCorners = sliceTable.SliceCorners(slice);

            gwPlugin.Run(magnitudeImage, sliceCorners, slice);

            magnitudeImage *= finalImageScaler;

            const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(slice);
            FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

            Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

            ShortMatrix shortImage(rotatedImage.shape());
            shortImage = cast<short>(rotatedImage);

            std::stringstream dicomFileName;
            dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";

            const ImageCorners imageCorners(sliceCorners, sliceOrientation);
            const GEDicom::MR::ImagePointer dicomImage = dicomSeries.NewImage(shortImage, imageNumber, imageCorners);

            // Take into account the reference phase if this is an integrated reference scan
            // by subtracting the firstImageAcquisitionPhase index from the current phase index.
            const int imagePhaseIndex = phase - firstImageAcquisitionPhase;
            const int currentBValue = (isT2Phase ? 0 : bValues[currentBValueIndex]);
            FillDiffusionDicomFields(dicomImage, diffusionImageTypesPerPhase[imagePhaseIndex], currentBValue, slice, numBValues);

            // If this is not a T2 Phase and the tensor vectors object was loaded from
            // tensor.dat then fill the tensor vector values in the DICOM image.
            if(!isT2Phase && tensorVectorsPtr)
            {
                FillTensorVectorDicomFields(dicomImage, tensorVectorsPtr->TensorVector(currentDirectionIndex));
            }

            dicomImage->Save(dicomFileName.str());

            ++imageNumber;
        } // End Slice Loop
    } // End Phase Loop
    
    return true;
}

void GERecon::Epi::FillDiffusionDicomFields(const GEDicom::MR::ImagePointer& dicomImage, const GERecon::DiffusionImageType diffusionImageType, const int bValue, const int currentSliceIndex, const int numBValues)
{
    // Multi-B Scans include a b-value bias factor added to each b-value. Replicate this product behavior here.
    const int bValueBiasFactor = ((numBValues > 1) ? 1000000000 : 0);

    std::stringstream bValueDicomString;
    bValueDicomString << (bValue + bValueBiasFactor) << "\\0\\0\\0";
    dicomImage->Insert<GEDicom::IntegerString>(0x0043, 0x1039, bValueDicomString.str());
    
    dicomImage->Insert<GEDicom::SignedShort>(0x0043, 0x1030, PrivateDicomUtils::TranslateDiffusionImageTypeValue(diffusionImageType));

    // 1 based geometric index (these indices repeat for each b-value/diffusion direction in the scan)
    dicomImage->Insert<GEDicom::UnsignedInteger>(0x0020, 0x9057, static_cast<unsigned int>((currentSliceIndex+1))); 
            
    if(bValueBiasFactor > 0)
    {
        // Product adds a bias factor to the bvalue tag for scans with more than one
        // b-value. If the bValueBiasFactor is greater than zero, then add that tag
        // to the DICOM header.
        dicomImage->Insert<GEDicom::IntegerString>(0x0043, 0x107f, bValueBiasFactor);
    }
}

void GERecon::Epi::FillTensorVectorDicomFields(const GEDicom::MR::ImagePointer& dicomImage, const MDArray::FloatTriple& tensorVectorComponents)
{
    dicomImage->Insert<GEDicom::DecimalString>(0x0019, 0x10bb, tensorVectorComponents(0));
    dicomImage->Insert<GEDicom::DecimalString>(0x0019, 0x10bc, tensorVectorComponents(1));
    dicomImage->Insert<GEDicom::DecimalString>(0x0019, 0x10bd, tensorVectorComponents(2));
}
