// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <iostream>
#include <string>

#include <boost/filesystem.hpp>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Asset/Calibration.h>
#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>

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

#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/LxControlSource.h>
#include <Orchestra/Epi/VariableRampSamplingKernel.h>
#include <Orchestra/Epi/VariableRampSamplingPlugin.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "EpiRecon.h"

using namespace GERecon;
using namespace GERecon::Legacy;
using namespace MDArray;
using namespace GERecon::Epi;

bool GERecon::Epi::EpiRecon(const GERecon::Legacy::PfilePointer& pfile) // parasoft-suppress  METRICS-22 "Procedural rehearsal code handles many cases in a single int main pipeline"
{
    // Setup a pfile object. Also, extract the Pfile directory from the full Pfile path passed
    // into this function. The Pfile directory will be used as a location to save DICOM images.
    // It will also be used to determine if a vrgf.dat or vrgf_kernels.dat file can be found in
    // the Pfile directory
    const boost::filesystem::path pFileDirectory = pfile->File().parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "CppRehearsalDicoms"; // parasoft-suppress  OPT-20 "Procedural Rehearsal Code"

    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(pfile); // parasoft-suppress  OPT-20 "Procedural Rehearsal Code"

    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfile->CreateOrchestraProcessingControl<Epi::LxControlSource>();    
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");    
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const int kissoffViews = processingControl->Value<int>("TransformKissoffViews") + processingControl->Value<int>("FinalImageKissoffViews");
    const bool assetPhaseCorrectionOptimization = processingControl->Value<bool>("AssetPhaseCorrectionOptimizationEnabled");
    const int numSlices = pfile->SliceCount();

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
    boost::shared_ptr<PhaseCorrectionReferenceFile> phaseCorrectionReferenceFile;
    nnpcPluginPtr = boost::make_shared<Epi::NearestNeighborStaticPhaseCorrectionPlugin>(*processingControl);
    phaseCorrectionReferenceFile = boost::make_shared<GERecon::Epi::PhaseCorrectionReferenceFile>(*processingControl, PhaseCorrectionReferenceFile::ReadMode);
    FloatVector linearCoefficients(acquiredYRes);
    FloatVector constantCoefficients(acquiredYRes);
   
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
    
    int imageNumber = 0;
    for(int slice = 0; slice < numSlices; ++slice)
    {
        const SliceCorners sliceCornersForAsset = sliceTable.AcquiredSliceCorners(slice);
        std::cout << "Current Slice: " << slice << std::endl;

        if(!isAssetScan)
        {
            // Use a sum of squares channel combiner
            channelCombiner.Reset();
        }

        for(int channel = 0; channel < numChannels; ++channel)
        {
            const int echo = 0;
            ComplexFloatMatrix rawData = pfile->KSpaceData<float>(slice, echo, channel);

            const int channelIndex = (assetPhaseCorrectionOptimization ? phaseCorrectionReferenceFile->MaximumSnrChannel(slice) : channel);
            linearCoefficients = phaseCorrectionReferenceFile->LinearCoefficients(slice, channelIndex);
            constantCoefficients = phaseCorrectionReferenceFile->ConstantCoefficients(slice, channelIndex);
            nnpcPluginPtr->ApplyPhaseCorrection(rawData, linearCoefficients, constantCoefficients);

            const ComplexFloatMatrix kSpaceToTransform = (rampSamplingPlugin ? rampSamplingPlugin->ApplyRampSampling(rawData) : rawData);

            if(isAssetScan && !homodyneEnabled)
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
            asset->Unalias(unaliasedHighPassData, unaliasedLowPassData, homodyneHighPassFilteredData, homodyneLowPassFilteredData, slice, sliceCornersForAsset);
            accumulatedChannels = unaliasedHighPassData;
            Cartesian2D::Homodyne::ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
        }
        else if(isAssetScan && !homodyneEnabled)
        {   
            asset->Unalias(accumulatedChannels, aliasedData, slice, sliceCornersForAsset);
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
            finalImage(Range::all(), Range(0,kissoffViews-1)) = 0.0f;
            finalImage(Range::all(), Range(finalImage.extent(secondDim)-kissoffViews, finalImage.extent(secondDim)-1)) = 0.0f;
        }

        const SliceCorners sliceCorners = sliceTable.SliceCorners(slice);
        gwPlugin.Run(finalImage, sliceCorners, slice);

        finalImage *= finalImageScaler;

        const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(slice);
        FloatMatrix rotatedImage = RotateTranspose::Apply<float>(finalImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

        Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

        ShortMatrix shortImage(rotatedImage.shape());
        shortImage = cast<short>(rotatedImage);

        std::stringstream dicomFileName;
        dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";
        const ImageCorners imageCorners(sliceCorners, sliceOrientation);
        dicomSeries.SaveImage(dicomFileName.str(), shortImage, slice, imageCorners);

        ++imageNumber;
    } // End Slice Loop
    
    return true;
}
