// Copyright 2018 General Electric Company. All rights reserved.

#include <boost/math/constants/constants.hpp>

#include <MDArray/MDArray.h>
#include <MDArray/Fourier.h>

#include <Hdf5/Snap.h>

#include <Orchestra/Acquisition/ControlPacket.h>
#include <Orchestra/Acquisition/ControlTypes.h>
#include <Orchestra/Acquisition/FrameControl.h>

#include <Orchestra/Acquisition/Core/ArchiveStorage.h>

#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>

#include <Orchestra/Cartesian2D/Homodyne.h>
#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian2D/RampSamplingKernel.h>
#include <Orchestra/Cartesian2D/RampSamplingPlugin.h>

#include <Orchestra/Common/ExamData.h>
#include <Orchestra/Common/ExamDataStorageInfo.h>
#include <Orchestra/Common/ExamStorageDirectory.h>
#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ProgramOptions.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Common/SliceInfoTable.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Epi/LxControlSource.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/RowFlipParameters.h>
#include <Orchestra/Epi/RowFlipPlugin.h>
#include <Orchestra/Epi/SelfNavDynamicPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/VariableRampSamplingKernel.h>
#include <Orchestra/Epi/VariableRampSamplingPlugin.h>

#include <Orchestra/Epi/Diffusion/DynamicPhaseCorrectionManager.h>

#include <Orchestra/Gradwarp/GradwarpFileSource.h>
#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/CommonTypes.h>
#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/LxDownloadData.h>

#include "EpiScanArchiveRecon.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Epi;

void GERecon::Epi::EpiScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive) // parasoft-suppress  METRICS-22 "Procedural rehearsal code"
{
    // Create a trace object to log messages to the console for this reeharsal
    GERecon::Trace trace("EpiScanArchiveRehearsal");
   
    // Initialize a directory to store the dicom image outputs from this rehearsal
    boost::filesystem::path scanArchiveFullPath = scanArchive->Path();
    const boost::filesystem::path dicomDirectory = scanArchiveFullPath.parent_path() / "dicomImages";
    if(!boost::filesystem::exists(dicomDirectory))
    {
        boost::filesystem::create_directory(dicomDirectory);
    }

    // Create an Epi LxControlSource object to interpret recon parameters. The LxControlSource object is used
    // to interpret parameters from the pool header and store the parameters in a processing control object.
    // The processing control object can be used to 
    const GERecon::Legacy::LxDownloadDataPointer downloadData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchive->LoadDownloadData());
    const boost::shared_ptr<GERecon::Epi::LxControlSource> controlSource = boost::make_shared<GERecon::Epi::LxControlSource>(downloadData);
    Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

    // Initialize DicomSeries object to use when creating dicom images
    const Legacy::DicomSeries dicomSeries(downloadData); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Extract parameters from processing control 
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const int extraFramesBottom = processingControl->Value<int>("ExtraFramesBottom");
    const int numAcquisitionsPerRepetition = processingControl->Value<int>("NumAcquisitionsPerRepetition");
    const bool variableKernelRampSamplingInterpolation = processingControl->Value<bool>("FastRampSamplingEnabled");
    const bool rampSamplingInterpolation = processingControl->Value<bool>("RampSamplingEnabled");
    const int phaseCorrectionType = processingControl->ValueStrict<Epi::LxControlSource::PhaseCorrectionType>("PhaseCorrectionType");
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const bool assetEnabled = processingControl->ValueStrict<bool>("Asset");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int kissoffViews = processingControl->Value<int>("TransformKissoffViews") + processingControl->Value<int>("FinalImageKissoffViews");
    const bool assetPhaseCorrectionOptimizationEnabled = processingControl->ValueStrict<bool>("AssetPhaseCorrectionOptimizationEnabled");
    const int totalNumReferenceViews = extraFramesTop + extraFramesBottom;
    const int numPcCoefficients = processingControl->Value<int>("NumPcCoefficients");

    // Prior to initializing external file sources, extract all files saved in the scan archive
    // Because files contained in the scan archive are extracted to the locations of their corresponding
    // GERecon::Path's, set the GERecon::Path locations prior to loading the saved files
    GERecon::Path::InputAppData(scanArchiveFullPath.parent_path());
    GERecon::Path::InputExamData(scanArchiveFullPath.parent_path());
    GERecon::Path::ScannerConfig(scanArchiveFullPath.parent_path());
    scanArchive->LoadSavedFiles();

    // Initialize the row flip plugin that is used to apply the EPI row flip. Whether or not
    // to flip a given readout is controlled by the rowflip.param file which is included as one
    // of the saved files in the scan archive.
    const int totalNumberOfAcquiredViews = totalNumReferenceViews + acquiredYRes;
    const RowFlipParametersPointer rowFlipParams = boost::make_shared<RowFlipParameters>(totalNumberOfAcquiredViews);
    RowFlipPlugin rowFlipPlugin(rowFlipParams, *processingControl); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Initialize phase correction plugin. Either static or dynamic phase correction is supported 
    // by this rehearsal pipeline. If dynamic phase correction is used, the dynamic phase correction
    // manager is used to coordinate passing data for a particular slice from one temporal phase
    // to the next temporal phase
    boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin> staticPhaseCorrectionPlugin;
    boost::shared_ptr<Epi::SelfNavDynamicPhaseCorrectionPlugin> dynamicPhaseCorrectionPlugin;
    boost::shared_ptr<Epi::Diffusion::DynamicPhaseCorrectionManager> dynamicPhaseCorrectionManager;
    boost::shared_ptr<PhaseCorrectionReferenceFile> phaseCorrectionReferenceFile;
    phaseCorrectionReferenceFile = boost::make_shared<GERecon::Epi::PhaseCorrectionReferenceFile>(*processingControl, PhaseCorrectionReferenceFile::ReadMode);
    std::vector<unsigned int> slicePhaseIndexCounters;
    if(phaseCorrectionType == Epi::LxControlSource::SelfNavigatedDynamic)
    {
        dynamicPhaseCorrectionPlugin = boost::make_shared<Epi::SelfNavDynamicPhaseCorrectionPlugin>(*processingControl);
        dynamicPhaseCorrectionManager = Epi::Diffusion::DynamicPhaseCorrectionManager::Instance(*processingControl);

        // Initialize vector used to keep track of which phase index each slice is on
        slicePhaseIndexCounters.resize(sliceTable.GeometricSliceLocations(), 0);
    }
    else if(phaseCorrectionType == Epi::LxControlSource::NearestNeighborStatic)
    {
        staticPhaseCorrectionPlugin = boost::make_shared<Epi::NearestNeighborStaticPhaseCorrectionPlugin>(*processingControl);
    }
    else
    {
        trace.ConsoleMsg("Unsupported phase correction type. Phase correction will be skipped!");
    }

    // Initialize ramp sampling file and plugin. The ramp sampling plugin interpolates data sampled
    // on the gradient ramps to linearly spaced samples in frequency space
    boost::shared_ptr<RampSamplingPlugin> rampSamplingPlugin;
    if(variableKernelRampSamplingInterpolation)
    {
        // If variable kernel ramp sampling interpolation is enabled then a truncated sinc function
        // is used to interpolate to linearly sampled points in frequency space. The truncated sinc
        // function length varies from one sample location to the next
        const boost::shared_ptr<GERecon::Epi::VariableRampSamplingKernel> variableRampSamplingKernel = boost::make_shared<VariableRampSamplingKernel>(*processingControl);
        rampSamplingPlugin = boost::make_shared<VariableRampSamplingPlugin>(variableRampSamplingKernel, *processingControl);
    }
    else if(rampSamplingInterpolation)
    {
        // If the variable kernel option is not enabled then a sinc function which uses all available
        // readout points is used to interpolate to linearly sampled frequency space
        const boost::shared_ptr<GERecon::RampSamplingKernel> rampSamplingKernel = boost::make_shared<RampSamplingKernel>(*processingControl);
        rampSamplingPlugin = boost::make_shared<RampSamplingPlugin>(rampSamplingKernel, *processingControl);
    }
    else
    {
        trace.ConsoleMsg("Ramp sampling interpolation is disabled for this recon");
    }
    
    // Setup a kSpace transformer object. This object handles kSpace apodization as well as partial fourier (if enabled)
    Cartesian2D::KSpaceTransformer transformer(*processingControl); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Initialize ASSET, if enabled 
    ComplexFloatCube homodyneHighPassFilteredDataForAsset;
    ComplexFloatCube homodyneLowPassFilteredDataForAsset;
    ComplexFloatMatrix unaliasedHighPassData;
    ComplexFloatMatrix unaliasedLowPassData;
    ComplexFloatCube aliasedImagesForAsset;
    boost::shared_ptr<Asset::Worker> assetWorker;
    if(assetEnabled)
    {
        // Initialize ASSET calibration
        const boost::shared_ptr<Asset::CalibrationFile> assetCalibration = boost::make_shared<Asset::CalibrationFile>(*processingControl);
        assetCalibration->TopDirectory(GERecon::Path::InputExamData());
        GERecon::ExamStorageDirectory examStorageDirectory;
        std::vector<boost::filesystem::path> assetCalibrationFiles = examStorageDirectory.FindPotentialFiles(assetCalibration->StorageInfo());

        if(assetCalibrationFiles.size() > 0)
        {
            const boost::filesystem::path assetCalibrationFilePath = assetCalibrationFiles[0];
            assetCalibration->Load(assetCalibrationFilePath);
            assetWorker = boost::make_shared<Asset::Worker>(assetCalibration, *processingControl);
        }
        else
        {
            trace.ConsoleMsg("Could not find AssetCalibration file!");
            return;
        }

        // Initialize Workspace
        if(homodyneEnabled)
        {
            // Initialize ASSET + Homodyne workspace
            homodyneHighPassFilteredDataForAsset.resize(reconXRes, reconYRes, numChannels);
            homodyneLowPassFilteredDataForAsset.resize(reconXRes, reconYRes, numChannels);
            unaliasedHighPassData.resize(reconXRes, reconYRes);
            unaliasedLowPassData.resize(reconXRes, reconYRes);
        }
        else
        {
            // Initialize ASSET only workspace
            aliasedImagesForAsset.resize(reconXRes, reconYRes, numChannels);
        }
    }

    // Initialize SumOfSquares channel combiner to use if SumOfSquares channel combination is enabled
    SumOfSquares channelCombiner(processingControl->ValueStrict<FloatVector>("ChannelWeights")); 

    // Initialize Gradwarp. The gradwarp file source will use the gw_coils.dat file stored in the ScanArchive
    // to initialize the gradwarp plugin
    const GradwarpFileSourcePtr gradwarpFileSource = boost::make_shared<GradwarpFileSource>();
    GradwarpPlugin gradwarpPlugin(*processingControl, GERecon::TwoDGradwarp, gradwarpFileSource);
    
    // Create MDArray workspaces
    ComplexFloatMatrix accumulatedChannels(reconXRes, reconYRes);
    ComplexFloatMatrix transformedData(reconXRes, reconYRes);

    // Match product scaling - for partial ky scans apply the following scalar
    const float finalImageScaler = homodyneEnabled ? 256.0f / (reconXRes * reconXRes) : 1.0f;

    // Load the archive storage which contains all acquisition data held in the archive
    GERecon::Acquisition::ArchiveStoragePointer archiveStorage = GERecon::Acquisition::ArchiveStorage::Create(scanArchive);
    
    // Determine how many control (DAB) packets are contained in the storage
    const size_t numControls = archiveStorage->AvailableControlCount();

    // Initialize an acquisition (pass) counter. This counter is incremented every time an acquisition (pass) done
    // packet is received from the PSD. Also initialize a phase index counter to keep track of which temporal
    // phase of the scan we're on. Each time the acquisitionPassCounter reaches the numAcquisitionsPerRepetition
    // value, the phaseIndex is incremented and the acquisitionPassCounter is reset to 0
    int acquisitionPassCounter = 0;
    int phaseIndex = 0;

    // Create an HDF5 file to save debug data to. This debug file will be enabled if the Snap Key
    // specified in this constructor is included on the command line for this rehearsal application:
    //      --GEHdf5.Snap.Key EpiScanArchiveReherasalDebug
    const GEHdf5::Snap debugFile("EpiScanArchiveRehearsal.h5", "EpiScanArchiveReherasalDebug");

    // Loop over all control packets in the archive. Some control packets are scan control packets
    // which may indicate the end of an acquisition (pass) or the end of the scan. Other control
    // packets are frame control packets which describe the raw frame (or view) data they're 
    // associated with. All control packets and associated frame data are stored in the archive
    // in the order they're acquired.
    for(size_t controlPacketIndex = 0; controlPacketIndex < numControls; ++controlPacketIndex)
    {
        const GERecon::Acquisition::FrameControlPointer controlPacketAndFrameData = archiveStorage->NextFrameControl();
        
        // Variables extracted from acquisition packets
        int sliceIndexInAcq = 0;
        short viewIncrement = 0;
        bool reconstructDataAssociatedWithThisPacket = false;

        if(controlPacketAndFrameData->Control().Opcode() == Acquisition::ScanControlOpcode)
        {
            const GERecon::Acquisition::ScanControlPacket scanControlPacket = controlPacketAndFrameData->Control().Packet().As<GERecon::Acquisition::ScanControlPacket>();

            // For diffusion scans, the same volume may be acquired multiple times for the various diffusion passes (T2, B-Value/Diffusion Direction). 
            // Thus, when this counter reaches the number of acquisitions per repetition, reset the counter back to 0.
            ++acquisitionPassCounter;
            if(acquisitionPassCounter == numAcquisitionsPerRepetition)
            {
                acquisitionPassCounter = 0;
                ++phaseIndex;
            }

            trace.ConsoleMsg("Received scan control packet, scanControlBitMask: %d, acquisitionCounter is now: %d", scanControlPacket.scanControl, acquisitionPassCounter);
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::EpiOpcode)
        {
            const Acquisition::HyperFrameControlPacket hyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::HyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(hyperFramePacket.sliceNumH, hyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(hyperFramePacket.viewSkipH, hyperFramePacket.viewSkipL));
            reconstructDataAssociatedWithThisPacket = true;
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::DiffusionHyperScanOpcode)
        {
            const Acquisition::DiffusionHyperFrameControlPacket diffusionHyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::DiffusionHyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(diffusionHyperFramePacket.sliceNumH, diffusionHyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(diffusionHyperFramePacket.viewSkipH, diffusionHyperFramePacket.viewSkipL));
            reconstructDataAssociatedWithThisPacket = true;
        }

        // Use the slice table to determine which geometric slice index this data corresponds to    
        const int geometricSliceIndex = sliceTable.GeometricSliceNumber(acquisitionPassCounter, sliceIndexInAcq);

        if(reconstructDataAssociatedWithThisPacket)
        {
            trace.ConsoleMsg("Reconstructing geometric slice %d of phase %d", geometricSliceIndex, phaseIndex);
            std::stringstream phaseNameString;
            std::stringstream sliceNameString;
            phaseNameString << "phase_" << std::setw(3) << std::setfill('0') << phaseIndex;
            sliceNameString << "slice_" << std::setw(3) << std::setfill('0') << geometricSliceIndex;

            // Extract acquisition data associated with this hyper frame packet.
            // Note that this rehearsal only supports single shot EPI data. Each hyper frame control packet
            // is associated with a single shot of data. Since this rehearsal only supports single shot,
            // each control packet is therefore associated with a single shot, single slice of data. Multi-shot
            // support may be added in the future. For example, for a 2 shot scan, each slice of data needs
            // is composed of 2 EPI shots. Thus, 2 control packets and their associated data must be used to
            // obtain the kSpace for a single slice
            ComplexFloatCube allRawData = controlPacketAndFrameData->Data();
            debugFile.Write(allRawData, phaseNameString.str() + "/" + sliceNameString.str() + "/allRawData");
                
            if(viewIncrement < 0)
            {
                // If the view increment in this packet is less than zero, that indicates that the EPI echo train was
                // acquired in a bottom-up fashion. Thus, the first view in the echo train is the "bottom" view in
                // the kSpace matrix. To orient this data, flip the data along the second dimension (phase encode dimension
                // such that phase encode index 0 corresponds to the "top" kSpace view and the last phase encode index
                // corresponds to the "bottom" kSpace view.
                // Note that this rehearsal only supports single shot acquisitions. For single shot acquisitions, all slice
                // data is accompanied by a single control packet; thus, additional sorting/bookkeeping of acquisition data
                // is not required.
                for(int channelIndex = 0; channelIndex < allRawData.extent(secondDim); ++channelIndex)
                {
                    MDArray::ComplexFloatMatrix currentChannel = allRawData(Range::all(), channelIndex, Range::all());
                    MDArray::Flip(currentChannel, secondDim);
                }
            }

            // If additional non-phase-encoded reference views are acquired, they are either at the top or bottom of the acquired
            // data. Consider these reference views and create a reference to them here, if they are acquired
            ComplexFloatCube rawImageData(allRawData.extent(firstDim), acquiredYRes, numChannels);
            const Range imageViewRange(extraFramesTop, extraFramesTop + acquiredYRes - 1);
            ComplexFloatCube rawReferenceData;
            Range referenceViewRange;
            if(totalNumReferenceViews > 0)
            {
                // If reference views are acquired, allocate space for them
                rawReferenceData.resize(allRawData.extent(firstDim), totalNumReferenceViews, numChannels);

                // Also create a range object representing the view-indices containing the reference views
                referenceViewRange = extraFramesTop > 0 ? Range(fromStart, extraFramesTop - 1) : (extraFramesBottom > 0) ? Range(acquiredYRes, extraFramesBottom + acquiredYRes - 1) : Range();
            }

            for(int channel = 0; channel < numChannels; ++channel)
            {
                rawImageData(Range::all(), Range::all(), channel) = allRawData(Range::all(), channel, imageViewRange);
                if(rawReferenceData.size() > 0)
                {
                    rawReferenceData(Range::all(), Range::all(), channel) = allRawData(Range::all(), channel, referenceViewRange);
                }
            }

            debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/rawImageData");
            debugFile.Write(rawReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/rawReferenceData");
            
            // RowFlip for EPI raw data. Flip both the raw image data and raw dynamic reference data
            for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
            {
                ComplexFloatMatrix currentRawImageData = rawImageData(Range::all(), Range::all(), channelIndex);
                rowFlipPlugin.ApplyImageDataRowFlip(currentRawImageData);

                if(rawReferenceData.size() > 0)
                {
                    ComplexFloatMatrix currentDynamicReferenceData = rawReferenceData(Range::all(), Range::all(), channelIndex);
                    rowFlipPlugin.ApplyReferenceDataRowFlip(currentDynamicReferenceData);
                }
            }
            debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/rowFlippedRawImageData");
            debugFile.Write(rawReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/rowFlippedRawReferenceData");

            if(staticPhaseCorrectionPlugin)
            {
                // Apply static phase correction
                for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
                {
                    const int channelToUseForPhaseCorrection = assetPhaseCorrectionOptimizationEnabled ? phaseCorrectionReferenceFile->MaximumSnrChannel(geometricSliceIndex) : channelIndex;
                    const FloatVector constantCoefficientsCurrentChannel = phaseCorrectionReferenceFile->ConstantCoefficients(geometricSliceIndex, channelToUseForPhaseCorrection);
                    const FloatVector linearCoefficientsCurrentChannel = phaseCorrectionReferenceFile->LinearCoefficients(geometricSliceIndex, channelToUseForPhaseCorrection);
                    ComplexFloatMatrix dataToCorrect = rawImageData(Range::all(), Range::all(), channelIndex);
                    staticPhaseCorrectionPlugin->ApplyPhaseCorrection(dataToCorrect, linearCoefficientsCurrentChannel, constantCoefficientsCurrentChannel);
                }
            }
            else if(dynamicPhaseCorrectionPlugin)
            {
                // Retrieve the static reference coefficients to use
                FloatMatrix linearCoefficients(rawImageData.extent(secondDim), numChannels);
                FloatMatrix constantCoefficients(rawImageData.extent(secondDim), numChannels);
                for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
                {
                    const int channelToUseForPhaseCorrection = assetPhaseCorrectionOptimizationEnabled ? phaseCorrectionReferenceFile->MaximumSnrChannel(geometricSliceIndex) : channelIndex;
                    constantCoefficients(Range::all(), channelIndex) = phaseCorrectionReferenceFile->ConstantCoefficients(geometricSliceIndex, channelToUseForPhaseCorrection);
                    linearCoefficients(Range::all(), channelIndex) = phaseCorrectionReferenceFile->LinearCoefficients(geometricSliceIndex, channelToUseForPhaseCorrection);
                }
                FloatCube pcCoefficients(acquiredYRes, numChannels, numPcCoefficients);
                pcCoefficients = phaseCorrectionReferenceFile->PcCoefficients(geometricSliceIndex);

                // Chop to center transforms
                rawImageData(Range(fromStart,toEnd,2), Range::all(), Range::all()) *= -1.0f;
                rawReferenceData(Range(fromStart,toEnd,2), Range::all(), Range::all()) *= -1.0f;

                Fourier::Ifft(rawImageData, MDArray::firstDim);
                Fourier::Ifft(rawReferenceData, MDArray::firstDim);

                rawImageData /= static_cast<float>(rawImageData.extent(MDArray::firstDim));
                rawReferenceData /= static_cast<float>(rawReferenceData.extent(MDArray::firstDim));

                ComplexFloatCube baselineImage;
                ComplexFloatCube baselineDynamicReferenceData;
                FloatMatrix constantPhasePerDynamicReferenceView;
                std::vector<MDArray::TinyVector<float, 4> > dynamicCoefficients(numChannels, MDArray::TinyVector<float, 4>(0.0f, 0.0f, 0.0f, 0.0f));
                const unsigned int currentSlicePhaseIndex = slicePhaseIndexCounters[geometricSliceIndex];
                if(currentSlicePhaseIndex == 0)
                {
                    // Copy baseline data to a new MDArray and commit that MDArray to be stored
                    // for use in all remaining temporal phases of this reconstruction
                    baselineImage.resize(rawImageData.shape());
                    baselineDynamicReferenceData.resize(rawReferenceData.shape());
                    baselineImage = rawImageData;
                    baselineDynamicReferenceData = rawReferenceData;
                    dynamicPhaseCorrectionManager->CommitBaselineData(geometricSliceIndex, baselineImage, baselineDynamicReferenceData);
                    
                    // Initialize dynamic phase correction input data from previous to 0.0 for the first (baseline) phase
                    // After the baseline phase, this data will come from the previous temporal phase for each slice
                    constantPhasePerDynamicReferenceView.resize(rawReferenceData.extent(secondDim), numChannels);
                    constantPhasePerDynamicReferenceView = 0.0f;
                    // The dynamicCoefficients from the previous phase were initialized to 0.0 above
                }
                else
                {
                    // Retrieve baseline data for the current slice
                    baselineImage.reference(dynamicPhaseCorrectionManager->BaselineImage(geometricSliceIndex));
                    baselineDynamicReferenceData.reference(dynamicPhaseCorrectionManager->BaselineDynamicReferenceViews(geometricSliceIndex));

                    // Retrieve data from previous temporal phase for the current geometric slice. This data
                    // is input data to the dynamic phase correction algorithm
                    const int previousPhaseIndex = currentSlicePhaseIndex - 1;
                    constantPhasePerDynamicReferenceView.reference(dynamicPhaseCorrectionManager->ConstantPhasePerDynamicReferenceView(previousPhaseIndex, geometricSliceIndex));
                    dynamicCoefficients = dynamicPhaseCorrectionManager->DynamicCoefficients(previousPhaseIndex, geometricSliceIndex);
                }

                const int shotIndex = 0; //Only support single shot
                dynamicPhaseCorrectionPlugin->ApplyPhaseCorrection(constantPhasePerDynamicReferenceView, dynamicCoefficients, rawImageData, rawReferenceData, baselineImage, baselineDynamicReferenceData, pcCoefficients, shotIndex);

                // Commit data from the current temporal phase for use in the next temporal phase
                dynamicPhaseCorrectionManager->CommitDynamicData(currentSlicePhaseIndex, geometricSliceIndex, dynamicCoefficients, constantPhasePerDynamicReferenceView);

                // Increment slice phase index counter
                slicePhaseIndexCounters[geometricSliceIndex] = slicePhaseIndexCounters[geometricSliceIndex] + 1;

                // Transform phase corrected data back to kx,ky for apodization filtering and 2D kSpace transform
                Fourier::Fft(rawImageData, MDArray::firstDim);
                rawImageData(Range(fromStart, toEnd, 2), Range::all(), Range::all()) *= -1.0f;
            }
            debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/phaseCorrectedRawImageData");

            // VRGF
            ComplexFloatCube kSpaceToTransform = rawImageData;
            if(rampSamplingPlugin)
            {
                kSpaceToTransform.resize(rampSamplingPlugin->InterpolatedXRes(), rawImageData.extent(secondDim), rawImageData.extent(thirdDim));
                for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
                {
                    kSpaceToTransform(Range::all(), Range::all(), channelIndex) = rampSamplingPlugin->ApplyRampSampling(rawImageData(Range::all(), Range::all(), channelIndex));
                }
            }
            debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/phaseCorrectedRawImageData");

            // Loop over all channels and transform the data
            for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
            {
                const ComplexFloatMatrix currentChannelDataToTransform = kSpaceToTransform(Range::all(), Range::all(), channelIndex);

                if(assetEnabled && !homodyneEnabled) // parasoft-suppress  METRICS-23 "rehearsal code handles many cases in a single int main pipeline"
                {
                    ComplexFloatMatrix aliasedChannelData = aliasedImagesForAsset(Range::all(), Range::all(), channelIndex);
                    transformer.Apply(aliasedChannelData, currentChannelDataToTransform);
                    debugFile.Write(aliasedChannelData, phaseNameString.str() + "/" + sliceNameString.str() + "/aliasedChannelData");
                }
                else if(homodyneEnabled && assetEnabled) // parasoft-suppress  METRICS-23 "Procedural rehearsal code handles many cases in a single int main pipeline"
                {
                    ComplexFloatMatrix highPassImage = homodyneHighPassFilteredDataForAsset(Range::all(), Range::all(), channelIndex);
                    ComplexFloatMatrix lowPassImage = homodyneLowPassFilteredDataForAsset(Range::all(), Range::all(), channelIndex);
                    transformer.Apply(highPassImage, lowPassImage, currentChannelDataToTransform);
                }
                else if(!assetEnabled) // parasoft-suppress  METRICS-23 "Procedural, for-loop rehearsal code"
                {
                    // Handles the homodyne only, zerofilling, or full kx,ky cases
                    transformer.Apply(transformedData, currentChannelDataToTransform);
                    channelCombiner.Accumulate(transformedData, channelIndex);
                }
            }

            // Compute the channel combined image
            const GERecon::SliceCorners sliceCornersForAsset = sliceTable.AcquiredSliceCorners(geometricSliceIndex);
            if(assetEnabled && homodyneEnabled)
            {
                debugFile.Write(homodyneHighPassFilteredDataForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/homodyneHighPassFilteredDataForAsset");
                debugFile.Write(homodyneLowPassFilteredDataForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/homodyneLowPassFilteredDataForAsset");
                assetWorker->Unalias(unaliasedHighPassData, unaliasedLowPassData, homodyneHighPassFilteredDataForAsset, homodyneLowPassFilteredDataForAsset, geometricSliceIndex, sliceCornersForAsset);
                accumulatedChannels = unaliasedHighPassData;
                Cartesian2D::Homodyne::ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
            }
            else if(assetEnabled && !homodyneEnabled)
            {   
                debugFile.Write(aliasedImagesForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/aliasedImagesForAsset");
                assetWorker->Unalias(accumulatedChannels, aliasedImagesForAsset, geometricSliceIndex, sliceCornersForAsset);
            }
            else if(!assetEnabled)
            {
                accumulatedChannels = channelCombiner.GetCombinedImage();
                channelCombiner.Reset();
            }
            debugFile.Write(accumulatedChannels, phaseNameString.str() + "/" + sliceNameString.str() + "/accumulatedChannels");

            FinalizeImage(accumulatedChannels, kissoffViews, sliceTable, gradwarpPlugin, dicomSeries, finalImageScaler, phaseIndex, geometricSliceIndex, dicomDirectory);
        }
    }
}

void GERecon::Epi::FinalizeImage(const MDArray::ComplexFloatMatrix& imageToFinalize, 
                                 const int kissoffViews, 
                                 const GERecon::SliceInfoTable& sliceTable, 
                                 GERecon::GradwarpPlugin& gradwarpPlugin,
                                 const GERecon::Legacy::DicomSeries& dicomSeries,
                                 const float finalImageScaler,
                                 const int phaseIndex,
                                 const int geometricSliceIndex,
                                 const boost::filesystem::path& dicomDirectory)
{
    // Create a magnitude image and zero out any kissoff views
    FloatMatrix magnitudeImage(imageToFinalize.shape());
    magnitudeImage = MDArray::abs(imageToFinalize);
    if(kissoffViews > 0)
    {
        magnitudeImage(Range::all(), Range(fromStart,kissoffViews-1)) = 0.0f;
        magnitudeImage(Range::all(), Range(magnitudeImage.extent(secondDim)-kissoffViews, magnitudeImage.extent(secondDim)-1)) = 0.0f;
    }

    const SliceCorners sliceCorners = sliceTable.AcquiredSliceCorners(geometricSliceIndex);
    gradwarpPlugin.Run(magnitudeImage, sliceCorners, geometricSliceIndex);
    magnitudeImage *= finalImageScaler;

    const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(geometricSliceIndex);
    FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

    Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

    ShortMatrix shortImage(rotatedImage.shape());
    shortImage = cast<short>(rotatedImage);

    const int imageNumber = sliceTable.GeometricSliceLocations() * phaseIndex + geometricSliceIndex;
    std::stringstream dicomFileName;
    dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(5) << std::setfill('0') << imageNumber << ".dcm";
    const ImageCorners imageCorners(sliceCorners, sliceOrientation);
    dicomSeries.SaveImage(dicomFileName.str(), shortImage, imageNumber, imageCorners);
}
