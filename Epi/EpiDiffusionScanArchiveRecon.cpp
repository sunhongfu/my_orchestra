// Copyright 2018 General Electric Company. All rights reserved.

#include <MDArray/MDArray.h>
#include <MDArray/Fourier.h>

#include <Orchestra/Acquisition/ControlPacket.h>
#include <Orchestra/Acquisition/ControlTypes.h>
#include <Orchestra/Acquisition/FrameControl.h>

#include <Orchestra/Acquisition/Core/ArchiveStorage.h>

#include <Orchestra/Arc/Arc.h>
#include <Orchestra/Arc/CalibrationManager.h>

#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>

#include <Orchestra/Calibration/Common/RawFile.h>
#include <Orchestra/Calibration/Common/RawFileInfo.h>

#include <Orchestra/Cartesian2D/Homodyne.h>
#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian2D/RampSamplingKernel.h>
#include <Orchestra/Cartesian2D/RampSamplingPlugin.h>

#include <Orchestra/Common/ExamData.h>
#include <Orchestra/Common/ExamDataStorageInfo.h>
#include <Orchestra/Common/ExamStorageDirectory.h>
#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Common/SliceInfoTable.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Epi/MultibandWorker.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/RowFlipParameters.h>
#include <Orchestra/Epi/RowFlipPlugin.h>
#include <Orchestra/Epi/SelfNavDynamicPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/VariableRampSamplingKernel.h>
#include <Orchestra/Epi/VariableRampSamplingPlugin.h>

#include <Orchestra/Epi/Diffusion/DynamicPhaseCorrectionManager.h>
#include <Orchestra/Epi/Diffusion/HoecReconCorrection.h>
#include <Orchestra/Epi/Diffusion/TensorVectors.h>

#include <Orchestra/Epi/DistortionCorrection/B0Correction.h>
#include <Orchestra/Epi/DistortionCorrection/RpgParameters.h>

#include <Orchestra/Gradwarp/GradwarpFileSource.h>
#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/CommonTypes.h>
#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/LxDownloadData.h>

#include <Orchestra/MoCo/MoCo.h>
#include <Orchestra/MoCo/MoCoParameters.h>
#include <Orchestra/MoCo/MoCoParametersTuner.h>

#include "EpiDiffusionRecon.h"
#include "EpiDiffusionScanArchiveRecon.h"
#include "EpiReferenceScanRecon.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Epi;

void GERecon::Epi::EpiDiffusionScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive) // parasoft-suppress  METRICS-22 "Procedural rehearsal code"
{
    // Create an Epi LxControlSource object to interpret recon parameters. The LxControlSource object is used
    // to interpret parameters from the pool header and store the parameters in a processing control object.
    // The processing control object can be used to 
    const GERecon::Legacy::LxDownloadDataPointer downloadData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchive->LoadDownloadData());
    const boost::shared_ptr<GERecon::Epi::LxControlSource> controlSource = boost::make_shared<GERecon::Epi::LxControlSource>(downloadData);
    const Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

    // Extract parameters from processing control 
    const int numT2Images = processingControl->Value<int>("NumT2ImagesPerLocation");

    // Prior to initializing external file sources, extract all files saved in the scan archive
    // Because files contained in the scan archive are extracted to the locations of their corresponding
    // GERecon::Path's, set the GERecon::Path locations prior to loading the saved files
    boost::filesystem::path scanArchiveFullPath = scanArchive->Path();
    GERecon::Path::InputAppData(scanArchiveFullPath.parent_path());
    GERecon::Path::InputExamData(scanArchiveFullPath.parent_path());
    GERecon::Path::ScannerConfig(scanArchiveFullPath.parent_path());
    scanArchive->LoadSavedFiles();

    // If it is an integrated ref scan, reconstruct the reference scan
    if(processingControl->Value<bool>("IntegratedReferenceScan") && !boost::filesystem::exists(scanArchiveFullPath.parent_path() / "ref.h5"))
    {
        EpiReferenceScanRecon(scanArchive);
    }

    // Get the t2 volume and b-value volume from the DAB packets.
    // Each of this volume will be one acquisition worth of data necessary for reconstruction
    MDArray::Array<float,5> t2Volume;
    MDArray::Array<float,6> bValueVolume;
    StoreVolumes(t2Volume, bValueVolume, scanArchive, processingControl);

    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const int numSlices = sliceTable.GeometricSliceLocations();

    // Initialize DicomSeries object to use when creating dicom images
    const Legacy::DicomSeries dicomSeries(downloadData); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Initialize a directory to store the dicom image outputs from this rehearsal
    const boost::filesystem::path dicomDirectory = scanArchiveFullPath.parent_path() / "dicomImages";
    if(!boost::filesystem::exists(dicomDirectory))
    {
        boost::filesystem::create_directory(dicomDirectory);
    }

    // These matrices are generated in T2SegmentRecon and used in DiffusionSegmentRecon
    MDArray::FloatCube displacementMap(reconXRes,reconYRes, numSlices);
    displacementMap = 0.0f;
    MDArray::FloatCube referenceImageForMoCo(reconXRes,reconYRes, numSlices);
    displacementMap = 0.0f;

    if(numT2Images > 0)
    {
        T2SegmentRecon(displacementMap, referenceImageForMoCo, t2Volume, processingControl, dicomDirectory, dicomSeries);
    }
    DiffusionSegmentRecon(bValueVolume, displacementMap, referenceImageForMoCo, processingControl, dicomDirectory, dicomSeries);
}

void Epi::StoreVolumes(MDArray::Array<float,5>& t2Volume, // parasoft-suppress  METRICS-22 "Required for procedual code"
                       MDArray::Array<float,6>& bValueVolume, 
                       const GERecon::ScanArchivePointer& scanArchive,
                       const Control::ProcessingControlPointer& processingControl)
{
    // Create a trace object to log messages to the console for this rehearsal
    GERecon::Trace trace("Decode control packets");

    const int numAcquisitionsPerRepetition = processingControl->Value<int>("NumAcquisitionsPerRepetition");
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const int numSlices = sliceTable.GeometricSliceLocations();

    // Resize all matrices and fill data
    const int acquiredXRes = processingControl->Value<int>("AcquiredXRes");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const int extraFramesBottom = processingControl->Value<int>("ExtraFramesBottom");
    const int totalNumReferenceViews = extraFramesTop + extraFramesBottom;
    const int totalNumberOfAcquiredViews = totalNumReferenceViews + acquiredYRes;
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int numT2Images = processingControl->Value<int>("NumT2ImagesPerLocation");
    const int t2Nex = processingControl->Value<int>("DiffusionT2NumberOfNexes");
    const int totalNumT2 = numT2Images * t2Nex * numSlices;
    const int totalNumRef = processingControl->Value<bool>("IntegratedReferenceScan") ? numSlices : 0;
    if(numT2Images > 0)
    {
        t2Volume.resize(reconXRes,reconYRes,t2Nex,numT2Images,numSlices);
        t2Volume = 0.0f;
    }

    const int numBValues = processingControl->Value<int>("NumBValues");
    const int numDiffusionDirections = processingControl->Value<int>("NumDiffusionDirections");
    const std::vector<int> diffusionNexTable = processingControl->Value<std::vector<int> >("DiffusionNexTable");
    int maxDiffusionNex = 0;
    for(size_t nex = 0; nex < diffusionNexTable.size(); ++nex)
    {
        if(maxDiffusionNex < diffusionNexTable[nex])
        {
            maxDiffusionNex = diffusionNexTable[nex];
        }
    }
    
    bValueVolume.resize(reconXRes,reconYRes,maxDiffusionNex,numDiffusionDirections,numBValues,numSlices);
    bValueVolume = 0.0f;

    // Load the archive storage which contains all acquisition data held in the archive
    GERecon::Acquisition::ArchiveStoragePointer archiveStorage = GERecon::Acquisition::ArchiveStorage::Create(scanArchive);

    // Determine how many control (DAB) packets are contained in the storage
    const size_t numControls = archiveStorage->AvailableControlCount();

    // Initialize an acquisition (pass) counter. This counter is incremented every time an acquisition (pass) done
    // packet is received from the PSD. Each time the acquisitionPassCounter reaches the numAcquisitionsPerRepetition
    // value, the acquisitionPassCounter is reset to 0
    int acquisitionPassCounter = 0;
    
    // Vector to track number of phases per slice
    std::vector<unsigned int> numPhasesPerSlice;
    numPhasesPerSlice.resize(sliceTable.GeometricSliceLocations(),0);

    ComplexFloatCube allRawData(acquiredXRes,totalNumberOfAcquiredViews,numChannels);
    FloatCube channelCombinedData;
    std::vector<int> unaliasedGeometricSliceNumbers;
    const bool isRpgVolumeAcquired = processingControl->Value<bool>("IsRpgVolumeAcquired");

    // Class intializations for Epi k-space to image space
    const RowFlipParametersPointer rowFlipParams = boost::make_shared<RowFlipParameters>(totalNumberOfAcquiredViews);
    RowFlipPlugin rowFlipPlugin(rowFlipParams, *processingControl); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Phase correction in EPI Diffusion pipeline
    boost::shared_ptr<PhaseCorrectionReferenceFile> phaseCorrectionReferenceFile;
    phaseCorrectionReferenceFile = boost::make_shared<GERecon::Epi::PhaseCorrectionReferenceFile>(*processingControl, PhaseCorrectionReferenceFile::ReadMode);

    // Initialize phase correction plugin. Either static or dynamic phase correction is supported 
    // by this rehearsal pipeline. If dynamic phase correction is used, the dynamic phase correction
    // manager is used to coordinate passing data for a particular slice from one temporal phase
    // to the next temporal phase
    boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin> staticPhaseCorrectionPlugin;
    boost::shared_ptr<Epi::SelfNavDynamicPhaseCorrectionPlugin> dynamicPhaseCorrectionPlugin;
    boost::shared_ptr<Epi::Diffusion::DynamicPhaseCorrectionManager> dynamicPhaseCorrectionManager;

    std::vector<unsigned int> slicePhaseIndexCounters;
    const int phaseCorrectionType = processingControl->ValueStrict<Epi::LxControlSource::PhaseCorrectionType>("PhaseCorrectionType");
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
    const bool variableKernelRampSamplingInterpolation = processingControl->Value<bool>("FastRampSamplingEnabled");
    const bool rampSamplingInterpolation = processingControl->Value<bool>("RampSamplingEnabled");
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

    // Create an HDF5 file to save debug data to. This debug file will be enabled if the Snap Key
    // specified in this constructor is included on the command line for this rehearsal application:
    //      --GEHdf5.Snap.Key EpiDiffusionScanArchiveRehearsalDebug
    const GEHdf5::Snap debugFile("EpiDiffusionScanArchiveRehearsal.h5", "EpiDiffusionScanArchiveRehearsalDebug");
    int packetCount = 0;
    for(size_t controlPacketIndex = 0; controlPacketIndex < numControls; ++controlPacketIndex)
    {
        const GERecon::Acquisition::FrameControlPointer controlPacketAndFrameData = archiveStorage->NextFrameControl();

        // Variables extracted from acquisition packets
        int sliceIndexInAcq = 0;
        int echoIndex = 0;
        short viewIncrement = 0;
        bool isT2Frame = false;
        bool isDiffusionFrame = false;
        int t2Index = 0;
        int bValueIndex = 0;
        int diffDirIndex = 0;
        int imageNexIndex = 0;
        int geometricSliceIndex = 0;
        if(controlPacketAndFrameData->Control().Opcode() == Acquisition::ScanControlOpcode)
        {
            const GERecon::Acquisition::ScanControlPacket scanControlPacket = controlPacketAndFrameData->Control().Packet().As<GERecon::Acquisition::ScanControlPacket>();

            // For diffusion scans, the same volume may be acquired multiple times for the various diffusion passes (T2, B-Value/Diffusion Direction). 
            // Thus, when this counter reaches the number of acquisitions per repetition, reset the counter back to 0.
            ++acquisitionPassCounter;
            if(acquisitionPassCounter == numAcquisitionsPerRepetition)
            {
                acquisitionPassCounter = 0;
            }

            trace.ConsoleMsg("Received scan control packet, scanControlBitMask: %d, acquisitionCounter is now: %d", scanControlPacket.scanControl, acquisitionPassCounter);
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::EpiOpcode)
        {
            const Acquisition::HyperFrameControlPacket hyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::HyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(hyperFramePacket.sliceNumH, hyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(hyperFramePacket.viewSkipH, hyperFramePacket.viewSkipL));
            echoIndex = hyperFramePacket.echoNum;
            imageNexIndex = echoIndex;
            packetCount++;
            int frameType = 0;
            GetPacketInfo(frameType, t2Index, bValueIndex, diffDirIndex, packetCount, numSlices, totalNumRef, totalNumT2, numBValues, numDiffusionDirections, diffusionNexTable);
            isT2Frame = (1 == frameType);
            isDiffusionFrame = (2 == frameType);
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::DiffusionHyperScanOpcode)
        {
            const Acquisition::DiffusionHyperFrameControlPacket diffusionHyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::DiffusionHyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(diffusionHyperFramePacket.sliceNumH, diffusionHyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(diffusionHyperFramePacket.viewSkipH, diffusionHyperFramePacket.viewSkipL));
            echoIndex = diffusionHyperFramePacket.echoNum;
            isT2Frame = (1 == diffusionHyperFramePacket.frameType);
            isDiffusionFrame = (2 == diffusionHyperFramePacket.frameType);
            t2Index = diffusionHyperFramePacket.instanceIndex;
            bValueIndex = diffusionHyperFramePacket.bValueIndex;
            diffDirIndex = Acquisition::GetPacketValue(diffusionHyperFramePacket.diffDirIndexH, diffusionHyperFramePacket.diffDirIndexL);
            imageNexIndex = echoIndex;
        }
        // Use the slice table to determine which geometric slice index this data corresponds to
        if(isT2Frame || isDiffusionFrame)
        {
            geometricSliceIndex = sliceTable.GeometricSliceNumber(acquisitionPassCounter, sliceIndexInAcq);

            // The API returns data with acqXRes x numChannels x acqYRes as dimensions. So, we need to tranpose the second and third dimensions to get acqXRes x acqYRes x numChannels
            const ComplexFloatCube rawDataFromApi = controlPacketAndFrameData->Data();
            allRawData = rawDataFromApi.transpose(firstDim,thirdDim,secondDim);
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
                for(int channelIndex = 0; channelIndex < allRawData.extent(thirdDim); ++channelIndex)
                {
                    MDArray::ComplexFloatMatrix currentChannel = allRawData(Range::all(), Range::all(), channelIndex);
                    MDArray::Flip(currentChannel, secondDim);
                }
            }
            const bool isRpgVolume = (isT2Frame && isRpgVolumeAcquired && t2Index == 0);
            DiffusionCommonRecon(channelCombinedData, unaliasedGeometricSliceNumbers, slicePhaseIndexCounters, allRawData, rowFlipPlugin, phaseCorrectionReferenceFile, staticPhaseCorrectionPlugin,
                dynamicPhaseCorrectionPlugin, dynamicPhaseCorrectionManager, rampSamplingPlugin, processingControl, debugFile, imageNexIndex, geometricSliceIndex, numPhasesPerSlice[geometricSliceIndex], isRpgVolume);
            numPhasesPerSlice[geometricSliceIndex] = numPhasesPerSlice[geometricSliceIndex] + 1;
            // Multiband unaliasing can potentially more than 1 slice worth of data. Store it in either the T2 volume or the B-value volume.
            for(int unaliasedSliceIndex = 0; unaliasedSliceIndex < boost::numeric_cast<int>(unaliasedGeometricSliceNumbers.size()); ++unaliasedSliceIndex)
            {
                const int geometricSliceIndexFromHyperBand = unaliasedGeometricSliceNumbers[unaliasedSliceIndex];
                if(geometricSliceIndexFromHyperBand < numSlices)
                {
                    if(isT2Frame)
                    {
                        t2Volume(Range::all(),Range::all(),imageNexIndex,t2Index,geometricSliceIndexFromHyperBand) = channelCombinedData(Range::all(),Range::all(),unaliasedSliceIndex);
                    }
                    if(isDiffusionFrame)
                    {
                        bValueVolume(Range::all(),Range::all(),imageNexIndex,diffDirIndex,bValueIndex,geometricSliceIndexFromHyperBand) = channelCombinedData(Range::all(),Range::all(),unaliasedSliceIndex);
                    }
                }
            }
        }
    }
}

void Epi::T2SegmentRecon(MDArray::FloatCube& displacementMap, // parasoft-suppress  METRICS-22 "Required for procedual code"
                         MDArray::FloatCube& referenceVolumeForMoCo,
                         const MDArray::Array<float,5>& t2Volume,
                         const Control::ProcessingControlPointer& processingControl,
                         const boost::filesystem::path& dicomDirectory,
                         const Legacy::DicomSeries& dicomSeries)
{
    // Create a trace object to log messages to the console for this rehearsal
    GERecon::Trace trace("T2 segment recon");
    
    // Get slice table for post processing
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    // Distortion correction parameters
    const bool isRpgVolumeAcquired = processingControl->ValueStrict<bool>("IsRpgVolumeAcquired");
    const bool isB0CorrectionEnabled = processingControl->ValueStrict<bool>("B0CorrectionEnabled");
    const bool isRigidCorrectionEnabled = processingControl->ValueStrict<bool>("RigidCorrectionEnabled");
    const bool isAffineCorrectionEnabled = processingControl->ValueStrict<bool>("AffineCorrectionEnabled");
    const int kissoffViews = processingControl->Value<int>("FinalImageKissoffViews");

    // Initialize MoCo class here with rigid correction parameters as T2 images will always use Rigid correction
    const boost::filesystem::path moCoConfigParametersConfigFile = "OrchestraMoCoOxEpiDiffusion.geopts";
    const boost::filesystem::path itkConfigParametersConfigFile = "OrchestraMoCo3pRigid.geopts";
    GERecon::MoCo::MoCoParameters moCoParameters(2);
    GERecon::MoCo::MoCoParametersTuner::PopulateWithConfigFile(moCoParameters, moCoConfigParametersConfigFile, itkConfigParametersConfigFile);
    GERecon::MoCo::MoCo<float,2> moCo(moCoParameters);

    // Initialize Gradwarp. The gradwarp file source will use the gw_coils.dat file stored in the ScanArchive
    // to initialize the gradwarp plugin
    const GradwarpFileSourcePtr gradwarpFileSource = boost::make_shared<GradwarpFileSource>();
    GradwarpPlugin gradwarpPlugin(*processingControl, GERecon::TwoDGradwarp, gradwarpFileSource);

    // Match product scaling - for partial ky scans apply the following scalar
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const float finalImageScaler = homodyneEnabled ? 256.0f / (reconXRes * reconXRes) : 1.0f;

    // For Distortion correction
    MDArray::FloatMatrix reverseImage(reconXRes,reconYRes);
    MDArray::FloatMatrix forwardImage(reverseImage.shape());
    MDArray::FloatMatrix reverseImageSmooth(reverseImage.shape());
    MDArray::FloatMatrix forwardImageSmooth(reverseImage.shape());
    const float fermiRadInput = processingControl->ValueStrict<float>("FermiRadius")/2;
    const float fermiWidthInput = processingControl->ValueStrict<float>("FermiWidth");
    const float fermiEccInput = GERecon::Epi::DistortionCorrection::RpgParameters::FermiEccentricityForImageSmoothing;
    DistortionCorrection::B0Correction b0Correction;
    const DistortionCorrection::RpgParameters rpgParameters;

    const int numT2Images = t2Volume.extent(fourthDim);
    const int numSlices = t2Volume.extent(fifthDim);

    // Workspaces
    MDArray::FloatMatrix imageData(reconXRes,reconYRes);
    MDArray::FloatCube nexCombinedData(reconXRes,reconYRes,numT2Images);
    FloatMatrix uncorrectedMagnitudeImage(reconXRes,reconYRes); 
    FloatMatrix magnitudeImage(reconXRes,reconYRes);
    
    // MoCo API accepts only 3D arrays for uncorrected and corrected images
    MDArray::FloatCube uncorrectedImages;
    MDArray::FloatCube motionCorrectedImages;

    const int numBValues = processingControl->Value<int>("NumBValues");
    int imageNumber = 0;

    const GEHdf5::Snap debugFile("T2SegmentReconScanArchiveRehearsal.h5", "EpiDiffusionScanArchiveRehearsalDebug");
    const int firstNexIndex = 0;
    const int rpgT2Index = 0;
    const int fpgT2Index = 1;
    for(int slice = 0; slice < numSlices; ++slice)
    {
        trace.ConsoleMsg("Reconstructing geometric slice %d", slice);
        std::stringstream sliceNameString;
        sliceNameString << "slice_" << std::setw(3) << std::setfill('0') << slice;
        if(isRpgVolumeAcquired)
        {
            reverseImage = t2Volume(Range::all(),Range::all(),firstNexIndex,rpgT2Index,slice);
            forwardImage = t2Volume(Range::all(),Range::all(),firstNexIndex,fpgT2Index,slice);
        }
        debugFile.Write(reverseImage, sliceNameString.str() + "/reverseImage");
        debugFile.Write(forwardImage, sliceNameString.str() + "/forwardImage");
        nexCombinedData = 0.0f;
        for(int t2Index = 0; t2Index < numT2Images; ++t2Index)
        {
            std::stringstream t2IndexNameString;
            t2IndexNameString << "t2Index_" << std::setw(3) << std::setfill('0') << t2Index;

            const int nunNex = (isRpgVolumeAcquired && (t2Index == 0)) ? 1 : t2Volume.extent(thirdDim); // Reverse polarity is always single NEX
            // Do NEX combination followed by image sending
            for(int nex = 0; nex < nunNex; ++nex)
            {
                std::stringstream nexNameString;
                nexNameString << "nex_" << std::setw(3) << std::setfill('0') << nex;
                imageData = t2Volume(Range::all(),Range::all(),nex,t2Index,slice);
                nexCombinedData(Range::all(),Range::all(),t2Index) += imageData;
                debugFile.Write(imageData, sliceNameString.str() + "/" + t2IndexNameString.str() + "/" + nexNameString.str() + "/imageData");
            }

            // Assure nex count is not zero, and then scale
            if(nunNex != 0)
            {
                nexCombinedData(Range::all(),Range::all(),t2Index) /= (static_cast<float>(nunNex));
            }
            debugFile.Write(nexCombinedData, sliceNameString.str() + "/" + t2IndexNameString.str() + "/nexCombinedData");
        }

        // Distortion Correction. Field Map estimation
        if(isB0CorrectionEnabled)
        {
            referenceVolumeForMoCo(Range::all(),Range::all(),slice) = forwardImage;
            reverseImageSmooth = reverseImage;
            forwardImageSmooth = forwardImage;

            GERecon::Epi::DistortionCorrection::B0Correction::FermiSmooth(reverseImageSmooth, fermiRadInput, fermiWidthInput, fermiEccInput);
            GERecon::Epi::DistortionCorrection::B0Correction::FermiSmooth(forwardImageSmooth, fermiRadInput, fermiWidthInput, fermiEccInput);

            debugFile.Write(reverseImageSmooth, sliceNameString.str() + "/reverseImageSmooth");
            debugFile.Write(forwardImageSmooth, sliceNameString.str() + "/forwardImageSmooth");

            trace.ConsoleMsg(" Distortion Correction: Launch field map computation for Slice: %d", slice);
            MDArray::FloatMatrix displacementMapPerSlice = displacementMap(Range::all(),Range::all(),slice);
            b0Correction.CalculateDisplacementMap(displacementMapPerSlice, forwardImageSmooth, reverseImageSmooth, rpgParameters);

            debugFile.Write(displacementMapPerSlice, sliceNameString.str() + "/displacementMapPerSlice");

            // The computed (raw) B0 maps suffer from some high-frequency artefacts in PE-direction. Hence a Fermi Filter is applied to 
            // suppress those artefacts. Using the combination of a Fermi radius Xres/4 with an eccentricity of 2 results in a smoothing in PE direction only.
            const float fermiRadOutput = static_cast<float>(displacementMap.extent(firstDim))/4;
            const float fermiWidthOutput = GERecon::Epi::DistortionCorrection::RpgParameters::FermiWidthForFieldMapSmoothing;
            const float fermiEccOutput = GERecon::Epi::DistortionCorrection::RpgParameters::FermiEccentricityForFieldMapSmoothing;
            GERecon::Epi::DistortionCorrection::B0Correction::FermiSmooth(displacementMapPerSlice, fermiRadOutput, fermiWidthOutput, fermiEccOutput);
            debugFile.Write(displacementMapPerSlice, sliceNameString.str() + "/displacementMapPerSliceAfterSmoothing");
        }

        const int startingT2Index = isB0CorrectionEnabled ? 1 : 0;
        for(int t2Index = startingT2Index; t2Index < numT2Images; ++t2Index)
        {
            std::stringstream t2IndexNameString;
            t2IndexNameString << "t2Index_" << std::setw(3) << std::setfill('0') << t2Index;
            uncorrectedMagnitudeImage = nexCombinedData(Range::all(),Range::all(),t2Index);
            magnitudeImage = uncorrectedMagnitudeImage;
            if(isRigidCorrectionEnabled || isAffineCorrectionEnabled)
            {
                // Motion correction will be followed by Rpg correction
                if(t2Index > 1)
                {
                    const int thirdDimSizeForMoCoAPI = 1;
                    if(uncorrectedImages.extent(firstDim) != reconXRes)
                    {
                        uncorrectedImages.resize(reconXRes,reconYRes,thirdDimSizeForMoCoAPI);
                    }
                    uncorrectedImages(Range::all(),Range::all(),0) = uncorrectedMagnitudeImage;
                    if(!MDArray::CompareShape(motionCorrectedImages, uncorrectedImages))
                    {
                        motionCorrectedImages.resize(uncorrectedImages.shape());
                    }
                    motionCorrectedImages = 0.0f;
                    moCo.ApplyMoCo(motionCorrectedImages, uncorrectedImages, referenceVolumeForMoCo(Range::all(),Range::all(),slice));
                    uncorrectedMagnitudeImage = motionCorrectedImages(Range::all(),Range::all(),0);
                }
            }

            debugFile.Write(uncorrectedMagnitudeImage, sliceNameString.str() + "/" + t2IndexNameString.str() + "/motioncorrected");

            if(isB0CorrectionEnabled)
            {
                const MDArray::FloatMatrix displacementMapPerSlice = displacementMap(Range::all(),Range::all(),slice);
                DistortionCorrection::B0Correction::ApplyDisplacementMap(magnitudeImage, uncorrectedMagnitudeImage, displacementMapPerSlice, true, Cubic);
            }

            debugFile.Write(magnitudeImage, sliceNameString.str() + "/" + t2IndexNameString.str() + "/rpgcorrected");

            if(kissoffViews > 0)
            {
                magnitudeImage(Range::all(), Range(0, kissoffViews - 1)) = 0.0f;
                magnitudeImage(Range::all(), Range(magnitudeImage.extent(secondDim) - kissoffViews, magnitudeImage.extent(secondDim) - 1)) = 0.0f;
            }

            const SliceCorners sliceCorners = sliceTable.SliceCorners(slice);
            gradwarpPlugin.Run(magnitudeImage, sliceCorners, slice);

            magnitudeImage *= finalImageScaler;

            const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(slice);
            FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

            Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

            ShortMatrix shortImage(rotatedImage.shape());
            shortImage = cast<short>(rotatedImage);

            // T2 images first, then followed by diffusion images
            imageNumber = t2Index*numSlices + slice;
            std::stringstream dicomFileName;
            dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";

            const ImageCorners imageCorners(sliceCorners, sliceOrientation);
            const GEDicom::MR::ImagePointer dicomImage = dicomSeries.NewImage(shortImage, imageNumber, imageCorners);

            FillDiffusionDicomFields(dicomImage, GERecon::DiffusionT2Image, 0, slice, numBValues);

            dicomImage->Save(dicomFileName.str());
        }
    }
}

void Epi::DiffusionSegmentRecon(const MDArray::Array<float,6>& diffusionVolume, // parasoft-suppress  METRICS-22 "Required for procedual code"
                                const MDArray::FloatCube& displacementMap,
                                const MDArray::FloatCube& referenceVolumeForMoCo,
                                const Control::ProcessingControlPointer& processingControl,
                                const boost::filesystem::path& dicomDirectory,
                                const Legacy::DicomSeries& dicomSeries)
{
    // Create a trace object to log messages to the console for this rehearsal
    GERecon::Trace trace("Diffusion Segment Recon");

    // Setup HOEC Correction
    Epi::Diffusion::HoecReconCorrection realtimeFieldAdjustment(*processingControl);
    const bool isIntegratedRefScan = processingControl->Value<bool>("IntegratedReferenceScan");

    // Distortion correction parameters
    const bool isB0CorrectionEnabled = processingControl->Value<bool>("B0CorrectionEnabled");
    const bool isRigidCorrectionEnabled = processingControl->Value<bool>("RigidCorrectionEnabled");
    const bool isAffineCorrectionEnabled = processingControl->Value<bool>("AffineCorrectionEnabled");
    const int kissoffViews = processingControl->Value<int>("FinalImageKissoffViews");
    
    // Initialize MoCo class for diffusion images
    const boost::filesystem::path moCoConfigParametersConfigFile = "OrchestraMoCoOxEpiDiffusion.geopts";
    const boost::filesystem::path itkConfigParametersConfigFile = isRigidCorrectionEnabled ? "OrchestraMoCo3pRigid.geopts" : "OrchestraMoCo3pAffine.geopts";
    GERecon::MoCo::MoCoParameters moCoParameters(2);
    GERecon::MoCo::MoCoParametersTuner::PopulateWithConfigFile(moCoParameters, moCoConfigParametersConfigFile, itkConfigParametersConfigFile);
    GERecon::MoCo::MoCo<float,2> moCo(moCoParameters);

    // Initialize Gradwarp. The gradwarp file source will use the gw_coils.dat file stored in the ScanArchive
    // to initialize the gradwarp plugin
    const GradwarpFileSourcePtr gradwarpFileSource = boost::make_shared<GradwarpFileSource>();
    GradwarpPlugin gradwarpPlugin(*processingControl, GERecon::TwoDGradwarp, gradwarpFileSource);

    // If tensor.dat exists in the Pfile directory and the number of diffusions
    // directions is greater than or equal to 6 then initialize a TensorVectors
    // object to read the tensor vector components from tensor.dat.
    const int tensorFileNumber = processingControl->Value<int>("TensorFileNumber");
    boost::shared_ptr<Epi::Diffusion::TensorVectors> tensorVectorsPtr;
    const int numDiffusionDirections = diffusionVolume.extent(fourthDim);
    if(numDiffusionDirections >= 6)
    {
        tensorVectorsPtr = boost::make_shared<Epi::Diffusion::TensorVectors>(numDiffusionDirections, tensorFileNumber);
    }
    
    // Match product scaling - for partial ky scans apply the following scalar
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const float finalImageScaler = homodyneEnabled ? 256.0f / (reconXRes * reconXRes) : 1.0f;

    // Workspaces and variables for final image sending
    MDArray::FloatMatrix imageData(reconXRes,reconYRes);
    MDArray::FloatMatrix nexCombinedData(reconXRes,reconYRes);
    MDArray::FloatMatrix uncorrectedMagnitudeImage(reconXRes,reconYRes);
    MDArray::FloatMatrix magnitudeImage(reconXRes,reconYRes); 
    MDArray::FloatCube uncorrectedImages;
    MDArray::FloatCube motionCorrectedImages;
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    const int numT2Images = processingControl->Value<int>("NumT2ImagesPerLocation");
    int imageNumber = 0;
    const std::vector<int> diffusionNexTable = processingControl->Value<std::vector<int> >("DiffusionNexTable"); // parasoft-suppress OPT-20 "Procedural, just for rehearsal"
    const int numBValues = diffusionVolume.extent(fifthDim);
    const int numSlices = diffusionVolume.extent(sixthDim);
    const GEHdf5::Snap debugFile("DiffusionSegmentReconScanArchiveRehearsal.h5", "EpiDiffusionScanArchiveRehearsalDebug");
    for(int slice = 0; slice < numSlices; ++slice)
    {
        std::stringstream sliceNameString;
        sliceNameString << "slice_" << std::setw(3) << std::setfill('0') << slice;
        for(int bValueIndex = 0; bValueIndex < numBValues; ++bValueIndex)
        {
            std::stringstream bValueNameString;
            bValueNameString << "bValue_" << std::setw(3) << std::setfill('0') << bValueIndex;
            for(int diffusionDirIndex = 0; diffusionDirIndex < numDiffusionDirections; ++diffusionDirIndex)
            {
                std::stringstream diffDirNameString;
                diffDirNameString << "diffDir_" << std::setw(3) << std::setfill('0') << diffusionDirIndex;

                trace.ConsoleMsg("Reconstructing geometric slice: %d, bvalue: %d, diffusion direction index: %d", slice, bValueIndex, diffusionDirIndex);

                // Do reconstruction and NEX combination
                const int numNex = diffusionNexTable[bValueIndex];
                nexCombinedData = 0.0f;
                for(int nex = 0; nex < numNex; ++nex)
                {
                    std::stringstream nexNameString;
                    nexNameString << "nex_" << std::setw(3) << std::setfill('0') << nex;

                    imageData = diffusionVolume(Range::all(),Range::all(),nex,diffusionDirIndex,bValueIndex,slice);
                    nexCombinedData += imageData;

                    debugFile.Write(imageData, sliceNameString.str() + "/" + bValueNameString.str() + "/" + diffDirNameString.str() + "/" + nexNameString.str() + "/imageData");
                }

                // Divide by number of Nex
                if(numNex != 0)
                {
                    nexCombinedData /= static_cast<float>(numNex);
                }

                debugFile.Write(nexCombinedData, sliceNameString.str() + "/" + bValueNameString.str() + "/" + diffDirNameString.str() + "/nexCombinedData");

                // HOEC
                uncorrectedMagnitudeImage = nexCombinedData;
                const int phaseIndex = numT2Images + bValueIndex*numDiffusionDirections + diffusionDirIndex;
                const int phaseIndexToUse = isIntegratedRefScan ? phaseIndex+1 : phaseIndex;
                realtimeFieldAdjustment.Apply(uncorrectedMagnitudeImage, phaseIndexToUse, slice);

                debugFile.Write(uncorrectedMagnitudeImage, sliceNameString.str() + "/" + bValueNameString.str() + "/" + diffDirNameString.str() + "/hoecOutput");

                magnitudeImage = uncorrectedMagnitudeImage;
                if(isRigidCorrectionEnabled || isAffineCorrectionEnabled)
                {
                    // Motion correction will be followed by Rpg correction
                    const int thirdDimSizeForMoCoAPI = 1;
                    if(uncorrectedImages.extent(firstDim) != reconXRes)
                    {
                        uncorrectedImages.resize(reconXRes,reconYRes,thirdDimSizeForMoCoAPI);
                    }
                    uncorrectedImages(Range::all(),Range::all(),0) = uncorrectedMagnitudeImage;
                    if(!MDArray::CompareShape(motionCorrectedImages, uncorrectedImages))
                    {
                        motionCorrectedImages.resize(uncorrectedImages.shape());
                    }
                    motionCorrectedImages = 0.0f;
                    moCo.ApplyMoCo(motionCorrectedImages, uncorrectedImages, referenceVolumeForMoCo(Range::all(),Range::all(),slice));
                    uncorrectedMagnitudeImage = motionCorrectedImages(Range::all(),Range::all(),0);
                }

                // Distortion Correction
                if(isB0CorrectionEnabled)
                {
                    const MDArray::FloatMatrix displacementMapPerSlice = displacementMap(Range::all(),Range::all(),slice);
                    DistortionCorrection::B0Correction::ApplyDisplacementMap(magnitudeImage, uncorrectedMagnitudeImage, displacementMapPerSlice, true, Cubic);
                }

                debugFile.Write(magnitudeImage, sliceNameString.str() + "/" + bValueNameString.str() + "/" + diffDirNameString.str() + "/distCorrImage");

                if(kissoffViews > 0)
                {
                    magnitudeImage(Range::all(), Range(0,kissoffViews-1)) = 0.0f;
                    magnitudeImage(Range::all(), Range(magnitudeImage.extent(secondDim)-kissoffViews, magnitudeImage.extent(secondDim)-1)) = 0.0f;
                }

                const SliceCorners sliceCorners = sliceTable.SliceCorners(slice);
                gradwarpPlugin.Run(magnitudeImage, sliceCorners, slice);

                magnitudeImage *= finalImageScaler;

                debugFile.Write(magnitudeImage, sliceNameString.str() + "/" + bValueNameString.str() + "/gradwarpedImage");

                const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(slice);
                FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

                Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

                ShortMatrix shortImage(rotatedImage.shape());
                shortImage = cast<short>(rotatedImage);

                imageNumber = numT2Images*numSlices + bValueIndex*numDiffusionDirections*numSlices + slice*numDiffusionDirections + diffusionDirIndex;
                std::stringstream dicomFileName;
                dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";

                const ImageCorners imageCorners(sliceCorners, sliceOrientation);
                const GEDicom::MR::ImagePointer dicomImage = dicomSeries.NewImage(shortImage, imageNumber, imageCorners);

                FillDiffusionDicomFields(dicomImage, GERecon::DiffusionCombinedImage, bValueIndex, slice, numBValues);

                // If this is not a T2 Phase and the tensor vectors object was loaded from
                // tensor.dat then fill the tensor vector values in the DICOM image.
                if(tensorVectorsPtr)
                {
                    FillTensorVectorDicomFields(dicomImage, tensorVectorsPtr->TensorVector(0));
                }

                dicomImage->Save(dicomFileName.str());
            }
        }
    }
}

void Epi::DiffusionCommonRecon(MDArray::FloatCube& imageData,  // parasoft-suppress  METRICS-22 "Required for procedual code"
                               std::vector<int>& unaliasedGeometricSliceNumbers,
                               std::vector<unsigned int>& slicePhaseIndexCounters,
                               const MDArray::ComplexFloatCube& allRawData,
                               RowFlipPlugin& rowFlipPlugin,
                               const boost::shared_ptr<PhaseCorrectionReferenceFile>& phaseCorrectionReferenceFile,
                               const boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin>& staticPhaseCorrectionPlugin,
                               const boost::shared_ptr<Epi::SelfNavDynamicPhaseCorrectionPlugin>& dynamicPhaseCorrectionPlugin,
                               const boost::shared_ptr<Epi::Diffusion::DynamicPhaseCorrectionManager>& dynamicPhaseCorrectionManager,
                               const boost::shared_ptr<RampSamplingPlugin>& rampSamplingPlugin,
                               const Control::ProcessingControlPointer& processingControl,
                               const GEHdf5::Snap& debugFile,
                               const int nexIndex,
                               const int slice,
                               const int phase,
                               const bool isRpgVolume)
{
    // Create a trace object to log messages to the console for this rehearsal
    GERecon::Trace trace("DiffusionRehearsalCommon");
    trace.ConsoleMsg("Processing slice index: %d, phase index: %d", slice, phase);
    std::stringstream phaseNameString;
    std::stringstream sliceNameString;
    phaseNameString << "phase_" << std::setw(3) << std::setfill('0') << phase;
    sliceNameString << "slice_" << std::setw(3) << std::setfill('0') << slice;
    debugFile.Write(allRawData, phaseNameString.str() + "/" + sliceNameString.str() + "/allRawData");

    // First, initialize everything we need
    // Initialize the row flip plugin that is used to apply the EPI row flip. Whether or not
    // to flip a given readout is controlled by the rowflip.param file which is included as one
    // of the saved files in the scan archive.
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const int extraFramesBottom = processingControl->Value<int>("ExtraFramesBottom");
    const int acquiredYRes = processingControl->Value<int>("AcquiredYRes");
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    // Setup a kSpace transformer object. This object handles kSpace apodization as well as partial fourier (if enabled)
    Cartesian2D::KSpaceTransformer transformer(*processingControl); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Initialize ASSET, if enabled 
    ComplexFloatCube homodyneHighPassFilteredDataForAsset;
    ComplexFloatCube homodyneLowPassFilteredDataForAsset;
    ComplexFloatMatrix unaliasedHighPassData;
    ComplexFloatMatrix unaliasedLowPassData;
    ComplexFloatCube aliasedImagesForAsset;
    boost::shared_ptr<Asset::Worker> assetWorker;
    const bool homodyneEnabled = processingControl->ValueStrict<bool>("HalfNex");
    const bool assetEnabled = processingControl->ValueStrict<bool>("Asset");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int numPcCoefficients = processingControl->Value<int>("NumPcCoefficients");

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

    // Initialize Multiband calibration (if needed)
    const bool multibandEnabled = processingControl->Value<bool>("MultibandEnabled");
    boost::shared_ptr<GERecon::Calibration::RawFile> calibrationFile;
    if(multibandEnabled)
    {
        // Initialize multiband calibration
        calibrationFile = boost::make_shared<GERecon::Calibration::RawFile>(*processingControl);
        calibrationFile->TopDirectory(GERecon::Path::InputExamData());
        GERecon::ExamStorageDirectory examStorageDirectory;
        std::vector<boost::filesystem::path> rawFileCalibrationFiles = examStorageDirectory.FindPotentialFiles(calibrationFile->StorageInfo());

        if(rawFileCalibrationFiles.size() > 0)
        {
            const boost::filesystem::path rawFileCalibrationFilePath = rawFileCalibrationFiles[0];
            calibrationFile->Load(rawFileCalibrationFilePath);
        }
        else
        {
            trace.ConsoleMsg("Could not find multiband calibration file!");
            return;
        }
    }

    // Set up variables for Multiband processing
    MDArray::Array<std::complex<float>, 4> multibandUnaliasedSliceData;

    // Initialize SumOfSquares channel combiner to use if SumOfSquares channel combination is enabled
    SumOfSquares channelCombiner(processingControl->ValueStrict<FloatVector>("ChannelWeights")); 

    // If additional non-phase-encoded reference views are acquired, they are either at the top or bottom of the acquired
    // data. Consider these reference views and create a reference to them here, if they are acquired
    ComplexFloatCube rawImageData(allRawData.extent(firstDim), acquiredYRes, numChannels);
    rawImageData = allRawData(Range::all(), Range(extraFramesTop, extraFramesTop + acquiredYRes - 1), Range::all());
    ComplexFloatCube rawReferenceData;
    const int totalNumReferenceViews = extraFramesTop + extraFramesBottom;
    if(totalNumReferenceViews > 0)
    {
        // If reference views are acquired, allocate space for them
        rawReferenceData.resize(allRawData.extent(firstDim), totalNumReferenceViews, numChannels);
        if(extraFramesTop > 0)
        {
            rawReferenceData = allRawData(Range::all(), Range(fromStart, extraFramesTop - 1), Range::all());
        }
        else if(extraFramesBottom > 0)
        {
            rawReferenceData = allRawData(Range::all(), Range(acquiredYRes, extraFramesBottom + acquiredYRes - 1), Range::all());

        }
        debugFile.Write(rawReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/rawReferenceData");
    }
    
    debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/rawImageData");

    // RowFlip for EPI raw data. Flip both the raw image data and raw dynamic reference data
    for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
    {
        ComplexFloatMatrix currentRawImageData = rawImageData(Range::all(), Range::all(), channelIndex);
        rowFlipPlugin.ApplyImageDataRowFlip(currentRawImageData);

        if(rawReferenceData.size() > 0)
        {
            ComplexFloatMatrix currentDynamicReferenceData = rawReferenceData(Range::all(), Range::all(), channelIndex);
            rowFlipPlugin.ApplyReferenceDataRowFlip(currentDynamicReferenceData);
            debugFile.Write(rawReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/rowFlippedRawReferenceData");
        }
    }
    debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/rowFlippedRawImageData");

    // Phase correction in EPI Diffusion pipeline
    const bool assetPhaseCorrectionOptimizationEnabled = processingControl->ValueStrict<bool>("AssetPhaseCorrectionOptimizationEnabled");
    if(staticPhaseCorrectionPlugin)
    {
        // Apply static phase correction
        for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
        {
            const int channelToUseForPhaseCorrection = assetPhaseCorrectionOptimizationEnabled ? phaseCorrectionReferenceFile->MaximumSnrChannel(slice) : channelIndex;
            const FloatVector constantCoefficientsCurrentChannel = phaseCorrectionReferenceFile->ConstantCoefficients(slice, channelToUseForPhaseCorrection);
            const FloatVector linearCoefficientsCurrentChannel = phaseCorrectionReferenceFile->LinearCoefficients(slice, channelToUseForPhaseCorrection);
            ComplexFloatMatrix dataToCorrect = rawImageData(Range::all(), Range::all(), channelIndex);
            staticPhaseCorrectionPlugin->ApplyPhaseCorrection(dataToCorrect, linearCoefficientsCurrentChannel, constantCoefficientsCurrentChannel);
        }
    }
    else if(dynamicPhaseCorrectionPlugin)
    {
        const bool rfChoppingEnabled = processingControl->ValueStrict<bool>("DynamicPCUnchopAlternateNexs");
        const bool inputDataChoppedRelativeToBaseline = (rfChoppingEnabled && (nexIndex % 2 == 1));

        if(inputDataChoppedRelativeToBaseline)
        {
            // If the input data is chopped relative to the baseline data, chop
            // the raw input data to make it consistent with the baseline kSpace data
            rawReferenceData *= -1.0f;
            rawImageData *= -1.0f;
        }

        // Retrieve the static reference coefficients to use
        FloatMatrix linearCoefficients(rawImageData.extent(secondDim), numChannels);
        FloatMatrix constantCoefficients(rawImageData.extent(secondDim), numChannels);
        for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
        {
            const int channelToUseForPhaseCorrection = assetPhaseCorrectionOptimizationEnabled ? phaseCorrectionReferenceFile->MaximumSnrChannel(slice) : channelIndex;
            constantCoefficients(Range::all(), channelIndex) = phaseCorrectionReferenceFile->ConstantCoefficients(slice, channelToUseForPhaseCorrection);
            linearCoefficients(Range::all(), channelIndex) = phaseCorrectionReferenceFile->LinearCoefficients(slice, channelToUseForPhaseCorrection);
        }

        FloatCube pcCoefficients(acquiredYRes, numChannels, numPcCoefficients);
        pcCoefficients = phaseCorrectionReferenceFile->PcCoefficients(slice);

        debugFile.Write(constantCoefficients, phaseNameString.str() + "/" + sliceNameString.str() + "/constantCoefficients");
        debugFile.Write(linearCoefficients, phaseNameString.str() + "/" + sliceNameString.str() + "/linearCoefficients");

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
        const unsigned int currentSlicePhaseIndex = slicePhaseIndexCounters[slice];
        if(currentSlicePhaseIndex == 0)
        {
            // Copy baseline data to a new MDArray and commit that MDArray to be stored
            // for use in all remaining temporal phases of this reconstruction
            baselineImage.resize(rawImageData.shape());
            baselineDynamicReferenceData.resize(rawReferenceData.shape());
            baselineImage = rawImageData;
            baselineDynamicReferenceData = rawReferenceData;
            dynamicPhaseCorrectionManager->CommitBaselineData(slice, baselineImage, baselineDynamicReferenceData);

            // Initialize dynamic phase correction input data from previous to 0.0 for the first (baseline) phase
            // After the baseline phase, this data will come from the previous temporal phase for each slice
            constantPhasePerDynamicReferenceView.resize(rawReferenceData.extent(secondDim), numChannels);
            constantPhasePerDynamicReferenceView = 0.0f;
            // The dynamicCoefficients from the previous phase were initialized to 0.0 above
        }
        else
        {
            // Retrieve baseline data for the current slice
            baselineImage.reference(dynamicPhaseCorrectionManager->BaselineImage(slice));
            baselineDynamicReferenceData.reference(dynamicPhaseCorrectionManager->BaselineDynamicReferenceViews(slice));

            // Retrieve data from previous temporal phase for the current geometric slice. This data
            // is input data to the dynamic phase correction algorithm
            const int previousPhaseIndex = currentSlicePhaseIndex - 1;
            constantPhasePerDynamicReferenceView.reference(dynamicPhaseCorrectionManager->ConstantPhasePerDynamicReferenceView(previousPhaseIndex, slice));
            dynamicCoefficients = dynamicPhaseCorrectionManager->DynamicCoefficients(previousPhaseIndex, slice);
        }
        debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/transformedRowFlippedImageData");
        debugFile.Write(rawReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/transformedRowFlippedReferenceData");
        debugFile.Write(baselineImage, phaseNameString.str() + "/" + sliceNameString.str() + "/transformedBaselineImageData");
        debugFile.Write(baselineDynamicReferenceData, phaseNameString.str() + "/" + sliceNameString.str() + "/transformedBaselineReferenceData");
        debugFile.Write(constantPhasePerDynamicReferenceView, phaseNameString.str() + "/" + sliceNameString.str() + "/inputConstantPhasePerView");
        const int transformYRes = processingControl->Value<int>("TransformYRes");
        debugFile.Write(SelfNavDynamicPhaseCorrectionBase::B0DriftInPixels(dynamicCoefficients[0](0), transformYRes), phaseNameString.str() + "/" + sliceNameString.str() + "/inputB0DynamicCoefInPixels");
        debugFile.Write(dynamicCoefficients[0], phaseNameString.str() + "/" + sliceNameString.str() + "/inputDynamicCoefficientDeltas");
        
        const int shotIndex = 0;
        dynamicPhaseCorrectionPlugin->ApplyPhaseCorrection(constantPhasePerDynamicReferenceView, dynamicCoefficients, rawImageData, rawReferenceData, baselineImage, baselineDynamicReferenceData, pcCoefficients, shotIndex);

        debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/outputPhaseCorrectedData");

        // Commit data from the current temporal phase for use in the next temporal phase
        dynamicPhaseCorrectionManager->CommitDynamicData(currentSlicePhaseIndex, slice, dynamicCoefficients, constantPhasePerDynamicReferenceView);

        // Increment slice phase index counter
        slicePhaseIndexCounters[slice] = slicePhaseIndexCounters[slice] + 1;

        // Transform phase corrected data back to kx,ky for apodization filtering and 2D kSpace transform
        Fourier::Fft(rawImageData, MDArray::firstDim);
        rawImageData(Range(fromStart, toEnd, 2), Range::all(), Range::all()) *= -1.0f;

        debugFile.Write(SelfNavDynamicPhaseCorrectionBase::B0DriftInPixels(dynamicCoefficients[0](0), transformYRes), phaseNameString.str() + "/" + sliceNameString.str() + "/outputB0DynamicCoefInPixels");
        debugFile.Write(constantPhasePerDynamicReferenceView, phaseNameString.str() + "/" + sliceNameString.str() + "/outputConstantPhasePerView");
        debugFile.Write(dynamicCoefficients[0], phaseNameString.str() + "/" + sliceNameString.str() + "/outputDynamicCoefficientDeltas");
        debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/finalTransformedOutput");

    }
    debugFile.Write(rawImageData, phaseNameString.str() + "/" + sliceNameString.str() + "/phaseCorrectedRawImageData");

    const bool flipDataInPhaseEncodeDirection = processingControl->ValueStrict<bool>("PhaseFlip");
    const bool flipToUse = isRpgVolume ? (!flipDataInPhaseEncodeDirection) : flipDataInPhaseEncodeDirection; 
    const int logicalSlice = sliceTable.AcquiredSliceInfo(slice).second;

    // Multiband processing
    if(multibandEnabled) //for Multiband scans
    {
        // Need to handle RPG and FPG calibration if Rpg volume is acquired.
        // Create pointers for Multiband Calibration and Multiband Worker
        GERecon::Arc::CalibrationPtr calPtr;
        GERecon::Arc::CalibrationPtr rpgCalPtr;
        GERecon::Epi::MultibandCalibrationProcessor multibandCalibrationProcessor(*processingControl);
        Epi::MultibandWorker multibandWorker(*processingControl);

        const bool storeCalibration = isRpgVolume ? false : true;
        GERecon::Epi::MultibandCalibration(calPtr, multibandCalibrationProcessor, calibrationFile, logicalSlice, storeCalibration, flipToUse);
        
        if(isRpgVolume)
        {
            GERecon::Epi::MultibandCalibration(rpgCalPtr, multibandCalibrationProcessor, calibrationFile, logicalSlice, storeCalibration, flipToUse);
        }

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

        debugFile.Write(resampledKSpaceToUnalias, phaseNameString.str() + "/" + sliceNameString.str() + "/resampledKSpaceToUnalias");

        // Perform Multiband Unaliasing
        if(isRpgVolume)
        {
            multibandWorker.UnaliasSliceData(multibandUnaliasedSliceData, unaliasedGeometricSliceNumbers, resampledKSpaceToUnalias, rpgCalPtr, logicalSlice);
        }
        else
        {
            multibandWorker.UnaliasSliceData(multibandUnaliasedSliceData, unaliasedGeometricSliceNumbers, resampledKSpaceToUnalias, calPtr, logicalSlice);
        }
        debugFile.Write(multibandUnaliasedSliceData, phaseNameString.str() + "/" + sliceNameString.str() + "/multibandUnaliasedSliceData");
    }
    else // for Non-Multiband scans
    {
        const int numUnaliasedSlicesForNonMultiband = 1;
        unaliasedGeometricSliceNumbers.resize(numUnaliasedSlicesForNonMultiband);
        unaliasedGeometricSliceNumbers[0] = slice;
    }

    // Create MDArray workspaces
    ComplexFloatMatrix accumulatedChannels(reconXRes, reconYRes);
    ComplexFloatMatrix transformedData(reconXRes, reconYRes);
    // Resize the array once
    if(imageData.extent(thirdDim) != boost::numeric_cast<int>(unaliasedGeometricSliceNumbers.size()))
    {
        imageData.resize(reconXRes, reconYRes, boost::numeric_cast<int>(unaliasedGeometricSliceNumbers.size()));
    }
    const int transformKissoffViews = processingControl->ValueStrict<int>("TransformKissoffViews");
    for(int unaliasedSliceIndex = 0; unaliasedSliceIndex < boost::numeric_cast<int>(unaliasedGeometricSliceNumbers.size()); ++unaliasedSliceIndex)
    {
        const int totalgeometricslices = sliceTable.GeometricSliceLocations();
        const int geometricSliceNumber = unaliasedGeometricSliceNumbers[unaliasedSliceIndex];

        if(geometricSliceNumber < totalgeometricslices) // Discarding extended slices for Multiband scans
        {
            std::stringstream unaliasedSliceNameString;
            unaliasedSliceNameString << "slice_" << std::setw(3) << std::setfill('0') << geometricSliceNumber;

            if(!assetEnabled)
            {
                // Use a sum of squares channel combiner
                channelCombiner.Reset();
            }

            // Loop over all channels and transform the data
            for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
            {
                const ComplexFloatMatrix kSpaceToTransform = multibandEnabled ? multibandUnaliasedSliceData(Range::all(), Range::all(), unaliasedSliceIndex, channelIndex): (rampSamplingPlugin ? rampSamplingPlugin->ApplyRampSampling(rawImageData(Range::all(), Range::all(), channelIndex)) : rawImageData(Range::all(), Range::all(), channelIndex));

                if(assetEnabled && !homodyneEnabled) // parasoft-suppress  METRICS-23 "rehearsal code handles many cases in a single int main pipeline"
                {
                    ComplexFloatMatrix aliasedChannelData = aliasedImagesForAsset(Range::all(), Range::all(), channelIndex);
                    transformer.Apply(aliasedChannelData, kSpaceToTransform, flipToUse);
                    debugFile.Write(aliasedChannelData, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/aliasedChannelData");
                }
                else if(homodyneEnabled && assetEnabled) // parasoft-suppress  METRICS-23 "Procedural rehearsal code handles many cases in a single int main pipeline"
                {
                    ComplexFloatMatrix highPassImage = homodyneHighPassFilteredDataForAsset(Range::all(), Range::all(), channelIndex);
                    ComplexFloatMatrix lowPassImage = homodyneLowPassFilteredDataForAsset(Range::all(), Range::all(), channelIndex);
                    transformer.Apply(highPassImage, lowPassImage, kSpaceToTransform, flipToUse);
                    debugFile.Write(highPassImage, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/aliasedHighPassImage");
                    debugFile.Write(lowPassImage, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/aliasedLowPassImage");
                }
                else if(!assetEnabled) // parasoft-suppress  METRICS-23 "Procedural, for-loop rehearsal code"
                {
                    // Handles the homodyne only, zerofilling, or full kx,ky cases
                    transformer.Apply(transformedData, kSpaceToTransform, flipToUse);
                    DiscardKissoffViews(transformedData,transformKissoffViews);
                    channelCombiner.Accumulate(transformedData, channelIndex);
                }
            }

            // Compute the channel combined image
            const GERecon::SliceCorners sliceCornersForAsset = sliceTable.AcquiredSliceCorners(geometricSliceNumber);
            if(assetEnabled && homodyneEnabled)
            {
                debugFile.Write(homodyneHighPassFilteredDataForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/homodyneHighPassFilteredDataForAsset");
                debugFile.Write(homodyneLowPassFilteredDataForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/homodyneLowPassFilteredDataForAsset");
                assetWorker->Unalias(unaliasedHighPassData, unaliasedLowPassData, homodyneHighPassFilteredDataForAsset, homodyneLowPassFilteredDataForAsset, geometricSliceNumber, sliceCornersForAsset);
                accumulatedChannels = unaliasedHighPassData;
                Cartesian2D::Homodyne::ApplyPhaseCorrection(accumulatedChannels, unaliasedLowPassData);
                DiscardKissoffViews(accumulatedChannels,transformKissoffViews);
            }
            else if(assetEnabled && !homodyneEnabled)
            {   
                debugFile.Write(aliasedImagesForAsset, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/aliasedImagesForAsset");
                assetWorker->Unalias(accumulatedChannels, aliasedImagesForAsset, geometricSliceNumber, sliceCornersForAsset);
                DiscardKissoffViews(accumulatedChannels,transformKissoffViews);
            }
            else if(!assetEnabled)
            {
                accumulatedChannels = channelCombiner.GetCombinedImage();
                channelCombiner.Reset();
            }
            debugFile.Write(accumulatedChannels, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/accumulatedChannels");
            imageData(Range::all(),Range::all(),unaliasedSliceIndex) = abs(accumulatedChannels);
            debugFile.Write(imageData, phaseNameString.str() + "/" + sliceNameString.str() + "/" + unaliasedSliceNameString.str() + "/imageData");
        }
    }
}

void GERecon::Epi::DiscardKissoffViews(const MDArray::ComplexFloatMatrix& imageData,
                                       const int transformKissoffViews)
{
    if(transformKissoffViews <= imageData.extent(secondDim) && transformKissoffViews > 0)
    {
        imageData(Range::all(), Range(0,transformKissoffViews-1)) = 0;
        imageData(Range::all(), Range(imageData.extent(secondDim) - transformKissoffViews, imageData.extent(secondDim)-1)) = 0;
    }
}

void GERecon::Epi::MultibandCalibration(GERecon::Arc::CalibrationPtr& calibrationPtr,
                                        GERecon::Epi::MultibandCalibrationProcessor& multibandCalibrationProcessor,
                                        const boost::shared_ptr<GERecon::Calibration::RawFile>& calibrationFile,
                                        const int logicalSliceNumber,
                                        const bool storeCalibration,
                                        const bool flipDataInPhaseEncodeDirection)
{
    
    // Prepare a single slice of calibration data and its associated parameters so that we can call Calibration() later
    const float calScanCenter = calibrationFile->Info().ScanCenter();
    const MDArray::Array<std::complex<float>, 4> imageSpaceToCalibrateWith = multibandCalibrationProcessor.CalibrationImageSpace(calibrationFile->ImageSpaceData(GERecon::Calibration::SurfaceImageSpaceAllPassSets), 
        calibrationFile->Info().FirstSliceCorners(),
        calibrationFile->Info().LastSliceCorners(),
        calScanCenter,
        logicalSliceNumber,
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
    const bool distribute = false;
    const int leftBlanks = multibandCalibrationProcessor.XCalPointsToSkipLeft();
    const int rightBlanks = multibandCalibrationProcessor.XCalPointsToSkipRight();
    const int unacceleratedDimKernelSize = multibandCalibrationProcessor.XKernelSize();
    const int acceleratedDimKernelSize = multibandCalibrationProcessor.YKernelSize();
    const int maxCalXSize = kSpaceToCalibrateWith.extent(MDArray::firstDim);
    const int maxCalYSize = kSpaceToCalibrateWith.extent(MDArray::secondDim);
    const int maxCalZSize = kSpaceToCalibrateWith.extent(MDArray::thirdDim);

    // The following function call will run ARC calibration and store the calibration in the calibration manager
    const int indexToStoreCalibration = logicalSliceNumber;
    calibrationPtr = storeCalibration ? Arc::CalibrationManager::Instance()->Calibration(kSpaceToCalibrateWith, calAcquiredLocations, undersampledDataAcquiredLocations, undersampledDataDimensionSizes,
        yAcceleration, zAcceleration, xPeakPosition, yPeakPosition, zPeakPosition, indexToStoreCalibration, distribute, leftBlanks, rightBlanks, unacceleratedDimKernelSize, acceleratedDimKernelSize,
        maxCalXSize, maxCalYSize, maxCalZSize) : 
    Arc::Calibrate(kSpaceToCalibrateWith, calAcquiredLocations, undersampledDataAcquiredLocations, undersampledDataDimensionSizes, yAcceleration, zAcceleration, xPeakPosition, yPeakPosition, zPeakPosition,
        distribute, leftBlanks, rightBlanks, unacceleratedDimKernelSize, acceleratedDimKernelSize, maxCalXSize, maxCalYSize, maxCalZSize);
}

bool GERecon::Epi::GetPacketInfo(int& frameType,
                                 int& t2Index,
                                 int& bValueIndex,
                                 int& diffDirIndex,
                                 const int packetCount,
                                 const int numSlices,
                                 const int totalNumRef,
                                 const int totalNumT2,
                                 const int numBValues,
                                 const int numDiffusionDirections,
                                 const std::vector<int>& diffusionNexTable)
{
    t2Index = 0; // t2Index > 0 only for tensor Diffusion (DiffusionHyperScanOpcode)

    if(packetCount <= totalNumRef)
    {
        frameType = 0;
        bValueIndex = 0;
        diffDirIndex = 0;
    }
    else if(packetCount <= totalNumRef + totalNumT2)
    {
        frameType = 1;
        bValueIndex = 0;
        diffDirIndex = 0;
    }
    else
    {
        frameType = 2;
        
        // Specify the diffusion direction index (0-based)
        diffDirIndex = (packetCount - 1) % numDiffusionDirections;

        // Specify the bValue index (0-based)
        int totalNumTillCurrentB = totalNumRef + totalNumT2;
        for(int i=0; i < numBValues; i++)
        {
            totalNumTillCurrentB += numSlices * diffusionNexTable[i];
            if(packetCount <= totalNumTillCurrentB)
            {
                bValueIndex = i;
                i = numBValues;
            }
        }
        if(packetCount > totalNumTillCurrentB)
        {
            return false;
        }
    }

    return true;
}