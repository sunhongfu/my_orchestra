// Copyright 2018 General Electric Company.  All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <iostream>
#include <string>

#include <MDArray/Fourier.h>
#include <MDArray/MDArray.h>
#include <MDArray/Utils.h>

#include <Orchestra/Acquisition/ControlPacket.h>
#include <Orchestra/Acquisition/ControlTypes.h>
#include <Orchestra/Acquisition/FrameControl.h>

#include <Orchestra/Acquisition/Core/ArchiveStorage.h>

#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Common/SliceInfoTable.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Epi/LxControlSource.h>
#include <Orchestra/Epi/PhaseCorrectionReferenceFile.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrectionPlugin.h>
#include <Orchestra/Epi/NearestNeighborStaticPhaseCorrection.h>
#include <Orchestra/Epi/RowFlipPlugin.h>


#include <Orchestra/Legacy/CommonTypes.h>
#include <Orchestra/Legacy/LxDownloadData.h>
#include <Orchestra/Legacy/Pfile.h>

#include "EpiReferenceScanRecon.h"

using namespace GERecon;
using namespace GERecon::Epi;
using namespace MDArray;

/**
* EPI Reference Scan
*
* @author Eric Printz
*/
bool GERecon::Epi::EpiReferenceScanRecon(const GERecon::Legacy::PfilePointer& pfile)
{    
    // Create a ProcessingControl from the Legacy Pfile.
    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfile->CreateOrchestraProcessingControl<Epi::LxControlSource>();

    // Pull info directly out of ProcessingControl using the Value function.
    // Values are checked and validated on creation.
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const float rationalScalingFactor = processingControl->Value<float>("RationalScalingFactor");
    const unsigned int numShots = processingControl->Value<unsigned int>("NumberOfShots");
    const bool doConstantPhaseAlignment = processingControl->Value<bool>("PhaseAlignmentEnabled");
    const MDArray::FloatVector channelWeights = processingControl->ValueStrict<MDArray::FloatVector>("NormalizedChannelWeights");
    const int middleKSpaceView = processingControl->ValueStrict<int>("ViewToUseForMaxChannelSnrComputation");

    // Pull slice and channel count from pfile
    const size_t numSlices = pfile->SliceCount();
    const size_t numChannels = pfile->ChannelCount();

    const NearestNeighborStaticPhaseCorrectionPlugin nnPlugin(*processingControl);

    PhaseCorrectionReferenceFile epiRef(*processingControl, PhaseCorrectionReferenceFile::WriteMode);
    FloatVector linearCoefs(acqYRes);
    FloatVector constantCoefs(acqYRes);

    // Create a static phase correction file object to compute and write the phase correction maps to an HDF5 file
    IntVector maximumSnrChannels(boost::numeric_cast<int>(numSlices)); //vector containing the zero based channel index of the channel with max SNR for each slice
    FloatVector channelSnr(boost::numeric_cast<int>(numChannels));

    for(unsigned int slice = 0; slice < numSlices; ++slice)
    {           
        const unsigned int currentEcho = 0;
        channelSnr = 0;

        FloatMatrix constantCoefficientsAllChannels(acqYRes, boost::numeric_cast<int>(numChannels));
        for(unsigned int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
        {
            std::cout << "Processing Slice " << slice << " : Channel " << currentChannel << std::endl;

            // Extract all raw data from the pfile for the current slice/echo/channel
            ComplexFloatMatrix allRawData = pfile->KSpaceData<float>(slice, currentEcho, currentChannel);

            // If dynamic reference views are acquired, do not use them when computing the static reference coefficients
            // The following line discards the dynamic reference views, if acquired, from the raw data matrix
            ComplexFloatMatrix rawData = allRawData(Range::all(), Range(extraFramesTop, extraFramesTop + acqYRes - 1));

            // FFT each readout to x,ky space. Chop to center the transform
            rawData(Range(fromStart,toEnd,2), Range::all()) *= -1.0f;
            Fourier::Ifft(rawData, firstDim); 
            rawData(Range(fromStart,toEnd,2), Range::all()) *= -1.0f;
            rawData *= rationalScalingFactor;

            nnPlugin.ComputeNearestNeighborStaticCoefficients(linearCoefs, constantCoefs, rawData, numShots);

            epiRef.LinearCoefficients(linearCoefs, slice, currentChannel);
            epiRef.ConstantCoefficients(constantCoefs, slice, currentChannel);
            constantCoefficientsAllChannels(Range::all(), currentChannel) = constantCoefs;

            // SNR = mean(sqrt(mag. squared)) of middle view in kspace
            channelSnr(currentChannel) = mean(abs(rawData(Range::all(), middleKSpaceView))) / channelWeights(currentChannel);
        }

        if(doConstantPhaseAlignment)
        {
            NearestNeighborStaticPhaseCorrection::PhaseAlignment(constantCoefficientsAllChannels);
            for(unsigned int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
            {
                const FloatVector constantCoefs = constantCoefficientsAllChannels(Range::all(), currentChannel);
                epiRef.ConstantCoefficients(constantCoefs, slice, currentChannel);
            }
        }

        // The channel with maximum SNR here is the index of the channel with the greatest value
        maximumSnrChannels(slice) = FindMaxPosition(channelSnr);
    }

    //Save ref.h5
    epiRef.MaximumSnrChannels(maximumSnrChannels);
    epiRef.Write();

    return true;
}

bool GERecon::Epi::EpiReferenceScanRecon(const GERecon::ScanArchivePointer& scanArchive) // parasoft-suppress  METRICS-22 "Required for procedual code"
{
    // Create a trace object to log messages to the console for this rehearsal
    GERecon::Trace trace("EpiReferenceScan-ScanArchiveRehearsal");

    // Create an Epi LxControlSource object to interpret recon parameters. The LxControlSource object is used
    // to interpret parameters from the pool header and store the parameters in a processing control object.
    // The processing control object can be used to 
    const GERecon::Legacy::LxDownloadDataPointer downloadData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchive->LoadDownloadData());
    const boost::shared_ptr<GERecon::Epi::LxControlSource> controlSource = boost::make_shared<GERecon::Epi::LxControlSource>(downloadData);
    const Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

    // Prior to initializing external file sources, extract all files saved in the scan archive
    // Because files contained in the scan archive are extracted to the locations of their corresponding
    // GERecon::Path's, set the GERecon::Path locations prior to loading the saved files
    const boost::filesystem::path scanArchiveFullPath = scanArchive->Path();
    GERecon::Path::InputAppData(scanArchiveFullPath.parent_path());
    GERecon::Path::InputExamData(scanArchiveFullPath.parent_path());
    GERecon::Path::ScannerConfig(scanArchiveFullPath.parent_path());
    scanArchive->LoadSavedFiles();

    // Pull info directly out of ProcessingControl using the Value function.
    // Values are checked and validated on creation.
    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int extraFramesTop = processingControl->Value<int>("ExtraFramesTop");
    const int extraFramesBottom = processingControl->Value<int>("ExtraFramesBottom");
    const float rationalScalingFactor = processingControl->Value<float>("RationalScalingFactor");
    const unsigned int numShots = processingControl->Value<unsigned int>("NumberOfShots");
    const bool doHighOrderPhaseCorrection = processingControl->Value<bool>("ApplyHighOrderPhaseCorrection");
    const int numCoefficients = doHighOrderPhaseCorrection ? 4 : 2;
    const bool doConstantPhaseAlignment = processingControl->Value<bool>("PhaseAlignmentEnabled");
    const MDArray::FloatVector channelWeights = processingControl->ValueStrict<MDArray::FloatVector>("NormalizedChannelWeights");
    const int middleKSpaceView = processingControl->ValueStrict<int>("ViewToUseForMaxChannelSnrComputation");
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const bool multibandEnabled = processingControl->Value<bool>("MultibandEnabled");
    const int numSlices = multibandEnabled ? processingControl->ValueStrict<int>("MultibandNumAcquiredSlices") : sliceTable.GeometricSliceLocations();
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int numAcquisitionsPerRepetition = processingControl->Value<int>("NumAcquisitionsPerRepetition");
    const int totalNumReferenceViews = extraFramesTop + extraFramesBottom;
    const int totalNumberOfAcquiredViews = totalNumReferenceViews + acqYRes;
    const RowFlipParametersPointer rowFlipParams = boost::make_shared<RowFlipParameters>(totalNumberOfAcquiredViews);
    RowFlipPlugin rowFlipPlugin(rowFlipParams, *processingControl); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    const NearestNeighborStaticPhaseCorrectionPlugin nnPlugin(*processingControl);

    PhaseCorrectionReferenceFile epiRef(*processingControl, PhaseCorrectionReferenceFile::WriteMode);
    FloatVector linearCoefs(acqYRes);
    FloatVector constantCoefs(acqYRes);
    FloatMatrix pcCoefs(acqYRes, numCoefficients);
    FloatCube pcCoefsAllChannel(acqYRes, numCoefficients, numChannels);

    // Create a static phase correction file object to compute and write the phase correction maps to an HDF5 file
    IntVector maximumSnrChannels(numSlices); //vector containing the zero based channel index of the channel with max SNR for each slice
    maximumSnrChannels = 0;
    FloatVector channelSnr(numChannels);

    // Load the archive storage which contains all acquisition data held in the archive
    GERecon::Acquisition::ArchiveStoragePointer archiveStorage = GERecon::Acquisition::ArchiveStorage::Create(scanArchive);

    // Determine how many control (DAB) packets are contained in the storage
    const size_t numControls = archiveStorage->AvailableControlCount();

    // Initialize an acquisition (pass) counter. This counter is incremented every time an acquisition (pass) done
    // packet is received from the PSD. Each time the acquisitionPassCounter reaches the numAcquisitionsPerRepetition value
    int acquisitionPassCounter = 0;

    // Diffusion hyperscan packet (opcode = DiffusionHyperScanOpcode) has all the information that
    // recon needes to sort data into integrate references, T2, BValues, Diffusion Directions, etc...).
    // The EpiOpcode does not have fields to identify the various diffusion passes. Thus, need to sort
    // manually. Integrated references are always scanned first.
    bool refPassComplete = false;

    // Loop over all the control packets and do reconstruction only if it is a reference scan
    ComplexFloatCube allRawData(acqXRes,totalNumberOfAcquiredViews,numChannels);
    for(size_t controlPacketIndex = 0; (controlPacketIndex < numControls) && (!refPassComplete); ++controlPacketIndex)
    {
        const GERecon::Acquisition::FrameControlPointer controlPacketAndFrameData = archiveStorage->NextFrameControl();

        // Variables extracted from acquisition packets
        int sliceIndexInAcq = 0;
        short viewIncrement = 0;
        bool isRefFrame = false; // This is used only if there is a integrated ref scan

        if(controlPacketAndFrameData->Control().Opcode() == Acquisition::ScanControlOpcode)
        {
            // For diffusion scans, the same volume may be acquired multiple times for the various diffusion passes (T2, B-Value/Diffusion Direction). 
            // Thus, when this counter reaches the number of acquisitions per repetition, reset the counter back to 0.
            ++acquisitionPassCounter;
            if(acquisitionPassCounter == numAcquisitionsPerRepetition)
            {
                acquisitionPassCounter = 0;
                refPassComplete = true;
            }

            trace.ConsoleMsg("Received scan control packet");
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::EpiOpcode)
        {
            const Acquisition::HyperFrameControlPacket hyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::HyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(hyperFramePacket.sliceNumH, hyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(hyperFramePacket.viewSkipH, hyperFramePacket.viewSkipL));
            isRefFrame = 1;
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::DiffusionHyperScanOpcode)
        {
            const Acquisition::DiffusionHyperFrameControlPacket diffusionHyperFramePacket = controlPacketAndFrameData->Control().Packet().As<Acquisition::DiffusionHyperFrameControlPacket>();
            sliceIndexInAcq = Acquisition::GetPacketValue(diffusionHyperFramePacket.sliceNumH, diffusionHyperFramePacket.sliceNumL);
            viewIncrement = static_cast<short>(Acquisition::GetPacketValue(diffusionHyperFramePacket.viewSkipH, diffusionHyperFramePacket.viewSkipL));
            isRefFrame = (0 == diffusionHyperFramePacket.frameType);
        }
        else
        {
            trace.ConsoleMsg("Ignoring unsupported packet with opcode: %d",controlPacketAndFrameData->Control().Opcode());
            isRefFrame = false;
        }
        if(isRefFrame)
        {
            // Use the slice table to determine which geometric slice index this data corresponds to
            const int geometricSliceIndex = sliceTable.GeometricSliceNumber(acquisitionPassCounter, sliceIndexInAcq);
            const ComplexFloatCube rawDataFromApi = controlPacketAndFrameData->Data();

            // The API returns data with acqXRes x numChannels x acqYRes as dimensions. So, we need to tranpose the second and third dimensions to get acqXRes x acqYRes x numChannels
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

            ComplexFloatCube rawImageData = allRawData(Range::all(), Range(extraFramesTop, extraFramesTop + acqYRes - 1), Range::all());

            // RowFlip for EPI raw data.
            for(int channelIndex = 0; channelIndex < numChannels; ++channelIndex)
            {
                ComplexFloatMatrix currentRawImageData = rawImageData(Range::all(), Range::all(), channelIndex);
                rowFlipPlugin.ApplyImageDataRowFlip(currentRawImageData);
            }

            FloatMatrix constantCoefficientsAllChannels(acqYRes, numChannels);
            for(unsigned int currentChannel = 0; currentChannel < boost::numeric_cast<unsigned int>(numChannels); ++currentChannel)
            {
                std::cout << "Processing Slice " << geometricSliceIndex << " : Channel " << currentChannel << std::endl;

                ComplexFloatMatrix rawData = rawImageData(Range::all(),Range::all(),currentChannel);

                // FFT each readout to x,ky space. Chop to center the transform
                rawData(Range(fromStart,toEnd,2), Range::all()) *= -1.0f;
                Fourier::Ifft(rawData, firstDim); 
                rawData(Range(fromStart,toEnd,2), Range::all()) *= -1.0f;
                rawData *= rationalScalingFactor;

                if( doHighOrderPhaseCorrection )
                {
                    nnPlugin.ComputeNearestNeighborStaticCoefficients(pcCoefs, rawData, numShots);
                }
                else
                {
                    nnPlugin.ComputeNearestNeighborStaticCoefficients(linearCoefs, constantCoefs, rawData, numShots);
                    pcCoefs(Range::all(), firstDim) = linearCoefs;
                    pcCoefs(Range::all(), secondDim) = constantCoefs;
                }

                epiRef.PcCoefficients(pcCoefs, geometricSliceIndex, currentChannel);
                constantCoefficientsAllChannels(Range::all(), currentChannel) = pcCoefs(Range::all(), numCoefficients - 1);
                pcCoefsAllChannel(Range::all(), Range::all(), currentChannel) = pcCoefs;

                // SNR = mean(sqrt(mag. squared)) of middle view in kspace
                channelSnr(currentChannel) = mean(abs(rawData(Range::all(), middleKSpaceView))) / channelWeights(currentChannel);
            }

            if(doConstantPhaseAlignment)
            {
                NearestNeighborStaticPhaseCorrection::PhaseAlignment(constantCoefficientsAllChannels);
                for(unsigned int currentChannel = 0; currentChannel < boost::numeric_cast<unsigned int>(numChannels); ++currentChannel)
                {
                    pcCoefs(Range::all(), Range(0, numCoefficients-2)) = pcCoefsAllChannel(Range::all(), Range(0, numCoefficients-2), currentChannel);
                    pcCoefs(Range::all(), numCoefficients - 1) = constantCoefficientsAllChannels(Range::all(), currentChannel);
                    epiRef.PcCoefficients(pcCoefs, geometricSliceIndex, currentChannel);
                }
            }

            // The channel with maximum SNR here is the index of the channel with the greatest value
            maximumSnrChannels(geometricSliceIndex) = FindMaxPosition(channelSnr);

        }
    }

    //Save ref.h5
    epiRef.MaximumSnrChannels(maximumSnrChannels);
    epiRef.Write();
    std::cout << "ref.h5 file was written out" << std::endl;

    return true;

}
