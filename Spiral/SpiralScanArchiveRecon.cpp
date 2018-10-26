// Copyright 2018 General Electric Company.  All rights reserved.

#include <Hdf5/Snap.h>

#include <Orchestra/Acquisition/ControlPacket.h>
#include <Orchestra/Acquisition/ControlTypes.h>
#include <Orchestra/Acquisition/Core/ArchiveStorage.h>
#include <Orchestra/Acquisition/DataTypes.h>
#include <Orchestra/Acquisition/FrameControl.h>

#include <Orchestra/Cartesian3D/ZTransformer.h>
#include <Orchestra/Common/ExamData.h>
#include <Orchestra/Common/ExamDataStorageInfo.h>
#include <Orchestra/Common/ExamStorageDirectory.h>
#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ProgramOptions.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Control/ProcessingControl.h>

// Algorithms
#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>
#include <Orchestra/Core/Windows.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Spiral/LxControlSource.h>
#include <Orchestra/Spiral/SpiralWorker.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/LxDownloadData.h>

#include "Rehearsal.h"
#include "SpiralUtil.h"

using namespace GERecon;
using namespace GERecon::Spiral;
using namespace GERecon::Legacy;
using namespace MDArray;

void GERecon::SpiralScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive) // parasoft suppress item METRICS-22 reason "High line count is required for algorithm."
{
    // Create a trace object to log messages to the console for this reeharsal
    GERecon::Trace trace("SpiralScanArchiveRehearsal");
   
    // Initialize a directory to store the dicom image outputs from this rehearsal
    const boost::filesystem::path scanArchiveFullPath = scanArchive->Path();

    // Create a Spiral LxControlSource object to interpret recon parameters. The LxControlSource object is used
    // to interpret parameters from the pool header and store the parameters in a processing control object.
    // The processing control object can be used to 
    const GERecon::Legacy::LxDownloadDataPointer downloadData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchive->LoadDownloadData());
    const boost::shared_ptr<GERecon::Spiral::LxControlSource> controlSource = boost::make_shared<GERecon::Spiral::LxControlSource>(downloadData);
    const Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

    // Initialize DicomSeries object to use when creating dicom images
    const Legacy::DicomSeries dicomSeries(downloadData); // parasoft-suppress  OPT-20 "Procedural rehearsal code"

    // Extract parameters from processing control 
    const bool is3DASL = processingControl->Value<bool>("Is3DASL");
    const GERecon::SliceInfoTable sliceTable = processingControl->ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const int numAcquisitions = processingControl->Value<int>("NumAcquisitions");
    const int numPointsPerArm = processingControl->Value<int>("NumPointsPerArm");
    const int numArms = processingControl->Value<int>("NumArms");
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int numTotalSlices = processingControl->Value<int>("NumSlices");
    const int numSlices = (numTotalSlices + (numAcquisitions-1))/numAcquisitions;
    const int kissoffs = processingControl->Value<int>("Kissoffs");
    const int numEchoes = processingControl->Value<int>("NumEchoes");
    const int imageSize = processingControl->Value<int>("TransformXRes");
    const bool isChopDataInZ = processingControl->Value<bool>("ChopZ");
    float scalingFactor = processingControl->Value<float>("FirstAcqScalingFactor");
    trace.ConsoleMsg("This is a %s scan", is3DASL ? "3DASL" : "Spiral2D");
    trace.ConsoleMsg("numAcquisitions: %d", numAcquisitions);
    trace.ConsoleMsg("numPointsPerArm: %d", numPointsPerArm);
    trace.ConsoleMsg("numArms: %d", numArms);
    trace.ConsoleMsg("numChannels: %d", numChannels);
    trace.ConsoleMsg("numSlices: %d", numSlices);
    trace.ConsoleMsg("kissoffs: %d", kissoffs);
    trace.ConsoleMsg("numEchoes: %d", numEchoes);
    trace.ConsoleMsg("ChopDataInZ is %s", isChopDataInZ ? "true" : "false");
    trace.ConsoleMsg("B0MapEnabled is %s", isChopDataInZ ? "true" : "false");
    trace.ConsoleMsg("FirstAcqScalingFactor = %f", scalingFactor);

    // Set the GERecon::Path locations prior to loading the saved files
    GERecon::Path::SetAllInputPaths(scanArchiveFullPath.parent_path() / "ScanArchiveFiles");
    scanArchive->LoadSavedFiles();

    // Load the archive storage which contains all acquisition data held in the archive
    GERecon::Acquisition::ArchiveStoragePointer archiveStorage = GERecon::Acquisition::ArchiveStorage::Create(scanArchive);
    
    // Determine how many control (DAB) packets are contained in the storage
    const size_t numControls = archiveStorage->AvailableControlCount();

    // Initialize variables
    int acquisitionCounter = 0;
    int viewValue = 0;
    int frameType = 0;
    int armIndex = 0;
    int echoIndex = 0;
    int sliceIndex = 0;
    // Assume nexCounter was same for all the passes to align with solo.
    // Actually nex may be different for PW & PD in 3DASL (improve opportunity) 
    int nexCounter = 0;
    // Slice number in last acquisition may be different from others and not equal to numSlices per acquisition
    int sliceCounter = 0;
    int geometricSliceIndex = 0;
    int frameCombineOperation = 0;
    
    const float maxPixelValue = static_cast<float>(std::numeric_limits<short>().max());
    const float minPixelValue = 0.0f;

    int numProcessedSlices = numSlices-2 * kissoffs;
    const int zStart = kissoffs;
    const int zEnd = numSlices - kissoffs - 1;
    ComplexFloatCube armStackData(numPointsPerArm, numArms, numSlices); // [Kx, Ky, Z]
    ComplexFloatCube transformedData(numPointsPerArm, numArms, numSlices); // [Kx, Ky, Kz]
    ComplexFloat5D acquiredRawData(numPointsPerArm, numArms, numChannels, numSlices, numEchoes); // Raw data loaded from ScanArchive Per acquisition
    ComplexFloat4D processedRawData(numPointsPerArm, numArms, numChannels, numProcessedSlices);  // Spiral data ready for Spiral recon Per echo per acquisition
    FloatCube imageData(imageSize, imageSize, numProcessedSlices);

    // Create an HDF5 file to save debug data to. This debug file will be enabled if the Snap Key
    // specified in this constructor is included on the command line for this rehearsal application:
    //      --GEHdf5.Snap.Key SpiralScanArchiveReherasalDebug
    const GEHdf5::Snap debugFile("SpiralScanArchiveRehearsal.h5", "SpiralScanArchiveReherasalDebug");
    std::stringstream outputNameString;

    // Create a Z-transform object to complete 1D IFFT along the Z-direction (slice)
    GERecon::Cartesian3D::ZTransformer zTransformer(*processingControl);

    // Initialize Gradwarp. The gradwarp file source will use the gw_coils.dat file stored in the ScanArchive
    // to initialize the gradwarp plugin
    const GERecon::GradwarpMode gradwarpMode = processingControl->ValueStrict<GERecon::GradwarpMode>("GradwarpMode");
    const GradwarpFileSourcePtr gradwarpFileSource = boost::make_shared<GradwarpFileSource>();
    GradwarpPlugin gradwarp(*processingControl, gradwarpMode, gradwarpFileSource);

    // Create a SpiralWorker
    const Spiral::SpiralWorker spiral;

    // Keep count of image numbers
    int imageNumber = 1;
    
    // Loop over all control packets in the archive. Some control packets are scan control packets
    // which may indicate the end of an acquisition (pass) or the end of the scan. Other control
    // packets are frame control packets which describe the raw frame (or view) data they're 
    // associated with. All control packets and associated frame data are stored in the archive
    // in the order they're acquired.
    for(size_t controlPacketIndex = 0; controlPacketIndex < numControls; ++controlPacketIndex)
    {
        const GERecon::Acquisition::FrameControlPointer controlPacketAndFrameData = archiveStorage->NextFrameControl();

        if(controlPacketAndFrameData->Control().Opcode() == Acquisition::ProgrammableOpcode)
        {
            const Acquisition::ProgrammableControlPacket framePacket = controlPacketAndFrameData->Control().Packet().As<GERecon::Acquisition::ProgrammableControlPacket>();

            // Only include ImageFrames
            viewValue = Acquisition::GetPacketValue(framePacket.viewNumH, framePacket.viewNumL);
            frameType = viewValue == 0 ? GERecon::Acquisition::BaselineFrame : GERecon:: Acquisition::ImageFrame;

            if( frameType == GERecon::Acquisition::ImageFrame )
            {
                // If packet view number == 0, then this is a baseline view, else correct the baseline view number
                armIndex = viewValue == 0 ? 0 : viewValue - 1;
                echoIndex = framePacket.echoNum;

                sliceIndex = Acquisition::GetPacketValue(framePacket.sliceNumH, framePacket.sliceNumL);
                if( (numArms-1 == armIndex) && (0 == echoIndex) && (0 == sliceIndex) && (0 == acquisitionCounter))
                {
                    ++nexCounter;
                }
                
                const ComplexFloatCube frameRawData = controlPacketAndFrameData->Data();
                frameCombineOperation = framePacket.operation;
                if(GERecon::Acquisition::StoreFrame == frameCombineOperation)
                {
                    if( (numArms-1 == armIndex) && (0 == echoIndex))
                    {
                        ++sliceCounter;
                    }
                    acquiredRawData(Range::all(), armIndex, Range::all(), sliceIndex, echoIndex) = frameRawData(Range::all(), Range::all(), 0);
                }
                else if(GERecon::Acquisition::FramePlusSource == frameCombineOperation)
                {
                    acquiredRawData(Range::all(), armIndex, Range::all(), sliceIndex, echoIndex) += frameRawData(Range::all(), Range::all(), 0);
                }
                else if(GERecon::Acquisition::FrameMinusSource == frameCombineOperation) // parasoft-suppress  METRICS-23 "Nesting required"
                {
                    acquiredRawData(Range::all(), armIndex, Range::all(), sliceIndex, echoIndex) =
                        frameRawData(Range::all(), Range::all(), 0) - acquiredRawData(Range::all(), armIndex, Range::all(), sliceIndex, echoIndex);
                }
                else if(GERecon::Acquisition::SourceMinusFrame == frameCombineOperation)
                {
                    acquiredRawData(Range::all(), armIndex, Range::all(), sliceIndex, echoIndex) -= frameRawData(Range::all(), Range::all(), 0);
                }
            }
        }
        else if(controlPacketAndFrameData->Control().Opcode() == Acquisition::ScanControlOpcode)
        {
            trace.ConsoleMsg("Received scan control packet, acquisition %d is done.", acquisitionCounter);

            outputNameString.str("");
            outputNameString << "acq" << std::setw(3) << std::setfill('0') << acquisitionCounter;
            debugFile.Write(acquiredRawData, outputNameString.str() + "/acquiredRawData");

            for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
            {
                outputNameString << "/echo" << std::setw(3) << std::setfill('0') << currentEcho;
                if( is3DASL )
                {
                    std::stringstream subNameString;
                    for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel) // parasoft-suppress  METRICS-23 "Nesting required"
                    {
                        subNameString.str("");
                        subNameString << "acq" << std::setw(3) << std::setfill('0') << acquisitionCounter;
                        subNameString << "/echo" << std::setw(3) << std::setfill('0') << currentEcho;
                        subNameString << "/channel" << std::setw(3) << std::setfill('0') << currentChannel;

                        armStackData = acquiredRawData(Range::all(), Range::all(), currentChannel, Range::all(), currentEcho);
                        debugFile.Write(armStackData, subNameString.str() + "/armStackDataRaw");

                        if(isChopDataInZ)
                        {
                            armStackData = MDArray::where(tensor::k % 2 == 1, -armStackData, armStackData);
                        }
                        debugFile.Write(armStackData, subNameString.str() + "/armStackData");

                        // Inverse FFT along Z
                        zTransformer.Apply(transformedData, armStackData);
                        debugFile.Write(transformedData, subNameString.str() + "/transformedData");

                        processedRawData(Range::all(), Range::all(), currentChannel, Range::all()) = transformedData(Range::all(), Range::all(), Range(zStart, zEnd));
                    }
                    debugFile.Write(processedRawData, outputNameString.str() + "/processedRawData");

                    if( acquisitionCounter > 0 ) // parasoft-suppress  METRICS-23 "Nesting required"
                    {
                        scalingFactor = 1.0f;
                    }
                    // Spiral reconstruction
                    SpiralUtil::SpiralRecon(imageData, processedRawData, processingControl, scanArchiveFullPath.parent_path().string(), scalingFactor/nexCounter);
                }
                else // Spiral2D
                {
                    // Spiral reconstruction
                    SpiralUtil::SpiralRecon(imageData, acquiredRawData(Range::all(), Range::all(), Range::all(), Range::all(), currentEcho), processingControl, scanArchiveFullPath.parent_path().string(), scalingFactor/nexCounter);
                    numProcessedSlices = sliceCounter;
                }

                debugFile.Write(imageData, outputNameString.str() + "/imageData");

                for(int currentSlice = 0; currentSlice < numProcessedSlices; ++currentSlice)
                {
                    outputNameString.str("");
                    outputNameString << "acq" << std::setw(3) << std::setfill('0') << acquisitionCounter;
                    outputNameString << "/echo" << std::setw(3) << std::setfill('0') << currentEcho;
                    outputNameString << "/slice" << std::setw(3) << std::setfill('0') << currentSlice;

                    FloatMatrix sliceImage = imageData(Range::all(), Range::all(), currentSlice);

                    // Get physical slice number and its orientation and corner points for current acquired slice
                    geometricSliceIndex = sliceTable.GeometricSliceNumber(acquisitionCounter, currentSlice);
                    const SliceCorners sliceCorners = sliceTable.AcquiredSliceCorners(geometricSliceIndex);

                    // Perform Gradwarp
                    gradwarp.Run(sliceImage, sliceCorners, geometricSliceIndex);
                    debugFile.Write(sliceImage, outputNameString.str() + "/imageGradwarp");

                    // Clip the image
                    Clipper::Apply(sliceImage, maxPixelValue, minPixelValue);
                    debugFile.Write(sliceImage, outputNameString.str() + "/imageClip");

                    // Rotate and transpose image
                    const SliceOrientation sliceOrientation = sliceTable.SliceOrientation(geometricSliceIndex);
                    FloatMatrix rotatedImage = RotateTranspose::Apply<float>(sliceImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());
                    debugFile.Write(rotatedImage, outputNameString.str() + "/imageRotate");

                    // Apply circular mask
                    spiral.ApplyCircularMask(rotatedImage);
                    debugFile.Write(rotatedImage, outputNameString.str() + "/imageCircular");

                    // Create storage for final image represented as a short matrix - what is sent to host.
                    ShortMatrix shortImage(rotatedImage.shape());

                    // Cast final image from float to short
                    shortImage = MDArray::cast<short>(rotatedImage);

                    // Save Image
                    std::ostringstream fileName;
                    fileName << "Image_Slice" << geometricSliceIndex << "_Echo" << currentEcho << "_Acq" << acquisitionCounter << ".dcm";
                    const ImageCorners imageCorners(sliceCorners, sliceOrientation);
                    dicomSeries.SaveImage(fileName.str(), shortImage, imageNumber, imageCorners);

                    ++imageNumber;
                } // end of slice loop
            } // end of echo loop

            ++acquisitionCounter;
            sliceCounter = 0;
            acquiredRawData = 0;
        } // opcode
    } // end of control loop
}
