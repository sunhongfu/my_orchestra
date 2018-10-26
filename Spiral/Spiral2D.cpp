// Copyright 2018 General Electric Company.  All rights reserved.

#include <Dicom/MR/Image.h>

#include <MDArray/Utils.h>

#include <Hdf5/Snap.h>

#include <Orchestra/Control/ProcessingControl.h>

// Algorithms
#include <Orchestra/Common/GradwarpMode.h>
#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/Windows.h>

#include <Orchestra/Gradwarp/GradwarpFileSource.h>
#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/Pfile.h>
#include <Orchestra/Legacy/SliceEntry.h>

#include <Orchestra/Spiral/SpiralWorker.h>

#include "Rehearsal.h"
#include "SpiralUtil.h"

using namespace GERecon;
using namespace GERecon::Spiral;
using namespace GERecon::Legacy;
using namespace MDArray;

void GERecon::Spiral2D(const Control::ProcessingControlPointer& processingControl, Legacy::Pfile& pfile) // parasoft-suppress  METRICS-22 "Procedural Rehearsal Code"
{
    const std::string dataDirectory = boost::filesystem::path(pfile.File()).parent_path().string();

    // Pull legacy header from Pfile object, as well as Series and Exam info.
    const Legacy::ExamDataTypeStruct& examData = boost::const_pointer_cast<const LxDownloadData>(pfile.DownloadData())->ExamData();
    const Legacy::SeriesDataTypeStruct& seriesData = boost::const_pointer_cast<const LxDownloadData>(pfile.DownloadData())->SeriesData();

    // Create an HDF5 file to save debug data to. This debug file will be enabled if the Snap Key
    // specified in this constructor is included on the command line for this rehearsal application:
    //      --GEHdf5.Snap.Key Spiral2DReherasalDebug
    const GEHdf5::Snap debugFile("Spiral2DRehearsal.h5", "Spiral2DReherasalDebug");
    std::stringstream outputNameString;

    // Acquisition Parameters.
    const int numSlices = processingControl->Value<int>("NumSlices");
    const int numEchoes = processingControl->Value<int>("NumEchoes");
    const int numChannels = processingControl->Value<int>("NumChannels");
    const int acqXRes = processingControl->Value<int>("NumPointsPerArm");
    const int acqYRes = processingControl->Value<int>("NumArms");
    const int imageSize = processingControl->Value<int>("TransformXRes");
    const float scalingFactor = processingControl->Value<float>("FirstAcqScalingFactor"); 
    const SliceInfoTable& sliceTable = processingControl->ValueStrict<SliceInfoTable>("SliceTable");

    // Image gradwarp plugin
    const GERecon::GradwarpMode gradwarpMode = processingControl->ValueStrict<GERecon::GradwarpMode>("GradwarpMode");
    GradwarpPlugin gradwarp(*processingControl, gradwarpMode, GERecon::XRMBGradient);

    // Create a SpiralWorker
    const Spiral::SpiralWorker spiral;

    // Keep count of image numbers
    int imageNumber = 1;

    for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
    {
        ComplexFloat4D rawData(acqXRes, acqYRes, numChannels, numSlices);
        FloatCube imageData(imageSize, imageSize, numSlices);
        for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
        {
            for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
            {
                // Extract an echo for a single channel of data from the pfile
                rawData(Range::all(), Range::all(), currentChannel, currentSlice) = pfile.KSpaceData<float>(currentSlice, currentEcho, currentChannel);
            }
        }
        
        // Spiral reconstruction
        SpiralUtil::SpiralRecon(imageData, rawData, processingControl, dataDirectory, scalingFactor);
        debugFile.Write(imageData, "imageData");

        for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
        {
            FloatMatrix sliceImage = imageData(Range::all(), Range::all(), currentSlice);
            outputNameString.str("");
            outputNameString << "echo" << std::setw(3) << std::setfill('0') << currentEcho;
            outputNameString << "_slice" << std::setw(3) << std::setfill('0') << currentSlice;

            // Get physical slice number and its orientation and corner points for current acquired slice
            const int sliceNumber = sliceTable.GeometricSliceNumber(currentSlice);
            const SliceCorners& sliceCorners = sliceTable.SliceCorners(sliceNumber);
            const SliceOrientation& sliceOrientation = sliceTable.SliceOrientation(sliceNumber);

            debugFile.Write(sliceImage, outputNameString.str() + "/imageBeforeGradwarp");

            // Perform Gradwarp
            gradwarp.Run(sliceImage, sliceCorners, sliceNumber);
            debugFile.Write(sliceImage, outputNameString.str() + "/imageAfterGradwarp");
            
            // Clip the image
            const float maxPixelValue = static_cast<float>(std::numeric_limits<short>().max()); // Max of short 32767
            const float minPixelValue = 0.0f;
            Clipper::Apply(sliceImage, maxPixelValue, minPixelValue);
            debugFile.Write(sliceImage, outputNameString.str() + "/imageClip");

            // Rotate and transpose image
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
            fileName << "Image_Slice" << currentSlice << "_Echo" << currentEcho << ".dcm";
            
            GEDicom::MR::Image dcmImage(shortImage, imageNumber, seriesData.se_no, examData.ex_no, "MR Orchestra", "Spiral 2D");
            dcmImage.Save(fileName.str());

            ++imageNumber;
        }
    }
}
