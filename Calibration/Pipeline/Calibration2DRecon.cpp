// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include <Dicom/MR/Image.h>
#include <Dicom/Instance.h>
#include <Dicom/Entity.h>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Asset/Calibration.h>

#include <Orchestra/Calibration/Common/RawFile.h>

#include <Orchestra/Cartesian2D/KSpaceTransformer.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Common/SliceOrientation.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/LxControlSource.h>
#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CalibrationRehearsalUtils.h"
#include "Calibration2DRecon.h"

using namespace MDArray;
using namespace GEDicom;
using namespace GERecon;
using namespace GERecon::Legacy;

void GERecon::Calibration::Calibration2DRecon(const std::string& pfileName, const GERecon::GradientType gradientType) // parasoft-suppress  METRICS-22 "Procedural rehearsal code handles many cases in a single int main pipeline"
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("Calibration2D");

    trace.ConsoleMsg("Running Orchestra Calibration2D Rehearsal");

    // Setup a pfile object. Also, extract the pfile directory from the full pfile path passed in to this function.
    const boost::filesystem::path pFileDirectory = boost::filesystem::path(pfileName).parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "RehearsalDicoms";
    boost::shared_ptr<Pfile> pfilePointer = Pfile::Create(pfileName);

    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(pfilePointer);

    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfilePointer->CreateOrchestraProcessingControl<GERecon::Legacy::LxControlSource>();    
    const int numSurfaceCoilChannels = processingControl->Value<int>("NumChannels");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const float pureScale = processingControl->Value<float>("PureScale");
    const float fermiWidth = processingControl->Value<float>("FermiWidth");
    const float alternateRadius = 9.0f; // Alternate fermi radius is always 9.0f
    const float fermiEccentricity = processingControl->Value<float>("FermiEccentricity");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");
    const GERecon::SliceInfoTable& sliceTable = processingControl->Value<SliceInfoTable>("SliceTable");

    // Create two transformer objects. One with the standard apodization filter settings
    // and one with alternate apodization filter settings.
    // These will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformerStandardApodization(*processingControl);
    Cartesian2D::KSpaceTransformer transformerAlternateApodization(*processingControl, fermiWidth, alternateRadius, fermiEccentricity, true);

    // Storage for transformed image data, thus the image sizes.
    ComplexFloatMatrix imageData(imageXRes, imageYRes);

    // Create channel combiner object that will do the channel combining work in channel loop.
    SumOfSquares channelCombiner(channelWeights);

    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, gradientType);

    const int numSlices = pfilePointer->SliceCount();

    const MDArray::Array<std::complex<float>, 4> surfaceCoilVolumeStandardApodization(imageXRes, imageYRes, numSurfaceCoilChannels, numSlices);
    const MDArray::Array<std::complex<float>, 4> surfaceCoilVolumeAlternateApodization(imageXRes, imageYRes, numSurfaceCoilChannels, numSlices);

    const ComplexFloatCube surfaceCoilChannelCombinedVolume(imageXRes, imageYRes, numSlices);
    const ComplexFloatCube referenceCoilVolumeForPureCal(imageXRes, imageYRes, numSlices);

    // There will always be one surface coil phase. If this calibration
    // was done with a PURE compatible coil then a second volume (body) coil
    // phase will be acquired
    const int numPhases = pfilePointer->PhaseCount();

    // Initialize raw calibration file. Remove file if it currently exists. 
    // This file may be used by reconstruction algorithms that require access 
    // to raw calibration kSpace or imageSpace
    GERecon::Calibration::RawFile rawCalibrationFile(*processingControl, GERecon::Calibration::TwoD);
    const boost::filesystem::path rawCalFilePath = pFileDirectory / "RawCalibration.h5";
    boost::filesystem::remove(rawCalFilePath);
    rawCalibrationFile.UseWritePaths(rawCalFilePath, rawCalFilePath);

    // Write meta-data to the raw calibration file
    rawCalibrationFile.WriteHeader();

    // Initialize matrix used to store a volume of calibration kSpace 
    // that is to be written to the raw calibration file after the slice
    // loop is complete.
    CalibrationData calibrationKSpace;

    // Initialize image number to zero and increment with each image that is generated
    int imageNumber = 0;

    for(int phase = 0; phase < numPhases; ++phase)
    {
        const int numChannelsInThisPass = (phase > 0) ? 1 : numSurfaceCoilChannels;

        for(int slice = 0; slice < numSlices; ++slice)
        {
            channelCombiner.Reset();

            // Raw calibration kSpace is stored in acquisition order. The slice loop is looping
            // over geometric slice locations. Thus, the acquired slice number corresponding to
            // the current geometric slice is needed below when data is copied into the 
            // calibration kSpace array. The calibration kSpace array is written to the raw 
            // calibration file as a single volume after the slice loop is complete.
            const int currentAcquiredSliceNumber = sliceTable.AcquiredSliceInfo(slice).second;

            ComplexFloatMatrix combinedImage;
            if(numChannelsInThisPass == 1)
            {
                // This is a volume (body) coil pass
                // No need for channel combine, this is single channel data
                // Transform and set combinedData matrix to this single channel's
                // complex image space. The single channel complex image space will
                // be used to create dicom images. It will also be passed to the PURE
                // calibration function which will convert the images to magnitude images
                // and run calibration processing.
                const int echo = 0;
                const int channel = 0;
                const ComplexFloatMatrix kSpace = pfilePointer->KSpaceData<float>(slice, echo, channel, phase);

                // Ensure calibrationKSpace matrix is allocated prior to copying kSpace data into it.
                MDArray::AssureShape(calibrationKSpace, TinyVector<int,4>(kSpace.extent(firstDim), kSpace.extent(secondDim), numChannelsInThisPass, numSlices));
                calibrationKSpace(Range::all(), Range::all(), channel, currentAcquiredSliceNumber) = kSpace;

                transformerStandardApodization.Apply(imageData, kSpace);

                // Scale by rdb_hdr_pure_scale
                imageData *= pureScale;
                referenceCoilVolumeForPureCal(Range::all(), Range::all(), slice) = imageData;
                combinedImage.reference(imageData);
            }
            else
            {
                // This is a surface coil pass.
                // Create two sets of individual channel images. One set with the standard
                // apodization filter settings from the raw header and one set with an 
                // alternate set of apodizationfilter settings. Hold on to these two sets
                // of images so they can be passed to ASSET calibration processing.
                // Also, run a sum of squares channel combine to produce complex channel
                // combined data that will be used to generate dicoms for this scan. The
                // complex channel combined data will also be passed to the PURE calibration
                // routine (if applicable). The PURE calibration function will convert these
                // to magnitude images and run PURE calibration processing.
                // Finally, hold on to the calibration kSpace and image space such that it
                // can be saved to the raw calibration file as a single volume after the
                // slice loop.
                channelCombiner.Reset();
                for(int channel = 0; channel < numSurfaceCoilChannels; ++channel)
                {
                    // 2D Calibration scans do not contain more than one echo
                    const int echo = 0;
                    const ComplexFloatMatrix kSpace = pfilePointer->KSpaceData<float>(slice, echo, channel, phase);

                    // Ensure calibrationKSpace matrix is allocated prior to copying kSpace data into it.
                    MDArray::AssureShape(calibrationKSpace, TinyVector<int,4>(kSpace.extent(firstDim), kSpace.extent(secondDim), numSurfaceCoilChannels, numSlices));
                    calibrationKSpace(Range::all(), Range::all(), channel, currentAcquiredSliceNumber) = kSpace;

                    transformerStandardApodization.Apply(imageData, kSpace);
                    surfaceCoilVolumeStandardApodization(Range::all(), Range::all(), channel, slice) = imageData;
                    channelCombiner.Accumulate(imageData, channel);

                    ComplexFloatMatrix alternateApodizedImage = surfaceCoilVolumeAlternateApodization(Range::all(), Range::all(), channel, slice);
                    transformerAlternateApodization.Apply(alternateApodizedImage, kSpace);
                }
                
                const ComplexFloatMatrix currentCombinedImage = channelCombiner.GetCombinedImage();
                surfaceCoilChannelCombinedVolume(Range::all(), Range::all(), slice) = currentCombinedImage;
                combinedImage.reference(currentCombinedImage);
            }
            
            // Create storage for final image represented as a float matrix
            FloatMatrix magnitudeImage(combinedImage.shape());

            // Convert complex data to specified image type.
            MDArray::ComplexToReal(magnitudeImage, combinedImage, MDArray::MagnitudeData);

            // Get information for current slice
            const SliceOrientation& sliceOrientation = pfilePointer->Orientation(slice);
            const SliceCorners& sliceCorners = pfilePointer->Corners(slice);
            const SliceCorners& acquiredCorners = pfilePointer->AcquiredCorners(slice);

            // Perform Gradwarp
            gradwarp.Run(magnitudeImage, slice, sliceCorners, acquiredCorners);

            // Rotate/Transpose image accordingly
            FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

            // Clip the image
            Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

            // Create storage for final image represented as a short matrix - what is sent to host.
            ShortMatrix finalImage(rotatedImage.shape());

            // Cast final image from float to short
            finalImage = MDArray::cast<short>(rotatedImage);

            std::stringstream dicomFileName;
            dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";

            const ImageCorners imageCorners(sliceCorners, sliceOrientation);
            const GEDicom::MR::ImagePointer dicomImage = dicomSeries.NewImage(finalImage, imageNumber, imageCorners);
            dicomImage->Save(dicomFileName.str());

            ++imageNumber;
        }

        // Save the volumes of kSpace and image space data to the raw calibration file for the pass that was just completed
        const int currentPassSetIndex = 0; // 2D cal always has only 1 pass-set
        const GERecon::Calibration::ReceiverDataType kSpaceType = (numChannelsInThisPass == 1 ? VolumeKSpace : SurfaceKSpace);
        rawCalibrationFile.AddData(calibrationKSpace, kSpaceType, currentPassSetIndex);

        const GERecon::Calibration::ReceiverDataType imageSpaceType = (numChannelsInThisPass == 1 ? VolumeImageSpace : SurfaceImageSpaceAllPassSets);
        if(imageSpaceType == VolumeImageSpace)
        {
            // Volume image space was copied into the array referenceCoilVolumeForPureCal. The data
            // was also scaled by pureScale. Copy this data into a CalibrationData array and divide
            // out the pureScale scalar that was applied. This allows for raw complex image space
            // to be saved to the raw calibration file.
            const CalibrationData volumeImageSpace(referenceCoilVolumeForPureCal.extent(firstDim), referenceCoilVolumeForPureCal.extent(secondDim), numChannelsInThisPass, numSlices);
            volumeImageSpace(Range::all(), Range::all(), 0, Range::all()) = referenceCoilVolumeForPureCal / pureScale;
            rawCalibrationFile.AddData(volumeImageSpace, imageSpaceType, currentPassSetIndex);
        }
        else
        {
            // Surface image space was copied to the array surfaceCoilVolumeStandardApodization in the loop above
            // No additional scaling was applied above, thus, this array can be written directly to the
            // raw calibration file
            rawCalibrationFile.AddData(surfaceCoilVolumeStandardApodization, imageSpaceType, currentPassSetIndex);
        }
    }

    if(processingControl->Value<bool>("AssetCalibration"))
    {
        ProcessAssetCalibration(surfaceCoilVolumeStandardApodization, surfaceCoilVolumeAlternateApodization, *processingControl, pFileDirectory);
    }

    if(processingControl->Value<bool>("PureCalibration"))
    {
        const float pureScale = processingControl->Value<float>("PureScale");
        ProcessPureCalibration(surfaceCoilChannelCombinedVolume, referenceCoilVolumeForPureCal, *processingControl, gradwarp, pFileDirectory, pureScale);
    }
}
