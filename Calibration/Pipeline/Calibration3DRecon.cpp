// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <complex>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include <Dicom/MR/Image.h>
#include <Dicom/Instance.h>
#include <Dicom/Entity.h>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Asset/Calibration.h>

#include <Orchestra/Calibration/3D/LxControlSource.h>
#include <Orchestra/Calibration/3D/PassSetMap.h>

#include <Orchestra/Calibration/Common/RawFile.h>

#include <Orchestra/Cartesian2D/KSpaceTransformer.h>
#include <Orchestra/Cartesian3D/ZTransformer.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Common/SliceOrientation.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CalibrationRehearsalUtils.h"
#include "Calibration3DRecon.h"

using namespace MDArray;
using namespace GEDicom;
using namespace GERecon;
using namespace GERecon::Legacy;

void GERecon::Calibration::Calibration3DRecon(const std::string& pfileName, const GERecon::GradientType gradientType) // parasoft-suppress  METRICS-22 "Procedural rehearsal code handles many cases in a single int main pipeline"
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("Calibration3D");

    trace.ConsoleMsg("Running Orchestra Calibration3D Rehearsal");

    // Setup a pfile object. Also, extract the pfile directory from the full pfile path passed in to this function.
    const boost::filesystem::path pFileDirectory = boost::filesystem::path(pfileName).parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "RehearsalDicoms";
    boost::shared_ptr<Pfile> pfilePointer = Pfile::Create(pfileName);

    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(pfilePointer);

    const boost::shared_ptr<Control::ProcessingControl> processingControl = pfilePointer->CreateOrchestraProcessingControl<GERecon::Calibration3D::LxControlSource>();    
    const int numSurfaceCoilChannels = processingControl->Value<int>("NumChannels");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int acqZRes = processingControl->Value<int>("AcquiredZRes");
    const int kissoffs = processingControl->Value<int>("Kissoffs");
    const float pureScale = processingControl->Value<float>("PureScale");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");
    const Calibration3D::PassSetMap& passSetMap = processingControl->ValueStrict<Calibration3D::PassSetMap>("LxPassSetMap");
    const GERecon::SliceInfoTable& sliceTable = processingControl->Value<GERecon::SliceInfoTable>("SliceTable");

    // Create a transformer (2D) object to complete the K-Space to image-space transformation.
    // This will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformer(*processingControl);

    // Create a ztransformer object to complete 1D IFFT along the Z-direction (slice)
    Cartesian3D::ZTransformer zTransformer(*processingControl);

    // Storage for transformed image data, thus the image sizes.
    ComplexFloatMatrix imageData(imageXRes, imageYRes);
    
    // Temporary storage for conducting Z-Transform 
    ComplexFloatCube dataKxKyKz(acqXRes, acqYRes, acqZRes); 

    // Create channel combiner object that will do the channel combining work in channel loop.
    SumOfSquares channelCombiner(channelWeights);

    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, gradientType);

    // Each pass of a 3D Calibration Scan appears as a separate acquisition
    // containing unique slice locations. Thus, the slice count in the pfile
    // (which indicates the count of unique slice locations in the pfile) 
    // is a count of all slices in the scan, not all unique slices in the scan
    // To account for this, divide numSlices by the total number of acquisitions
    // (passes) in the scan.
    const int numUniqueSliceLocations = pfilePointer->SliceCount() / pfilePointer->AcqCount();
    const int totalNumPassSets = passSetMap.NumPassSets();

    // Allocate space for surface coil pass set combined volume
    // and initialize normalization scalar which is used to scale
    // the pass set combined volume based on the weights applied
    // to data from each surface coil pass. This combined surface
    // coil data is used for Asset calibration
    const int numReconstructedSlices = numUniqueSliceLocations - (2 * kissoffs);
    Array<std::complex<float>, 4> surfaceCoilPassSetCombinedVolume(imageXRes, imageYRes, numSurfaceCoilChannels, numReconstructedSlices);
    surfaceCoilPassSetCombinedVolume = 0.0f;
    float surfaceCoilNormalizationScalar = 0.0f;

    // Allocate space for pure pass set combined data. Pure calibration
    // data is combined from all pass sets that contain both a surface
    // coil acquisition and a volume coil acquisition
    float pureCalibrationNormalizationScalar = 0.0f;
    Array<std::complex<float>, 4> pureCalSurfaceCoilPassSetCombinedVolume(imageXRes, imageYRes, numSurfaceCoilChannels, numReconstructedSlices);
    Array<std::complex<float>, 3> pureCalVolumeCoilPassSetCombinedVolume(imageXRes, imageYRes, numReconstructedSlices);
    pureCalSurfaceCoilPassSetCombinedVolume = 0.0f;
    pureCalVolumeCoilPassSetCombinedVolume = 0.0f;

    // Initialize raw calibration file. Remove file if it currently exists. 
    // This file may be used by reconstruction algorithms that require access 
    // to raw calibration kSpace or imageSpace
    GERecon::Calibration::RawFile rawCalibrationFile(*processingControl, GERecon::Calibration::ThreeD);
    const boost::filesystem::path rawCalFilePath = pFileDirectory / "RawCalibration.h5";
    boost::filesystem::remove(rawCalFilePath);
    rawCalibrationFile.UseWritePaths(rawCalFilePath, rawCalFilePath);

    // Write meta-data to the raw calibration file
    rawCalibrationFile.WriteHeader();

    // Initialize matrix used to store a volume of calibration kSpace 
    // that is to be written to the raw calibration file after the slice
    // loop is complete.
    CalibrationData calibrationKSpace;

    // Keep a flat count of the current pass index. This pass index will be
    // multiplied by the number of slices in a pass to determine the slice
    // index in the pfile.
    int flatPassIndex = 0;

    for(int passSetIndex = 0; passSetIndex < totalNumPassSets; ++passSetIndex)
    {
        const int numPassesInPassSet = passSetMap.NumPassesInPassSet(passSetIndex);

        // Pfile data has been accumulated for each nex, but it has not been divided
        // by the total nex count. The total nex count can be different for each 
        // pass in a 3D Cal acquisition. Thus, determine the nex count for the 
        // current pass set here and use this value to scale the raw kSpace data
        // for all passes in this pass set.
        const int nexCount = passSetMap.PassSetNexCount(passSetIndex);

        // All of the currently supported 3D Calibration paradigms
        // contain a surface pass in all pass sets that are acquired.
        surfaceCoilNormalizationScalar += passSetMap.PassSetWeight(passSetIndex);
        
        // Pure normalization scalar Pass set combined data from pass sets
        // with both a surface and volume pass is used for pure calibration.
        // Check if this pass set contains both a surface and a volume pass
        // and accumulate normalization scalar.
        if(passSetMap.PassSetType(passSetIndex) == GERecon::Calibration::SurfaceAndVolume)
        {
            pureCalibrationNormalizationScalar += passSetMap.PassSetWeight(passSetIndex);
        }

        for(int passIndex = 0; passIndex < numPassesInPassSet; ++passIndex)
        {
            trace.ConsoleMsg("Processing Pass Set: %d, Pass %d of %d", passSetIndex, passIndex, numPassesInPassSet);
            const int echo = 0; // 3D Calibration scans do not contain more than one echo

            // Load K-space data
            const int numChannels = passSetMap.IsVolumePass(flatPassIndex) ? 1 : numSurfaceCoilChannels; 
            ComplexFloat4D dataKxKyZChan(acqXRes, acqYRes, numUniqueSliceLocations, numChannels); 
            dataKxKyZChan = 0; 
            if(pfilePointer->IsZEncoded()) // Apply 1D-FFT along the Z-direction (slice) if needed
            {
                for(int channel = 0; channel < numChannels; ++channel)
                {
                    for(int kz = 0; kz < acqZRes; ++kz)
                    {
                        dataKxKyKz(Range::all(), Range::all(), kz) = pfilePointer->KSpaceData<float>
                            (Legacy::Pfile::PassSlicePair(flatPassIndex, kz), echo, channel);
                    }
                    // Apply 1D-IFFT along the Z-direction (slice)
                    ComplexFloatCube data(dataKxKyZChan(Range::all(), Range::all(), Range::all(), channel));
                    zTransformer.Apply(data, dataKxKyKz); 
                }   
            }
            else
            {
                for(int channel = 0; channel < numChannels; ++channel)
                {
                    for(int slice = 0; slice < numUniqueSliceLocations; ++slice)
                    {
                        dataKxKyZChan(Range::all(), Range::all(), slice, channel) = pfilePointer->KSpaceData<float>(Legacy::Pfile::PassSlicePair(flatPassIndex, slice), echo, channel);
                    }
                }
            }

            // Account for the kissoff slices on either end of the acquired slab
            // Pfile data for 3D Calibration scans already has the z-transform
            // applied to it. Thus, the kissoff slice data is no longer needed
            // Simply ignore the kissoff slices here.
            for(int slice = kissoffs; slice < (numUniqueSliceLocations - kissoffs); ++slice)
            {
                channelCombiner.Reset();

                // This loop is looping over the slice, or 'z', index of the dataKyKyZChan array that
                // was populated above. Each z index in this matrix corresponds to an acquired slice
                // index. Since kissoff slices are discarded, adjust for kissoff slices in the acquired
                // slice number here.
                // Also, lookup the geometric slice number corresponding to the current acquired slice
                // number. Since the slice table does not include kissoff slices, the kissoff-adjusted
                // acquired slice number is used as an input to the slice table API. 
                const int currentAcquiredSliceNumber = slice-kissoffs;
                const int currentGeometricSliceNumber = sliceTable.GeometricSliceNumber(currentAcquiredSliceNumber);

                ComplexFloatMatrix combinedImage;
                if(passSetMap.IsVolumePass(flatPassIndex))
                {
                    // No need for channel combine, this is single channel data
                    // Transform and set combinedData matrix to this single channel's
                    // complex image space
                    
                    const int channel = 0;
                    ComplexFloatMatrix kSpace = dataKxKyZChan(Range::all(), Range::all(), slice, channel); 

                    // Nex Scaling
                    kSpace /= static_cast<float>(nexCount);

                    // Ensure calibrationKSpace matrix is allocated prior to copying kSpace data into it.
                    const int numChannelsInThisPass = 1;
                    MDArray::AssureShape(calibrationKSpace, TinyVector<int,4>(kSpace.extent(firstDim), kSpace.extent(secondDim), numChannelsInThisPass, numReconstructedSlices));
                    calibrationKSpace(Range::all(), Range::all(), channel, currentAcquiredSliceNumber) = kSpace;

                    transformer.Apply(imageData, kSpace);

                    // Scale by rdb_hdr_pure_scale
                    imageData *= pureScale;

                    pureCalVolumeCoilPassSetCombinedVolume(Range::all(), Range::all(), currentGeometricSliceNumber) += (imageData * passSetMap.PassSetWeight(passSetIndex));

                    combinedImage.reference(imageData);
                }
                else
                {
                    for(int channel = 0; channel < numSurfaceCoilChannels; ++channel)
                    {
                        ComplexFloatMatrix kSpace = dataKxKyZChan(Range::all(), Range::all(), slice, channel);
                        
                        // Nex Scaling
                        kSpace /= static_cast<float>(nexCount);

                        // Ensure calibrationKSpace matrix is allocated prior to copying kSpace data into it.
                        MDArray::AssureShape(calibrationKSpace, TinyVector<int,4>(kSpace.extent(firstDim), kSpace.extent(secondDim), numSurfaceCoilChannels, numReconstructedSlices));
                        calibrationKSpace(Range::all(), Range::all(), channel, currentAcquiredSliceNumber) = kSpace;

                        transformer.Apply(imageData, kSpace);

                        channelCombiner.Accumulate(imageData, channel);

                        surfaceCoilPassSetCombinedVolume(Range::all(), Range::all(), channel, currentGeometricSliceNumber) += (imageData * passSetMap.PassSetWeight(passSetIndex));

                        // If the current surface pass is from a pass set that contains both a volume and
                        // a surface pass then include this data in the pure calibration pass set combined
                        // surface volume
                        if(passSetMap.PassSetType(passSetIndex) == GERecon::Calibration::SurfaceAndVolume) // parasoft-suppress  METRICS-23 "Procedural Rehearsal Code"
                        {
                            pureCalSurfaceCoilPassSetCombinedVolume(Range::all(), Range::all(), channel, currentGeometricSliceNumber) += (imageData * passSetMap.PassSetWeight(passSetIndex));
                        }
                    }

                    combinedImage.reference(channelCombiner.GetCombinedImage());
                }

                // Create storage for final image represented as a float matrix
                FloatMatrix magnitudeImage(combinedImage.shape());

                // Convert complex data to specified image type.
                MDArray::ComplexToReal(magnitudeImage, combinedImage, MDArray::MagnitudeData);

                // Get information for current slice
                const SliceOrientation& sliceOrientation = sliceTable.SliceOrientation(currentGeometricSliceNumber);
                const SliceCorners& sliceCorners = sliceTable.SliceCorners(currentGeometricSliceNumber);
                const SliceCorners& acquiredCorners = sliceTable.AcquiredSliceCorners(currentGeometricSliceNumber);

                // Perform Gradwarp
                gradwarp.Run(magnitudeImage, slice, sliceCorners, acquiredCorners);

                // Rotate/Transpose image and corner points accordingly
                FloatMatrix rotatedImage = RotateTranspose::Apply<float>(magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());

                // Clip the image
                Clipper::Apply(rotatedImage, GERecon::MagnitudeImage);

                // Create storage for final image represented as a short matrix - what is sent to host.
                ShortMatrix finalImage(rotatedImage.shape());

                // Cast final image from float to short
                finalImage = MDArray::cast<short>(rotatedImage);

                const unsigned int imageNumber = currentGeometricSliceNumber + flatPassIndex * numReconstructedSlices;
                std::stringstream dicomFileName;
                dicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";

                const ImageCorners imageCorners(sliceCorners, sliceOrientation);
                const GEDicom::MR::ImagePointer dicomImage = dicomSeries.NewImage(finalImage, imageNumber, imageCorners);
                dicomImage->Save(dicomFileName.str());
            }

            // Save the volumes of kSpace and image space data to the raw calibration file for the pass that was just completed
            const GERecon::Calibration::ReceiverDataType kSpaceType = (passSetMap.IsVolumePass(flatPassIndex) ? VolumeKSpace : SurfaceKSpace);
            rawCalibrationFile.AddData(calibrationKSpace, kSpaceType, passSetIndex);

            ++flatPassIndex;
        }
    }

    // Data from multiple pass sets is combined in the loop above. Each pass set can be weighted
    // differently. Normalize the pass set combined data here to account for the varying pass set weights
    surfaceCoilPassSetCombinedVolume /= surfaceCoilNormalizationScalar;
    pureCalSurfaceCoilPassSetCombinedVolume /= pureCalibrationNormalizationScalar;
    pureCalVolumeCoilPassSetCombinedVolume /= pureCalibrationNormalizationScalar;

    // Unchop and save image space data to the raw calibration file
    CalibrationData unchoppedCalData(surfaceCoilPassSetCombinedVolume.shape()); 
    unchoppedCalData = UnchopCalImageSpaceData(surfaceCoilPassSetCombinedVolume); 
    rawCalibrationFile.AddData(unchoppedCalData, SurfaceImageSpaceAllPassSets);

    unchoppedCalData.resize(pureCalSurfaceCoilPassSetCombinedVolume.shape());
    unchoppedCalData = UnchopCalImageSpaceData(pureCalSurfaceCoilPassSetCombinedVolume); 
    rawCalibrationFile.AddData(unchoppedCalData, SurfaceImageSpacePairedPassSets);

    // CalibrationData is saved as a 4D array. The pureCalVolumeCoilPassSetCombinedVolume array
    // is a 3D array. Create a 4D array with channel dimension size = 1 to save this data
    // to the raw calibration file
    const int numVolumeCoilChannels = 1; // Only 1 volume coil channel is currently supported
    CalibrationData volumeImageSpace(pureCalVolumeCoilPassSetCombinedVolume.extent(firstDim), pureCalVolumeCoilPassSetCombinedVolume.extent(secondDim), 
                                           numVolumeCoilChannels, pureCalVolumeCoilPassSetCombinedVolume.extent(thirdDim));
    volumeImageSpace(Range::all(), Range::all(), 0, Range::all()) = pureCalVolumeCoilPassSetCombinedVolume;
    unchoppedCalData.resize(volumeImageSpace.shape());
    unchoppedCalData = UnchopCalImageSpaceData(volumeImageSpace); 
    rawCalibrationFile.AddData(unchoppedCalData, VolumeImageSpace);

    // Once all pass sets have been combined, process the pass set combined surface coil and
    // volume coil data (if applicable) to generate Asset, Pure, and Mcsi Calibration Files
    if(processingControl->Value<bool>("AssetCalibration"))
    {
        trace.ConsoleMsg("Running Asset Calibration Processing");

        Array<std::complex<float>, 4> surfaceCoilPassSetCombinedVolumeAlternateApodization(surfaceCoilPassSetCombinedVolume.shape());
        surfaceCoilPassSetCombinedVolumeAlternateApodization = surfaceCoilPassSetCombinedVolume; 

        ComplexFloatCube alternateApodizedImages;
        const float originalFermiRadius = processingControl->Value<float>("FermiRadius");
        const float originalFermiWidth = processingControl->Value<float>("FermiWidth");
        for(int slice = 0; slice < surfaceCoilPassSetCombinedVolume.extent(fourthDim); ++slice)
        {
            ComplexFloatCube dataToReapodize = surfaceCoilPassSetCombinedVolumeAlternateApodization(Range::all(), Range::all(), Range::all(), slice);
            GERecon::Asset::Calibration::GenerateAlternatelyApodizedCalData(alternateApodizedImages, dataToReapodize, originalFermiRadius, originalFermiWidth);
            dataToReapodize = alternateApodizedImages;
        }

        ProcessAssetCalibration(surfaceCoilPassSetCombinedVolume, surfaceCoilPassSetCombinedVolumeAlternateApodization, *processingControl, pFileDirectory);
    }

    if(processingControl->Value<bool>("PureCalibration"))
    {
        trace.ConsoleMsg("Running Pure Calibration Processing");

        const ComplexFloatCube channelCombinedComplexSurfaceCoilDataForPure(pureCalSurfaceCoilPassSetCombinedVolume.extent(firstDim), pureCalSurfaceCoilPassSetCombinedVolume.extent(secondDim), pureCalSurfaceCoilPassSetCombinedVolume.extent(fourthDim));
        for (int geometricSliceNumber = 0; geometricSliceNumber < channelCombinedComplexSurfaceCoilDataForPure.extent(thirdDim); ++geometricSliceNumber)
        {
            // Combine Channels for the surface Coil for this slice
            channelCombiner.Reset();
            for(int channel = 0; channel < numSurfaceCoilChannels; ++channel)
            {
                const ComplexFloatMatrix channelImage = pureCalSurfaceCoilPassSetCombinedVolume(Range::all(), Range::all(), channel, geometricSliceNumber);
                channelCombiner.Accumulate(channelImage, channel);
            }            
            const ComplexFloatMatrix channelCombinedImage = channelCombiner.GetCombinedImage();
            channelCombinedComplexSurfaceCoilDataForPure(Range::all(), Range::all(), geometricSliceNumber) = channelCombinedImage;
        }

        ProcessPureCalibration(channelCombinedComplexSurfaceCoilDataForPure, pureCalVolumeCoilPassSetCombinedVolume, *processingControl, gradwarp, pFileDirectory);
    }
}
