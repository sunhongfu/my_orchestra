// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <MDArray/Utils.h>

#include <Orchestra/Asset/Calibration.h>
#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/CalibrationInfo.h>

#include <Orchestra/Calibration/Common/Types.h>

#include <Orchestra/Common/SliceInfoTable.h>

#include <Orchestra/Core/ChannelCombiner.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Legacy/Pfile.h>

#include <Orchestra/Pure/Pure1/PureCalibrationFile.h>
#include <Orchestra/Pure/Pure1/PureCalibrationInfo.h>
#include <Orchestra/Pure/Pure1/PurePlugin.h>

#include <Orchestra/Spectro/MCSI/CalibrationFile.h>

#include "CalibrationRehearsalUtils.h"

using namespace GERecon;
using namespace MDArray;


void GERecon::Calibration::ProcessAssetCalibration(const MDArray::Array<std::complex<float>, 4>& surfaceCoilVolume, const MDArray::Array<std::complex<float>, 4>& alternatelyApodizedVolume, 
                                                   const GERecon::Control::ProcessingControl& processingControl, const boost::filesystem::path& outputDirectory)
{
    // Run asset calibration processing. As part of asset calibration
    // processing, an Mcsi calibration file is also generated
    // Asset calibration processing is run on two different data sets,
    // each with different 2D apodization parameters. First process
    // the datasets reconstructed above using apodization parameters
    // from processing control (apodization is applied inside of the
    // 2D KSpaceTransformer). Then, re-apodize the data reconstructed
    // above and run asset calibration one more time on this re-apodized
    // data.
    const int numAssetCalApodizationTypes = 2;

    FloatCube virtualCoilVolume(surfaceCoilVolume.extent(firstDim), surfaceCoilVolume.extent(secondDim), surfaceCoilVolume.extent(fourthDim));
    Array<std::complex<float>, 4> sensitivityRatiosVolume(surfaceCoilVolume.shape());
    Array<float, 4> magnitudeChannelDataVolume(surfaceCoilVolume.shape());
    Asset::CalibrationInfo calInfo = AssetCalibrationInfo(processingControl);

    Asset::CalibrationFile assetCalFile(processingControl.Value<unsigned int>("ExamNumber"), processingControl.Value<unsigned int>("CoilConfigUID"), processingControl.ValueStrict<int64_t>("CalibrationUID"));
    const boost::filesystem::path assetCalFilePath = outputDirectory / "AssetCalibration.h5";
    assetCalFile.UseWritePaths(assetCalFilePath, assetCalFilePath);
    assetCalFile.WriteHeader(calInfo);
        
    GERecon::Spectro::MCSI::CalibrationFile mcsiCalFile(processingControl.Value<unsigned int>("ExamNumber"), processingControl.Value<unsigned int>("CoilConfigUID"), processingControl.ValueStrict<int64_t>("CalibrationUID"));
    const boost::filesystem::path mcsiCalFilePath = outputDirectory / "McsiCalibration.h5";
    mcsiCalFile.UseWritePaths(mcsiCalFilePath, mcsiCalFilePath);
    mcsiCalFile.WriteHeader(processingControl);

    const GERecon::SliceInfoTable& sliceTable = processingControl.ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    for(int apodizationType = 0; apodizationType < numAssetCalApodizationTypes; ++apodizationType)
    {
        const GERecon::Asset::Calibration::Type currentType = (apodizationType == 0 ? GERecon::Asset::Calibration::Regular : GERecon::Asset::Calibration::AlternateApodization);

        Array<std::complex<float>, 4> surfaceCoilVolumeForCurrentCalType = surfaceCoilVolume;
        if(currentType == GERecon::Asset::Calibration::AlternateApodization)
        {
            // In case of alternate apodization, use the alternately apodized surface coil volume
            surfaceCoilVolumeForCurrentCalType.reference(alternatelyApodizedVolume);
        }

        for(int slice = 0; slice < surfaceCoilVolumeForCurrentCalType.extent(fourthDim); ++slice)
        {
            // Input Data
            const ComplexFloatCube passSetCombinedDataCurrentSlice = surfaceCoilVolumeForCurrentCalType(Range::all(), Range::all(), Range::all(), slice);

            // Output Data
            FloatMatrix virtualCoil = virtualCoilVolume(Range::all(), Range::all(), slice);
            ComplexFloatCube sensitivityRatios = sensitivityRatiosVolume(Range::all(), Range::all(), Range::all(), slice);
            FloatCube magnitudeChannelData = magnitudeChannelDataVolume(Range::all(), Range::all(), Range::all(), slice);
                
            const SliceOrientation& orientation = sliceTable.SliceOrientation(slice);
            AssetAndMcsiCalSliceProcessing(virtualCoil, sensitivityRatios, magnitudeChannelData, passSetCombinedDataCurrentSlice, orientation, slice, processingControl);
        }

        // Normalize calibration data and calculate mean, max, and ratio mean.
        float calMean = 0.0f;
        float calMax = 0.0f;
        float ratioMean = 0.0f;

        Asset::Calibration::Normalize(calMean, calMax, ratioMean, virtualCoilVolume, sensitivityRatiosVolume);

        // Write information specific to calibration type
        calInfo.CalibrationMean(calMean);
        calInfo.CalibrationMax(calMax);
        calInfo.SensitivityRatiosMean(ratioMean);

        assetCalFile.Write(calInfo, sensitivityRatiosVolume, virtualCoilVolume, currentType);

        // MCSI Cal Uses Alternate Apodization
        if(currentType == GERecon::Asset::Calibration::AlternateApodization)
        {
            mcsiCalFile.Write(calInfo.FirstSliceCornerPoints(), calInfo.LastSliceCornerPoints(), calInfo.ScanCenter(), magnitudeChannelDataVolume);
        }
    }
}

Asset::CalibrationInfo GERecon::Calibration::AssetCalibrationInfo(const GERecon::Control::ProcessingControl& processingControl)
{
    const GERecon::SliceInfoTable& sliceTable = processingControl.ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    Asset::CalibrationInfo calibrationInfo;
    calibrationInfo.Exam(processingControl.Value<unsigned int>("ExamNumber"));
    calibrationInfo.CoilID(processingControl.Value<unsigned int>("CoilConfigUID"));
    calibrationInfo.Channels(processingControl.Value<int>("NumChannels"));
    calibrationInfo.XRes(processingControl.Value<int>("TransformXRes"));
    calibrationInfo.YRes(processingControl.Value<int>("TransformYRes"));
    calibrationInfo.Landmark(processingControl.Value<float>("Landmark"));
    calibrationInfo.PatientEntry(processingControl.Value<int>("PatientEntry"));
    calibrationInfo.PatientPosition(processingControl.Value<int>("PatientPosition"));
    calibrationInfo.ScanCenter(processingControl.Value<float>("ScanCenter"));
    calibrationInfo.FirstSliceCornerPoints(sliceTable.FirstSliceCorners());
    calibrationInfo.LastSliceCornerPoints(sliceTable.LastSliceCorners());
    calibrationInfo.Slices(sliceTable.SliceCount(0));

    return calibrationInfo;
}

PureCalibrationInfo GERecon::Calibration::CalibrationInfoForPure(const GERecon::Control::ProcessingControl& processingControl)
{
    const GERecon::SliceInfoTable& sliceTable = processingControl.ValueStrict<GERecon::SliceInfoTable>("SliceTable");

    PureCalibrationInfo calibrationInfo;

    calibrationInfo.Exam(processingControl.Value<unsigned int>("ExamNumber"));
    calibrationInfo.CoilID(processingControl.Value<unsigned int>("CoilConfigUID"));
    calibrationInfo.PatientPosition(processingControl.Value<int>("PatientPosition"));
    calibrationInfo.PatientEntry(processingControl.Value<int>("PatientEntry"));
    calibrationInfo.Landmark(processingControl.Value<float>("Landmark"));
    calibrationInfo.XRes(processingControl.Value<int>("TransformXRes"));
    calibrationInfo.YRes(processingControl.Value<int>("TransformYRes"));
    calibrationInfo.ScanCenter(processingControl.Value<float>("ScanCenter"));
    calibrationInfo.FermiRadius(processingControl.Value<float>("FermiRadius"));
    calibrationInfo.FermiWidth(processingControl.Value<float>("FermiWidth"));
    calibrationInfo.FirstSliceCornerPoints(sliceTable.FirstSliceCorners());
    calibrationInfo.LastSliceCornerPoints(sliceTable.LastSliceCorners());
    calibrationInfo.Slices(sliceTable.SliceCount(0));

    return calibrationInfo;
}

void GERecon::Calibration::AssetAndMcsiCalSliceProcessing(MDArray::FloatMatrix& virtualCoil, MDArray::ComplexFloatCube& sensitivityRatios, MDArray::FloatCube& magnitudeChannelData, 
                                                     const MDArray::ComplexFloatCube& channelData, const SliceOrientation& orientation, const int sliceNumber,
                                                     const GERecon::Control::ProcessingControl& processingControl)
{

    ComplexFloatCube normalizedChannelDataForAsset(channelData.shape());
    ComplexFloatCube normalizedChannelDataForMcsi(channelData.shape());    
    FloatMatrix rotateTransposeWorkspace(channelData.extent(firstDim), channelData.extent(secondDim));

    // Create magnitude data (weighted by normalized channedl weights). This magnitude data
    // is ultimately included in the Mcsi calibration file
    const FloatVector& normalizedChannelWeights = processingControl.ValueStrict<FloatVector>("NormalizedChannelWeights");
    normalizedChannelDataForMcsi = channelData;
    ChannelCombiner::NormalizeChannelData(normalizedChannelDataForMcsi, normalizedChannelWeights);
    MDArray::ComplexToMagnitude(magnitudeChannelData, normalizedChannelDataForMcsi);

    for(int channel = 0; channel < magnitudeChannelData.extent(thirdDim); ++channel)
    {
        rotateTransposeWorkspace = magnitudeChannelData(Range::all(), Range::all(), channel);
        FloatMatrix rotatedTransposeOutput = magnitudeChannelData(Range::all(), Range::all(), channel);
        RotateTranspose::Apply(rotatedTransposeOutput, rotateTransposeWorkspace, orientation.RotationType(), orientation.TransposeType());
    }

    // Normalize channel data based on channel weights. This normalized channel data is used
    // for computing data for the asset calibration file
    const FloatVector& channelWeights = processingControl.ValueStrict<FloatVector>("ChannelWeights");
    normalizedChannelDataForAsset = channelData;
    ChannelCombiner::NormalizeChannelData(normalizedChannelDataForAsset, channelWeights);

    // Compute reference "virtual coil" image
    GERecon::Asset::Calibration::ComputeVirtualCoil(virtualCoil, normalizedChannelDataForAsset);
        
    // Compute sensitivity ratios from virtual coil and channel data
    GERecon::Asset::Calibration::ComputeSensitivityRatios(sensitivityRatios, virtualCoil, normalizedChannelDataForAsset, sliceNumber);
}

void GERecon::Calibration::ProcessPureCalibration(const MDArray::Array<std::complex<float>, 3>& channelCombinedComplexSurfaceCoilVolume, const MDArray::ComplexFloatCube& complexReferenceCoilVolume,
                                                  const GERecon::Control::ProcessingControl& processingControl, GERecon::GradwarpPlugin& gradwarp, const boost::filesystem::path& outputDirectory,
                                                  const float refererenceCoilScalar)
{
    // Run pure calibration processing
    PureCalibrationFile pureCalFile(processingControl);
    const PureCalibrationInfo pureCalibrationInfo = CalibrationInfoForPure(processingControl);
    const boost::filesystem::path pureCalFilePath = outputDirectory / "PureCalibration.h5";
    pureCalFile.UseWritePaths(pureCalFilePath, pureCalFilePath);

    pureCalFile.WriteHeader(pureCalibrationInfo);

    const GERecon::SliceInfoTable& sliceTable = processingControl.ValueStrict<GERecon::SliceInfoTable>("SliceTable");
    const int numGeometricSlices = channelCombinedComplexSurfaceCoilVolume.extent(thirdDim);
    
    // Perform Gradwarp correction on each slice image.        
    FloatCube magnitudeVolumeCoilDataForPure(complexReferenceCoilVolume.extent(firstDim), complexReferenceCoilVolume.extent(secondDim), complexReferenceCoilVolume.extent(thirdDim));
    FloatCube channelCombinedMagnitudeSurfaceCoilDataForPure(channelCombinedComplexSurfaceCoilVolume.extent(firstDim), channelCombinedComplexSurfaceCoilVolume.extent(secondDim), channelCombinedComplexSurfaceCoilVolume.extent(thirdDim));
    for (int geometricSliceNumber = 0; geometricSliceNumber < numGeometricSlices; ++geometricSliceNumber)
    {  
        const ComplexFloatMatrix channelCombinedImage = channelCombinedComplexSurfaceCoilVolume(Range::all(), Range::all(), geometricSliceNumber);
        FloatMatrix surfaceSliceImage = channelCombinedMagnitudeSurfaceCoilDataForPure(Range::all(), Range::all(), geometricSliceNumber);
        MDArray::ComplexToReal(surfaceSliceImage, channelCombinedImage, MDArray::MagnitudeData);

        FloatMatrix referenceSliceImage = magnitudeVolumeCoilDataForPure(Range::all(), Range::all(), geometricSliceNumber);
        const ComplexFloatMatrix passSetCombinedImage = complexReferenceCoilVolume(Range::all(), Range::all(), geometricSliceNumber);
        MDArray::ComplexToReal(referenceSliceImage, passSetCombinedImage, MDArray::MagnitudeData);

        // Scale reference coil data
        referenceSliceImage *= refererenceCoilScalar;

        const GERecon::SliceCorners& sliceCorners = sliceTable.SliceCorners(geometricSliceNumber);

        gradwarp.Run(surfaceSliceImage,   sliceCorners, geometricSliceNumber);
        gradwarp.Run(referenceSliceImage, sliceCorners, geometricSliceNumber);
    }
    
    // Calculate Pure sensitivity ratios using PurePlugin
    FloatCube sensitivityRatios;
    PurePlugin::ProcessCalibration(sensitivityRatios, channelCombinedMagnitudeSurfaceCoilDataForPure, magnitudeVolumeCoilDataForPure, processingControl);

    // Write Pure sentivity data to file.
    pureCalFile.Write(pureCalibrationInfo, sensitivityRatios, channelCombinedMagnitudeSurfaceCoilDataForPure, magnitudeVolumeCoilDataForPure);
}

GERecon::Calibration::CalibrationData GERecon::Calibration::UnchopCalImageSpaceData(const GERecon::Calibration::CalibrationData& calImageSpaceData)
{
    GERecon::Calibration::CalibrationData unchoppedCalImageSpaceData(calImageSpaceData.shape());
    unchoppedCalImageSpaceData = calImageSpaceData; 
    unchoppedCalImageSpaceData = MDArray::where(tensor::l % 2 == 1, 
        -unchoppedCalImageSpaceData, unchoppedCalImageSpaceData);
    unchoppedCalImageSpaceData = MDArray::where(tensor::i % 2 == tensor::j % 2, 
        unchoppedCalImageSpaceData, -unchoppedCalImageSpaceData);
    return unchoppedCalImageSpaceData; 
}
