// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/filesystem.hpp>

#include <MDArray/MDArray.h>

#include <Orchestra/Asset/CalibrationInfo.h>
#include <Orchestra/Calibration/Common/Types.h>
#include <Orchestra/Common/SliceOrientation.h>
#include <Orchestra/Control/ProcessingControl.h>
#include <Orchestra/Gradwarp/GradwarpPlugin.h>
#include <Orchestra/Pure/Pure1/PureCalibrationInfo.h>

namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
    }

    namespace Calibration
    {
        /**
         * Run Asset Calibration Processing
         */
        void ProcessAssetCalibration(const MDArray::Array<std::complex<float>, 4>& surfaceCoilVolume,  const MDArray::Array<std::complex<float>, 4>& alternatelyApodizedVolume, 
                                     const GERecon::Control::ProcessingControl& processingControl, const boost::filesystem::path& outputDirectory);

        /**
         * Slice by slice portion of Asset Calibration processing
         */
        void AssetAndMcsiCalSliceProcessing(MDArray::FloatMatrix& virtualCoil, MDArray::ComplexFloatCube& sensitivityRatios, MDArray::FloatCube& magnitudeChannelData, const MDArray::ComplexFloatCube& channelData,
                                            const GERecon::SliceOrientation& orientation, const int sliceNumber, const GERecon::Control::ProcessingControl& processingControl);

        /**
         * Fill Calibration Info Structure for ASSET
         */
        Asset::CalibrationInfo AssetCalibrationInfo(const GERecon::Control::ProcessingControl& processingControl);

        /**
         * Run Pure Calibration Processing
         */
        void ProcessPureCalibration(const MDArray::Array<std::complex<float>, 3>& channelCombinedComplexSurfaceCoilVolume, const MDArray::ComplexFloatCube& complexReferenceCoilVolume,
                                    const GERecon::Control::ProcessingControl& processingControl, GERecon::GradwarpPlugin& gradwarp, const boost::filesystem::path& outputDirectory,
                                    const float refererenceCoilScalar = 1.0f);

        /**
         * Fill Calibration Info Structure for PURE
         */
        PureCalibrationInfo CalibrationInfoForPure(const GERecon::Control::ProcessingControl& processingControl);

        /**
         * Unchop calibration data in image space: 
         * <i> Elements with odd indices along the fourth dimension (Z-Axis) will be multiplied by -1; 
         * <ii> Elements with 'unmatched' indices along the first and the second dimensions (i.e., one index is odd while the other is even) will be multiplied by -1. 
         * @param calibration data to be unchopped
         * @return resulting unchopped calibration data 
         */
        GERecon::Calibration::CalibrationData UnchopCalImageSpaceData(const CalibrationData& calImageSpaceData); 
    }
}


