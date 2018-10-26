// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <MDArray/MDArray.h>

#include <Orchestra/Control/ProcessingControl.h>

namespace GERecon
{
    namespace SpiralUtil
    {
        /**
         * Run a simple Spiral Recon of 4D raw data. No gradwarp and image clip/rotation/transpose
         * @param imageData         Reconstructed image output (imageSize x imageSize x numSlices)
         * @param rawData           Input kSpace data (numPointsPerArm x numArms x numChannels x numSlices)
         * @param processingControl
         * @param dataDirectory     Directory for saved trajectory file
         * @param scalingFactor
         */
        void SpiralRecon(MDArray::FloatCube& imageData,
            const MDArray::ComplexFloat4D& rawData,
            const Control::ProcessingControlPointer& processingControl,
            const std::string& dataDirectory,
            const float scalingFactor);
    }
}
