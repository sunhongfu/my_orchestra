// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>

#include <MDArray/MDArray.h>

#include <Orchestra/Control/ProcessingControl.h>

namespace GERecon
{
    namespace Spectro
    {
        /**
         * Entry point to the Multi-Channel Multi-Voxel Rehearsal Pipeline
         * @param pfileName
         */
        void MultiChannelMultiVoxelRecon(const std::string& pfileName);
    }
}
