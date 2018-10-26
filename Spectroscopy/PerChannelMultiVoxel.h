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
         * Entry point to the Per-Channel Multi-Voxel Rehearsal Pipeline
         * @param pfileName
         */
        void PerChannelMultiVoxelRecon(const std::string& pfileName);
    }
}
