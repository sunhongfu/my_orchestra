// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

namespace GERecon
{
    namespace Spectro
    {
        enum PipelineOptions
        {
            MultiChannelMultiVoxelReconPipeline = 76,
            PerChannelMultiVoxelPipeline = 81,
            SingleVoxelReconPipeline = 88,
            InvalidPipeline = 999
        };

        void AddOptions();

        std::string PfilePath();

        PipelineOptions PipelineToRun();
    }
}
