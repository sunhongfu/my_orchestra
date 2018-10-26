// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include <Orchestra/Common/GradwarpMode.h>

namespace GERecon
{
    namespace Calibration
    {
        enum PipelineOptions
        {
            Calibration2DPipeline = 632,
            Calibration3DPipeline = 854,
            InvalidPipeline = 999
        };

        void AddOptions();

        std::string PfilePath();

        PipelineOptions PipelineToRun();

        GERecon::GradientType GradientType();
    }
}
