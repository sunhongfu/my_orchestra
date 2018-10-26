// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <Orchestra/Common/GradwarpMode.h>

namespace GERecon
{
    namespace Calibration
    {
        /**
         * The main recon function
         */
        void Calibration3DRecon(const std::string& pfileName, const GERecon::GradientType gradientType);
    }
}

