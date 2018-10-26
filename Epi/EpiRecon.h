// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/shared_ptr.hpp>

namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    namespace Epi
    {
        /* The main recon function */
        bool EpiRecon(const GERecon::Legacy::PfilePointer& pfile);
    }
}

