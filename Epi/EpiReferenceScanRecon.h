// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/shared_ptr.hpp>

namespace GERecon
{
    class ScanArchive;
    typedef boost::shared_ptr<ScanArchive> ScanArchivePointer;

    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    namespace Epi
    {
        /**
         * Epi Reference Scan Pipeline
         *
         * @param pfile  Pfile containing raw data and header information
         * @return true if application completed successfully
         */
        bool EpiReferenceScanRecon(const GERecon::Legacy::PfilePointer& pfile);

        /**
         * Epi Reference Scan recon from scan archive
         * Currently only Diffusion hyper packets are supported
         *
         * @param scanArchive
         */
        bool EpiReferenceScanRecon(const GERecon::ScanArchivePointer& scanArchive);

    }
}
