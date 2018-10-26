// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/filesystem/path.hpp>

#include <boost/shared_ptr.hpp>

#include <MDArray/MDArray.h>

namespace GERecon
{
    class SliceInfoTable;
    class GradwarpPlugin;
    class ScanArchive;
    typedef boost::shared_ptr<ScanArchive> ScanArchivePointer;

    namespace Legacy
    {
        class DicomSeries;
    }

    namespace Epi
    {
        /**
         * Reconstruct the given EPI ScanArchive.
         *
         * @param scanArchive
         */
        void EpiScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive);

        /**
         * Finalize a single image
         */
        void FinalizeImage(const MDArray::ComplexFloatMatrix& imageToFinalize, 
                           const int kissoffViews, 
                           const GERecon::SliceInfoTable& sliceTable, 
                           GERecon::GradwarpPlugin& gradwarpPlugin,
                           const GERecon::Legacy::DicomSeries& dicomSeries,
                           const float finalImageScaler,
                           const int phaseIndex,
                           const int geometricSliceIndex,
                           const boost::filesystem::path& dicomDirectory);
    }
}
