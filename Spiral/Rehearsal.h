// Copyright 2018 General Electric Company.  All rights reserved.

#pragma once

#include <MDArray/MDArray.h>

#include <Orchestra/Control/ProcessingControl.h>
#include <Orchestra/Legacy/LegacyForwardDeclarations.h>
#include <Orchestra/Legacy/LegacyImageDB.h>
#include <Orchestra/Spiral/SpiralWorker.h>

/**
 * This header defines a list of functions that act as simple
 * recon pipeline that can be called from the main method
 * in the corresponding source .cpp file.
 *
 * @author Alan Thompson
 */
namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
    }
    namespace Cartesian3D
    {
        class ZTransformer;
    }

    /**
     * Spiral 2D pipeline to recon one Pfile worth of data.
     * This function uses the Pfile to pull out parameters for the reconstruction.
     * The function is intended to show the use of the Spiral and other
     * common Orchestra algorithms and building blocks.
     *
     * @param pfileName name of pfile to load and run from
     */
    void Spiral2D(const Control::ProcessingControlPointer& processingControl, Legacy::Pfile& pfile);

    /**
     * Reconstruct the given Spiral ScanArchive.
     *
     * @param scanArchive
     */
    void SpiralScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive);
}
