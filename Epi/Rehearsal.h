// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>

#include <boost/filesystem/path.hpp>

namespace GERecon
{
    namespace Epi
    {
        // Enumeration of pipeline options. Each assigned a random enum integer value.
        enum PipelineOptions
        {
            EpiReferenceScanReconPipeline = 125,
            EpiReconPipeline = 247,
            EpiMultiPhaseReconPipeline = 423,
            EpiDiffusionReconPipeline = 634,
            EpiScanArchivePipeline = 636,
            InvalidPipeline = 999
        };

        /**
         * Add known command line options to list of available options
         */
        void AddOptions();

        /**
         * Retrieve input file from command line
         */
        boost::filesystem::path InputFilePath();

        /**
         * Select and run a pfile rehearsal pipeline. If the pipeline command line option
         * is available, the pipeline to run is based on this option. If the pipeline command
         * line option is not available, the pipeline to run is determined by inspecting the
         * raw header contained in the input pfile
         */
        void RunPfileRehearsalPipeline(const boost::filesystem::path& inputFilePath);

        /**
         * Select and run a scan archive rehearsal pipeline. If the pipeline command line option
         * is available, the pipeline to run is based on this option. If the pipeline command
         * line option is not available, the pipeline to run is determined by inspecting the
         * raw header contained in the input scan archive
         */
        void RunScanArchiveRehearsalPipeline(const boost::filesystem::path& inputFilePath);

        /**
         * Set input/output paths to the input file path
         *
         * @param inputFilePath
         */
        void SetPathsToInputFileDirectory(const boost::filesystem::path& inputFilePath);
    }
}
