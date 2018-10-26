// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <cstdlib>
#include <exception>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Common/ReconException.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/ScanArchive.h>

#include <Orchestra/Legacy/Pfile.h>

#include "EpiDiffusionRecon.h"
#include "EpiDiffusionScanArchiveRecon.h"
#include "EpiScanArchiveRecon.h"
#include "EpiMultiPhaseRecon.h"
#include "EpiRecon.h"
#include "EpiReferenceScanRecon.h"
#include "Rehearsal.h"

using namespace GERecon::Epi;

/**
 * EPI Rehearsal Recon Pipelines
 * 
 * @author Eric Printz
 */
int main(const int argc, const char* const argv[])
{
    try
    {
        GESystem::ProgramOptions().SetupCommandLine(argc, argv);
        AddOptions();

        const boost::filesystem::path inputFilePath = InputFilePath();
        SetPathsToInputFileDirectory(inputFilePath);

        // Determine the EPI pipeline to run.
        if(GERecon::ScanArchive::IsArchiveFilePath(inputFilePath))
        {
            // If this is a scan archive, run the scan archive pipeline
            RunScanArchiveRehearsalPipeline(inputFilePath);
        }
        else
        {
            // The input data must be a pfile. Determine the pipeline to run
            // based on the input file and the command line. If the user specified
            // a specific pipeline on the command line, that pipeline will be used.
            // If a pipeline was not specified on the command line, the raw header
            // will be used to determine the pipeline to run
            RunPfileRehearsalPipeline(inputFilePath);
        }

        std::cout << "EPI Recon Complete" << std::endl;
    }
    catch( std::exception& e )
    {
        std::cout << "Runtime Exception! " << e.what() << std::endl;
        std::cin.ignore();
    }
    catch( ... )
    {
        std::cout << "Unknown Runtime Exception!" << std::endl;
        std::cin.ignore();
    }

    return 0;
}

// Add valid command line options
void GERecon::Epi::AddOptions()
{
    boost::program_options::options_description options;
    options.add_options()
        ("file", boost::program_options::value<std::string>(), "Specify input file (Pfile or ScanArchive) to run. The pfile option should be omitted from the command line if this option is used.")
        ("pfile", boost::program_options::value<std::string>(), "Specify pfile to run.")
        ("pipeline", boost::program_options::value<std::string>(), "Specify pipeline to run.");

    GESystem::ProgramOptions().AddOptions(options);
}

// Get the Pfile path from command line
boost::filesystem::path GERecon::Epi::InputFilePath()
{
    // Check if the command line has a "--pfile" option
    const boost::optional<std::string> pfileOption = GESystem::ProgramOptions().Get<std::string>("pfile");
    const boost::optional<std::string> fileOption = GESystem::ProgramOptions().Get<std::string>("file");
    boost::filesystem::path inputFilePath;

    // First check for the file option
    if(fileOption)
    {
        inputFilePath = *fileOption;
    }
    else if(pfileOption)
    {
        inputFilePath = *pfileOption;
    }
    else
    {
        GESystem::ProgramOptions().PrintUsage();
        throw GERecon::Exception(__SOURCE__, "No input file specified! Use '--file' on the command line.");
    }

    if(!boost::filesystem::exists(inputFilePath))
    {
        throw GERecon::Exception(__SOURCE__, "Pfile [%s] doesn't exist!", *pfileOption);
    }

    return inputFilePath;
}

void GERecon::Epi::SetPathsToInputFileDirectory(const boost::filesystem::path& inputFilePath)
{
    // Set the InputAppData directory to be the same as location of Pfile supplied on command line.
    // The ramp sampling and phase correction plugins will look for data such as ref.dat and vrgf.dat
    // in the InputAppData recon path.
    const boost::filesystem::path inputFileDirectory = inputFilePath.parent_path();
    GERecon::Path::InputAppData(inputFileDirectory);

    // Set the OutputAppData directory to the Pfile directory. This will 
    // cause ref.h5 to be saved to the Pfile directory when a reference scan
    // recon is run
    GERecon::Path::OutputAppData(inputFileDirectory);
}

// Get the pipeline to run from the command line
void GERecon::Epi::RunPfileRehearsalPipeline(const boost::filesystem::path& inputFilePath)
{
    // Create a pfile object based on the input file path
    const GERecon::Legacy::PfilePointer pfile = Legacy::Pfile::Create(inputFilePath);

    // Check if the command line has a "--pipeline" option
    const boost::optional<std::string> pipelineOption = GESystem::ProgramOptions().Get<std::string>("pipeline");
    const bool pipelineOptionExists = pipelineOption ? true : false;
    if( (pipelineOptionExists && pipelineOption->compare("EpiDiffusionRecon") == 0) || (!pipelineOptionExists && pfile->IsDiffusionEpi()) )
    {
        // If the pipeline option is specified and is EpiDiffusionRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiDiffusionRecon
        EpiDiffusionRecon(pfile);
    }
    else if( (pipelineOptionExists && pipelineOption->compare("EpiMultiPhaseRecon") == 0) || (!pipelineOptionExists && (pfile->IsMultiPhaseEpi() || pfile->IsFunctionalMri())) )
    {
        // If the pipeline option is specified and is EpiMultiPhaseRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiMultiPhaseRecon
        EpiMultiPhaseRecon(pfile);
    }
    else if( (pipelineOptionExists && pipelineOption->compare("EpiReferenceScanRecon") == 0) || (!pipelineOptionExists && pfile->IsEpiRefScan()) )
    {
        // If the pipeline option is specified and is EpiReferenceScanRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiReferenceScanRecon
        EpiReferenceScanRecon(pfile);
    }
    else if( (pipelineOptionExists && pipelineOption->compare("EpiRecon") == 0) || (!pipelineOptionExists && pfile->IsEpi()) )
    {
        // If the pipeline option is specified and is EpiRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiRecon
        EpiRecon(pfile);
    }
    else
    {
        std::cout << "Could not determine EPI rehearsal pipeline to run. Use the --pipeline option to force a given EPI pipeline to run" << std::endl;
        GESystem::ProgramOptions().PrintUsage();
    }
}

// Get the pipeline to run from the command line
void GERecon::Epi::RunScanArchiveRehearsalPipeline(const boost::filesystem::path& inputFilePath)
{
    // If this is a scan archive, run the scan archive pipeline
    const GERecon::ScanArchivePointer scanArchive = GERecon::ScanArchive::Create(inputFilePath, GESystem::Archive::LoadMode);

    const GERecon::Legacy::LxDownloadDataPointer downloadData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchive->LoadDownloadData());

    // Check if the command line has a "--pipeline" option
    const boost::optional<std::string> pipelineOption = GESystem::ProgramOptions().Get<std::string>("pipeline");
    const bool pipelineOptionExists = pipelineOption ? true : false;
    if( (pipelineOptionExists && pipelineOption->compare("EpiDiffusionRecon") == 0) || (!pipelineOptionExists && downloadData->IsDiffusionEpi()) )
    {
        // If the pipeline option is specified and is EpiDiffusionRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiDiffusionRecon
        GERecon::Epi::EpiDiffusionScanArchiveRecon(scanArchive);
    }
    else if( (pipelineOptionExists && pipelineOption->compare("EpiRecon") == 0) || (!pipelineOptionExists && downloadData->IsEpi()) )
    {
        // If the pipeline option is specified and is EpiRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiRecon
        GERecon::Epi::EpiScanArchiveRecon(scanArchive);
    }
    else if( (pipelineOptionExists && pipelineOption->compare("EpiReferenceScanRecon") == 0) || (!pipelineOptionExists && downloadData->IsEpiRefScan()) )
    {
        // If the pipeline option is specified and is EpiReferenceScanRecon -or- if the pipeline option is not specified and the raw header indicates this is EpiReferenceScanRecon
        GERecon::Epi::EpiReferenceScanRecon(scanArchive);
    }
    else
    {
        std::cout << "Could not determine EPI rehearsal pipeline to run. Use the --pipeline option to force a given EPI pipeline to run" << std::endl;
        GESystem::ProgramOptions().PrintUsage();
    }
}