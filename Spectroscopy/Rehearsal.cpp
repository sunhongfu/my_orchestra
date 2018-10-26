// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <cstdlib>
#include <exception>
#include <iostream>

#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Common/ReconException.h>
#include <Orchestra/Common/ReconPaths.h>

#include "Rehearsal.h"
#include "MultiChannelMultiVoxel.h"
#include "PerChannelMultiVoxel.h"
#include "SingleVoxel.h"

using namespace GERecon::Spectro;

int main(const int argc, const char* const argv[])
{
    try
    {
        GESystem::ProgramOptions().SetupCommandLine(argc, argv);
        AddOptions();

        const std::string pfilePath = PfilePath();
        const GERecon::Spectro::PipelineOptions pipelineOption = PipelineToRun();

        switch(pipelineOption)
        {
        case MultiChannelMultiVoxelReconPipeline:
            GERecon::Spectro::MultiChannelMultiVoxelRecon(pfilePath);
            break;
        case PerChannelMultiVoxelPipeline:
            GERecon::Spectro::PerChannelMultiVoxelRecon(pfilePath);
            break;
        case SingleVoxelReconPipeline:
            GERecon::Spectro::Rehearsal::SingleVoxelRecon(pfilePath);
            break;
        default:
            // Invalid pipeline option displayed in PipelineToRun function
            // Give the user an opportunity to read the message and continue
            std::cout << "Press enter to continue" << std::endl;
            std::cin.ignore();
        }

        std::cout << "Spectro Recon Complete" << std::endl;
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
void GERecon::Spectro::AddOptions()
{
    boost::program_options::options_description options;
    options.add_options()
        ("pfile", boost::program_options::value<std::string>(), "Specify pfile to run.")
        ("pipeline", boost::program_options::value<std::string>(), "Specify pipeline to run [MultiChannelMultiVoxel, PerChannelMultiVoxel, SingleVoxel]");

    GESystem::ProgramOptions().AddOptions(options);
}

// Get the Pfile path from command line
std::string GERecon::Spectro::PfilePath()
{
    // Check if the command line has a "--pfile" option
    const boost::optional<std::string> pfileOption = GESystem::ProgramOptions().Get<std::string>("pfile");

    if(!pfileOption)
    {
        GESystem::ProgramOptions().PrintUsage();

        throw GERecon::Exception(__SOURCE__, "No input Pfile specified! Use '--pfile' on command line.");
    }

    if(!boost::filesystem::exists(*pfileOption))
    {
        throw GERecon::Exception(__SOURCE__, "Pfile [%s] doesn't exist!", *pfileOption);
    }

    return *pfileOption;
}

// Get the pipeline to run from the command line
GERecon::Spectro::PipelineOptions GERecon::Spectro::PipelineToRun()
{
    // Check if the command line has a "--pipeline" option
    const boost::optional<std::string> pipelineOption = GESystem::ProgramOptions().Get<std::string>("pipeline");

    if(!pipelineOption)
    {
        // Pipeline not specified, print usage and return InvalidPipeline
        std::cout << "No pipeline specified." << std::endl;
        std::cout << "Specify --pipeline [MultiChannelMultiVoxel, PerChannelMultiVoxel, SingleVoxel] on the command line" << std::endl;
        return GERecon::Spectro::InvalidPipeline;
    }

    std::string pipelineOptionString = *pipelineOption;

    if(pipelineOptionString.compare("MultiChannelMultiVoxel") == 0)
    {
        return GERecon::Spectro::MultiChannelMultiVoxelReconPipeline;
    }
    else if(pipelineOptionString.compare("PerChannelMultiVoxel") == 0)
    {
        return GERecon::Spectro::PerChannelMultiVoxelPipeline;
    }
    else if(pipelineOptionString.compare("SingleVoxel") == 0)
    {
        return GERecon::Spectro::SingleVoxelReconPipeline;
    }
    else
    {
        std::cout << "Invalid pipeline specified." << std::endl;
        std::cout << "Specify --pipeline [MultiChannelMultiVoxel, PerChannelMultiVoxel, SingleVoxel] on the command line" << std::endl;
        return GERecon::Spectro::InvalidPipeline;
    }
}
