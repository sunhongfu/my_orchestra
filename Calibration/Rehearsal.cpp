// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <exception>
#include <iostream>
#include <cstdlib>

#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Common/ReconException.h>
#include <Orchestra/Common/ReconPaths.h>

#include "Pipeline/Calibration2DRecon.h"
#include "Pipeline/Calibration3DRecon.h"
#include "Rehearsal.h"

using namespace GERecon::Calibration;

/**
 * Calibration Rehearsal Recon Pipelines
 */
int main(const int argc, const char* const argv[])
{
    try
    {
        GESystem::ProgramOptions().SetupCommandLine(argc, argv);
        AddOptions();

        const std::string pfilePath = PfilePath();
        const GERecon::Calibration::PipelineOptions pipelineOption = PipelineToRun();
        const GERecon::GradientType gradientType = GradientType();

        switch(pipelineOption)
        {
        case Calibration2DPipeline:
            Calibration2DRecon(pfilePath, gradientType);
            break;
        case Calibration3DPipeline:
            Calibration3DRecon(pfilePath, gradientType);
            break;
        default:
            // Invalid pipeline option was displayed in PipelineToRun function
            // Give the user an opportunity to read the message and continue
            std::cout << "Press enter to continue" << std::endl;
            std::cin.ignore();
            break;
        }

        std::cout << "Calibration3D Recon Complete" << std::endl;
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
void GERecon::Calibration::AddOptions()
{
    boost::program_options::options_description options;
    options.add_options()
        ("pfile", boost::program_options::value<std::string>(), "Specify pfile to run.")
        ("pipeline", boost::program_options::value<std::string>(), "Specify pipeline to run. [Calibration2DRecon, Calibration3DRecon]")
        ("gradientType", boost::program_options::value<std::string>(), "Specify XRMB (standard) or XRMW (wide) gradient. Defaults to XRMB.");
        
    GESystem::ProgramOptions().AddOptions(options);
}

// Get the Pfile path from command line
std::string GERecon::Calibration::PfilePath()
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
GERecon::Calibration::PipelineOptions GERecon::Calibration::PipelineToRun()
{
    // Check if the command line has a "--pipeline" option
    const boost::optional<std::string> pipelineOption = GESystem::ProgramOptions().Get<std::string>("pipeline");

    if(!pipelineOption)
    {
        // Pipeline not specified, print usage and return InvalidPipeline
        std::cout << "No pipeline specified." << std::endl;
        std::cout << "Specify --pipeline [Calibration2DRecon, Calibration3DRecon] on the command line" << std::endl;
        return GERecon::Calibration::InvalidPipeline;
    }

    const std::string pipelineOptionString = *pipelineOption;

    if(pipelineOptionString.compare("Calibration2DRecon") == 0)
    {
        return GERecon::Calibration::Calibration2DPipeline;
    }
    else if(pipelineOptionString.compare("Calibration3DRecon") == 0)
    {
        return GERecon::Calibration::Calibration3DPipeline;
    }
    else
    {
        std::cout << "Invalid pipeline specified." << std::endl;
        std::cout << "Specify --pipeline [Calibration2DRecon, Calibration3DRecon] on the command line" << std::endl;
        return GERecon::Calibration::InvalidPipeline;
    }
}

GERecon::GradientType GERecon::Calibration::GradientType()
{
    // Check if the command line has a "--pipeline" option
    const boost::optional<std::string> gradientTypeOption = GESystem::ProgramOptions().Get<std::string>("gradientType");

    if(!gradientTypeOption)
    {
        return GERecon::XRMBGradient;
    }

    const std::string gradientTypeOptionString = *gradientTypeOption;

    if(gradientTypeOptionString.compare("XRMW") == 0)
    {
        return GERecon::XRMWGradient;
    }
    else
    {
        return GERecon::XRMBGradient;
    }
}
