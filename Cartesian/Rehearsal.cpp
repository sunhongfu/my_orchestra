// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <exception>
#include <iostream>

#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Common/ReconTrace.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CommandLine.h"
#include "Rehearsal.h"

using namespace GERecon;
    
/*****************************************************************
 ** Main function that calls the specific recon pipeline to run **
 ******************************************************************/
int main(const int argc, const char* const argv[])
{
    GESystem::ProgramOptions().SetupCommandLine(argc, argv);

    // Read Pfile from command line
    const boost::filesystem::path pfilePath = CommandLine::PfilePath();
    const Legacy::PfilePointer pfile = Legacy::Pfile::Create(pfilePath, Legacy::Pfile::AllAvailableAcquisitions, AnonymizationPolicy(AnonymizationPolicy::None));
    try
    {
        if(pfile->IsZEncoded())
        {
            RunSimple3D(pfile);
        }
        else
        {
            RunSimple2D(pfile);
        }
        std::cout << std::endl << "Application completed successfully. Hit Enter to quit..." << std::endl;
        std::cin.get();
        return 0;
    }
    catch( std::exception& e )
    {
        std::cout << "Runtime Exception! " << e.what() << std::endl;
        std::cout << "Hit Enter to quit...";
        std::cin.get();
    }
    catch( ... )
    {
        std::cout << "Unknown Runtime Exception!" << std::endl;
        std::cout << "Hit Enter to quit...";
        std::cin.get();
    }

    return -1;
}

std::string GERecon::GenerateFileName(const int exam, const int series, const int image,
                                      const std::string& extension, const int examChars, 
                                      const int seriesChars, const int imageChars)
{
    std::ostringstream fileName;

    fileName << std::setfill('0') 
            << "Exam" << std::setw(examChars) << exam 
            << "_Series" << std::setw(seriesChars) << series 
            << "_Image" << std::setw(imageChars) << image 
            << extension;

    return fileName.str();
}
