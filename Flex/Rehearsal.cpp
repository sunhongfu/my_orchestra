// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <iostream>
#include <exception>

#include <System/Utilities/ProgramOptions.h>

#include "CommandLine.h"
#include "Rehearsal.h"

using namespace GERecon;
    
/*****************************************************************
 ** Main function that calls the specific recon pipeline to run **
 ******************************************************************/
int main(const int argc, const char* const argv[])
{
    try
    {
        GESystem::ProgramOptions().SetupCommandLine(argc, argv);

        // Two available APIs.  RunFlexSimple() is easier to understand and follow, but more memory demanding.
        if( CommandLine::RunSimple() )
        {
            RunFlexSimple();
        }
        else
        {
            RunFlex();
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
