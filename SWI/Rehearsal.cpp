// Copyright 2011 General Electric Company. All rights reserved.

#include <iostream>
#include <exception>

#include <System/Utilities/ProgramOptions.h>

#include "Rehearsal.h"


// Include this to avoid having to type fully qualified names
using namespace GERecon;
    
/*****************************************************************
 ** Main function that calls the specific recon pipeline to run **
 ******************************************************************/
int main(const int argc, const char* const argv[])
{
    GESystem::ProgramOptions().SetupCommandLine(argc, argv);

    try
    {

        RunSWI();

        std::cout << std::endl << "Application completed successfully" << std::endl;
        return 0;
    }
    catch( std::exception& e )
    {
        std::cout << "Runtime Exception! " << e.what() << std::endl;
    }
    catch( ... )
    {
        std::cout << "Unknown Runtime Exception!" << std::endl;
    }

    return -1;
}


