// Copyright 2018 General Electric Company. All rights reserved.

#include <exception>
#include <iostream>

#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Legacy/Pfile.h>

#include <Orchestra/Spiral/LxControlSource.h>
#include <Orchestra/Spiral/SpiralWorker.h>

#include "Rehearsal.h"

using namespace GERecon;
using namespace GERecon::Spiral;
using namespace GERecon::Legacy;

/*****************************************************************
 ** Main function that calls the specific recon pipeline to run **
 ******************************************************************/
int main(const int argc, const char* const argv[])
{
    try
    {        
        GESystem::ProgramOptions().SetupCommandLine(argc,argv);
        boost::program_options::options_description options;
        options.add_options()("file", boost::program_options::value<std::string>(), "Specify file.");
        const GESystem::ProgramOptions progOptions;
        progOptions.AddOptions(options);
        if(!progOptions.Get<std::string>("file"))
        {
            throw std::runtime_error("No input file specified. Use '--file'");
        }

        boost::filesystem::path inputFilePath = *progOptions.Get<std::string>("file");;

        // Determine the Spiral pipeline to run
        if(GERecon::ScanArchive::IsArchiveFilePath(inputFilePath))
        {
            const GERecon::ScanArchivePointer scanArchive = GERecon::ScanArchive::Create(inputFilePath, GESystem::Archive::LoadMode);
            SpiralScanArchiveRecon(scanArchive);

            std::cout << std::endl << "ScanArchive Application completed successfully. Hit Enter to quit..." << std::endl;
            std::cin.get();
            return 0;
        }
        else
        {
            const boost::shared_ptr<Spiral::LxControlSource> spiralControlSource = boost::make_shared<Spiral::LxControlSource>();
            boost::shared_ptr<Pfile> pFile = Pfile::Create(inputFilePath, spiralControlSource, Pfile::AllAcquisitions);
            const Control::ProcessingControlPointer reconControl = pFile->CreateOrchestraProcessingControl();

            const bool is3DASL = reconControl->Value<bool>("Is3DASL");
            if( is3DASL )
            {
                std::cout << "3DASL rehearsal only support ScanArchive input." << std::endl;
                std::cout << "Hit Enter to quit...";
                std::cin.get();
                return 0;
            }
            else
            {
                Spiral2D(reconControl, *pFile );

                std::cout << std::endl << "Spiral2D Application completed successfully. Hit Enter to quit..." << std::endl;
                std::cin.get();
                return 0;
            }
        }
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
