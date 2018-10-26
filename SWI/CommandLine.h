// Copyright 2014 General Electric Company. All rights reserved.

#pragma once

#include <string>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>

/**
 * Forwared declare the GEDicom::Network class and GEDicom::NetworkPointer typedef
 */
namespace GEDicom
{
    class Network;
    typedef boost::shared_ptr<Network> NetworkPointer;
}

namespace GERecon
{
    /**
     * Class that contains utilties for parsing parameters/values/flags from
     * the command line for usage in simple programs. The class requires the
     * GESystem::ProgramOptions to be initialized after main(...):
     * Example:
     * 
     *   int main(const int argc, const char* const argv[])
     *   {
     *       GESystem::ProgramOptions().SetupCommandLine(argc, argv);
     *      
     *       // code...
     *
     *       return 0;
     *   }
     *
     * @author Matt Bingen
     */
    class CommandLine
    {
    public:

        /**
         * Get the Pfile path specified on the command line. If it is not set
         * or does not exist, the function will throw an exception.
         *
         * Usage:
         *   --pfile </path/to/pfile>
         */
        static boost::filesystem::path PfilePath();

        /**
         * Flag to output spectral bins as h5 file If it is not set
         * the boost::optional will be empty.
         *
         * Usage:
         *   --outputBins 0/1
         */
        //static boost::optional<int> OutputSpectralBins();

        /**
         * Get the name and path of the config file
         *
         * Usage:
         *   --config <description>
         */
        static boost::optional<std::string> ConfigFile();

        /**
         * Get a DICOM network from parameters passed on the command line. If all
         * parameters are not set or the network cannot be create an empty pointer
         * will be returned.
         *
         * Usage:
         *   --ip <ip address of peer> --port <port #> --peer <peer AE title> --title <local AE title>
         *
         * Example:
         *   --ip 3.7.25.18 --port 4006 --peer t18 --title ese
         */
         // static GEDicom::NetworkPointer DicomNetwork();

    private:

        /**
         * Constructor - do not allow.
         */
        CommandLine();
    };
}
