// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>

/**
 * Forward declare the GEDicom::Network class and GEDicom::NetworkPointer typedef
 */
namespace GEDicom
{
    class Network;
    typedef boost::shared_ptr<Network> NetworkPointer;
}

namespace GERecon
{
    /**
     * Class that contains utilities for parsing parameters/values/flags from
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
     * @author Jason Darby
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
         * Get the series number specified on the command line. If it is not set
         * the boost::optional will be empty.
         *
         * Usage:
         *   --series <series#>
         */
        static boost::optional<int> SeriesNumber();

        /**
         * Get the series description specified on the command line. If it is not set
         * or does not exist, an empty string is returned.
         *
         * Usage:
         *   --description <description>
         */
        static boost::optional<std::string> SeriesDescription();

        /**
         * Get a DICOM network from parameters passed on the command line. If all
         * parameters are not set or the network cannot be created, an empty pointer
         * will be returned.
         *
         * Usage:
         *   --ip <ip address of peer> --port <port #> --peer <peer AE title> --title <local AE title>
         *
         * Example:
         *   --ip 3.7.25.18 --port 4006 --peer t18 --title ese
         */
        static GEDicom::NetworkPointer DicomNetwork();

    private:

        /**
         * Constructor - do not allow.
         */
        CommandLine();
    };
}
