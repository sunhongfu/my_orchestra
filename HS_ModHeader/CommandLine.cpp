// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>

#include <System/Utilities/ProgramOptions.h>

#include <Dicom/Network.h>

#include <Orchestra/Common/ReconException.h>

#include "CommandLine.h"

// testing the version control here.

using namespace GERecon;

boost::filesystem::path CommandLine::PfilePath()
{
    // Create option to get Pfile path. Additional options can be registered and loaded here.
    boost::program_options::options_description options;

    options.add_options()
        ("pfile", boost::program_options::value<std::string>(), "Specify pfile to run.");

    const GESystem::ProgramOptions programOptions;
    programOptions.AddOptions(options);

    // Check if the command line has a "--pfile" option
    const boost::optional<std::string> pfileOption = programOptions.Get<std::string>("pfile");

    if(!pfileOption)
    {
        throw GERecon::Exception(__SOURCE__, "No input Pfile specified! Use '--pfile' on command line.");
    }

    // Get path from string option
    const boost::filesystem::path pfilePath = *pfileOption;

    if(!boost::filesystem::exists(pfilePath))
    {
        throw GERecon::Exception(__SOURCE__, "Pfile [%s] doesn't exist!", pfilePath.string());
    }

    return pfilePath;
}

boost::optional<int> CommandLine::SeriesNumber()
{
    boost::program_options::options_description options;

    options.add_options()
        ("series", boost::program_options::value<int>(), "Series number to create images into");

    const GESystem::ProgramOptions programOptions;
    programOptions.AddOptions(options);

    return programOptions.Get<int>("series");
}

boost::optional<std::string> CommandLine::SeriesDescription()
{
    boost::program_options::options_description options;

    options.add_options()
        ("description", boost::program_options::value<std::string>(), "Series description");

    const GESystem::ProgramOptions programOptions;
    programOptions.AddOptions(options);

    return programOptions.Get<std::string>("description");
}

GEDicom::NetworkPointer CommandLine::DicomNetwork()
{
    // Create option to get Dicom Network info.
    boost::program_options::options_description options;

    options.add_options()
        ("ip", boost::program_options::value<std::string>(), "Peer IP Address")
        ("port", boost::program_options::value<unsigned short>(), "Peer Port")
        ("peer", boost::program_options::value<std::string>(), "Peer AE Title")
        ("title", boost::program_options::value<std::string>(), "Local/Host AE Title");

    const GESystem::ProgramOptions programOptions;
    programOptions.AddOptions(options);

    // Check if the command line options have been specified
    const boost::optional<std::string> ip = programOptions.Get<std::string>("ip");
    const boost::optional<unsigned short> port = programOptions.Get<unsigned short>("port");
    const boost::optional<std::string> peer = programOptions.Get<std::string>("peer");
    const boost::optional<std::string> title = programOptions.Get<std::string>("title");

    if(ip && port && peer && title)
    {
        std::cout << "Sending DICOM Image to:" << std::endl;
        std::cout << "IP: " << *ip << " Port: " << *port << " Peer: " << *peer << " Title: " << *title << std::endl;
        return boost::make_shared<GEDicom::Network>(*ip, *port, *title, *peer, true);
    }

    return boost::shared_ptr<GEDicom::Network>();
}
