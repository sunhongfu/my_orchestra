#include <boost/filesystem.hpp>  
#include <System/Utilities/ProgramOptions.h>

#include <Orchestra/Legacy/PfileReader.h>  
//#include <Orchestra/Legacy/LegacyRdbm.h>
#include <Orchestra/Legacy/LegacyImageDB.h>

#include <iostream>
#include <string>
//
//#include <MDArray/Fourier.h>
//#include <MDArray/Utils.h>
//
//#include <Orchestra/Common/ImageCorners.h>
//#include <Orchestra/Common/ReconPaths.h>
//#include <Orchestra/Common/SliceInfoTable.h>
//#include <Orchestra/Common/SliceOrientation.h>
//
//#include <Orchestra/Control/ProcessingControl.h>
//
//#include <Orchestra/Core/Clipper.h>
//#include <Orchestra/Core/RotateTranspose.h>
//#include <Orchestra/Core/SumOfSquares.h>
//
//#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CommandLine.h"
//#include "Rehearsal.h"
//#include <boost/hof/lambda.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>

// From GERecon::RunSimple2D() or main(...)
using namespace GERecon;
using namespace MDArray;


int main(const int argc, const char* const argv[])
{
    GESystem::ProgramOptions().SetupCommandLine(argc, argv);
    
    // Read Pfile from command line
    const boost::filesystem::path pfilePath = CommandLine::PfilePath();
    
    // Create PfileReader to read the original raw header. Do not anonymize
    // any fields to leave it up to the user of the new modified Pfile.
    const AnonymizationPolicy noAnonymizationPolicy(AnonymizationPolicy::None);
    Legacy::PfileReader pfileReader(pfilePath, false, noAnonymizationPolicy, 25.002f);
    
    
    // Read POOL Header
    Legacy::PoolHeaderStruct header;
    pfileReader.ReadHeader(header);
    std::cout << "rdb_hdr_image.psdname is " << header.rdb_hdr_image.psdname << std::endl;
    std::cout << "rdb_hdr_image.psd_iname is " << header.rdb_hdr_image.psd_iname << std::endl;
    std::cout << "rdb_hdr_image.psdnameannot is " << header.rdb_hdr_image.psdnameannot << std::endl;

//    std::cout << "rdb_hdr_image.proj_name is " << header.rdb_hdr_image.proj_name << std::endl;
    std::cout << "rdb_hdr_image.cname is " << header.rdb_hdr_image.cname << std::endl;
    std::cout << "rdb_hdr_image.GEcname is " << header.rdb_hdr_image.GEcname << std::endl;
    std::cout << "rdb_hdr_series.prtcl is " << header.rdb_hdr_series.prtcl << std::endl;
    std::cout << "rdb_hdr_series.se_desc is " << header.rdb_hdr_series.se_desc << std::endl;


    std::cout << "rdb_hdr_rec.rdb_hdr_asset is " << header.rdb_hdr_rec.rdb_hdr_asset << std::endl;
    std::cout << "rdb_hdr_rec.rdb_hdr_asset_R is " << header.rdb_hdr_rec.rdb_hdr_asset_R << std::endl;
    std::cout << "rdb_hdr_rec.rdb_hdr_asset_slice_R is " << header.rdb_hdr_rec.rdb_hdr_asset_slice_R << std::endl;
    std::cout << "rdb_hdr_rec.rdb_hdr_image.ti is " << header.rdb_hdr_image.ti << std::endl;
    std::cout << "rdb_hdr_rec.rdb_hdr_kacq_uid is " << header.rdb_hdr_rec.rdb_hdr_kacq_uid << std::endl; // kacq uid, e.g. kacq_yz.txt.1104145905
    std::cout << "rdb_hdr_rec.rdb_hdr_user16 is " << header.rdb_hdr_rec.rdb_hdr_user16 << std::endl;
    // this is the uniPolar_flag, 1 means uni, otherwiese bi, for the efgre3d based sequences
    
    
    
    

    
    return 0;
    
}
