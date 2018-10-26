#include <boost/filesystem.hpp>  
#include <System/Utilities/ProgramOptions.h>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>

#include <Orchestra/Legacy/PfileReader.h>  
#include <Orchestra/Legacy/LegacyRdbm.h>  

#include <iostream>
#include <string>

#include <MDArray/Fourier.h>
#include <MDArray/Utils.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceInfoTable.h>
#include <Orchestra/Common/SliceOrientation.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>
#include <Orchestra/Core/SumOfSquares.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include "CommandLine.h"
#include "Rehearsal.h"
//#include <boost/hof/lambda.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>

// From GERecon::RunSimple2D() or main(...)
using namespace GERecon;
using namespace MDArray;

struct Print
{
    template<typename T>
    void operator()(T& t) const
    {
        std::cout << t << std::endl;
    }
};

int main(const int argc, const char* const argv[])
{
    GESystem::ProgramOptions().SetupCommandLine(argc, argv);
    
    // Read Pfile from command line
    const boost::filesystem::path pfilePath = CommandLine::PfilePath();
    
    // Create PfileReader to read the original raw header. Do not anonymize
    // any fields to leave it up to the user of the new modified Pfile.
    const AnonymizationPolicy noAnonymizationPolicy(AnonymizationPolicy::None);
    Legacy::PfileReader pfileReader(pfilePath, false, noAnonymizationPolicy);
    
    // Read POOL Header
    Legacy::PoolHeaderStruct header;
    pfileReader.ReadHeader(header);
    
//    std::cout << "rdb_hdr_rec.rdb_hdr_asset is " << header.rdb_hdr_rec << std::endl;

//    BOOST_HOF_STATIC_LAMBDA_FUNCTION(Print) = [](const auto& x)
//    {
//        std::cout << x << std::endl;
//    };


//    boost::fusion::for_each<_RDB_HEADER_REC>(header.rdb_hdr_rec, Print());
    
    // Set the value to 0 so the Pfile class will know how to correctly parse the header
    //header.rdb_hdr_rec.rdb_hdr_asset = 7;
    //header.rdb_hdr_rec.rdb_hdr_asset_R = 1.0;
    //header.rdb_hdr_rec.rdb_hdr_asset_slice_R = 1.0;
    //header.rdb_hdr_rec.rdb_hdr_channel_combine_method = 1.0;
    //header.rdb_hdr_rec.rdb_hdr_channel_combine_filter_width = 0;
    //header.rdb_hdr_rec.rdb_hdr_channel_combine_filter_beta = 0;
    
    // header.rdb_hdr_series.assetcal_serno = 1;
    // header.rdb_hdr_rec.rdb_hdr_SliceAsset = 1;
    // _SERIESDATATYPE.
    //                     short int assetcal_serno;
    //                     short int assetcal_scnno;
    //                     short int asset_cal_type;
    //                     char asset_appl[12];
    
    // _MRIMAGEDATATYPE
    //                     float SliceAsset;
    //                     float PhaseAsset;
    
    
    
    return 0;
    
}
