#include <boost/filesystem.hpp>  
#include <System/Utilities/ProgramOptions.h>
  
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
Legacy::PfileReader pfileReader(pfilePath, false, noAnonymizationPolicy);  
  
// Read POOL Header  
Legacy::PoolHeaderStruct header;  
pfileReader.ReadHeader(header);  

  
// Set the value to 0 so the Pfile class will know how to correctly parse the header  
header.rdb_hdr_rec.rdb_hdr_asset = 7;  
header.rdb_hdr_rec.rdb_hdr_asset_R = 1.0;  
header.rdb_hdr_rec.rdb_hdr_asset_slice_R = 1.0;
header.rdb_hdr_rec.rdb_hdr_channel_combine_method = 1.0;
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
  
// Output Pfile with user19 set to 0 to /path/to/orig/Pxxxxx.7.mod  
boost::filesystem::path newFile = pfilePath;  
const std::string newName = newFile.leaf().string() + ".mod";  
newFile.remove_leaf() /= newName;  
  
// Write out new Pfile  
pfileReader.WritePfile(newFile, &header, sizeof(Legacy::PoolHeaderStruct));  
  
// Create Pfile object for reading data in standard format  
const Legacy::PfilePointer pfile = Legacy::Pfile::Create(newFile, Legacy::Pfile::AllAvailableAcquisitions, AnonymizationPolicy(AnonymizationPolicy::None));  
  
// Continue normal recon pipeline...  
return 0;

}

