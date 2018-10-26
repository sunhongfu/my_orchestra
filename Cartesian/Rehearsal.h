// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <sstream>
#include <string>

#include <boost/shared_ptr.hpp>

/**
 * This header defines a list of functions that act as simple
 * recon pipelines (rehearsals) that can be called from the main
 * method in the corresponding source .cpp file.
 *
 * Define any new pipelines here and implement in a new
 * .cpp file. A typical use case would be to copy one
 * of the existing pipelines and modify it for development.
 *
 * The file contains a few helper functions useful for basic
 * pipeline creation and control.
 *
 * Also, note that everything is nested in the GERecon namespace.
 * This is the convention that is seen throughout all Orchestra
 * code. Namespaces allow for components/classes to be scoped
 * appropriately. If it lives in Orchestra, it's probably nested
 * somewhere in the GERecon namespace. Example: GERecon::Cartesian3D
 *
 * @author Jason Darby
 */
namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    /**
     * Simple 2D pipeline to recon one Pfile worth of data.
     * This function uses the Pfile to generate an Orchestra
     * style ProcessingControl instead of the Lx raw header
     * to pull out parameters for the reconstruction.
     *
     * Limitations: 
     *  No parallel imaging (Asset or Arc)
     *  No retrospective phase correction     
     */
    void RunSimple2D(const Legacy::PfilePointer& pfile);

    /**
     * Simple 3D pipeline to recon one Pfile worth of data.
     * This function uses the Pfile to generate an Orchestra
     * style ProcessingControl instead of the Lx raw header
     * to pull out parameters for the reconstruction.
     *
     * Limitations: 
     * No Asset support
     */
    void RunSimple3D(const Legacy::PfilePointer& pfile);

    /**
     * Get the image number for 2D based on slice/echo/phase index.
     * Image numbering scheme:
     * P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
     * P1S0E0, P1S0E1, ... PnSnEn
     */
    int ImageNumber2D (const int phase, const int slice, const int echo, const Legacy::PfilePointer& pfile);

    /**
     * Get the image number for 3D based on pass/slice/echo information.
     * Image numbering scheme (P = Phase; S = Slice; E = Echo):
     * P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
     * P1S0E0, P1S0E1, ... PnSnEn
     */
    int ImageNumber3D(const int pass, const int slice, const int echo, const Legacy::PfilePointer& pfile);

    /**
     * Convenience function for generating a filename for an output file.
     * The default generated string is of the form:
     * 
     * ExamXXXXX_SeriesYY_ImageZZZ.dcm
     * 
     * Where XXXXX, YY, ZZZ represent the numeric information specified.
     * The length of each numeric set of characters can be adjusted as
     * can the extension of the file.
     */
    std::string GenerateFileName(const int exam, const int series, const int image,
                                 const std::string& extension = ".dcm", const int examChars = 5, 
                                 const int seriesChars = 2, const int imageChars = 3);
}
