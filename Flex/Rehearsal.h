// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include <MDArray/MDArray.h>

#include <Orchestra/Common/ReconTrace.h>
#include <Orchestra/Gradwarp/GradwarpPlugin.h>

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
 * somewhere in the GERecon namespace. Example: GERecon::Flex
 *
 * @author Matt Bingen, Kang Wang
 */
namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    /**
     * Simple Flex pipeline to recon one Pfile worth of data.
     * With an image space 5D array as input argument,
     * it only takes one function call to do the 2-pt Dixon
     * water-fat separation. Simpler to use, but more memory
     * demanding.
     */
    void RunFlexSimple();

    /**
     * Flex pipeline to recon one Pfile worth of data,
     * and generate water and fat.
     * This function includes a few steps to do the entire
     * water-fat separation. It is less memory demanding,
     * but does have some repeated computations.
     */
    void RunFlex();

    /**
     * Generate water and fat DICOM images
     */
    void GenerateWaterFatDicom(const MDArray::Array<float, 4>& imageDataWaterFat,
                               const Legacy::PfilePointer& pfile,
                               GradwarpPlugin& gradwarp);
}
