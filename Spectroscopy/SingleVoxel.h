// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <string>

#include <Hdf5/File.h>

#include <MDArray/MDArray.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Spectro/SingleVoxel/SingleVoxelQuant.h>
#include <Orchestra/Spectro/SingleVoxel/SingleVoxelWorker.h>
#include <Orchestra/Spectro/SingleVoxel/SingleVoxelParameters.h>

namespace GERecon
{
    namespace Spectro
    {
        namespace Rehearsal
        {
            /**
             * Entry point to the single voxel rehearsal pipeline
             *
             * @param pfilePath
             */
            void SingleVoxelRecon(const std::string& pfilePath);
            
            void ProcessFid(MDArray::ComplexFloatVector& cumulativePhaseCorrectionVector, MDArray::ComplexFloatVector& fidToProcess, const GERecon::Spectro::SingleVoxel::Options& processingParameters, const int transformSize);

            void GenerateChannelImages(const MDArray::ComplexFloatVector& combinedSignalFid, const MDArray::ComplexFloatVector& combinedReferenceFid, 
                                       const float lowestPpmToPlot, const float ppmRangeToPlot, const GERecon::Control::ProcessingControl& processingControl,
                                       const GEHdf5::File& svDebugFile, const int channelIndex);
        }
    }
}
