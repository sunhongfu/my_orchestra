// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <vector>

#include <MDArray/MDArray.h>

#include <Orchestra/Arc/Types.h>

#include <Orchestra/Calibration/Common/RawFile.h>

#include <Orchestra/Epi/MultibandCalibrationProcessor.h>

namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    namespace Epi
    {
        /* The main recon function */
        bool EpiMultiPhaseRecon(const GERecon::Legacy::PfilePointer& pfile);

        /**
         * typedef defining the type used to propagate dynamic coefficients for all slices
         * of all channels from phase to phase.
         * The dynamic phase correction plugin automatically updates these coefficients with
         * each call. Thus, the only action for the rehearsal pipeline is to allocate space
         * for these coefficients. The coefficient details are included below for informational
         * purpose only:
         * An individual set of dynamic coefficients is held in a TinyVector of size 4 and 
         * type float. This represents four coefficients:
         * linear phase slope along phase direction, global constant image phase, linear phase
         * coefficient delta along readout direction, constant phase coefficient delta along
         * readout direction.
         * A set of these four coefficients is held for each slice and channel in the scan. These
         * deltas are propagated from phase to phase. Thus, the coefficients are held for each
         * slice and channel in a vector of vectors.
         */
        typedef std::vector<std::vector<MDArray::TinyVector<float, 4> > > DynamicCoefficientsType;

        /**
         * Build the Multiband Calibration 
         * @param calibrationFile                        Raw calibration file 
         * @param slice                                  Acquired slice number
         * @param multibandCalibrationProcessor          Multiband Calibration Processor
         * @param flipDataInPhaseEncodeDirection         Flip data in phase encode direction
         * @return pointer to Arc calibration            Pointer to Arc calibration data
         */
        GERecon::Arc::CalibrationPtr CalibrationProcessorMultiband(const boost::shared_ptr<GERecon::Calibration::RawFile>& calibrationfile,int slice, GERecon::Epi::MultibandCalibrationProcessor& multibandCalibrationProcessor, const bool flipDataInPhaseEncodeDirection);
    }
}

