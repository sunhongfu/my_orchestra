// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <boost/filesystem/path.hpp>

#include <boost/shared_ptr.hpp>

#include <Hdf5/Snap.h>

#include <MDArray/MDArray.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Epi/MultibandCalibrationProcessor.h>

namespace GERecon
{
    class SliceInfoTable;
    class GradwarpPlugin;
    class ScanArchive;
    typedef boost::shared_ptr<ScanArchive> ScanArchivePointer;
    class RampSamplingPlugin;

    namespace Legacy
    {
        class DicomSeries;
    }

    namespace Calibration
    {
        class RawFile;
    }

    namespace Epi
    {
        class MultibandCalibrationProcessor;
        class NearestNeighborStaticPhaseCorrectionPlugin;
        class PhaseCorrectionReferenceFile;
        class RowFlipPlugin;
        class SelfNavDynamicPhaseCorrectionPlugin;

        namespace Diffusion
        {
            class DynamicPhaseCorrectionManager;
        }

        /**
         * Reconstruct the given EPI Diffusion ScanArchive.
         *
         * @param scanArchive
         */
        void EpiDiffusionScanArchiveRecon(const GERecon::ScanArchivePointer& scanArchive);

        /**
         * Decodes all control packets in scan archive and stores the T2 and B-value volumes
         *
         * @param t2Volume
         * @param bValueVolume
         * @param scanArchive
         * @param processingControl
         */
        void StoreVolumes(MDArray::Array<float,5>& t2Volume,
                          MDArray::Array<float,6>& bValueVolume, 
                          const GERecon::ScanArchivePointer& scanArchive,
                          const Control::ProcessingControlPointer& processingControl);

        /**
         * DiffusionCommon recon. Common reconstruction for both T2 volumes and b-value volume
         *
         * @param imageData
         * @param unaliasedGeometricSliceNumbers
         * @param allRawData
         * @param rowFlipPlugin
         * @param phaseCorrectionReferenceFile
         * @param staticPhaseCorrectionPlugin
         * @param dynamicPhaseCorrectionPlugin
         * @param dynamicPhaseCorrectionManager
         * @param slicePhaseIndexCounters
         * @param rampSamplingPlugin
         * @param processingControl
         * @param nexIndex
         * @param slice
         * @param phase
         * @param isRpgVolume
         */
        void DiffusionCommonRecon(MDArray::FloatCube& imageData, 
                                  std::vector<int>& unaliasedGeometricSliceNumbers,
                                  std::vector<unsigned int>& slicePhaseIndexCounters,
                                  const MDArray::ComplexFloatCube& allRawData,
                                  Epi::RowFlipPlugin& rowFlipPlugin,
                                  const boost::shared_ptr<Epi::PhaseCorrectionReferenceFile>& phaseCorrectionReferenceFile,
                                  const boost::shared_ptr<Epi::NearestNeighborStaticPhaseCorrectionPlugin>& staticPhaseCorrectionPlugin,
                                  const boost::shared_ptr<Epi::SelfNavDynamicPhaseCorrectionPlugin>& dynamicPhaseCorrectionPlugin,
                                  const boost::shared_ptr<Epi::Diffusion::DynamicPhaseCorrectionManager>& dynamicPhaseCorrectionManager,
                                  const boost::shared_ptr<RampSamplingPlugin>& rampSamplingPlugin,
                                  const Control::ProcessingControlPointer& processingControl,
                                  const GEHdf5::Snap& debugFile,
                                  const int nexIndex,
                                  const int slice,
                                  const int phase,
                                  const bool isRpgVolume = false);

        /** 
         * T2 Segment recon
         *
         * Handles Nex'ing, HOEC, distortion correction and image sending for T2 segment of Diffusion
         *
         * @param displacementMap          Distortion correction field map
         * @param referenceVolumeForMoCo   referenceImages for Motion correction
         * @param t2Volume                 T2 Volume
         * @param processingControl
         * @param dicomDirectory
         * @param dicomSeries
         */
        void T2SegmentRecon(MDArray::FloatCube& displacementMap,
                            MDArray::FloatCube& referenceVolumeForMoCo,
                            const MDArray::Array<float,5>& t2Volume, 
                            const Control::ProcessingControlPointer& processingControl,
                            const boost::filesystem::path& dicomDirectory,
                            const Legacy::DicomSeries& dicomSeries);

        /**
         * Diffusion segment recon
         *
         * Handles Nex'ing, HOEC, distortion correction and image sending for diffusion segment of Diffusion
         *
         * @param t2Volume                 Diffusion Volume
         * @param displacementMap          Distortion correction field map
         * @param referenceVolumeForMoCo   referenceImages for Motion correction
         * @param processingControl
         * @param dicomDirectory
         * @param dicomSeries
         */
        void DiffusionSegmentRecon(const MDArray::Array<float,6>& diffusionVolume,
                                   const MDArray::FloatCube& displacementMap,
                                   const MDArray::FloatCube& referenceVolumeForMoCo,
                                   const Control::ProcessingControlPointer& processingControl,
                                   const boost::filesystem::path& dicomDirectory,
                                   const Legacy::DicomSeries& dicomSeries);


        /**
         * Discards kissoff views in transform.
         *
         * @param imageData
         * @param transformKissoffViews
         */
        void DiscardKissoffViews(const MDArray::ComplexFloatMatrix& imageData,
                                 const int transformKissoffViews);

        /**
         * Multiband calibration processing done using Arc calibration manager
         *
         * If this is a Rpg volume, calibration will not be stored using the calibration manager. Arc::Calibrate will be used instead.
         *
         * @param calibrationPtr
         * @param multibandCalibrationProcessor
         * @param calibrationFile
         * @param logicalSliceNumber
         * @param storeCalibration
         * @param flipDataInPhaseEncodeDirection
         */
        void MultibandCalibration(GERecon::Arc::CalibrationPtr& calibrationPtr,
                                   MultibandCalibrationProcessor& multibandCalibrationProcessor,
                                   const boost::shared_ptr<GERecon::Calibration::RawFile>& calibrationFile,
                                   const int logicalSliceNumber,
                                   const bool storeCalibration,
                                   const bool flipDataInPhaseEncodeDirection);

        /**
         * Get packet information for EpiOpcode
         *
         * Diffusion hyperscan packet (opcode = DiffusionHyperScanOpcode) has all the information that
         * recon needes to sort data into T2 images, BValues, Diffusion Directions, etc...).
         * The EpiOpcode does not have fields to identify the various diffusion passes. Thus, manual
         * bookkeeping implemented in this function is used to identify which diffusion pass corresponds
         * with the current control packet when EpiOpcode control packets are received.
         *
         * @param frameType     1: T2, 2: Diffusion
         * @param t2Index
         * @param bValueIndex
         * @param diffDirIndex
         * @param packetCount   count for EpiOpcode packets
         * @param numSlices
         * @param totalNumT2
         * @param numBValues
         * @param numDiffusionDirections
         * @param diffusionNexTable
         * @return false if received packets are more than expected, else true
         */
        bool GetPacketInfo(int& frameType,
                           int& t2Index,
                           int& bValueIndex,
                           int& diffDirIndex,
                           const int packetCount,
                           const int numSlices,
                           const int totalNumRef,
                           const int totalNumT2,
                           const int numBValues,
                           const int numDiffusionDirections,
                           const std::vector<int>& diffusionNexTable);
    }
}
