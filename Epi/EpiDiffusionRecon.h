// Copyright 2018 General Electric Company. All rights reserved.

#pragma once

#include <Dicom/MR/Types.h>

#include <MDArray/MDArray.h>

#include <Orchestra/Common/ImageType.h>

namespace GERecon
{
    namespace Legacy
    {
        class Pfile;
        typedef boost::shared_ptr<Pfile> PfilePointer;
    }

    namespace Epi
    {
        /**
         * Reconstruct the given EPI Diffusion Pfile. This recon looks for all secondary
         * inputs (ref.dat, vrgf.dat, AssetCalibration.h5) in the Pfile directory.
         * This recon reconstructs all acquired images (opening subsequent Pfiles based
         * on the first Pfile's name), but does not reconstruct combined diffusion images.
         * 
         * @param pfile
         */
        bool EpiDiffusionRecon(const GERecon::Legacy::PfilePointer& pfile);

        /**
         * Fill diffusion specific fields in the given DICOM image based on the provided
         * parameters.
         *
         * @param dicomImage
         * @param diffusionImageType
         * @param bValue
         * @param currentSliceIndex
         * @param numBValues
         */
        void FillDiffusionDicomFields(const GEDicom::MR::ImagePointer& dicomImage, const GERecon::DiffusionImageType diffusionImageType, const int bValue, const int currentSliceIndex, const int numBValues);

        /**
         * Fill tensor vector DICOM fields based on the given tensor vector components.
         *
         * @param dicomImage
         * @param tensorVectorComponents
         */
        void FillTensorVectorDicomFields(const GEDicom::MR::ImagePointer& dicomImage, const MDArray::FloatTriple& tensorVectorComponents);

    }
}

