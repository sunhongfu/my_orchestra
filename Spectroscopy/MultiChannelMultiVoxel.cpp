// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <Hdf5/File.h>

#include <MDArray/MDArray.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include <Orchestra/Spectro/MCSI/CalibrationFile.h>
#include <Orchestra/Spectro/MCSI/Filters.h>
#include <Orchestra/Spectro/MCSI/LxControlSource.h>
#include <Orchestra/Spectro/MCSI/ProcessAndCombine.h>
#include <Orchestra/Spectro/MCSI/Utilities.h>

#include <Orchestra/Spectro/MultiVoxel/Utilities.h>

#include "MultiChannelMultiVoxel.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Spectro::MCSI;
using namespace GERecon::Legacy;

void GERecon::Spectro::MultiChannelMultiVoxelRecon(const std::string& pfileName) // parasoft-suppress  METRICS-22 "Procedural Rehearsal Code"
{  
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("MultiChannelMultiVoxelRehearsal");

    trace.ConsoleMsg("Running Multi-Channel Multi-Voxel Spectroscopy Reconstruction");

    // Setup a pfile object
    const boost::filesystem::path pFileDirectory = boost::filesystem::path(pfileName).parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "cppDicoms";
    boost::shared_ptr<Pfile> pfilePointer = Pfile::Create(pfileName);
    const boost::filesystem::path pfileFullFilePath(pfileName);

    // Create processing control and extract parameters
    const boost::shared_ptr<Control::ProcessingControl> processingControlPtr = pfilePointer->CreateOrchestraProcessingControl<GERecon::Spectro::MCSI::LxControlSource>();
    const int acquiredXRes = processingControlPtr->Value<int>("AcquiredXRes");
    const int acquiredYRes = processingControlPtr->Value<int>("AcquiredYRes");
    const int acquiredZRes = processingControlPtr->Value<int>("AcquiredZRes");
    const int numChannels = processingControlPtr->Value<int>("NumChannels");
    const int numPointsInFID = processingControlPtr->Value<int>("AcquiredFidLength");
    const bool chopX = processingControlPtr->Value<bool>("ChopX");
    const GERecon::Rotation rotationFlag = processingControlPtr->Value<GERecon::Rotation>("RotateFlag");
    const GERecon::Transpose transposeFlag = processingControlPtr->Value<GERecon::Transpose>("TransposeFlag");
    // Additional parameters which may be of interest, but are not directly utilized in this example pipeline are:
    //                 processingControlPtr->Value<int>("SpectralBandwidthHz");
    //                 processingControlPtr->Value<int>("TransformXRes"); // Spatial TransformXRes
    //                 processingControlPtr->Value<int>("TransformYRes"); // Spatial TransformYRes
    //                 processingControlPtr->Value<int>("NumSlicesToReconstruct");
    //                 processingControlPtr->Value<int>("SpectralTransformSize");
    //                 processingControlPtr->Value<int>("Temperature");
    //                 processingControlPtr->Value<float>("ScanCenterFrequencyMHz");
    //                 processingControlPtr->Value<float>("XFieldOfView");
    //                 processingControlPtr->Value<float>("YFieldOfView");
    //                 processingControlPtr->Value<float>("ZFieldOfView");
    //                 processingControlPtr->Value<float>("VoxelCenterRLDirRAS");
    //                 processingControlPtr->Value<float>("VoxelCenterAPDirRAS");
    //                 processingControlPtr->Value<float>("VoxelCenterSIDirRAS");
    //                 processingControlPtr->Value<float>("VoxelWidthRLDirRAS");
    //                 processingControlPtr->Value<float>("VoxelWidthAPDirRAS");
    //                 processingControlPtr->Value<float>("VoxelWidthSIDirRAS");
    
    // Create DICOM Series object to save images into
    const Legacy::DicomSeries dicomSeries(pfilePointer->DownloadData(), Legacy::DicomSeries::SeriesTimes100);

    boost::shared_ptr<GERecon::Spectro::MCSI::CalibrationFile> mcsiCalibration = boost::make_shared<GERecon::Spectro::MCSI::CalibrationFile>(*processingControlPtr);
    boost::filesystem::path calibrationPath = pFileDirectory / "McsiCalibration.h5";
    mcsiCalibration->Load(calibrationPath);

    GERecon::Spectro::MCSI::ProcessAndCombine processAndCombine(*processingControlPtr, mcsiCalibration);
    GERecon::Spectro::MCSI::Utilities utilities(*processingControlPtr);    
    GERecon::Spectro::MCSI::Filters filters(*processingControlPtr);
    const FloatVector spectralApodizationFilter = filters.SpectralApodizationFilter();

    GEHdf5::File hdf("multiVoxelRehearsal.h5", GEHdf5::File::Overwrite);

    for(int channel = 0; channel < numChannels; ++channel)
    {
        trace.ConsoleMsg("Processing Channel: %d of %d", (channel+1), numChannels);

        // Read in multi-slice, multi-voxel data from Pfile
        // Multi-Voxel spectroscopy data is stored in a Pfile on a 
        // channel by channel basis using a flat view index. That is, 
        // the slice and echo dimensions remain constant. The view 
        // index is incremented with each acquired FID. Thus, to access 
        // the first acquired FID, use slice = 0, echo = 0, view = 0. 
        // To access the next acquired FID, use 
        // slice = 0, echo = 0, view = 1, and so on.
        // Consider a hypothetical 3x4x3 multi-voxel scan
        // (acquiredXRes = frequency res on UI = 3), (acquiredYRes =
        // phase res on UI = 4), (acquiredZRes = num slices on UI =
        // 3)
        // 
        // Add equation in here: flatViewIndex = xIndex +
        // yIndex*(acquiredXRes) + slice*(acquiredXRes*acquiredYRes);
        // 
        //  Voxel Indices  ----------  Pfile Indices
        //  x=0, y=0, z=0  ----------  slice=0, echo=0, view=0
        //  x=1, y=0, z=0  ----------  slice=0, echo=0, view=1
        //  x=2, y=0, z=0  ----------  slice=0, echo=0, view=2
        //  x=0, y=1, z=0  ----------  slice=0, echo=0, view=3
        //  x=1, y=1, z=0  ----------  slice=0, echo=0, view=4
        //  x=2, y=1, z=0  ----------  slice=0, echo=0, view=5
        //  x=0, y=2, z=0  ----------  slice=0, echo=0, view=6
        //  ...
        //  x=0, y=3, z=2  ----------  slice=0, echo=0, view=33
        //  x=1, y=3, z=2  ----------  slice=0, echo=0, view=34
        //  x=2, y=3, z=2  ----------  slice=0, echo=0, view=35
        const int slice = 0;
        const int echo = 0;

        // Allocate space to hold a single channel's acquired data (must re-allocate each time through this loop to
        // handle non-square acquisition matrices being rotated/transposed below)
        Array<std::complex<float>, 4> singleChannelSortedFids(numPointsInFID, acquiredXRes, acquiredYRes, acquiredZRes);
        ComplexFloatMatrix singleChannelAllFids = pfilePointer->KSpaceData<float>(slice, echo, channel);
        for(int slice = 0; slice < acquiredZRes; ++slice)
        {
            for(int yIndex = 0; yIndex < acquiredYRes; ++yIndex)
            {
                for(int xIndex = 0; xIndex < acquiredXRes; ++xIndex)
                {
                    const int flatViewIndex = xIndex + yIndex*acquiredXRes + slice*(acquiredXRes*acquiredYRes);

                    singleChannelSortedFids(Range::all(),xIndex,yIndex,slice) = singleChannelAllFids(Range::all(), flatViewIndex);
                }
            }
        }

        // Chop along x-direction
        // chopX indicates to chop along x when it is true.
        // This is a negated version of the PSD RHTYPCHP bit which
        // indicates that rf chopping was turned on
        if(chopX)
        {
            singleChannelSortedFids(Range::all(),Range(1,toEnd,2),Range::all(),Range::all()) = singleChannelSortedFids(Range::all(),Range(1,toEnd,2),Range::all(),Range::all()) * -1.0f;
        }
        
        // Baseline correction
        // Inspect the end of each FID to determine an average baseline
        // offset. Subtract the average baseline offset from the raw
        // acquired data.
        Utilities::ApplyBaselineCorrection(singleChannelSortedFids);
        
        // Orient
        // Rotate/Transpose the data according to the slice orientation
        // flags in the raw header for the pfile that is currently loaded.
        Utilities::ApplyRotateTranspose(singleChannelSortedFids, rotationFlag, transposeFlag);
        
        // Zerofill
        // Zerofill data such that the spatial x/y dimensions have the 
        // same size and are equal to a power of 2. Example: An 8x10 
        // acquisition will be zerofilled to 16x16 spatial voxels. The FID 
        // spectral dimensions is also zerofilled to twice the acquired 
        // FID length.
        MDArray::Array<std::complex<float>, 4> zerofilledMatrix;
        utilities.Zerofill(zerofilledMatrix, singleChannelSortedFids);

        // Spectral transform 
        // Apply the spectral transform along the FID dimension of the
        // zerofilled data. A spectral filter is applied along
        // the FID dimension prior to transforming.
        const int xDimensionIndex = 1;
        const int yDimensionIndex = 2;
        const int zDimensionIndex = 3;
        for(int z = 0; z < zerofilledMatrix.extent(zDimensionIndex); ++z)
        {
            for(int y = 0; y < zerofilledMatrix.extent(yDimensionIndex); ++y)
            {
                for(int x = 0; x < zerofilledMatrix.extent(xDimensionIndex); ++x)
                {
                    // First apply the spectral apodization filter to each FID (no need to apply filter to zerofilled portion of FID)
                    ComplexFloatVector nonZerofilledPortionOfDataToTransform = zerofilledMatrix(Range(0,spectralApodizationFilter.extent(firstDim)-1), x, y, z);
                    nonZerofilledPortionOfDataToTransform *= spectralApodizationFilter;

                    // Next do the spectral transform on this FID
                    ComplexFloatVector dataToTransformVector = zerofilledMatrix(Range::all(), x, y, z);
                    utilities.SpectralTransform(dataToTransformVector);
                }
            }
        }
    
        // Spatial Transform
        // Apply the spatial transform along the x,y,z spatial dimensions.
        // Note that this function does not apply any spatial apodization filters.
        Utilities::SpatialTransform(zerofilledMatrix);

        // Accumulate data
        // Accumulate data using the GE multi-channel combine algorithm.
        // Note that this algorithm is specifically designed to be used 
        // with the product pulse sequence and processing steps above. This
        // algorithm will extract a subset of the spectrum to work with. The
        // algorithm is not compatible with 7.0T scans.
        // Part of the processing done in this function is to transform back
        // to kz and adjust the z-transform such that it yields slices at
        // slice locations with the same spacing as the localizer scan. This
        // is accomplished by adding a conversion factor to the exponent of
        // the Fourier transform. The conversion factor is stored in
        // rhuser48 and equals:
        //
        // rhuser48=(LocalizerSliceThickness + SliceSpacing) / (CSISliceThickness)
        //
        // Including the conversion factor above in the z-transform exponent
        // results in each reconstructed slice 'n' being spaced by 
        // (LocalizerSliceThickness + SliceSpacing). The number of output slices 
        // can be obtained from the pfile's SliceCount() function.
        processAndCombine.ProcessAndAccumulate(zerofilledMatrix, channel);
    }
    
    // Retrieve Accumulated Data
    // Retrieve the processed and accumulated data. The data returned here
    // will be a 256 point spectrum from ~ -0.4ppm to 4.3ppm. The ppm range
    // returned for plotting by GE multi-channel accumulation function is
    // not customizable at this time.
    const Array<std::complex<float>, 4> accumulatedData = processAndCombine.AccumulatedData();
    
    // Postage Stamp Generation    
    // Determine the number of reconstructed slices to make a postage stamp
    // for. The number of reconstructed slices is not equal to the number of
    // acquired slices because of the z-transform adjustment described
    // above. The slice dimension of the accumulated data is the fourth
    // dimension (SpectralDimension x X x Y x Slice).
    // Postage stamps are a 2D representation of spectral data for multiple
    // voxels of a single slice. A postage stamp always consists of a 16 x
    // 16 grid of squares. Each square contains one point of a 256 point
    // spectrum for each spatial voxel in the slice. Thus, for an 8x8 set of
    // spatial locations, each of the 16 squares has a size of 8x8. The
    // resulting image size is 16*8 x 16*8 = 128x128. For a 16x16 set of
    // spatial locations the resulting image size is 16*16 x 16*16 =
    // 256x256.
    const int numReconStructedSlices = accumulatedData.extent(fourthDim);
       
    FloatCube dataToGeneratePostageStampFrom(accumulatedData.extent(firstDim), accumulatedData.extent(secondDim), accumulatedData.extent(thirdDim));
    int imageNumber = 0;
    for(int slice = 0; slice < numReconStructedSlices; ++slice)
    {
        trace.ConsoleMsg("Generating Postage Stamps for Slice: %d of %d", (slice+1), numReconStructedSlices);
        ComplexFloatCube currentSlice = accumulatedData(Range::all(), Range::all(), Range::all(), slice);

        // Generate real and stamps
        FloatMatrix postageStamp;
        dataToGeneratePostageStampFrom = real(currentSlice);
        GERecon::Spectro::MultiVoxel::Utilities::GeneratePostageStamp(postageStamp, dataToGeneratePostageStampFrom);
        Clipper::Apply(postageStamp, GERecon::RealImage);

        // X x Y x [Real/Imag]
        ShortMatrix outputStamp(postageStamp.extent(firstDim), postageStamp.extent(secondDim));
        outputStamp = MDArray::cast<short>(postageStamp);

        const GERecon::SliceOrientation& sliceOrientation = pfilePointer->Orientation(slice);
        const GERecon::SliceCorners& sliceCorners = pfilePointer->Corners(slice);
        const GERecon::ImageCorners imageCorners(sliceCorners, sliceOrientation);
                
        std::stringstream realDicomFileName;
        realDicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";
        dicomSeries.SaveImage(realDicomFileName.str(), outputStamp, imageNumber, imageCorners);
        imageNumber = imageNumber + 1;

        dataToGeneratePostageStampFrom = imag(currentSlice);
        GERecon::Spectro::MultiVoxel::Utilities::GeneratePostageStamp(postageStamp, dataToGeneratePostageStampFrom);
        Clipper::Apply(postageStamp, GERecon::ImaginaryImage);
        outputStamp = MDArray::cast<short>(postageStamp);

        std::stringstream imagDicomFileName;
        imagDicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";
        dicomSeries.SaveImage(imagDicomFileName.str(), outputStamp, imageNumber, imageCorners);
        imageNumber = imageNumber + 1;
    }
}
