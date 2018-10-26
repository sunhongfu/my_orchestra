// Copyright 2018 General Electric Company. All rights reserved.
// GE Proprietary and Confidential Information. Only to be distributed with
// permission from GE. Resulting outputs are not for diagnostic purposes.

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <MDArray/MDArray.h>

#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ReconPaths.h>

#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>

#include <Orchestra/Spectro/MultiVoxel/Utilities.h>

#include "PerChannelMultiVoxel.h"

using namespace MDArray;
using namespace GERecon;
using namespace GERecon::Spectro::MultiVoxel;
using namespace GERecon::Legacy;

void GERecon::Spectro::PerChannelMultiVoxelRecon(const std::string& pfileName) // parasoft-suppress  METRICS-22 "Procedural Rehearsal Code"
{  
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out. In a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("PerChannelMultiVoxelRehearsal");

    trace.ConsoleMsg("Running Multi-Channel Multi-Voxel Spectroscopy Reconstruction");

    // Setup a pfile object
    const boost::filesystem::path pFileDirectory = boost::filesystem::path(pfileName).parent_path();
    const boost::filesystem::path dicomDirectory = pFileDirectory / "cppDicoms";
    boost::shared_ptr<Pfile> pfilePointer = Pfile::Create(pfileName);
    const boost::filesystem::path pfileFullFilePath(pfileName);

    // Create processing control and extract parameters
    const boost::shared_ptr<Control::ProcessingControl> processingControlPtr = pfilePointer->CreateOrchestraProcessingControl<GERecon::Spectro::MultiVoxel::LxControlSource>();
    const int acquiredXRes = processingControlPtr->Value<int>("AcquiredXRes");
    const int acquiredYRes = processingControlPtr->Value<int>("AcquiredYRes");
    const int acquiredZRes = processingControlPtr->Value<int>("AcquiredZRes");
    const int numChannels = processingControlPtr->Value<int>("NumChannels");
    const int numPointsInFID = processingControlPtr->Value<int>("AcquiredFidLength");
    const GERecon::Rotation rotationFlag = processingControlPtr->Value<GERecon::Rotation>("RotateFlag");
    const GERecon::Transpose transposeFlag = processingControlPtr->Value<GERecon::Transpose>("TransposeFlag");
    const NegateFinalSpectrumFlags negateFinalSpectrumFlag = processingControlPtr->ValueStrict<NegateFinalSpectrumFlags>("NegateFinalSpectrumFlag");
    const float zTransformConversionFactor = processingControlPtr->Value<float>("ZTransformLocationConversionFactor");
    const int numSlicesToReconstruct = processingControlPtr->Value<int>("NumSlicesToReconstruct");
    const int transformXRes = processingControlPtr->Value<int>("TransformXRes"); // Spatial TransformXRes
    const int transformYRes = processingControlPtr->Value<int>("TransformYRes"); // Spatial TransformYRes
    const float scalingFactor = processingControlPtr->Value<float>("ScalingFactor");
    const float temperature = processingControlPtr->Value<float>("Temperature");
    // Additional parameters which may be of interest, but are not directly utilized in this example pipeline are:
    //                 processingControlPtr->Value<int>("SpectralBandwidthHz");
    //                 processingControlPtr->Value<float>("ScanCenterFrequencyMHz");
    //                 processingControlPtr->Value<bool>("ChopX");
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

    // Allocate space to hold raw acquired data. This matrix is filled with acquisition data
    // each time through the channel loop below.
    Array<std::complex<float>, 4> singleChannelSortedFids(numPointsInFID, acquiredXRes, acquiredYRes, acquiredZRes);

    // Allocate space for the z-transform
    ComplexFloatMatrix transformWorkspace(numPointsInFID, acquiredZRes);
    MDArray::Array<std::complex<float>, 4> dataAfterZDirectionTransform(numPointsInFID, acquiredXRes, acquiredYRes, numSlicesToReconstruct);
    
    // Create a spatial apodization vector to be applied along the z dimension of our data
    GERecon::Spectro::MultiVoxel::Filters multiVoxelFilters(*processingControlPtr);
    FloatVector zTransformApodizationFilter(acquiredZRes);
    multiVoxelFilters.ZDirectionFermiWindow(zTransformApodizationFilter);

    // Create an xy spatial domain filter
    ComplexFloatMatrix xySpatialFilter(acquiredXRes, acquiredYRes);
    multiVoxelFilters.CombinedApodizationAndPhaseShiftFilterXYPlane(xySpatialFilter);

    // Create our spectral filter
    FloatVector spectralFilter(numPointsInFID);
    multiVoxelFilters.SpectralApodizationWindow(spectralFilter);

    //// Channel Loop
    // This rehearsal reconstructs each channel individually. The rehearsal does
    // not implement a channel combine and is used for single channel 1.5T
    // PROSE exams in product.
    int imageNumber = 0;
    for(int channel = 0; channel < numChannels; ++channel)
    {
        //// Read in multi-slice, multi-voxel data from pfile
        // Multi-Voxel spectroscopy data is stored in a pfile on a 
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
        // * Voxel Indices  ----------  Pfile Indices
        // * x=0, y=0, z=0  ----------  slice=0, echo=0, view=0
        // * x=1, y=0, z=0  ----------  slice=0, echo=0, view=1
        // * x=2, y=0, z=0  ----------  slice=0, echo=0, view=2
        // * x=0, y=1, z=0  ----------  slice=0, echo=0, view=3
        // * x=1, y=1, z=0  ----------  slice=0, echo=0, view=4
        // * x=2, y=1, z=0  ----------  slice=0, echo=0, view=5
        // * x=0, y=2, z=0  ----------  slice=0, echo=0, view=6
        // * ...
        // * x=0, y=3, z=2  ----------  slice=0, echo=0, view=33
        // * x=1, y=3, z=2  ----------  slice=0, echo=0, view=34
        // * x=2, y=3, z=2  ----------  slice=0, echo=0, view=35
        //     
        const int slice = 0;
        const int echo = 0;
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
        
        //// Z-direction Transform
        // The Multi-Voxel Spectroscopy Z direction transform yields slices at
        // slice locations with the same spacing as the localizer scan. This
        // is accomplished by adding a conversion factor to the exponent of
        // the Fourier transform. The conversion factor is stored in
        // rhuser48 and equals:
        //
        // $$rhuser48=\frac{(LocalizerSliceThickness + SliceSpacing)}{(CSISliceThickness)}$$
        //
        // Transform with conversion factor:
        //
        // $$f(n)=\sum\limits_{k=0}^{n-1}F(k)e^{-i\frac{2\pi}{N}kn(rhuser48)}$$
        //
        // Using the z-transform above results in each transform location
        // 'n' being spaced by (LocalizerSliceThickness + SliceSpacing). The
        // number of output slices can be obtained from the Pfile's
        // sliceCount.        
        if(acquiredZRes > 1)
        {
            for(int slice = 0; slice < numSlicesToReconstruct; ++slice)
            {
                const float locationIncrementToObtainFromTransform = slice - (numSlicesToReconstruct - 1.0f)/2.0f;
            
                for(int x = 0; x < singleChannelSortedFids.extent(secondDim); ++x)
                {
                    for(int y = 0; y < singleChannelSortedFids.extent(thirdDim); ++y)
                    {
                        // Apply z-direction apodization filter and format data into a 2D ComplexFloatMatrix
                        // with the FID dimension in the first dimension and the acquired Z dimension in the
                        // second dimension
                        for(int z = 0; z < acquiredZRes; ++z)
                        {
                            transformWorkspace(Range::all(), z) = singleChannelSortedFids(Range::all(), x, y, z) * zTransformApodizationFilter(z);
                        }

                        const ComplexFloatVector voxelAfterZTransformCurrentReconSlice = Utilities::ReformattedTransform(transformWorkspace, zTransformConversionFactor, locationIncrementToObtainFromTransform, GERecon::Spectro::MultiVoxel::ForwardReformattedTransform);
                        dataAfterZDirectionTransform(Range::all(),x,y,slice) = voxelAfterZTransformCurrentReconSlice;
                    }
                }
            }
        }
        else
        {
            // This is single slice data, no transform processing to do
            dataAfterZDirectionTransform = singleChannelSortedFids;
        }
        
        for(int reconstructedSlice = 0; reconstructedSlice < dataAfterZDirectionTransform.extent(fourthDim); ++reconstructedSlice)
        {
            //// Spectral and X/Y Spatial Direction Transforms
            // This function will apply
            // the spectral and XY spatial direction transform (the Z-direction
            // transform was applied above).
            // This function will also apply baseline correction to remove any DC
            // offset in the FIDs and subtract the phase of the first point in
            // each spectrum such that the first point in each spectrum has a
            // phase angle of zero degrees.
            ComplexFloatCube sliceDataToTransform = dataAfterZDirectionTransform(Range::all(), Range::all(), Range::all(), reconstructedSlice);
            ComplexFloatCube transformedSliceData = Utilities::Transform(sliceDataToTransform, xySpatialFilter, spectralFilter, transformXRes, transformYRes);
            
            //// Extract ROI
            // Extract the region of interest to include in postage stamps.
            // This function will interpolate the spectrum such that 256
            // points are extracted between ~4.3ppm and 0.4ppm. The extracted
            // region will be adjusted based on the temperature setting passed 
            // to the function.
            const int numPointsToPlot = 256; // postage stamp format always contains 256 points
            ComplexFloatCube extractedSpectra(numPointsToPlot, transformedSliceData.extent(secondDim), transformedSliceData.extent(thirdDim));
            Utilities::ExtractSpectralRegionOfInterest(extractedSpectra, transformedSliceData, temperature);
            extractedSpectra *= scalingFactor;
            
            //// Orient
            // Orient the spectrums according to the rotate and transpose values from processing control.
            // Also, if the final spectrum is to be negated, apply that operation. The negate operation
            // is dictated by the negateFinalSpectrumFlag which may have these three possible values:
            //    NeverNegate
            //    NegateIfMetaboliteRegionMeanIsNegative
            //    NegateIfLargestPeakIsNegative
            // If the input spatial X/Y dimensions are non-square this function will zerofill the data such
            // that a square matrix is returned. The square matrix can be given to the GeneratePostageStamp
            // function to generate postage stamp images.
            ComplexFloatCube squareOrientedData = Utilities::Orient(extractedSpectra, negateFinalSpectrumFlag, rotationFlag, transposeFlag);

            //// Postage Stamp Generation    
            // Postage stamps are a 2D representation of spectral data for multiple
            // voxels of a single slice. A postage stamp always consists of a 16 x
            // 16 grid of squares. Each square contains one point of a 256 point
            // spectrum for each spatial voxel in the slice. Thus, for an 8x8 set of
            // spatial locations, each of the 16 squares has a size of 8x8. The
            // resulting image size is 16*8 x 16*8 = 128x128. For a 16x16 set of
            // spatial locations the resulting image size is 16*16 x 16*16 =
            // 256x256.
            const GERecon::SliceOrientation& sliceOrientation = pfilePointer->Orientation(slice);
            const GERecon::SliceCorners& sliceCorners = pfilePointer->Corners(slice);
            const GERecon::ImageCorners imageCorners(sliceCorners, sliceOrientation);

            FloatMatrix postageStamp;
            Utilities::GeneratePostageStamp(postageStamp, real(squareOrientedData));
            Clipper::Apply(postageStamp, GERecon::RealImage);

            // X x Y x [Real/Imag]
            ShortMatrix outputStamp(postageStamp.extent(firstDim), postageStamp.extent(secondDim));
            outputStamp = MDArray::cast<short>(postageStamp);

            std::stringstream realDicomFileName;
            realDicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";
            dicomSeries.SaveImage(realDicomFileName.str(), outputStamp, imageNumber, imageCorners);
            imageNumber = imageNumber + 1;

            Utilities::GeneratePostageStamp(postageStamp, imag(squareOrientedData));
            Clipper::Apply(postageStamp, GERecon::ImaginaryImage);
            outputStamp = MDArray::cast<short>(postageStamp);

            std::stringstream imagDicomFileName;
            imagDicomFileName << dicomDirectory.string() << "/Image_" << std::setw(3) << std::setfill('0') << imageNumber << ".dcm";
            dicomSeries.SaveImage(imagDicomFileName.str(), outputStamp, imageNumber, imageCorners);
            imageNumber = imageNumber + 1;
        }
    }
}
