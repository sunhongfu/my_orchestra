// Copyright 2011 General Electric Company. All rights reserved.
#include <MDArray/Utils.h>
#include <MDArray/Fourier.h>

#include <Dicom/Network.h>
#include <Dicom/MR/Image.h>
#include <Hdf5/Snap.h>

#include <Orchestra/Arc/Arc.h>
#include <Orchestra/Asset/Calibration.h>
#include <Orchestra/Asset/CalibrationFile.h>
#include <Orchestra/Asset/Worker.h>


#include <Orchestra/Cartesian2D/KSpaceTransformer.h>

#include <Orchestra/Common/ImageCorners.h>

#include <Orchestra/Common/AnonymizationPolicy.h>
#include <Orchestra/Common/ReconPaths.h>
#include <Orchestra/Common/SliceOrientation.h>
#include <Orchestra/Common/SliceCorners.h>
#include <Orchestra/Common/ReconTrace.h>

#include <Orchestra/Control/ProcessingControl.h>
#include <Orchestra/Control/CommonTypes.h>

#include <Orchestra/Core/SumOfSquares.h>
#include <Orchestra/Core/Clipper.h>
#include <Orchestra/Core/Windows.h>
#include <Orchestra/Core/RotateTranspose.h>

#include <Orchestra/Gradwarp/GradwarpPlugin.h>

#include <Orchestra/Legacy/DicomSeries.h>
#include <Orchestra/Legacy/Pfile.h>
#include <Orchestra/Legacy/LxDownloadData.h>
#include <Orchestra/Legacy/LegacyImageDB.h>
#include <Orchestra/Legacy/LegacyRdbm.h>


#include <boost/filesystem.hpp>

#include "SWI.h"
#include "Rehearsal.h"
#include "CommandLine.h"

#if defined(__APPLE__)
#include <dispatch/dispatch.h>
#endif

// Include this to avoid having to type fully qualified names
using namespace GERecon;
using namespace MDArray;

/**
 * Function that performs a conventional phase-based susceptibility-weighted
 * reconstruction from multi-echo 3D gradient echo images
 *
 * Supports ARC as well as ASSET
 *
 * @author Kevin Koch, Medical College of Wisconsin, 2014
 * updated by Robin Karr, Medical College of Wisconsin, 2016
 *
 */
void GERecon::RunSWI()
{
    // Create a Trace buffer to log messages. In this environment the messages will go to standard
    // out, in a product environment messages will be logged to the appropriate file.
    GERecon::Trace trace("OxSWIRecon");
    
    trace.ConsoleMsg("Running Orchestra Phase-Based SWI Rehearsal");
    
    // Local configuration management, as read from configuration file input on command line
    const boost::optional<std::string> configFileName = CommandLine::ConfigFile();
    
    int anonFlag = 0 ;
    
    int inputVals[100] ;
    
    std::ifstream in (configFileName->c_str());
    
    if ( in.is_open() ) {
        std::string line;
        
        int pCounter = 0 ;
        while ( getline ( in, line ) ) {
            std::string::size_type i = line.find_first_not_of ( " \t\n\v" );
            
            if ( i != std::string::npos && line[i] == '#' )
                continue;
            
            inputVals[pCounter] = atoi(line.c_str()) ;
            pCounter++ ;
            
        }
    }
    
    anonFlag = inputVals[0] ;
    
    
    
    // Read Pfile and fill important information.
    const boost::filesystem::path pfileName = CommandLine::PfilePath();
    
    AnonymizationPolicy myAnonPolicy ;
    if(anonFlag = 0) myAnonPolicy = AnonymizationPolicy(AnonymizationPolicy::None) ;
    else myAnonPolicy = AnonymizationPolicy(AnonymizationPolicy::Default) ;
    
    
    //const Legacy::PfilePointer pfile = Legacy::Pfile::Create(pfileName, Legacy::Pfile::AllAvailableAcquisitions, AnonymizationPolicy(AnonymizationPolicy::None)) ;
    const Legacy::PfilePointer pfile = Legacy::Pfile::Create(pfileName, Legacy::Pfile::AllAvailableAcquisitions, myAnonPolicy) ;
    
    
    const Control::ProcessingControlPointer processingControl = pfile->CreateOrchestraProcessingControl();
    
    const unsigned int examNumber = processingControl->Value<unsigned int>("ExamNumber");
    const unsigned int seriesNumber = processingControl->Value<unsigned int>("SeriesNumber");
    const std::string patientID = processingControl->Value<std::string>("PatientID");
    const std::string protocol = processingControl->Value<std::string>("Protocol");
    const int imageXRes = processingControl->Value<int>("ImageXRes");
    const int imageYRes = processingControl->Value<int>("ImageYRes");
    const int reconXRes = processingControl->Value<int>("TransformXRes");
    const int reconYRes = processingControl->Value<int>("TransformYRes");
    const FloatVector channelWeights = processingControl->Value<FloatVector>("ChannelWeights");
    
    const int acqXRes = processingControl->Value<int>("AcquiredXRes");
    const int acqYRes = processingControl->Value<int>("AcquiredYRes");
    const int acqZRes = processingControl->Value<int>("AcquiredZRes");
    
    // Pulled directly from Pfile API and not ProcessingControl since the Pfile class knows exactly
    // how many slices and acquisitions it has loaded and can produce raw data for.
    const int numPhases = pfile->PhaseCount();
    const int numSlices = pfile->SliceCount();
    int numEchoes = pfile->EchoCount();
    int numChannels = pfile->ChannelCount();
    
    
    // Determine scaling factor for 3D transform
    const float scale3d = processingControl->Value<float>("3DScalingFactor");
    const float winApodization3d = processingControl->Value<float>("3DWinApodization");
    const float winQ3d = processingControl->Value<float>("3DWinQ");
    const int winType3d = processingControl->Value<short>("3DWinType");
    
    FloatVector zFilter(numSlices);
    Windows::Apodization3D(zFilter, winType3d, winApodization3d, winQ3d);
    zFilter /= scale3d;
    
    const Legacy::LxDownloadData& dd = *pfile->DownloadData();
    const Legacy::RdbHeaderRec& myHeader = dd.RawHeader();

    float echoSpacing ;
    echoSpacing = myHeader.rdb_hdr_te - myHeader.rdb_hdr_te2 ;
    trace.ConsoleMsg("Echo spacing as determined from header: %f ",echoSpacing);
    
    // setup for Asset
    boost::shared_ptr<Asset::CalibrationFile> assetCalibration;
    const bool isAssetScan = processingControl->Value<bool>("Asset");
    if(isAssetScan)
    {
        // If this is an asset scan then an asset calibration file must exist next to the pfile and have a name of: AssetCalibration.h5
        const unsigned int examNumber = processingControl->Value<unsigned int>("ExamNumber");
        const unsigned int coilID = processingControl->Value<unsigned int>("CoilConfigUID");
        assetCalibration = boost::make_shared<Asset::CalibrationFile>(examNumber, coilID, Asset::Calibration::Regular, *processingControl);
        if( !assetCalibration->UseForReadingIfContentsMatch(GERecon::Path::InputAppData()/"AssetCalibration.h5", GERecon::Path::InputAppData()/"AssetCalibration.h5") )
        {
            throw GERecon::Exception(__SOURCE__, "Invalid Asset Calibration file!");
        }
        
        assetCalibration->ReadyForReading();
        
        trace.ConsoleMsg("Performing 1D ASSET -- Calibratin Loaded Successfully");
        
        
    }
    
    trace.ConsoleMsg("Parameters Loaded, Filters Contructed");
    
    // a convenience for MDArray operations
    const Range all = Range::all();
    
    // Create a transformer object to complete the K-Space to image-space transformation.
    // This will apply a Fermi/phase-shift filter, scale the data, and perform a 2D IFFT.
    Cartesian2D::KSpaceTransformer transformer(*processingControl);
    
    // Create channel combiner object that will do the channel combining work in channel loop.
    
    // Gradwarp Plugin
    GradwarpPlugin gradwarp(*processingControl, GERecon::TwoDGradwarp, GERecon::XRMBGradient);
    
    // allocate storage for all of the complex echo images
    MDArray::Array<std::complex<float>,4> echoImages(reconXRes,reconYRes,numSlices,numEchoes) ;
    
    // construct phase chopper for spatial phase processing
    FloatMatrix altGrid(reconXRes,reconYRes) ;
    for(int thisX=0; thisX<reconXRes;thisX++){
        for(int thisY=0; thisY<reconYRes;thisY++){
            altGrid(thisX,thisY) = pow(-1.0,thisX)*pow(-1.0,thisY+1);
        }
    }
    
    // Log Message
    trace.ConsoleMsg("Reconstructing Exam: %d, Series: %d", examNumber, seriesNumber);
    
    // setup arc parameters.  i don't use the kacq file -- but just define the strides in rhuser variable in my psd
    int arcPhStride = 1 ;
    int arcSlStride = 1 ;
    if( processingControl->Value<bool>("ArcScan") ){
        arcPhStride = (int)processingControl->Value<float>("UserValue12");
        arcSlStride = (int)processingControl->Value<float>("UserValue13");
        
        trace.ConsoleMsg("Running ARC. Strides [%d x %d]",arcPhStride, arcSlStride);
    }
    
    
    // ok -- process the echoes
    for(int phase = 0; phase < numPhases; ++phase)
    {
        // For Complex Coil Combiner
        Array<std::complex<float>,4> coilSensitivities ;
        if(!isAssetScan) coilSensitivities.resize(reconXRes,reconYRes,numSlices,numChannels) ;
        
        for(int currentEcho = 0; currentEcho < numEchoes; ++currentEcho)
        {
            trace.ConsoleMsg("Processing Echo [%d of %d]",currentEcho+1, numEchoes);
            
            MDArray::Array<std::complex<float>,4> kSpaceEcho(acqXRes,acqYRes,acqZRes,numChannels) ;
            
            // Get KSpace [kx,ky,kz,channel] data for this bin
            for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
            {
             
		   for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel){
                    
                    ComplexFloatMatrix kSpace = pfile->KSpaceData<float>(currentSlice, currentEcho, currentChannel, phase);
                    // Extract a single channel of data from the pfile
                    kSpaceEcho(all,all, currentSlice, currentChannel) = kSpace ;
                }
            }
            
            // If Arc is enabled, run and transform/filter in Z
            if( processingControl->Value<bool>("ArcScan") )
            {
                trace.ConsoleMsg("Running ARC");
                
                
                // ARC Processing
                Arc::Process(kSpaceEcho,true);
                
                // Filter across Z
                for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                {
                    for(int y = 0; y < acqYRes; ++y)
                    {
                        for(int x = 0; x < acqXRes; ++x)
                        {
                            kSpaceEcho(x, y, all, currentChannel) *= zFilter;
                        }
                    }
                }
                
                // TODO: handle slice zipping
                
                // Transform across Z
                Fourier::Ifft(kSpaceEcho, thirdDim);
            }
            
            
            
            trace.ConsoleMsg("In-Plane Processing of Volumetric Data");
#if defined(__APPLE__)
            dispatch_apply(numSlices, dispatch_get_global_queue(0, 0), ^(size_t currentSlice){
#else
#pragma omp parallel for // multi threaded distribution here
                for(int currentSlice = 0; currentSlice < numSlices; ++currentSlice)
                {
#endif
                    
                    
                    
                    // Get information for current slice
                    const SliceCorners& sliceCorners = pfile->Corners(currentSlice);
                    
                    //std::cout << "Corners: " << sliceCorners << std::endl ;
                    
                    boost::shared_ptr<Asset::Worker> asset;
                    if(isAssetScan){
                        asset.reset(new Asset::Worker(assetCalibration, *processingControl));
                    }
                    
                    // Create channel combiner object that will do the channel combining work in channel loop.
                    
                    Cartesian2D::KSpaceTransformer transformer(*processingControl);
                    
                    // Zero out channel combiner buffer for the next set of channels.
                    
                    ComplexFloatCube channelData(reconXRes, reconYRes, numChannels);
                    
                    // Storage for transformed image data, thus the image sizes.
                    ComplexFloatMatrix imageData(imageXRes, imageYRes);
                    
                    for(int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                    {
                        
                        // Extract a single channel of data from the pfile
                        ComplexFloatMatrix kSpace = kSpaceEcho(all,all,currentSlice,currentChannel);
                        
                        // Tranform to image space. Data will be zipped from acquired size.
                        transformer.Apply(imageData, kSpace);
                        
                        
                        channelData(all,all,currentChannel) = imageData ;
                        
                    }
                    
                    ComplexFloatMatrix accumulatedChannels(reconXRes, reconYRes);
                    
                    // run asset if necessary
                    if(isAssetScan){
                        asset->Unalias(accumulatedChannels, channelData, currentSlice, sliceCorners);
                        // Chop Correct to smooth phase computation
                        accumulatedChannels *= altGrid ;
                    }
                    else{
                        
                        // my complex coil combination routine
                        ComplexFloatMatrix thisPlane(reconXRes, reconYRes) ;
                        for (int currentChannel = 0; currentChannel < numChannels; ++currentChannel)
                        {
                            thisPlane = channelData(Range::all(),Range::all(),currentChannel);
                            if(channelWeights(currentChannel) == 0) thisPlane = 0.0;
                            else thisPlane /= channelWeights(currentChannel);
                            
                            // Chop Correct to smooth phase computation
                            channelData(all,all,currentChannel) = thisPlane*altGrid;
                        }
                        
                        if(currentEcho==0){ // compute coil sensitivity maps using first echo
                            coilSensitivities(all,all,currentSlice,all)  = GaussianBlurVolume(channelData,10) ;
                        }
                        accumulatedChannels = SWIComplexChannelCombine( // do the combination here
                                                                       channelData,coilSensitivities(all,all,currentSlice,all)) ;
                        
                    }
                    
                    
                    echoImages(all,all,currentSlice,currentEcho) = accumulatedChannels ;
                    
                    
                }
#if defined(__APPLE__)
                );
#endif
                
            }
                           
    }
                           
   // output the echo images
   trace.ConsoleMsg("Outputing Echo Data as DICOM Images");
   //debugSnap.Write(echoSpacing,"echoSpacing") ;
   //debugSnap.Write(echoImages,"echoImages") ;
   
   // Pull real and imaginary components
   MDArray::Array<float,4> realImages(reconXRes,reconYRes,numSlices,numEchoes) ;
   MDArray::Array<float,4> imagImages(reconXRes,reconYRes,numSlices,numEchoes) ;
   
   
   // Create storage for final mag image represented as a float matrix
   FloatMatrix realImage(reconXRes, reconYRes);
   FloatMatrix imagImage(reconXRes, reconYRes);
   
   realImages = real(echoImages) ;
   imagImages = imag(echoImages) ;
   
   // now scale them down to short values (lossy compression here !!)
   FloatVector maxVals(2) ;
   maxVals(0) = max(abs(realImages)) ;
   maxVals(1) = max(abs(imagImages)) ;
   
   float totalMax = max(maxVals) ;
   
   // Equivalently scale real and complex across all echoes
   realImages *= 32767.0/totalMax ;
   imagImages *= 32767.0/totalMax ;
   
   
   // Create DICOM Series object to save images into
   const Legacy::DicomSeries realDicomSeries(501,pfile);
   const Legacy::DicomSeries imagDicomSeries(502,pfile);
   
   boost::filesystem::create_directory("realDicom") ;
   boost::filesystem::create_directory("imagDicom") ;
   
   for(int thisEcho =0; thisEcho < numEchoes; thisEcho++)
   {
       for(int currentSlice = 0; currentSlice < numSlices; currentSlice++)
       {
           // Get information for current slice
           const SliceOrientation& sliceOrientation = pfile->Orientation(currentSlice);
           const SliceCorners& sliceCorners = pfile->Corners(currentSlice);
           
           realImage = realImages(all,all,currentSlice,thisEcho) ;
           imagImage = imagImages(all,all,currentSlice,thisEcho) ;
           
           // Perform Gradwarp
           gradwarp.Run(realImage, sliceCorners, currentSlice);
           gradwarp.Run(imagImage, sliceCorners, currentSlice);
           
           
           FloatMatrix rotatedRealImage = RotateTranspose::Apply<float>(realImage,
                                                                        sliceOrientation.RotationType(), sliceOrientation.
                                                                        TransposeType());
           
           FloatMatrix rotatedImagImage = RotateTranspose::Apply<float>(imagImage,
                                                                        sliceOrientation.RotationType(), sliceOrientation.TransposeType());
           
           // Create storage for final image represented as a short matrix - what is sent to host.
           ShortMatrix finalRealImage(realImages.extent(0),realImages.extent(1));
           ShortMatrix finalImagImage(imagImages.extent(0),imagImages.extent(1));
           
           
           // Cast final image from float to short
           finalRealImage = MDArray::cast<short>(rotatedRealImage);
           finalImagImage = MDArray::cast<short>(rotatedImagImage);
           
           
           // Rotate/Transpose image and corner points accordingly
           const ImageCorners imageCorners(sliceCorners, sliceOrientation);
           
           const int imageNumber = currentSlice; // TODO: Generate!
           
           // Save Real Image
           std::ostringstream strm2;
           strm2 << "realDicom/RealSlice" << currentSlice << "Echo" << thisEcho << ".dcm";
           
           
           const GEDicom::MR::ImagePointer dicom2 = realDicomSeries.NewImage(finalRealImage, imageNumber, imageCorners);
           // Add custom annotation if specified on command line.
           // Note, boost::optional<> types only get inserted if set.
           dicom2->Insert<GEDicom::LongString>(0x0008, 0x103E, "Real Image"); // Series description, showed in image browser
           dicom2->Insert<GEDicom::LongString>(0x0008, 0x1090, "Real Image"); // Manufacturers model name, showed in annotation
           dicom2->Insert<GEDicom::IntegerString>(0x0018, 0x0086, thisEcho+1); // Record the echo number
           dicom2->Insert<GEDicom::LongString>(0x0028, 0x0107, echoSpacing); // Record the echo spacing in the "LargestPixelValue" DICOM Tag
           
           dicom2->Save(strm2.str());
           
           // Save Imaginary Image
           std::ostringstream strm3;
           strm3 << "imagDicom/ImagSlice" << currentSlice << "Echo" << thisEcho << ".dcm";
           
           
           const GEDicom::MR::ImagePointer dicom3 = imagDicomSeries.NewImage(finalImagImage, imageNumber, imageCorners);
           // Add custom annotation if specified on command line.
           // Note, boost::optional<> types only get inserted if set.
           dicom3->Insert<GEDicom::LongString>(0x0008, 0x103E, "Imag Image"); // Series description, showed in image browser
           dicom3->Insert<GEDicom::LongString>(0x0008, 0x1090, "Imag Image"); // Manufacturers model name, showed in annotation
           dicom3->Insert<GEDicom::IntegerString>(0x0018, 0x0086, thisEcho+1); // Record the echo number
           dicom3->Insert<GEDicom::LongString>(0x0028, 0x0107, echoSpacing); // Record the echo spacing in the "LargestPixelValue" DICOM Tag
           
           dicom3->Save(strm3.str());
       }
   }
}

                           
ComplexFloatCube GERecon::GaussianBlurVolume(ComplexFloatCube thisImage,float sigmaValue){
   int xNum = thisImage.extent(0) ;
   int yNum = thisImage.extent(1) ;
   int zNum = thisImage.extent(2) ;
   
   ComplexFloatMatrix GaussianKernel(xNum,yNum) ;
   
   for(int thisX=0; thisX<xNum; thisX++){
       for(int thisY=0; thisY<yNum; thisY++){
           GaussianKernel(thisX,thisY) =
           exp(-1.0*((float)((thisX-xNum/2)*(thisX-xNum/2))+
                     (float)((thisY-yNum/2)*(thisY-yNum/2)))/(2*sigmaValue*sigmaValue)) ;
       }
   }
   
   const float invSum = (1.0f / static_cast<float>(sum(real(GaussianKernel))));
   GaussianKernel *= invSum;
   
   // Move filter to k-space
   Fourier::Fft(GaussianKernel);
   
   ComplexFloatCube blurredOutput(xNum,yNum,zNum) ;
   
   ComplexFloatMatrix thisPlane(xNum,yNum);
   
   float scale = 1.0/(xNum*yNum) ;
   
   for(int thisZ=0; thisZ< zNum; thisZ++){
       
       thisPlane = thisImage(Range::all(),Range::all(),thisZ) ;
       
       // Move to k-space
       Fourier::Fft(thisPlane);
       
       
       // Convolve
       thisPlane *= GaussianKernel ;
       
       // Back to image space
       Fourier::Ifft(thisPlane) ;
       Fourier::IfftShift(thisPlane) ;
       
       // store in output cube
       blurredOutput(Range::all(),Range::all(),thisZ) = thisPlane*scale ;
       
       
   }
   
   return blurredOutput ;
}

FloatCube GERecon::GaussianBlurVolume(FloatCube thisImage,float sigmaValue){
   int xNum = thisImage.extent(0) ;
   int yNum = thisImage.extent(1) ;
   int zNum = thisImage.extent(2) ;
   
   ComplexFloatMatrix GaussianKernel(xNum,yNum) ;
   
   for(int thisX=0; thisX<xNum; thisX++){
       for(int thisY=0; thisY<yNum; thisY++){
           GaussianKernel(thisX,thisY) =
           exp(-1.0*((float)((thisX-xNum/2)*(thisX-xNum/2))+
                     (float)((thisY-yNum/2)*(thisY-yNum/2)))/(2*sigmaValue*sigmaValue)) ;
       }
   }
   
   const float invSum = (1.0f / static_cast<float>(sum(real(GaussianKernel))));
   GaussianKernel *= invSum;
   
   // Move filter to k-space
   Fourier::Fft(GaussianKernel);
   
   FloatCube blurredOutput(xNum,yNum,zNum) ;
   
   ComplexFloatMatrix thisPlane(xNum,yNum);
   
   float scale = 1.0/(xNum*yNum) ;
   
   for(int thisZ=0; thisZ< zNum; thisZ++){
       
       thisPlane = thisImage(Range::all(),Range::all(),thisZ) ;
       
       
       // Move to k-space
       Fourier::Fft(thisPlane);
       
       
       // Convolve
       thisPlane *= GaussianKernel ;
       
       
       // Back to image space
       Fourier::Ifft(thisPlane) ;
       Fourier::IfftShift(thisPlane) ;
       
       // store in output cube
       blurredOutput(Range::all(),Range::all(),thisZ) = real(thisPlane)*scale ;
   }
   
   return blurredOutput ;
}



ComplexFloatMatrix GERecon::SWIComplexChannelCombine(ComplexFloatCube channelData, ComplexFloatCube coilSensitivities){
   
   int xNum = coilSensitivities.extent(0) ; 
   int yNum = coilSensitivities.extent(1) ; 
   int channelNum = coilSensitivities.extent(2) ; 
   
   // Apply coil noise weighting on all the source images
   
   ComplexFloatMatrix outputPlane(xNum,yNum) ; 
   outputPlane = 0.0 ;
   
   // compute normalizer from sensitivity map
   FloatMatrix coilNorm(xNum, yNum);
   coilNorm = 0.0f;
   for (int currentChannel = 0; currentChannel < channelNum; ++currentChannel)
   {
       coilNorm += ( real(coilSensitivities(Range::all(), Range::all(),currentChannel))*
                    real(coilSensitivities(Range::all(), Range::all(),currentChannel)) +
                    imag(coilSensitivities(Range::all(), Range::all(),currentChannel))*
                    imag(coilSensitivities(Range::all(), Range::all(),currentChannel)) );
       
       
       outputPlane += channelData(Range::all(),Range::all(),currentChannel)*conj(coilSensitivities(Range::all(),Range::all(),currentChannel)) ;
       
   }
   
   coilNorm = sqrt(coilNorm);
   
   outputPlane = outputPlane/coilNorm ;
   
   return outputPlane ;
   
}

void GERecon::ScaleAndClipVolume(FloatCube& inputVolume){
   float maxVal ;
   maxVal = max(inputVolume) ;
   
   inputVolume = 32767.0*inputVolume/maxVal ;
   
}
