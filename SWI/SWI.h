#include <MDArray/MDArray.h>
#include <Orchestra/Common/ReconTrace.h>

namespace GERecon
{

MDArray::ComplexFloatCube GaussianBlurVolume(MDArray::ComplexFloatCube thisImage,const float sigmaValue) ;
    


MDArray::ComplexFloatMatrix SWIComplexChannelCombine(MDArray::ComplexFloatCube channelData, MDArray::ComplexFloatCube coilSensitivities) ;
    
    MDArray::FloatCube GaussianBlurVolume(MDArray::FloatCube thisImage,float sigmaValue);
    
    void ScaleAndClipVolume(MDArray::FloatCube& inputVolume);
    
    MDArray::FloatCube ComputePhaseMaskConventional(
                                                    MDArray::Array<std::complex<float>,4> echoImages,const float kernelSize, const float maskVal, GERecon::Trace) ;

}
