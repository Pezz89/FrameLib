
#include "FrameLib_iFFT.h"
#include "FrameLib_MaxClass.h"

extern "C" int C74_EXPORT main(void)
{
    FrameLib_MaxClass_Expand<FrameLib_iFFT>::makeClass("fl.ifft~");
}
