
#include "FrameLib_Ticks.h"
#include "FrameLib_MaxClass.h"

extern "C" int C74_EXPORT main(void)
{
    FrameLib_MaxClass<FrameLib_Expand <FrameLib_Ticks> >::makeClass(CLASS_BOX, "fl.ticks~");
}