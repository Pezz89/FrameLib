
#include "FrameLib_FIRPhase.h"
#include "FrameLib_MaxClass.h"

extern "C" int C74_EXPORT main(void)
{
    FrameLib_MaxClass_Expand<FrameLib_FIRPhase>::makeClass(CLASS_BOX, "fl.firphase~");
}