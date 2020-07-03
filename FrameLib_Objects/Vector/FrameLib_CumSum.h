
#ifndef FRAMELIB_CUMSUM_H
#define FRAMELIB_CUMSUM_H

#include "FrameLib_DSP.h"

class FrameLib_CumSum final : public FrameLib_Processor
{
public:

    // Constructor

    FrameLib_CumSum(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    // Process
    
    void process() override;
};

#endif

