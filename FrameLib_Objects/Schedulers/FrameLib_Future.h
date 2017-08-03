
#ifndef FRAMELIB_FUTURE_H
#define FRAMELIB_FUTURE_H

#include "FrameLib_DSP.h"

class FrameLib_Future : public FrameLib_Scheduler
{
    // Parameter Enums and Info

    enum ParameterList { kTime, kUnits };
    enum Units { kSamples, kMS, kSeconds };

    struct ParameterInfo : public FrameLib_Parameters::Info { ParameterInfo(); };

public:
    
    // Constructor
    
    FrameLib_Future(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, void *owner);
    
    // Info
    
    std::string objectInfo(bool verbose);
    std::string inputInfo(unsigned long idx, bool verbose);
    std::string outputInfo(unsigned long idx, bool verbose);
    
private:

    // Schedule
    
    SchedulerInfo schedule(bool newFrame, bool noAdvance);
    
    // Data
    
    static ParameterInfo sParamInfo;
};

#endif
