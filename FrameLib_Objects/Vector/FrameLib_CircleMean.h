
#ifndef FRAMELIBE_CCIRCLEMEAN_H
#define FRAMELIB_CIRCLEMEAN_H

#include "FrameLib_DSP.h"
#include "../../FrameLib_Dependencies/SpectralProcessor.hpp"
#include "../../FrameLib_Dependencies/Statistics.hpp"
#include <numeric>
#include <iterator>
#include <algorithm>
#include <math.h>

class FrameLib_CircleMean final : public FrameLib_Processor
{
	using EdgeMode = spectral_processor<double, FrameLib_DSP::Allocator>::EdgeMode;
	enum ParameterList { kRangeMax };

	struct ParameterInfo : public FrameLib_Parameters::Info { ParameterInfo(); };
public:

    // Constructor

    FrameLib_CircleMean(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    // Process
	static ParameterInfo sParamInfo;

    void process() override;
};

#endif