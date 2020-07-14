
#ifndef FRAMELIB_YIN_H
#define FRAMELIB_YIN_H

#include "FrameLib_DSP.h"
#include "../../FrameLib_Dependencies/SpectralProcessor.hpp"
#include <numeric>
#include <iterator>
#include <algorithm>

class FrameLib_Yin final : public FrameLib_Processor
{
	using EdgeMode = spectral_processor<double, FrameLib_DSP::Allocator>::EdgeMode;
	enum ParameterList { kF0Min, kF0Max, kHarmoThresh };

	struct ParameterInfo : public FrameLib_Parameters::Info { ParameterInfo(); };
public:

    // Constructor

    FrameLib_Yin(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    // Process
    
    void process() override;
	void differenceFunction(const double *x, double * output, unsigned int N, double tauMax);

	spectral_processor<double, FrameLib_DSP::Allocator> mProcessor;
};

#endif

