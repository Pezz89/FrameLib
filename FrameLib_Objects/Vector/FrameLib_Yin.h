
#ifndef FRAMELIB_YIN_H
#define FRAMELIB_YIN_H

#include "FrameLib_DSP.h"
#include "../../FrameLib_Dependencies/SpectralProcessor.hpp"
#include <numeric>
#include <iterator>
#include <algorithm>
#include <math.h>

class FrameLib_Yin final : public FrameLib_Processor
{
	using EdgeMode = spectral_processor<double, FrameLib_DSP::Allocator>::EdgeMode;
	enum ParameterList { kF0Min, kF0Max, kHarmoThresh };

	struct ParameterInfo : public FrameLib_Parameters::Info { ParameterInfo(); };
public:

    // Constructor

    FrameLib_Yin(FrameLib_Context context, const FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    // Process
    
    void process() override;
	void differenceFunction(const double *x, double * output, unsigned int N, unsigned int tauMax);
	void differenceFunction_slow(const double * x, double * output, unsigned int N, unsigned int tau_max);
	void cmndf(double *df, double *output, unsigned int N);
	//void getPitch(double *cmndf, double *df, double *f, double *harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th, double * tau_out);
	void getPitch(double *cmndf, double *df, double *f, double *harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th);

	spectral_processor<double, FrameLib_DSP::Allocator> mProcessor;

	static ParameterInfo sParamInfo;

	template<typename T>
	void parabolic_interp(T f, T harm, double * cmndf, unsigned int tau);
};

#endif

template<typename T>
inline void FrameLib_Yin::parabolic_interp(T f, T harm, double * cmndf, unsigned int tau)
{
	*f = 1 / 2. * (cmndf[tau - 1] - cmndf[tau + 1]) / (cmndf[tau - 1] - 2 * cmndf[tau] + cmndf[tau + 1]) + tau;
	*harm = cmndf[tau] - 1 / 4. * (cmndf[tau - 1] - cmndf[tau + 1]) * (*f - tau);
}
