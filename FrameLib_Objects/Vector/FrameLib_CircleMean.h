
#ifndef FRAMELIB_CIRCLEMEAN_H
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

    FrameLib_CircleMean(FrameLib_Context context, const FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    // Process
	static ParameterInfo sParamInfo;
	double to_positive_angle(double angle)
	{
		angle = fmod(angle, 360.0);
		if (angle < 0.0) angle += 360.0;
		return angle;
	}
	double last_avg;

    void process() override;
	template<typename T> T sign(T t) { return t > T(0) ? T(1) : T(-1); };
	template<typename T> T py_mod(T n, T a) {
		T r = fmod(n, a);
		if (r * sign(a) < T(0)) r += a;
		return r;
	};
	template<typename T> T distance(double len_my_list, T idx_1, T idx_2) {
		T ret = 0.0;
		double distances[2] = { py_mod(idx_1 - idx_2, len_my_list), py_mod(idx_2 - idx_1, len_my_list) };
		int selection = statMinPosition(distances, 2);	
		if (selection == 1) {
			ret = -distances[selection];
		}
		else {
			ret = distances[selection];
		}
		return ret;
	};
};

#endif