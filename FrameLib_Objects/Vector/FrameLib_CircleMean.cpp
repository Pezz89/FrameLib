#include "FrameLib_CircleMean.h"

/*
Implementation of the YIN fundamental frequency (f0) estimation algorithm, as described in:
[1] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.                                                                                                                                                           3 of America, 111(4), 1917-1930.⏎
*/

// Constructor

FrameLib_CircleMean::FrameLib_CircleMean(FrameLib_Context context, const FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, &sParamInfo, 1, 2)
{
	mParameters.addDouble(kRangeMax, "RangeMax", 360.0, 0);
	mParameters.setMin(0.0);
	mParameters.setInstantiation();
	mParameters.set(serialisedParameters);
}

FrameLib_CircleMean::ParameterInfo FrameLib_CircleMean::sParamInfo;

FrameLib_CircleMean::ParameterInfo::ParameterInfo()
{

}

// Info

std::string FrameLib_CircleMean::objectInfo(bool verbose)
{
    return formatInfo("Output the cummulative sum of the input frame",
                   "Output the cummulative sum of the input frame", verbose);
}

std::string FrameLib_CircleMean::inputInfo(unsigned long idx, bool verbose)
{
    return formatInfo("Input Frames", "Input Frames", verbose);
}

std::string FrameLib_CircleMean::outputInfo(unsigned long idx, bool verbose)
{
    return "cumulative sum of the frame";
}

// Process
void FrameLib_CircleMean::process()
{
	unsigned long sizeIn, sizeOut2, sizeOut;
    const double *input = getInput(0, &sizeIn);
    requestOutputSize(0, sizeIn);
	requestOutputSize(1, 1);
    allocateOutputs();

    double *output = getOutput(0, &sizeOut);
	double *output2 = getOutput(1, &sizeOut2);
	
	auto in_angle = alloc<double>(sizeIn);
	auto cos_angle = alloc<double>(sizeIn);
	auto sin_angle = alloc<double>(sizeIn);

	auto diffs = alloc<double>(2);
	auto absDiffs = alloc<double>(2);
	double avg;
	int argMinInd;

	auto range = mParameters.getValue(kRangeMax);

    if (output)
    {
		for (unsigned int i = 0; i < sizeIn; i++) {
			in_angle[i] = input[i] * (360.0 / range);
			in_angle[i] = in_angle[i] * (M_PI / 180);
			cos_angle[i] = cos(in_angle[i]);
			sin_angle[i] = sin(in_angle[i]);
		}
		avg = atan2(statMean(sin_angle, sizeIn), statMean(cos_angle, sizeIn));
		avg *= (180 / M_PI);
		avg = avg / (360 / range);

		for (unsigned int i = 0; i < sizeIn; i++) {
			diffs[0] = avg - input[i];
			diffs[1] = avg - (input[i] - range);
			for (int j = 0; j < 2; j++) {
				absDiffs[j] = abs(diffs[j]);
			}

			argMinInd = (int) statMinPosition(absDiffs, 2);
			*output2 = avg;
			output[i] = -diffs[argMinInd];
		}
    }
	dealloc(diffs);
	dealloc(absDiffs);

}
