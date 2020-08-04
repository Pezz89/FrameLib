#include "FrameLib_CircleMean.h"

/*
Implementation of the YIN fundamental frequency (f0) estimation algorithm, as described in:
[1] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.                                                                                                                                                           3 of America, 111(4), 1917-1930.⏎
*/

// Constructor

FrameLib_CircleMean::FrameLib_CircleMean(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, &sParamInfo, 1, 2)
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
	auto in_180 = alloc<double>(sizeIn);
	auto in_90 = alloc<double>(sizeIn);
	auto in_270 = alloc<double>(sizeIn);
	auto diffs = alloc<double>(4);
	auto absDiffs = alloc<double>(4);
	double avg, avg90, avg180, avg270, diff0, diff90, diff180, diff270;
	int argMinInd;

	auto range = mParameters.getValue(kRangeMax);
	for (int i = 0; i < sizeIn; i++) {
		in_180[i] = fmod(fmod((input[i] - (range*0.5)), range) + range, range);
		in_90[i] = fmod(fmod((input[i] - (range*0.33)), range) + range, range);
		in_270[i] = fmod(fmod((input[i] - (range*0.66)), range) + range, range);
	}
    if (output)
    {
		avg180 = fmod(fmod(statMean(in_180, sizeIn) + (range*0.5), range) + range, range);
		avg90 = fmod(fmod(statMean(in_90, sizeIn) + (range*0.33), range) + range, range);
		avg270 = fmod(fmod(statMean(in_270, sizeIn) + (range*0.66), range) + range, range);
		avg = statMean(input, sizeIn);
		for (int i = 0; i < sizeIn; i++) {
			diffs[0] = avg - input[i];
			diffs[1] = avg90 - in_90[i];
			diffs[2] = avg180 - in_180[i];
			diffs[3] = avg270 - in_270[i];
			for (int j = 0; j < 4; j++) {
				absDiffs[j] = abs(diffs[j]);
			}

			argMinInd = (int) statArgMin(absDiffs, 4);
			*output2 = (double)argMinInd;
			output[i] = diffs[argMinInd];
		}
    }
	dealloc(diffs);
	dealloc(absDiffs);
	dealloc(in_180);
	dealloc(in_90);
	dealloc(in_270);
}
