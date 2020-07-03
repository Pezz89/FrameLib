#include "FrameLib_Yin.h"


// Constructor

FrameLib_Yin::FrameLib_Yin(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, nullptr, 1, 1)
{
	mParameters.addDouble(kF0Min, "f0Min", 0);
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addEnum(kF0Max, "f0Max", 1);
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addEnum(kHarmoThresh, "HarmoThresh", 2);
	mParameters.setClip(0.0, 1.0);
	mParameters.setInstantiation();

	mParameters.set(serialisedParameters);
}

// Info

std::string FrameLib_Yin::objectInfo(bool verbose)
{
    return formatInfo("Output the cummulative sum of the input frame",
                   "Output the cummulative sum of the input frame", verbose);
}

std::string FrameLib_Yin::inputInfo(unsigned long idx, bool verbose)
{
    return formatInfo("Input Frames", "Input Frames", verbose);
}

std::string FrameLib_Yin::outputInfo(unsigned long idx, bool verbose)
{
    return "cumulative sum of the frame";
}

// Process
void FrameLib_Yin::process()
{

	const double tauMin = floor(mParameters.getValue(kF0Min) / mSamplingRate);
	const double tauMax = floor(mParameters.getValue(kF0Max) / mSamplingRate);
	
	unsigned long sizeIn;
    const double *input = getInput(0, &sizeIn);

    unsigned long sizeOut = sizeIn;

    requestOutputSize(0, sizeOut);
    allocateOutputs();

    double *output = getOutput(0, &sizeOut);

    if (output)
    {
		output[0] = input[0];
        for (unsigned long i=1; i < sizeOut; i++)
        {
			output[i] = input[i] + input[i - 1];
        }
    }
}

const double * FrameLib_Yin::differenceFunction(double * x, unsigned int N, double tauMax)
{
	tauMax = std::min(tauMax, static_cast<double>(N));

	return nullptr;
}

