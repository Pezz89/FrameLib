#include "FrameLib_CumSum.h"

// Constructor

FrameLib_CumSum::FrameLib_CumSum(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, nullptr, 1, 1)
{}

// Info

std::string FrameLib_CumSum::objectInfo(bool verbose)
{
    return formatInfo("Output the cummulative sum of the input frame",
                   "Output the cummulative sum of the input frame", verbose);
}

std::string FrameLib_CumSum::inputInfo(unsigned long idx, bool verbose)
{
    return formatInfo("Input Frames", "Input Frames", verbose);
}

std::string FrameLib_CumSum::outputInfo(unsigned long idx, bool verbose)
{
    return "cumulative sum of the frame";
}

// Process
void FrameLib_CumSum::process()
{
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

