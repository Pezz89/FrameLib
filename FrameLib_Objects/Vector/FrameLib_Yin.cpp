#include "FrameLib_Yin.h"


// Constructor

FrameLib_Yin::FrameLib_Yin(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, nullptr, 1, 1), mProcessor(*this)
{
	mParameters.addDouble(kF0Min, "f0Min");
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addDouble(kF0Max, "f0Max");
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addDouble(kHarmoThresh, "HarmoThresh");
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

    unsigned long sizeOut = mProcessor.convolved_size(sizeIn, sizeIn, EdgeMode::kEdgeLinear);;

    requestOutputSize(0, sizeOut);
    allocateOutputs();

    double *output = getOutput(0, &sizeOut);



    if (output)
    {
		this->differenceFunction(input, output, sizeIn, tauMax);
    }
}

void FrameLib_Yin::differenceFunction(const double * x, double * output, unsigned int N, double tauMax)
{
	tauMax = std::min(tauMax, static_cast<double>(N));
	auto x_cumsum = std::make_unique<double[]>(N+1);
	x_cumsum[0] = 0;
	std::transform(x, x + N, x_cumsum.get() + 1, [](const double val) {return std::pow(val, 2.0);});
	std::partial_sum(x_cumsum.get()+1, x_cumsum.get() + 1 + N, x_cumsum.get() + 1);

	
	unsigned long convSize = mProcessor.convolved_size(N, N, EdgeMode::kEdgeLinear);
	auto conv = std::make_unique<double[]>(convSize);
	double x2[10] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
	mProcessor.convolve(output, { x, N }, { x2, N }, EdgeMode::kEdgeLinear);
}

