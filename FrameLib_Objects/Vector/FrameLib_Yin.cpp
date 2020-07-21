#include "FrameLib_Yin.h"


// Constructor

FrameLib_Yin::FrameLib_Yin(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, &sParamInfo, 1, 2), mProcessor(*this)
{
	mParameters.addDouble(kF0Min, "f0Min", 0.0, 0);
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addDouble(kF0Max, "f0Max", floor(mSamplingRate*0.5), 1);
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addDouble(kHarmoThresh, "HarmoThresh", 0.0, 2);
	mParameters.setClip(0.0, 1.0);
	mParameters.setInstantiation();

	mParameters.set(serialisedParameters);

}

FrameLib_Yin::ParameterInfo FrameLib_Yin::sParamInfo;

FrameLib_Yin::ParameterInfo::ParameterInfo()
{
	add("Sets the maximum output length. The output length will be M + N - 1 where M and N are the sizes of the two inputs respectively");
	add("Sets the type of input expected / output produced.");
	add("Sets the type of input expected / output produced.");
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
	unsigned long sizeIn, sizeOut;
    const double *input = getInput(0, &sizeIn);
	auto buffer = std::make_unique<double[]>(sizeIn);

	const unsigned int tauMin = (unsigned int) floor(mSamplingRate / mParameters.getValue(kF0Max));
	const unsigned int tauMax = (unsigned int) std::min(floor(mSamplingRate / mParameters.getValue(kF0Min)), (double) sizeIn);

    requestOutputSize(0, 1);
	requestOutputSize(1, 1);
    allocateOutputs();
    double *output = getOutput(0, &sizeOut);
	double *harmonicity = getOutput(1, &sizeOut);
	

    if (output)
    {
		this->differenceFunction(input, buffer.get(), sizeIn, tauMin, tauMax);
		this->cmndf(buffer.get(), buffer.get(), sizeIn);
		this->getPitch(buffer.get(), output, harmonicity, tauMin, tauMax, mParameters.getValue(kHarmoThresh));

		if (*output > 0.0) {
			output[0] = mSamplingRate / output[0];
		}
		else {
			double* min_ind = std::min_element(buffer.get(), buffer.get() + sizeIn);
		}
    }
}

void FrameLib_Yin::differenceFunction(const double * x, double * output, unsigned int N, unsigned int tauMin, unsigned int tauMax)
{
	tauMax = std::min(tauMax, N);
	auto x_cumsum = std::make_unique<double[]>(N+1);
	x_cumsum[0] = x[0];
	std::transform(x, x + N, x_cumsum.get() + 1, [](const double val) {return std::pow(val, 2.0);});
	std::partial_sum(x_cumsum.get()+1, x_cumsum.get() + 1 + N, x_cumsum.get() + 1);

	unsigned long convSize = mProcessor.convolved_size(N, N, EdgeMode::kEdgeLinear);
	auto conv = std::make_unique<double[]>(convSize);
	// TODO: This is wasting memory, would be better to iterate over x in reverse or in covolution using iterators
	auto x_rev = std::make_unique<double[]>(N);
	std::reverse_copy(x, x + N, x_rev.get());
	mProcessor.convolve(conv.get(), { x, N }, { x_rev.get(), N }, EdgeMode::kEdgeLinear);

	for (unsigned int i = tauMin; i < tauMax; i++) {
		output[i] = x_cumsum[N - i] + x_cumsum[N] - x_cumsum[i] - 2.0 * conv[(N - 1) + i];
	}
}

void FrameLib_Yin::cmndf(double * df, double * output, unsigned int N)
{
	double rolling_sum = 0.0;
	for (int i = 1; i < N; i++) {
		rolling_sum += df[i];
		output[i] = df[i] * i / rolling_sum;
	}
	output[0] = 1.0;
}

void FrameLib_Yin::getPitch(double * cmndf, double * f, double * harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th)
{
	unsigned int tau = tau_min;
	f[0] = 0.0;
	// Find the first peak that is above the harmonicity threshold and is beyond the minimum frequency
	while (tau < tau_max) {
		if (cmndf[tau] < harmo_th) {
			while ((tau + 1 < tau_max) && (cmndf[tau + 1] < cmndf[tau])) {
				tau++;
			}
			break;
		}
		tau++;
	}

	// Parabolic interpolation requires a sample on either side of the current tau value 
	if (tau > tau_min && tau < (tau_max - 1)) {
		// Interpolate 3 samples around estimate to increase accuracy of f0 and harmonicity values
		*f = 1 / 2. * (cmndf[tau - 1] - cmndf[tau + 1]) / (cmndf[tau - 1] - 2 * cmndf[tau] + cmndf[tau + 1]) + tau;
		*harm = cmndf[tau] - 1 / 4. * (cmndf[tau - 1] - cmndf[tau + 1]) * (*f - tau);
	}
	else if (tau == tau_max) {
		// If no f0 was found, calculate the harmonicity only
		tau = *std::min_element(cmndf+tau_min, cmndf + tau_max);
		// Parabolic interpolation requires a sample on either side of the current tau value
		if (tau > tau_min && tau < (tau_max - 1)) {
			*harm = cmndf[tau] - 1 / 4. * (cmndf[tau - 1] - cmndf[tau + 1]) * (*f - tau);
		}
		else {
			harm[0] = cmndf[tau];
		}
	}
	else {
		f[0] = tau;
		harm[0] = cmndf[tau];
	}

}
