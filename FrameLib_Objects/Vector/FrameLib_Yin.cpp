﻿#include "FrameLib_Yin.h"

/*
Implementation of the YIN fundamental frequency (f0) estimation algorithm, as described in:
[1] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.                                                                                                                                                           3 of America, 111(4), 1917-1930.⏎
*/

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
	auto df = std::make_unique<double[]>(sizeIn);
	auto cmndf = std::make_unique<double[]>(sizeIn);

	const unsigned int tauMin = (unsigned int) floor(mSamplingRate / mParameters.getValue(kF0Max));
	const unsigned int tauMax = (unsigned int) std::min(floor(mSamplingRate / mParameters.getValue(kF0Min)), (double) sizeIn);

    requestOutputSize(0, 1);
	requestOutputSize(1, 1);
    allocateOutputs();
    double *output = getOutput(0, &sizeOut);
	double *harmonicity = getOutput(1, &sizeOut);
	

    if (output)
    {
		this->differenceFunction(input, df.get(), sizeIn, tauMin, tauMax);
		this->cmndf(df.get(), cmndf.get(), sizeIn);
		this->getPitch(cmndf.get(), df.get(), output, harmonicity, tauMin, tauMax, mParameters.getValue(kHarmoThresh));
    }
}

void FrameLib_Yin::differenceFunction(const double * x, double * output, unsigned int N, unsigned int tauMin, unsigned int tauMax)
{
	/*
	Implement the "Difference function" as expressed in equation 6 of [1]
	*/
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
	/*
	Converts the difference function to a "Cumulative Normalised Difference Function" as expressed in equation 8 of [1]
	*/
	double rolling_sum = 0.0;
	for (int i = 1; i < N; i++) {
		rolling_sum += df[i];
		output[i] = df[i] * i / rolling_sum;
	}
	output[0] = 1.0;
}

void FrameLib_Yin::getPitch(double * cmndf, double * df, double * f, double * harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th)
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
	
	if (tau == tau_max) {
		// If no f0 was found, calculate the harmonicity only
		tau = *std::min_element(df + tau_min, df + tau_max);
	// Parabolic interpolation requires a sample on either side of the current tau value 
	if (tau > tau_min && tau < (tau_max - 1)) {
		// Interpolate 2 samples around estimate to increase accuracy of f0 and harmonicity values
		// Use the difference function as opposed to cmndf to ensure unbiased interpolation, as per [1]
		*f = mSamplingRate / (1 / 2. * (df[tau - 1] - df[tau + 1]) / (df[tau - 1] - 2 * df[tau] + df[tau + 1]) + tau);
		*harm = std::max(df[tau] - 1 / 4. * (df[tau - 1] - df[tau + 1]) * (*f - tau), 0.0);
	}
	else {
		f[0] = mSamplingRate / tau;
		harm[0] = df[tau];
	}
	if (tau == tau_max) {
		// Check for nans as a result of DC signals
		isnan(df[tau]) ? harm[0] = 1.0 : harm[0] = df[tau];
		f[0] = 0.0;
	}
}
