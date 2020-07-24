#include "FrameLib_Yin.h"

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

	const unsigned int tau_min = (unsigned int)floor(mSamplingRate / mParameters.getValue(kF0Max));
	const unsigned int tau_max = (unsigned int)std::min(floor(mSamplingRate / mParameters.getValue(kF0Min)), (double)sizeIn);

    requestOutputSize(0, 1);
	requestOutputSize(1, 1);
    allocateOutputs();
    double *output = getOutput(0, &sizeOut);
	double *harmonicity = getOutput(1, &sizeOut);
	

    if (output)
    {
		auto df = alloc<double>(sizeIn);
		auto cmndf = alloc<double>(sizeIn);
		this->differenceFunction_slow(input, df, sizeIn, tau_max);
		this->cmndf(df, cmndf, sizeIn);
		this->getPitch(cmndf, df, output, harmonicity, tau_min, tau_max, mParameters.getValue(kHarmoThresh));
		dealloc(df);
		dealloc(cmndf);
    }

}

void FrameLib_Yin::differenceFunction(const double * x, double * df, unsigned int N, unsigned int tauMax)
{
	/*
	Implement the "Difference function" as expressed in equation 6 of [1]
	*/
	tauMax = std::min(tauMax, N);
	auto x_cumsum = alloc<double>(N+1);
	x_cumsum[0] = x[0];
	std::transform(x, x + N, x_cumsum + 1, [](const double val) {return std::pow(val, 2.0);});
	std::partial_sum(x_cumsum+1, x_cumsum + 1 + N, x_cumsum + 1);

	unsigned long convSize = mProcessor.convolved_size(N, N, EdgeMode::kEdgeLinear);
	auto conv = alloc<double>(convSize);

	// TODO: This is wasting memory, would be better to iterate over x in reverse or in covolution using iterators
	auto x_rev = alloc<double>(N);
	std::reverse_copy(x, x + N, x_rev);
	mProcessor.convolve(conv, { x, N }, { x_rev, N }, EdgeMode::kEdgeLinear);

	for (unsigned int i = 0; i < tauMax; i++) {
		df[i] = x_cumsum[N - i] + x_cumsum[N] - x_cumsum[i] - 2.0 * conv[(N - 1) + i];
	}
	dealloc(x_cumsum);
	dealloc(conv);
	dealloc(x_rev);
}

void FrameLib_Yin::differenceFunction_slow(const double * x, double * output, unsigned int N, unsigned int tauMax)
{
	tauMax = std::min(tauMax, N);
	unsigned int j, tau;
	double tmp;
	for (tau = 0; tau < tauMax; tau++) {
		output[tau] = 0.;
	}
	for (tau = 1; tau < tauMax; tau++) {
		for (j = 0; j < N-tauMax; j++) {
			tmp = x[j] - x[j + tau];
			output[tau] += tmp * tmp;
		}
	}
}

void FrameLib_Yin::cmndf(double * df, double * cmndf, unsigned int tau_max)
{
	/*
	Converts the difference function to a "Cumulative Normalised Difference Function" as expressed in equation 8 of [1]
	*/
	double rolling_sum = 0.0;
	for (unsigned int i = 1; i < tau_max; i++) {
		rolling_sum += df[i];
		cmndf[i] = df[i] * i / rolling_sum;
	}
	cmndf[0] = 1.0;
}

void FrameLib_Yin::getPitch(double * cmndf, double * df, double * f, double * harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th)
{
	unsigned int tau = tau_min;
	f[0] = 0.0;
	// Find the first peak that is above the harmonicity threshold and is beyond the minimum frequency
	while (tau < tau_max) {
		if (cmndf[tau] < (1.0 - harmo_th)) {
			while ((tau + 1 < tau_max) && (cmndf[tau + 1] < cmndf[tau])) {
				tau++;
			}
			break;
		}
		tau++;
	}
	
	bool is_voiced = true;
	if (tau == tau_max) {
		// If no f0 was found, calculate the harmonicity only
		tau = *std::min_element(cmndf + tau_min, cmndf + tau_max);
		is_voiced = false;
	}
	// Parabolic interpolation requires a sample on either side of the current tau value 
	if (tau > tau_min && tau < (tau_max - 1)) {
		is_voiced = true;
		// Interpolate 2 samples around estimate to increase accuracy of f0 and harmonicity values
		// Use the difference function as opposed to cmndf for frequency to ensure unbiased interpolation, as per [1]

		//a = (x1 + x3 - 2 * x2) / 2;
		//b = (x3 - x1) / 2;

		*f = 1 / 2. * (cmndf[tau - 1] - cmndf[tau + 1]) / (cmndf[tau - 1] - 2 * cmndf[tau] + cmndf[tau + 1]) + tau;
		// cnmdf used for harmonicity
		*harm = cmndf[tau];
		*f = mSamplingRate / *f;
	}
	else {
		f[0] = mSamplingRate / tau;
		harm[0] = cmndf[tau];
	}
	if (is_voiced == false) {
		// Check for nans as a result of DC signals
		isnan(cmndf[tau]) ? harm[0] = 1.0 : harm[0] = cmndf[tau];
		f[0] = 0.0;
	}
}
