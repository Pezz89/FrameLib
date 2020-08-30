#include "FrameLib_Yin.h"

/*
Implementation of the YIN fundamental frequency (f0) estimation algorithm, as described in:
[1] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.                                                                                                                                                           3 of America, 111(4), 1917-1930.⏎
*/

// Constructor

FrameLib_Yin::FrameLib_Yin(FrameLib_Context context, const FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, &sParamInfo, 1, 4), mProcessor(*this)
{
	mParameters.addDouble(kF0Min, "f0Min", 0.0, 0);
	// Maximum frequency must be one sample less than the nyquist rate, so that all frequencies can be properly interpolated
	mParameters.setClip(0.0, mSamplingRate*((mSamplingRate/floor(mSamplingRate*0.5))-1.0));
	mParameters.setInstantiation();
	mParameters.addDouble(kF0Max, "f0Max", mSamplingRate*((mSamplingRate / floor(mSamplingRate*0.5)) - 1.0), 1);
	mParameters.setClip(0.0, floor(mSamplingRate*0.5));
	mParameters.setInstantiation();
	mParameters.addDouble(kHarmoThresh, "HarmoThresh", 0.0, 2);
	mParameters.setMin(0.0);
	mParameters.setInstantiation();

	mParameters.set(serialisedParameters);

}

FrameLib_Yin::ParameterInfo FrameLib_Yin::sParamInfo;

FrameLib_Yin::ParameterInfo::ParameterInfo()
{
	add("Sets the minimum fundamental frequency that can be detected by YIN.");
	add("Sets the maximum fundamental frequency that can be detected by YIN.");
	add("Sets the harmonicity threshold. Harmonicities above this threshold are classified as inharmonic, and the minimum of the CMNDF is taken.");
}

// Info

std::string FrameLib_Yin::objectInfo(bool verbose)
{
    return formatInfo("Implementation of the YIN fundamental frequency estimation algorithm", "Implementation of the YIN fundamental frequency estimation algorithm", verbose);
}

std::string FrameLib_Yin::inputInfo(unsigned long idx, bool verbose)
{
    return formatInfo("Input Frames", "Input Frames", verbose);
}

std::string FrameLib_Yin::outputInfo(unsigned long idx, bool verbose)
{
	switch (idx) {
		case 0:
			return formatInfo("Fundamental frequency (hz) frame", "Fundamental frequency (hz) frame", idx, verbose);
			break;
		case 1:
			return formatInfo("Harmonicity", "Harmonicity", idx, verbose);
			break;
		case 2:
			return formatInfo("Cumulative normalised distribution function", "Cumulative normalised distribution function", idx, verbose);
			break;
		case 3:
			return formatInfo("Tau", "Tau", idx, verbose);
			break;
		default:
			return "";
			break;
	}
		
}

// Process
void FrameLib_Yin::process()
{
	unsigned long sizeIn, sizeOut, sizeCMNDF;
    const double *input = getInput(0, &sizeIn);
	// double * tau_out;

	// Add the extra sample removed during parameter clipping to allow all frequencies to be interpolated
	const unsigned int tau_min = (unsigned int)floor(mSamplingRate / mParameters.getValue(kF0Max));
	const unsigned int tau_max = (unsigned int)floor(mSamplingRate / mParameters.getValue(kF0Min));

    requestOutputSize(0, 1);
	requestOutputSize(1, 1);
	//requestOutputSize(2, tau_max + 1);
	//requestOutputSize(3, 1);
    allocateOutputs();
    double *output = getOutput(0, &sizeOut);
	double *harmonicity = getOutput(1, &sizeOut);

	// Output cmndf function and tau to plot and debug in real-time
	//double *cmndf = getOutput(2, &sizeCMNDF);
	//double *tau_out = getOutput(3, &sizeOut);

    if (output)
    {
		auto w_len = sizeIn - tau_max;
		auto df = alloc<double>(tau_max+1);
		auto cmndf = alloc<double>(tau_max+1);
		this->differenceFunction_slow(input, df, w_len, tau_max);
		this->cmndf(df, cmndf, tau_max);
		// this->getPitch(cmndf, df, output, harmonicity, tau_min, tau_max, mParameters.getValue(kHarmoThresh), tau_out);
		this->getPitch(cmndf, df, output, harmonicity, tau_min, tau_max, mParameters.getValue(kHarmoThresh));
		dealloc(df);
		dealloc(cmndf);
    }

}

void FrameLib_Yin::differenceFunction(const double * x, double * df, unsigned int N, unsigned int tauMax)
{
	/*
	Implement the "Difference function" as expressed in equation 6 of [1]
	Fast FFT based implementation based off of python implementation. Resulting function is similar,but not identical to that implemented by DeChevigne
	*/
	auto x_cumsum = alloc<double>(N+1);
	x_cumsum[0] = x[0];
	std::transform(x, x + N, x_cumsum + 1, [](const double val) {return std::pow(val, 2.0);});
	std::partial_sum(x_cumsum+1, x_cumsum + 1 + N, x_cumsum + 1);

	uintptr_t convSize = mProcessor.convolved_size(N, N, EdgeMode::kEdgeLinear);
	auto conv = alloc<double>(convSize);

	// TODO: Replace with correlation function
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

void FrameLib_Yin::differenceFunction_slow(const double * x, double * output, unsigned int N, unsigned int tau_max)
{
	/*
	Implement the "Difference function" as expressed in equation 6 of [1], copying the method implemented by Alaine Dechevigne in his Matlab implementation.
	*/
	unsigned int m, t;
	double tmp;
	for (t = 0; t < tau_max+1; t++) {
		output[t] = 0.;
	}
	for (m = tau_max; m < tau_max+N; m++) {
		for (t = 0; t < tau_max+1; t++) {
			tmp = x[m-t] - x[m];
			// TODO: What's the best way to square a double in C++/FrameLib?
			tmp = tmp * tmp;
			output[t] += tmp;
		}
	}
	for (t = 0; t < tau_max + 1; t++) {
		output[t] /= N;
		output[t] /= 2.0;
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
		cmndf[i] = df[i] * (double) i / rolling_sum;
	}
	cmndf[0] = 1.0;
}

void FrameLib_Yin::getPitch(double * cmndf, double * df, double * f, double * harm, const unsigned int tau_min, const unsigned int tau_max, double harmo_th)
{
	unsigned int tau = tau_min;
	unsigned int best_tau = tau_min;
	// Variables to store interpolation intermediate values
	double a, b, shift, best_f;
	f[0] = 0.0;
	best_f = 0.0;
	// TODO: Replace with double maximum (how's that number found in C++?)
	double best_harm = 10000000.0;
	// Find the first peak that is above the harmonicity threshold and is at or beyond the minimum frequency
	while (tau < tau_max) {
		// If the index is a local minima...
		if ((cmndf[tau - 1] > cmndf[tau]) && (cmndf[tau + 1] > cmndf[tau])) {
			// Interpolate 2 samples around estimate to increase accuracy of f0 and harmonicity values
			// Use the difference function as opposed to cmndf for frequency to ensure unbiased interpolation, as per [1]
			a = 0.5 * (cmndf[tau - 1] + cmndf[tau + 1] - 2.0 * cmndf[tau]);
			b = 0.5 * (cmndf[tau + 1] - cmndf[tau - 1]);
			// TODO: Is this clipping correct? Or should the CMNDF function never go above 1.0 or below 0.0 even with interpolation?
			*harm = std::max((cmndf[tau] - b * b / (4.0 * a)), 0.0); // value of interpolated minimum, or 0.0 if interpolation undershoots
			if (*harm < best_harm) {
				best_harm = *harm;
				shift = -b / (2.0 * a);											// offset of interpolated minimum re current sample
				best_f = mSamplingRate / ((double)tau + shift);
				best_tau = tau;
			}
			if (*harm < harmo_th) {
				unsigned int c = std::min(tau * 2, tau_max);
				tau++;
				// If a harmonic is found, check that there isn't a lower one in the next octave, which is more likely to be the actual fundamental
				while (tau < c-1) {
					if ((cmndf[tau - 1] > cmndf[tau]) && (cmndf[tau + 1] > cmndf[tau])) {
						// Interpolate 2 samples around estimate to increase accuracy of f0 and harmonicity values
						// Use the difference function as opposed to cmndf for frequency to ensure unbiased interpolation, as per [1]
						a = 0.5 * (cmndf[tau - 1] + cmndf[tau + 1] - 2.0 * cmndf[tau]);
						b = 0.5 * (cmndf[tau + 1] - cmndf[tau - 1]);
						// TODO: Is this clipping correct? Or should the CMNDF function never go below 0.0 even with interpolation?
						*harm = std::max((cmndf[tau] - b * b / (4.0 * a)), 0.0); // value of interpolated minimum, or 0.0 if interpolation undershoots
						if (*harm < best_harm) {
							best_harm = *harm;
							shift = -b / (2.0 * a);											// offset of interpolated minimum re current sample
							best_f = mSamplingRate / ((double)tau + shift);
							best_tau = tau;
						}
					}
					tau++;
				}
				*f = best_f;
				break;
			}
		}
		tau++;
	}
	//*tau_out = (double)best_tau;
	if (tau == tau_max) {
		// Check for nans as a result of DC signals
		if (isnan(cmndf[best_tau])) {
			// TODO: calculate maximum aperiodicity value possible?
			*harm = 3.0;
			*f = 0.0;
			//*tau_out = 0.0;
			return;
		}
		*f = best_f;
		*harm = best_harm;
	}
}