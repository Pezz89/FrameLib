#ifndef FRAMELIB_CIRCLEDISTANCE_H
#define FRAMELIB_CIRCLEDISTANCE_H

#include "FrameLib_DSP.h"
#include "../../FrameLib_Dependencies/SpectralProcessor.hpp"
#include "../../FrameLib_Dependencies/Statistics.hpp"
#include <numeric>
#include <iterator>
#include <algorithm>
#include <math.h>

// Binary Operator


class FrameLib_CircleDistance final : public FrameLib_Processor
{
	// Parameter Enums and Info

	struct ParameterInfo : public FrameLib_Parameters::Info
	{
		ParameterInfo()
		{
			add("Sets the mode used when dealing with mismatched input lengths: "
				"wrap - the smaller input is read modulo against the larger input. "
				"shrink - the output length is set to that of the smaller input. "
				"pad_in - the smaller input is padded prior to calculation to match the larger input. "
				"pad_out - the output is padded to match the length of the larger input.");
			add("Sets which inputs trigger output.");
			add("Sets the value used for padding (for either pad_in or pad_out modes).");
		}
	};

	enum ParameterList { kMismatchMode, kTriggers, kPadding, kLength, kDCFilterLength};
	enum MismatchModes { kWrap, kShrink, kPadIn, kPadOut };
	enum TriggerModes { kBoth, kLeft, kRight };

public:

	// Constructor

	FrameLib_CircleDistance(FrameLib_Context context, const FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy)
		: FrameLib_Processor(context, proxy, getParameterInfo(), 3, 1)
	{
		mParameters.addEnum(kMismatchMode, "mismatch");
		mParameters.addEnumItem(kWrap, "wrap");
		mParameters.addEnumItem(kShrink, "shrink");
		mParameters.addEnumItem(kPadIn, "pad_in");
		mParameters.addEnumItem(kPadOut, "pad_out");
		mParameters.setInstantiation();

		mParameters.addEnum(kTriggers, "trigger_ins");
		mParameters.addEnumItem(kBoth, "both");
		mParameters.addEnumItem(kLeft, "left");
		mParameters.addEnumItem(kRight, "right");
		mParameters.setInstantiation();

		mParameters.addDouble(kPadding, "pad", 0.0);
		mParameters.setInstantiation();

		mParameters.addDouble(kLength, "length", 12.0);
		mParameters.setInstantiation();

		mParameters.addInt(kDCFilterLength, "DC_filter_length", 10);
		mParameters.setInstantiation();

		mParameters.set(serialisedParameters);

		mMismatchMode = static_cast<MismatchModes>(mParameters.getInt(kMismatchMode));
		mPadValue = mParameters.getValue(kPadding);

		history_length = mParameters.getInt(kDCFilterLength);
		wrap_history = alloc<double>(history_length);
		std::fill(wrap_history, wrap_history + history_length, 0.0);
		TriggerModes triggers = (TriggerModes)mParameters.getInt(kTriggers);

		if (triggers == kLeft)
			setInputMode(1, false, false, false);
		if (triggers == kRight)
			setInputMode(0, false, false, false);
	}
	~FrameLib_CircleDistance() {
		dealloc(wrap_history);
	}

	// Info

	std::string objectInfo(bool verbose) override
	{
		return formatInfo("#: Calculation is performed on pairs of values in turn. "
			"The output is a frame at least as long as the smaller of the two inputs. "
			"When inputs mismatch in length the result depends on the mismatch parameter. "
			"Either or both inputs may be set to trigger output.",
			"#.", getDescriptionString(), verbose);
	}

	std::string inputInfo(unsigned long idx, bool verbose) override { return idx ? "Right Operand" : "Left Operand"; }
	std::string outputInfo(unsigned long idx, bool verbose) override { return "Output"; }

private:
	double *wrap_history;
	size_t history_length;
	double offset = 0.0;
	double hist_mean = 0.0;
	double last_samp = 0.0;
	size_t i = 0;
	// Process
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

	double distance_unwrap(int len, double idx_1, double idx_2) {
		double wrapped_res = distance(len, idx_1, idx_2);
		double wrap_diff = wrapped_res - last_samp;

		if (wrap_diff > (len / 2.0)) {
			wrapped_res -= len;
		} 
		else if (wrap_diff < -(len / 2.0)) {
			wrapped_res += len;
		}
		last_samp = wrapped_res;
		wrap_history[i] = wrapped_res;
		hist_mean = statMean(wrap_history, history_length);
		wrapped_res -= hist_mean;
		//wrap_history[i] = wrapped_res;
		i++;
		i %= mParameters.getInt(kDCFilterLength);
		return wrapped_res;
	}

	void process() override
	{
		MismatchModes mode = mMismatchMode;

		unsigned long sizeIn1, sizeIn2, sizeIn3, sizeCommon, sizeOut;

		const double *input1 = getInput(0, &sizeIn1);
		const double *input2 = getInput(1, &sizeIn2);
		double defaultValue = mPadValue;

		// Get common size

		sizeCommon = std::min(sizeIn1, sizeIn2);

		// Calculate output size by mode

		switch (mode)
		{
		case kShrink:
			sizeOut = sizeCommon;
			break;
		default:
			sizeOut = std::max(sizeIn1, sizeIn2);
			if (mode == kWrap)
				sizeOut = sizeIn1 && sizeIn2 ? sizeOut : 0;
			break;
		}

		// Allocate output

		requestOutputSize(0, sizeOut);
		allocateOutputs();
		double *output = getOutput(0, &sizeOut);
		sizeCommon = std::min(sizeCommon, sizeOut);

		if (!sizeOut)
			return;

		// Do first part

		for (unsigned long i = 0; i < sizeCommon; i++)
			output[i] = distance_unwrap(mParameters.getValue(kLength) , input1[i], input2[i]);

		// Clean up if sizes don't match

		if (sizeIn1 != sizeIn2)
		{
			switch (mode)
			{
			case kShrink:
				break;

			case kWrap:

				if (sizeIn1 > sizeIn2)
				{
					if (sizeIn2 == 1)
					{
						double value = input2[0];
						for (unsigned long i = 1; i < sizeOut; i++)
							output[i] = distance_unwrap(mParameters.getValue(kLength), input1[i], value);
					}
					else
					{
						for (unsigned long i = sizeCommon; i < sizeOut;)
							for (unsigned long j = 0; j < sizeIn2 && i < sizeOut; i++, j++)
								output[i] = distance_unwrap(mParameters.getValue(kLength), input1[i], input2[j]);
					}
				}
				else
				{
					if (sizeIn1 == 1)
					{
						double value = input1[0];
						for (unsigned long i = 1; i < sizeOut; i++)
							output[i] = distance_unwrap(mParameters.getValue(kLength), value, input2[i]);
					}
					else
					{
						for (unsigned long i = sizeCommon; i < sizeOut;)
							for (unsigned long j = 0; j < sizeIn1 && i < sizeOut; i++, j++)
								output[i] = distance_unwrap(mParameters.getValue(kLength), input1[j], input2[i]);
					}
				}
				break;

			case kPadIn:

				if (sizeIn1 > sizeIn2)
				{
					for (unsigned long i = sizeCommon; i < sizeOut; i++)
						output[i] = distance_unwrap(mParameters.getValue(kLength), input1[i], defaultValue);
				}
				else
				{
					for (unsigned long i = sizeCommon; i < sizeOut; i++)
						output[i] = distance_unwrap(mParameters.getValue(kLength), defaultValue, input2[i]);
				}
				break;

			case kPadOut:

				for (unsigned long i = sizeCommon; i < sizeOut; i++)
					output[i] = defaultValue;
				break;
			}
		}
	}

	// Description (specialise to change description)

	const char *getDescriptionString() { return "Binary Operator - No operator info available"; }

	ParameterInfo *getParameterInfo()
	{
		static ParameterInfo info;
		return &info;
	}

	// Data

	double mPadValue;
	MismatchModes mMismatchMode;
};

#endif

