
#include "FrameLib_MedianFilter.h"
#include "FrameLib_Sort_Functions.h"
#include <algorithm>

// Constructor

FrameLib_MedianFilter::FrameLib_MedianFilter(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy) : FrameLib_Processor(context, proxy, &sParamInfo, 2, 1)
{
    mParameters.addInt(kWidth, "width", 1, 0);
    mParameters.setMin(1);
    
    mParameters.addDouble(kPadding, "pad", 0.0, 1);
    
    mParameters.addEnum(kMode, "mode", 3);
    mParameters.addEnumItem(kPad, "pad");
    mParameters.addEnumItem(kWrap, "wrap");
    mParameters.addEnumItem(kFold, "fold");
    
    mParameters.set(serialisedParameters);
    
    setParameterInput(1);
}

// Info

std::string FrameLib_MedianFilter::objectInfo(bool verbose)
{
    return formatInfo("Median filter an input frame: The output is the same size as the input. "
                   "Each output value is the median of the area surrounding the input value. "
                   "The width of the area, and the edge behaviours are controllable.",
                   "Median filter an input frame.", verbose);
}

std::string FrameLib_MedianFilter::inputInfo(unsigned long idx, bool verbose)
{
    if (idx)
        return parameterInputInfo(verbose);
    else
        return "Input Frames";
}

std::string FrameLib_MedianFilter::outputInfo(unsigned long idx, bool verbose)
{
    return "Median Filtered Frames";
}

// Parameter Info

FrameLib_MedianFilter::ParameterInfo FrameLib_MedianFilter::sParamInfo;

FrameLib_MedianFilter::ParameterInfo::ParameterInfo()
{
    add("Sets the width of the median filtering in samples.");
    add("Sets the padding value.");
    add("Sets the mode that controls the edge behaviour: "
        "pad - the edges are treated as though infinitely padded with the padding value. "
        "wrap - the edges are treated as though the frame is wrapped cyclically. "
        "fold - the edges are treated as though they fold over (suitable for spectral purposes).");
}

// Helpers

double FrameLib_MedianFilter::insertMedian(double *temp, unsigned long *indices, double value, long index, long width)
{
    long current, insert, gap;
    
    // Insert actual value
    
    temp[index] = value;
    
    // Find sort position
    
    current = -10000000;
    
    for (long i = 0; i < width; i++)
    {
        if (indices[i] == index)
        {
            current = i;
            break;
        }
    }
    
    // Search right
    
    for (insert = current, gap = std::max(1L, (width - current) >> 1); gap; gap >>= 1)
    {
        for (long i = insert + gap; i < width; i += gap)
        {
            if (temp[indices[i]] >= value)
                break;
            else
                insert = i;
        }
    }
    
    // Swaps
    
    for (long i = current; i < insert; i++)
        indices[i] = indices[i + 1];
    indices[insert] = index;
    
    // Search left
    
    for (current = insert, gap = std::max(1L, current >> 1); gap; gap >>= 1)
    {
        for (long i = insert - gap; i >= 0; i -= gap)
        {
            if (temp[indices[i]] <= value)
                break;
            else
                insert = i;
        }
    }
    
    // Swaps
    
    for (long i = current; i > insert; i--)
        indices[i] = indices[i - 1];
    indices[insert] = index;
    
    return temp[indices[width >> 1]];
}

double FrameLib_MedianFilter::median(double *temp, unsigned long *indices, long width)
{
    sortIndicesAscending(indices, temp, width);
    return temp[indices[width >> 1]];
}

double FrameLib_MedianFilter::getPad(const double *input, long index, long sizeIn, double padValue)
{
    return (index >= 0 && index < sizeIn) ? input[index] : padValue;
}

double FrameLib_MedianFilter::getWrap(const  double *input, long index, long sizeIn)
{
    index %= sizeIn;
    index = index < 0 ? index + sizeIn : index;
    
    return input[index];
}

double FrameLib_MedianFilter::getFold(const double *input, long index, long sizeIn)
{
    index = std::abs(index) % ((sizeIn - 1) * 2);
    index = index > (sizeIn - 1) ? ((sizeIn - 1) * 2) - index : index;
    
    return input[index];
}

// Process

void FrameLib_MedianFilter::process()
{
    // Get Input
    
    unsigned long sizeIn, sizeOut;
    long width = mParameters.getInt(kWidth);
    const double *input = getInput(0, &sizeIn);
    double padValue = mParameters.getValue(kPadding);
    Modes mode = static_cast<Modes>(mParameters.getInt(kMode));
    
    requestOutputSize(0, sizeIn);
    allocateOutputs();
    
    double *output = getOutput(0, &sizeOut);    
    double *temp = alloc<double>(width);
    unsigned long *indices = alloc<unsigned long>(width);
    
    // Do filtering
    
    if (sizeOut)
    {
        switch (mode)
        {
            case kWrap:
                for (long i = 0; i < width; i++)
                    temp[i] = getWrap(input, i - (width >> 1), sizeIn);
                output[0] = median(temp, indices, width);
                for (long i = 1; i < static_cast<long>(sizeIn); i++)
                {
                    double newValue = getWrap(input, i + (width >> 1), sizeIn);
                    output[i] = insertMedian(temp, indices, newValue, (i - 1) % width, width);
                }
                break;
                
            case kPad:
                for (long i = 0; i < width; i++)
                    temp[i] = getPad(input, i - (width >> 1), sizeIn, padValue);
                output[0] = median(temp, indices, width);
                for (long i = 1; i < static_cast<long>(sizeIn); i++)
                {
                    double newValue = getPad(input, i + (width >> 1), sizeIn, padValue);
                    output[i] = insertMedian(temp, indices, newValue, (i - 1) % width, width);
                }
                break;
                
            case kFold:
                for (long i = 0; i < width; i++)
                    temp[i] = getFold(input, i - (width >> 1), sizeIn);
                output[0] = median(temp, indices, width);
                for (long i = 1; i < static_cast<long>(sizeIn); i++)
                {
                    double newValue = getFold(input, i + (width >> 1), sizeIn);
                    output[i] = insertMedian(temp, indices, newValue, (i - 1) % width, width);
                }
                break;
        }
    }
    
    dealloc(temp);
    dealloc(indices);
}
