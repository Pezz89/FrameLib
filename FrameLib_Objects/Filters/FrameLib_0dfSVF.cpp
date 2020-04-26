
#include "FrameLib_0dfSVF.h"

SVF::ModeType SVF::sModes
{{
    Mode("lpf", &SVF::lpf),
    Mode("bpf", &SVF::bpf),
    Mode("hpf", &SVF::hpf)
}};

SVF::ParamType SVF::sParameters
{{
    Param("freq", 500.0, Min(0.0)),
    Param("reson", 500.0, Clip(0.0, 1.0))
}};

void SVF::reset()
{
    s1 = 0.0;
    s2 = 0.0;
    lp = 0.0;
    bp = 0.0;
    hp = 0.0;
}

void SVF::updateCoefficients(double freq, double reson, double samplingRate)
{
    double srConst = 0.5 / samplingRate;

    r = std::min(std::max(1.0 - reson, 0.005), 1.0);
    g = ((2.0 * samplingRate) * tan((freq * twopi()) * srConst) * srConst);
}

double SVF::process(double x)
{
    // Compute highpass then bandpass  by applying 1st integrator to highpass output and update state
    
    hp = (x - (2.0 * r * s1) - (g * s1) - s2) / (1.0 + (2.0 * r * g) + (g * g));
    bp = g * hp + s1;
    lp = g * bp + s2;
    
    s1 = g * hp + bp;
    s2 = g * bp + lp;
    
    return 0.0;
}

double SVF::hpf(double x)
{
    return hp;
}

double SVF::bpf(double x)
{
    return bp;
}

double SVF::lpf(double x)
{
    return lp;
}
