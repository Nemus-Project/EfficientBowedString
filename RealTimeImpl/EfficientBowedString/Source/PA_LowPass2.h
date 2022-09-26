#pragma once

#include <cmath>
#include <stdlib.h>
#include <stdio.h>

// ------------------------------------------------------------------------------------------
// Lowpass for 2X resampling. Property of Physical Audio Ltd.
// ------------------------------------------------------------------------------------------

class PA_LowPass2
{
    
public:
    
    PA_LowPass2()
    {
        double aF2coeffs[16] = {-2.67732545660392e-15, 2.16964477543704,     -4.20828990548919e-15, 1.76504194938032,
                                -1.75297444707951e-15, 0.685302479580568,     5.73987936236921e-16, 0.134084268528723,
                                 6.49025897555478e-16, 0.0128506295130287,    9.02436423742187e-17, 0.000540383728671063,
                                -1.31607067733054e-17, 7.64925129563502e-06, -9.60372445660347e-19, 1.54896336661535e-08};
        
        double bF2coeffs[16] = { 0.00140807425559308,   0.0105605569169481,   0.0492825989457579,   0.160168446573713,
                                 0.384404271776912,     0.704741164924338,    1.00677309274905,     1.13261972934269,
                                 1.00677309274905,      0.704741164924338,    0.384404271776912,    0.160168446573713,
                                 0.0492825989457579,    0.0105605569169481,   0.00140807425559308,  8.80046409745677e-05};
        
        for (int i = 0; i < 16; ++i)
        {
            m_aF[i] = aF2coeffs[i];
            m_bF[i] = bF2coeffs[i];
        }
        
        clear();
    }
    
    
    void  clear() noexcept
    {
        for (int i = 0; i < 16; ++i)
        {
            m_yold[i] = 0.0;
            m_xold[i] = 0.0;
        }
    }
    
    
    float update(float input) noexcept
    {
        double b0    = 8.80046409745677e-05;
        double outst = (double)input * b0;
        
        for (int q = 0; q < 16; ++q)
        {
            outst  += m_xold[q] * m_bF[q] - m_yold[q] * m_aF[q];
        }
        
        for (int q = 15; q > 0; --q)
        {
            m_xold[q] = m_xold[q-1];
        }
        
        for (int q = 15; q > 0; --q)
        {
            m_yold[q] = m_yold[q-1];
        }
        
        m_xold[0] = (double)input;
        m_yold[0] = outst;
        
        double output = 0.0;
        
        if (std::fabs(outst) > 1.0e-6)
        {
            output = outst;
        }
        
        float fout = (float)output;
        
        return fout;
    }
    
    
private:
    
    double m_aF[16],   m_bF[16];
    double m_xold[16], m_yold[16];
    
};


















