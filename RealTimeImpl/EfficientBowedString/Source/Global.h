/*
  ==============================================================================

    Global.h
    Created: 25 Apr 2022 12:37:28pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once
#include "../eigen/Eigen/Eigen"
#define RUN_ALL  // define this macro if you want to run all methods (reference, optimised matrix and optimised vector)
#define TIME_DOMAIN_STRING 0

namespace Global
{
    namespace Strings
    {
        static struct String
        {
            String(int aId,
                juce::String aName,
                float aRadius,
                float aDensity, 
                float aTension, 
                float aYoungMod, 
                float aLength)
            {
                mId = aId;
                mName = aName;
                mRadius = aRadius;
                mDensity = aDensity;
                mTension = aTension;
                mYoungMod = aYoungMod;
                mLength = aLength;
            }

            int mId;
            juce::String mName;
            float mRadius;
            float mDensity;
            float mTension;
            float mYoungMod;
            float mLength;
        };

        static String* kpCelloA3 = new String(1,
            "CelloA3",
            (float)3.75e-04,
            (float)3.7575e3, 
            153.f, 
            (float)25e9, 
            0.69f);
        static String* kpCelloD3 = new String(2,
            "CelloD3",
            (float)4.4e-04,
            (float)4.1104e3,
            102.6f,
            (float)25e9,
            0.69f);
        static String* kpCelloG2 = new String(3,
            "CelloG2",
            (float)6.05e-04,
            (float)5.3570e3,
            112.67f,
            (float)8.6e9,
            0.69f);
        static String* kpCelloC2 = new String(4,
            "CelloC2",
            (float)7.2e-04,
            (float)1.3017e4,
            172.74f,
            (float)22.4e9,
            0.69f);
    }
    
    static double cubicInterpolation (double* xVec, int l, double alpha)
    {
        return xVec[l - 1] * (alpha * (alpha - 1) * (alpha - 2)) / -6.0
        + xVec[l] * ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2.0
        + xVec[l + 1] * (alpha * (alpha + 1) * (alpha - 2)) / -2.0
        + xVec[l + 2] * (alpha * (alpha + 1) * (alpha - 1)) / 6.0;
    }

    static void cubicExtrapolation (double* xVec, int l, double alpha, double val)
    {
        xVec[l - 1] = xVec[l - 1] + val * (alpha * (alpha - 1) * (alpha - 2)) / -6.0;
        xVec[l] = xVec[l] + val * ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2.0;
        xVec[l + 1] = xVec[l + 1] + val * (alpha * (alpha + 1) * (alpha - 2)) / -2.0;
        xVec[l + 2] = xVec[l + 2] + val * (alpha * (alpha + 1) * (alpha - 1)) / 6.0;

    }

    static float limitOutput (float x)
    {
        if (x > 1.0)
            return 1.0;
        else if (x < -1.0)
            return -1.0;
        else
            return x;
    }
};
