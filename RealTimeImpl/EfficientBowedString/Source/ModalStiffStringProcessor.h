/*
  ==============================================================================

    ModalStiffStringView.h
    Created: 04/05/2022
    Author:  Riccardo Russo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Global.h"

class ModalStiffStringProcessor
{
public:
    //==========================================================================
    ModalStiffStringProcessor(double aSampleRate, Global::Strings::String* apString);
    ~ModalStiffStringProcessor();

    //==========================================================================
    //Set the time sampling step, to be called inside the PrepareToPlay
    void SetTimeStep(double aTimeStep);

    /*
    Play or pause the sound. If the sound is paused the state is not computed, 
    but the string is not reset.
    */
    void SetPlayState(bool aPlayState);

    /*
    Resets the string states, setting each oscillator to zero.
    If the PlayState is true it is set to false
    */
    void ResetStringStates();

    //Recomputes the mode for the input location at runtime
    void SetInputPos(float aNewPos);

    //Recomputes the mode for the output location at runtime
    void SetReadPos(float aNewPos);

    //Sets the gain to be multiplied to the output value
    void SetGain(float aGain);

    //Sets the bowing pressure Fb at runtime
    void SetBowPressure(float aPressure);

    //Sets the bowing speed Vb at runtime
    void SetBowSpeed(float aSpeed);

    /*
    Change the string being played.This stops the 
    playback and recomputes all the string matrices
    */
    void SetString(Global::Strings::String* apString);

    /*
    Calculates the next string state. 
    To be called for each sample inside the audio process
    */
    void ComputeState();

    //Returns the output value at the output location
    float ReadOutput();

    //Return the modes number
    int GetModesNumber();

    /*
    Take in input a vector (by reference!), empty it and fill it with the mode
    shapes at the requested location for each oscillator. Useful for visualizing
    the string state at a certain location. The location is expressed in fraction
    of string. Therefore, position = aLocationPerc * string length
    */
    void GetModesAtLocation(std::vector<float>& aModesArray, float aLocationFrac);

    /*
    Return the entire string state, i.e. the current state of each oscillator.
    Useful for visualization purposes.
    */
    std::vector<float> GetStringState();

private:
    //==========================================================================
    Global::Strings::String* mpString;

    //PlayState
    std::atomic<bool> mPlayState{ false };
    std::atomic<float> mGain{ 0.f };

    //String params
    float mRadius{ 0.f };
    float mDensity{ 0.f };
    float mTension{ 0.f };
    float mArea{ 0.f };
    float mLinDensity{ 0.f };
    float mC{ 0.f };
    float mYoungMod{ 0.f };
    float mInertia{ 0.f };
    float mK{ 0.f };
    float mLength{ 0.f };
    float mExcitPos{ 0.f };
    float mReadPos{ 0.f };
    std::vector<float> mDampCoeffs;

    //==========================================================================
    //String states
    std::vector<std::vector<float>> mStates;
    std::vector<float*> mpStatesPtrs;

    //==========================================================================
    //Bow params
    std::atomic<float> mFb;
    std::atomic<float> mVb;
    float mA{ 0.f };

    //==========================================================================
    //FDS & Modal params
    int mOversamplingFactor{ 0 };
    double mTimeStep{ 0.0 };
    int mModesNumber{ 0 };
    std::vector<float> mEigenFreqs;

    std::vector<std::vector<float>> mModesIn;
    std::atomic<float*> mpModesInCurr;
    std::atomic<float*> mpModesInNew;

    std::vector<std::vector<float>> mModesOut;
    std::atomic<float*> mpModesOutCurr;
    std::atomic<float*> mpModesOutNew;

    std::vector<int> mT11;
    std::vector<float> mT12;
    std::vector<float> mT21;
    std::vector<int> mT22;

    std::vector<float> mSchurComp;

    std::vector<int> mB11;
    std::vector<float> mB12;
    std::vector<float> mB21;
    std::vector<float> mB22;

    std::vector<float> mZeta2;
    std::vector<float> mB1;
    std::vector<float> mB2;

    std::vector<float> mZ1;
    std::vector<float> mInvAv2;
    std::vector<float> mInvAv1;

    std::vector<float> mY2;
    std::vector<float> mZ2;
    std::vector<float> mInvAb2;
    std::vector<float> mInvAb1;

    //==========================================================================
    //Utility Functions
    float ComputeEigenFreq(int aModeNumber);
    float ComputeMode(float aPos, int aModeNumber);
    float ComputeDampCoeff(float aFreq);

    void RecomputeModesNumber();
    void RecomputeEigenFreqs();
    void InitializeInModes();
    void InitializeOutModes();
    void RecomputeInModes();
    void RecomputeOutModes();
    void RecomputeDampProfile();

    void InitializeStates();
    void ResetMatrices();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalStiffStringProcessor)
};
