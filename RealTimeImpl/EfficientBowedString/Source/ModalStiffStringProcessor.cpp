#include "ModalStiffStringProcessor.h"

ModalStiffStringProcessor::ModalStiffStringProcessor (double aSampleRate, Global::Strings::String* apString)
{
    auto vPi = juce::MathConstants<float>::pi;
    mOversamplingFactor = 1;

    mTimeStep = 1.0 / (aSampleRate * mOversamplingFactor);

    mpString = apString;
    mRadius = mpString->mRadius;
    mDensity = mpString->mDensity;
    mTension = mpString->mTension;
    mLength = mpString->mLength;

    mArea = vPi * mRadius * mRadius;
    mLinDensity = mDensity * mArea;
    mInertia = (vPi * mRadius * mRadius * mRadius * mRadius) / 4;
    mK = sqrt(mYoungMod * mInertia / (mLinDensity * mLength * mLength * mLength * mLength));
    mC = sqrt(mTension / mLinDensity);

    mA = 100.f;

    RecomputeModesNumber();
    RecomputeEigenFreqs();
    InitializeInModes();
    InitializeOutModes();
    RecomputeDampProfile();

    InitializeStates();
    ResetMatrices();
}

ModalStiffStringProcessor::~ModalStiffStringProcessor()
{
}

//==========================================================================
void ModalStiffStringProcessor::SetTimeStep(double aTimeStep)
{
    bool vCurrPlayState = mPlayState;
    if (vCurrPlayState)
    {
        mPlayState.store(false);
    }
    mTimeStep = aTimeStep;
    if (vCurrPlayState)
    {
        mPlayState.store(true);
    }
}

void ModalStiffStringProcessor::SetPlayState(bool aPlayState)
{
    mPlayState.store(aPlayState);
}

void ModalStiffStringProcessor::ResetStringStates()
{
    if (mPlayState.load())
    {
        mPlayState.store(false);
    }
    std::fill(mStates[0].begin(), mStates[0].end(), 0);
    std::fill(mStates[1].begin(), mStates[1].end(), 0);

    mPreviousSample = 0.0;
}

void ModalStiffStringProcessor::SetInputPos(float aNewPos)
{
    //Position is in normalized percentage of string length
    if (aNewPos >= 0 && aNewPos <= 1)
    {
        mExcitPos = aNewPos * mLength;
    }
    else
    {
        jassertfalse;
    }
    RecomputeInModes();
}

void ModalStiffStringProcessor::SetReadPos(float aNewPos)
{
    //Position is in normalized percentage of string length
    if (aNewPos >= 0 && aNewPos <= 1)
    {
        mReadPos = aNewPos * mLength;
    }
    else
    {
        jassertfalse;
    }
    RecomputeOutModes();
}

void ModalStiffStringProcessor::SetGain(float aGain)
{
    mGain.store(aGain);
}

void ModalStiffStringProcessor::SetBowPressure(float aPressure)
{
    mFb.store(aPressure);
}

void ModalStiffStringProcessor::SetBowSpeed(float aSpeed)
{
    mVb.store(aSpeed);
}

void ModalStiffStringProcessor::SetString(Global::Strings::String* apString)
{
    auto vPi = juce::MathConstants<float>::pi;

    mpString = apString;
    mRadius = mpString->mRadius;
    mDensity = mpString->mDensity;
    mTension = mpString->mTension;
    mLength = mpString->mLength;

    mArea = vPi * mRadius * mRadius;
    mLinDensity = mDensity * mArea;
    mInertia = (vPi * mRadius * mRadius * mRadius * mRadius) / 4;
    mK = sqrt(mYoungMod * mInertia / (mLinDensity * mLength * mLength * mLength * mLength));
    mC = sqrt(mTension / mLinDensity);

    ResetStringStates();
    RecomputeModesNumber();
    RecomputeEigenFreqs();
    InitializeInModes();
    InitializeOutModes();
    RecomputeDampProfile();
    InitializeStates();
    ResetMatrices();
}

void ModalStiffStringProcessor::ComputeState()
{
    if (mPlayState.load())
    {
        //for (int vOS = 0; vOS < mOversamplingFactor; ++vOS)
        //{
            //Computing input projection
            float vZeta1 = 0.f;
            for (int i = 0; i < mModesNumber; ++i)
            {
                vZeta1 += mpModesInCurr.load()[i] * mpStatesPtrs[0][i + mModesNumber];
            }

            //Computing bow input
            float vEta = vZeta1 - mVb.load();
            float vD = sqrt(2 * mA) * exp(-mA * vEta * vEta + 0.5);
            float vLambda = vD * (1 - 2 * mA * vEta * vEta);

            float vVt1 = 0.f;
            float vVt2 = 0.f;

            //Computing known terms
            for (int i = 0; i < mModesNumber; ++i)
            {
                float vZeta2 = mpModesInCurr.load()[i] * vZeta1;

                //Notice that the first half of zeta in the matlab code is made of zeroes, 
                //so there is no point of computing multiplications by it
                float vB1 = mB11[i] * mpStatesPtrs[0][i] + mB12[i] * mpStatesPtrs[0][i + mModesNumber];
                float vB2 = mB21[i] * mpStatesPtrs[0][i] + mB22[i] * mpStatesPtrs[0][i + mModesNumber] +
                    vZeta2 * 0.5f * mTimeStep * mFb.load() * (vLambda - 2 * vD) +
                    mTimeStep * mFb.load() * vD * mpModesInCurr.load()[i] * mVb.load();

                //Computing T^-1*a (see overleaf notes)
                float vZ1 = 0.5f * mTimeStep * mFb.load() * vLambda * mpModesInCurr.load()[i];
                mInvAv2[i] = (1 / mSchurComp[i]) * vZ1;
                mInvAv1[i] = -mT11[i] * mT12[i] * mInvAv2[i];

                //Computing T^-1*[j1;j1] (see overleaf notes)
                float vY2 = mT11[i] * vB1;
                float vZ2 = vB2 - mT21[i] * vY2;
                mInvAb2[i] = (1 / mSchurComp[i]) * vZ2;
                mInvAb1[i] = vY2 - mT11[i] * mT12[i] * mInvAb2[i];

                vVt1 += mpModesInCurr.load()[i] * mInvAv2[i];
                vVt2 += mpModesInCurr.load()[i] * mInvAb2[i];
            }

            float vCoeff = 1 / (1 + vVt1);

            for (int i = 0; i < mModesNumber; ++i)
            {
                mpStatesPtrs[1][i] = mInvAb1[i] - vCoeff * mInvAv1[i] * vVt2;
                mpStatesPtrs[1][i + mModesNumber] = mInvAb2[i] - vCoeff * mInvAv2[i] * vVt2;
            }

            //Pointers switch
            auto vpStatePointer = mpStatesPtrs[0];
            mpStatesPtrs[0] = mpStatesPtrs[1];
            mpStatesPtrs[1] = vpStatePointer;
        //}
    }
}

float ModalStiffStringProcessor::ReadOutput()
{
    float vOutputValue = 0.f;
    if (mPlayState.load())
    {
        for (int i = 0; i < mModesNumber; ++i)
        {
            vOutputValue += mpModesOutCurr.load()[i] * mpStatesPtrs[0][i];
        }
    }

    float vDiffOutput = mGain.load() * (vOutputValue - mPreviousSample) / mTimeStep;
    mPreviousSample = vOutputValue;
    return vDiffOutput;
}

int ModalStiffStringProcessor::GetModesNumber()
{
    return mModesNumber;
}

void ModalStiffStringProcessor::GetModesAtLocation(std::vector<float>& aModesArray, float aLocationFrac)
{
    aModesArray.resize(mModesNumber, 0.f);
    auto vPos = aLocationFrac * mLength;
    for (int i = 0; i < mModesNumber; ++i)
    {
        auto vMode = ComputeMode(vPos, i + 1);
        aModesArray[i] = vMode;
    }
}

std::vector<float> ModalStiffStringProcessor::GetStringState()
{
    std::vector<float> vState(mModesNumber);
    std::copy(mStates[0].begin(), mStates[0].begin() + mModesNumber, vState.begin());
    return vState;
}

//==========================================================================
float ModalStiffStringProcessor::ComputeEigenFreq(int aModeNumber)
{
    auto vN = aModeNumber * juce::MathConstants<float>::pi / mLength;
    return sqrt((mTension / mLinDensity) * vN * vN + (mYoungMod * mInertia / mLinDensity) * vN * vN * vN * vN);
}

float ModalStiffStringProcessor::ComputeMode(float aPos, int aModeNumber)
{
    return sqrt(2 / mLength) * sin(aModeNumber * juce::MathConstants<float>::pi * aPos / mLength);
}

float ModalStiffStringProcessor::ComputeDampCoeff(float aFreq)
{
    auto vPi = juce::MathConstants<float>::pi;
    float vRhoAir = 1.225f;
    float vMuAir = (float)1.619e-5;
    auto vD0 = -2 * vRhoAir * vMuAir / (mDensity * mRadius * mRadius);
    auto vD1 = -2 * vRhoAir * sqrt(2 * vMuAir) / (mDensity * mRadius);
    auto vD2 = static_cast<float>(-1 / 18000);
    auto vD3 = -0.003f * mYoungMod * mDensity * vPi * vPi * mRadius * mRadius * mRadius * mRadius * mRadius * mRadius / (4 * mTension * mTension);
    return vD0 + vD1 * sqrt(aFreq) + vD2 * aFreq + vD3 * aFreq * aFreq * aFreq;
}

void ModalStiffStringProcessor::RecomputeModesNumber()
{
    int vModesNumber = 1;
    float vLimitFreq = 20e3 * 2 * juce::MathConstants<float>::pi;
    while (true)
    {
        auto vFreq = ComputeEigenFreq(vModesNumber);
        if (vFreq > vLimitFreq) 
        {
            --vModesNumber;
            break;
        }
        ++vModesNumber;
    }
    mModesNumber = vModesNumber;    
}

void ModalStiffStringProcessor::RecomputeEigenFreqs()
{
    mEigenFreqs.resize(mModesNumber);
    for (int i = 0; i < mModesNumber; ++i)
    {
        mEigenFreqs[i] = ComputeEigenFreq(i + 1);
    }
}

void ModalStiffStringProcessor::InitializeInModes()
{
    mModesIn.resize(2);
    mModesIn = std::vector<std::vector<float>>(2, std::vector<float>(mModesNumber, 0));

    mpModesInCurr.store(&mModesIn[0][0]);
    mpModesInNew.store(&mModesIn[1][0]);

    RecomputeInModes();
}

void ModalStiffStringProcessor::InitializeOutModes()
{
    mModesOut.resize(2);
    mModesOut = std::vector<std::vector<float>>(2, std::vector<float>(mModesNumber, 0));

    mpModesOutCurr.store(&mModesOut[0][0]);
    mpModesOutNew.store(&mModesOut[1][0]);

    RecomputeOutModes();
}

void ModalStiffStringProcessor::RecomputeInModes()
{
    //Computing new modes offline on another thread
    for (int i = 0; i < mModesNumber; ++i)
    {
        mpModesInNew.load()[i] = ComputeMode(mExcitPos, i + 1);
    }
    //Atomic pointers switch allows to change position online
    auto vpModesPtr = mpModesInCurr.load();
    mpModesInCurr.store(mpModesInNew.load());
    mpModesInNew.store(vpModesPtr);
}

void ModalStiffStringProcessor::RecomputeOutModes()
{
    //Computing new modes offline on another thread
    for (int i = 0; i < mModesNumber; ++i)
    {
        mpModesOutNew.load()[i] = ComputeMode(mReadPos, i + 1);
    }
    //Atomic pointers switch allows to change position online
    auto vpModesPtr = mpModesOutCurr.load();
    mpModesOutCurr.store(mpModesOutNew.load());
    mpModesOutNew.store(vpModesPtr);
}

void ModalStiffStringProcessor::RecomputeDampProfile()
{
    mDampCoeffs.resize(mModesNumber);
    for (int i = 0; i < mModesNumber; ++i)
    {
        auto vFreq = mEigenFreqs[i];
        mDampCoeffs[i] = - ComputeDampCoeff(vFreq);
    }
}

void ModalStiffStringProcessor::InitializeStates()
{
    // Initialise xVectors
    mStates.resize(2);
    mpStatesPtrs.resize(2);
    // initialise states container with two vectors of 0s
    mStates = std::vector<std::vector<float>>(2, std::vector<float>(mModesNumber * 2, 0));
    mpStatesPtrs = std::vector<float*>(2, nullptr);
    // initialise pointers to state vectors
    for (int i = 0; i < 2; ++i)
    {
        mpStatesPtrs[i] = &mStates[i][0];
    }
}

void ModalStiffStringProcessor::ResetMatrices()
{
    mT11.resize(mModesNumber);
    mT12.resize(mModesNumber);
    mT21.resize(mModesNumber);
    mT22.resize(mModesNumber);
    mSchurComp.resize(mModesNumber);

    mB11.resize(mModesNumber);
    mB12.resize(mModesNumber);
    mB21.resize(mModesNumber);
    mB22.resize(mModesNumber);

    mZeta2.resize(mModesNumber);
    mB1.resize(mModesNumber);
    mB2.resize(mModesNumber);

    mZ1.resize(mModesNumber);
    mInvAv2.resize(mModesNumber);
    mInvAv1.resize(mModesNumber);

    mY2.resize(mModesNumber);
    mZ2.resize(mModesNumber);
    mInvAb2.resize(mModesNumber);
    mInvAb1.resize(mModesNumber);

    for (int i = 0; i < mModesNumber; ++i)
    {
        mT11[i] = 1;
        mT12[i] = - 0.5f * mTimeStep;
        mT21[i] = - 0.5f * mTimeStep * (-mEigenFreqs[i] * mEigenFreqs[i]);
        mT22[i] = 1 + 0.5 * mTimeStep * mDampCoeffs[i];

        mSchurComp[i] = mT22[i] - mT21[i] * (mT11[i] * mT12[i]);

        mB11[i] = 1;
        mB12[i] = 0.5f * mTimeStep;
        mB21[i] = 0.5f * mTimeStep * (-mEigenFreqs[i] * mEigenFreqs[i]);
        mB22[i] = 1 - 0.5 * mTimeStep * mDampCoeffs[i];
    }
}
