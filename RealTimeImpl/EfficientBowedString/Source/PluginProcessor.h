/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Global.h"
#include "ModalStiffStringProcessor.h"
#include "PA_LowPass2.h"

//==============================================================================
/**
*/
class FastBowedStringAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    FastBowedStringAudioProcessor();
    ~FastBowedStringAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    //==============================================================================
    std::shared_ptr<ModalStiffStringProcessor> GetModalStringProcessor();

private:
    //==============================================================================
    double mSampleRate{ 0.0 }; // sample rate
    int mBlockSize{ 0 };
    

    std::shared_ptr<ModalStiffStringProcessor> mpModalStiffStringProcessor;

    std::unique_ptr<PA_LowPass2> mpLPFilter;
    
    // Current sample (debugging purposes only)
    unsigned long curSample = 0;
    
    // Current buffer (debugging purposes only)
    unsigned long curBuffer = 0;

    // Cumulative time for calculating the reference, optimised matrix, or optimised vector method (debugging purposes only)
    double cumulativeTimePerBufferRef = 0;
    double cumulativeTimePerBufferOpt = 0;
    double cumulativeTimePerBufferOptVec = 0;

    double mCumulativeTimePerBufferMod{ 0.0 };
   

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (FastBowedStringAudioProcessor)
};
