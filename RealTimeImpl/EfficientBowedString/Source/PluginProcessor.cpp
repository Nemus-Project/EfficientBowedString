/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
FastBowedStringAudioProcessor::FastBowedStringAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
}

FastBowedStringAudioProcessor::~FastBowedStringAudioProcessor()
{
}

//==============================================================================
const juce::String FastBowedStringAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool FastBowedStringAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool FastBowedStringAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool FastBowedStringAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double FastBowedStringAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int FastBowedStringAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int FastBowedStringAudioProcessor::getCurrentProgram()
{
    return 0;
}

void FastBowedStringAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String FastBowedStringAudioProcessor::getProgramName (int index)
{
    return {};
}

void FastBowedStringAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void FastBowedStringAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    

    if (!mpModalStiffStringProcessor)
    {
        mpModalStiffStringProcessor = std::make_shared<ModalStiffStringProcessor>(sampleRate, Global::Strings::kpCelloG2.get());
    }
    else if (mSampleRate != sampleRate)
    {
        mpModalStiffStringProcessor->SetTimeStep(1.0 / sampleRate);
    }   


    // save samplerate and block size
    mSampleRate = sampleRate;
    mBlockSize = samplesPerBlock;

    if (!mpLPFilter)
    {
        mpLPFilter.reset(new PA_LowPass2());
    }
}

void FastBowedStringAudioProcessor::releaseResources()
{
    /*delete(Global::Strings::kpCelloA3);
    delete(Global::Strings::kpCelloD3);
    delete(Global::Strings::kpCelloG2);
    delete(Global::Strings::kpCelloC2);
    delete(Global::Strings::kpBassG2);
    delete(Global::Strings::kpBassD2);
    delete(Global::Strings::kpBassA1);
    delete(Global::Strings::kpBassE1);*/
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool FastBowedStringAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void FastBowedStringAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    //// In case we have more outputs than inputs, this code clears any output
    //// channels that didn't contain input data, (because these aren't
    //// guaranteed to be empty - they may contain garbage).
    //// This is here to avoid people getting screaming feedback
    //// when they first compile a plugin, but obviously you don't need to keep
    //// this code if your algorithm always overwrites all the output channels.
    //for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
    //    buffer.clear (i, 0, buffer.getNumSamples());

    float vOutput = 0.0;

    // Get the current time
    //double vNow = Time::getMillisecondCounterHiRes();

    for (int i = 0; i < buffer.getNumSamples(); ++i)
    {
        mpModalStiffStringProcessor->ComputeState();
        //vOutput = mpLPFilter->update(mpModalStiffStringProcessor->ReadOutput());
        vOutput = mpModalStiffStringProcessor->ReadOutput();
        //DBG(Global::limitOutput(vOutput));
        //jassert(vOutput <= 1 && vOutput >= -1);
        for (int channel = 0; channel < totalNumOutputChannels; ++channel)
        {
            auto vpBuffer = buffer.getWritePointer(channel);
            vpBuffer[i] = Global::limitOutput(vOutput);
        }
    }

    /*mCumulativeTimePerBufferMod = (Time::getMillisecondCounterHiRes() - vNow);
    float vRealTime = (1000 / mSampleRate) * mBlockSize;

    Logger::getCurrentLogger()->outputDebugString("3: Modal: " + String(mCumulativeTimePerBufferMod/vRealTime));*/

}

//==============================================================================
bool FastBowedStringAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* FastBowedStringAudioProcessor::createEditor()
{
    return new FastBowedStringAudioProcessorEditor (*this);
}

//==============================================================================
void FastBowedStringAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void FastBowedStringAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

std::shared_ptr<ModalStiffStringProcessor> FastBowedStringAudioProcessor::GetModalStringProcessor()
{
    return mpModalStiffStringProcessor;
}


//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new FastBowedStringAudioProcessor();
}
