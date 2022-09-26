/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
FastBowedStringAudioProcessorEditor::FastBowedStringAudioProcessorEditor (FastBowedStringAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{   
    mpModalStiffString = std::make_unique<ModalStiffStringView>();
    addAndMakeVisible(*mpModalStiffString);
    mpModalStiffString->SetProcessor(p.GetModalStringProcessor());

    // Refresh the graphics at a rate of 15 Hz
    startTimerHz (15);

    setSize (800, 600);
}

FastBowedStringAudioProcessorEditor::~FastBowedStringAudioProcessorEditor()
{
}

//==============================================================================
void FastBowedStringAudioProcessorEditor::paint (juce::Graphics& g)
{
    
}

void FastBowedStringAudioProcessorEditor::resized()
{
    mpModalStiffString->setBounds(getLocalBounds());
}

void FastBowedStringAudioProcessorEditor::timerCallback()
{
    repaint();
}
