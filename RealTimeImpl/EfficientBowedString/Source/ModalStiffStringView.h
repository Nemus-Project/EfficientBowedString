/*
  ==============================================================================

    ModalStiffStringView.h
    Created: 04/05/2022
    Author:  Riccardo Russo

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "ModalStiffStringProcessor.h"

class ModalStiffStringView
    : public juce::Component
    , public juce::Button::Listener
    , public juce::Slider::Listener
    , public juce::ComboBox::Listener
{
public:
    ModalStiffStringView();
    ~ModalStiffStringView();

    //==========================================================================
    // juce::Component
    void paint(juce::Graphics&) override;
    void resized() override;

    //==========================================================================
    void buttonClicked(juce::Button* apButton) override;
    void sliderValueChanged(juce::Slider* apSlider) override;
    void comboBoxChanged(juce::ComboBox* comboBoxThatHasChanged) override;

    //==========================================================================
    void SetProcessor(std::shared_ptr<ModalStiffStringProcessor> apProcessor);

private:
    std::shared_ptr<ModalStiffStringProcessor> mpStiffStringProcessor;

    bool mPlayState{ false };

    juce::TextButton mPlayButton;
    juce::TextButton mResetButton;

    juce::Slider mGainSlider;
    juce::Label mGainLabel;

    juce::Slider mInputPosSlider;
    juce::Slider mReadPosSlider;
    juce::Label mInputPosLabel;
    juce::Label mReadPosLabel;

    juce::Slider mBowPressureSlider;
    juce::Label mBowPressureLabel;

    juce::Slider mBowSpeedSlider;
    juce::Label mBowSpeedLabel;

    juce::ComboBox mStringChoiceBox;

    int mVisualizationPoints{ 0 };
    std::vector<std::vector<float>> mVisualizationModes;
    std::vector<float> mInputModes;
    std::vector<float> mReadModes;

    int mStringModesNumber{ 0 };
    void SetVisualizationModes();
    juce::Path VisualiseState(juce::Graphics& g, std::vector<float> aStringState);
    juce::Path DrawInPointer(juce::Graphics& g, juce::Path aStringPath);
    juce::Path DrawOutPointer(juce::Graphics& g, juce::Path aStringPath);

    /*
    * Found: https://forum.juce.com/t/path-intersectsline-line-float-line-float-tolerance-point-float-intersection/32283
    */
    bool IntersectsPath(const juce::Path& p, juce::Line<float> line, juce::Array<juce::Point<float>>& possibleIntersections);
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModalStiffStringView)
};

