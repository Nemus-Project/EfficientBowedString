#include "ModalStiffStringView.h"

//==========================================================================
ModalStiffStringView::ModalStiffStringView()
{
	mGainInitialValue = 0.f;
	mInputPosInitialValue = 0.833f;
	mReadPosInitialValue = 0.33f;
	mBowPressureInitialValue = 5.f;
	mBowSpeedInitialValue = 0.2f;

	mPlayButton.addListener(this);
	mResetButton.addListener(this);
	mGainSlider.addListener(this);
	mInputPosSlider.addListener(this);
	mReadPosSlider.addListener(this);
	mBowPressureSlider.addListener(this);
	mBowSpeedSlider.addListener(this);
	mStringChoiceBox.addListener(this);

	mGainSlider.setRange(0.0, 1, 0.001);
	mGainSlider.setValue(mGainInitialValue, juce::sendNotification);
	mGainSlider.setSliderStyle(juce::Slider::SliderStyle::RotaryVerticalDrag);
	mGainSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 15);

	addAndMakeVisible(mGainLabel);
	mGainLabel.setText("Gain", juce::dontSendNotification);
	mGainLabel.attachToComponent(&mGainSlider, false);
	mGainLabel.setJustificationType(juce::Justification::centred);

	mInputPosSlider.setRange(0.0, 100.0, 0.01);
	mInputPosSlider.setValue(mInputPosInitialValue * 100.0, juce::sendNotification);
	mInputPosSlider.setSliderStyle(juce::Slider::SliderStyle::LinearHorizontal);
	mInputPosSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 15);
	mInputPosSlider.setTextValueSuffix(" %");

	mReadPosSlider.setRange(0.0, 100.0, 0.01);
	mReadPosSlider.setValue(mReadPosInitialValue * 100.0, juce::sendNotification);
	mReadPosSlider.setSliderStyle(juce::Slider::SliderStyle::LinearHorizontal);
	mReadPosSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 15);
	mReadPosSlider.setTextValueSuffix(" %");

	addAndMakeVisible(mInputPosLabel);
	mInputPosLabel.setText("Input Pos", juce::dontSendNotification);
	mInputPosLabel.attachToComponent(&mInputPosSlider, false);
	mInputPosLabel.setJustificationType(juce::Justification::centred);

	addAndMakeVisible(mReadPosLabel);
	mReadPosLabel.setText("Read Pos", juce::dontSendNotification);
	mReadPosLabel.attachToComponent(&mReadPosSlider, false);
	mReadPosLabel.setJustificationType(juce::Justification::centred);

	mBowPressureSlider.setRange(0.0, 100.0, 0.001);
	mBowPressureSlider.setValue(mBowPressureInitialValue, juce::sendNotification);
	mBowPressureSlider.setSliderStyle(juce::Slider::SliderStyle::LinearHorizontal);
	mBowPressureSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 15);

	addAndMakeVisible(mBowPressureLabel);
	mBowPressureLabel.setText("Bow Pressure", juce::dontSendNotification);
	mBowPressureLabel.attachToComponent(&mBowPressureSlider, false);
	mBowPressureLabel.setJustificationType(juce::Justification::centred);

	mBowSpeedSlider.setRange(0.0, 2.0, 0.0001);
	mBowSpeedSlider.setValue(mBowSpeedInitialValue, juce::sendNotification);
	mBowSpeedSlider.setSliderStyle(juce::Slider::SliderStyle::LinearHorizontal);
	mBowSpeedSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 50, 15);
	mBowSpeedSlider.setTextValueSuffix("m/s");

	addAndMakeVisible(mBowSpeedLabel);
	mBowSpeedLabel.setText("Bow Speed", juce::dontSendNotification);
	mBowSpeedLabel.attachToComponent(&mBowSpeedSlider, false);
	mBowSpeedLabel.setJustificationType(juce::Justification::centred);

	addAndMakeVisible(mStringChoiceBox);
	mStringChoiceBox.addItem(Global::Strings::kpCelloA3->mName, 1);
	mStringChoiceBox.addItem(Global::Strings::kpCelloD3->mName, 2);
	mStringChoiceBox.addItem(Global::Strings::kpCelloG2->mName, 3);
	mStringChoiceBox.addItem(Global::Strings::kpCelloC2->mName, 4);
	mStringChoiceBox.addItem(Global::Strings::kpBassG2->mName, 5);
	mStringChoiceBox.addItem(Global::Strings::kpBassD2->mName, 6);
	mStringChoiceBox.addItem(Global::Strings::kpBassA1->mName, 7);
	mStringChoiceBox.addItem(Global::Strings::kpBassE1->mName, 8);
	mStringChoiceBox.setSelectedId(3, juce::sendNotification);
}

ModalStiffStringView::~ModalStiffStringView()
{
}

//==========================================================================
void ModalStiffStringView::paint(juce::Graphics& g)
{

	std::vector<float> vCurStringState = mpStiffStringProcessor->GetStringState();

	juce::Path vStringPath = VisualiseState(g, vCurStringState);
    g.setColour (Colours::cyan);
    g.strokePath (vStringPath, PathStrokeType(2.0f));
	g.fillPath(DrawInPointer(g, vStringPath));
	g.fillPath(DrawOutPointer(g, vStringPath));
}

juce::Path ModalStiffStringView::VisualiseState (juce::Graphics& g, std::vector<float> aStringState)
{
	//amplifies string movement for visualization purposes (string displacement is too tiny)
    double vVisualScaling = 5.0;
    
    // String-boundaries are in the vertical middle of the component
    double vStringBoundaries = (getHeight() / 2.0) - (getHeight() / 4.0);
    
    // Initialise path
    juce::Path vStringPath;
    
    // Start path
	vStringPath.startNewSubPath (0, vStringBoundaries);
    
    double vSpacing = getWidth() / static_cast<double>(mVisualizationPoints);
    double vX = vSpacing;
    
    std::vector<float> vPointStates (mVisualizationPoints, 0);

	for (int l = 0; l < mVisualizationPoints; ++l)
	{
		for (int m = 0; m < mStringModesNumber; ++m)
		{
			vPointStates[l] += aStringState[m] * mVisualizationModes[l][m];
		}
	}
        

    for (int l = 0; l < mVisualizationPoints; l++)
    {
        // Needs to be -u, because a positive u would visually go down
        float vNewY = -vPointStates[l] * vVisualScaling * getHeight() + vStringBoundaries;
        
        // if we get NAN values, make sure that we don't get an exception
        if (isnan(vNewY))
			vNewY = 0;
        
		vStringPath.lineTo (vX, vNewY);
		vX += vSpacing;
    }

	//juce::Path vInPointerPath = DrawInPointer(g, vStringPath);
	//juce::Path vOutPointerPath = DrawOutPointer(g, vStringPath);
	//vStringPath.addPath(vInPointerPath);
	//vStringPath.addPath(vOutPointerPath);

    return vStringPath;
}

juce::Path ModalStiffStringView::DrawInPointer(juce::Graphics& g, juce::Path aStringPath)
{
	// String-boundaries are in the vertical middle of the component
	double vStringBoundaries = (getHeight() / 2.0) - (getHeight() / 4.0);
	auto vInputPos = getWidth() * mInputPosSlider.getValue() / 100.f;

	juce::Line<float> vInLine(vInputPos, getHeight(), vInputPos, 0.f);
	juce::Array<juce::Point<float>> vPossibleIntersections;

	bool vIntersects = IntersectsPath(aStringPath, vInLine, vPossibleIntersections);

	//auto vInPoint = aStringPath.getPointAlongPath(vInputPos);

	float vPointerSize = 10.f;
	juce::Path vInPointerPath;
	//juce::Rectangle<float> vInRectangle(vInPoint.getX() - vPointerSize / 2, vInPoint.getY() - vPointerSize / 2, vPointerSize, vPointerSize);
	juce::Rectangle<float> vInRectangle(
		vPossibleIntersections[0].getX() - vPointerSize / 2, 
		vPossibleIntersections[0].getY() - vPointerSize / 2, 
		vPointerSize, vPointerSize);

	vInPointerPath.addEllipse(vInRectangle);

	return vInPointerPath;
}

juce::Path ModalStiffStringView::DrawOutPointer(juce::Graphics& g, juce::Path aStringPath)
{
	// String-boundaries are in the vertical middle of the component
	double vStringBoundaries = (getHeight() / 2.0) - (getHeight() / 4.0);
	auto vReadPos = getWidth() * mReadPosSlider.getValue() / 100.f;

	juce::Line<float> vInLine(vReadPos, getHeight(), vReadPos, 0.f);
	juce::Array<juce::Point<float>> vPossibleIntersections;

	bool vIntersects = IntersectsPath(aStringPath, vInLine, vPossibleIntersections);

	//auto vReadPoint = aStringPath.getPointAlongPath(vReadPos);

	float vPointerSize = 10.f;
	//juce::Rectangle<float> vOutRectangle(vReadPoint.getX() - vPointerSize / 2, vReadPoint.getY() - vPointerSize / 2, vPointerSize, vPointerSize);
	juce::Rectangle<float> vOutRectangle(
		vPossibleIntersections[0].getX() - vPointerSize / 2, 
		vPossibleIntersections[0].getY() - vPointerSize / 2, 
		vPointerSize, vPointerSize);
	juce::AffineTransform vTransform = juce::AffineTransform::rotation(
		-juce::MathConstants<float>::pi / 4.f, 
		vPossibleIntersections[0].getX(), vPossibleIntersections[0].getY());

	juce::Path vReadPointerPath;
	vReadPointerPath.addRectangle(vOutRectangle);
	vReadPointerPath.applyTransform(vTransform);

	return vReadPointerPath;
}

void ModalStiffStringView::resized()
{
	int vButtonsWidth = 100;
	int vButtonHeigth = 30;
	int vGainSliderDims = 100;
	int vSpacing = 50;

	mPlayButton.setButtonText("PLAY/PAUSE");
	mPlayButton.setBounds(getWidth() / 2 - getWidth() / 4 - vButtonsWidth / 2, getHeight() - getHeight() / 4 + vGainSliderDims / 2 - vButtonHeigth, vButtonsWidth, vButtonHeigth);
	mPlayButton.setClickingTogglesState(true);

	mResetButton.setBounds(getWidth() / 2 + getWidth() / 4 - vButtonsWidth / 2, getHeight() - getHeight() / 4 + vGainSliderDims / 2 - vButtonHeigth, vButtonsWidth, vButtonHeigth);
	mResetButton.setButtonText("RESET STATE");
	mResetButton.setClickingTogglesState(false);

	mInputPosSlider.setBounds(getWidth() / 2 - getWidth() / 4 - vButtonsWidth / 2, getHeight() - getHeight() / 4 - vSpacing - vButtonHeigth, vButtonsWidth, vButtonHeigth);
	mReadPosSlider.setBounds(getWidth() / 2 + getWidth() / 4 - vButtonsWidth / 2, getHeight() - getHeight() / 4 - vSpacing - vButtonHeigth, vButtonsWidth, vButtonHeigth);

	mBowPressureSlider.setBounds(getWidth() / 2 - getWidth() / 8 - vButtonsWidth / 2, getHeight() - getHeight() / 4 - vSpacing - vButtonHeigth, vButtonsWidth, vButtonHeigth);
	mBowSpeedSlider.setBounds(getWidth() / 2 + getWidth() / 8 - vButtonsWidth / 2, getHeight() - getHeight() / 4 - vSpacing - vButtonHeigth, vButtonsWidth, vButtonHeigth);

	mGainSlider.setBounds(getWidth() / 2 - vGainSliderDims / 2, getHeight() - (getHeight() / 5) * 2, vGainSliderDims, vGainSliderDims);

	mStringChoiceBox.setBounds(getWidth() / 2 - vButtonsWidth / 2, getHeight() - getHeight() / 4 + vGainSliderDims / 2 - vButtonHeigth, vButtonsWidth, vButtonHeigth);

	addAndMakeVisible(mPlayButton);
	addAndMakeVisible(mResetButton);
	addAndMakeVisible(mGainSlider);
	addAndMakeVisible(mInputPosSlider);
	addAndMakeVisible(mReadPosSlider);
	addAndMakeVisible(mBowPressureSlider);
	addAndMakeVisible(mBowSpeedSlider);
}

void ModalStiffStringView::buttonClicked(juce::Button* apButton)
{
	if (apButton == &mPlayButton)
	{
		mpStiffStringProcessor->SetPlayState(apButton->getToggleState());
	}
	else if (apButton == &mResetButton)
	{
		if (mPlayButton.getToggleState())
		{
			mPlayButton.setToggleState(false, juce::sendNotification);
		}
		mpStiffStringProcessor->ResetStringStates();

		mGainSlider.setValue(mGainInitialValue, juce::sendNotification);
		mInputPosSlider.setValue(mInputPosInitialValue * 100.f, juce::sendNotification);
		mReadPosSlider.setValue(mReadPosInitialValue * 100.f, juce::sendNotification);
		mBowPressureSlider.setValue(mBowPressureInitialValue, juce::sendNotification);
		mBowSpeedSlider.setValue(mBowSpeedInitialValue, juce::sendNotification);
	}
}

void ModalStiffStringView::sliderValueChanged(juce::Slider* apSlider)
{
	if (apSlider == &mGainSlider)
	{
		mpStiffStringProcessor->SetGain(static_cast<float>(mGainSlider.getValue()));
	}
	else if (apSlider == &mInputPosSlider)
	{
		//Value is in percentage if string length
		auto vValue = juce::jlimit<float>(0.f, 1.f, mInputPosSlider.getValue() / 100.0);
		mpStiffStringProcessor->SetInputPos(vValue);
		mpStiffStringProcessor->GetModesAtLocation(mInputModes, vValue);
	}
	else if (apSlider == &mReadPosSlider)
	{
		//Value is in percentage if string length
		auto vValue = juce::jlimit<float>(0.f, 1.f, mReadPosSlider.getValue() / 100.0);
		mpStiffStringProcessor->SetReadPos(vValue);
		mpStiffStringProcessor->GetModesAtLocation(mReadModes, vValue);
	}
	else if (apSlider == &mBowPressureSlider)
	{
		mpStiffStringProcessor->SetBowPressure(static_cast<float>(apSlider->getValue()));
	}
	else if (apSlider == &mBowSpeedSlider)
	{
		mpStiffStringProcessor->SetBowSpeed(static_cast<float>(apSlider->getValue()));
	}
}

void ModalStiffStringView::comboBoxChanged(juce::ComboBox* comboBoxThatHasChanged)
{
	if (comboBoxThatHasChanged == &mStringChoiceBox)
	{ 
		if (mPlayButton.getToggleState())
		{
			mPlayButton.setToggleState(false, juce::sendNotification);
		}
		if (mStringChoiceBox.getSelectedId() == Global::Strings::kpCelloA3->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpCelloA3.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpCelloD3->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpCelloD3.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpCelloG2->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpCelloG2.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpCelloC2->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpCelloC2.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpBassG2->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpBassG2.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpBassD2->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpBassD2.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpBassA1->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpBassA1.get());
		}
		else if (mStringChoiceBox.getSelectedId() == Global::Strings::kpBassE1->mId)
		{
			mpStiffStringProcessor->SetString(Global::Strings::kpBassE1.get());
		}
		mStringModesNumber = mpStiffStringProcessor->GetModesNumber();
		SetVisualizationModes();
	}
}

void ModalStiffStringView::SetProcessor(std::shared_ptr<ModalStiffStringProcessor> apProcessor)
{
	jassert(apProcessor);
	if (apProcessor)
	{
		mpStiffStringProcessor = apProcessor;
	}
	mStringModesNumber = mpStiffStringProcessor->GetModesNumber();
	SetVisualizationModes();
}

void ModalStiffStringView::SetVisualizationModes()
{
	mVisualizationPoints = 101;
	float vStep = 1.f / (mVisualizationPoints - 1); //position expressed in normalized percentage
	mVisualizationModes = std::vector<std::vector<float>>(mVisualizationPoints, std::vector<float>(mStringModesNumber, 0));
	for (int i = 0; i < mVisualizationPoints; ++i)
	{
		mpStiffStringProcessor->GetModesAtLocation(mVisualizationModes[i], vStep * i);
	}

	mInputModes.resize(mStringModesNumber, 0.f);
	mReadModes.resize(mStringModesNumber, 0.f);
}

bool ModalStiffStringView::IntersectsPath(const juce::Path& p, juce::Line<float> line, juce::Array<juce::Point<float>>& possibleIntersections)
{
	juce::PathFlatteningIterator i(p, juce::AffineTransform(), juce::Path::defaultToleranceForTesting);

	juce::Point<float> intersectionPoint;
	while (i.next())
	{
		if (line.intersects(juce::Line<float>(i.x1, i.y1, i.x2, i.y2), intersectionPoint))
		{
			possibleIntersections.add(intersectionPoint);
		}
	}

	if (possibleIntersections.isEmpty())
		return false;

	return true;
}