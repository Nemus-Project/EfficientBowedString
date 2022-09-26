/*
  ==============================================================================

    Bowed1DWaveFirstOrder.h
    Created: 25 Apr 2022 12:20:12pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Global.h"

//==============================================================================
/*
*/
class Bowed1DWaveFirstOrder : public juce::Component
{
public:
    Bowed1DWaveFirstOrder (double k);
    ~Bowed1DWaveFirstOrder() override;

    void paint (juce::Graphics&) override;
    
    // Function to draw the state of the system
    Path visualiseState (Graphics& g, double visualScaling, double* x, double offset);
    
    void calculateFirstOrderRef(); // Reference first order system calculation
    void calculateFirstOrderOptVec(); // Optimised first order system calculation with vectors
    void calculateFirstOrderOpt(); // Optimised first order system calculation

    float getOutput (float outRatio) { return xVec[1][N + (int)floor(outRatio * N)]; };
    
    // Function to check whether the optimised version is equal to the reference (within machine precision (here considered to be < 1e-10))
    double getDiffSum();
    
private:
    // Recalculate the zeta vector
    void recalculateZeta();
    
    // Time step
    double k;
    
    // Scheme parameters (length, wave speed and grid spacing)
    double L, c, h;
    
    // Number of grid points
    int N;
    
    // Length of x vector (2*N-1 for first-order, 2 * (N-1) for modal)
    int NN;
    
    // Bowing variables
    double a;   // free parameter
    double xB;  // Bowing location (as a ratio of the length) (can be made mouse-controlled)
    double vB;  // Bowing velocity (in m/s)
    double Fb;  // Bowing force (in m^2/s^2) (?)
    double eta; // relative velocity between the string and bow (in m/s)
    
    double lambda, d; // noniterative factors
    double outPos; // output location (as a ratio of the length)
        
    // Vectors and matrices (eigen library)
    Eigen::VectorXd xNext, x, xNextRef, xRef, xPaint, xRefPaint;
    
    Eigen::MatrixXd T;
    Eigen::SparseMatrix<double, Eigen::RowMajor> I, J, Tinv, zetaZetaT, TzzT;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Amat, Bmat, Apre, Bpre, Ainv;
    Eigen::VectorXd zetaTinv, TinvZeta, b, bx;
    Eigen::SparseVector<double> zeta;
    
    // C++ vector equivalents of the above
    std::vector<std::vector<double>> xStates;
    std::vector<double*> xVec;
    std::vector<std::vector<double>> BpreVec, zetaZetaTVec, TinvVec, TzzTVec, AinvVec;
    std::vector<double> bxVec, zetaVec, TinvZetaVec;
    
    // Variables used to
    int zetaStartIdx, zetaEndIdx;
    
    // zeta^T * T^{-1} * zeta
    double zTz;
    
    // Sum of the difference between the refence and the optimised states (used for debugging purposes only)
    double diffsum;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Bowed1DWaveFirstOrder)
};
