/*
  ==============================================================================

    Bowed1DWaveFirstOrder.cpp
    Created: 25 Apr 2022 12:20:12pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Bowed1DWaveFirstOrder.h"

//==============================================================================
Bowed1DWaveFirstOrder::Bowed1DWaveFirstOrder (double k) : k (k)
{
    L = 0.7;
    c = 300;
    h = c * k; // using h to calculate number of modes
    N = floor (L / h);
    
    h = L / N; // recalculation of h for the non-modal scheme
    
    // calculate number of states in first order system
#ifdef MODAL
    NN = 2 * (N - 1);
#else
    NN = 2 * N - 1;
#endif
    
    a = 100;
    xB = 0.633 * L;
    vB = 0.2;
    Fb = 5;

    outPos = 0.33 * L;

    // Initialise xVectors
    xStates.resize (2);
    xVec.resize (2);
    // initialise states container with two vectors of 0s
    xStates = std::vector<std::vector<double>> (2,
                                               std::vector<double> (NN, 0));
    // initialise pointers to state vectors
    for (int i = 0; i < 2; ++i)
        xVec[i] = &xStates[i][0];
    
    
    // Initialise x for optimised algorithm
    using namespace Eigen;
    
    xNext = VectorXd (NN);
    xNext.setZero();

    x = VectorXd (NN);
    x.setZero();
      
    // Initialise x for reference algorithm
    xNextRef = VectorXd (NN);
    xNextRef.setZero();

    xRef = VectorXd (NN);
    xRef.setZero();

    I = SparseMatrix<double, RowMajor> (NN, NN);
    I.setIdentity();
    J = SparseMatrix<double, RowMajor> (NN, NN);
    J.setZero();
    
    for (int i = 0; i < N-1; ++i)
    {
        // top right quadrant
        J.coeffRef(i,N+i) += c / h;
        J.coeffRef(i+1,N+i) += -c / h;
        
        // bottom left quadrant
        J.coeffRef(N+i, i) += -c / h;
        J.coeffRef(N+i, i+1) += c / h;

    }

    using namespace Eigen;
    
    T = I / k - J / 2.0;
    
    Tinv = T.inverse().sparseView();
    
//    Tinv.pruned();
    Tinv.prune (1e-8); //test the prune value and whether it makes it faster..
    
    TzzT = SparseMatrix<double> (NN, NN);

    // Vector forms
    TinvVec = std::vector<std::vector<double>> (NN,
                                               std::vector<double> (NN, 0));
    for (int i = 0; i < NN; ++i)
        for (int j = 0; j < NN; ++j)
            TinvVec[i][j] = Tinv.coeff(i, j);
    TinvZetaVec = std::vector<double> (NN, 0);
    
    Apre = I / k - J / 2;
    Bpre = I / k + J / 2;
    
    Amat = SparseMatrix<double, RowMajor> (NN, NN);
    Bmat = SparseMatrix<double, RowMajor> (NN, NN);
    Amat.setZero();
    Bmat.setZero();

    zeta = SparseVector<double> (NN);
    zeta.setZero();
    
    bx = SparseVector<double> (NN);
    bx.setZero();

    bxVec = std::vector<double> (NN, 0);
    zetaVec = std::vector<double> (NN, 0);
    
    TzzTVec  = std::vector<std::vector<double>> (NN,
                                                std::vector<double> (NN, 0));
    AinvVec  = std::vector<std::vector<double>> (NN,
                                                std::vector<double> (NN, 0));

#ifdef MODAL
    
#else
    zeta.coeffRef(N + (int)floor(xB * N / L)) = 1.0 / h;
    zetaVec[N + (int)floor(xB * N / L)]= 1.0 / h;
    recalculateZeta(); // if done in the loop this can be excluded here
#endif
    
    BpreVec = std::vector<std::vector<double>> (NN,
                                               std::vector<double> (NN, 0));
    
    for (int i = 0; i < NN; ++i)
        for (int j = 0; j < NN; ++j)
            BpreVec[i][j] = Bpre.coeff(i, j);

}

Bowed1DWaveFirstOrder::~Bowed1DWaveFirstOrder()
{
}

void Bowed1DWaveFirstOrder::paint (juce::Graphics& g)
{
    // clear the background
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    
#ifdef RUN_ALL
    // draw the state of the reference solution
    g.setColour (Colours::green);
    g.strokePath (visualiseState (g, 50000, &xRef.coeffRef (0), 0.2 * getHeight()), PathStrokeType(2.0f));

    // draw the state of the optimised matrix form
    g.setColour (Colours::yellow);
    g.strokePath (visualiseState (g, 50000, &x.coeffRef (0), -0.2 * getHeight()), PathStrokeType(2.0f));
    
    // draw the state of the vector form
    g.setColour (Colours::cyan);
    g.strokePath (visualiseState (g, 50000, xVec[1], 0), PathStrokeType(2.0f));
#else
    // only draw the state of the vector form
    g.setColour (Colours::cyan);
    g.strokePath (visualiseState (g, 50000, xVec[1], 0), PathStrokeType(2.0f));

#endif
    
}


Path Bowed1DWaveFirstOrder::visualiseState (Graphics& g, double visualScaling, double* x, double offset)
{
    // String-boundaries are in the vertical middle of the component with a given offset
    double stringBoundaries = getHeight() / 2.0 + offset;
    
    // initialise path
    Path stringPath;
    
    // start path
    stringPath.startNewSubPath (0, 0 * visualScaling + stringBoundaries);
    
    double spacing = getWidth() / static_cast<double>(N);
    double xLoc = spacing;
    
    for (int l = 0; l < N; l++)
    {
        // Needs to be -x, because a positive x would visually go down
        float newY = -x[l+N] * visualScaling + stringBoundaries;
        
        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newY))
            newY = 0;
        
        stringPath.lineTo (xLoc, newY);
        xLoc += spacing;
    }
    
    stringPath.lineTo (xLoc, stringBoundaries);
    return stringPath;
}


void Bowed1DWaveFirstOrder::recalculateZeta()
{
    // Identify where the non-zero values of zeta are
    bool zetaFlag = false;
    for (int i = 0; i < NN; ++i)
    {
        if (zetaVec[i] != 0)
        {
            zetaFlag = true;
            zetaStartIdx = i;
        }
        if (zetaFlag && zetaVec[i] == 0)
        {
            zetaFlag = false;
            zetaEndIdx = i;
        }
    }
    
    // Get the zeta * zeta^T matrix
    zetaZetaT = (zeta * zeta.transpose()).pruned();
    
    zetaZetaTVec = std::vector<std::vector<double>> (NN, std::vector<double> (NN, 0));
    for (int i = 0; i < NN; ++i)
        for (int j = 0; j < NN; ++j)
            zetaZetaTVec[i][j] = zetaZetaT.coeff (i, j);
    
    // Calculate T^{-1}zeta
    TinvZeta = Tinv * zeta;
   
    // Get it in c++ vector form
    for (int i = 0; i < NN; ++i)
    {
        TinvZetaVec[i] = 0;
        for (int j = 0; j < NN; ++j)
            TinvZetaVec[i] += TinvVec[i][j] * zetaVec[j];
    }
    
    /// Sherman-Morrison
    
    // Calculate zeta^T * T^{-1} * zeta
    zTz = 0;
    for (int i = 0; i < NN; ++i)
        zTz += zetaVec[i] * TinvZetaVec[i];
    
    // Calculate T^{-1} * zeta * zeta^T * T^{-1}
    TzzT = (TinvZeta * zeta.transpose() * Tinv).pruned();
    
    // Get it in c++ vector form
    for (int i = 0; i < NN; ++i)
        for (int j = 0; j < NN; ++j)
            TzzTVec[i][j] = TzzT.coeff (i, j);

}

// Reference solution (using matrix inversion)
void Bowed1DWaveFirstOrder::calculateFirstOrderRef()
{
    
    double bowLoc = xB * N / L; // xB can be made user-controlled
    eta = h * 1.0 / h * xRef.data()[N + (int)floor(bowLoc)] - vB;

    lambda = sqrt(2.0*a) * (1.0 - 2.0 * a * eta * eta) * exp(-a * eta * eta + 0.5);
    d = sqrt(2.0 * a) * exp(-a * eta * eta + 0.5);

    Bmat = Bpre + (Fb * h * (0.5 * lambda - d) * zetaZetaT);
    
    /// Linear system solve ///
    using namespace Eigen;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;

    // Matrix to invert
    Amat = Apre + (Fb * h * 0.5 * lambda * zetaZetaT);
    
    // Right hand side
    b = Bmat * xRef + Fb * zeta * d * vB;

    // For the below, see SparseLU documentation on eigen.com
    Amat.makeCompressed();
    solver.analyzePattern(Amat);
    solver.factorize(Amat);
    
    // Solve the system for x^{n+1}
    xNextRef = solver.solve(b);
       
    // Update states here
    xRef = xNextRef;
}

// Sherman-Morrison using matrices
void Bowed1DWaveFirstOrder::calculateFirstOrderOpt()
{
    double bowLoc = xB * N / L; // If xB is made user-controlled, recalculateZeta() should be called either every sample, or every time xB is changed

    // Relative velocity between bow and string
    eta = h * 1.0 / h * x.data()[N + (int)floor(bowLoc)] - vB;
    
    // Non-iterative coefficients
    lambda = sqrt(2.0*a) * (1.0 - 2.0 * a * eta * eta) * exp(-a * eta * eta + 0.5);
    d = sqrt(2.0 * a) * exp(-a * eta * eta + 0.5);

    // Sherman-Morrison
    double invDiv = 1.0 + Fb * h * lambda * 0.5 * zTz;
    double divTerm = (Fb * h * lambda * 0.5) / invDiv;

    Ainv = Tinv - (TzzT * divTerm).pruned();
    
    // B matrix
    Bmat = Bpre + (Fb * h * (0.5 * lambda - d) * zetaZetaT);
    
    // Calculate x^{n+1}
    xNext = Ainv * (Bmat * x + Fb * zeta * d * vB);

    // Update states here
    x = xNext;
    
}

// Sherman-Morrison optimised (only using c++ vectors)
void Bowed1DWaveFirstOrder::calculateFirstOrderOptVec()
{
    double bowLoc = xB * N / L; // If xB is made user-controlled, recalculateZeta() should be called either every sample, or every time xB is changed
    
    // Relative velocity between bow and string
    eta = h * 1.0 / h * xVec[1][N + (int)floor(bowLoc)] - vB; // should include zeta here as well if interpolation is used
    
    // Non-iterative coefficients
    lambda = sqrt(2.0*a) * (1.0 - 2.0 * a * eta * eta) * exp(-a * eta * eta + 0.5);
    d = sqrt(2.0 * a) * exp(-a * eta * eta + 0.5);
    
    // Calculate A^{-1} (Sherman-Morrison)
    double invDiv = 1.0 + Fb * h * lambda * 0.5 * zTz;
    double divTerm = (Fb * h * lambda * 0.5) / invDiv;
    for (int i = 0; i < NN; ++i)
        for (int j = 0; j < NN; ++j)
            AinvVec[i][j] = TinvVec[i][j] - TzzTVec[i][j] * divTerm;
    
    
    
    /// Prepare the RHS of the linear system (i.e., B * x + Fb * zeta * d * vB) named bxVec here
    
    // Non-zero values due to I/k on the diagonal
    for (int i = 0; i < NN; ++i)
        bxVec[i] = BpreVec[i][i] * xVec[1][i]; // Overwrite (=) bxVec here and add (+=) in the operations below
    
    // Non-zero values due to J/2
    for (int i = 0; i < N; ++i) // top-right quadrant of BpreVec
        for (int j = std::max(N, N-1 + i); j <= std::min(NN-1, N+i); ++j)
            bxVec[i] += BpreVec[i][j] * xVec[1][j];

    for (int i = N; i < NN; ++i) // bottom-left quadrant of BpreVec
        for (int j = i - N; j <= i - (N-1); ++j)
            bxVec[i] += BpreVec[i][j] * xVec[1][j];
    
    // Add effect of bow term. If no interpolation is used, this loop is just one iteration.
    for (int i = zetaStartIdx; i < zetaEndIdx; ++i)
    {
        bxVec[i] += Fb * zetaVec[i] * d * vB; // this is assuming that vB doesn't change (otherwise we need "\mu_+ vB")
        
        // Non-iterative term
        for (int j = zetaStartIdx; j < zetaEndIdx; ++j)
            bxVec[i] += Fb * h * (0.5 * lambda - d) * zetaZetaTVec[i][j] * xVec[1][j];
    }
    
    // Calculate x^{n+1} by multiplying bxVec by A^{-1}
    for (int i = 0; i < NN; ++i) // if not modal, otherwise i < N-1
    {
        xVec[0][i] = 0;
        for (int j = 0; j < NN; ++j)
            xVec[0][i] += AinvVec[i][j] * bxVec[j];
    }
        
    // Pointer switch (update states)
    double* xTmp = xVec[1];
    xVec[1] = xVec[0];
    xVec[0] = xTmp;
    
}

double Bowed1DWaveFirstOrder::getDiffSum()
{
   diffsum = 0;
   // using xVec[1] here for the vector version because the pointer switch has been done already
   for (int i = 0; i < NN; ++i)
       diffsum += (xVec[1][i] - xNextRef.coeff (i));
    
   return diffsum;
}
