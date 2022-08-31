%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Iterative Bowed String 
%    Modal system - Non Iterative
%         Riccardo Russo
%       University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
realTimeDraw = true;
play = false;

osFac = 1;          %Oversampling factor
SR = 44100*osFac;   %Sample rate [Hz]
T = 0.2;            %Time of Simulation [s]

L = 1;                % String length [m]
excitPos = 0.68*L;

% radius = 0.00029;
% rho = 8e3;             % String Density [kg/m] 
% T0 = 40;               % Tension [N] 
% A = pi*radius^2;
% rA = rho*A;
% c = sqrt(T0/rA);
c = 400;

Fb = 20;            %Bow pressure scaled by rho A [force/linear density]
a = 100;            %Bow free parameter
vb = 0.2;           %Bow velocity [m/s]

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Derived Parameters
k = 1/SR;
timeSamples = floor(T*SR);                       %Length of sim [samples]
timeVec = (1:floor(timeSamples))*k;              %Time vector [sec]
timeVecUS = (1:floor(timeSamples/osFac))*k*osFac;%Time vector after undersampling

bowVel = zeros(1,timeSamples);
bowVel(:) = vb;

h = 1.00*c*k;
N = floor(L/h);
% h = 1.00*L/N;

vectorLength = N - 1; %Using Dirichlet condition

outPoint = floor(0.5*N);

%%%%% Dirac Delta approx & J coefficient
excitSample = excitPos/h;
excitPosFloor = floor(excitSample);
beta = excitSample - excitPosFloor;
delta = zeros(vectorLength,1);
delta(excitPosFloor) = (1-beta)/h;
delta(excitPosFloor+1) = beta/h;

delta = sparse(delta);

%%%%% Offline initialization of matrices
I = speye(2*vectorLength);

%%%%% Computing Eigenvectors & Eigenfreqs
eigenFreqs = ComputeEigenFreq(c, h, N, (1:vectorLength).');
eigenFreqsMat = sparse(diag(eigenFreqs.^2));
omegaMat = [eigenFreqsMat,zeros(vectorLength);zeros(vectorLength),speye(vectorLength)];

R = ComputeEigenVector(N, (1:vectorLength).', (1:vectorLength));
Rt = R.';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices
x = [zeros(vectorLength,1) ; zeros(vectorLength,1)];

zeta = [zeros(vectorLength,1) ; Rt*delta];
zetaTR = zeta.';
zetaOutProd = zeta*zetaTR;

J = [zeros(vectorLength),speye(vectorLength);-speye(vectorLength),zeros(vectorLength)];
JomegaMat = J*omegaMat;

%%%%% Initializing outputs
OutPlay = zeros(1,floor(timeSamples/osFac));
Out = zeros(1,timeSamples);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
tic
for i = 1:timeSamples

    eta = h*zetaTR*x - bowVel(i);
    d = sqrt(2*a)*exp(-a*eta^2 + 0.5);
    lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2);

    A = I + 0.5*k*h*Fb*lambda*zetaOutProd - 0.5*k*JomegaMat;
    B = (I + 0.5*k*h*Fb*(lambda - 2*d)*zetaOutProd + 0.5*k*JomegaMat)*x + k*Fb*d*zeta*bowVel(i);
    xNext = A\B;

    uNext = R*xNext(1:vectorLength);
    Out(i)= uNext(outPoint);

    if ~mod(i,osFac)
        index = i/osFac;
        OutPlay(index) = uNext(outPoint);
    end

    x = xNext;

    if realTimeDraw
        plot(R*xNext(1:vectorLength));
        ylim([-1e-4,1e-4]);
        drawnow
        hold off
    end
end
toc

if play soundsc(OutPlay,SR); end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
figure(1)
plot(timeVec,Out);
title('String Displacement [m]');
legend('String Displacement at Readout Position');
xlabel('Time [s]');

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function omega = ComputeEigenFreq(c, h, N, index)
    omega = 2*c*sin(index*pi/2/N)/h;
end

function eigenVec = ComputeEigenVector(N, index1, index2)
    eigenVec = sqrt(2/N)*sin(index1*index2*pi/N);
end