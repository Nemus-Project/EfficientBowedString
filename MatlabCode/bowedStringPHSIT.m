%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Iterative Bowed String 
%       Ist order system
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
T = 0.3;            %Time of Simulation [s]

L = 0.7;                % String length [m]

excitPos = 0.633*L;
outPos = 0.33*L;

% radius = 0.00029;
% rho = 8e3;             % String Density [kg/m] 
% T0 = 40;               % Tension [N] 
% A = pi*radius^2;
% rA = rho*A;
% c = sqrt(T0/rA);
c = 150;

Fb = 30;            %Bow pressure scaled by rho A [force/linear density]
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

vectorLength = N - 1; %Using Dirichlet condition
outPoint = floor(outPos/h);

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
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices
x = [zeros(vectorLength,1) ; zeros(vectorLength,1)];

zeta = [zeros(vectorLength,1) ; delta];
zetaTR = zeta.';
zetaOutProd = zeta*zetaTR;

J = [zeros(vectorLength),speye(vectorLength);c^2*Dxx,zeros(vectorLength)];

tol = 0; tolThresh = 1e-9;

%%%%% Initializing outputs
OutPlay = zeros(1,floor(timeSamples/osFac));
Out = zeros(1,timeSamples);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
coeff1 = Fb*sqrt(2*a);
itersAverage = 0;
tic
for i = 1:timeSamples
    %Newton-Raphson solver
    tol = 1; iters = 1;
    r = x; %this because in this point of the code x = xNext
    while tol>tolThresh && iters < 1000
        iters = iters + 1;
        eta = h*zetaTR*r - bowVel(i);
        expCoeff = exp(-a*eta^2 + 0.5);
        f = (2/k)*r - (2/k)*x - J*r + zeta*coeff1*eta*expCoeff;
        fp = (2/k)*I - J + h*coeff1*(1 - 2*a*eta^2)*expCoeff*zetaOutProd;
        rNext = r - fp \ f;
        tol = max(abs(r - rNext));
        r = rNext;
    end

    xNext = 2*r - x;

    itersAverage = itersAverage + iters;

    uNext = xNext(1:vectorLength);
    Out(i)= uNext(outPoint);

    x = xNext;

    if realTimeDraw
        plot(xNext(1:vectorLength));
        ylim([-1e-4,1e-4]);
        drawnow
    %     pause(0.05);
        hold off
    end
end
toc

itersAverage = itersAverage/timeSamples

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Retrieving undersampled output for comparison with other schemes
finalOSFac = 1;
if finalOSFac>osFac disp("Undersampling Error."); return; end

OutPlay = zeros(1,floor(timeSamples/(osFac/finalOSFac)));

for i=1:size(Out,2)
    if ~mod(i,osFac) || mod(i,osFac) == osFac/finalOSFac
        index = i/(osFac/finalOSFac);
        OutPlay(index) = Out(i);
    end
end

if play soundsc(OutPlay,SR); end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
figure(1)
plot(timeVec,Out);
title('String Displacement [m]');
legend('String Displacement at Readout Position');
xlabel('Time [s]');
