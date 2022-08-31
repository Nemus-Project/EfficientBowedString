%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Iterative Bowed String 
% Ist order system - Non Iterative
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
T = 1;            %Time of Simulation [s]

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

Fb = 10;            %Bow pressure scaled by rho A [force/linear density]
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

    A = I + 0.5*k*h*Fb*lambda*zetaOutProd - 0.5*k*J;
    B = (I + 0.5*k*h*Fb*(lambda - 2*d)*zetaOutProd + 0.5*k*J)*x + k*Fb*d*zeta*bowVel(i);
    xNext = A\B;

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