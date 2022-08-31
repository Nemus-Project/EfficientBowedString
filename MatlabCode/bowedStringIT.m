%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Iterative Bowed String
%        2nd order system
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

Fb = 25;            %Bow pressure scaled by rho A [force/linear density]
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

existenceCond = (2*c/L)/(-Fb*(-11));
h = (1/existenceCond)*c*k;
if c*k/h > 1
    h = 1.0*c*k;
end

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
deltaTR = delta.';
deltaOutProd = delta*deltaTR;

%%%%% Offline initialization of matrices
I = speye(vectorLength);
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States
uNext = zeros(vectorLength,1); 
uPrev = zeros(vectorLength,1);

%The correct initialization  considers uPrev as u^0, u as u^1, therefore 
%u needs to be initialized taking into account the bow velocity.
u = -0.5*k^2*delta*Fb*sqrt(2*a)*(-bowVel(1))*exp(-a*(-bowVel(1))^2+ 0.5); 

tol = 0; tolThresh = 1e-9;

%%%%% Initializing outputs
Out = zeros(1,timeSamples);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
coeff1 = Fb*sqrt(2*a);
itersAverage = 0;
tic
for i = 1:timeSamples
    
    %Newton-Raphson solver
    b = -(2*I + c^2*k^2*Dxx)*u + uPrev*2;
    tol = 1; iters = 1;
    r = u - uPrev; %this because in this point of the code u = uNext
    while tol>tolThresh && iters < 1000
        iters   = iters + 1;
        eta = (0.5*h/k)*deltaTR*r - bowVel(i);
        expCoeff = exp(-a*eta^2 + 0.5);
        f = r + b + delta*k^2*coeff1*eta*expCoeff;
        fp = I + 0.5*h*k*coeff1*(1 - 2*a*eta^2)*expCoeff*deltaOutProd;
        rNext = r - fp \ f;
        tol = max(abs(r - rNext));
        r = rNext;
    end

    itersAverage = itersAverage + iters;

    uNext = r + uPrev;
    
    Out(i) = uNext(outPoint);

    uPrev = u;
    u = uNext;

    if realTimeDraw
        plot(uNext);
        ylabel("u(t,x_o) (m)");
        ylim([-1e-4,1e-4]);
        xlim([0,N]);
        drawnow
    %     pause(0.05);
        hold off
    end
end
toc

itersAverage = itersAverage/timeSamples
    
if play soundsc(OutPlay,SR); end

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

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
figure(1)
plot(timeVec,Out);
title('String Displacement [m]');
legend('String Displacement at Readout Position');
xlabel('Time [s]');