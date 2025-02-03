%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Non Iterative Bowed stiff String FDTD
% Ist order system - Non Iterative
%         Riccardo Russo
%       University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
realTimeDraw = false;
play = false;
dampOn = true;
smSolver = true;

osFac = 1;          %Oversampling factor
SR = 44100*osFac;   %Sample rate [Hz]
T = 1;            %Time of Simulation [s]
k = 1/SR;
timeSamples = floor(T*SR);                       %Length of sim [samples]
timeVec = (1:floor(timeSamples))*k;              %Time vector [sec]
timeVecUS = (1:floor(timeSamples/osFac))*k*osFac;%Time vector after undersampling

%%%%% Define variables of the system
radius = 0.00029;
rho = 8e3;             % String Density [kg/m] 
T0 = 40;               % Tension [N] 
E = 174e9;            % Young modulus [Pa]
L = 0.7;              %length in meters
             
Area = pi*radius^2;      % Area of string section
rA = rho*Area;
I0 = (pi*radius^4)/ 4;   % Moment of Inertia
K =sqrt(E*I0/rA);    % Stiffness parameter

c = sqrt(T0/rA); 

Fb = 10*ones(1,timeSamples);            %Bow pressure scaled by rho A [force/linear density]
a = 100;            %Bow free parameter
vb = 0.2;           %Bow velocity [m/s]

bowVel = zeros(1,timeSamples);
bowVel(:) = vb;

excitPos = 0.633*L;
outPos = 0.33*L;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Damping
omega1 = 0*2*pi;
omega2 = 1000*2*pi;
T601 = 15; T602 = 10;
param1 = (-c^2 + sqrt(c^4 + 4*K^2*omega1^2))/(2*K^2);
param2 = (-c^2 + sqrt(c^4 + 4*K^2*omega2^2))/(2*K^2);

sigma0 = (6*log(10)/(param2 - param1))*(param2/T601 - param1/T602); 
sigma1 = (6*log(10)/(param2 - param1))*(-1/T601 + 1/T602);

if ~dampOn sigma0 = 0; sigma1 = 0; end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Compute grid spacing from variables
epsilon = 0;          %Deviation from stability condition. eps \in [0,1)
%h = sqrt((c^2*k^2 + 4*sigma1*k + sqrt((c^2*k^2 + 4*sigma1*k)^2+((16*K^2*k^2))))/2)*(1 + epsilon);
h = 1.0*c*k;   %linear part is unconditionally stable, I can use a smaller grid!
N = floor(L/h);
h = L/N;
nPoints = N - 1;
lambda = c*k/h;

outPoint = floor(outPos/h);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Offline initialization of matrices
Id = speye(nPoints);
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;
Dxxxx = Dxx.'*Dxx;
I = speye(2*nPoints);


%%%%% Dirac Delta approx & J coefficient
excitSample = excitPos/h;
excitPosFloor = floor(excitSample);
beta = excitSample - excitPosFloor;
delta = zeros(nPoints,1);
delta(excitPosFloor) = (1-beta)/h;
delta(excitPosFloor+1) = beta/h;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices
x = [zeros(nPoints,1) ; zeros(nPoints,1)];

zeta = [zeros(nPoints,1) ; delta];
zetaTR = zeta.';
zetaOutProd = zeta*zetaTR;

J = [zeros(nPoints),Id;
    c^2*Dxx - K*Dxxxx, -2*(sigma0*Id - sigma1*Dxx)];

A1 = I - 0.5*k*J;
A11 = A1(1:nPoints,1:nPoints);
A12 = A1(1:nPoints,nPoints + 1:end);
A21 = A1(nPoints + 1:end,1:nPoints);
A22 = A1(nPoints + 1:end,nPoints + 1:end);

%Computation of Shur Complement of A
invA11 = A11^-1; %Notice that A11=eye so invA11=A11
shurComp = A22 - A21*(invA11*A12);
invShurComp = shurComp^-1;

B1 = sparse(I + 0.5*k*J);

%%%%% Initializing outputs
OutPlay = zeros(1,floor(timeSamples/osFac));
Out = zeros(1,timeSamples);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
tic
for i = 1:timeSamples

    zeta1 = h*zetaTR*x;
    eta = zeta1 - bowVel(i);
    d = sqrt(2*a)*exp(-a*eta^2 + 0.5);
    lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2);

    if smSolver
%     %Calculating known terms
        zeta2 = zeta*zeta1;
     
        B2 = B1*x;
        B = B2 + zeta2*0.5*k*Fb(i)*(lambda - 2*d) + k*Fb(i)*d*zeta*bowVel(i);
        
        b1 = B(1:nPoints); b2 = B(nPoints + 1:end);
        
        %Sherman morrison solver
        v = (0.5*k*h*Fb(i)*lambda)*zeta;
        v1 = v(1:nPoints); v2 = v(nPoints + 1:end);
        
        invAv2 = invShurComp*v2; invAv1 = - invA11*A12*invAv2;
        invAv = [invAv1;invAv2];
        
        y2 = A11*b1; z2 = b2 - A21*y2;
        invAb2 = invShurComp*z2; invAb1 = y2 - invA11*A12*invAb2;
        invAb = [invAb1;invAb2];
        
        vt1 = zetaTR*invAv; vt2 = zetaTR*invAb;
        
        xNext = invAb - (1/(1+vt1))*invAv*vt2;
    else
        A = I + 0.5*k*h*Fb(i)*lambda*zetaOutProd - 0.5*k*J;
        B = (I + 0.5*k*h*Fb(i)*(lambda - 2*d)*zetaOutProd + 0.5*k*J)*x + k*Fb(i)*d*zeta*bowVel(i);
        xNext = A\B;
    end

    uNext = xNext(1:nPoints);
    Out(i)= uNext(outPoint);

    x = xNext;

    if realTimeDraw
        plot(xNext(1:nPoints));
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
% plot(timeVecUS,OutPlay);
plot(timeVec,Out);
title('String Displacement [m]');
legend('String Displacement at Readout Position');
xlabel('Time [s]');


% figure(2)
% plot(timeVec,BowPower);
% title('Energy Balance [W]');
% xlabel('Time [s]');