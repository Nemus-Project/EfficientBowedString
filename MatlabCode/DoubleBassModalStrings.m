%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Double Bass Bowed Stiff String 
%    Modal system - Non Iterative
%         Riccardo Russo
%       University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;

%if set true the system is solved with Sherman Morrison formula, otherwise with backslash
smSolver = true;

%sets if to use the improved friction model from desvages
desvagesFriction = false;

stringToPlay = 3;   %0=G2, 1=D2, 2=A1, 3=E1

%sets if to let the string free to vibrate at the end or to stop it
freeVib = false;

saveAudio = false;
if saveAudio
    cond = '_Stop';
    if freeVib
        cond = '_Free';
    end
switch stringToPlay
    case 0
        string = 'G2';
    case 1
        string = 'D2';
    case 2
        string = 'A1';
    case 3
        string = 'E1';
end
    fileName = strcat('../Sounds/DoubleBass/Notes/',string,cond,'.wav');
end

osFac = 2;          %Oversampling factor
SR = 44100*osFac;   %Sample rate [Hz]
T = 5;              %Time of Simulation [s]

%%%%% Time Parameters
k = 1/SR;
timeSamples = floor(T*SR);                       %Length of sim [samples]
timeVec = (1:floor(timeSamples))*k;              %Time vector [sec]
timeVecUS = (1:floor(timeSamples/osFac))*k*osFac;%Time vector after undersampling

%%%%% String Length & Fretting Position
baseLength = 1.06;  %String base length [m]
frettingPos = 1;
                %27/32; %3rdm
                %64/81; 3rdM
                %8/9;   %1 semitone
                %4/5;   %4th
                %2/3;   %5th
                %1/2^(nSemitone/12) %semitones

L = baseLength*frettingPos;           

Fb = zeros(1,timeSamples);

%%%%% Double bass Strings
switch stringToPlay
    case 0
        % G2           
        radius = 5.18e-04;                  % string radius [m]
        rho = 7.8532e+03;                     % string Density [kg/m] 
        T0 = 285.53;                           % tension [N] 
        A = pi*radius^2;
        rA = rho*A;
        E = 200e9;                          % young modulus [Pa]
        Inertia = (pi*radius^4)/ 4;        % moment of inertia    
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.833;
        outPos = 0.33*L;
    
        if freeVib
            startFb = 20; maxFb = 20; endFb = 0;
        else
            maxFb = 35;
                    %20; %(5th)
        end
    case 1
        % D2           
        radius = 6.99e-04;
        rho = 7.8437e+03;                     
        T0 = 291.54;                         
        A = pi*radius^2;
        rA = rho*A;
        E = 200e9;                           
        Inertia = (pi*radius^4)/ 4;          
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.6;%0.833;
        outPos = 0.33*L;
        
        if freeVib
            startFb = 30; maxFb = 30; endFb = 0;
        else
            maxFb = 35;
                    %20; %(5th)
        end
    case 2
        % A1           
        radius = 9.5e-04;
        rho = 7.8158e+03;                          
        T0 = 301.35;                           
        A = pi*radius^2;
        rA = rho*A;
        E = 200e9;                          
        Inertia = (pi*radius^4)/ 4;             
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.6;%0.833;
        outPos = 0.33*L;

        if freeVib
            startFb = 10; maxFb = 10; endFb = 0;
        else
            maxFb = 20;
                    %20;%(5th)
        end
    case 3
        % E1 
        radius = 12.86e-04;
        rho = 7.8375e+03;                          
        T0 = 310.65;                           
        A = pi*radius^2;
        rA = rho*A;
        E = 200e9;                          
        Inertia = (pi*radius^4)/ 4;            
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

%         excitPos = 0.733*L; %G2 C2
%         outPos = 0.53*L;
        excitPos = 0.833;
        outPos = 0.33*L;

        if freeVib
            startFb = 10; maxFb = 10; endFb = 0;
        else
            maxFb = 35;
                    %20;%(5th)
        end
end

a = 100;            %Bow free parameter
muD = 0.3;          %Desvages friction parameter


%%%%% Bow Speed & Pressure
bowVel = zeros(1,timeSamples);

% %Linear ramp
bowRampLength = 2*timeSamples/T;
maxVb = 0.2; startVb = 0.1; endVb = 0.0;
timeFrac = 2;
bowVel(1:floor(1*bowRampLength/timeFrac)) = linspace(startVb,maxVb,floor(1*bowRampLength/timeFrac));
bowVel(floor(1*bowRampLength/timeFrac)+1:floor(2*bowRampLength/timeFrac)) = maxVb;
bowVel(floor(2*bowRampLength/timeFrac)+1:floor(3*bowRampLength/timeFrac))= linspace(maxVb,endVb,bowRampLength/timeFrac);

if freeVib
    Fb(1:floor(1*bowRampLength/5)) = linspace(startFb,maxFb,floor(1*bowRampLength/5));
    Fb(floor(1*bowRampLength/5)+1:floor(2*bowRampLength/5)) = maxFb;
    Fb(floor(2*bowRampLength/5)+1:floor(3*bowRampLength/5))= linspace(maxFb,endFb,bowRampLength/5);
else
    Fb(:) = maxFb;
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and eigenvectors
modesNumber = 1;
while true
    if ComputeEigenFreq(T0, E,Inertia,L,rA,modesNumber)>(20e3*2*pi)
        modesNumber = modesNumber-1;
        break; 
    end
    modesNumber = modesNumber+1;
end

eigenFreqs = ComputeEigenFreq(T0, E,Inertia,L,rA,(1:modesNumber).');

modesNumber = length(eigenFreqs) ;

modesIn = zeros(modesNumber,1);
modesOut = zeros(modesNumber,1);

for i = 1:modesNumber
    modesIn(i) = sqrt(2/L)*sin(i*pi*excitPos/L);
    modesOut(i) = sqrt(2/L)*sin(i*pi*outPos/L);
end

waveNumsSq = (-c^2 + sqrt(c^4 + 4*K^2*eigenFreqs.^2))/(2*K^2);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing Damping Parameters
sigma0 = zeros(modesNumber,1);
for i = 1:modesNumber
    sigma0(i) = - ComputeDamp(rho,radius,E,T0,eigenFreqs(i));
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices
I = speye(2*modesNumber);

x = [zeros(modesNumber,1) ; zeros(modesNumber,1)];

zeta = [zeros(modesNumber,1) ; modesIn];
zetaTR = zeta.';
zetaOutProd = zeta*zetaTR;

JOmega = [zeros(modesNumber),eye(modesNumber);-diag(eigenFreqs.^2),zeros(modesNumber)];

C = [zeros(modesNumber),zeros(modesNumber);zeros(modesNumber),diag(sigma0)];

%%%%% Initializing outputs
OutPlay = zeros(1,floor(timeSamples/osFac));
Out = zeros(1,timeSamples);

%%%%% Initializing Sherman Morrison Algorithm
%Offline computation of part of matrix A and extracting diagonal components
%into vectors for speeding up computation
A = full(I - 0.5*k*JOmega + 0.5*k*C);
A11 = diag(A(1:modesNumber,1:modesNumber));
A12 = diag(A(1:modesNumber,modesNumber + 1:end));
A21 = diag(A(modesNumber + 1:end,1:modesNumber));
A22 = diag(A(modesNumber + 1:end,modesNumber + 1:end));

%Computation of Shur Complement of A
invA11 = (1./A11); %Notice that A11=eye so invA11=A11
shurComp = A22 - A21.*(invA11.*A12);
invShurComp = 1./shurComp;

B1 = sparse(I + 0.5*k*JOmega - 0.5*k*C);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
tic
for i = 1:timeSamples
    
    %calculating bow input
    zeta1 = zetaTR*x;
    eta = zeta1 - bowVel(i);
    if desvagesFriction
        %Desvages friction
        d = sqrt(2*a)*exp(-a*eta^2 + 0.5) + 2*muD*atan(eta/0.02)/pi/eta;
        lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2) + 2*muD*50/pi/(2500*eta^2 + 1);
    else
        %Bilbao friction
        d = sqrt(2*a)*exp(-a*eta^2 + 0.5);
        lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2);
    end
    
    if smSolver
%     %Calculating known terms
        zeta2 = zeta*zeta1;
     
        B2 = B1*x;
        B = B2 + zeta2*0.5*k*Fb(i)*(lambda - 2*d) + k*Fb(i)*d*zeta*bowVel(i);
        
        b1 = B(1:modesNumber); b2 = B(modesNumber + 1:end);
        
        %Sherman morrison solver
        v = (0.5*k*Fb(i)*lambda)*zeta;
        v1 = v(1:modesNumber); v2 = v(modesNumber + 1:end);
        
        z1 = v2;
        invAv2 = invShurComp.*z1; invAv1 = - invA11.*A12.*invAv2;
        invAv = [invAv1;invAv2];
        
        y2 = A11.*b1; z2 = b2 - A21.*y2;
        invAb2 = invShurComp.*z2; invAb1 = y2 - invA11.*A12.*invAb2;
        invAb = [invAb1;invAb2];
        
        vt1 = zetaTR*invAv; vt2 = zetaTR*invAb;
        
        xNext = invAb - (1/(1+vt1))*invAv*vt2;
    else
        A = I + 0.5*k*Fb(i)*lambda*zetaOutProd - 0.5*k*JOmega;
        B = (I + 0.5*k*Fb(i)*(lambda - 2*d)*zetaOutProd + 0.5*k*JOmega - k*C)*x + k*Fb(i)*d*zeta*bowVel(i);

        xNext = A\B;
    end
    
    uNext = (modesOut.')*xNext(1:modesNumber);
    Out(i)= uNext;

    x = xNext;
end
realTimeFrac = toc/T

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Retrieving undersampled output

%Calculating derivative for better sound output
OutDiff = diff(Out);

finalOSFac = 1;
if finalOSFac>osFac disp("Undersampling Error."); return; end

OutPlay = zeros(1,floor(timeSamples/(osFac/finalOSFac)));

%lowpass(Out,20000,SR);
for i=1:size(OutDiff,2)
    if ~mod(i,osFac) || mod(i,osFac) == osFac/finalOSFac
        index = i/(osFac/finalOSFac);
        OutPlay(index) = OutDiff(i);
    end
end

%if play soundsc(OutPlay,SR/osFac*finalOSFac); end
if play soundsc(OutDiff,SR); end
if saveAudio
    audiowrite(fileName,OutPlay/max(abs(OutPlay)),SR/osFac*finalOSFac);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
fontSize = 18;
lineWidth = 1.5;

figure
plot(timeVec*1000,Out)
ylim([min(Out),max(Out)]);
xlabel('Time [s]');
ylabel("u(t,x_o) [m]");

hold on
plot(timeVec*1000,bowVel*1e-4);
hold on
plot(timeVec*1000,Fb*1e-6);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function omega = ComputeEigenFreq(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
end

function sigma = ComputeDamp(rho,r,E,T0,omega)
    %Desvages
    rhoAir = 1.225;
    muAir = 1.619e-5;
    d0 = -2*rhoAir*muAir/(rho*r^2);
    d1 = -2*rhoAir*sqrt(2*muAir)/(rho*r);
    d2 = -1/18000;
    d3 = -0.003*E*rho*pi^2*r^6/(4*T0^2); 

    sigma = d0 + d1*sqrt(omega) + d2*omega + d3*omega^3;
    
    %Issanchou
%     rA = rho*pi*r^2;
%     nu = omega/2/pi;
%     R = 2*pi*rhoAir + 2*pi*2*r*sqrt(pi*rhoAir*muAir*nu);
%     Qa = R/(2*pi*rA*nu);
%     I = (pi*r)/ 4;
%     Qve = (4*pi^2*E*I*nu^2*4.5e-3)/T0^2;
%     Qdisl = 2.02e-4;
% 
%     sigma = pi*nu*(Qa + Qdisl + Qve);
end