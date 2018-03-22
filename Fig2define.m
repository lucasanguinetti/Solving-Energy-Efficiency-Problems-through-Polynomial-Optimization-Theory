%This Matlab script is called by the script Fig2main.m to generate Figure 2 in the article:
%
%Andrea Pizzo, Alessio Zappone and Luca Sanguinetti, "Solving Energy Efficiency Problems
%through Polynomial Optimization Theory," IEEE Signal Processing Letters, Submitted to.
%
%This is version 1.0 (Last edited: 2018-22-03)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%

% *********************************************
%% Simulation parameters
% *********************************************

%Maximal number of antennas and users considered when plotting the EE.
Mmax = 600;
Kmax = 40;

%Select the values of log2(1+gamma) that should be considered
rateValue = 2;
%Extract the current gamma value
gammaval = 2^rateValue - 1;

%Spectral resources
Bw = 2e7;   % system bandwidth
T = 1/(Bw); %Symbol time (based on 20 MHz)

%Propagation parameters
alphaval = 3.76; %Pathloss exponent
tau = 400; %Length of coherence block (in symbols)
omegaval_db = 130;
omegaval = 10^(-omegaval_db/10);

%Hardware characterization
eta = 0.39; %Power amplifier efficiency

%SNR pilot sequence
SNRp_db = 5;
SNRp = 10^(SNRp_db/10);

%SNR payload sequence
SNR_db = 0;
SNR = 10^(SNR_db/10);
sigma2 = 1e-20/T; % noise variance (in watt)
rhoval = SNR/sigma2;

%Define the coverage area (as a square with wrap-around)
squareLength = 1000; % in m

%BS density
lambdaBS = 5; % in BSs/km^2

%Average number of base stations per setup
averageNumberofBSs = lambdaBS*(squareLength/1e3)^2;

% Energy consumption parameters (according to [8])
%
%[8] E. Björnson, L. Sanguinetti, and M. Kountouris, “Deploying dense
%networks for maximal energy efficiency: Small cells meet massive
%MIMO,” IEEE Journal on Selected Areas in Communications, vol. 34,
%no. 4, pp. 832–847, April 2016.
%
P_FIX = 10;
P_LO = .2;
P_BS = .4;
P_UE = .2;
P_COD = .1e-9;
P_DEC = .8e-9;
P_BT = .25e-9;
L_BS = 75e9;
L_UE = 3e9;
% auxiliary energy consumption parameters
aux = 3/T/tau/L_BS;
ucal = rhoval/eta*omegaval*gamma(alphaval/2+1)/(pi*lambdaBS)^(alphaval/2);
C0cal = (P_FIX + P_LO);
C1cal = (P_UE + aux*2);
C1barcal = C1cal + ucal/tau;
C3cal = 0;
D0cal = (P_BS);
D1cal = (aux*(tau + 2));
D2cal = 0;
Acal = ((P_COD + P_DEC + P_BT)/lambdaBS);

% *********************************************
%% Optimization parameters
% *********************************************

% B1bar and B2bar coefficients
b11 = (4/(alphaval-2)^2 + 1/(alphaval-1) + 2/(alphaval-2))*gammaval/tau;
b12 = 1/(alphaval-1)*gammaval/tau;
b10 = 2/(alphaval-2)/SNR*gammaval/tau;
b21 = (1 + 2/(alphaval-2))*(1 + 1/SNR)*gammaval;
b20 = (1 + 1/SNR)/SNR*gammaval;
rvalue = log2(1+gammaval)*Bw;
psivalue = rvalue*Acal + ucal;

%ASE: f(K,M)
f10 = -b20*rvalue;
f11 = rvalue;
f20 = -(b21+b10)*rvalue;
f21 = -b12*rvalue;
f30 = -b11*rvalue;

%APC: g(K,M)
g00 = -b20*C0cal;
g10 = -(b21*C0cal + b20*C1barcal + b20*psivalue);
g01 = C0cal - b20*D0cal;
g20 = -(b21*C1barcal + b21*psivalue + b10*psivalue);
g11 = C1barcal - b21*D0cal - b20*D1cal + psivalue;
g02 = D0cal;
g21 = -(b21*D1cal + b12*psivalue);
g12 = D1cal;
g30 = -b11*psivalue;

%Constraints: h(K,M) and q(K,M)
h00 = tau*b10 + b20;
h10 = tau*b11 + b21;
h01 = tau*b12 - 1;
%
q00 = -b20;
q10 = -b21 - b10;
q01 = 1;
q11 = -b12;
q20 = -b11;

%scaling up parameters for avoiding numerical problems
sK = 40;
sM = 600;
%
f10bar = f10*sK;
f11bar = f11*(sM*sK);
f20bar = f20*(sK^2);
f21bar = f21*(sM*sK^2);
f30bar = f30*(sK^3);
%
g10bar = g10*sK;
g01bar = g01*sM;
g20bar = g20*(sK^2);
g11bar = g11*(sM*sK);
g02bar = g02*(sM^2);
g21bar = g21*(sM*sK^2);
g12bar = g12*(sM^2*sK);
g30bar = g30*(sK^3);
%
h10bar = h10*sK;
h01bar = h01*sM;
q10bar = q10*sK;
q01bar = q01*sM;
q11bar = q11*(sM*sK);
q20bar = q20*(sK^2);
