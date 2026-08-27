% Membrane Parameters
% Area
Aa=1300*2/(28.7*2.51); %?m^2
Ab=1.51*Aa; %?m^2
% Cellular Height
Hi0=28.7;  %?m 
delta=(Aa+Ab)/2;

% Ca2+-activated-K/Cl channels
KCaCC=0.26+100; % ?M
KCaKC=0.3+100;  % ?M
eta1=1.49;  % Hill Coefficient
eta2=1.7;    % Hill Coefficient

% NaK
r=1.305e6;    % M^-3 s^-1
alpha1=0.641; % M^-1

% Thermodynamical Constants
R=8.3144621;    % J mol^-1 K^-1
T=310;          % K
F=96485.3365;   % C mol^-1  
RTF=1e3*R*T/F;

% Steady State Values for pH
pHl=6.81;
pHi=6.91;
pHe=7.41;
kn = 2.6e4;             % mM^-1 s^-1  
kp = 11; 

kn=kn*0.012;
kp=kp*0.012;

% Nkkc1
a1 = 157.55;       % /s
a2 = 2.0096*1e7;   % /M^4.s
a3 = 1.0306;       % /s
a4 = 1.3852*1e6;   % /M^4.s


% Steady State Concentrations
Cle=102.6;      %mM
Nae=140.2;      %mM
Ke=5.3;             %mM
He=1e3*10^(-pHe);   %mM
HCO3e=40;       %mM
%HCO3e=20;
CO2e=1.6;
CO2e2=kn*HCO3e*He/kp;

Nal0=118.7;      %mM
Kl0=5.6;         %mM
Cll0=Nal0+Kl0;    %mM
Hl=1e3*10^(-pHl);%mM
HCO3l=Kl0+Nal0-Cll0+Hl; %mM
CO2l=11.6;      %mM

Na0=25;         %mM
K0=120;          %mM
Cl0=50;       %mM
H0=1e3*10^(-pHi);   %mM
HCO30=10;           %mM

% Membrane potential SS
Va0=-50.24;        %mV
Vb0=-62.8;      %mV
Vt0=Va0-Vb0;    %mV

% Anion Exchangers (AE2/AE4/NHE1)
ap_Ae2=0.17;
ap_Nhe1=1;
kp_Nhe1=1e-3;
kp_Ae2=6.85;

KCl=5.6;
KB=1e-4;
KNa=15;
KH=4.5e-4;

k1=1.92e-2;
k2=1.341840733667157e-5;


% Membrane Permeabilities to water
alpha2 = 2.7e-2;     % decresing alpha2 ensures water flow remains constant with increase Cai
b1 = alpha2*0.01194/(2*(K0 + Na0 - Nal0 - Kl0 + H0));
b2 = 7*b1;
b3=(0.94*b1);
b1=(b1*74.4);

% Ca=50e-3;
% PCO2=(0.197e4);
% CO20=(PCO2*(CO2l+CO2e)-kn*HCO30*H0)/(2*PCO2-kp);
% xl=(b2/b1)*(2*(Na0+K0+H0)+CO20-(Nae+Ke+Cle+HCO3e))-(2*(Nal0+Kl0-Na0-K0-H0)-CO20)
% s=1e-3;
% alpha_Nkcc1=(0.02279)*0.28;     % amol/micro-metre^3

