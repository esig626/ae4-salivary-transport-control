clear 
close all
format longg
clc

global alpha_NaK alpha_Nkcc1 GtNa GtK GCaCC GCaKC GNhe1 GAe2 GAe4
global RTF s PCaCC PCaKC xl b1 b2 b3 GBB PCO2 Ca 

run Parameters

load('Par.mat')

Ca=50e-3;
PCO2=(0.197e4);
CO20=(PCO2*(CO2l+CO2e)-kn*HCO30*H0)/(2*PCO2-kp);
xl=(b2/b1)*(2*(Na0+K0+H0)+CO20-(Nae+Ke+Cle+HCO3e))-(2*(Nal0+Kl0-Na0-K0-H0)-CO20);
s=1e-3;
alpha_Nkcc1=(0.02279)*0.28;     % amol/micro-metre^3

% Ca2+ Activated Channels
VCaCC=RTF*log(Cll0/Cl0); % mV
PCaCC=1/(1+(KCaCC/Ca)^eta1);
vCaCC=PCaCC*(Va0+VCaCC)/F;

VCaKC=RTF*log(Ke/K0); % mV
PCaKC=1/(1+(KCaKC/Ca)^eta2);
vCaKC=PCaKC*(Vb0-VCaKC)/F;

% Tight Junction
VtNa=RTF*log(Nal0/Nae);
vtNa=(Vt0-VtNa)/F;

VtK=RTF*log(Kl0/Ke);
vtK=(Vt0-VtK)/F;

% Osmolarities
qa=b1*(2*(Nal0+Kl0-Na0-K0-H0)-CO20+xl);
qb=b2*(2*(Na0+K0+H0)+CO20-(Nae+Ke+Cle+HCO3e));
qt=b3*(2*(Nal0+Kl0)+xl-(Nae+Ke+Cle+HCO3e));
qtot0=(qa+qt);

% NaK-ATPase
vNaK=(r*(Ke*s)^2*(Na0*s)^3/((Ke*s)^2+alpha1*(Na0*s)^3));

% Nkcc1
Nkcc1=alpha_Nkcc1*((a1-a2*(Na0*s)*(K0*s)*(Cl0*s)^2)/(a3+a4*(Na0*s)*(K0*s)*(Cl0*s)^2));

% Ae4
vAe4=(k1*Cle*HCO30^2*(Na0)-(k2*118.7219)*Cl0*21^2*(Nae));
% Ae2
vAe2=((Cle/(Cle+KCl))*(HCO30/(HCO30+KB))-(Cl0/(Cl0+KCl))*(21/(21+KB)));

% Nhe1
vNhe1=((Nae/(Nae+KNa))*(H0/(KH+H0))-(Na0/(Na0+KNa))*(He/(KH+He)));

% Bicarbonate Buffer
vBB=(kp*CO20-kn*HCO30*H0);

% Membrane Permeability to CO2
JCO2=(PCO2*(1.99997*CO20-CO2l-CO2e));

% Find Conductivity and Densities of the system's channels and
% transporters.

% Tight Junction Na current
GtNa=qtot0*Nal0/vtNa;
tNa=GtNa*vtNa;

% Tight Junction K current
GtK=qtot0*Kl0/vtK;
tK=GtK*vtK;

% Apical Ca2+ activated Cl channel 
GCaCC=-qtot0*Cll0/vCaCC;
CaCC=GCaCC*vCaCC;

% Bicarbonate Buffer
GBB=JCO2/vBB;
BB=GBB*vBB;

% Sodium Proton Antiporter
GNhe1=BB/vNhe1;
Nhe1=GNhe1*vNhe1;

% Sodium Potassium ATPase
alpha_NaK=((tNa+tK)-Nkcc1)/(3*vNaK);
NaK=alpha_NaK*vNaK;

% Ca2+ activated K+ channel
GCaKC=(Nkcc1+2*(tNa+tK))/(3*vCaKC);
CaKC=GCaKC*vCaKC;

% Anion exchanger 4
GAe4=(Nkcc1-3*NaK+Nhe1)/(vAe4);
Ae4=(GAe4)*vAe4;

%GAe2=-(2*Nkcc1+Ae4+CaCC)/(vAe2);
%Ae2=vAe2*GAe2;

% Anion exchanger 2
GAe2=(BB-2*Ae4)/(vAe2); 
Ae2=vAe2*GAe2;

% Densities and conductances
conductances=[alpha_Nkcc1; PCO2; GtNa; GtK; GCaCC; GBB; GNhe1; alpha_NaK; GCaKC; GAe4; GAe2];

% Fluxes of the system
Fluxes=[Ae4; Ae2; Nkcc1; NaK; JCO2; BB; -CaCC; CaKC];

% Equations of the system
dx(1)=tNa-qtot0*Nal0;
dx(2)=tK-qtot0*Kl0;
dx(3)=-CaCC-qtot0*Cll0;
dx(4)=qb-qa;
dx(5)=(Nkcc1-3*NaK+Nhe1-Ae4-dx(4)*Na0)/Hi0;
dx(6)=(Nkcc1+2*NaK-CaKC-dx(4)*K0)/Hi0;
dx(7)=(2*Nkcc1+Ae2+Ae4+CaCC-dx(4)*Cl0)/Hi0;
dx(8)=(BB-2*Ae4-Ae2-dx(4)*HCO30)/Hi0;
dx(9)=(BB-Nhe1-dx(4)*H0)/Hi0;
dx(10)=(JCO2-BB-dx(4)*CO20)/Hi0;
dx(11)=-CaCC-(tNa+tK);
dx(12)=-NaK-CaKC+(tNa+tK);
dx=dx';


save('Conductances.mat','alpha_Nkcc1' ,...
    'PCO2', 'GtNa' , 'GtK' , 'GCaCC' , 'GBB' ,...
    'GNhe1',  'alpha_NaK',  'GCaKC', 'GAe4', 'GAe2')


IC = [118.700000000000,5.60000000000000,28.7000000000000,25,120,50,10];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

[t,x]=ode15s(@(t,x) Salivary(t,x,50e-3,par),[0, 30],IC,options);

Nal =x(:,1);
Kl  =x(:,2);
Hi = x(:,3);
vol =delta*x(:,3)*10^(-15); %(picoLitres)
Na  =x(:,4);
K   =x(:,5);
Cl  =x(:,6);
HCO3=x(:,7);
%H   =x(:,8);
Cll = Nal+Kl;
% pH = log(1000./H)./(log(2)+log(5));

Xi = 2439.50353087137;
HH=-(Na+K-HCO3-Cl-Xi./Hi);

figure(1)
subplot(3,3,1)
plot(t,Na,'linewidth',3)
xlabel('t (sec)')
ylabel('[Na^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,2)
plot(t,K,'linewidth',3)
xlabel('t (sec)')
ylabel('[K^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,3)
plot(t,Cl,'linewidth',3)
xlabel('t (sec)')
ylabel('[Cl^-]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,4)
plot(t,HCO3,'linewidth',3)
xlabel('t (sec)')
ylabel('[HCO_3^-]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,5)
plot(t,Nal,'linewidth',3)
xlabel('t (sec)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
ylabel('[Na^+]_l (mM)')
box off
hold on
subplot(3,3,6)
plot(t,Kl,'linewidth',3)
xlabel('t (sec)')
ylabel('[K^+]_l (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,7)
plot(t,Cll,'linewidth',3)
xlabel('t (sec)')
ylabel('[Cl^-]_l (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,8)
plot(t,vol,'linewidth',3)
xlabel('t (sec)')
ylabel('Volume (\mu m^3)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
subplot(3,3,9)
% plot(t,H,'linewidth',3)
% hold on
plot(t,HH ,'linewidth',3)
xlabel('t (sec)')
ylabel('H^+')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off

