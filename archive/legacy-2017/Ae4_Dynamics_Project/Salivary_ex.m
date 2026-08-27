close all; clear; clc; format longg; load('Par.mat')
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
%% Full Model
%close all;
clc
g4 = par.g4;
g2 = par.g2*0;
IC = [118.700000000000,5.60000000000000,28.7000000000000,25,120,50,10];
[t,x]=ode15s(@(t,x) Salivary(t,x,50e-3,par,g4,g2),[0, 200],IC,options);


nal =x(:,1);
kl  =x(:,2);
hi = x(:,3);
vol =x(:,3) * par.delta; %(picoLitres)
na  =x(:,4);
k   =x(:,5);
cl  =x(:,6);
hco=x(:,7);

cll = nal+kl;

h=-(na+k-hco-cl-par.Xi./hi);
pH = log(1000./h)./(log(2)+log(5));
co2 = (par.pco2*par.co2le - par.gb*par.kn.*hco.*h)/(par.a);

figure(1)
subplot(3,3,1)
plot(t,na,'linewidth',3)
xlabel('t (sec)')
ylabel('[Na^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,2)
plot(t,k,'linewidth',3)
xlabel('t (sec)')
ylabel('[K^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,3)
plot(t,cl,'linewidth',3)
xlabel('t (sec)')
ylabel('[Cl^-]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,4)
plot(t,hco,'linewidth',3)
xlabel('t (sec)')
ylabel('[HCO_3^-]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,5)
plot(t,nal,'linewidth',3)
xlabel('t (sec)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
ylabel('[Na^+]_l (mM)')
box off
hold on
subplot(3,3,6)
plot(t,kl,'linewidth',3)
xlabel('t (sec)')
ylabel('[K^+]_l (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(3,3,7)
plot(t,cl,'linewidth',3)
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
hold on
subplot(3,3,9)
plot(t,pH,'linewidth',3)
xlabel('t (sec)')
ylabel('pH')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on

j4=g4*(par.k1*par.cle*hco.^2.*na-par.k2.*cl*21^2*par.nae);
jnk1=par.aNkcc1*((par.a1-par.a2*na.*k.*cl.^2) ...
                            ./(par.a3+par.a4.*na.*k.*cl.^2));
j2=g2*((par.cle/(par.cle+par.KCl))*(hco./(hco+par.KB))...
                          -(cl./(cl+par.KCl))*(21/(21+par.KB)));

qa = par.b1*(2*(nal+kl-na-k-h)-co2+par.xl);
qt = par.b3*(2*(nal+kl)+par.xl-par.Ie);
qT = qa + qt;

figure(2)                      
subplot(2,2,1)
plot(t,j4,'linewidth',3)
xlabel('t (sec)')
ylabel('Ae4')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(2,2,2)
plot(t,j2,'linewidth',3)
xlabel('t (sec)')
ylabel('Ae2')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(2,2,3)
plot(t,jnk1,'linewidth',3)
xlabel('t (sec)')
ylabel('Nk1')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on
subplot(2,2,4)
plot(t,qT/qT(1,1),'linewidth',3)
xlabel('t (sec)')
ylabel('FFR')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on

%%
IC2 = [28.7000000000000,25,120,50,10];
[t1,x1]=ode15s(@(t,x) Salivary2(t,x,50e-3,par),[0, 200],IC2,options);
close all;
His = x1(:,1);
vols =x1(:,1) * par.delta; %(picoLitres)
Nas  =x1(:,2);
Ks   =x1(:,3);
Cls  =x1(:,4);
HCO3s=x1(:,5);

Hs=-(Nas+Ks-HCO3s-Cls-par.Xi./His);
pHs = log(1000./Hs)./(log(2)+log(5));

figure(2)
subplot(3,3,1)
plot(t1,Nas,'linewidth',3)
xlabel('t (sec)')
ylabel('[Na^+]_i (mM)')
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
box off
hold on
subplot(3,3,2)
plot(t1,Ks,'linewidth',3)
xlabel('t (sec)')
ylabel('[K^+]_i (mM)')
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
box off
hold on
subplot(3,3,3)
plot(t1,Cls,'linewidth',3)
xlabel('t (sec)')
ylabel('[Cl^-]_i (mM)')
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
box off
hold on
subplot(3,3,4)
plot(t1,HCO3s,'linewidth',3)
xlabel('t (sec)')
ylabel('[HCO_3^-]_i (mM)')
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
box off
hold on
subplot(3,3,8)
plot(t1,vols,'linewidth',3)
xlabel('t (sec)')
ylabel('Volume (\mu m^3)')
ax=gca;
hold on
ax.LineWidth=2;
ax.FontSize=20;
box off
subplot(3,3,9)
plot(t1,pHs,'linewidth',3)
xlabel('t (sec)')
ylabel('pH')
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
box off

