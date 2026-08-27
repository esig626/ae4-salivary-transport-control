clear; close all; clc; format longg;




par.gb = 0.199652558943136;
par.pco2= 1970;
par.kp = 0.132;
par.kn = 312;
par.Xi = 2439.50353087137;
par.a = par.pco2*1.99997 - par.kp*par.gb;
par.co2le = 13.2;
par.aNaK = 0.00138020583534274;
par.r = 1305000;
par.ke = 5.3;
par.alpha1 = 0.641;
par.s = 0.001;
par.RTF = 26.7137302360965;
par.KCaCC = 100.26;
par.KCaKC = 100.3;
par.eta1 = 1.49;
par.eta2 = 1.7;
par.F = 96485.3365;
par.gtna = 357.863312079523;
par.aNkcc1 = 0.0063812;
par.gtk = 25.893178360928;
par.gcl = 20468119.8347071;
par.gk=87495178.2340436;
par.g1 = 0.674737308094733;
par.g2 = 0.399084872489615;
par.g4 = 18.9502607221547;
par.b1 = 0.000579346121973713;
par.b2 = 5.45083716910752e-05;
par.b3 = 7.31969562708724e-06;
par.xl = 48.8001357236744;
par.k1 = 0.0192;
par.k2 = 0.00159305881398359;
par.cle = 102.6;
par.nae = 140.2;
par.KCl = 5.6;
par.KB = 0.0001;
par.KNa = 15;
par.KH = 0.00045;
par.he = 3.8904514499428e-05;
par.a4 = 1.3852e-06;
par.a2 = 2.0096e-05;
par.a1 = 157.55;
par.a3 = 1.0306;
par.kn = 312;
par.Ie = 288.1;
par.delta = 4.52961672473868e-14;
save('Par.mat')