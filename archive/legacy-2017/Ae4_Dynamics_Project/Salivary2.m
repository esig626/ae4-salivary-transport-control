function dx=Salivary2(t,x,ca,par,g4)

if t > 100
    ca = 550e-3;
end


hi=x(1);
na=x(2);
k=x(3);
cl=x(4);
hco=x(5);



nal = 118.7;
kl = 5.6;
cll = nal + kl;

h = -(na+k-hco-cl-par.Xi./hi);
co2 = (par.pco2*par.co2le - par.gb*par.kn*hco*h)/(par.a);

% NaK-ATPase
jnak=par.aNaK*(par.r*par.ke^2*na^3*par.s^5....
                    /(par.s^2*(par.ke^2+par.alpha1*na^3*par.s)));
vcl=par.RTF*log(cll/cl);
Pcl=1/(1+(par.KCaCC/ca)^par.eta1);
vk=par.RTF*log(par.ke/k); 
Pk=1/(1+(par.KCaKC/ca)^par.eta2);
vtna=par.RTF*log(nal/par.nae);
vtk=par.RTF*log(kl/par.ke);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
va=-((par.gtna+par.gtk)*(par.F * jnak)-(par.gk * Pk)...
    *(par.gtna*vtna+par.gtk*vtk+(par.gtna+par.gtk)*vk)...
    +(par.gcl * Pcl*vcl)*((par.gk * Pk)+(par.gtna+par.gtk)))...
    /((par.gtna+par.gtk)*((par.gk * Pk)+(par.gcl * Pcl))...
    +(par.gk * Pk)*(par.gcl * Pcl));

vb=-(((par.gtna+par.gtk)+(par.gcl * Pcl))*((par.F * jnak)...
    -(par.gk * Pk)*vk)...
    +(par.gcl * Pcl)*(par.gtna*vtna+par.gtk*vtk)...
    +(par.gtna+par.gtk)*(par.gcl * Pcl * vcl))...
    /((par.gtna+par.gtk)*(par.gk * Pk)...
    +(par.gcl * Pcl)*((par.gk * Pk)+(par.gtna+par.gtk)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ca2+ Activated Channels
jcl=par.gcl*Pcl*(va+vcl)/par.F;
jk=par.gk*Pk*(vb-vk)/par.F;

% Osmolarities
qa = par.b1*(2*(nal+kl-na-k-h)-co2+par.xl);
qb = par.b2*(2*(na+k+h)+co2-par.Ie);

% Nkcc1
jnk1=par.aNkcc1*((par.a1-par.a2*na*k*cl^2)/(par.a3+par.a4*na*k*cl^2));

% Ae4
j4=par.g4*(par.k1*par.cle*hco^2*na-par.k2*cl*21^2*par.nae);

j2=par.g2*((par.cle/(par.cle+par.KCl))*(hco/(hco+par.KB))...
                                      -(cl/(cl+par.KCl))*(21/(21+par.KB)));
% Nhe1
j1=par.g1*((par.nae/(par.nae+par.KNa))*(h/(par.KH+h))...
                              -(na/(na+par.KNa))*(par.he/(par.KH+par.he)));

jB=par.gb*(par.kp*co2-par.kn*hco*h);


dx(1)=qb-qa;
dx(2)=(jnk1-3*jnak+j1-j4-(qb-qa)*na)/hi;
dx(3)=(jnk1+2*jnak-jk-(qb-qa)*k)/hi;
dx(4)=(2*jnk1+j2+j4+jcl-(qb-qa)*cl)/hi;
%dx(6)=(3*jnk1-3*jnak+j2+j1+jcl-dx(3)*cl)/hi;
dx(5)=(jB-2*j4-j2-(qb-qa)*hco)/hi;
%dx(5)=0;
dx=1000*dx';
end