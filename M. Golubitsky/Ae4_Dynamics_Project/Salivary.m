function dx=Salivary(t,x,ca,par,g4,g2)

if t > 100
    ca = 550e-3;
end

nal=x(1);
kl=x(2);
hi=x(3);
na=x(4);
k=x(5);
cl=x(6);
hco=x(7);

h = -(na+k-hco-cl-par.Xi./hi);
co2 = (par.pco2*par.co2le - par.gb*par.kn*hco*h)/(par.a);
cll = nal + kl;
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
vt=va-vb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ca2+ Activated Channels
jcl=par.gcl*Pcl*(va+vcl)/par.F;
jk=par.gk*Pk*(vb-vk)/par.F;

% Tight Junction
jtna=par.gtna*(vt-vtna)/par.F;
jtk=par.gtk*(vt-vtk)/par.F;


% Osmolarities
qa = par.b1*(2*(nal+kl-na-k-h)-co2+par.xl);
qb = par.b2*(2*(na+k+h)+co2-par.Ie);
qt = par.b3*(2*(nal+kl)+par.xl-par.Ie);
qT = qa + qt;


% Nkcc1
jnk1=par.aNkcc1*((par.a1-par.a2*na*k*cl^2)/(par.a3+par.a4*na*k*cl^2));

% Ae4
j4=g4*(par.k1*par.cle*hco^2*na-par.k2*cl*21^2*par.nae);

j2=g2*((par.cle/(par.cle+par.KCl))*(hco/(hco+par.KB))...
                                      -(cl/(cl+par.KCl))*(21/(21+par.KB)));
% Nhe1
j1=par.g1*((par.nae/(par.nae+par.KNa))*(h/(par.KH+h))...
                              -(na/(na+par.KNa))*(par.he/(par.KH+par.he)));

jB=par.gb*(par.kp*co2-par.kn*hco*h);


dx(1)=jtna-qT*nal;
%dx(1)=0;
dx(2)=jtk-qT*kl;
dx(3)=qb-qa;
dx(4)=(jnk1-3*jnak+j1-j4-dx(3)*na)/hi;
dx(5)=(jnk1+2*jnak-jk-dx(3)*k)/hi;
dx(6)=(2*jnk1+j2+j4+jcl-dx(3)*cl)/hi;
%dx(6)=(3*jnk1-3*jnak+j2+j1+jcl-dx(3)*cl)/hi;
dx(7)=(jB-2*j4-j2-dx(3)*hco)/hi;
dx=10000*dx';
end