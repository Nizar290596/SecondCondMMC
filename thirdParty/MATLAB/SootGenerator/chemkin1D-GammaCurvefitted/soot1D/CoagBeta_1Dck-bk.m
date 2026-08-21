
%calculate coagulation coefficient
function A_b = CoagBeta(numC1,numC2,numH1,numH2)
%numC1 = 3548;
%numC2 = 2791100;
parDensity1 = 1800; %[kg/m^3]
parDensity2 = 1800; %[kg/m^3]
Avog = 6.023e23; %[moleculemole]
k = 1.38064852e-23; %[m2 kg s^-2 K^-1]

%K_coag = gamma * beta * ni * nj 
%[molecule collision/(m^3-s)] = ([1][m^3/s][molecule/m^3][molecule/m^3])
%[mole collision/(m^3-s)] = [1][m^3/s] Avog [mole/m^3]Avog [mole/m^3]/ Avog
%[kmole collision/(m^3-s)] = 1/kAvog [m^3/s-kmole] kAvog [kmole/m^3] kAvog [kmole/m^3]

%                         = kAvog * gamma * beta * Xi * Xj  

%                         = gamma * A_b * T^0.5 * Xi * Xj

%A_b = 1e3 Avog * beta / T^0.5 , this is reasonably contant wrt T, 
%calculated by T=1500
%beta = min(B_free_molecule_regimn,B_continum)

parM1 = (12 * numC1 + numH1)/ Avog / 1000; %[kg] 
parM2 = (12 * numC2 + numH2) / Avog / 1000; %[kg] 

%molecule volume scaling factor ------%%%%-----!!!!!!!!!
%Vfactor = 1;


%assuming sphere soot
%parR1 = ((3*parM1*Vfactor)/(4*parDensity1*pi))^(1/3); %[m]
%parR2 = ((3*parM2*Vfactor)/(4*parDensity2*pi))^(1/3); %[m]

parR1 = ((3*parM1)/(4*parDensity1*pi))^(1/3); %[m]
parR2 = ((3*parM2)/(4*parDensity2*pi))^(1/3); %[m]


%parR1 = 4.2798e-9/2; %[m]
%parR2 = 5.64796e-9/2;


Tnorm = 1500; %[K]

%B_free = (8 pi k T)^0.5 * (1/m1+1/m2)^0.5 * (r1+r2)^2
B_free = (8*pi*k*Tnorm)^0.5 * (1/parM1+1/parM2)^0.5 * (parR1+parR2)^2; 

%B_con = C* (2kT)/(3mu) * (1/r1+1/r2) * (r1+r2)
%mu = mu @1500 
mu = 5.4e-5; %[kg/m-s]
% C = 1 + Kn*(1.257+0.4*exp(-0.55/(Kn/2)))
%Kn = lamda/(0.5(r1+r2))
%lamda = 2.37e-10 * T [m]

lamda = 2.37e-10 * Tnorm; %[m]
Kn = lamda/(0.5*(parR1+parR2));
C = 1 + Kn*(1.257+0.4*exp(-0.55/(Kn/2)));
B_con = C* (2*k*Tnorm)/(3*mu) * (1/parR1+1/parR2) * (parR1+parR2);

beta = min(B_free,B_con)/sqrt(Tnorm);
A_b = 1e3 * Avog * beta;



