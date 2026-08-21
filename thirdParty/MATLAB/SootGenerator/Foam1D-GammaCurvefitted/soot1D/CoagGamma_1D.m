

function Gamma = CoagGamma(numC1,numC2,numH1,numH2)

%numC1=1e4;numC2=1e4;numH1=1e4;numH2=2;

parDensity1 = 1800; %[kg/m^3]
parDensity2 = 1800; %[kg/m^3]
Avog = 6.023e23; %[moleculemole]
k = 1.38064852e-23; %[m2 kg s^-2 K^-1]

%molecule volume scaling factor 
Vfactor = 1;

Tnorm = 1500; %[K]

parM1 = (12 * numC1 + numH1)/ Avog / 1000; %[kg] 
parM2 = (12 * numC2 + numH2) / Avog / 1000; %[kg] 

%r1 = ((3*parM1*Vfactor)/(4*parDensity1*pi))^(1/3); %[m]
%r2 = ((3*parM2*Vfactor)/(4*parDensity2*pi))^(1/3); %[m]

r1 = ((3*parM1)/(4*parDensity1*pi))^(1/3); %[m]
r2 = ((3*parM2)/(4*parDensity2*pi))^(1/3); %[m]

HCratio1 = numH1/numC1;
HCratio2 = numH2/numC2;



if(HCratio1+HCratio2)<0.1 % 0+0
    %C-C
    r_eq_molecule = 1.122*3.298e-10; %[m]
    Eeq_molecule = -71.4*1.38e-23;%-71.4*1.38e-23;%-0.2326e3/6.023e23; %[J]
    V_molecule1 = 20.58e-30;%14e-30;%147.58e-30; %[m3]
    V_molecule2 = 20.58e-30;%14e-30;%147.58e-30; %[m3]

elseif(HCratio1+HCratio2)<1.0 % 0+1
    %C-benzene
    r_eq_molecule = 3.3898e-10; %[m]
    Eeq_molecule = -4.2612e-21; %[J]
    V_molecule1 = 20.58e-30; %[m3]
    V_molecule2 = 147.58e-30; %[m3]
else %1+1
    %benzene
    r_eq_molecule = 3.77e-10; %[m]
    Eeq_molecule = -13.8e3/6.023e23; %[J]
    V_molecule1 = 147.58e-30; %[m3]
    V_molecule2 = 147.58e-30; %[m3]
end

phi = 1;
for(n = 1:500)

    D = r1+r2 + (0.01e-10)*n ;

    K = Eeq_molecule * pi^2 / (Vfactor^2 * V_molecule1 * V_molecule2 * D);
    %A = -2 * r_eq_molecule^6 * Eeq_molecule * pi^2 * (1/V_molecule^2);


    %-----------------

    one = D + r2 + r1;
    two = D + r2 - r1;
    thr = D - r2 + r1;
    fou = D - r2 - r1;

    first1 = -1/105*one^-5 + -(r1+r2)/21*one^-6 + -2*r1*r2/7*one^-7;
    first2 = +1/105*two^-5 + +(r2-r1)/21*two^-6 + -2*r1*r2/7*two^-7;
    first3 = +1/105*thr^-5 + +(r1-r2)/21*thr^-6 + -2*r1*r2/7*thr^-7;
    first4 = -1/105*fou^-5 + +(r1+r2)/21*fou^-6 + -2*r1*r2/7*fou^-7;
    first = r_eq_molecule^12/360 *(first1 + first2 + first3 + first4);

    second1 = log( (D^2-(r2+r1)^2)/(D^2-(r2-r1)^2));
    second2 = 2*r1*r2/(D^2-(r2+r1)^2);
    second3 = 2*r1*r2/(D^2-(r1-r2)^2);

    second = r_eq_molecule^6*2*D/6 * (second1 + second2 + second3);

    LJp = K * (first + second);

        if(LJp<phi)
            phi = LJp;
        end

end

phi = abs(phi/Tnorm/k);
Gamma = 1-(1+phi)*exp(-phi);
