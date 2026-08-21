
 phi0 = [-4.98994189453472e-20];
 k = 1.38064852e-23; %[m2 kg s^-2 K^-1]

phi = abs(phi0/1500/k);
Gamma = 1-(1+phi)*exp(-phi);

ii = 1;
for Temi=50:220
    Ti = Temi*10;
    phi = abs(phi0/Ti/k);
    GammaPoints(ii,1) = 1-(1+phi)*exp(-phi);
    
    TT(ii,1) = Ti;
    
    ii = ii + 1;
end
X=TT;
Y=GammaPoints;

% Initialize the coefficients of the function.
X0=[0 0 0]';
% Calculate the new coefficients using LSQNONLIN.
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1000);
x =lsqnonlin(@fit_diff,X0,[],[],options,X,Y);

TTT = [300:10:2200];
Y_new = x(1).*TTT.^(x(2)).*exp(x(3)./TTT);

plot(X,Y,'+r',TTT,Y_new,'b')




%plot total rate.

x1 = [9.968614e+12    0.500000    0.000000     ];

x2 = [1.707068e+19    -1.432163    -656.822727    ];

ii = 1;
for Ti = 300:10:2200
    y1(ii) = x1(1)*Ti^x1(2)*exp(0);
    
    y2(ii) = x2(1)*Ti^x2(2)*exp(x2(3)/8.314*4.184/Ti);
    
    TT(ii) = Ti;
    
    ii=ii+1;
end

plot(TT,y1); hold on
plot(TT,y2,'k');
