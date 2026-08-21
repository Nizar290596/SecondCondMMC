function diff = fit_diff(x,X,Y)
% This function is called by lsqnonlin.
% x is a vector which contains the coefficients of the
% equation.  X and Y are the option data sets that were
% passed to lsqnonlin.
A=x(1);
B=x(2);
C=x(3);

diff = A.*X.^(B).*exp(C./X) - Y;
