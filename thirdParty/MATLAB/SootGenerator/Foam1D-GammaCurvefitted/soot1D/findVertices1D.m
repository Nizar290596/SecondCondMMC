%nCp  - product C 
%nC  - grid of C 
function [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC)

nCintvl = 1;
if nCp <= nC(end) 
   
    %find interval in nC
    while (nCp>nC(nCintvl))
        nCintvl = nCintvl +1;  
    end        
        pC1 = nCintvl-1;
        pC2 = nCintvl; 
        pC1coeff = (nC(pC2)-nCp)/(nC(pC2)-nC(pC1));
        pC2coeff = (nCp-nC(pC1))/(nC(pC2)-nC(pC1));
else
        pC1 = 0;
        pC2 = 0; 
        pC1coeff = 0;
        pC2coeff = 0;
end

