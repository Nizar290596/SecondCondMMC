
if reactantAi{ri}(1)
    RType{1} = 'Ai';
end
if reactantAj{ri}(1)
    RType{1} = 'Aj';
end
if reactantAi{ri}(2)
    RType{2} = 'Ai';
end
if reactantAj{ri}(2)
    RType{2} = 'Aj';
end

LPi = 1;
nCpTemp = 0;%nC(3) + nC(11);% nC(2)+nC(8);
nHpTemp = 0;%1 + nC(11)-1;% nC(8)-1+2-1
%reactants are both Ai Aj   
    % substract normal species nC nH for reactant
        for j = 1:size(products{ri},1)
            Rindx = find(strcmpi(SpName,products{ri}{j,3}),1);
            if isempty(Rindx)
                %is the product Ait1 or Ajt1?
                isAit1 = strcmpi('Ait1',products{ri}{j,3});
                isAjt1 = strcmpi('Ajt1',products{ri}{j,3});
            else
                factorNum = products{ri}{j,2};
                if max(strcmp(SpElements{Rindx,2},'C'))
                nCpTemp = nCpTemp - factorNum * str2double(SpElements{Rindx,3}...
                    {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
                end
                if max(strcmp(SpElements{Rindx,2},'H'))
                nHpTemp = nHpTemp - factorNum * str2double(SpElements{Rindx,3}...
                    {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
                end
            end
        end
        
        
    %calc lump species nC nH
    LoopClow = zeros(2,1);LoopCup = zeros(2,1);
    for ReN = 1:2
        if strcmp(reactants{ri}{ReN,3},strcat(RType{ReN},'t1'))
            LoopClow(ReN) = 1; LoopCup(ReN) = size(nC,1);
        else
            Expression = strcat(RType{ReN},...
                '(?<loopStyC>[><_])(?<startC>\d+)(?<loopStyH>[><_])(?<startH>\d+)');
            LoopParameters = regexp(reactants{ri}{ReN,3},Expression,'names');
            %loop boundary for C
            if LoopParameters.loopStyC == '>'
                LoopClow(ReN) = str2double(LoopParameters.startC)+1; LoopCup(ReN) = size(nC,1);
            elseif LoopParameters.loopStyC == '<'
                LoopClow(ReN) = 1; LoopCup(ReN) = str2double(LoopParameters.startC)-1;
            elseif LoopParameters.loopStyC == '_'
                LoopClow(ReN) = str2double(LoopParameters.startC); LoopCup(ReN) = LoopClow(ReN);
            end

            
            
        end     
    
    end
    
    
%end        
   
%LoopCup = [1;1];%testing

%calculate vertex of products and coefficients
for lpCRe1 = LoopClow(1) : LoopCup(1)
    for lpCRe2 = LoopClow(2) : lpCRe1%LoopCup(2)
            
                nCp = nCpTemp;
                nHp = nHpTemp;
                %sum reactant1 nC nH
                factorNum = reactants{ri}{1,2};
                nCp = nCp + factorNum * nC(lpCRe1);
                if strcmp(RType{1},'Ai')
                    nHp = nHp + factorNum * nH;
                elseif strcmp(RType{1},'Aj')
                    nHp = nHp + factorNum * nHaj;
                end
                %sum reactant2 nC nH
                factorNum = reactants{ri}{2,2};
                nCp = nCp + factorNum * nC(lpCRe2);
                if strcmp(RType{2},'Ai')
                    nHp = nHp + factorNum * nH;
                elseif strcmp(RType{2},'Aj')
                    nHp = nHp + factorNum * nHaj;
                end

            if nCp > nC(1) && nCp <= nC(end)
                %locate p in 2-D sectional interval for seperation
                if isAit1 %is product Ait1
                    [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                    PrdtType = 'Ai';
                    Hcompensate = nHp - nH*(pC1coeff+pC2coeff);
                    %Asp = strcat('Ai_',num2str(nCintvl-1),'_',num2str(nHintvl-1))
                elseif isAjt1 %is product Ajt1
                    [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                    PrdtType = 'Aj';
                    Hcompensate = nHp - nHaj*(pC1coeff+pC2coeff);
                end
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff];
                
                buildReaction2_1Dck;
            else
                %LPi
                %{
            elseif nCp<=nC(1)
                if isAit1
                    [pC,pH,nCintvl,nHintvl] = findVerticesHC(nCp,nHp,nC,nH);
                    PrdtType = {'Ai'};
                elseif isAjt1
                    [pC,pH,nCintvl,nHintvl] = findVerticesHC(nCp,nHp,nC,nHaj);
                    PrdtType = {'Aj'};
                end
                %[pC,pH] = findVerticesHC(nCp,nHp,nC,nH);
                %nucleation, only form particles at larger edge of C
                c = (nHp-nCp/nC(1)*pH(2))/(pH(3)-pH(2));
                b = nCp/nC(1)-c;
                ABCD = [ 0 b c 0];
                LumpC = [nCintvl-1,nCintvl,nCintvl,nCintvl-1];
                LumpH = [nHintvl,nHintvl,nHintvl+1,nHintvl+1];
                buildReaction2;%build reaction
                %}
            end

            if nCp <= nC(end) && (nCp >= nHp)
                %calculate coeffs------------------------------------------------------
                if SecCoeffs(ri,1) == -1 %for coagulation
                    if size(reactants{ri},1)==2 %coag must have 2 reactant
                        ReactantNC(1) = nC(lpCRe1);
                        ReactantNC(2) = nC(lpCRe2); 
                        ReactantNH(1) = nH; %nC(lpHRe1);
                        ReactantNH(2) = nH; %nC(lpHRe2); 
                    end

                    coeffA = 1e3*CoagBeta_1Dck(ReactantNC(1),ReactantNC(2),ReactantNH(1),ReactantNH(2));
                    %sticking efficiency
                    
                    GammaArrConst{ri}(LPi,1:4) = CoagGamma_1Dck_CFit(ReactantNC(1),ReactantNC(2),ReactantNH(1),ReactantNH(2));
                    
                    %Gamma@1500K
                    GammaArrConst{ri}(LPi,5) = CoagGamma_1Dck(ReactantNC(1),ReactantNC(2),ReactantNH(1),ReactantNH(2));

                    x = GammaArrConst{ri}(LPi,:);
                    %Gamma of fitted curve@1500K 
                    GammaArrConst{ri}(LPi,6) = x(1)*1500^(x(2))*exp(x(3)/1500);   

                    FoamArrCoeffs{ri}(LPi,1) = coeffA;
                    FoamArrCoeffs{ri}(LPi,2) = 0.5;
                    FoamArrCoeffs{ri}(LPi,3) = 0;
                    FoamArrCoeffs{ri}(LPi,4) = 1;
                    FoamArrCoeffs{ri}(LPi,5) = 1;
                    FoamArrCoeffs{ri}(LPi,6) = 1;
                    FoamArrCoeffs{ri}(LPi,7) = 2;
                    clear ReactantNC;

                else
                %convert/calculate unit of coefficients
                %convert unit of A
                    
                    FoamArrCoeffs{ri}(LPi,1) = 1 * SecCoeffs(ri,1);
                    

                    FoamArrCoeffs{ri}(LPi,2) = SecCoeffs(ri,2);
                    %E(cal/mol)
                    FoamArrCoeffs{ri}(LPi,3) = SecCoeffs(ri,3);
                    FoamArrCoeffs{ri}(LPi,4) = SecCoeffs(ri,4);
                    FoamArrCoeffs{ri}(LPi,5) = SecCoeffs(ri,5);

                    ReactantNC(1) = nC(lpCRe1);
                    ReactantNC(2) = nC(lpCRe2);

                    if SecCoeffs(ri,6) == 1
                        FoamArrCoeffs{ri}(LPi,6) = sum(ReactantNC)/2;
                    elseif SecCoeffs(ri,6) == 2
                        FoamArrCoeffs{ri}(LPi,6) = max(ReactantNC);
                    elseif SecCoeffs(ri,6) == 3
                        FoamArrCoeffs{ri}(LPi,6) = sum(ReactantNC);
                    end
                    clear ReactantNC;

                    FoamArrCoeffs{ri}(LPi,7) = SecCoeffs(ri,7);
                end
            end
            %done calculate coeffs------------------------------------------------
            
            LPi = LPi + 1;
            
        
    end
end
clear RType;

    %}
