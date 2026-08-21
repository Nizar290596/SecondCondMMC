
Aindex = find([reactantAi{ri};reactantAj{ri}]==1);
if size(reactants{ri},1)==2
    if Aindex == 1
        takeReactant = 2;
        skipReactant = 1;
        RType = 'Ai';
    end
    if Aindex == 3
        takeReactant = 2;
        skipReactant = 1;
        RType = 'Aj';
    end
    if Aindex == 2
        takeReactant = 1;
        skipReactant = 2;
        RType = 'Ai';
    end
    if Aindex == 4
        takeReactant = 1;
        skipReactant = 2;
        RType = 'Aj';
    end
elseif size(reactants{ri},1)==1
    if Aindex == 1
        skipReactant = 1;
        RType = 'Ai';
    end
    if Aindex == 2
        skipReactant = 1;
        RType = 'Aj';
    end
end

LPi = 1;
nCpTemp = 0;%nC(3) + nC(11);% nC(2)+nC(8);
nHpTemp = 0;%1 + nC(11)-1;% nC(8)-1+2-1
if size(reactants{ri},1)==2
    % add normal species nC nH for reactant
    Rindx = find(strcmpi(SpName,reactants{ri}{takeReactant,3}),1);
    factorNum = reactants{ri}{takeReactant,2};
    if max(strcmp(SpElements{Rindx,2},'C'))
    nCpTemp = nCpTemp + factorNum * str2double(SpElements{Rindx,3}...
        {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
    end
    if max(strcmp(SpElements{Rindx,2},'H'))
    nHpTemp = nHpTemp + factorNum * str2double(SpElements{Rindx,3}...
        {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
    end
end   
    
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
    %reactants{ri}{skipReactant,3}
    if strcmp(reactants{ri}{skipReactant,3},strcat(RType,'t1'))
        LoopClow = 1; LoopCup = size(nC,1);
    else
        Expression = strcat(RType,...
            '(?<loopStyC>[><_])(?<startC>\d+)(?<loopStyH>[><_])(?<startH>\d+)');
        LoopParameters = regexp(reactants{ri}{skipReactant,3},Expression,'names');
        %loop boundary for C
        if LoopParameters.loopStyC == '>'
            LoopClow = str2double(LoopParameters.startC)+1; LoopCup = size(nC,1);
        elseif LoopParameters.loopStyC == '<'
            LoopClow = 1; LoopCup = str2double(LoopParameters.startC)-1;
        elseif LoopParameters.loopStyC == '_'
            LoopClow = str2double(LoopParameters.startC); LoopCup = LoopClow;
        end
          
    end   
%end        
   
%calculate vertex of products and coefficients
for lpC = LoopClow : LoopCup
    
        nCp = nCpTemp;
        nHp = nHpTemp;
        factorNum = reactants{ri}{skipReactant,2};
        nCp = nCp + factorNum * nC(lpC);
        if strcmp(RType,'Ai')
            nHp = nHp + factorNum * nH;
        elseif strcmp(RType,'Aj')
            nHp = nHp + factorNum * nHaj;
        end
         
    if (nCp > nC(1)) && (nCp <= nC(end))
        %locate p in 2-D sectional interval for seperation
        %for nCp>nC1, (nHp <= nCp) is in Ai polygon,(nHp <= (nCp-1)) is in Aj polygon,
        if isAit1 && SecCoeffs(ri,7) ~= -1 %is product Ait1? and surface oxydation
                [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                PrdtType = 'Ai';
                Hcompensate = nHp - nH*(pC1coeff+pC2coeff);
            
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff];

                buildReaction1_1D;
           
            %Asp = strcat('Ai_',num2str(nCintvl-1),'_',num2str(nHintvl-1))
        elseif isAjt1 && SecCoeffs(ri,7) ~= -1 %is product Ajt1 and surface oxydation
                
                [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                PrdtType = 'Aj';
                Hcompensate = nHp - nHaj*(pC1coeff+pC2coeff);
            
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff];

                buildReaction1_1D;
        end
        
        
    %for fragementation, start from sec3  
        if isAit1 && SecCoeffs(ri,7) == -1 && (nCp > nC(2))%is product Ait1? and fragmentation
             
                %splite into 2 particles of half size
                nCp = nCp/2;
                [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                PrdtType = 'Ai';
                Hcompensate = nHp - nH*(pC1coeff+pC2coeff)*2;
               
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff] *2 ; %fragementation, number *2

                buildReaction1_1D;
            
            
            %Asp = strcat('Ai_',num2str(nCintvl-1),'_',num2str(nHintvl-1))
        elseif isAjt1 && SecCoeffs(ri,7) == -1 && (nCp > nC(2))%is product Ajt1 and fragmentation
                %splite into 2 particles with half size
                nCp = nCp/2;
                [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
                PrdtType = 'Aj';
                Hcompensate = nHp - nHaj*(pC1coeff+pC2coeff)*2;
                
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff] * 2;%fragementation, number *2

                buildReaction1_1D;
            
        end
        
        
        
    elseif nCp<=nC(1) 
            
        if isAit1 && SecCoeffs(ri,7) ~= -1
            
                pC1=0;pC2=1;pC1coeff=0;pC2coeff=nCp/nC(1);
                PrdtType = {'Ai'};
                Hcompensate = nHp - nH*pC2coeff;
                
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff];

                buildReaction1_1D;
            
            
        elseif isAjt1 && SecCoeffs(ri,7) ~= -1
            
                pC1=0;pC2=1;pC1coeff=0;pC2coeff=nCp/nC(1);
                PrdtType = {'Aj'};
                Hcompensate = nHp - nHaj*pC2coeff;
                
                LumpC = [pC1,pC2];
                AB = [pC1coeff,pC2coeff];

                buildReaction1_1D;

        end
        
    end
    
    
    %calculate Ahrr coeffs------------------------------------------------------     
    if nCp <= nC(end) && (nCp >= nHp)
        %convert/calculate unit of coefficients
        %convert unit of A
        if SecCoeffs(ri,1) == -1 %for coagulation
            if size(reactants{ri},1)==2
                Rindx = find(strcmpi(SpName,reactants{ri}{takeReactant,3}),1);
                ReactantNC(takeReactant) = str2double(SpElements{Rindx,3}{1,find(strcmp(SpElements{Rindx,2},'C'))+1});
                ReactantNC(skipReactant) = nC(lpC); 
                
                ReactantNH(takeReactant) = str2double(SpElements{Rindx,3}{1,find(strcmp(SpElements{Rindx,2},'H'))+1});
                ReactantNH(skipReactant) = nH; 
            end
            
            coeffA = CoagBeta(ReactantNC(1),ReactantNC(2),ReactantNH(1),ReactantNH(2));
            Gamma = CoagGamma(ReactantNC(1),ReactantNC(2),ReactantNH(1),ReactantNH(2));
            
            FoamArrCoeffs{ri}(LPi,1) = coeffA;
            FoamArrCoeffs{ri}(LPi,2) = 0.5;
            FoamArrCoeffs{ri}(LPi,3) = 0;
            FoamArrCoeffs{ri}(LPi,4) = 1;
            FoamArrCoeffs{ri}(LPi,5) = 1;
            FoamArrCoeffs{ri}(LPi,6) = Gamma;
            FoamArrCoeffs{ri}(LPi,7) = 2;
            clear ReactantNC;
      
        else
            if size(reactants{ri},1)==1
                FoamArrCoeffs{ri}(LPi,1) = 1 * SecCoeffs(ri,1);
            elseif size(reactants{ri},1) == 2
                FoamArrCoeffs{ri}(LPi,1) = 1e-3 * SecCoeffs(ri,1);
            elseif size(reactants{ri},1) == 3
                FoamArrCoeffs{ri}(LPi,1) = 1e-6 * SecCoeffs(ri,1);
            end

            FoamArrCoeffs{ri}(LPi,2) = SecCoeffs(ri,2);
            %E(cal/mol) -> T(K)
            FoamArrCoeffs{ri}(LPi,3) = SecCoeffs(ri,3)*4.184/8.3144621;
            FoamArrCoeffs{ri}(LPi,4) = SecCoeffs(ri,4);
            FoamArrCoeffs{ri}(LPi,5) = SecCoeffs(ri,5);

            if size(reactants{ri},1)==2 && contains(reactants{ri}{takeReactant,3},'C')
                Rindx = find(strcmpi(SpName,reactants{ri}{takeReactant,3}),1);
                ReactantNC(takeReactant) = str2double(SpElements{Rindx,3}{1,find(strcmp(SpElements{Rindx,2},'C'))+1});
                ReactantNC(skipReactant) = nC(lpC);
            else %size(reactants{ri},1)==1
                ReactantNC(skipReactant) = nC(lpC);
            end

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

   

clear RType;

    %}
