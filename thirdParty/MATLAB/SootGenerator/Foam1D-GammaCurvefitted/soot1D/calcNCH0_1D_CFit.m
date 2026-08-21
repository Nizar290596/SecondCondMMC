
LPi = 1;
nCp = 0;%nC(3) + nC(11);% nC(2)+nC(8);
nHp = 0;%1 + nC(11)-1;% nC(8)-1+2-1
        for j = 1:size(reactants{ri},1)
            Rindx = find(strcmpi(SpName,reactants{ri}{j,3}),1);
            factorNum = reactants{ri}{j,2};
            nCp = nCp + ...
                factorNum * str2double(SpElements{Rindx,3}...
                {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
            nHp = nHp + ...
                factorNum * str2double(SpElements{Rindx,3}...
                {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
        end
        
        for j = 1:size(products{ri},1)
            Rindx = find(strcmpi(SpName,products{ri}{j,3}),1);
            if isempty(Rindx)
                %is the product Ait1 or Ajt1?
                isAit1 = strcmpi('Ait1',products{ri}{j,3});
                isAjt1 = strcmpi('Ajt1',products{ri}{j,3});
            else
                
                
                factorNum = products{ri}{j,2};
                if max(strcmp(SpElements{Rindx,2},'C'))
                nCp = nCp - factorNum * str2double(SpElements{Rindx,3}...
                    {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
                end
                if max(strcmp(SpElements{Rindx,2},'H'))
                nHp = nHp - factorNum * str2double(SpElements{Rindx,3}...
                    {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
                end
                
            end
        end

    if nCp > nC(1)
        %locate p in 2-D sectional interval for seperation
        if isAit1 %is product Ait1
            [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
            PrdtType = 'Ai';
            Hcompensate = nHp - nH*(pC1coeff+pC2coeff);
            
        elseif isAjt1 %is product Ajt1
            [pC1,pC2,pC1coeff,pC2coeff] = findVertices1D(nCp,nC);
            PrdtType = 'Aj';
            Hcompensate = nHp - nHaj*(pC1coeff+pC2coeff);
        end
        
        LumpC = [pC1,pC2];
        AB = [pC1coeff,pC2coeff];

        buildReaction0_1D;
        %ExpandReactions{ri} = buildReaction(reactions{ri},ABCD,LumpC,LumpH);
    else
        if isAit1
            pC1=0;pC2=1;pC1coeff=0;pC2coeff=nCp/nC(1);
            PrdtType = {'Ai'};
            Hcompensate = nHp - nH*pC2coeff;
        elseif isAjt1
            pC1=0;pC2=1;pC1coeff=0;pC2coeff=nCp/nC(1);
            PrdtType = {'Aj'};
            Hcompensate = nHp - nHaj*pC2coeff;
        end
       
        LumpC = [pC1,pC2];
        AB = [pC1coeff,pC2coeff];
        buildReaction0_1D;%build reaction        
    end
    
    
    if SecCoeffs(ri,1) == -1 %for coagulation
        
        Rindx = find(strcmpi(SpName,reactants{ri}{1,3}),1);
        reactantCn1 = str2double(SpElements{Rindx,3}...
            {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
        reactantHn1 = str2double(SpElements{Rindx,3}...
            {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
        
        Rindx = find(strcmpi(SpName,reactants{ri}{2,3}),1);
        reactantCn2 = str2double(SpElements{Rindx,3}...
            {1,find(strcmp(SpElements{Rindx,2},'C'))+1});
        reactantHn2 = str2double(SpElements{Rindx,3}...
            {1,find(strcmp(SpElements{Rindx,2},'H'))+1});
        
        coeffA = CoagBeta_1D(reactantCn1,reactantCn2,reactantHn1,reactantHn2);
        
        GammaArrConst{ri}(1:4) = CoagGamma_1D_CFit(reactantCn1,reactantCn2,reactantHn1,reactantHn2);
        
        GammaArrConst{ri}(5) = CoagGamma_1D(reactantCn1,reactantCn2,reactantHn1,reactantHn2);

        x = GammaArrConst{ri};
        
        GammaArrConst{ri}(6) = x(1)*1500^(x(2))*exp(x(3)/1500);     
        
        FoamArrCoeffs{ri}(1) = coeffA;
        FoamArrCoeffs{ri}(2) = 0.5;
        FoamArrCoeffs{ri}(3) = 0;
        FoamArrCoeffs{ri}(4) = 1;
        FoamArrCoeffs{ri}(5) = 1;
        FoamArrCoeffs{ri}(6) = 1;
        FoamArrCoeffs{ri}(7) = 2;
        clear reactantCn1;
        clear reactantCn2;

    else %for soot reactions
        %convert/calculate unit of coefficients
        %convert unit of A
        if size(reactants{ri},1)==1
            FoamArrCoeffs{ri}(1) = 1 * SecCoeffs(ri,1);
        elseif size(reactants{ri},1) == 2
            FoamArrCoeffs{ri}(1) = 1e-3 * SecCoeffs(ri,1);
        elseif size(reactants{ri},1) == 3
            FoamArrCoeffs{ri}(1) = 1e-6 * SecCoeffs(ri,1);
        end

        FoamArrCoeffs{ri}(2) = SecCoeffs(ri,2);
        %E(cal/mol) -> T(K)
        FoamArrCoeffs{ri}(3) = SecCoeffs(ri,3)*4.184/8.3144621;
        FoamArrCoeffs{ri}(4) = SecCoeffs(ri,4);
        FoamArrCoeffs{ri}(5) = SecCoeffs(ri,5);

    for j = 1:2
        Rindx = find(strcmpi(SpName,reactants{ri}{j,3}),1);
        ReactantNC(j) = str2double(SpElements{Rindx,3}{1,find(strcmp(SpElements{Rindx,2},'C'))+1});
    end
        if SecCoeffs(ri,6) == 1
            FoamArrCoeffs{ri}(6) = sum(ReactantNC)/2;
        elseif SecCoeffs(ri,6) == 2
            FoamArrCoeffs{ri}(6) = max(ReactantNC);
        elseif SecCoeffs(ri,6) == 3
            FoamArrCoeffs{ri}(6) = sum(ReactantNC)/2;
        end
        clear ReactantNC;

        FoamArrCoeffs{ri}(7) = SecCoeffs(ri,7);
    end

    
    clear RType;
