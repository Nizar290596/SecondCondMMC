clear all
%chemk is the raw data of rations and coefficients
%SpNames is a list of species names
%rectns is reactions 
%ArrhenCoef is the Ahrrenius coefficients
%reactionType is the type of each reacions and the related coefficients

elementsList = cell(5,1);
elementsList = {'O','H','C','N','Ar'};

%read data of reactions
fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/gasChem.txt';
%fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/NAPSchemkin/sootChemNAPS.inp';
%data = readtable(fileDict,'TextType','string');
fid = fopen(fileDict);
%data = textscan(FID,'%s');
iRow = 1;
while (~feof(fid)) 
    myData(iRow,:) = textscan(fid,'%s','delimiter', '\n');
    iRow = iRow + 1;
end
fclose(fid);
myData{1} = strtrim(myData{1});

comSty = '!';
chemkData= cell(size(myData{1},1),1);
%read data and delte comments
for i=(1:size(myData{1},1))
    if contains(myData{1}{i},comSty)
        chemkData{i,:} = extractBefore(myData{1}{i},comSty);
    else
        chemkData{i,:} = myData{1}{i};
    end
end

%raw data chemk
chemkData = chemkData(~cellfun('isempty',chemkData));
cleanData = chemkData;

for i = (1:size(cleanData,1))
    if contains(cleanData{i},'SPECIES')
        SpStartLine = i+1;
    end
    
    if contains(cleanData{i},'REACTIONS')
        SpEndLine = i-1;
        RnStartLine = i+1;
    end      
end

%{
SpRaw = cleanData(SpStartLine:SpEndLine); % for chemkin format species
for i = 1:size(SpRaw,1)
    Sptemp{i} = split(SpRaw{i});
end
SpData = {Sptemp{1}{:}};
for i = 2:size(Sptemp,2)
    SpData = { SpData{:} Sptemp{i}{:}};
end
SpData=SpData';                          % for chemikn format species
%}


SpData = split(cleanData(SpStartLine:SpEndLine));%for john format species
%extract species names
SpeciesNames = SpData(:,1);
for SpN = 1: size(SpeciesNames,1)
    if contains(SpeciesNames{SpN},'=')
        SpeciesNames{SpN} = extractBefore(SpeciesNames{SpN},'=');
    end  
end                      %for john format species


%reaction data is stored in chemk
chemk = cleanData(RnStartLine:end);

Nchemk = size(chemk,1);
rectns = cell(Nchemk,1);
ArrhenCoef = zeros(Nchemk,3);

%seperate reactions and their coefficienrs
for i = (1:Nchemk)
    if contains(chemk{i},'=')
    %rectns{i} = extractBefore(chemk{i},"   "); 
    %ArrhenCoef(i,:) = str2num(extractAfter(chemk{i},"   "));
    splitTemp = split(chemk{i});
    splitTemp = splitTemp(~cellfun('isempty',splitTemp));
    splitTempSize = size(splitTemp,1);
    rectns{i} = extractBefore(chemk{i},splitTemp(splitTempSize-2));
    %rectns{i} = extractBefore(chemk{i},"   "); 
    ArrhenCoef(i,:) = str2num(extractAfter(chemk{i},rectns{i}));
    clear splitTemp
    end
end

%clean the ' ' in reactions
for i = (1:Nchemk)
    if ~isempty(rectns{i})
        rectns{i} = erase(rectns{i},' ');
    end
end

%determine reaction type
reactionType = cell(Nchemk,2);
for i = (1:Nchemk)
    if ~isempty(rectns{i})
        %reversibility
        if contains(rectns{i},'<=>')
            reactionType{i,1} = 'reversible';
        elseif contains(rectns{i},'<=') || contains(rectns{i},'=>')
            reactionType{i,1} = 'irreversible';
        else
            reactionType{i,1} = 'reversible';
        end
        
        %third body ractions
        if contains(rectns{i},'(+M)')
            if contains(chemk{i+2},'TROE')
                reactionType{i,2} = 'ArrheniusTroeFallOffReaction';
            else
                reactionType{i,2} = 'ArrheniusLindemannFallOffReaction';
            end
        elseif contains(rectns{i},'+M')
            reactionType{i,2} = 'thirdBodyArrheniusReaction';
        else
            reactionType{i,2} = 'ArrheniusReaction';
        end      
    end        
end

%read M coefficients
for i = (1:Nchemk)
    if ~isempty(reactionType{i,2})
        if strcmp(reactionType{i,2},'thirdBodyArrheniusReaction')
            if isempty(reactionType{i+1,2})
                reactionType{i+1,3} = split(chemk(i+1),'/');
                reactionType{i+1,3} = erase(reactionType{i+1,3},' ');
            else
                reactionType{i,3} = 'all ones';
            end
        end
        
        if strcmp(reactionType{i,2},'ArrheniusLindemannFallOffReaction')
            reactionType{i+1,3} = split(extractBetween(chemk(i+1),'/','/'));
            reactionType{i+1,3} = erase(reactionType{i+1,3},' ');
            reactionType{i+1,3} = reactionType{i+1,3}(~cellfun('isempty',reactionType{i+1,3}));
            ArrhenCoef(i+1,1)  = str2num(reactionType{i+1,3}{1});
            ArrhenCoef(i+1,2)  = str2num(reactionType{i+1,3}{2});
            ArrhenCoef(i+1,3)  = str2num(reactionType{i+1,3}{3});
            
            reactionType{i+2,3} = split(chemk(i+2),'/');
            reactionType{i+2,3} = erase(reactionType{i+2,3},' ');
            
        elseif strcmp(reactionType{i,2},'ArrheniusTroeFallOffReaction')
            reactionType{i+1,3} = split(extractBetween(chemk(i+1),'/','/'));
            reactionType{i+1,3} = erase(reactionType{i+1,3},' ');
            reactionType{i+1,3} = reactionType{i+1,3}(~cellfun('isempty',reactionType{i+1,3}));
            ArrhenCoef(i+1,1)  = str2num(reactionType{i+1,3}{1});
            ArrhenCoef(i+1,2)  = str2num(reactionType{i+1,3}{2});
            ArrhenCoef(i+1,3)  = str2num(reactionType{i+1,3}{3});
            
            reactionType{i+2,3} = split(extractBetween(chemk(i+2),'/','/'));
            reactionType{i+2,3} = erase(reactionType{i+2,3},' ');
            reactionType{i+2,3} = reactionType{i+2,3}(~cellfun('isempty',reactionType{i+2,3}));
            ArrhenCoef(i+2,1)  = str2num(reactionType{i+2,3}{1});
            ArrhenCoef(i+2,2)  = str2num(reactionType{i+2,3}{2});
            ArrhenCoef(i+2,3)  = str2num(reactionType{i+2,3}{3});
            if size(reactionType{i+2,3},1) == 4
                ArrhenCoef(i+2,4) = str2num(reactionType{i+2,3}{4});
            else
                ArrhenCoef(i+2,4) = 1e15;
            end
            
            reactionType{i+3,3} = split(chemk(i+3),'/');
            reactionType{i+3,3} = erase(reactionType{i+3,3},' ');
        end
    end
end

%convert Ea to Ta (K) = Ea/R ;Ea (cal/mole)
% A's unit: (cm^3/mole)^(r-1) * s^(-1); r-order of reaction

%distinguish among uni, bi and ter reaction
%erase '<' and '>' then extract left side of '=' to reactionLeft
%extract stoichiometric coefficients
reactionsCleanedTemp = rectns;%cell(Nchemk,1);
reactionsLeft = cell(Nchemk,2);
reactionsLeftSp = cell(Nchemk,1);
reactionsLeftSpCoef = cell(Nchemk,1);
for i = (1:Nchemk)
    if ~isempty(rectns{i})
        reactionsCleanedTemp{i} = erase(reactionsCleanedTemp{i},' ');
        reactionsCleanedTemp{i} = erase(reactionsCleanedTemp{i},'>');
        reactionsCleanedTemp{i} = erase(reactionsCleanedTemp{i},'<');
        reactionsCleanedTemp{i} = replace(reactionsCleanedTemp{i},'(+M)','+M');
        reactionsLeft{i,1} = extractBefore(reactionsCleanedTemp{i},'=');
        reactionsLeftSp{i} = split(reactionsLeft{i},'+');
        for si = (1: size(reactionsLeftSp{i},1))
            reactionsLeftSpCoef{i}{si} = sscanf(reactionsLeftSp{i}{si},'%f');
            if isempty(reactionsLeftSpCoef{i}{si})
                reactionsLeftSpCoef{i}{si} = 1;
            end
        end
            
    end
end

%distinguish uni, bi ortem reaction by 1 2 3
for i = (1:Nchemk)
    if ~isempty(reactionsLeft{i,1})
        reactionsLeft{i,2} = sum([reactionsLeftSpCoef{i}{:}]); 
    end    
end

%convert to foam-format unit
FoamArrhenCoef = ArrhenCoef;
for i = (1:Nchemk)
    if ~isempty(rectns{i,1})
        if strcmp(reactionType{i,2},'ArrheniusReaction')||...
                strcmp(reactionType{i,2},'thirdBodyArrheniusReaction')
            FoamArrhenCoef(i,1) = FoamArrhenCoef(i,1)*(1e-3)^(reactionsLeft{i,2}-1);
            FoamArrhenCoef(i,3) = FoamArrhenCoef(i,3)*4.184/8.3144621;
        elseif strcmp(reactionType{i,2},'ArrheniusLindemannFallOffReaction')||...
                strcmp(reactionType{i,2},'ArrheniusTroeFallOffReaction')
            FoamArrhenCoef(i,1) = FoamArrhenCoef(i,1)*(1e-3)^(reactionsLeft{i,2}-2);
            FoamArrhenCoef(i+1,1) = FoamArrhenCoef(i+1,1)*(1e-3)^(reactionsLeft{i,2}-1);
            FoamArrhenCoef(i,3) = FoamArrhenCoef(i,3)*4.184/8.3144621;
            FoamArrhenCoef(i+1,3) = FoamArrhenCoef(i+1,3)*4.184/8.3144621;
            
        end
    end    
end



%format reaction in FOAM format
FOAMreactions = rectns;%cell(Nchemk,1);
for i = (1:Nchemk)
    if ~isempty(rectns{i})
        FOAMreactions{i} = erase(FOAMreactions{i},' ');
        FOAMreactions{i} = erase(FOAMreactions{i},'>');
        FOAMreactions{i} = erase(FOAMreactions{i},'<');
        FOAMreactions{i} = erase(FOAMreactions{i},'(+M)');
        FOAMreactions{i} = erase(FOAMreactions{i},'+M');
        FOAMreactions{i} = replace(FOAMreactions{i},'=',' = ');
        FOAMreactions{i} = replace(FOAMreactions{i},'+',' + ');
    end
end

nC = 20;
nH = 1;
secSpind = size(SpeciesNames,1);
for nCi = 1:nC
    for nHi = 1:nH
        secSpind = secSpind +1;
        SpeciesNames{secSpind} = ...
            strcat('Ai_',num2str(nCi,'%02.f'),'_',num2str(nHi,'%02.f'),'_01_1');
    end
end

for nCi = 1:nC
    for nHi = 1:nH
        secSpind = secSpind +1;
        SpeciesNames{secSpind} = ...
            strcat('Aj_',num2str(nCi,'%02.f'),'_',num2str(nHi,'%02.f'),'_01_1');
    end
end



