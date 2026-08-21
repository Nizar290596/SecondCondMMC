
clear all

%read data of reactions
%fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D-GammaCurvefitted/soot1D/sootMechComplete.txt';
fileDict = '/Users/zhijiehuo/Dropbox/SootGenerator/Foam1D-GammaCurvefitted/soot1D/sootMechCompleteZhijie-case94.txt';
%fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/testing/sootMechTest.txt';
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
%read data and delete comments
for i=(1:size(myData{1},1))
    if contains(myData{1}{i},comSty)
        chemkData{i,:} = extractBefore(myData{1}{i},comSty);
    else
        chemkData{i,:} = myData{1}{i};
    end
end

%raw data chemk
cleanData = chemkData(~cellfun('isempty',chemkData));
%cleanData = chemkData;

for i = (1:size(cleanData,1))
    if contains(cleanData{i},'SPECIES')
        SpStartLine = i+1;
    end
    
    if contains(cleanData{i},'SECTIONS')
        SpEndLine = i-1;
        SecStartLine = i+1;
    end
    
    if contains(cleanData{i},'REACTIONS')
        SecEndLine = i-1;
        RnStartLine = i+1;
    end      
end


% for chemkin format species
%{
SpRaw = cleanData(SpStartLine:SpEndLine); 
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
load('/Users/zhijiehuo/Dropbox/SootGenerator/SpeciesInfo');
SpeciesNames = cell(size(SpData,1),1);
SpeciesCompos = cell(size(SpData,1),1);
for SpN = 1: size(SpeciesNames,1)
    if contains(SpData{SpN,1},'=')
        SpeciesNames{SpN} = extractBefore(SpData{SpN,1},'=');
        SpeciesCompos{SpN} = extractAfter(SpData{SpN,1},'=');
    else       
        SpeciesNames{SpN} = SpData{SpN,1};
        SpeciesCompos{SpN} = SpData{SpN,1};
    end   
end

SectionData = split(cleanData(SecStartLine:SecEndLine));%for sectional species
%extract section nC nH
nC = str2double(SectionData(:,2));
nH = 2;
nHaj = 1;

%HCratio = [0 1];
%nH = floor(nC*HCratio); %floor the value to avoid decimal element
%nH(:,1) = nH(:,1) + 2;% nH=2 for H/C=0 for Ait1
%nHaj = nH - 1; %for Ajt1



%read sectional reactions
SecRawReactionData = (cleanData(RnStartLine:end));
NsecRn = size(SecRawReactionData,1);
reactions = cell(NsecRn,1);
reactants = cell(NsecRn,1);
reactantAi = cell(NsecRn,1);
reactantAj = cell(NsecRn,1);
products = cell(NsecRn,1);
SecCoeffs = zeros(NsecRn,7);

%seperate reactions and their coefficienrs
for i = (1:NsecRn)
    if contains(SecRawReactionData{i},'=')
        TempBin = (split(SecRawReactionData{i}));
        reactions{i} = extractBefore(SecRawReactionData{i},TempBin{end-6}); 
        SecCoeffs(i,:) = str2num(extractAfter(SecRawReactionData{i},reactions{i}));
        reactions{i} = erase(reactions{i},' ');
        reactants{i} = split(extractBefore(reactions{i},'='),'+');
        products{i} = split(extractAfter(reactions{i},'=>'),'+');        
        clear TempBin;
    end
end

%read and seperate species and their coefficients
for i = 1: size(reactants,1)
    for j = 1:size(reactants{i},1)
        [tempRA1,tempRA2] = regexp(reactants{i}{j,1},'[a-zA-Z]\w*','match','once','split');
        if isempty(tempRA2{1})
            reactants{i}{j,2} = 1;
        else
            reactants{i}{j,2} = str2double(tempRA2{1});
        end
        reactants{i}{j,3} = strcat(tempRA1,tempRA2{2});
        reactantAi{i} = contains(reactants{i}(:,1),'Ai');
        reactantAj{i} = contains(reactants{i}(:,1),'Aj');

        clear tempRA1 tempRA2
    end
end

%read and seperate species and their coefficients
for i = 1: size(products,1)
    for j = 1:size(products{i},1)
        [tempRA1,tempRA2] = regexp(products{i}{j,1},'[a-zA-Z]\w*','match','once','split');
        if isempty(tempRA2{1})
            products{i}{j,2} = 1;
        else
            products{i}{j,2} = str2double(tempRA2{1});
        end
        products{i}{j,3} = strcat(tempRA1,tempRA2{2});
        clear tempRA1 tempRA2
    end
end



ExpandReactions = cell(NsecRn,1);
FoamArrCoeffs = cell(NsecRn,1);
GammaArrConst = cell(NsecRn,1);

for ri = 1:NsecRn
    %if reactant does not have Aj or Ai, no loop is required
    %if sum(contains(reactants{ri}(:,1),'Aj'))+sum(contains(reactants{ri}(:,1),'Ai'))==0
    if sum(reactantAi{ri})+sum(reactantAj{ri})==0
        calcNCH0_1D_CFit;
    end
    
    if sum(reactantAi{ri})+sum(reactantAj{ri})==1
        if sum(contains(products{ri}(:,1),'C'))
            calcNCH1wtOxC_1D_CFit;
        else
            calcNCH1_1D_CFit;
        end
    end
    
    if sum(reactantAi{ri})+sum(reactantAj{ri})==2
        calcNCH2_1D_CFit;
    end
     
end


%clean ExpandReactions and FoamCoeffs
for ri = 1:NsecRn
    FoamArrCoeffs{ri}(cellfun(@isempty,ExpandReactions{ri}),:) = [];
    ExpandReactions{ri}(cellfun(@isempty,ExpandReactions{ri})) = [];
end


