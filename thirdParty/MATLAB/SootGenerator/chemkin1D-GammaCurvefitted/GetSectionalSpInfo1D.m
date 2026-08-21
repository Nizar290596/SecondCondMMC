
clear all

%read data of reactions
fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/sootMechanism/sootMech.txt';
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
%load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D/SpeciesInfo1D');
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
HCratio = 0;
nH = 2;
nHaj = 1;

save('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D/SectionInfo1D','nC','nH','nHaj');
