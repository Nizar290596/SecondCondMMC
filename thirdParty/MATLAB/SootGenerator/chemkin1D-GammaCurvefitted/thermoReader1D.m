
GetSectionalSpInfo1D;

clear all

%read atomic weight list, data is from OpenFOAM
AWdict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/AtomicWeight';
fid = fopen(AWdict);
AWData = textscan(fid,'%s','delimiter','\n');
fclose(fid);
AtomW = split(strtrim(AWData{1}));


fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/thermo/thermo.in';
%fileDict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/NAPSchemkin20sections/thermo.in';
fid = fopen(fileDict);

iRow = 1;
while (~feof(fid)) 
    myData(iRow,:) = textscan(fid,'%s','delimiter', '\n','whitespace','');
    iRow = iRow + 1;
end
fclose(fid);

comSty = '!';
TempData= cell(size(myData{1},1),1);
%read data and delte comments
for i=(1:size(myData{1},1))
    if contains(myData{1}{i},comSty)
        TempData{i,:} = extractBefore(myData{1}{i},comSty);
    else
        TempData{i,:} = myData{1}{i};
    end
end

TempData = TempData(~cellfun('isempty',TempData));
cleanData = TempData;

defaultTemperature = sscanf(cleanData{2},'%f');

Nrow = size(cleanData,1);
SpName = cell((Nrow-3)/4,1);
Trange = cell((Nrow-3)/4,1);
TempTrange = cell((Nrow-3)/4,1);
SpElements = cell((Nrow-3)/4,4);
Coeffs = zeros(Nrow,5);
for i = 3:4:(size(cleanData,1)-4) %each species line as start
    SpName{(i+1)/4} = sscanf(cleanData{i}(1:18),'%[^ ]'); %species name
    SpElements{(i+1)/4} = sscanf(cleanData{i}(25:44),'%[^\n]'); %its elements
    if contains(SpElements{(i+1)/4},'0')
        SpElements{(i+1)/4} = extractBefore(SpElements{(i+1)/4},' 0');
    end
    SpElements{(i+1)/4} = erase(SpElements{(i+1)/4},' ');
    SpElements{(i+1)/4} = erase(SpElements{(i+1)/4},'.');
    [SpElements{(i+1)/4,2},SpElements{(i+1)/4,3}] = ...
        regexp(SpElements{(i+1)/4},'[A-Z]*','match','split');
    
    %read temperature range
    TempTrange{(i+1)/4} = sscanf(cleanData{i}(46:75),'%[0-9 .]');
    Trange{(i+1)/4} = regexp(TempTrange{(i+1)/4},'[0-9]*[.]?[0-9]*','match');

    %reading 14 coefficients
    indexE1 = strfind(cleanData{i+1}(1:75),'E');
    indexE2 = strfind(cleanData{i+2}(1:75),'E');
    indexE3 = strfind(cleanData{i+3}(1:75),'E');
    if isempty(indexE1)
        indexE1 = strfind(cleanData{i+1}(1:75),'e');
        indexE2 = strfind(cleanData{i+2}(1:75),'e');
        indexE3 = strfind(cleanData{i+3}(1:75),'e');
    end
    for ci = 1:5
        Coeffs(i+1,ci) = sscanf(cleanData{i+1}(indexE1(ci)-11:indexE1(ci)+3),'%f');
        Coeffs(i+2,ci) = sscanf(cleanData{i+2}(indexE2(ci)-11:indexE2(ci)+3),'%f');  
    end
    for ci = 1:4
        Coeffs(i+3,ci) = sscanf(cleanData{i+3}(indexE3(ci)-11:indexE3(ci)+3),'%f');
    end
   
end

%calculate molecular weight
SpElements(:,4) = {0};
for i = 1:size(SpElements,1)
    for Spi = 1:size(SpElements{i,2},2)
        
        %find element index in atom list
        AtWtIndex = find(strcmpi(AtomW(:,1),SpElements{i,2}{Spi}),1);

        SpElements{i,4} = SpElements{i,4} + str2num(SpElements{i,3}{Spi+1})...
                            *str2num(AtomW{AtWtIndex,2});
    end      

end

%save Species composition information for later use
save('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/foam1D/SpeciesInfo1D','SpName','SpElements')


%find(strcmp(AtomW(:,1),'S'),1)

