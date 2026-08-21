%read atomic weight list, data is from OpenFOAM
AWdict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/AtomicWeight';
fid = fopen(AWdict);
AWData = textscan(fid,'%s','delimiter','\n');
fclose(fid);
AtomW = split(strtrim(AWData{1}));
%thermo-Foam writer

nH=2;
nHj=1;
fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D/converted/Section','w');

%nC = 20;
%nH = 2;
i = 154; %A4 index - pyrene- thermodynamic data scale as WtAi/WtA4 * coeffA4
wtA4 = SpElements{i,4};
load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D/SectionInfo1D');
load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D/SpeciesInfo1D');

%For Ai
for nCi = 1:size(nC,1)
        %calculate scaling factor
        Awt = nC(nCi)*str2double(AtomW{10,2}) + nH;
        sclFcA4 = Awt/wtA4;
        
        fprintf(fid,'%s    ',strcat('BINi-',num2str(nCi,'%02.f')));
        

end

        fprintf(fid,'   \n');         

%for Aj
for nCi = 1:size(nC,1)
        Awt = nC(nCi)*str2double(AtomW{10,2})+1;
        sclFcA4 = Awt/wtA4;
        
        %print specie{}
        %fprintf(fid,'%s',strcat('BINj-',num2str(nCi,'%02.f')));
        fprintf(fid,'%s    ',strcat('BINj-',num2str(nCi,'%02.f')));
        



end

fclose('all');

