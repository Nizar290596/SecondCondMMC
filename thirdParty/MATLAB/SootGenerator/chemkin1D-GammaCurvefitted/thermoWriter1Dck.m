%read atomic weight list, data is from OpenFOAM
AWdict = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/AtomicWeight';
fid = fopen(AWdict);
AWData = textscan(fid,'%s','delimiter','\n');
fclose(fid);
AtomW = split(strtrim(AWData{1}));
%thermo-Foam writer

nH=2;
nHj=1;
fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D-GammaCurvefitted/converted/thermo1Dck80','w');

%nC = 20;
%nH = 2;
i = 154; %A4 index - pyrene- thermodynamic data scale as WtAi/WtA4 * coeffA4
wtA4 = SpElements{i,4};
%load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D/SectionInfo1D');
nC=load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D-GammaCurvefitted/sectionC80');
load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/chemkin1D/SpeciesInfo1D');

%For Ai
for nCi = 1:size(nC,1)
        %calculate scaling factor
        Awt = nC(nCi)*str2double(AtomW{10,2}) + nH;
        sclFcA4 = Awt/wtA4;
        
        fprintf(fid,'%s                 ',strcat('BINi-',num2str(nCi,'%02.f')));
        
        %fprintf(fid,'%s  %d','C',nC(nCi));         
        fprintf(fid,'%s   %d','H',nH);         
        fprintf(fid,'               %s','G');         
        
        fprintf(fid,'   %s','200.000');         
        fprintf(fid,'  %s             1','3000.0000'); 
        fprintf(fid,'%s   %d','C',nC(nCi)); 
        
        fprintf(fid,'\n %.8e %.8e %.8e %.8e %.8e    2\n %.8e %.8e',...
            sclFcA4*Coeffs(4*i,1),sclFcA4*Coeffs(4*i,2),sclFcA4*Coeffs(4*i,3),sclFcA4*Coeffs(4*i,4),...
            sclFcA4*Coeffs(4*i,5),sclFcA4*Coeffs(4*i+1,1),sclFcA4*Coeffs(4*i+1,2)); %end of highCp
        

        fprintf(fid,'%.8e %.8e %.8e    3\n %.8e %.8e %.8e %.8e                   4\n',...
            sclFcA4*Coeffs(4*i+1,3),sclFcA4*Coeffs(4*i+1,4),sclFcA4*Coeffs(4*i+1,5),sclFcA4*Coeffs(4*i+2,1),...
            sclFcA4*Coeffs(4*i+2,2),sclFcA4*Coeffs(4*i+2,3),sclFcA4*Coeffs(4*i+2,4)); %end of highC        
 

end

        fprintf(fid,'   \n');         

%for Aj
for nCi = 1:size(nC,1)
        Awt = nC(nCi)*str2double(AtomW{10,2})+1;
        sclFcA4 = Awt/wtA4;
        
        %print specie{}
        %fprintf(fid,'%s',strcat('BINj-',num2str(nCi,'%02.f')));
        fprintf(fid,'%s                 ',strcat('BINj-',num2str(nCi,'%02.f')));
        
        fprintf(fid,'%s   %d','H',nHj);         
        fprintf(fid,'               %s','G');         
        
        fprintf(fid,'   %s','200.000');         
        fprintf(fid,'   %s            1','3000.0000'); 
        fprintf(fid,'%s   %d','C',nC(nCi));         


        fprintf(fid,'\n %.8e %.8e %.8e %.8e %.8e    2\n %.8e %.8e',...
            sclFcA4*Coeffs(4*i,1),sclFcA4*Coeffs(4*i,2),sclFcA4*Coeffs(4*i,3),sclFcA4*Coeffs(4*i,4),...
            sclFcA4*Coeffs(4*i,5),sclFcA4*Coeffs(4*i+1,1),sclFcA4*Coeffs(4*i+1,2)); %end of highCp

        fprintf(fid,'%.8e %.8e %.8e    3\n %.8e %.8e %.8e %.8e                   4\n',...
            sclFcA4*Coeffs(4*i+1,3),sclFcA4*Coeffs(4*i+1,4),sclFcA4*Coeffs(4*i+1,5),sclFcA4*Coeffs(4*i+2,1),...
            sclFcA4*Coeffs(4*i+2,2),sclFcA4*Coeffs(4*i+2,3),sclFcA4*Coeffs(4*i+2,4)); %end of highCp       

        

end

fclose('all');

