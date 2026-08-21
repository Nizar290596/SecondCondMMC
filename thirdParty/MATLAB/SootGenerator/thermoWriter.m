
%thermo-Foam writer

fid = fopen(...
    '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/converted/gasFOAMthermo','w');
for i = 1:size(SpName,1)
    %print specie{}
    fprintf(fid,'%s\n{\n',SpName{i});
    fprintf(fid,'    specie\n    {\n');
    fprintf(fid,'        molWeight        %f;\n',SpElements{i,4});    
    fprintf(fid,'    }\n'); %end of species
    
    %print thermodynamics{}
    fprintf(fid,'    thermodynamics\n    {\n');
    %lower the temperature range to avoid OF warning
    if str2double(Trange{i}{1}) < 290
        fprintf(fid,'        Tlow        %f;\n',str2double(Trange{i}{1})); 
    else
        fprintf(fid,'        Tlow        %f;\n',290);         
    end
    fprintf(fid,'        Thigh       %f;\n',str2double(Trange{i}{2})); 
    if size(Trange{i},2)==3
        fprintf(fid,'        Tcommon     %f;\n',str2num(Trange{i}{3})); 
    else
        fprintf(fid,'        Tcommon     %f;\n',(defaultTemperature(2)));         
    end
    
    fprintf(fid,'        highCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
        Coeffs(4*i,1),Coeffs(4*i,2),Coeffs(4*i,3),Coeffs(4*i,4),...
        Coeffs(4*i,5),Coeffs(4*i+1,1),Coeffs(4*i+1,2)); %end of highCp
    
    fprintf(fid,'        lowCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
        Coeffs(4*i+1,3),Coeffs(4*i+1,4),Coeffs(4*i+1,5),Coeffs(4*i+2,1),...
        Coeffs(4*i+2,2),Coeffs(4*i+2,3),Coeffs(4*i+2,4)); %end of highCp
    
    fprintf(fid,'    }\n'); %end of thermodynamics
    
    %print transport{}
    fprintf(fid,'    transport\n    {\n');
    fprintf(fid,'        As        %f;\n', 1.67212e-06);    
    fprintf(fid,'        Ts        %f;\n', 170.672);        
    
    fprintf(fid,'    }\n'); %end of transport

    %print elements{}
    fprintf(fid,'    elements\n    {\n');
    
    for SpInx = 1:size(SpElements{i,2},2)
        fprintf(fid,'        %s        %d;\n',...
            SpElements{i,2}{SpInx},str2num(SpElements{i,3}{SpInx+1}));    
    end
    fprintf(fid,'    }\n'); %end of elements

    fprintf(fid,'}\n\n');

end

%nC = 20;
%nH = 2;
i = 154; %A4 index - pyrene- thermodynamic data scale as WtAi/WtA4 * coeffA4
wtA4 = SpElements{i,4};
load('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/SectionInfo');

%For Ai
for nCi = 1:size(nC,1)
    for nHi = 1:size(nH,2)
        %calculate scaling factor
        Awt = (nC(nCi)*str2double(AtomW{10,2})+nH(nCi,nHi)*str2double(AtomW{3,2}));
        sclFcA4 = Awt/wtA4;
        
        %print specie{}
        fprintf(fid,'%s\n{\n',strcat('Ai_',num2str(nCi,'%02.f'),'_',num2str(nHi,'%02.f'),'_01_1'));
        fprintf(fid,'    specie\n    {\n');
        fprintf(fid,'        molWeight        %f;\n',Awt);    
        fprintf(fid,'    }\n'); %end of species
        
        %print thermodynamics{}
        fprintf(fid,'    thermodynamics\n    {\n');
        fprintf(fid,'        Tlow        %f;\n',290);         
        fprintf(fid,'        Thigh       %f;\n',str2double(Trange{i}{2})); 
        if size(Trange{i},2)==3
            fprintf(fid,'        Tcommon     %f;\n',str2num(Trange{i}{3})); 
        else
            fprintf(fid,'        Tcommon     %f;\n',(defaultTemperature(2)));         
        end
        
        fprintf(fid,'        highCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
            sclFcA4*Coeffs(4*i,1),sclFcA4*Coeffs(4*i,2),sclFcA4*Coeffs(4*i,3),sclFcA4*Coeffs(4*i,4),...
            sclFcA4*Coeffs(4*i,5),sclFcA4*Coeffs(4*i+1,1),sclFcA4*Coeffs(4*i+1,2)); %end of highCp

        fprintf(fid,'        lowCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
            sclFcA4*Coeffs(4*i+1,3),sclFcA4*Coeffs(4*i+1,4),sclFcA4*Coeffs(4*i+1,5),sclFcA4*Coeffs(4*i+2,1),...
            sclFcA4*Coeffs(4*i+2,2),sclFcA4*Coeffs(4*i+2,3),sclFcA4*Coeffs(4*i+2,4)); %end of highCp

        fprintf(fid,'    }\n'); %end of thermodynamics
        
        %print transport{}
        fprintf(fid,'    transport\n    {\n');
        fprintf(fid,'        As        %f;\n', 1.67212e-06);    
        fprintf(fid,'        Ts        %f;\n', 170.672);        

        fprintf(fid,'    }\n'); %end of transport

        %print elements{}
        fprintf(fid,'    elements\n    {\n');

        fprintf(fid,'        %s        %d;\n','C',nC(nCi));    
        fprintf(fid,'        %s        %d;\n','H',nH(nCi,nHi));    
        
        fprintf(fid,'    }\n'); %end of elements

        fprintf(fid,'}\n\n');

    end
end

%for Aj
for nCi = 1:size(nC,1)
    for nHi = 1:size(nHaj,2)
        Awt = (nC(nCi)*str2double(AtomW{10,2})+nHaj(nCi,nHi)*str2double(AtomW{3,2}));
        sclFcA4 = Awt/wtA4;
        
        %print specie{}
        fprintf(fid,'%s\n{\n',strcat('Aj_',num2str(nCi,'%02.f'),'_',num2str(nHi,'%02.f'),'_01_1'));
        fprintf(fid,'    specie\n    {\n');
        fprintf(fid,'        molWeight        %f;\n',Awt);    
        fprintf(fid,'    }\n'); %end of species
        
        %print thermodynamics{}
        fprintf(fid,'    thermodynamics\n    {\n');
        fprintf(fid,'        Tlow        %f;\n',290);         
        fprintf(fid,'        Thigh       %f;\n',str2double(Trange{i}{2})); 
        if size(Trange{i},2)==3
            fprintf(fid,'        Tcommon     %f;\n',str2num(Trange{i}{3})); 
        else
            fprintf(fid,'        Tcommon     %f;\n',(defaultTemperature(2)));         
        end

        fprintf(fid,'        highCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
            sclFcA4*Coeffs(4*i,1),sclFcA4*Coeffs(4*i,2),sclFcA4*Coeffs(4*i,3),sclFcA4*Coeffs(4*i,4),...
            sclFcA4*Coeffs(4*i,5),sclFcA4*Coeffs(4*i+1,1),sclFcA4*Coeffs(4*i+1,2)); %end of highCp

        fprintf(fid,'        lowCpCoeffs      ( %.8e %.8e %.8e %.8e %.8e %.8e %.8e );\n',...
            sclFcA4*Coeffs(4*i+1,3),sclFcA4*Coeffs(4*i+1,4),sclFcA4*Coeffs(4*i+1,5),sclFcA4*Coeffs(4*i+2,1),...
            sclFcA4*Coeffs(4*i+2,2),sclFcA4*Coeffs(4*i+2,3),sclFcA4*Coeffs(4*i+2,4)); %end of highCp

        fprintf(fid,'    }\n'); %end of thermodynamics
        
        %print transport{}
        fprintf(fid,'    transport\n    {\n');
        fprintf(fid,'        As        %f;\n', 1.67212e-06);    
        fprintf(fid,'        Ts        %f;\n', 170.672);        

        fprintf(fid,'    }\n'); %end of transport

        %print elements{}
        fprintf(fid,'    elements\n    {\n');

        fprintf(fid,'        %s        %d;\n','C',nC(nCi));    
        fprintf(fid,'        %s        %d;\n','H',nHaj(nCi,nHi));    
        
        fprintf(fid,'    }\n'); %end of elements

        fprintf(fid,'}\n\n');

    end
end

fclose('all');

