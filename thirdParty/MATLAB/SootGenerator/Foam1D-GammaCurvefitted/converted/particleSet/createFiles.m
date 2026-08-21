clear all

for i=1:20
    sootSp{i,1} = strcat('YAi_',num2str(i,'%02.f'),'_01_01_1');
end
for i=1:20
    sootSp{i+20,1} = strcat('YAj_',num2str(i,'%02.f'),'_01_01_1');
end

locDir = '/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D-GammaCurvefitted/converted/particleSet/';

for i = 1: size(sootSp,1)
    
    fname = strcat(locDir,sootSp{i});
    
    fid = fopen(fname,'w');
    
    fprintf(fid,'%s\n','/*--------------------------------*- C++ -*----------------------------------*\');
    fprintf(fid,'%s\n','| =========                 |                                                 |');
    fprintf(fid,'%s\n','| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |');
    fprintf(fid,'%s\n','|  \\    /   O peration     | Version:  5.0                                   |');
    fprintf(fid,'%s\n','|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |');
    fprintf(fid,'%s\n','|    \\/     M anipulation  |                                                 |');
    fprintf(fid,'%s\n','\*---------------------------------------------------------------------------*/');
    %fprintf(fid,'%s\n','');
    %fprintf(fid,'%s\n','');
    %fprintf(fid,'%s\n','');
    
    
   
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    
    fprintf(fid,'FoamFile\n{\n');
    
    fprintf(fid,'    version     2.0;\n');
    fprintf(fid,'    format      ascii;\n');
    fprintf(fid,'    class       scalarField;\n');
    fprintf(fid,'    location    "0.306/lagrangian/popeParticlesSet_1";\n');
    fprintf(fid,'    object      %s;\n',sootSp{i});
    fprintf(fid,'}\n\n');
    
    fprintf(fid,'%s\n','// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //');

    
    fprintf(fid,'933943{0}\n');
    
    fprintf(fid,'%s\n','// ************************************************************************* //');
    
    
    fclose(fid);
    
    
end

%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D-GammaCurvefitted/converted/particleSet/','w');
