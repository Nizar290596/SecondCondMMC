%write Expanded reactions

%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D-GammaCurvefitted/converted/Tvar/secFOAM3V3B1D-Tvar','w');
fid = fopen('/Users/zhijiehuo/Dropbox/SootGenerator/Foam1D-GammaCurvefitted/changeKinetics/secFOAM3V3B1Dck-T-Zhijiecase944444','w');
%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/converted/FraFixGamma1/3VsecFOAMfix1','w');
%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/testing/secFOAM','w');
%gamma = 1;

%starting number of sectional reaction
rnNo = 243;
for EXRi = 1:size(ExpandReactions,1)
    for EXRj = 1:size(ExpandReactions{EXRi},1)
        fprintf(fid,'    un-named-reaction-%d\n    {\n',rnNo);
        fprintf(fid,'        type            ');
        fprintf(fid,'irreversibleArrheniusReaction;\n');
        fprintf(fid,'        reaction        ');
        fprintf(fid,'"%s";\n',ExpandReactions{EXRi}{EXRj}{1});


        fprintf(fid,'        A            ');
        if FoamArrCoeffs{EXRi}(EXRj,7) == 2
            fprintf(fid,'%e;\n',FoamArrCoeffs{EXRi}(EXRj,1)*GammaArrConst{EXRi}(EXRj,1));
            fprintf(fid,'        beta            ');
            fprintf(fid,'%f;\n',FoamArrCoeffs{EXRi}(EXRj,2)+GammaArrConst{EXRi}(EXRj,2));
            fprintf(fid,'        Ta            ');
            fprintf(fid,'%f;\n',GammaArrConst{EXRi}(EXRj,3)); %Ta(K)
            fprintf(fid,'    }\n');            
        else
            % sootA = A*(mC)^k
            sootA = (FoamArrCoeffs{EXRi}(EXRj,5) * FoamArrCoeffs{EXRi}(EXRj,6))^(FoamArrCoeffs{EXRi}(EXRj,4));
            
            if sootA == 0 %nucleation m = 0
                fprintf(fid,'%e;\n',FoamArrCoeffs{EXRi}(EXRj,1));
            else
                fprintf(fid,'%e;\n',FoamArrCoeffs{EXRi}(EXRj,1)*sootA);
            end
            
            fprintf(fid,'        beta            ');
            fprintf(fid,'%f;\n',FoamArrCoeffs{EXRi}(EXRj,2));
            fprintf(fid,'        Ta            ');
            fprintf(fid,'%f;\n',FoamArrCoeffs{EXRi}(EXRj,3));
            fprintf(fid,'    }\n');  
            
        end
        


        rnNo = rnNo +1;
    end
end


fclose('all');
