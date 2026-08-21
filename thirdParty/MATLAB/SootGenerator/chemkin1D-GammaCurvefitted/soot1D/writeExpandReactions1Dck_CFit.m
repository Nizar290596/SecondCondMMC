%write Expanded reactions

fid = fopen('/Users/zhijiehuo/Dropbox/SootGenerator/chemkin1D-GammaCurvefitted/converted/secFOAM3V3B1Dck-T-Zhijiecase129','w');
%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/converted/FraFixGamma1/3VsecFOAMfix1','w');
%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/testing/secFOAM','w');
%gamma = 1;

%starting number of sectional reaction
rnNo = 243;
for EXRi = 1:size(ExpandReactions,1)
    for EXRj = 1:size(ExpandReactions{EXRi},1)
        %fprintf(fid,'    un-named-reaction-%d\n    {\n',rnNo);
        %fprintf(fid,'        type            ');
        %fprintf(fid,'irreversibleArrheniusReaction;\n');
        %fprintf(fid,'        reaction        ');
        fprintf(fid,'%s         ',ExpandReactions{EXRi}{EXRj}{1});


        %fprintf(fid,'        A            ');
        if FoamArrCoeffs{EXRi}(EXRj,7) == 2
            %coag = Beta * Gamma
            %GammaArrConst{EXRi}(EXRj,1)=GammaArrConst{EXRi}(EXRj,5);%0.3;
            %GammaArrConst{EXRi}(EXRj,2)=0;
            %GammaArrConst{EXRi}(EXRj,3)=0;
            
            fprintf(fid,'%e    ',FoamArrCoeffs{EXRi}(EXRj,1)*GammaArrConst{EXRi}(EXRj,1));

            fprintf(fid,'%f    ',FoamArrCoeffs{EXRi}(EXRj,2)+GammaArrConst{EXRi}(EXRj,2));
            fprintf(fid,'%f    ',GammaArrConst{EXRi}(EXRj,3)*8.314/4.184); %Ta(K) -> Ea(cal/mol)
            fprintf(fid,'\n');        
        
        else
            % sootA = (m*C)^k
            sootA = (FoamArrCoeffs{EXRi}(EXRj,5) * FoamArrCoeffs{EXRi}(EXRj,6))^(FoamArrCoeffs{EXRi}(EXRj,4));
            
            if sootA == 0 %nucleation m = 0
                fprintf(fid,'%e    ',FoamArrCoeffs{EXRi}(EXRj,1));
            else
                fprintf(fid,'%e    ',FoamArrCoeffs{EXRi}(EXRj,1)*sootA);
            end
            
            fprintf(fid,'%f    ',FoamArrCoeffs{EXRi}(EXRj,2));
            fprintf(fid,'%f    ',FoamArrCoeffs{EXRi}(EXRj,3));
            fprintf(fid,'\n');
            
        end

        rnNo = rnNo +1;
    end
end


fclose('all');
