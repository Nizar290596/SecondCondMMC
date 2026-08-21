
%start writing

fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/Foam1D/converted/gasFOAM','w');
%fid = fopen('/Users/zhijiehuo/Dropbox/DelftFlameResults/FOAMformat/NAPSchemkin20sections/sootFOAMNAPS','w');
%print elements
fprintf(fid, 'elements\n%d\n',size(elementsList,2));
fprintf(fid,'(\n');
for i=1:size(elementsList,2)
    fprintf(fid,'%s\n',elementsList{i});
end
fprintf(fid,')\n;\n\n');

SpNo = size(SpeciesNames,1);
%print species
fprintf(fid, 'species\n%d\n',SpNo);
fprintf(fid,'(\n');
for i=1:SpNo
    fprintf(fid,'%s\n',SpeciesNames{i});
end
fprintf(fid,')\n;\n\n');


%print reactions
fprintf(fid, 'reactions\n{\n');

rnNo = 0;
for i=1:size(rectns,1)
    if ~isempty(FOAMreactions{i})
        fprintf(fid,'    un-named-reaction-%d\n    {\n',rnNo);
        rnNo = rnNo + 1;
        fprintf(fid,'        type            ');
        fprintf(fid,'%s;\n',strcat(reactionType{i,1},reactionType{i,2}));
        fprintf(fid,'        reaction        ');
        fprintf(fid,'"%s";\n',FOAMreactions{i}); 
        
        %print A beta Ta
        if strcmp(reactionType{i,2},'ArrheniusReaction') ...
              || strcmp(reactionType{i,2},'thirdBodyArrheniusReaction')

            fprintf(fid,'        A            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i,1));
            fprintf(fid,'        beta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,2));
            fprintf(fid,'        Ta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,3));
        end

        %print third body coeffs
        if strcmp(reactionType{i,2},'thirdBodyArrheniusReaction')
            fprintf(fid,'        coeffs            \n%d\n(\n',SpNo);
            %3rd body coeffs r all 1
            if isempty(reactionType{i+1,3})
                for nSp = 1:SpNo
                    fprintf(fid,'(%s 1) ',SpeciesNames{nSp});
                end
            %find 3rd body coeff that are not 1 in reactionType{i+1,3)
            else 
                for nSp = 1:SpNo
                    if sum(strcmp(reactionType{i+1,3},SpeciesNames{nSp}))
                        index = find(strcmp(reactionType{i+1,3},SpeciesNames{nSp}));
                        fprintf(fid,'(%s %s) ',SpeciesNames{nSp},reactionType{i+1,3}{index+1});
                    else
                        fprintf(fid,'(%s 1) ',SpeciesNames{nSp});
                    end
                    
                end
            end
            fprintf(fid,'\n)\n;\n');%end of coeffs 
        end
        
        %print lindemann 
        if strcmp(reactionType{i,2},'ArrheniusLindemannFallOffReaction')
            fprintf(fid,'        k0\n        {\n'); %print k0
            fprintf(fid,'            A            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i+1,1));
            fprintf(fid,'            beta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+1,2));
            fprintf(fid,'            Ta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+1,3));
            fprintf(fid,'        }\n');    
            
            fprintf(fid,'        kInf\n        {\n'); %print kf
            fprintf(fid,'            A            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i,1));
            fprintf(fid,'            beta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,2));
            fprintf(fid,'            Ta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,3));
            fprintf(fid,'        }\n');
            fprintf(fid,'        F\n        {\n'); %print F, empty
            fprintf(fid,'        }\n');

            fprintf(fid,'        thirdBodyEfficiencies\n');
            fprintf(fid,'        {\n');
            fprintf(fid,'        coeffs            \n%d\n(\n',SpNo);
            
            for nSp = 1:SpNo
                    if sum(strcmp(reactionType{i+2,3},SpeciesNames{nSp}))
                        index = find(strcmp(reactionType{i+2,3},SpeciesNames{nSp}));
                        fprintf(fid,'(%s %s) ',SpeciesNames{nSp},reactionType{i+2,3}{index+1});
                    else
                        fprintf(fid,'(%s 1) ',SpeciesNames{nSp});
                    end
            end            
            fprintf(fid,'\n)\n;\n');%end of coeffs 
            fprintf(fid,'        }\n');%end of thirdBodyEfficiencies
        end
        
        %print Troe 
        if strcmp(reactionType{i,2},'ArrheniusTroeFallOffReaction')
            fprintf(fid,'        k0\n        {\n'); %print k0
            fprintf(fid,'            A            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i+1,1));
            fprintf(fid,'            beta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+1,2));
            fprintf(fid,'            Ta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+1,3));
            fprintf(fid,'        }\n');    
            
            fprintf(fid,'        kInf\n        {\n'); %print kf
            fprintf(fid,'            A            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i,1));
            fprintf(fid,'            beta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,2));
            fprintf(fid,'            Ta            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i,3));
            fprintf(fid,'        }\n');
            fprintf(fid,'        F\n        {\n'); %print F, 4p
            fprintf(fid,'            alpha            ');
            fprintf(fid,'%e;\n',FoamArrhenCoef(i+2,1));
            fprintf(fid,'            Tsss            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+2,2));
            fprintf(fid,'            Ts            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+2,3));
            fprintf(fid,'            Tss            ');
            fprintf(fid,'%f;\n',FoamArrhenCoef(i+2,4));
            fprintf(fid,'        }\n');

            fprintf(fid,'        thirdBodyEfficiencies\n');
            fprintf(fid,'        {\n');
            fprintf(fid,'        coeffs            \n%d\n(\n',SpNo);
            
            for nSp = 1:SpNo
                    if sum(strcmp(reactionType{i+3,3},SpeciesNames{nSp}))
                        index = find(strcmp(reactionType{i+3,3},SpeciesNames{nSp}));
                        fprintf(fid,'(%s %s) ',SpeciesNames{nSp},reactionType{i+3,3}{index+1});
                    else
                        fprintf(fid,'(%s 1) ',SpeciesNames{nSp});
                    end
            end            
            fprintf(fid,'\n)\n;\n');%end of coeffs 
            fprintf(fid,'        }\n');%end of thirdBodyEfficiencies
        end
        
        fprintf(fid,'    }\n');%end of this reaction
    end
end


fprintf(fid,'}\n\n');%end of reaction section

fclose(fid);


