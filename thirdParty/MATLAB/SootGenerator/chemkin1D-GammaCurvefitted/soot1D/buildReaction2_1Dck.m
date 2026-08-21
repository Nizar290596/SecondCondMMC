
%build reaction by: ABCD,LumpC,LumpH

%buildReaction(reactions{ri},ABCD,LumpC,LumpH)
fmPrsn = '%.10f';%'%.6g';

%construct reactants
    if reactants{ri}{1,2}~=1
        Re1factor = strcat(num2str(reactants{ri}{1,2}),{' '});
    else
        Re1factor = ' ';
    end
    if reactants{ri}{2,2}~=1
        Re2factor = strcat(num2str(reactants{ri}{2,2}),{' '});
    else
        Re2factor = ' ';
    end
    
    leftLumpSp1 = strcat(RType{1},'_',num2str(lpCRe1,'%02.f'),'_01_01_1');
    leftLumpSp2 = strcat(RType{2},'_',num2str(lpCRe2,'%02.f'),'_01_01_1');
    LEFT = strcat(Re1factor,leftLumpSp1,...
        {' + '},Re2factor,leftLumpSp2);
    
%construct products
RIGHT = extractAfter(reactions{ri},'=>');

for i = 1:2
vertex{i} = strcat(PrdtType,'_',num2str(LumpC(i),'%02.f'),'_01_01_1');
end

LumpSp={' '};
writtenSp = 0;
for li = 1:2
    if AB(li)~=0
        if writtenSp>=1
            LumpSp = strcat(LumpSp,{'+'});
        end
        
        if AB(li)~=1
        LumpSp = strcat(LumpSp,num2str(AB(li),fmPrsn),{' '},vertex{li});
        else
        LumpSp = strcat(LumpSp,{' '},vertex{li}); 
        end
        writtenSp = writtenSp + 1;
        
    end
end

if Hcompensate ~= 0 
    if Hcompensate ~= 2
        LumpSp = strcat(LumpSp,{'+'},num2str(Hcompensate/2),' H2');
    else
        LumpSp = strcat(LumpSp,{'+'},' H2');
    end
end

proRight = strrep(RIGHT,strcat(PrdtType,'t1'),LumpSp);
proRight = strrep(proRight,'+',' + ');


%expand and store
ExpandReactions{ri}{LPi,1} = strcat(LEFT,{' => '},proRight);
clear LEFT
clear RIGHT
