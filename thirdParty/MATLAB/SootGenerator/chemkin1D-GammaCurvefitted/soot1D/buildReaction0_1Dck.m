
%build reaction by: ABCD,LumpC,LumpH

%buildReaction(reactions{ri},ABCD,LumpC,LumpH)
fmPrsn = '%.10f';%'%.6g';

LEFT = extractBefore(reactions{ri},'=');

LEFT = strrep(LEFT,'+',' + ');

RIGHT = extractAfter(reactions{ri},'=>');

for i = 1:2
vertex{i} = strcat(PrdtType,'_',num2str(LumpC(i),'%02.f'),'_01_01_1');
end

LumpSp={' '};
writtenSp = 0;
for li = 1:2
    if AB(li)~=0
        if writtenSp>=1
            LumpSp = strcat(LumpSp,{' + '});
        end
        
        LumpSp = strcat(LumpSp,num2str(AB(li),fmPrsn),{' '},vertex{li});
        writtenSp = writtenSp + 1;
    end
end

if Hcompensate ~= 0 
    LumpSp = strcat(LumpSp,{' + '},num2str(Hcompensate/2),' H2');
end

proRight = strrep(RIGHT,strcat(PrdtType,'t1'),LumpSp);

ExpandReactions{ri}{LPi,1} = strcat(LEFT,{' => '},proRight);
clear LEFT
clear RIGHT
