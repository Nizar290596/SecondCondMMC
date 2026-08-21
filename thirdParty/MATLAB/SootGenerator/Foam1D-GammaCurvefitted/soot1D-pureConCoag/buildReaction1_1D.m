
%build reaction by: ABCD,LumpC,LumpH

%buildReaction(reactions{ri},ABCD,LumpC,LumpH)
fmPrsn = '%.6g';

%construct reactants
%LEFT = extractBefore(reactions{ri},'=');
if size(reactants{ri},1)==2
    if reactants{ri}{skipReactant,2}~=1
        skRfactor = strcat(num2str(reactants{ri}{skipReactant,2}),{' '});
    else
        skRfactor = ' ';
    end
    if reactants{ri}{takeReactant,2}~=1
        tkRfactor = strcat(num2str(reactants{ri}{takeReactant,2}),{' '});
    else
        tkRfactor = ' ';
    end
    
    leftLumpSp = strcat(RType,'_',num2str(lpC,'%02.f'),'_01_01_1');
    LEFT = strcat(skRfactor,leftLumpSp,...
        {' + '},tkRfactor,reactants{ri}{takeReactant,3});
    
elseif size(reactants{ri},1)==1
    if reactants{ri}{skipReactant,2}~=1
        skRfactor = strcat(num2str(reactants{ri}{skipReactant,2}),{' '});
    else
        skRfactor = ' ';
    end
    leftLumpSp = strcat(RType,'_',num2str(lpC,'%02.f'),'_01_01_1');
    LEFT = strcat(skRfactor,leftLumpSp);
end

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
            LumpSp = strcat(LumpSp,{' + '});
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
ExpandReactions{ri}{LPi,1} = strcat(LEFT,{' = '},proRight);
clear LEFT
clear RIGHT
