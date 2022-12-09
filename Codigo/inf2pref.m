function [resposta]=inf2pref(reg)

load ind
reg

resposta=cell(1);

i=1;
for j=1:length(reg)

    asteristicos=find(reg{i,j}=='*')
    
    if isempty(asteristicos)
        
    else
%         resposta=
        base=1;
        for k=1:length(asteristicos)
            reg{i,j}(1,base:(asteristicos(k)-1))
            base=asteristicos(k)+1
        end
        reg{i,j}(1,base:length(reg{i,j}))
    end
    
    
end

end