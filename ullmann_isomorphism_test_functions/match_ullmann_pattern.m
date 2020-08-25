function [ matched, ptrn ] = match_ullmann_pattern( A_molecule, labels_molecule, crosslinks_molecule )
global pattern;
matched = {};
ptrn = [];

if length(A_molecule)<2
    for j = 1 : length(pattern)
        if length(pattern(j).adj)<2
            if isequal(pattern(j).lbl, labels_molecule) && isequal(pattern(j).lbl, labels_molecule)
                ptrn = [ptrn; j];
            end
        end
    end
    
    
else
    for j = 1 : length(pattern)
        match_int = {};
        if length(pattern(j).adj)>=2
            
            match_int = ullmann_pattern( pattern(j).adj, pattern(j).lbl, pattern(j).crsl, A_molecule, labels_molecule, crosslinks_molecule );
        end
        if ~isempty(match_int)
            matched = [matched; match_int];
            for k = 1 : length(match_int)
                ptrn = [ptrn; j];
            end
        end
    end
    
end

end





