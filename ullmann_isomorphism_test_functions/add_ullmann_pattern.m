function [ matched, ptrn ] = add_ullmann_pattern( A_molecule, labels_molecule, crosslinks_molecule, pattern, j )
matched = {};
ptrn = [];

if length(A_molecule)>5
        match_int = {};
        if length(pattern.adj)>=2
            match_int = ullmann_pattern( pattern.adj, pattern.lbl, pattern.crsl, A_molecule, labels_molecule, crosslinks_molecule );
        end
        if ~isempty(match_int)
            matched = [matched; match_int];
            for k = 1 : length(match_int)
                ptrn = [ptrn; j];
            end
        end
    
end

end





