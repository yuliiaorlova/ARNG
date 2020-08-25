function [ idx_prod ] = add_molecule( A, labels, crosslinks )
% this function checks whether a newly created molecule is already in
% the list of species or not. Additionally, if molecule is not in the list of Species yet,
% it finds patterns in a molecule and saves it to Species
global worth_updating;
global Species;
idx_prod = [];
if length(A)<3
    for i = 1 : length(Species)
        if length(Species(i))<=3 && isequal(Species(i).lbl, labels)
            idx_prod = i;
            break
        end
    end
    if isempty(idx_prod)
        Species(end + 1).adj = A;
        Species(end).lbl = labels;
        Species(end).crsl = crosslinks;
        % here we need to match patterns
        [ matched, ptrn ] = match_ullmann_pattern( A, labels, crosslinks ); % find patterns in newly created molecule
        Species(end).match = matched;
        Species(end).ptrn = ptrn;
        Species(end).stuff = 1;
        idx_prod = length(Species);
        worth_updating = 1;
        
    end
else
    stuff = get_stuff(A, labels);
    % matching molecule
    
    for i = 1 : length(Species)
        
        if length(Species(i).adj)>=3
            
            if isequal(stuff, Species(i).stuff) %n1==n2 && O1==O2 && H1==H2 && C1==C2
                g1 = sparse(Species(i).adj);
                g2 = sparse(A);
                labels_1 = Species(i).lbl;
                labels_2 = labels;
                crosslinks_1 = Species(i).crsl;
                crosslinks_2 = crosslinks;
                [indicator, map] = graphisomorphism(g1,g2);
                if indicator
                    if isequal(sort(labels_1), sort(labels_2)) && isequal(crosslinks_1, crosslinks_2(map))
                        idx_prod = i;
                    end
                end
            end
        end
    end
end
%[ idx_prod ] = match_ullmann_molecule( A, labels, stuff ); % main check: if the molecule is already in the list
if isempty(idx_prod)
    Species(end + 1).adj = A;
    Species(end).lbl = labels;
    Species(end).crsl = crosslinks;
    [ matched, ptrn ] = match_ullmann_pattern( A, labels, crosslinks ); % find patterns in newly created molecule
    Species(end).match = matched;
    Species(end).ptrn = ptrn;
    Species(end).stuff = get_stuff(A, labels);
    idx_prod = length(Species);
    worth_updating = 1;
end


end

