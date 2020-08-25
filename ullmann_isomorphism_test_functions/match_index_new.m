function [ match_all ] = match_index_new( M, A_pattern, A_molecule, labels_pattern, labels_molecule )
%MATCH_INDEX Summary of this function goes here
match_all = {};
for i = 1 : size(M,1) % count pattern
    for j = 1 : size(M,2) % count molecule
        if M(i, j)~= 0
            matched = zeros(length(A_pattern),2);
            root_pattern = i;
            root_molecule = j;
            matched(:, 1) = 1:length(A_pattern); 
            matched(1, 2) = root_molecule;
            order = b_first_search( A_pattern, root_pattern);
            
            for k = 2 : length(order)
                
                candidates_ = find(M(matched(order(k),1),:)~=0);
                for l = 1 : length(candidates_) % neighbors of the atom to add are already identified
                    % find is atom k is connected to the previously identif
                    neighbors_candidate_mol = neighbors(A_molecule, candidates_(l));
                    neighbors_candidate_pattern = intersect(neighbors(A_pattern, matched(order(k), 1)), matched(order(1:k),1));
                    %                     if ~isempty(intersect( neighbors_candidate_mol, matched(find(matched(:,1)==intersect(neighbors(A_pattern, matched(order(k), 1)), matched(1:order(k),1))), 2))) && isempty(intersect(candidates_(l), matched(:,2)))
                    %                         matched(k,2) = candidates_(l);
                    %                     end
                    if ~isempty(intersect( neighbors_candidate_mol, matched(intersect(neighbors(A_pattern, matched(order(k), 1)), matched(:,1)), 2))) && isempty(intersect(candidates_(l), matched(:,2)))
                        matched(order(k),2) = candidates_(l);
                    end
                end
            end
            
            if isempty(find(matched(:,2)==0))
                % try to validate connections between carbons
                c_atom_pattern = find(ismember(labels_pattern, 'C'));
                single_bond = 0;
                double_bond = 0;
                if length(c_atom_pattern)>1
                    [~, b] = intersect(matched(:,1), c_atom_pattern);
                    single_bond = length(find(tril(A_molecule(matched(b,2), matched(b,2))) == 1)) == length(find(tril(A_pattern(matched(b,1), matched(b,1))) == 1));
                    double_bond = length(find(tril(A_molecule(matched(b,2), matched(b,2))) == 2)) == length(find(tril(A_pattern(matched(b,1), matched(b,1))) == 2));
                    
                end
                if single_bond && double_bond || length(c_atom_pattern)<=1
                    indicator = 0;
                    if isempty(match_all)
                        match_all = [match_all; matched];
                        indicator = 1;
                    else
                        for k = 1 : length(match_all)
                            common_atoms = intersect(match_all{k}(:,2), matched);
                            uncommon_atoms = setdiff(matched(:,2), common_atoms);
                            if length(common_atoms) == length(A_pattern)
                                indicator = 1; % if indicator is 1 it means that this match is already there
                                % break
                            else if length(uncommon_atoms) == 1
                                    if ismember(labels_molecule(uncommon_atoms), 'H')
                                       indicator = 1; 
                                    end
                                end
                            end
                        end
                        
                        if ~indicator
                            match_all = [match_all; matched];
                        end
                    end
                end
            end
        end
    end
end

