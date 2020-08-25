function [ match_all ] = ullmann_pattern( A_pattern, labels_pattern, crosslinks_pattern, A_molecule, labels_molecule, crosslinks_molecule )
% ULLMANN_PATTERN finds patterns in a molecule

M = initialize_M_pattern(A_pattern, labels_pattern, crosslinks_pattern, A_molecule, labels_molecule, crosslinks_molecule);
match_all = {};
% update matrix till match is found or one of the rows gets all zeros: it
% means that there is no match for graphs
[M, dont_match_1] = first_update_M(M, A_pattern, A_molecule );
% refine matrix M for each non_zero entry
if ~dont_match_1
    % find candidates
    [M, dont_match_2] = second_update_M(M, A_pattern, A_molecule );
    if ~dont_match_2
        [ match_all ] = match_index_new( M, A_pattern, A_molecule, labels_pattern, labels_molecule ); % these labels are needed only for excluding symmetric hydrogens
    end
end
end

