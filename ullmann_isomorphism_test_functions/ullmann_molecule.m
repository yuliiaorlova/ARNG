function [ indicator ] = ullmann_molecule( A_pattern, labels_pattern, A_molecule, labels_molecule )
%ULLMANN_NOLECULE checks if two molecular graphs are isomorphic
indicator = 0;
M = initialize_M(A_pattern, labels_pattern, A_molecule, labels_molecule);
% update matrix till match is found or one of the rows gets all zeros: it
% means that there is no match for graphs
[M, dont_match_1] = first_update_M(M, A_pattern, A_molecule );
% refine matrix M for each non_zero entry
if ~dont_match_1
    % find candidates
    [M, dont_match_2] = second_update_M(M, A_pattern, A_molecule );
    if ~dont_match_2
        indicator = 1; % if the algorithm didn't reach this part, it means that there is no isomorphism
    end
end


end

