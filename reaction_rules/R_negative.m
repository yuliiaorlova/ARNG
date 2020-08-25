function  candidate_product = R_negative( reactant, product, candidate,  size_diff )
global Species;
% reactant is smaller than product
candidate_product = [];
idx_prod = 0;

for k = 1 : length(candidate.patterns)
    
    % applying transformation to the molecule
    labels_temp = Species(candidate.index).lbl;
    add_matrix = zeros(abs(size_diff), abs(size_diff));
    match_temp = Species(candidate.index).match{candidate.patterns(k)};
    crosslinks_temp = Species(candidate.index).crsl;
    
    for l = 1 : abs(size_diff) % extend matrix
        match_temp(end+1, :) = [length(reactant.adj)+l  length(Species(candidate.index).adj)+l];
        labels_temp = [labels_temp; product.lbl{length(reactant.adj)+l} ];
        crosslinks_temp{end+1} = product.crsl{length(reactant.adj)+l};
    end
    A_temp = blkdiag(Species(candidate.index).adj, add_matrix);
    reactant_1 = blkdiag(reactant.adj, add_matrix);
    
    diag1 = diag(A_temp(match_temp(:, 2), match_temp(:, 2)));
    diag2 = diag(reactant_1(match_temp(:,1), match_temp(:,1)));
    
    A_temp(match_temp(:, 2), match_temp(:, 2)) = product.adj(match_temp(:,1), match_temp(:,1)) + diag((diag1 - diag2));
    for i = 1 : length(match_temp(:, 2))
        if ~isempty(product.crsl{match_temp(i, 1)}) % add only if there is something to add
            if product.crsl{match_temp(i, 1)}> 0 % if positive - add crosslink
                crosslinks_temp{match_temp(i, 2)} = [crosslinks_temp{match_temp(i, 2)} product.crsl{match_temp(i, 1)}];
            else % if negative - make corresponding crosslink vanish
                [~, b] = ismember(abs(product.crsl{match_temp(i, 1)}), crosslinks_temp{match_temp(i, 2)});
                crosslinks_temp{match_temp(i, 2)}(b) = [];
                if isempty(crosslinks_temp{match_temp(i, 2)})
                    crosslinks_temp{match_temp(i, 2)} = [];
                end
            end
        end
        %crosslinks_temp{match_temp(i, 2)} = [crosslinks_temp{match_temp(i, 2)} product.crsl{match_temp(i, 1)}];
    end
    
    idx_prod = add_molecule(A_temp, labels_temp, crosslinks_temp);
    candidate_product = [candidate_product;  idx_prod];
    
    
    clear A_temp reactant_1;
    clear labels_temp match_temp;
end


end