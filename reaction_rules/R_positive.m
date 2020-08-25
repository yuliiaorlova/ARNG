function candidate_product = R_positive( reactant, product, candidate, size_diff )
% reactant is bigger than product
global Species;

idx_prod = 0;
add_matrix = zeros(size_diff, size_diff);
product_1 = blkdiag(product.adj, add_matrix);
crosslinks_1 = product.crsl;
for i = 1 : size_diff
    crosslinks_1{end+1} = [];
end
candidate_product = [];


for k = 1 : length(candidate.patterns)
    
    % applying transformation to the molecule
    A_temp = Species(candidate.index).adj;
    labels_temp = Species(candidate.index).lbl;
    crosslinks_temp = Species(candidate.index).crsl;
    diag1 = diag(A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)));
    diag2 = diag(reactant.adj(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)));
    
    A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)) = product_1(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)) + diag((diag1 - diag2));
    % should include the same for crosslinks
    for i = 1 : length(Species(candidate.index).match{candidate.patterns(k)}(:, 2))
        if ~isempty(crosslinks_1{Species(candidate.index).match{candidate.patterns(k)}(i,1)}) % add only if there is something to add
            if crosslinks_1{Species(candidate.index).match{candidate.patterns(k)}(i,1)}> 0 % if positive - add crosslink
                crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} = [crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} crosslinks_1{Species(candidate.index).match{candidate.patterns(k)}(i,1)}];
            else % if negative - make corresponding crosslink vanish
                [~, b] = ismember(abs(crosslinks_1{Species(candidate.index).match{candidate.patterns(k)}(i,1)}), crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)});
                crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)}(b) = [];
                if isempty(crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)})
                    crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} = [];
                end
            end
        end
    end
    z = find(sum(A_temp(:,:)) == 0);
    A_temp(z, :) = [];
    A_temp(:, z) = [];
    labels_temp(z) = [];
    crosslinks_temp(z) = [];
    
    idx_prod = add_molecule(A_temp, labels_temp, crosslinks_temp);
    candidate_product = [candidate_product; idx_prod];
    
    clear A_temp;
    clear labels_temp;
    clear ind;
    clear matched_temp;
end


end


