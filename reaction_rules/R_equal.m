function candidate_product =  R_equal( reactant, product, candidate )
global Species;

candidate_product = [];
idx_prod = 0;

for k = 1 : length(candidate.patterns)
    
    labels_temp = Species(candidate.index).lbl;
    A_temp = Species(candidate.index).adj;
    crosslinks_temp = Species(candidate.index).crsl;
    
    % manipulations with crosslink-pattern overlap:
    diag1 = diag(A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)));
    diag2 = diag(reactant.adj(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)));
    
    A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)) = product.adj(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)) + diag((diag1 - diag2));
    labels_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2)) = product.lbl(Species(candidate.index).match{candidate.patterns(k)}(:,1));
    for i = 1 : length(Species(candidate.index).match{candidate.patterns(k)}(:, 2))
        if ~isempty(product.crsl{Species(candidate.index).match{candidate.patterns(k)}(i,1)}) % add only if there is something to add
            if product.crsl{Species(candidate.index).match{candidate.patterns(k)}(i,1)}> 0 % if positive - add crosslink
                crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} = [crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} product.crsl{Species(candidate.index).match{candidate.patterns(k)}(i,1)}];
            else % if negative - make corresponding crosslink vanish
                [~, b] = ismember(abs(product.crsl{Species(candidate.index).match{candidate.patterns(k)}(i,1)}), crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)});
                crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)}(b) = [];
                 if isempty(crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)})
                    crosslinks_temp{Species(candidate.index).match{candidate.patterns(k)}(i, 2)} = [];
                end
            end
        end
    end
    
    idx_prod  = add_molecule(A_temp, labels_temp, crosslinks_temp);
    candidate_product = [candidate_product; idx_prod];
    
    clear A_temp diag1 diag2;
    clear labels_temp;
end


end

