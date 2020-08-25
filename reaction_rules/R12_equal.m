function candidate_product = R12_equal( R12, candidate )
global Species;

candidate_product = [];

for k = 1 : length(candidate.patterns)
    
    product = R12.product1;
    reactant = R12.reactant1;
    labels_temp = Species(candidate.index).lbl;
    A_temp = Species(candidate.index).adj;
    crosslinks_temp = Species(candidate.index).crsl;
    
    diag1 = diag(A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)));
    diag2 = diag(reactant.adj(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)));
    
    A_temp(Species(candidate.index).match{candidate.patterns(k)}(:, 2), Species(candidate.index).match{candidate.patterns(k)}(:, 2)) = product.adj(Species(candidate.index).match{candidate.patterns(k)}(:,1), Species(candidate.index).match{candidate.patterns(k)}(:,1)) + diag((diag1 - diag2));
    
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
    
    [prod] = find_conn_comp(A_temp);
    
    A_temp_prod1 = A_temp(prod{1}, prod{1});
    labels_temp_prod1 = labels_temp(prod{1});
    crosslinks_temp_prod1 = crosslinks_temp(prod{1});
    
    A_temp_prod2 = A_temp(prod{2}, prod{2});
    labels_temp_prod2 = labels_temp(prod{2});
    crosslinks_temp_prod2 = crosslinks_temp(prod{2});
    
    idx_prod1  = add_molecule(A_temp_prod1, labels_temp_prod1, crosslinks_temp_prod1 );
    idx_prod1_ = idx_prod1;
    
    if isequal(labels_temp_prod2, {'Co'})
       idx_prod2_ = 1; 
    else
    idx_prod2 = add_molecule(A_temp_prod2, labels_temp_prod2, crosslinks_temp_prod2 );
    idx_prod2_ = idx_prod2;
    end
    
    candidate_product = [candidate_product; idx_prod1_ idx_prod2_];
    
    clear A_temp;
    clear A_temp_prod1;
    clear A_temp_prod2;
    clear labels_temp;
    clear labels_temp_prod1;
    clear labels_temp_prod2;
    clear idx_prod1;
    clear idx_prod2;
    
end

end

