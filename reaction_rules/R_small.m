function candidate_product = R_small( product1, candidate )

candidate_product = [];
% match species from the list
% check if the product is already in the list
idx_prod = add_molecule(product1.adj, product1.lbl, product1.crsl); % idx_prod will always be non-empty in this case
candidate_product = [candidate_product; idx_prod];


end


