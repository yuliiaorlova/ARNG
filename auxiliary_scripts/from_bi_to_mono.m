%% full projected reaction network

% take care of small heavily connected species
patterns_nr = [1,2,3,4,5,6,38];
nr_full = [];
for i = 1 : length(Species_full)
   
    if ~isempty(intersect(Species_full(i).ptrn, patterns_nr))
        nr_full = [nr_full; i]; 
    end
    
end

nr_small = [];
for i = 1 : length(Species)

    if ~isempty(intersect(Species(i).ptrn, patterns_nr))
        nr_small = [nr_small; i];
    end

end

%%
A_mono = sparse(length(Species_full), length(Species_full));

for i = 1 : length(R11_candidates)
A_mono(sub2ind(size(A_mono),R11_candidates(i).reactant1, R11_candidates(i).product1)) = R11(i).type;
end

for i = 1 : length(R12_candidates)
A_mono(sub2ind(size(A_mono),R12_candidates(i).reactant1, R12_candidates(i).product1)) = R12(i).type;
A_mono(sub2ind(size(A_mono),R12_candidates(i).reactant1, R12_candidates(i).product2)) = R12(i).type;
end

for i = 1 : length(R21_candidates)
A_mono(sub2ind(size(A_mono),R21_candidates(i).reactant1, R21_candidates(i).product1)) = R21(i).type;
end

for i = 1 : length(R22_candidates)
A_mono(sub2ind(size(A_mono),R22_candidates(i).reactant1, R22_candidates(i).product1)) = R22(i).type;
A_mono(sub2ind(size(A_mono),R22_candidates(i).reactant2, R22_candidates(i).product2)) = R22(i).type;

end

for i = 1 : length(R23_candidates)
A_mono(sub2ind(size(A_mono),R23_candidates(i).reactant1, R23_candidates(i).product1)) = R23(i).type;
A_mono(sub2ind(size(A_mono),R23_candidates(i).reactant2, R23_candidates(i).product2)) = R23(i).type;
end

A_mono_1 = full(A_mono);
%A_mono_1 = A_mono_1 + A_mono_1';
A_mono_1(nr_full, :) = [];
A_mono_1(:, nr_full) = [];
A_mono_1(find(A_mono_1>0)) = 1;
Mol_1 = Molecules_full;
Mol_1(nr_full) = [];
di_tri(nr_full) = [];
%%
% dlmwrite('A_monopartite_full_31.10.txt',A_mono_1','delimiter',' ');
% graphtogml('A_monopartite_full_31.10.gml', A_mono_1', Mol_1);

%% try to make reduced reaction network

A_mono_small = full(A_mono);
%A_mono_small = A_mono_small + A_mono_small';
species_to_remove = setdiff(1:length(Species_full), ids_1)';
%%
to_remove = species_to_remove;
A_mono_small(to_remove,:) = [];
A_mono_small(:,to_remove) = [];
A_mono_small(find(A_mono_small>0)) = 1;

A_mono_small(nr_small,:) = [];
A_mono_small(:, nr_small) = [];
Mol_small = Molecules_full;
Mol_small(to_remove) = [];
Mol_small(nr_small) = [];
%%
Species_small = Species;
Species_small(nr_small) = [];
%%
% crsl_file = zeros(length(Species_small), 2);
% for i = 1 : length(Species_small)
%     crsl_file(i, 2) = length(find(~cellfun(@isempty,Species_small(i).crsl)));
% end
% 
% crsl_file(:, 1) = [1:length(crsl_file)];
% dlmwrite('crsl_A_monopartite_small_31_10.txt',crsl_file,'delimiter',' ');
%%
 dlmwrite('A_monopartite_full_14_08.txt',A_mono_small','delimiter',' ');
 graphtogml('A_monopartite_full_14_98.gml', A_mono_small', Mol_small);% di_tri_small);
%% clear out_deg in_deg degree
% 
% for i =1 : size(A_mono_small,1)
%     out_deg(i) = sum(A_mono_small(i, :));
%     in_deg(i) = sum(A_mono_small(:,i));
%     
% end
% 
% out_deg = full(out_deg);
% in_deg = full(in_deg);
% degree = [out_deg' in_deg'];
% hanging_in_nodes = find(in_deg==0)';
%%

comp_mat = find_conn_comp(A_mono_small);

% dlmwrite('A_monopartite_small.txt',A_mono_small,'delimiter',' ');
% graphtogml('A_mono_small.gml', tril(A_mono_small), Mol_small);
% 
% %% find a spanning tree of A_mono_small
% 
% components = find_conn_comp(A_mono_small);
% 
% %%
% A_mono_small_conn = A_mono_small(components{1}, components{1});
% %%
% Mol_small_conn = Mol_small(components{1});
% %clear components
% %%
% sp_tree = zeros(size(A_mono_small_conn, 1), size(A_mono_small_conn, 1));
% %%
% for i = 2 : size(A_mono_small_conn, 1)
% 
%     [dist,P]=dijkstra(A_mono_small_conn,1,i);
%     edges = zeros(length(P)-1, 2);
%     edges(:,1) = P(1:end-1);
%     edges(:,2) = P(2: end);
%     sp_tree(sub2ind(size(A_mono_small_conn), edges(:,1), edges(:,2))) = 1;
%     clear dist P
% end
% %%
% graphtogml('sp_tree.gml', sp_tree', Mol_small_conn);
%% 

% for i = 1 : length(R23_candidates)
% if ~isempty(intersect(R23_candidates(i).product1, 780))
% disp(i)
% end
% end
% for i = 1 : length(R23_candidates)
% if ~isempty(intersect(R23_candidates(i).product2, 1096))
% disp(i)
% end
% end

% for i = 1 : length(R23_candidates)
% if ~isempty(intersect(R23_candidates(i).product3, 780))
% disp(i)
% end
% end
%%
% for i =1 : length(Species_full)
%     [ matched, ptrn ] = match_ullmann_pattern( Species_full(i).adj, Species_full(i).lbl, Species_full(i).crsl ); % find patterns in newly created molecule
%     Species_full(i).match1 = matched;
%     Species_full(i).ptrn1 = ptrn;
% end
% %%
% for i =1 : length(Species_full)
%     if ~isequal(Species_full(i).ptrn, Species_full(i).ptrn1)
%         disp(i)
%         graph2mol(Species_full(i).adj, Species_full(i).lbl, i)
%     end
% end
