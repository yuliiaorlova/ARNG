% remove di and tri radicals
A_full = A;
radicals = [8, 9, 11, 17, 22, 23, 24, 32, 42]; % all radical patterns
ids_1 = [1: length(Species)]';
di_tri_radicals = [];
%% finad all di and tri radicals
for i =1 : length(Species)
    
    if sum(histc(Species(i).ptrn, intersect(Species(i).ptrn, radicals)))>1
        di_tri_radicals = [di_tri_radicals i];
    end
    
end

%% remove incoming and outgoing reaction to and from di and tri radicals
di_tri_radicals_reactions = [];
for i = 1 : length(di_tri_radicals)

    ingoing = find(A(:,di_tri_radicals(i))>0);
    outgoing = find(A(di_tri_radicals(i), :)>0);
    di_tri_radicals_reactions = [di_tri_radicals_reactions; ingoing];
    di_tri_radicals_reactions = [di_tri_radicals_reactions; outgoing'];
    clear ingoing outgoing
end
di_tri_radicals_reactions = unique(di_tri_radicals_reactions);
% combine ids of di- tri- radicals nodes and their surrounding reaction in one array
di_tri_radicals_neighbors = [di_tri_radicals'; di_tri_radicals_reactions]; % stopped here

% remove this array from the reaction network

A(di_tri_radicals_neighbors, :) = [];
A(:,di_tri_radicals_neighbors) = [];
ids_1(di_tri_radicals) = [];
% Species(di_tri_radicals) = [];
% 
% weight(di_tri_radicals_neighbors) = [];
% type(di_tri_radicals_neighbors) = [];

% Molecules = [ ];
% Reactions   = [ ];
% for i = 1 : length(type)
%     if type(i) == 1
%         Molecules(end+1) = i;
%     else
%         Reactions(end + 1) = i;
%     end
%     
% end

%% check out if there are still hanging nodes here

for i =1 : length(ids_1)
    out_deg(i) = sum(A(i, :));
    in_deg(i) = sum(A(:,i));

end
out_deg = full(out_deg);
in_deg = full(in_deg);
degree = [out_deg' in_deg'];
%% species that would not be reached without di and tri radicals (can be safely removed)

hanging_in_nodes = find(in_deg==0)'; % find all hanging nodes
hanging_in_nodes(1) = []; % remove EL from there
hanging_species_1 = ids_1(hanging_in_nodes);
%% find all the outgoing reactions from the 1st generation of the hanging nodes 

hanging_in_nodes_neighbors = [];
for i = 1 : length(hanging_in_nodes)
    outgoing = find(A(hanging_in_nodes(i), :)>0);
    hanging_in_nodes_neighbors = [hanging_in_nodes_neighbors; outgoing'];
    clear outgoing
end
hanging_in_nodes_neighbors = unique(hanging_in_nodes_neighbors);
hanging_in_nodes_neighbors = [hanging_in_nodes; hanging_in_nodes_neighbors];
%

A(hanging_in_nodes_neighbors, :) = [];
A(:,hanging_in_nodes_neighbors) = [];

ids_1(hanging_in_nodes) = [];

clear out_deg in_deg degree
%% find 2nd generation of hanging nodes: make sure their removal makes sense


for i =1 : length(ids_1)
    out_deg(i) = sum(A(i, :));
    in_deg(i) = sum(A(:,i));

end

out_deg = full(out_deg);
in_deg = full(in_deg);
degree = [out_deg' in_deg'];
%% species that would not be reached without di and tri radicals (can be safely removed)

hanging_in_nodes = find(in_deg==0)'; % find all hanging nodes
hanging_in_nodes(1) = []; % remove EL from there
hanging_species_2 = ids_1(hanging_in_nodes);
%% find all the outgoing reactions from the 1st generation of the hanging nodes 

hanging_in_nodes_neighbors = [];
for i = 1 : length(hanging_in_nodes)
    outgoing = find(A(hanging_in_nodes(i), :)>0);
    hanging_in_nodes_neighbors = [hanging_in_nodes_neighbors; outgoing'];
    clear outgoing
end
hanging_in_nodes_neighbors = unique(hanging_in_nodes_neighbors);
hanging_in_nodes_neighbors = [hanging_in_nodes; hanging_in_nodes_neighbors];
%%

A(hanging_in_nodes_neighbors, :) = [];
A(:,hanging_in_nodes_neighbors) = [];

ids_1(hanging_in_nodes) = [];
clear out_deg in_deg degree

%% this was the last generation
% next time should make this provess in a loop (while hanging_in_nodes is empty)
hanging_species = [hanging_species_1; hanging_species_2; hanging_species_3];

% Species(hanging_in_nodes) = [];
% 
% weight(hanging_in_nodes_neighbors) = [];
% type(hanging_in_nodes_neighbors) = [];
% 
% Molecules = [ ];
% Reactions   = [ ];
% for i = 1 : length(type)
%     if type(i) == 1
%         Molecules(end+1) = i;
%     else
%         Reactions(end + 1) = i;
%     end
%     
% end

