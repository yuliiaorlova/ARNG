%% prepare all unwanted species
% load('08.02_all_reactions_NEW_CL_with_RN.mat')
global A Species Molecules Reactions type weight 
A_full = A;
Species_full = Species;
Molecules_full = Molecules;
Reactions_full = Reactions; 
type_full = type;

radicals = [8, 9, 11, 19, 20, 21, 29, 39]; % all radical patterns
ids_1 = [1: length(Species)]';
di_tri = zeros(length(Species), 1);
hanging_di_tri = zeros(length(Species), 1);
unwanted_species = [];
% unwanted_species_1 = [];
% find all di and tri radicals
%%
% for i =1 : length(Species_full)
%     
%     if sum(histc(Species_full(i).ptrn, intersect(Species_full(i).ptrn, radicals)))>1
%         unwanted_species_1 = [unwanted_species_1 i];
%     end
%     
% end
%%
for i =1 : length(Species_full)
    
    if sum(diag(Species_full(i).adj))>1
        unwanted_species = [unwanted_species i];
    end
    
end

%%
di_tri(unwanted_species) = 1; 

%% 1) remove incoming and outgoing reaction to and from unwanted species
unwanted_reactions = [];
for i = 1 : length(unwanted_species)

    incoming = find(A(:,unwanted_species(i))>0);
    outgoing = find(A(unwanted_species(i), :)>0);
    unwanted_reactions = [unwanted_reactions; incoming];
    unwanted_reactions = [unwanted_reactions; outgoing'];
    clear ingoing outgoing
end
%%
unwanted_reactions = unique(unwanted_reactions);
% combine ids of unwanted nodes and their surrounding reaction in one array
unwanted_species_and_reactions = [unwanted_species'; unwanted_reactions]; 

% remove this array from the reaction network

A(unwanted_species_and_reactions, :) = [];
A(:,unwanted_species_and_reactions) = [];
ids_1(unwanted_species) = []; % keeps species IDs from the initial full list
Species(unwanted_species) = [];

weight(unwanted_species_and_reactions) = [];
type(unwanted_species_and_reactions) = [];



%% 2) check out if there are still hanging nodes here (make this in a while loop)

% initialize loop
hanging_in_species = [];
for i = 1 : length(Species)
  
    in_deg(i) = sum(A(:,i));

end
in_deg = full(in_deg);

% species that would not be reached without di and tri radicals (can be safely removed)

hanging_in_species = find(in_deg==0)'; % find all hanging nodes
%%
hanging_in_species(1) = []; % remove EL from there
count = 0;
di_tri(ids_1(hanging_in_species)) = 2;

%%
while ~isempty(hanging_in_species)
    
    count = count + 1;
    disp(count);
    
    
    hanging_species_1 = ids_1(hanging_in_species);
    % find all the outgoing reactions from the i-th generation of the hanging nodes
    
    hanging_in_reactions = [];
    hanging_in_species_and_reactions = [];
    for i = 1 : length(hanging_in_species)
        outgoing = find(A(hanging_in_species(i), :)>0);
        hanging_in_reactions = [hanging_in_reactions; outgoing'];
        clear outgoing
    end
    hanging_in_reactions = unique(hanging_in_reactions);
    hanging_in_species_and_reactions = [hanging_in_species; hanging_in_reactions];
    
    A(hanging_in_species_and_reactions, :) = [];
    A(:,hanging_in_species_and_reactions) = [];
    ids_1(hanging_in_species) = [];
    Species(hanging_in_species) = [];
    weight(hanging_in_species_and_reactions) = [];
    type(hanging_in_species_and_reactions) = [];
    
    clear in_deg
    
    % initialize loop
    for i = 1 : length(Species)
        
        in_deg(i) = sum(A(:,i));
        
    end
    in_deg = full(in_deg);
    % species that would not be reached without di and tri radicals (can be safely removed)
    hanging_in_species = [];
    hanging_in_species = find(in_deg==0)'; % find all hanging nodes
    hanging_in_species(1) = []; % remove EL from there
    if ~isempty(hanging_in_species)
        di_tri(ids_1(hanging_in_species)) = 2+count;
    end
end

%% reassign all molecules and raction lists
Molecules = [ ];
Reactions   = [ ];
for i = 1 : length(type)
    if type(i) == 1
        Molecules(end+1) = i;
    else
        Reactions(end + 1) = i;
    end
    
end

%% make sure the reaction network remains connected
A_ = A + A'; % make A undirected

[S, ~] = graphconncomp(sparse(A_), 'Directed', false);

if S == 1
    disp('A is reduced and connected')
end

clear S unwanted_species unwanted_reactions unwanted_species_and_reactions
clear di_tri hanging_di_tri hanging_in_species ids_1 in_deg A_

if length(Species) == length(Species_full)
        
        clear A_full Species_full Molecules_full Reactions_full type_full
end
