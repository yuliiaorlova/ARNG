
% estimate how many symmetric reactions are in R22_candidates and R23_candidates

symm_rules_22 = [];
symm_rules_23 = [];

n_rections_R11 = 0;
for i = 1 : length(R11_candidates) % no symmetric reactions
    n_rections_R11 = n_rections_R11 + length(R11_candidates(i).reactant1);
end

n_rections_R12 = 0;
for i = 1 : length(R12_candidates) % no symmetric reactions
    n_rections_R12 = n_rections_R12 + length(R12_candidates(i).reactant1);
end

n_rections_R21 = 0;
for i = 1 : length(R21_candidates) % no symmetric reactions
    n_rections_R21 = n_rections_R21 + length(R21_candidates(i).reactant1)* length(R21_candidates(i).reactant2);
end

n_rections_R22_reduced = 0;
for i = 1 : length(R22_candidates)
%     if ismember(i, symm_rules_22)
%         n_rections_R22_reduced = n_rections_R22_reduced + (size(combnk(1:length(R22_candidates(i).reactant1),2), 1) + length(R22_candidates(i).reactant2));
%     else
        n_rections_R22_reduced = n_rections_R22_reduced + length(R22_candidates(i).reactant1)* length(R22_candidates(i).reactant2);
%    end
end
n_rections_R23_reduced = 0;
if ~isempty(R23(1).reactant1)
for i = 1 : length(R23_candidates)
        n_rections_R23_reduced = n_rections_R23_reduced + length(R23_candidates(i).reactant1)* length(R23_candidates(i).reactant2);
end
end
% total size of the bipartite reaction network
NNS = length(Species) + n_rections_R11 + n_rections_R12 + n_rections_R21 + n_rections_R22_reduced + n_rections_R23_reduced;
%%

global  Reactions Molecules weight type y0  A
NS = length(Species);
weight = zeros(NNS, 1);
type = zeros(NNS, 1);

type(1:NS) = 1;
%%
% initalize concentrations
c_ba_0= 3.7106; % 3.89 bis-allyl per molecule (ref. Mallegol)
c_O2_0=2.18e-3;
c_ini_0=1e-2; % was 1e-5

weight(1) = c_ini_0;
weight(2) = c_O2_0;
weight(3) = c_ba_0;

%%
% generate indices for the nodes in the reaction network
LA = length(Species);

list_i_all = [];
list_j_all = [];
global pointer
pointer = 0;

% first order rules: R1 -> P1
for i =1  : length(R11_candidates)
    if ~isempty(R11_candidates(i).reactant1)&& ~isempty(R11_candidates(i).product1)
        n_r_temp = length(R11_candidates(i).reactant1);
        list = zeros(3, n_r_temp);
        list(1, :) = repmat( (R11_candidates(i).reactant1(:)), [1, 1] );
        list(2, :) = repmat( (R11_candidates(i).product1(:)), [1, 1] );
        list(3, :) = (LA+1:LA+n_r_temp);
        
        weight(LA+1 : LA + n_r_temp) = R11(i).k;
        type(LA+1 : LA + n_r_temp) = R11(i).type;
        LA = LA+n_r_temp;
        
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(3,:));
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(3,:), list(2,:));
        
        clear n_r_temp idx_temp_1 idx_temp_2 list
    end
end


%%
% first order rules: R1 -> P1 + P2
for i =1  : length(R12_candidates)
    if ~isempty(R12_candidates(i).reactant1)&& ~isempty(R12_candidates(i).product1) && ~isempty(R12_candidates(i).product2)
        n_r_temp = length(R12_candidates(i).reactant1);
        list = zeros(4, n_r_temp);
        list(1, :) = repmat( (R12_candidates(i).reactant1(:)), [1, 1] );
        list(2, :) = repmat( (R12_candidates(i).product1(:)), [1, 1] );
        list(3, :) = repmat( (R12_candidates(i).product2(:)), [1, 1] );
        list(4, :) = (LA+1:LA+n_r_temp);
        
        weight(LA+1 : LA + n_r_temp) = R12(i).k;
        type(LA+1 : LA + n_r_temp) = R12(i).type;
        LA = LA+n_r_temp;
        
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(4,:));
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(4,:), list(2,:));
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(4,:), list(3,:));
        
        clear n_r_temp idx_temp_1 idx_temp_2 idx_temp_3 list
    end
end


%%
% second order rules: R1 + R2 -> P1

for i =1  : length(R21_candidates)
    if ~isempty(R21_candidates(i).reactant1)&& ~isempty(R21_candidates(i).reactant2)&& ~isempty(R21_candidates(i).product1)
        n_r_temp = length(R21_candidates(i).reactant1)*length(R21_candidates(i).reactant2);
        list = zeros(4, n_r_temp);
        list(1, :) = repmat( (R21_candidates(i).reactant1(:)), [length(R21_candidates(i).reactant2), 1] );
        y_temp = repmat(R21_candidates(i).reactant2(:), 1, length(R21_candidates(i).reactant1));
        y_temp = y_temp';
        list(2, :) = y_temp(:)';
        clear y_temp
        list(3, :) = repmat( (R21_candidates(i).product1(:)), [length(R21_candidates(i).reactant2), 1] );
        list(4, :) = (LA+1:LA+n_r_temp);
        

        weight(LA+1 : LA + n_r_temp) = R21(i).k;
        type(LA+1 : LA + n_r_temp) = R21(i).type;
        LA = LA+n_r_temp;
        
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(4,:));
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(2,:), list(4,:));
        [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(4,:), list(3,:));
        

        clear n_r_temp  idx_temp_1 idx_temp_2 idx_temp_3 list
    end
end

%%
% second order rules: R1 + R2 -> P1 + P2
for i = 1 : length(R22_candidates)
    if ~isempty(R22_candidates(i).reactant1)&& ~isempty(R22_candidates(i).reactant2)&& ~isempty(R22_candidates(i).product1)&& ~isempty(R22_candidates(i).product2)
        if ~ismember(i, symm_rules_22)
            n_r_temp = length(R22_candidates(i).reactant1)*length(R22_candidates(i).reactant2);
            list = zeros(5, n_r_temp);
            list(1, :) = repmat( (R22_candidates(i).reactant1(:)), [length(R22_candidates(i).reactant2), 1] );
            list(3, :) = repmat( (R22_candidates(i).product1(:)), [length(R22_candidates(i).product2), 1] );
            
            y_temp = repmat(R22_candidates(i).reactant2(:), 1, length(R22_candidates(i).reactant1));
            y_temp = y_temp';
            list(2, :) = y_temp(:)';
            clear y_temp
            
            y_temp = repmat(R22_candidates(i).product2(:), 1, length(R22_candidates(i).product1));
            y_temp = y_temp';
            list(4, :) = y_temp(:)';
            clear y_temp
            
            list(5, :) = (LA+1:LA+n_r_temp);
            
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(5,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(2,:), list(5,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(5,:), list(3,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(5,:), list(4,:));
            
            weight(LA+1 : LA + n_r_temp) = R22(i).k;
            type(LA+1 : LA + n_r_temp) = R22(i).type;
            LA = LA+n_r_temp;
            
            
        else
            n_r_temp = size(combnk(1:length(R22_candidates(i).reactant1),2), 1) + length(R22_candidates(i).reactant2);
            list = zeros(5, n_r_temp);
            y_temp = combnk(R22_candidates(i).reactant1(:), 2);
            y_temp_1 = zeros(length(R22_candidates(i).reactant1), 2);
            y_temp_1(:, 1) = R22_candidates(i).reactant1;
            y_temp_1(:, 2) = R22_candidates(i).reactant1;
            y_temp_2 = [y_temp; y_temp_1];
            list(1, :) =  y_temp_2(:, 2); % reactant1
            list(2, :) =  y_temp_2(:, 1); % reactant2
            clear y_temp y_temp_1 y_temp_2
            
            y_temp = combnk(R22_candidates(i).product1(:), 2);
            y_temp_1 = zeros(length(R22_candidates(i).product1), 2);
            y_temp_1(:, 1) = R22_candidates(i).product1;
            y_temp_1(:, 2) = R22_candidates(i).product1;
            y_temp_2 = [y_temp; y_temp_1];
            list(3, :) =  y_temp_2(:, 2); % reactant1
            list(4, :) =  y_temp_2(:, 1); % reactant2
            clear y_temp y_temp_1 y_temp_2
            
            
            list(5, :) = (LA+1:LA+n_r_temp);
            
            
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(5,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(2,:), list(5,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(5,:), list(3,:));
            [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(5,:), list(4,:));
            
            
            weight(LA+1 : LA + n_r_temp) = R22(i).k;
            type(LA+1 : LA + n_r_temp) = R22(i).type;
            LA = LA+n_r_temp;
            
        end
        clear n_r_temp list idx_temp_1 idx_temp_2 idx_temp_3 idx_temp_4
    end
end
%%
% second order rules: R1 + R2 -> P1 + P2 + P3
if ~isempty(R23(1).reactant1)
    for i = 1 : length(R23_candidates)
        if ~isempty(R23_candidates(i).reactant1)&& ~isempty(R23_candidates(i).reactant2)&& ~isempty(R23_candidates(i).product1)&& ~isempty(R23_candidates(i).product2)
            if ~ismember(i, symm_rules_23)
                n_r_temp = length(R23_candidates(i).reactant1)*length(R23_candidates(i).reactant2);
                list = zeros(6, n_r_temp);
                list(1, :) = repmat( (R23_candidates(i).reactant1(:)), [length(R23_candidates(i).reactant2), 1] );
                list(3, :) = repmat( (R23_candidates(i).product1(:)), [length(R23_candidates(i).product2), 1] );
                
                y_temp = repmat(R23_candidates(i).reactant2(:), 1, length(R23_candidates(i).reactant1));
                y_temp = y_temp';
                list(2, :) = y_temp(:)';
                clear y_temp
                
                y_temp = repmat(R23_candidates(i).product2(:), 1, length(R23_candidates(i).product1));
                y_temp = y_temp';
                list(4, :) = y_temp(:)';
                clear y_temp
                list(5, :) = R23_candidates(i).product3;
                
                list(6, :) = (LA+1:LA+n_r_temp);
                
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(6,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(2,:), list(6,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(3,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(4,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(5,:));
                
                weight(LA+1 : LA + n_r_temp) = R23(i).k;
                type(LA+1 : LA + n_r_temp) = R23(i).type;
                LA = LA+n_r_temp;
            else
                n_r_temp = size(combnk(1:length(R23_candidates(i).reactant1),2),1) + length(R23_candidates(i).reactant2);
                list = zeros(6, n_r_temp);
                
                y_temp = combnk(R23_candidates(i).reactant1(:), 2);
                y_temp_1 = zeros(length(R23_candidates(i).reactant1), 2);
                y_temp_1(:, 1) = R23_candidates(i).reactant1;
                y_temp_1(:, 2) = R23_candidates(i).reactant1;
                y_temp_2 = [y_temp; y_temp_1];
                list(1, :) =  y_temp_2(:, 2); % reactant1
                list(2, :) =  y_temp_2(:, 1); % reactant2
                clear y_temp y_temp_1 y_temp_2
                
                y_temp = combnk(R23_candidates(i).product1(:), 2);
                y_temp_1 = zeros(length(R23_candidates(i).product1), 2);
                y_temp_1(:, 1) = R23_candidates(i).product1;
                y_temp_1(:, 2) = R23_candidates(i).product1;
                y_temp_2 = [y_temp; y_temp_1];
                list(3, :) =  y_temp_2(:, 2); % reactant1
                list(4, :) =  y_temp_2(:, 1); % reactant2
                clear y_temp y_temp_1 y_temp_2
                
                list(5, :) = R23_candidates(i).product3;
                
                list(6, :) = (LA+1:LA+n_r_temp);
                
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(1,:), list(6,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(2,:), list(6,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(3,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(4,:));
                [list_i_all,list_j_all] = make_list(list_i_all,list_j_all,list(6,:), list(5,:));
                
                
                weight(LA+1 : LA + n_r_temp) = R23(i).k;
                type(LA+1 : LA + n_r_temp) = R23(i).type;
                LA = LA+n_r_temp;
            end
            clear n_r_temp list idx_temp_1 idx_temp_2 idx_temp_3 idx_temp_4 idx_temp_5
        end
    end
end

%%
Molecules = [ ];
Reactions   = [ ];
for i = 1 : length(type)
    if type(i) == 1
        Molecules(end+1) = i;
    else
        Reactions(end + 1) = i;
    end
    
end

% initialize y0 - initial concentrations needed to solve ODEs
y0 = zeros(length(Molecules),1);
y0 = weight(find(type==1));

% assign values to the adjacency matrix of the reaction network
A=sparse(list_i_all,list_j_all,1,NNS,NNS);

clear list_i_all list_j_all n_rections* pointer symm_rules_* LA

%% remove di and tri radicals from the reaction network
run('remove_unwanted_species.m'); 

%% save output
save('output/RN_EL.mat')

%% run next step
run('main_Rate_Parameters.m')
