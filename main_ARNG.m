addpath('auxiliary_scripts/')
addpath('functions/')
addpath('input_molecular_graphs/');
addpath('input_patterns/');
addpath('patterns/')
addpath('reaction_rules/')
addpath('ullmann_isomorphism_test_functions');
warning('off','all')

% initialize patterns

global number_of_reaction_families;
number_of_reaction_families = 28;

global pattern;
pattern = pattern_struct();

global pattern_non_reactive;
pattern_non_reactive = pattern_non_reactive_struct();

% initialize reaction rules
% R11 - first order reaction rules with one reactant and one product
% R12 - first order reaction rules with one reactant and two products
% R21 - second order reaction rules with two reactants and one product
% R22 - second order reaction rules with two reactant and two products
% R23 - second order reaction rules with two reactant and three products

global R11 R12 R21 R22 R23
%              TAKES                   GIVES
R11 = struct('reactant1', [], 'product1', [], 'k',[], 'type', []); % first order rule: R1 -> P1
%              TAKES                   GIVES               GIVES
R12 = struct('reactant1', [], 'product1', [], 'k',[], 'type', []); % first order rule: R1 -> P1 + P2;
%              TAKES                   TAKES               GIVES
R21 = struct('reactant1', [], 'reactant2', [], 'product1', [], 'k',[], 'type', []); % second order rule: R1 + R2 -> P1;
%              TAKES                   TAKES               GIVES               GIVES
R22 = struct('reactant1', [], 'reactant2', [], 'product1', [],'product2', [], 'k',[], 'type', []); % second order rule: R1 + R2 -> P1 + P2;
%              TAKES                   TAKES               GIVES               GIVES 
R23 = struct('reactant1', [], 'reactant2', [], 'product1', [],'product2', [], 'product3', [], 'k',[], 'type', []); % second order rule: R1 + R2 -> P1 + P2;

rules();

R11_candidates(1:length(R11)) = struct('reactant1', [], 'product1', []); % first order rule: R1 -> P1
R12_candidates(1:length(R12)) = struct('reactant1', [], 'product1', [], 'product2', []); % first order rule: R1 -> P1 + P2;
R21_candidates(1:length(R21)) = struct('reactant1', [], 'reactant2', [], 'product1', []); % second order rule: R1 + R2 -> P1;
R22_candidates(1:length(R22)) = struct('reactant1', [], 'reactant2', [], 'product1', [],'product2', []); % second order rule: R1 + R2 -> P1 + P2;
R23_candidates(1:length(R23)) = struct('reactant1', [], 'reactant2', [], 'product1', [],'product2', [], 'product3', []); % second order rule: R1 + R2 -> P1 + P2;

% initialize species
global Species;
Species = struct('adj', [], 'lbl', [], 'match', [], 'ptrn', [], 'stuff', [], 'crsl', []);

% species which are not defined explicitly: Initiator, Initiator radical
% and Metal drier (n) and (n+1) - two states of oxidation

A_molecule = importdata('adjacency_EL.txt'); % ethyl linoleate
labels_molecule = importdata('atoms_EL.txt');


Species(1).adj = 1;
Species(1).lbl = 'Co2';
Species(1).match = [];
Species(1).ptrn = 1;
Species(1).stuff = 1;
Species(1).crsl = [];
 
% Oxygen
Species(2).adj = 1;
Species(2).lbl = 'O2';
Species(2).match = [];
Species(2).ptrn = 4;
Species(2).stuff = 1;
Species(2).crsl = [];

% EL molecules
Species(3).adj = A_molecule;
Species(3).lbl = labels_molecule;
Species(3).crsl = cell(length(Species(3).lbl),1);
stuff = get_stuff(Species(3).adj, Species(3).lbl);
[matched_temp, pattern_temp] = match_ullmann_pattern( A_molecule, labels_molecule,  Species(3).crsl);
Species(3).match = matched_temp;
Species(3).ptrn = pattern_temp;
Species(3).stuff = stuff;
Species(3).crsl = cell(length(Species(3).lbl),1);

% H2O
Species(4).adj = 1;
Species(4).lbl = 'H2O';
Species(4).match = [];
Species(4).ptrn = 6;
Species(4).stuff = 1;
Species(4).crsl = [];

clear matched_temp A_molecule labels_molecule stuff
clear pattern_temp

warning('off','all')

global worth_updating;
worth_updating = 1;
count = 0;
previous_cycle(1) = 0;
%%
% automatic application of reaction rules to the input and newly generated
% molecules
tic
while worth_updating
    
    count = count + 1;
    disp('cycle #')
    disp(count)
    worth_updating = 0;
    n_species = length(Species);
    previous_cycle(count+1) = length(Species);
    str = sprintf('from: %d to: %d', previous_cycle(count), n_species);
    disp(str)
    
    % first order rules: R1 -> P1
    for n = (previous_cycle(count)+1) : n_species
        if ~isempty(R11(1).reactant1)
            for i = 1 : length(R11)
                % found candidate
                patterns_number = find(Species(n).ptrn == R11(i).reactant1.id);
                product11 = [];
                if ~isempty(patterns_number)
                    candidate.index = n;
                    candidate.patterns = patterns_number;
                    if length(R11(i).reactant1.adj)<2
                        product11 = R_small( R11(i).product1, candidate );
                    else
                        
                        size_diff = length(R11(i).reactant1.adj) - length(R11(i).product1.adj);
                        
                        if size_diff>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                            product11 = R_positive( R11(i).reactant1, R11(i).product1, candidate, size_diff );
                            
                        else if size_diff<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                product11 = R_negative( R11(i).reactant1, R11(i).product1, candidate, size_diff );
                                
                                
                            else if size_diff==0 % reactant has the same number of atoms as product
                                    product11 = R_equal(  R11(i).reactant1, R11(i).product1, candidate);
                                    
                                end
                            end
                        end
                    end
                   for k = 1 : length(product11) 
                   R11_candidates(i).reactant1 = [R11_candidates(i).reactant1; n];
                   R11_candidates(i).product1 = [R11_candidates(i).product1; product11(k)]; 
                   end
                end
                clear patterns_numebr size_diff candidate product11
            end
            
        end

        % first order rules:  R1 -> P1 + P2
        
        if ~isempty(R12(1).reactant1)
            for i = 1 : length(R12)
                % reactant and product have to be the same length
                % this type of reaction makes change only in connectivity
                % such as the adjacency matirx of the reactant can be split into
                % two connected components
                patterns_number = find(Species(n).ptrn == R12(i).reactant1.id);
                product12 = [];
                if ~isempty(patterns_number)
                    candidate.index = n;
                    candidate.patterns = patterns_number;                    
                    product12 = R12_equal( R12(i), candidate);
                    
                    for k = 1 : size(product12,1)
                        R12_candidates(i).reactant1 = [R12_candidates(i).reactant1; n];
                        R12_candidates(i).product1 = [R12_candidates(i).product1; product12(k, 1)];
                        R12_candidates(i).product2 = [R12_candidates(i).product2; product12(k, 2)];
                    end
                end
                clear  patterns_number candidate product12
            end
        end
        
        % second order rules:  R1 + R2 -> P1
        
        if ~isempty(R21(1).reactant1) % only reactant_1 undergoes transformation, reactant_2 is always a small species
            for i = 1 : length(R21)
                product21 = [];
                patterns_number = find(Species(n).ptrn == R21(i).reactant1.id);
                if ~isempty(patterns_number)
                    candidate.index = n;
                    candidate.patterns = patterns_number;
                    
                    if length(R21(i).reactant1.adj) <2
                        product21 = R_small( R21(i).product1, candidate);
                    else
                        size_diff = length(R21(i).reactant1.adj) - length(R21(i).product1.adj);
                        if size_diff>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                            product21 = R_positive( R21(i).reactant1, R21(i).product1, candidate, size_diff);
                        else if size_diff<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                product21 = R_negative(  R21(i).reactant1, R21(i).product1, candidate, size_diff );
                            else if size_diff==0 % reactant has the same number of atoms as product
                                    product21 = R_equal(  R21(i).reactant1, R21(i).product1, candidate );
                                end
                            end
                        end
                    end
                    for k = 1 : length(product21)
                        R21_candidates(i).reactant1 = [R21_candidates(i).reactant1; n];
                        R21_candidates(i).product1 = [R21_candidates(i).product1; product21(k)];                   
                    end
                end
                reactant_2 = find(Species(n).ptrn == R21(i).reactant2.id);
                if ~isempty(reactant_2)
                    R21_candidates(i).reactant2 = [R21_candidates(i).reactant2; n];
                end
                clear candidate patterns_number size_diff reactant_2 product21
            end
            
        end
        
        % second order rules:  R1 + R2 -> P1 + P2
        
        if ~isempty(R22(1).reactant1)
            for i = 1 : length(R22)
                product22_1 = [];
                product22_2 = [];
                
                patterns_number1 = find(Species(n).ptrn == R22(i).reactant1.id);
                patterns_number2 = find(Species(n).ptrn == R22(i).reactant2.id);
                if ~isempty(patterns_number1)
                     candidate1.index = n;
                     candidate1.patterns = patterns_number1;
                    % firstly, take care of reactant1 and form products22_1:
                    if length(R22(i).reactant1.adj) <2
                        % consider all transformations to reactant1
                        product22_1= R_small( R22(i).product1, candidate1 );
                        
                    else if length(R22(i).reactant1.adj)>=2
                            size_diff1 = length(R22(i).reactant1.adj) - length(R22(i).product1.adj);
                            
                            if size_diff1>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                                product22_1 = R_positive( R22(i).reactant1, R22(i).product1, candidate1, size_diff1 );
                                
                            else if size_diff1<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                    product22_1 = R_negative(  R22(i).reactant1, R22(i).product1, candidate1, size_diff1 );
                                    
                                else if size_diff1 == 0
                                        product22_1 = R_equal( R22(i).reactant1, R22(i).product1, candidate1 );
                                        
                                    end
                                end
                            end
                        end
                    end
                    for k = 1 : length(product22_1)
                        R22_candidates(i).reactant1 = [R22_candidates(i).reactant1; n];
                        R22_candidates(i).product1 = [R22_candidates(i).product1; product22_1(k)];
                    end
                end
                
                
                
                
                if ~isempty(patterns_number2)
                    candidate2.index = n;
                    candidate2.patterns = patterns_number2;
                    if length(R22(i).reactant2.adj) <2
                        % consider all transformations to reactant1
                        product22_2 = R_small( R22(i).product2, candidate2 );
                        
                        
                    else if length(R22(i).reactant2.adj)>=2
                            size_diff2 = length(R22(i).reactant2.adj) - length(R22(i).product2.adj);
                            
                            if size_diff2>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                                product22_2 = R_positive( R22(i).reactant2, R22(i).product2, candidate2, size_diff2 );
                                
                            else if size_diff2<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                    product22_2 = R_negative( R22(i).reactant2, R22(i).product2, candidate2, size_diff2 );
                                    
                                else if size_diff2 == 0
                                        product22_2 = R_equal( R22(i).reactant2, R22(i).product2, candidate2 );
                                        
                                    end
                                end
                            end
                        end
                    end
                    for k = 1 : length(product22_2)
                        R22_candidates(i).reactant2 = [R22_candidates(i).reactant2; n];
                        R22_candidates(i).product2 = [R22_candidates(i).product2; product22_2(k)];
                    end
                    
                end
                

                clear candidate1 candidate2 patterns_number1 patterns_number2 size_diff1 size_diff2 product22_1 product22_2
            end
        end
        
        % second order rules:  R1 + R2 -> P1 + P2 + P3
        
        if ~isempty(R23(1).reactant1) % only reactant_1 and reactant_2 undergo transformations to product_1 and product_2, product_3 is a small species
            for i = 1 : length(R23)
                product23_1 = [];
                product23_2 = [];
                
                patterns_number1 = find(Species(n).ptrn == R23(i).reactant1.id);
                patterns_number2 = find(Species(n).ptrn == R23(i).reactant2.id);
                if ~isempty(patterns_number1)
                    candidate1.index = n;
                    candidate1.patterns = patterns_number1;
                    % firstly, take care of reactant1 and form products22_1:
                    if length(R23(i).reactant1.adj) <2
                        % consider all transformations to reactant1
                        product23_1= R_small( R23(i).product1, candidate1 );
                        
                    else if length(R23(i).reactant1.adj)>=2
                            size_diff1 = length(R23(i).reactant1.adj) - length(R23(i).product1.adj);
                            
                            if size_diff1>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                                product23_1 = R_positive( R23(i).reactant1, R23(i).product1, candidate1, size_diff1 );
                                
                            else if size_diff1<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                    product23_1 = R_negative(  R23(i).reactant1, R23(i).product1, candidate1, size_diff1 );
                                    
                                else if size_diff1 == 0
                                        product23_1 = R_equal( R23(i).reactant1, R23(i).product1, candidate1 );
                                        
                                    end
                                end
                            end
                        end
                    end
                    for k = 1 : length(product23_1)
                        R23_candidates(i).reactant1 = [R23_candidates(i).reactant1; n];
                        R23_candidates(i).product1 = [R23_candidates(i).product1; product23_1(k)];
                    end
                end
                
                if ~isempty(patterns_number2)
                    candidate2.index = n;
                    candidate2.patterns = patterns_number2;
                    if length(R23(i).reactant2.adj) <2
                        % consider all transformations to reactant1
                        product23_2 = R_small( R23(i).product2, candidate2 );
                        
                        
                    else if length(R23(i).reactant2.adj)>=2
                            size_diff2 = length(R23(i).reactant2.adj) - length(R23(i).product2.adj);
                            
                            if size_diff2>0 % reactant has more atoms than the product, we need to add n rows and columns filled with zeros
                                product23_2 = R_positive( R23(i).reactant2, R23(i).product2, candidate2, size_diff2 );
                                
                            else if size_diff2<0 % reactant has less atoms than the product, we need to add new rows and columns to the molecule
                                    product23_2 = R_negative( R23(i).reactant2, R23(i).product2, candidate2, size_diff2 );
                                    
                                else if size_diff2 == 0
                                        product23_2 = R_equal( R23(i).reactant2, R23(i).product2, candidate2 );
                                        
                                    end
                                end
                            end
                        end
                    end
                    for k = 1 : length(product23_2)
                        R23_candidates(i).reactant2 = [R23_candidates(i).reactant2; n];
                        R23_candidates(i).product2 = [R23_candidates(i).product2; product23_2(k)];
                    end
                end
                
                product_3 = find(Species(n).ptrn == R23(i).product3.id);
                if ~isempty(product_3)
                    R23_candidates(i).product3 = [R23_candidates(i).product3; n];
                end
                
                clear candidate1 candidate2 patterns_number1 patterns_number2 size_diff1 size_diff2 product23_1 product23_2 produc_3
            end
        end
        
    end
    
    disp('species #')
    disp(length(Species))
end
toc
clear count i n k previous_cycle str worth_updating product_3

%% save output
save('output/species_EL.mat');

%% run next step
run('main_Reaction_Network.m')

