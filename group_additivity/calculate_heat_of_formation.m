% clear
% clear global
% load('Species_17.07_no_di_tri.mat')
all_groups = [];
crsl_ref = cell(3,2);
crsl_ref(:,2) = num2cell([4; 5; 3]);
crsl_ref(:,1) = cellstr(['O'; 'O'; 'C']);
run('group_library.m');

%% group is a polyvalent atom with all its adjucent ligands

for n = 1 : length(Species)
% n = 3;   
     if length(Species(n).adj)>2
         if isempty(intersect('Co', Species(n).lbl))
            crsl = [];
            crsl_table = [];
            % check for crosslinks
            
            identifier = cell(length(Species(n).lbl), 1);
           
            for i = 1 : length(Species(n).lbl)
                
                if isequal('C', Species(n).lbl{i}) && ~isempty(find(Species(n).adj(i, :)==2))
                    if isequal('O', Species(n).lbl{(find(Species(n).adj(i, :)==2))})
                        if Species(n).adj(i,i) == 1
                            identifier{i} = 'CO.';
                        else
                            identifier{i} = 'CO';
                        end
                    else
                        identifier{i} = 'Cd';
                    end
                else if Species(n).adj(i,i) == 1 && isequal('C', Species(n).lbl{i})
                        identifier{i} = 'C.';
                    else if Species(n).adj(i,i) == 1 && isequal('O', Species(n).lbl{i})
                            identifier{i} = 'O.';
                        else
                            identifier{i} = [Species(n).lbl{i}];
                        end
                    end
                end
                
            end

            
            
            if ~isempty(find(cellfun(@isempty,Species(n).crsl)==0))
                
                crsl = find(cellfun(@isempty,Species(n).crsl)==0);
                crsl_table = cell(length(crsl), 2);
                
                
                for i = 1 : length(crsl)
                    type_crsl = Species(n).crsl{crsl(i)};
                    for j = 1 : length(type_crsl)
                        
                        atoms_to_add{j} = crsl_ref{find(ismember([crsl_ref{:,2}], type_crsl(j) )),1};
                    end
                    crsl_table{i,1} = crsl(i);
                    identifier(end +1 : end+length(atoms_to_add)) = atoms_to_add';
                    crsl_table{i,2} = [crsl_table{i,2}, length(identifier)- length(atoms_to_add)+1 : length(identifier)];
                    
                end
                clear atoms_to_add
                % take care of crosslink contributions
                nbrs_crsl = {};
                for i = 1 : length(crsl)
                    type_crsl = Species(n).crsl{crsl(i)};
                    
                    for j = 1 : length(type_crsl)
                       nbrs = {};
                        if type_crsl(j) == 4
                        nbrs_crsl{end+1} = 'O-';
                        nbrs{1} = 'O';
                        nbrs{2} = identifier{crsl_table{i,1}};
                        nbrs_unique = unique(nbrs);
                        multiplicity=cellfun(@(x) sum(ismember(nbrs,x)), nbrs_unique);
                        for k = 1 : length(nbrs_unique)
                            if multiplicity(k)>1
                                nbrs_crsl{end} = strcat(nbrs_crsl{end},'(',nbrs_unique{k},')',num2str(multiplicity(k)));
                            else
                                nbrs_crsl{end} = strcat(nbrs_crsl{end},'(',nbrs_unique{k},')');
                            end
                            
                        end
                        end
                        if type_crsl(j) == 5
                        nbrs_crsl{end+1} = 'E-';
                        nbrs{1} = 'C';
                        nbrs{2} = identifier{crsl_table{i,1}};
                        nbrs_unique = unique(nbrs);
                        multiplicity=cellfun(@(x) sum(ismember(nbrs,x)), nbrs_unique);
                        for k = 1 : length(nbrs_unique)
                            if multiplicity(k)>1
                                nbrs_crsl{end} = strcat(nbrs_crsl{end},'(',nbrs_unique{k},')',num2str(multiplicity(k)));
                            else
                                nbrs_crsl{end} = strcat(nbrs_crsl{end},'(',nbrs_unique{k},')');
                            end
                            
                        end
                        end
                       
                    end
                end
            end
           
            % nbrs_list = cell(length(Species(3).lbl), 1);
            
    for i = 1 : length(Species(n).lbl)
              % i = 2; 
                % select all the atoms that will contribute to the heat of formation
                if ~isequal('CO.', identifier{i}) && ~isequal('O.', identifier{i})
                    if isequal('C', Species(n).lbl{i}) || (isequal('O', Species(n).lbl{i}) && (isempty(find(Species(n).adj(i,:)==2))))
                        nbrs_list{i} = strcat(identifier{i}, '-');
                        nbrs = neighbors(Species(n).adj, i);
                        if ~isempty(crsl_table)
                            if ~isempty(intersect(i, [crsl_table{:,1}]))
                                nbrs = [nbrs, crsl_table{find([crsl_table{:,1}]==i),2}];
                            end
                        end
                        % first get rid of Cd in neighbors if the central atom is also Cd
                        % (each Cd will have one neighbor Cd, so the intersect will never be empty)
%                       % fixed on 24.12.19
                    
                    if isequal('Cd', identifier{i})
                        b = find(ismember(identifier(nbrs), 'Cd'));
                        if length(b) == 1
                            nbrs(b) = [];
                            clear b
                        else if length(b) > 1
                                bb = find(Species(n).adj(i,nbrs(b))==2);
                                nbrs(b(bb)) = [];
                            end
                        end

                    end
                        % also get rid of O from neighbors of C=O
                        if isequal('CO', identifier{i})
                            %isequal('O', Species(3).lbl{i}) && isempty(find(Species(3).adj(i,:)==2)))
                            b = find(ismember(identifier(nbrs), 'O'));
                            for j = 1 : length(b)
                                if ~isempty(find(Species(n).adj(nbrs(b(j)),:)==2))
                                    nbrs(b(j)) = [];
                                    break
                                end
                            end
                            clear b
                        end
                        
                        
                        % work with neighbors identifiers and their multiplicity
                        % MODIFICATION WITH T goes here
                        
                        identifier_T = identifier(nbrs);
                        if ~isempty(find(ismember(identifier_T, 'T')))
                        b = find(ismember(identifier_T, 'T'));
                        identifier_T{b} = 'C';
                        clear b
                        end
                        nbrs_unique = unique(identifier_T);
                        multiplicity=cellfun(@(x) sum(ismember(identifier_T,x)), nbrs_unique,'un',0);
                        for j = 1 : length(nbrs_unique)
                            if multiplicity{j}>1
                                nbrs_list{i} = strcat(nbrs_list{i},'(',nbrs_unique{j},')',num2str(multiplicity{j}));
                            else
                                nbrs_list{i} = strcat(nbrs_list{i},'(',nbrs_unique{j},')');
                            end
                            
                        end
                        
                    end
                end
   end
            
            % add contributions from crosslinks
            if ~isempty(find(cellfun(@isempty,Species(n).crsl)==0))
                for j = 1 : size(nbrs_crsl, 2)
                    nbrs_list{end + 1} = nbrs_crsl{j};
                end
            end
            
            clear multiplicity nbrs nbrs_unique
            %nbrs_list = nbrs_list';
            nbrs_list = nbrs_list(~cellfun('isempty',nbrs_list));
            multiplicity=cellfun(@(x) sum(ismember(nbrs_list,x)), unique(nbrs_list),'un',0);
            multiplicity = num2cell([multiplicity{:}]');
            groups = struct('group', unique(nbrs_list)', 'multiplicity', multiplicity);
            % find group in the group library
            enthalpy = 0;
            for i = 1 : length(groups)
                b = find(strcmp({groups_library.group}, groups(i).group)==1);
                enthalpy = enthalpy + groups_library(b).enthalpy_kcal * groups(i).multiplicity;
                groups(i).enthalpy = groups_library(b).enthalpy_kcal;
            end
            
            
            Species(n).groups = groups;
            Species(n).enthalpy = enthalpy;
            
%             all_groups =  [all_groups; unique(nbrs_list)'];
            clear nbrs_list multiplicity enthalpy groups
         end
     end
    clear identifier
end

%% O2 enthalpy -230 in kJ, -55 in kcal
% O2 - pattern4
for i = 1 : length(Species)
    % oxygen
    if ~isempty(intersect(4,Species(i).ptrn))
        Species(i).enthalpy = 0;
    end
    % OH.
    if ~isempty(intersect(38,Species(i).ptrn))
        Species(i).enthalpy = -55;
    end
end


%% all_groups = unique(all_groups);
for i = 1 : length(Species)
    if isempty(Species(i).enthalpy)
    
        Species(i).enthalpy = 0;
        
    end
end

clear all_groups b bb crsl crsl_ref crsl_table data groups_library identifier_T nbrs_crsl
