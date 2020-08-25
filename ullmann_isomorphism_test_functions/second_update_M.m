function [ M, dont_match ] = second_update_M( M, A_pattern, A_molecule )
%SECOND_M Summary of this function goes here
dont_match = 0;
candidates = cell(size(M,1), 1);
    for i = 1 : size(M,1)
        candidates{i,1} = find(M(i,:)==1);
    end
    
    % refine M matrix
    for i = 1 : length(candidates) % i is for pattern
        % choose i as a root
        %q = b_first_search(A_pattern, i);
        if length(candidates{i})>1
            
            for j = 1 : length(candidates{i}) % j is for molecule
                % set all entries of row i and column candidates{i}(j), keep M_(i,candidates{i}(j)) = 1
                M_ = M;
                M_(:,candidates{i}(j)) = 0;
                M_(i, :) = 0;
                M_(i,candidates{i}(j)) = 1;
                % matching step
                change  = 0; % variable which checks whether there were changes to the M matrix on one loop or not
                
                while ~isempty(change)
                    [M_, change ]= update_M(M_, A_pattern, A_molecule);
                    % check for empty rows
                    if empty_rows(M_)
                        M(i,candidates{i}(j)) = 0;
                        break
                        
                    end
                end
            end
        end
    end
    
   if empty_rows(M)
       dont_match = 1;
   end
    
end

