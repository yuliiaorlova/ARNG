function [ match_all ] = match_index( M, A_pattern, A_molecule, labels_pattern )
%MATCH_INDEX Summary of this function goes here
match_all = {};
for i = 1 : 1 %size(M,1) % atoms pattern
    
    for j = 1 : size(M,2) % atoms molecule
        
        if M(i, j)~= 0
            matched = zeros(length(A_pattern),2);
            root_pattern = i;
            root_molecule = j;
            matched(:, 1) = b_first_search( A_pattern, root_pattern);
            matched(1, 2) = root_molecule;
            
            for k = 2 : length(matched)
                
                candidates_ = find(M(matched(k,1),:)~=0);
                for l = 1 : length(candidates_)
                    [~, a2]= ismember(intersect(neighbors(A_pattern, matched(k, 1)), matched(1:k,1)), (matched(:,1)));
                    
                    if ~isempty(intersect(neighbors(A_molecule, candidates_(l)), matched(a2, 2))) && isempty(intersect(candidates_(l), matched(:,2)))
                        matched(k,2) = candidates_(l);
                    end
                end
            end
            if isempty(find(matched(:,2)==0))
                indicator = 0;
                if isempty(match_all)
                    match_all = [match_all; matched];
                    indicator = 1;
                else
                    for k = 1 : length(match_all)
                        if length(intersect(match_all{k}(:,2), matched)) == length(A_pattern) || length(intersect(match_all{k}(:,2), matched)) == length(A_pattern)-length(find(ismember(labels_pattern, 'H')))
                            indicator = 1; % if indicator is 1 it means that this match is already there
                            break
                        end
                    end
                end
                
                if ~indicator
                    match_all = [match_all; matched];
                end
            end
        end
    end
end
end

