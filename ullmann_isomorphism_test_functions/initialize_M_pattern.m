function [ M ] = initialize_M_pattern( A_pattern, labels_pattern, crosslinks_pattern, A_molecule, labels_molecule, crosslinks_molecule )

M = zeros(length(A_pattern), length(A_molecule));

for i = 1 : length(A_molecule) % i for molecule
    for j = 1 : length(A_pattern) % j for pattern
    
        if isequal(labels_molecule(i), labels_pattern(j))
            
            if length(neighbors(A_pattern, j))<= length(neighbors(A_molecule, i)) % = for <
            
                n1 = labels_pattern(neighbors(A_pattern, j)); % labels
                n2 = labels_molecule(neighbors(A_molecule, i));
                
                l1 = A_pattern(j, neighbors(A_pattern, j)); % links
                l2 = A_molecule(i, neighbors(A_molecule, i));
                
                n1_un = unique(n1);
                n2_un = unique(n2);
                
                l1_un = unique(l1);
                l2_un = unique(l2);
                
                if length(n1_un)<=length(n2_un) && length(l1_un)<=length(l2_un) && prod(ismember(n1_un, n2_un)) == 1 && prod(ismember(l1_un, l2_un)) == 1 && A_molecule(i,i) == A_pattern(j,j)% (for first two) = for < 
                    count1 = 0;
                    count2 = 0;
                    for k = 1 : length(n1_un)
                        if length(find(strcmp(n1, n1_un(k))))<= length(find(strcmp(n2, n1_un(k)))) % = for <
                            count1 = count1+1; % count labels
                        else 
                            break
                        end
                    end
                    for k = 1 : length(l1_un)
                        if length(find(l1== l1_un(k)))<= length(find(l2== l1_un(k))) % = for <
                            count2 = count2+1; % count labels
                        else
                            break
                        end
                    end
                    if count1 == length(n1_un) && count2 == length(l1_un)
                        if isequal(crosslinks_molecule(i), crosslinks_pattern(j) ) || isempty(crosslinks_molecule{i})&& isempty(crosslinks_pattern{j})
                            M(j,i) = 1;
                            
                        else if isempty(crosslinks_pattern{j}) &&  ~isempty(crosslinks_molecule{i}) 
                                M(j,i) = 1;
                                
                            else if ~isempty(crosslinks_pattern{j}) && ismember(crosslinks_pattern{j}, crosslinks_molecule{i})
                                    M(j,i) = 1;
                                end
                            end
                        end
                        
                    end
                end
                
            end
        end
    end
end


end

