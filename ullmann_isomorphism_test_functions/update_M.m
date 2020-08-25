function [ M, change ] = update_M( M, A_pattern, A_molecule )
change = [];

for i = 1 : size(M,1)
    for j = 1 : size(M,2)
    
        if M(i,j)
        
            n1 = neighbors(A_pattern, i);
            n2 = neighbors(A_molecule, j);
            l1 = A_pattern(i, n1);
            l2 = A_molecule(j, n2);
            
            l1_un = unique(l1);
            l2_un = unique(l2);
            
            % check links
            if length(l1_un)>length(l2_un) && prod(ismember(l1_un, l2_un)) == 0 % ~= to > , 
                M(j,i) = 0;
                change = [change; i];
                break
            end
            for k = 1 : length(l1_un)
                if length(find(l1== l1_un(k))) > length(find(l2== l1_un(k))) % change was here
                    M(j,i) = 0;
                    change = [change; i];
                    break
                    
                end
                
            end

            for k = 1 : length(n1)
            
                if sum(M(n1(k), n2))>=1 % check if at least one neighbor from n2 has 1 with neighbor n1(k)
                    
                    continue
                else 
                    M(i,j) = 0;
                    change = [change; i];
                end
            
            end
        
        end
    
    end
end

end

