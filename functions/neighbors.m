function [ neighbors ] = neighbors( A, root )
neighbors = find(A(root, :)>=1);
neighbors = neighbors(neighbors~=root);
end

