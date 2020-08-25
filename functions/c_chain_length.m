function [ c_length ] = c_chain_length( labels )
% c_chain_length calculates carbon chain length of the molecule
c = find(ismember(labels, 'C'));
c_length = length(c);

end

