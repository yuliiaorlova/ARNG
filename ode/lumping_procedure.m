function [ p_y ] = lumping_procedure( y )
%LUMPING_PROCEDURE lumps species according to patterns they have
% matrix p_y corresponds to lumped concentrations of all species
% column corresponds to patterms
% row corresponds to time step
global Species pattern
p_y = zeros(size(y,1), length(pattern));
for i = 1 : length(pattern)
    for j = 1 : size(y,2)
        if ~isempty(intersect(i, Species(j).ptrn))
           p_y(:, i) = p_y(:,i) + y(:,j);
        end
    end
    
end

