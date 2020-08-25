function [ M, dont_match ] = first_update_M( M, A_pattern, A_molecule )
%FIRST_UPDATE_M - initial update after initialization 

change  = 0; % variable which checks whether there were changes to the M matrix on one loop or not
count = 0;
dont_match = 0;
while ~isempty(change)
    count = count + 1; % count number of cycles
    [M, change ]= update_M(M, A_pattern, A_molecule);
    % check for empty rows
    if empty_rows(M)
        dont_match = 1;
        break
    end
end


end

