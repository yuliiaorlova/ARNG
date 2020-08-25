function [ pattern_list ] = pattern_non_reactive_struct(  )

% intialize structure for patters

pattern_list = struct('adj', [], 'lbl', [], 'id', [], 'crsl', []);

for i = 1 : 34
    
    adj_name =  sprintf('adjacency_non_r_pattern_%d.txt', i);
    lbl_name = sprintf('atoms_non_r_pattern_%d.txt', i);
    crsl_name = sprintf('crosslinks_non_r_pattern_%d.txt', i);
    
    A = importdata(adj_name);
    labels = importdata(lbl_name);
    crsl = importdata(crsl_name);
    
    pattern_list(i).adj = A;
    pattern_list(i).lbl = labels;
    pattern_list(i).id = i+100;
    pattern_list(i).crsl= cell(length(crsl),1);
    
    n_cl = find(crsl~=0 );
    for j = 1 : length(n_cl)
            pattern_list(i).crsl{n_cl(j)}(end+1) = crsl(n_cl(j));
    end
    
    clear adj_name lbl_name crsl_name
    clear A labels crsl n_cl 
end

end

