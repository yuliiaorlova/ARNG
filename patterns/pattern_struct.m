function [ pattern_list ] = pattern_struct(  )

% intialize structure for patters

pattern_list = struct('adj', [], 'lbl', [], 'id', [], 'crsl', []);


% 1 - Co(II) initiator
pattern_list(1).adj = 1;
pattern_list(1).lbl = 'Co2';
pattern_list(1).id = length(pattern_list);
pattern_list(1).crsl = [];

% 2 -Co(III)OO
pattern_list(end+1).adj = 1;
pattern_list(end).lbl = 'Co3OO';
pattern_list(end).id = length(pattern_list);
pattern_list(end).crsl = [];

% 3 - Co(III)OOH
pattern_list(end+1).adj = 1;
pattern_list(end).lbl = 'Co3OOH';
pattern_list(end).id = length(pattern_list);
pattern_list(end).crsl = [];

% 4 - O2
pattern_list(end + 1).adj = 1;
pattern_list(end).lbl = 'O2';
pattern_list(end).id = length(pattern_list);
pattern_list(end).crsl = [];

% 5 - Co(III)OH
pattern_list(end + 1).adj = 1;
pattern_list(end).lbl = 'Co3OH';
pattern_list(end).id = length(pattern_list);
pattern_list(end).crsl = [];

% 6 - H2O
pattern_list(end + 1).adj = 1;
pattern_list(end).lbl = 'H2O';
pattern_list(end).id = length(pattern_list);
pattern_list(end).crsl = [];


for i = 7 : 40
    
    adj_name =  sprintf('adjacency_pattern_%d.txt', i);
    lbl_name = sprintf('atoms_pattern_%d.txt', i);
    crsl_name = sprintf('crosslinks_pattern_%d.txt', i);
    
    A = importdata(adj_name);
    labels = importdata(lbl_name);
    crsl = importdata(crsl_name);
    
    
    
    pattern_list(end+1).adj = A;
    pattern_list(end).lbl = labels;
    pattern_list(end).id = length(pattern_list);
    pattern_list(i).crsl= cell(length(crsl),1);
    n_cl = find(crsl >0 );
    for j = 1 : length(n_cl)
            pattern_list(i).crsl{n_cl(j)}(end+1) = crsl(n_cl(j));
    end
    clear A;
    clear labels;
    clear adj_name;
    clear lbl_name crsl crsl_name
end



end

