function  recombination_C_C( in, out )
global pattern;
global pattern_non_reactive;
global R22;
global ks;
   
if in<=10
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(10); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(3); % R-...
    R22(end).k = ks(12);
    R22(end).type = length(R22) + 1;
end
if in<=10
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(10); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(21); % R-...
    R22(end).k = ks(12);
    R22(end).type = length(R22) + 1;
end
if in<=19    
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(19); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(16); % R-...
    R22(end).k = ks(12);
    R22(end).type = length(R22) + 1;
end
if in<=24
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(24); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(11); % R-...
    R22(end).k = ks(12);
    R22(end).type = length(R22) + 1;
end
if in<=25
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(25); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(31); % R-...
    R22(end).k = ks(12);
    R22(end).type = length(R22) + 1;
end
end

