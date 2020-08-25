function  recombination_O_C( in, out )
% recombination of O. with all carbon radicals (10, 19, 24, 25)
global pattern;
global pattern_non_reactive;
global R22;

   
if in>=8
    R22(end+1).reactant1 = pattern(in); % RO.
    R22(end).reactant2 = pattern(8); % C.
    R22(end).product1 = pattern(out); % R-O-...
    R22(end).product2 = pattern_non_reactive(12); % R-...
    R22(end).k = 0;
    R22(end).type = 18;
end

if in>=19
    R22(end+1).reactant1 = pattern(in); % RO.
    R22(end).reactant2 = pattern(19); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(14); % R-...
    R22(end).k = 0;
    R22(end).type = 18;
end
if in>=20
    R22(end+1).reactant1 = pattern(in); % RO.
    R22(end).reactant2 = pattern(20); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(28); % R-...
    R22(end).k = 0;
    R22(end).type = 18;
end
end

