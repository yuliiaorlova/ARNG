function  recombination_OO_C( in, out )
% recombination of OO. with all carbon radicals (10, 19, 24, 25)
global pattern;
global pattern_non_reactive;
global R22;

   
if in>=8
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(8); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(11); % R-...
    R22(end).k = 0;
    R22(end).type = 16;
end

if in>=19
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(19); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern_non_reactive(9); % R-...
    R22(end).k = 0;
    R22(end).type = 16;
end
if in>=20
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(20); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(27); % R-...
    R22(end).k = 0;
    R22(end).type = 16;
end
end

