function  recombination_C_OO( in, out )
% recombination of C. with all peroxy radicals (11, 26)
global pattern;
global R22;
   
if in>=9
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(9); % COO.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(14); % R-...
    R22(end).k = 0;
    R22(end).type = 16;

end
if in>=21
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(21); % COO.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(27); % R-...
    R22(end).k = 0;
    R22(end).type = 16;

end

end

