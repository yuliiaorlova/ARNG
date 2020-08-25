function  recombination_O_O( in, out )
% recombination of O. with all alkoxy radicals (13, 34, 39, 42, 45, 47, 48)
global pattern;
global R22;
   
   
if in>=11
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(11); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(14); % R-...
    R22(end).k = 0;
    R22(end).type = 17;
end
if in>=29
    R22(end+1).reactant1 = pattern(in); % R.
    R22(end).reactant2 = pattern(29); % C.
    R22(end).product1 = pattern(out); % R-...
    R22(end).product2 = pattern(27); % R-...
    R22(end).k = 0;
    R22(end).type = 17;
end

end

