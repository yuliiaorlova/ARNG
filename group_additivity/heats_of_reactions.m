load('25_09_no_di_tri.mat');
run('calculate_heat_of_formation.m');

% calculate reate parameters according to heats of reactions and
% temperature

T = 298;
value = zeros(length(weight), 1);
Hrxn = zeros(length(weight), 1);
H1 = zeros(length(Reactions),1);
H2 = zeros(length(Reactions),1);
for i = 1 : length(Reactions)
    % each i corresponds to individual reaction
    %i = 1;
    indicator = 0;
%     X = [' Reaction: ', num2str(i) ,'; reactants: ', sprintf('%d,' , find(A(:,length(Molecules)+i)>0)), '; products:', sprintf('%d,', find(A(length(Molecules)+i, :)>0)) ];
%     disp(X)
    reactants = find(A(:,length(Molecules)+i)>0);
    products = find(A(length(Molecules)+i, :)>0);
    for j = 1 : length(reactants)
        if Species(reactants(j)).enthalpy == 0 && reactants(j)~=2
            indicator = 1;
        end
        H1(i) = H1(i) + Species(reactants(j)).enthalpy*A(reactants(j),length(Molecules)+i);
    end
    for j = 1 : length(products)
        if Species(products(j)).enthalpy == 0 && products(j)~=2
            indicator = 1;
        end
        H2(i) = H2(i) + Species(products(j)).enthalpy*A(length(Molecules)+i, products(j));
    end
    if indicator
        Hrxn(length(Molecules)+i) = 0;
    else
    Hrxn(length(Molecules)+i) = H2(i) - H1(i);
    end
    
    value(length(Molecules)+i) = rate_constant( T, type(length(Molecules)+i), Hrxn(length(Molecules)+i) );
%     % reactants
%     find(A(:,length(Molecules)+i)>0)
%     % products
%     find(A(length(Molecules)+i, :)>0)
end

weight(length(Molecules)+1: length(weight)) = value(length(Molecules)+1: length(weight));
compare = [weight, value, Hrxn, type];
save('25.09_no_di_tri_RN_Hrxn.mat')
run('main_run.m');