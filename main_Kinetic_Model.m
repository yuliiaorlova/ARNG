addpath('scripts_to_make_plots/')
addpath('ode/');

global hexanal pentanal c_hexanal_equil c_pentanal_equil k_hexanal_transfer k_pentanal_transfer oxygen_containing_species diff

hexanal = [];
pentanal = [];
crosslinks = [14, 15, 16, 26, 27, 28, 39];

k_hexanal_transfer = 10;
k_pentanal_transfer = 10;
c_hexanal_equil = 0;
c_pentanal_equil = 0;


%% find hexanal and pentanal species
% hexanal and pentanal volatilize, thus need species tratment in the
% kinetic model

 for i = 1 : length(Species) 
     if c_chain_length(Species(i).lbl)== 5 && isempty(intersect(crosslinks, Species(i).ptrn)) && ~isempty(intersect(18, Species(i).ptrn))
         pentanal = [pentanal; i];
     else if c_chain_length(Species(i).lbl) == 6  && isempty(intersect(crosslinks, Species(i).ptrn)) && ~isempty(intersect(18, Species(i).ptrn))
             hexanal = [hexanal; i];             
         end
     end
 end

%% find oxygen containing species
% this block is needed to extract correct oxygen uptake value
oxygen_containing_species = zeros(length(Species), 1);
for i = 1 : length(Species)
    
    oxygen_containing_species(i) = sum(ismember(Species(i).lbl, 'O'));
    
    % take care of ethyl group with 2 oxygens
    if ~isempty(intersect(Species(i).ptrn, 40))
        oxygen_containing_species(i) = oxygen_containing_species(i) - 2;
    end
    
    % take care of crosslinks: ether
    if ~isempty(intersect(Species(i).ptrn, 15)) || ~isempty(intersect(Species(i).ptrn, 28))
        oxygen_containing_species(i) = oxygen_containing_species(i) + 0.5; % one oxygen per crosslink, half per monomer
    end
    
    % take care of crosslinks: peroxy
    if ~isempty(intersect(Species(i).ptrn, 14)) || ~isempty(intersect(Species(i).ptrn, 27))
        oxygen_containing_species(i) = oxygen_containing_species(i) + 1; % two oxygens per crosslink, one per monomer
    end
    
    oxygen_containing_species(i)= oxygen_containing_species(i)/2; % because one molecule is O2
end
oxygen_containing_species(2) = 1;


diff = 0;


[O2_upt, PV, p_y, hex_sum, pent_sum, y2, t2] = model_run_1(  );

save('output/Kinetic_model_EL.mat')

%%

% visualize pos and neg fluxes for functional groups
plot_functional_groups( t2, y2, O2_upt, PV, [], [] )
