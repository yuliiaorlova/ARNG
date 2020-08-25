%% need to make a library of enthalpy for groups
data = importdata('group_additivity_values.txt');
groups_library = struct('group', cell(length(data.data),1), 'enthalpy_kcal',  [], 'enthalpy_kJ',  []);
for i = 1 : length(data.data)

  groups_library(i).group = strrep(data.textdata{i+1},'''','');
  groups_library(i).enthalpy_kcal = data.data(i,1);
  groups_library(i).enthalpy_kJ = data.data(i,2);

end



