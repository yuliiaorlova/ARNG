function  plot_functional_groups( t2, y2, O2_upt, PV, hex_sum, pent_sum )
global Species
figure(11)
hFig = figure(11);
set(gcf,'PaperPositionMode','auto')
set(gcf,'color','w');
set(gca,'FontSize',16); 
set(hFig, 'Position', [0 0 1600 900])

hold on
f1 = subplot(3, 4, 1); % EL
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f2 = subplot(3, 4, 2); % O2
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f3 = subplot(3, 4, 3); % PV
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f4 = subplot(3, 4, 4); %  hexana and pentanal
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f5 = subplot(3, 4, 5); % crosslinks
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f6 = subplot(3, 4, 6); % hydroxides 
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f7 = subplot(3, 4, 7); % aldehydes
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f8 = subplot(3, 4, 8); % conjugated db water
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f9 = subplot(3, 4, 9); % carboxylic acid
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f10 = subplot(3, 4, 10); % radicals
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f11 = subplot(3, 4, 11); % water
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on
f12 = subplot(3, 4, 12); % epoxides
set(gca, 'XScale', 'log', 'FontSize', 16);
hold on

crosslinks = [14, 15, 16, 26, 27, 28, 39];
radicals = [8, 9, 11, 19, 20, 21, 29, 39];
c_ba_0= 3.7106;

p_y = lumping_procedure(y2);
p_y = p_y/c_ba_0;
p_y(:, crosslinks) = p_y(:, crosslinks)/2;
alkyl_cl = [16, 26];
ether_cl = [15, 28];
peroxy_cl = [14, 27, 39];

%1 EL
semilogx(f1, t2(:), p_y(:, 7), 'LineWidth',2);
xlim(f1, [min(t2) max(t2)])
title(f1, 'EL')
xlabel(f1, 'Time (hours)')
ylabel(f1, 'Concentration (mol/mol EL)')

if ~isempty(O2_upt)
    %2 O2
    semilogx(f2, t2(:), O2_upt(:), 'LineWidth',2);
    title(f2, 'Oxygen uptake')
    xlim(f2, [min(t2) max(t2)])
    xlabel(f2, 'Time (hours)')
    ylabel(f2, 'Concentration (mol/mol EL)')
end
%3 PV
% OOH = [10, 22];
% OO_cl = [14, 27, 39];
% peroxides = [2, 3, 9, 10, 14, 21, 22, 25, 27, 32, 39];
% pv_cl = [14, 27, 39];
% p_y(:, crosslinks) = p_y(:, crosslinks)/2;
% PV = sum(p_y(:, peroxides), 2);
semilogx(f3, t2(:), PV(:), 'LineWidth',2);
title(f3, 'Peroxide value')
xlim(f3, [min(t2) max(t2)])
xlabel(f3, 'Time (hours)')
ylabel(f3, 'Concentration (mol/mol EL)')

%4 hydroxides (alcohols)
semilogx(f4, t2(:), p_y(:, 13), 'LineWidth',2);
title(f4, 'Alcohols')
xlim(f4, [min(t2) max(t2)])
xlabel(f4, 'Time (hours)')
ylabel(f4, 'Concentration (mol/mol EL)')



% find all double bonds concentrations 

n_db = zeros(length(Species), 1);
for i = 1 : length(Species)

[I, ~] =  ind2sub(length(Species(i).adj), find(Species(i).adj==2));
if ~isempty(I)
n_db(i) = 0.5*(0.5*length(I) - sum(ismember(Species(i).lbl(I),'O')));
end

end
y2_db = zeros(length(t2), length(Species));
for i = 1 : length(n_db)

    y2_db(:, i) = y2(:, i)*n_db(i);

end
y2_db_sum = sum(y2_db, 2);


semilogx(f5, t2(:),y2_db_sum/c_ba_0, 'LineWidth',2);
title(f5, 'All double bonds')
xlim(f5, [min(t2) max(t2)])
xlabel(f5, 'Time (hours)')
ylabel(f5, 'Concentration (mol/mol EL)')

% %5 all crosslinks
% semilogx(f5, t2(:), sum(p_y(:, crosslinks),2));
% title(f5, 'All crosslinks')
% xlabel(f5, 'Time (hours)')
% ylabel(f5, 'Concentration (mol/mol EL)')
% hold on

%6 conjugated double bonds
semilogx(f6, t2(:), p_y(:, 12), 'LineWidth',2);
title(f6, 'Conj. double bonds')
xlim(f6, [min(t2) max(t2)])
xlabel(f6, 'Time (hours)')
ylabel(f6, 'Concentration (mol/mol EL)')

%7 aldehydes
semilogx(f7, t2(:), p_y(:, 18), 'LineWidth',2);
title(f7, 'Aldehydes')
xlim(f7, [min(t2) max(t2)])
xlabel(f7, 'Time (hours)')
ylabel(f7, 'Concentration (mol/mol EL)')


%8 crosslinks
semilogx(f8, t2(:), sum(p_y(:, alkyl_cl),2), t2(:), sum(p_y(:, ether_cl),2), t2(:), sum(p_y(:, peroxy_cl),2), 'LineWidth',2); %, t2(:), pent_sum(:));
title(f8, 'crosslinks')
xlim(f8, [min(t2) max(t2)])
xlabel(f8, 'Time (hours)')
ylabel(f8, 'Concentration (mol/mol EL)')



%9 carboxylic acid
semilogx(f9, t2(:), p_y(:, 23), 'LineWidth',2);
title(f9, 'Carboxylic acid')
xlim(f9, [min(t2) max(t2)])
xlabel(f9, 'Time (hours)')
ylabel(f9, 'Concentration (mol/mol EL)')


%10 radicals
semilogx(f10, t2(:), sum(p_y(:, radicals),2), 'LineWidth',2);
title(f10, 'Radicals')
xlim(f10, [min(t2) max(t2)])
xlabel(f10, 'Time (hours)')
ylabel(f10, 'Concentration (mol/mol EL)')


%11 water
semilogx(f11, t2(:), p_y(:, 6), 'LineWidth',2);
title(f11, 'Water')
xlim(f11, [min(t2) max(t2)])
xlabel(f11, 'Time (hours)')
ylabel(f11, 'Concentration (mol/mol EL)')


%12 epoxides
semilogx(f12, t2(:), p_y(:, 37), 'LineWidth',2);
title(f12, 'Epoxides')
xlim(f12, [min(t2) max(t2)])
xlabel(f12, 'Time (hours)')
ylabel(f12, 'Concentration (mol/mol EL)')



end

