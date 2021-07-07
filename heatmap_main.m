% HEATMAP C NO NOISE %

% filename = 'table_for_heatmap_c_no_noise.txt';
% tbl = readtable(filename, 'VariableNamingRule', 'preserve');
% 
% h = heatmap(tbl,'ci/cr', 'di/dr', 'ColorVariable', 'cf/cr', 'CellLabelColor', 'none');
% 
% ax = gca;
% ax.FontSize = 25;
% 
% h.XLabel = 'c_{i}/c_{r}';
% h.YLabel = 'd_{i}/d_{r}'; 
% h.Title = 'c_{f}/d_{r}';
% h.ColorLimits = [0.30 2.0367];
% h.Colormap = parula;

% HEATMAP D NO NOISE %
% 
filename = 'table_for_heatmap_d_noise.txt';
tbl = readtable(filename, 'VariableNamingRule', 'preserve');

h = heatmap(tbl,'ci/cr', 'di/dr', 'ColorVariable', 'df/dr', 'CellLabelColor', 'none');

ax = gca;
ax.FontSize = 25;

h.XLabel = 'c_{i}/c_{r}';
h.YLabel = 'd_{i}/d_{r}'; 
h.Title = 'd_{f}/d_{r}';
h.ColorLimits = [0.2500 20];
h.Colormap = summer;

