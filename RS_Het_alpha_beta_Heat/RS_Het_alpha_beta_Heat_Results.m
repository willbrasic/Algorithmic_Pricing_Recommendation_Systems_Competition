%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Q-learning and Sellers Q-learning Heat Map Results
% William Brasic 
% The University of Arizona
% wbrasic@arizona.edu 
% williambrasic.com
% January 2024
%
% This script obtains results for RS_Het_alpha_beta_Heat
% 
% Before executing script:
% 1. Ensure R is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clear;  

% Colors for plots
color_1 = "black";
color_2 = "cyan";
color_3 = "magenta";
color_4 = "#D95319";

% Array of colors for plots corresponding to omega value
omega_colors = [color_1; color_2; color_3; color_4];

% Do not show warnings
warning off all; 

% Numbers are rounded without scientific notation
format longG;    

% Number of episodes
E = 100;

% Run to compute results for
R = 1;

% Number of firms
n = 2;

% Platform profit weights
omega = [0, 4/5, 1];

% Platform fee
f = 0.2;

% Overall version
version = "RS_Het_alpha_beta_Heat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Equilibrium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logit version
logit_version =  "RS_Het";

% Add path to logit equilibrium solver scripts
addpath(strcat("LE_", logit_version));

% Solve for logit collusive equilibrium
run(strcat("LE_", logit_version, "\LE_Coll_", logit_version, ".m"));

% Solve for logit competitive equilibrium
run(strcat("LE_", logit_version, "\LE_Comp_", logit_version, ".m"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path to with platform files
for w = 1:length(omega)
    addpath(strcat(version, "\", version, "_omega_", num2str(w), "_Results"));
end

% File name for results for with platform files
p_sellers_file_names = cell(1, length(omega));
cs_file_names = cell(1, length(omega));
a_plat_file_names = cell(1, length(omega));
for w = 1:length(omega)
    p_sellers_file_names{w} = strcat(version, "\", ...
        version, "_omega_", num2str(w), "_Results\", version, "_omega_", num2str(w), "_Prices_Sellers.csv");
    cs_file_names{w} = strcat(version, "\", ...
        version, "_omega_", num2str(w), "_Results\", version, "_omega_", num2str(w), "_CS.csv");
    a_plat_file_names{w} = strcat(version, "\", ...
        version, "_omega_", num2str(w), "_Results\", version, "_omega_", num2str(w), "_Actions_Plat.csv");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' prices with platform
p_sellers = zeros(size(table2array(readtable(p_sellers_file_names{w})), 1), size(table2array(readtable(p_sellers_file_names{w})), 2), length(omega));
for w = 1:length(omega)
    p_sellers(:, :, w) = table2array(readtable(p_sellers_file_names{w}));
end

% Find unique (alpha, beta) pairs across all slices to ensure consistent size
alpha = p_sellers(:, 2, 1);
beta = p_sellers(:, 3, 1);
unique_combinations = unique([alpha, beta], "rows");

% Mean price across firms for given (alpha, beta) pair
num_rows = size(unique_combinations, 1);
p_sellers_mean = zeros(num_rows, 3, length(omega)); 
for w = 1:length(omega)
    % Extract p_sellers_mean_relative_w for current omega (w)
    p_sellers_w = p_sellers(:, :, w);
    
    % Extract relevant columns (ignore theta_2 and firm)
    alpha = p_sellers_w(:, 2);
    beta = p_sellers_w(:, 3);
    prices = p_sellers_w(:, 5);
    
    % Create unique combinations of alpha and beta
    unique_combinations = [alpha, beta];
    [unique_vals, ~, idx] = unique(unique_combinations, "rows");

    % Compute mean prices for each unique (alpha, beta) group
    mean_prices = accumarray(idx, prices, [], @mean);

    % Store results in preallocated array
    p_sellers_mean(:, :, w) = [unique_vals, mean_prices];
end

% Preallocate the results array for relative mean calculations
p_sellers_mean_relative = zeros(size(p_sellers_mean, 1), 3, length(omega));

for w = 1:length(omega)
    p_sellers_mean_relative(:, 1, w) = p_sellers_mean(:, 1, w); 
    p_sellers_mean_relative(:, 2, w) = p_sellers_mean(:, 2, w); 
    p_sellers_mean_relative(:, 3, w) = p_sellers_mean(:, 3, w) ./ mean(comp_p);
end


%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%

% Create a figure with a 2x2 tiled layout and compact spacing
figure;
t = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% Extract unique alpha and beta values 
unique_alpha = unique(cs_mean_relative(:, 1, 1));
unique_beta = unique(cs_mean_relative(:, 2, 1));

% Identify the index for omega == 4/5 and the indices for the other omegas
idx_omega_4_5 = find(omega == 4/5);
other_indices = find(omega ~= 4/5);

% Plot the ω = 4/5 heatmap in the first row spanning both columns
if ~isempty(idx_omega_4_5)
    ax = nexttile(t, [1 2]);
    w = idx_omega_4_5;
    
    cs_mean_relative_w = cs_mean_relative(:, :, w);
    alpha = cs_mean_relative_w(:, 1);
    beta = cs_mean_relative_w(:, 2);
    relative_cs = cs_mean_relative_w(:, 3);
    
    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = relative_cs(idx);
        end
    end
    
    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title("\omega = 4/5", "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Plot the other omega heatmaps in the second row
for k = 1:length(other_indices)
    w = other_indices(k);
    ax = nexttile(t);
    
    cs_mean_relative_w = cs_mean_relative(:, :, w);
    alpha = cs_mean_relative_w(:, 1);
    beta = cs_mean_relative_w(:, 2);
    relative_cs = cs_mean_relative_w(:, 3);
    
    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = relative_cs(idx);
        end
    end
    
    omega_label = num2str(omega(w));
    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title(strcat("\omega = ", omega_label), "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Add a colorbar below all tiles
h = colorbar("southoutside");
h.Label.String = "Relative Consumer Surplus";
h.Layout.Tile = "south";

% Save prices heat map
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_CS.png"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consumer Surplus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
% With Platform
%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs with platform
cs = zeros(size(table2array(readtable(cs_file_names{w})), 1), size(table2array(readtable(cs_file_names{w})), 2), length(omega));
for w = 1:length(omega)
    cs(:, :, w) = table2array(readtable(cs_file_names{w}));
end

% Find unique (alpha, beta) pairs across all slices to ensure consistent size
alpha = cs(:, 2, 1);
beta = cs(:, 3, 1);
unique_combinations = unique([alpha, beta], "rows");

% Preallocate the results array
num_rows = size(unique_combinations, 1);
cs_mean = zeros(num_rows, 3, length(omega));  % Columns: alpha, beta, mean_cs

for w = 1:length(omega)
    % Extract cs_mean_relative_w for current omega (w)
    cs_w = cs(:, :, w);
    
    % Extract relevant columns (ignore theta_2)
    alpha = cs_w(:, 2);
    beta = cs_w(:, 3);
    cs_now = cs_w(:, 4);
    
    % Create unique combinations of alpha and beta
    unique_combinations = [alpha, beta];
    [unique_vals, ~, idx] = unique(unique_combinations, "rows");

    % Compute mean cs for each unique (alpha, beta) group
    mean_cs = accumarray(idx, cs_now, [], @mean);

    % Store results in preallocated array
    cs_mean(:, :, w) = [unique_vals, mean_cs];
end

% Preallocate the results array for relative mean calculations
cs_mean_relative = zeros(size(cs_mean, 1), 3, length(omega));

for w = 1:length(omega)
    cs_mean_relative(:, 1, w) = cs_mean(:, 1, w); 
    cs_mean_relative(:, 2, w) = cs_mean(:, 2, w); 
    cs_mean_relative(:, 3, w) = cs_mean(:, 3, w) ./ comp_cs;
end


%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%

% Create a figure with a 2x2 tiled layout and compact spacing
figure;
t = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% Extract unique alpha and beta values 
unique_alpha = unique(cs_mean_relative(:, 1, 1));
unique_beta = unique(cs_mean_relative(:, 2, 1));

% Find index of omega = 4/5
idx_omega_4_5 = find(omega == 4/5);
other_indices = find(omega ~= 4/5);

% Plot omega = 4/5 in a centered position spanning both columns in first row
if ~isempty(idx_omega_4_5)
    ax = nexttile(t, [1 2]);
    w = idx_omega_4_5;

    cs_mean_relative_w = cs_mean_relative(:, :, w);
    alpha = cs_mean_relative_w(:, 1);
    beta = cs_mean_relative_w(:, 2);
    relative_cs = cs_mean_relative_w(:, 3);

    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = relative_cs(idx);
        end
    end

    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title("\omega = 4/5", "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Plot the other omega values in the second row
for k = 1:length(other_indices)
    w = other_indices(k);
    ax = nexttile(t);

    cs_mean_relative_w = cs_mean_relative(:, :, w);
    alpha = cs_mean_relative_w(:, 1);
    beta = cs_mean_relative_w(:, 2);
    relative_cs = cs_mean_relative_w(:, 3);

    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = relative_cs(idx);
        end
    end

    omega_label = num2str(omega(w));
    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title(strcat("\omega = ", omega_label), "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Add colorbar below all tiles
h = colorbar("southoutside");
h.Label.String = "Relative Consumer Surplus";
h.Layout.Tile = "south";

% Save cs heat map
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_CS.png"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Actions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a figure with a 2x2 tiled layout and compact spacing
figure;
t = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% Extract unique alpha and beta values 
unique_alpha = unique(a_plat(:, 1, 1));
unique_beta = unique(a_plat(:, 2, 1));

% Find index of omega = 4/5
idx_omega_4_5 = find(omega == 4/5);
other_indices = find(omega ~= 4/5);

% Plot omega = 4/5 in a centered position spanning both columns in first row
if ~isempty(idx_omega_4_5)
    ax = nexttile(t, [1 2]);
    w = idx_omega_4_5;

    a_plat_w = a_plat(:, :, w);
    alpha = a_plat_w(:, 1);
    beta = a_plat_w(:, 2);
    count = a_plat_w(:, 3);

    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = count(idx);
        end
    end

    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title("\omega = 4/5", "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Plot the other omega values in the second row
for k = 1:length(other_indices)
    w = other_indices(k);
    ax = nexttile(t);

    a_plat_w = a_plat(:, :, w);
    alpha = a_plat_w(:, 1);
    beta = a_plat_w(:, 2);
    count = a_plat_w(:, 3);

    heatmap_matrix = zeros(length(unique_beta), length(unique_alpha));
    for i = 1:length(unique_beta)
        for j = 1:length(unique_alpha)
            idx = (alpha == unique_alpha(j)) & (beta == unique_beta(i));
            heatmap_matrix(i, j) = count(idx);
        end
    end

    omega_label = num2str(omega(w));
    imagesc(unique_alpha, unique_beta, heatmap_matrix);
    xlabel("\alpha");
    ylabel("\beta");
    title(strcat("\omega = ", omega_label), "FontWeight", "normal");
    set(gca, "YDir", "normal");
    colormap jet;
end

% Add colorbar below all tiles
h = colorbar("southoutside");
h.Label.String = "Proportion of Action 2 Selections";
h.Layout.Tile = "south";

% Save actions heat map
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_Platform_Actions.png"));










