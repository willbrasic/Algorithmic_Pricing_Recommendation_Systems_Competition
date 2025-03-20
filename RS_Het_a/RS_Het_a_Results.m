%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Q-learning and Sellers Q-learning Varying a Results
% William Brasic 
% The University of Arizona
% wbrasic@arizona.edu 
% williambrasic.com
% January 2024
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
color_3 = "green";
color_4 = "magenta";

% Array of colors for plots corresponding to omega value
colors = [color_1; color_2; color_3; color_4];

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
version = "RS_Het_a";

% Add path
addpath(version + "\" + version + "_Results");
addpath("LE_" + version);
addpath("No_" + version + "\" + "No_" + version + "_Results");
addpath("LE_No_" + version);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File names for results
p_sellers_file_name = version + "\" + version + "_Results" + "\" + version + "_Prices_Sellers.csv";
d_sellers_file_name = version + "\" + version + "_Results" + "\" + version + "_Demand_Sellers.csv";
cs_file_name = version + "\" + version + "_Results" + "\" + version + "_CS.csv";
plat_actions_file_name = version + "\" + version + "_Results" + "\" + version + "_Actions_Plat.csv";

% File names for competitive equilibrium results
p_sellers_comp_file_name = version + "_Comp_Price.mat";
d_sellers_comp_file_name = version + "_Comp_Demand.mat";
cs_comp_file_name = version + "_Comp_CS.mat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File names for results
no_rs_p_sellers_file_name = "No_" + version + "\" + "No_" + version + "_Results" + "\No_" + version + "_Prices_Sellers.csv";
no_rs_d_sellers_file_name = "No_" + version + "\" + "No_" + version + "_Results" + "\No_" + version + "_Demand_Sellers.csv";
no_rs_cs_file_name = "No_" + version + "\" + "No_" + version + "_Results" + "\No_" +  version + "_CS.csv";

% File names for competitive equilibrium results
no_rs_p_sellers_comp_file_name = "No_" + version + "_Comp_Price.mat";
no_rs_d_sellers_comp_file_name = "No_" + version + "_Comp_Demand.mat";
no_rs_cs_comp_file_name = "No_" + version + "_Comp_CS.mat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS CS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs
cs = table2array(readtable(cs_file_name));

% Read in results competitive cs
cs_comp = load(cs_comp_file_name);
cs_comp = cs_comp.comp_cs_results;

% Make sure arrays have sime size to compute relative cs
cs_comp = repmat(cs_comp, length(omega), 1);

% Calculate ratio of cs to competitive outcome
cs_relative = [cs(:, 1:2), cs(:, 3) ./ cs_comp(:, 2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS CS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs
no_rs_cs = table2array(readtable(no_rs_cs_file_name));

% Read in results competitive cs
no_rs_cs_comp = load(no_rs_cs_comp_file_name);
no_rs_cs_comp = no_rs_cs_comp.comp_cs_results;

% Calculate ratio of cs to competitive outcome
no_rs_cs_relative = [no_rs_cs(:, 1), no_rs_cs(:, 2) ./ no_rs_cs_comp(:, 2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' prices
p_sellers = table2array(readtable(p_sellers_file_name));

% Read in results competitive sellers' prices
p_sellers_comp = load(p_sellers_comp_file_name);
p_sellers_comp = p_sellers_comp.comp_p_results;

% Mean price across sellers for each unique omega and a pair
unique_combinations = [p_sellers(:,1), p_sellers(:,2)];
[unique_vals, ~, idx] = unique(unique_combinations, 'rows');
p_sellers_mean = [unique_vals, accumarray(idx, p_sellers(:,4), [], @mean)];

% Mean competitive price across sellers
p_sellers_comp_mean = mean(p_sellers_comp(:, 2:3), 2);

% Make sure arrays have sime size to compute relative prices
p_sellers_comp_mean = repmat(p_sellers_comp_mean, length(omega), 1);

% Calculate ratio of average seller prices to competitive outcome
p_sellers_mean_relative = [p_sellers_mean(:, 1:2), p_sellers_mean(:, 3) ./ p_sellers_comp_mean];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS Prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' prices
no_rs_p_sellers = table2array(readtable(no_rs_p_sellers_file_name));

% Read in results competitive sellers' prices
no_rs_p_sellers_comp = load(no_rs_p_sellers_comp_file_name);
no_rs_p_sellers_comp = no_rs_p_sellers_comp.comp_p_results;

% Mean price across sellers for each unique a
unique_combinations = [no_rs_p_sellers(:,1)];
[unique_vals, ~, idx] = unique(unique_combinations, 'rows');
no_rs_p_sellers_mean = [unique_vals, accumarray(idx, no_rs_p_sellers(:, 3), [], @mean)];

% Mean competitive price across sellers
no_rs_p_sellers_comp_mean = mean(no_rs_p_sellers_comp(:, 2:3), 2);

% Calculate ratio of average seller prices to competitive outcome
no_rs_p_sellers_mean_relative = [no_rs_p_sellers_mean(:, 1), no_rs_p_sellers_mean(:, 2) ./ no_rs_p_sellers_comp_mean];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' demand
d_sellers = table2array(readtable(d_sellers_file_name));

% Read in results competitive sellers' demand
d_sellers_comp = load(d_sellers_comp_file_name);
d_sellers_comp = d_sellers_comp.comp_d_results;

% Total demand across sellers for each unique omega and a pair
unique_combinations = [d_sellers(:,1), d_sellers(:,2)];
[unique_vals, ~, idx] = unique(unique_combinations, 'rows');
D_sellers = [unique_vals, accumarray(idx, d_sellers(:,4))];

% Total competitive demand across sellers
D_sellers_comp = sum(d_sellers_comp(:, 2:3), 2);

% Make sure arrays have sime size to compute relative total demand
D_sellers_comp= repmat(D_sellers_comp, length(omega), 1);

% Calculate ratio of total seller demand to competitive outcome
D_sellers_relative = [D_sellers(:, 1:2), D_sellers(:, 3) ./ D_sellers_comp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS Demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' demand
no_rs_d_sellers = table2array(readtable(no_rs_d_sellers_file_name));

% Read in results competitive sellers' demand
no_rs_d_sellers_comp = load(no_rs_d_sellers_comp_file_name);
no_rs_d_sellers_comp = no_rs_d_sellers_comp.comp_d_results;

% Total demand across sellers for each unique a 
unique_combinations = [no_rs_d_sellers(:,1)];
[unique_vals, ~, idx] = unique(unique_combinations, 'rows');
no_rs_D_sellers = [unique_vals, accumarray(idx, no_rs_d_sellers(:,3))];

% Mean competitive demand across sellers
no_rs_D_sellers_comp = sum(no_rs_d_sellers_comp(:, 2:3), 2);

% Calculate ratio of total seller demand to competitive outcome
no_rs_D_sellers_relative = [no_rs_D_sellers(:, 1), no_rs_D_sellers(:, 2) ./ no_rs_D_sellers_comp];


%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%

% Create tiled layout with 2 rows and 3 columns
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% CS plot
nexttile;
hold on;
h_omega = gobjects(length(omega), 1);
for w = 1:length(omega)
    mask = cs_relative(:, 1) == omega(w);
    a_values = cs_relative(mask, 2);
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    h_omega(w) = plot(a_values, cs_relative(mask, 3), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
h_no_rs = plot(a_values, no_rs_cs_relative(:, 2), "--", "Color", colors(4), ...
        "LineWidth", 2, "DisplayName", "No RS");
hold off;
xlabel("a_{12} and a_{21}");
xtickangle(45);
xticks([1, 1.5, 2]); 
xticklabels({"1", "1.5", "2"});
ylabel("Relative Consumer Surplus");
grid on;
xlim([1, 2]); 
ylim([min(no_rs_cs_relative(:, 2) - 0.01), max(cs_relative(:, 3) + 0.01)]);

% Prices plot
nexttile;
hold on;
gobjects(length(omega), 1);
for w = 1:length(omega)
    mask = p_sellers_mean_relative(:, 1) == omega(w);
    a_values = p_sellers_mean_relative(mask, 2);
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    plot(a_values, p_sellers_mean_relative(mask, 3), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
plot(a_values, no_rs_p_sellers_mean_relative(:, 2), "--", "Color", colors(4), ...
        "LineWidth", 2, "DisplayName", "No RS");
hold off;
xlabel("a_{12} and a_{21}");
xtickangle(45);
xticks([1, 1.5, 2]); 
xticklabels({"1", "1.5", "2"});
ylabel("Relative Prices");
grid on;
xlim([1, 2]);  
ylim([min(p_sellers_mean_relative(:, 3) - 0.01), max(no_rs_p_sellers_mean_relative(:, 2) + 0.01)]);

% Demand plot
nexttile;
hold on;
for w = 1:length(omega)
    mask = D_sellers_relative(:, 1) == omega(w);
    a_values = D_sellers_relative(mask, 2);
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    plot(a_values, D_sellers_relative(mask, 3), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
plot(a_values, no_rs_D_sellers_relative(:, 2), "--", "Color", colors(4), ...
        "LineWidth", 2, "DisplayName", "No RS");
hold off;
xlabel("a_{12} and a_{21}");
xtickangle(45);
xticks([1, 1.5, 2]); 
xticklabels({"1", "1.5", "2"});
ylabel("Relative Total Output");
grid on;
xlim([1, 2]); 
ylim([min(no_rs_D_sellers_relative(:, 2) - 0.01), max(D_sellers_relative(:, 3) + 0.01)]);

% Global legend 
legend_labels = [h_omega; h_no_rs];
lgd = legend(legend_labels, "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save Prices Plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_Outcome_Variables.png"));


%%%%%%%%%%%%%%%%%%%%%
% With RS Actions
%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs with RS
a_plat = table2array(readtable(plat_actions_file_name));


%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%

% Create a tiled layout
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Omega values for legend
legend_labels = cell(1, size(length(unique(a_plat(:, 3))), 2)); 

% Plot demand learning curves
for w = 1:length(omega)
    nexttile;
    % Convert omega value to string for title
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    hold on;
    for plat_action = 1:length(unique(a_plat(:, 3)))
        mask = a_plat(:, 1) == omega(w) & a_plat(:, 3) == plat_action;
        plot(a_plat(mask, 2), a_plat(mask, 4), "--", "Color", colors(plat_action), "LineWidth", 2, ...
            "DisplayName", strcat("Action ", num2str(plat_action)));
        legend_labels{plat_action} = strcat("Action ", num2str(plat_action));
    end
    hold off;
    title(strcat("\omega = ", omega_label), "FontWeight", "normal");
    if w == 1
       ylabel("Proportion"); 
    end
    xtickangle(45);
    xlabel("a_{12} and a_{21}");
    grid on;
    ylim([0 1]); 
end

% Global legend 
lgd = legend(legend_labels, "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save CS Plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_Platform_CDF_Actions.png"));


