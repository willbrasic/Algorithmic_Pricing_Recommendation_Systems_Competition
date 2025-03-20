%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Q-learning and Sellers Q-learning theta gamma Results
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

% Type one price sensitivities
theta_1 = [0.9, 1.0, 1.1];

% Platform fee
f = 0.2;

% Overall version
version = "RS_Het_theta_gamma";

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
cs_file_name = version + "\" + version + "_Results" + "\" + version + "_CS.csv";
plat_actions_file_name = version + "\" + version + "_Results" + "\" + version + "_Actions_Plat.csv";

% File names for competitive equilibrium results
p_sellers_comp_file_name = version + "_Comp_Price.mat";
cs_comp_file_name = version + "_Comp_CS.mat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File names for results
no_rs_p_sellers_file_name = "No_" + version + "\" + "No_" + version + "_Results" + "\No_" + version + "_Prices_Sellers.csv";
no_rs_cs_file_name = "No_" + version + "\" + "No_" + version + "_Results" + "\No_" +  version + "_CS.csv";

% File names for competitive equilibrium results
no_rs_p_sellers_comp_file_name = "No_" + version + "_Comp_Price.mat";
no_rs_cs_comp_file_name = "No_" + version + "_Comp_CS.mat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' prices
p_sellers = table2array(readtable(p_sellers_file_name));

% Read in results competitive sellers' prices
p_sellers_comp = load(p_sellers_comp_file_name);
p_sellers_comp = p_sellers_comp.comp_p_results;

% Find unique (omega, gamma, theta_1) pairs across all slices to ensure consistent size
omega_values = p_sellers(:, 1);
gamma_values = p_sellers(:, 2);
theta_1_values = p_sellers(:, 3);

% Make sure arrays have sime size to compute relative prices
p_sellers_comp = repmat(p_sellers_comp(:, 3:4), length(omega), 1);

% Calculate ratio of average seller prices to competitive outcome
p_sellers_relative = [p_sellers(p_sellers(:, 4) == 1, 1:3), p_sellers(p_sellers(:, 4) == 1, 5) ./ p_sellers_comp(:, 1), p_sellers(p_sellers(:, 4) == 2, 5) ./ p_sellers_comp(:, 2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without RS Prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in results for sellers' prices
no_rs_p_sellers = table2array(readtable(no_rs_p_sellers_file_name));

% Read in results competitive sellers' prices
no_rs_p_sellers_comp = load(no_rs_p_sellers_comp_file_name);
no_rs_p_sellers_comp = no_rs_p_sellers_comp.comp_p_results;

% Calculate ratio of average seller prices to competitive outcome
no_rs_p_sellers_relative = [no_rs_p_sellers(no_rs_p_sellers(:, 3) == 1, 1:2), no_rs_p_sellers(no_rs_p_sellers(:, 3) == 1, 4)./ no_rs_p_sellers_comp(:, 3), no_rs_p_sellers(no_rs_p_sellers(:, 3) == 2, 4) ./ no_rs_p_sellers_comp(:, 4)];


%%%%%%%%%%%%%%%%%%%%%
% Prices Plot
%%%%%%%%%%%%%%%%%%%%%

% Create tiled layout with 2 rows and 3 columns
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Seller 1 prices plot
for theta = 1:length(theta_1)
    nexttile;
    hold on;
    h_omega = gobjects(length(omega), 1);
    for w = 1:length(omega)
        mask = p_sellers_relative(:, 1) == omega(w) & p_sellers_relative(:, 3) == theta_1(theta);
        gamma_values = p_sellers_relative(mask, 2);
        if omega(w) == 4/5
            omega_label = "4/5";
        else
            omega_label = num2str(omega(w));
        end
        h_omega(w) = plot(gamma_values, p_sellers_relative(mask, 4), "--", "Color", colors(w), ...
            "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
    end
    mask_no_rs = no_rs_p_sellers_relative(:, 2) == theta_1(theta);
    h_no_rs = plot(gamma_values, no_rs_p_sellers_relative(mask_no_rs, 3), "--", "Color", colors(4), ...
            "LineWidth", 2, "DisplayName", "No RS");
    hold off;
    theta_label = "0.9";
    if theta_1(theta) == 0.9
        ylabel("Relative Prices");
    elseif theta_1(theta) == 1
        theta_label = "1.0";
    else
        theta_label = theta_1(theta);
    end
    title("\theta = (" + theta_label + ", 1.0)", "FontWeight", "Normal");
    xtickangle(45);
    grid on;
    ylim([min(p_sellers_relative(:, 4) - 0.01), max(no_rs_p_sellers_relative(:, 3) + 0.01)])
end

% Seller 2 prices plot
for theta = 1:length(theta_1)
    nexttile;
    hold on;
    for w = 1:length(omega)
        mask = p_sellers_relative(:, 1) == omega(w) & p_sellers_relative(:, 3) == theta_1(theta);
        gamma_values = p_sellers_relative(mask, 2);
        if omega(w) == 4/5
            omega_label = "4/5";
        else
            omega_label = num2str(omega(w));
        end
        plot(gamma_values, p_sellers_relative(mask, 5), "--", "Color", colors(w), ...
            "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
    end
    mask_no_rs = no_rs_p_sellers_relative(:, 2) == theta_1(theta);
    plot(gamma_values, no_rs_p_sellers_relative(mask_no_rs, 4), "--", "Color", colors(4), ...
            "LineWidth", 2, "DisplayName", "No RS");
    hold off;
    if theta_1(theta) == 0.9
        ylabel("Relative Prices");
    end
    xtickangle(45);
    xlabel("\gamma");
    grid on;
    ylim([min(p_sellers_relative(:, 5) - 0.01), max(no_rs_p_sellers_relative(:, 4) + 0.01)])
end

% Global legend 
legend_labels = [h_omega; h_no_rs];
lgd = legend(legend_labels, "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save Prices Plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_Seller_Prices.png"));


%%%%%%%%%%%%%%%%%%%%%
% With RS CS
%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs with RS
cs = table2array(readtable(cs_file_name));

% Read in results competitive cs
cs_comp = load(cs_comp_file_name);
cs_comp = cs_comp.comp_cs_results;

% Find unique (omega, gamma, theta_1) pairs across all slices to ensure consistent size
omega_values = cs(:, 1);
gamma_values = cs(:, 2);
theta_1_values = cs(:, 3);

% Make sure arrays have sime size to compute relative cs
cs_comp = repmat(cs_comp(:, 3), length(omega), 1);

% Calculate ratio of cs to competitive outcome
cs_relative = [cs(:, 1:3), cs(:, 4) ./ cs_comp];


%%%%%%%%%%%%%%%%%%%%%
% Without RS CS
%%%%%%%%%%%%%%%%%%%%%

% Read in results for cs without RS
no_rs_cs = table2array(readtable(no_rs_cs_file_name));

% Read in results competitive cs
no_rs_cs_comp = load(no_rs_cs_comp_file_name);
no_rs_cs_comp = no_rs_cs_comp.comp_cs_results;

% Calculate ratio of cs to competitive outcome
no_rs_cs_relative = [no_rs_cs(:, 1:2), no_rs_cs(:, 3) ./ no_rs_cs_comp(:, 3)];


%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%

% Create tiled layout with 1 row and 3 columns
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% CS plot
for theta = 1:length(theta_1)
    nexttile;
    hold on;
    h_omega = gobjects(length(omega), 1);
    for w = 1:length(omega)
        mask = cs_relative(:, 1) == omega(w) & cs_relative(:, 3) == theta_1(theta);
        gamma_values = cs_relative(mask, 2);
        if omega(w) == 4/5
            omega_label = "4/5";
        else
            omega_label = num2str(omega(w));
        end
        h_omega(w) = plot(gamma_values, cs_relative(mask, 4), "--", "Color", colors(w), ...
            "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
    end
    mask_no_rs = no_rs_cs_relative(:, 2) == theta_1(theta);
    h_no_rs = plot(gamma_values, no_rs_cs_relative(mask_no_rs, 3), "--", "Color", colors(4), ...
            "LineWidth", 2, "DisplayName", "No RS");
    hold off;
    theta_label = "0.9";
    if theta_1(theta) == 0.9
        ylabel("Relative Consumer Surplus");
    elseif theta_1(theta) == 1
        theta_label = "1.0";
    else
        theta_label = theta_1(theta);
    end
    title("\theta = (" + theta_label + ", 1.0)", "FontWeight", "Normal");
    xtickangle(45);
    xlabel("\gamma");
    grid on;
    ylim([min(no_rs_cs_relative(:, 3) - 0.01), max(cs_relative(:, 4) + 0.01)])
end

% Global legend 
legend_labels = [h_omega; h_no_rs];
lgd = legend(legend_labels, "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save CS Plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_CS.png"));






