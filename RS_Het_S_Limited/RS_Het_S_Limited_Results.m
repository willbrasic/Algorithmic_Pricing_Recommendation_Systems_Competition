%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RS Het Results
% William Brasic 
% The University of Arizona
% wbrasic@arizona.edu 
% williambrasic.com
% January 2024
%
% This script obtains results for RS_Het_S_Limited.m
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
color_5 = "blue";
color_6 = "red";
colors = [color_1, color_2, color_3, color_4, color_5, color_6];

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

% Version of simulation
version = "RS_Het_S_Limited";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Logit Equilibrium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logit version
logit_version = strrep(version, "_S_Limited", "");

% Add path to logit equilibrium solver scripts
addpath(strcat("LE_", logit_version));

% Solve for logit collusive equilibrium
run(strcat("LE_", logit_version, "\LE_Coll_", logit_version, ".m"));

% Solve for logit competitive equilibrium
run(strcat("LE_", logit_version, "\LE_Comp_", logit_version, ".m"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path to with RS files
addpath(strcat(version));

% File name for final results averaged over episodes
results_final_file_name = strcat(version, "\", ...
    version, "_Results_Final_", num2str(R), ".csv");

% File name for learning curve results for sellers' prices
p_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Prices_Sellers_", num2str(R), ".mat");

% File name for learning curve results for sellers' revenues
rvn_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Revenues_Sellers_", num2str(R), ".mat");

% File name for learning curve results for sellers' demand
d_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Demand_Sellers_", num2str(R), ".mat");

% File name for learning curve results for sellers' profits
pi_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Profits_Sellers_", num2str(R), ".mat");

% File name for learning curve results for platform profits
pi_plat_lc_file_name = strcat(version, "\", ...
    version, "_LC_Profits_Plat_", num2str(R), ".mat");

% File name for learning curve results for consumer surplus
cs_lc_file_name = strcat(version, "\", ...
    version, "_LC_CS_", num2str(R), ".mat");

% File name for storing learning curve results for platform optimal action selection
plat_opt_action_lc_file_name = strcat(version, "\", ...
    version, "_LC_Plat_Optimal_Action_", num2str(R), ".mat");

% File name for storing learning curve results for platform action selection CDF
plat_cdf_actions_lc_file_name = strcat(version, "\", ...
    version, "_LC_Plat_CDF_Actions_", num2str(R), ".mat");

% File name for learning curve results for platform actions
plat_actions_lc_file_name = strcat(version, "\", ...
    version, "_LC_Plat_Actions_", num2str(R), ".mat");

% File name for time steps until convergence
converge_file_name = strcat(version, "\", ...
    version, "_Converge_", num2str(R), ".mat");

% File name for platform action selections
plat_actions_file_name = strcat(version, "\", ...
    version, "_Plat_Actions_", num2str(R), ".mat");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With RS Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in final results 
results_final = readtable(results_final_file_name);

% Mean convergence            
mean_converge = results_final.Mean_Convergence;

% Mean market and firm profit
mean_market_pi = results_final.Mean_Market_Profit;
mean_firm_pi = [results_final.Mean_Firm_1_Profit, results_final.Mean_Firm_2_Profit];

% Mean market and firm revenues
mean_market_rvn = results_final.Mean_Market_Revenue;
mean_firm_rvn = [results_final.Mean_Firm_1_Revenue, results_final.Mean_Firm_2_Revenue];

% Mean market and firm demand
mean_market_d = results_final.Mean_Market_Quantity;
mean_firm_d = [results_final.Mean_Firm_1_Quantity, results_final.Mean_Firm_2_Quantity];

% Mean market and firm price
mean_market_p = results_final.Mean_Market_Price;
mean_firm_p = [results_final.Mean_Firm_1_Price, results_final.Mean_Firm_2_Price];

% Mean platform profit
mean_platform_profit = results_final.Mean_Platform_Profit;

% Mean consumer surplus
mean_cs = results_final.Mean_CS;

% Obtain results across omega values
for w = 1:length(omega)
    % Average results across *E* episodes
    fprintf(1,"\n*********************************************************\n");
    fprintf(1,"* OVERAL RESULTS                  ***********************\n");
    fprintf(1,"* (Averaged across %4.0f episodes) ***********************\n",E);
    fprintf(1,"* (omega = %1.2f) ***********************\n",omega(w));
    fprintf(1,"*********************************************************\n");
    fprintf(1,"\nAverage number of time steps until convergence: %1.0f\n", mean_converge(w));
    fprintf(1,"\n                                 Firms                 \n");
    fprintf(1,"                       ----------------------------------\n");
    fprintf(1,"             Tot/Avg");
    fprintf(1,"            %1.0f", [1:n]')
    fprintf(1,"\n---------------------------------------------------------\n");
    fprintf(1,"\nProfits            %1.4f", mean_market_pi(w))
    fprintf(1,"      %1.4f", mean_firm_pi(w, :))
    fprintf(1,"\nRevnues            %1.4f", mean_market_rvn(w))
    fprintf(1,"      %1.4f", mean_firm_rvn(w, :))
    fprintf(1,"\nDemand             %1.4f", mean_market_d(w))
    fprintf(1,"      %1.4f", mean_firm_d(w, :))
    fprintf(1,"\nPrices             %1.4f", mean_market_p(w))
    fprintf(1,"      %1.4f", mean_firm_p(w, :))
    fprintf(1,"\nPlatform Profit    %1.4f", mean_platform_profit(w))
    fprintf(1,"\nCS                 %1.4f", mean_cs(w))
    fprintf(1,"\n---------------------------------------------------------\n");
    
    % Average percentage change from competitive outcome across *E* episodes
    fprintf(1,"\n*********************************************************\n");
    fprintf(1,"* PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ************\n");
    fprintf(1,"* (Averaged across %4.0f episodes)            ************\n",E);
    fprintf(1,"* (omega = %1.2f) ***********************\n",omega(w));
    fprintf(1,"*********************************************************\n");
    fprintf(1,"\nAverage number of time steps until convergence: %1.0f\n", mean_converge(w));
    fprintf(1,"\n                                 Firms                 \n");
    fprintf(1,"                       ----------------------------------\n");
    fprintf(1,"               Tot/Avg");
    fprintf(1,"        %1.0f", [1:n]')
    fprintf(1,"\n---------------------------------------------------------\n");
    fprintf(1,"\nProfits         %1.2f%%", 100 * (mean_market_pi(w) - sum(comp_pi))/sum(comp_pi))
    fprintf(1,"    %1.2f%%", 100 * (mean_firm_pi(w, :) - comp_pi')./comp_pi')
    fprintf(1,"\nRevenues        %1.2f%%", 100 * (mean_market_rvn(w) - sum(comp_rvn)) / sum(comp_rvn))
    fprintf(1,"    %1.2f%%", 100 * (mean_firm_rvn(w, :) - comp_rvn')./comp_rvn')
    fprintf(1,"\nDemand          %1.2f%%", 100 * (mean_market_d(w) - sum(comp_d))/sum(comp_d))
    fprintf(1,"    %1.2f%%", 100 * (mean_firm_d(w, :) - comp_d')./comp_d')
    fprintf(1,"\nPrices          %1.2f%%", 100 * (mean_market_p(w) - mean(comp_p))/mean(comp_p))
    fprintf(1,"      %1.2f%%", 100 * (mean_firm_p(w, :) - comp_p')./comp_p')
    fprintf(1,"\nCS              %1.2f%%", 100*(mean_cs(w) - comp_cs)/comp_cs)
    fprintf(1,"\n---------------------------------------------------------\n");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CS, Prices, and Demand LC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% CS
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for cs w/ RS
cs_lc_data = load(cs_lc_file_name);
cs_lc = cs_lc_data.cs_lc;

% Read in learning curve data for cs w/o RS
cs_lc_no_rs_data = load(cs_lc_no_rs_file_name);
cs_lc_no_rs = cs_lc_no_rs_data.cs_lc;

% Replace zero values for in learning curve cs with last non-zero element for each episode for platform case
for w = 1:length(omega)
    for e = 1:E
        % Get the time-series slice for episode and omega value
        cs_slice = cs_lc(:, e, w);

        % Find the last non-zero value
        last_non_zero_cs = find(cs_slice ~= 0, 1, "last");

        if ~isempty(last_non_zero_cs)
            % Replace all zero entries with the last non-zero value
            for t = 1:length(cs_slice)
                if cs_slice(t) == 0
                    cs_lc(t, e, w) = cs_slice(last_non_zero_cs);
                else
                    cs_lc(t, e, w) = cs_slice(t);
                end
            end
        else
            % If no non-zero value exists, keep the slice as is
            cs_lc(:, e, w) = cs_slice;
        end
    end
end

% Replace zero values for in learning curve cs with last non-zero element for each episode for no platform case
for e = 1:E
    % Get the time-series slice for seller i, episode e
    cs_slice = cs_lc_no_rs(:, e);

    % Find the last non-zero value
    last_non_zero_cs = find(cs_slice ~= 0, 1, "last");

    if ~isempty(last_non_zero_cs)
        % Replace all zero entries with the last non-zero value
        for t = 1:length(cs_slice)
            if cs_slice(t) == 0
                cs_lc_no_rs(t, e) = cs_slice(last_non_zero_cs);
            else
                cs_lc_no_rs(t, e) = cs_slice(t);
            end
        end
    else
        % If no non-zero value exists, keep the slice as is
        cs_lc_no_rs(:, e) = cs_slice;
    end
end

% Calculate quantiles across episodes
Q1_cs = quantile(cs_lc, 0.25, 2);
Q3_cs = quantile(cs_lc, 0.75, 2);
Q1_cs_no_rs = quantile(cs_lc_no_rs, 0.25, 2);
Q3_cs_no_rs = quantile(cs_lc_no_rs, 0.75, 2);

% Calculate median cs over episodes
cs_lc_med = median(cs_lc, 2);
cs_lc_no_rs_med = median(cs_lc_no_rs, 2);


%%%%%%%%%%%%%%%%%%%%%%
% Prices
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for sellers' prices w/ RS
p_sellers_lc_data = load(p_sellers_lc_file_name);
p_sellers_lc = p_sellers_lc_data.p_sellers_lc;

% Read in learning curve data for sellers' prices w/o RS
p_sellers_lc_no_rs_data = load(p_sellers_lc_no_rs_file_name);
p_sellers_lc_no_rs = p_sellers_lc_no_rs_data.p_sellers_lc;

% Replace zero values for in learning curve price with last non-zero element for each episode for platform case
for w = 1:length(omega)
    for e = 1:E
        for i = 1:n
            % Get the time-series slice for seller i, episode e, omega value w
            p_slice = p_sellers_lc(:, i, e, w);

            % Find the last non-zero value
            last_non_zero_p = find(p_slice ~= 0, 1, "last");

            if ~isempty(last_non_zero_p)
                % Replace all zero entries with the last non-zero value
                for t = 1:length(p_slice)
                    if p_slice(t) == 0
                        p_sellers_lc(t, i, e, w) = p_slice(last_non_zero_p);
                    else
                        p_sellers_lc(t, i, e, w) = p_slice(t);
                    end
                end
            else
                % If no non-zero value exists, keep the slice as is
                p_sellers_lc(:, i, e, w) = p_slice;
            end
        end
    end
end

% Replace zero values for in learning curve price with last non-zero element for each episode for no platform case
for e = 1:E
    for i = 1:n
        % Get the time-series slice for seller i, episode e
        p_slice = p_sellers_lc_no_rs(:, i, e);

        % Find the last non-zero value
        last_non_zero_p = find(p_slice ~= 0, 1, "last");

        if ~isempty(last_non_zero_p)
            % Replace all zero entries with the last non-zero value
            for t = 1:length(p_slice)
                if p_slice(t) == 0
                    p_sellers_lc_no_rs(t, i, e) = p_slice(last_non_zero_p);
                else
                    p_sellers_lc_no_rs(t, i, e) = p_slice(t);
                end
            end
        else
            % If no non-zero value exists, keep the slice as is
            p_sellers_lc_no_rs(:, i, e) = p_slice;
        end
    end
end

% Calculate mean prices over sellers
p_sellers_lc_avg = mean(p_sellers_lc, 2);
p_sellers_lc_no_rs_avg = mean(p_sellers_lc_no_rs, 2);

% Calculate quantiles across episodes
Q1_p_sellers = quantile(p_sellers_lc_avg, 0.25, 3);
Q3_p_sellers = quantile(p_sellers_lc_avg, 0.75, 3);
Q1_p_sellers_no_rs = quantile(p_sellers_lc_no_rs_avg, 0.25, 3);
Q3_p_sellers_no_rs = quantile(p_sellers_lc_no_rs_avg, 0.75, 3);

% Calculate mean prices over sellers and median over episodes
p_sellers_lc_avg_med = median(p_sellers_lc_avg, 3);
p_sellers_lc_no_rs_avg_med = median(p_sellers_lc_no_rs_avg, 3);


%%%%%%%%%%%%%%%%%%%%%%
% Demand
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for sellers' demand w/ RS
d_sellers_lc_data = load(d_sellers_lc_file_name);
d_sellers_lc = d_sellers_lc_data.d_sellers_lc;

% Read in learning curve data for sellers' demand w/o RS
d_sellers_lc_no_rs_data = load(d_sellers_lc_no_rs_file_name);
d_sellers_lc_no_rs = d_sellers_lc_no_rs_data.d_sellers_lc;

% Replace zero values for in learning curve demand with last non-zero element for each episode for platform case
for w = 1:length(omega)
    for e = 1:E
        for i = 1:n
            % Get the time-series slice for seller i, episode e, omega value w
            d_slice = d_sellers_lc(:, i, e, w);

            % Find the last non-zero value
            last_non_zero_d = find(d_slice ~= 0, 1, "last");

            if ~isempty(last_non_zero_d)
                % Replace all zero entries with the last non-zero value
                for t = 1:length(d_slice)
                    if d_slice(t) == 0
                        d_sellers_lc(t, i, e, w) = d_slice(last_non_zero_d);
                    else
                        d_sellers_lc(t, i, e, w) = d_slice(t);
                    end
                end
            else
                % If no non-zero value exists, keep the slice as is
                d_sellers_lc(:, i, e, w) = d_slice;
            end
        end
    end
end

% Replace zero values for in learning curve demand with last non-zero element for each episode for no platform case
for e = 1:E
    for i = 1:n
        % Get the time-series slice for seller i, episode e
        d_slice = d_sellers_lc_no_rs(:, i, e);

        % Find the last non-zero value
        last_non_zero_d = find(d_slice ~= 0, 1, "last");

        if ~isempty(last_non_zero_d)
            % Replace all zero entries with the last non-zero value
            for t = 1:length(d_slice)
                if d_slice(t) == 0
                    d_sellers_lc_no_rs(t, i, e) = d_slice(last_non_zero_d);
                else
                    d_sellers_lc_no_rs(t, i, e) = d_slice(t);
                end
            end
        else
            % If no non-zero value exists, keep the slice as is
            d_sellers_lc_no_rs(:, i, e) = d_slice;
        end
    end
end

% Calculate total demand over sellers
D_sellers_lc = sum(d_sellers_lc, 2);
D_sellers_lc_no_rs = sum(d_sellers_lc_no_rs, 2);

% Calculate quantiles across episodes
Q1_D_sellers = quantile(D_sellers_lc, 0.25, 3);
Q3_D_sellers = quantile(D_sellers_lc, 0.75, 3);
Q1_D_sellers_no_rs = quantile(D_sellers_lc_no_rs, 0.25, 3);
Q3_D_sellers_no_rs = quantile(D_sellers_lc_no_rs, 0.75, 3);

% Calculate median total demand over episodes
D_sellers_lc_med = median(D_sellers_lc, 3);
D_sellers_lc_no_rs_med = median(D_sellers_lc_no_rs, 3);


%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%

% Create tiled layout with 1 row and 3 columns
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Consumer Surplus plot
nexttile;
hold on;
h_omega = gobjects(length(omega), 1);
x = 1:length(cs_lc_med);
x_fill = [x, fliplr(x)];
for w = 1:length(omega)
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    y_fill = [Q1_cs(:, w)', fliplr(Q3_cs(:, w)')];
    fill(x_fill, y_fill, colors(w), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h_omega(w) = plot(cs_lc_med(:, w), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
y_fill = [Q1_cs_no_rs', fliplr(Q3_cs_no_rs')];
    fill(x_fill, y_fill, colors(4), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h_no_rs = plot(cs_lc_no_rs_med, "--", "Color", colors(4), "LineWidth", 2, ...
    "DisplayName", "No RS");
h_bertrand = yline(comp_cs, "-", "Color", colors(5), "LineWidth", 2, ...
    "DisplayName", "Bertrand-Nash");
h_collusive = yline(coll_cs, "-", "Color", colors(6), "LineWidth", 2, ...
    "DisplayName", "Joint-Collusive");
hold off;
ylabel("Consumer Surplus");
xtickangle(45);
xticks([1, 20, 40, 60, 80, 100]); 
xticklabels({"10", "200", "400", "600", "800", "1,000"});
xlabel("Time Period \times 1,000");
grid on;
ylim([coll_cs - 0.01, max(max(cs_lc_med)) + 0.01]);

% Prices plot
nexttile;
hold on;
x = 1:length(p_sellers_lc_avg_med);
x_fill = [x, fliplr(x)];
for w = 1:length(omega)
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    y_fill = [Q1_p_sellers(:, w)', fliplr(Q3_p_sellers(:, w)')];
    fill(x_fill, y_fill, colors(w), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(p_sellers_lc_avg_med(:, w), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
y_fill = [Q1_p_sellers_no_rs', fliplr(Q3_p_sellers_no_rs')];
    fill(x_fill, y_fill, colors(4), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(p_sellers_lc_no_rs_avg_med, "--", "Color", colors(4), "LineWidth", 2, ...
    "DisplayName", "No RS");
yline(mean(comp_p), "-", "Color", colors(5), "LineWidth", 2, ...
    "DisplayName", "Bertrand-Nash");
yline(mean(coll_p), "-", "Color", colors(6), "LineWidth", 2, ...
    "DisplayName", "Joint-Collusive");
hold off;
ylabel("Prices");
xtickangle(45);
xticks([1, 20, 40, 60, 80, 100]); 
xticklabels({"10", "200", "400", "600", "800", "1,000"});
xlabel("Time Period \times 1,000");
grid on;
ylim([min(min(p_sellers_lc_avg_med)) + 0.01, mean(coll_p) + 0.01]);

% Demand plot
nexttile;
hold on;
x = 1:length(D_sellers_lc_med);
x_fill = [x, fliplr(x)];
for w = 1:length(omega)
    if omega(w) == 4/5
        omega_label = "4/5";
    else
        omega_label = num2str(omega(w));
    end
    y_fill = [Q1_D_sellers(:, w)', fliplr(Q3_D_sellers(:, w)')];
    fill(x_fill, y_fill, colors(w), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(D_sellers_lc_med(:, w), "--", "Color", colors(w), ...
        "LineWidth", 2, "DisplayName", strcat("\omega = ", omega_label));
end
y_fill = [Q1_D_sellers_no_rs', fliplr(Q3_D_sellers_no_rs')];
    fill(x_fill, y_fill, colors(4), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(D_sellers_lc_no_rs_med, "--", "Color", colors(4), "LineWidth", 2, ...
    "DisplayName", "No RS");
yline(sum(comp_d), "-", "Color", colors(5), "LineWidth", 2, ...
    "DisplayName", "Bertrand-Nash");
yline(sum(coll_d), "-", "Color", colors(6), "LineWidth", 2, ...
    "DisplayName", "Joint-Collusive");
hold off;
ylabel("Total Output");
xtickangle(45);
xticks([1, 20, 40, 60, 80, 100]); 
xticklabels({"10", "200", "400", "600", "800", "1,000"});
xlabel("Time Period \times 1,000");
grid on;
ylim([sum(coll_d) - 0.01, max(max(D_sellers_lc_med)) + 0.01]);

% Global legend 
legend_labels = [h_omega; h_no_rs; h_bertrand; h_collusive];
lgd = legend(legend_labels, "Orientation", "horizontal", "NumColumns", 2);
lgd.Layout.Tile = "south";

% Save learning curve plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_LC.png"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Platform Actions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% CDF Actions LC
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for platform actions
a_plat_cdf_lc_data = load(plat_cdf_actions_lc_file_name);
a_plat_cdf_lc = a_plat_cdf_lc_data.a_plat_cdf_lc;

% Replace zero values for in learning curve platform actions with last non-zero element for each episode for platform case
for w = 1:length(omega)
    for e = 1:E
        for plat_action = 1:size(a_plat_cdf_lc, 2)
            % Get the time-series slice for platform action plat_action, episode e, omega value w
            a_plat_cdf_slice = a_plat_cdf_lc(:, plat_action, e, w);

            % Find the last non-zero value
            last_non_zero_d = find(a_plat_cdf_slice ~= 0, 1, "last");

            if ~isempty(last_non_zero_d)
                % Replace all zero entries with the last non-zero value
                for t = 1:length(a_plat_cdf_slice)
                    if a_plat_cdf_slice(t) == 0
                        a_plat_cdf_lc(t, plat_action, e, w) = a_plat_cdf_slice(last_non_zero_d);
                    else
                        a_plat_cdf_lc(t, plat_action, e, w) = a_plat_cdf_slice(t);
                    end
                end
            else
                % If no non-zero value exists, keep the slice as is
                a_plat_cdf_lc(:, plat_action, e, w) = a_plat_cdf_slice;
            end
        end
    end
end

% Calculate mean platform actions over episodes
a_plat_cdf_lc_avg = mean(a_plat_cdf_lc, 3);


%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%

% Create a tiled layout
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Omega values for legend
legend_labels = cell(1, size(a_plat_cdf_lc, 2)); 

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
    for plat_action = 1:size(a_plat_cdf_lc, 2)
        plot(a_plat_cdf_lc_avg(:, plat_action, w), "--", "Color", colors(plat_action), "LineWidth", 2, ...
            "DisplayName", strcat("Action ", num2str(plat_action)));
        legend_labels{plat_action} = strcat("Action ", num2str(plat_action));
    end
    hold off;
    title(strcat("\omega = ", omega_label), "FontWeight", "normal");
    if w == 1
       ylabel("Proportion"); 
    end
    xtickangle(45);
    xticks([1, 20, 40, 60, 80, 100]); 
    xticklabels({"10", "200", "400", "600", "800", "1,000"}); 
    xlabel("Time Period \times 1,000");
    grid on;
    ylim([0 1]); 
end

% Global legend 
lgd = legend(legend_labels, "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save platform action learning curve
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_LC_Platform_CDF_Actions.png"));


% %%%%%%%%%%%%%%%%%%%%%%
% % Optimal Actions LC
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Read in learning curve data for platform actions
% a_plat_opt_lc_data = load(plat_opt_action_lc_file_name);
% a_plat_opt_lc = a_plat_opt_lc_data.a_plat_opt_lc;
% 
% % Replace zero values for in learning curve platform actions with last non-zero element for each episode for platform case
% for w = 1:length(omega)
%     for e = 1:E
%         % Get the time-series slice for episode e, omega value w
%         a_plat_opt_slice = a_plat_opt_lc(:, e, w);
% 
%         % Find the last non-zero value
%         last_non_zero_a_plat = find(a_plat_opt_slice ~= 0, 1, "last");
% 
%         if ~isempty(last_non_zero_a_plat)
%             % Replace all zero entries with the last non-zero value
%             for t = 1:length(a_plat_opt_slice)
%                 if a_plat_opt_slice(t) == 0
%                     a_plat_opt_lc(t, e, w) = a_plat_opt_slice(last_non_zero_a_plat);
%                 else
%                     a_plat_opt_lc(t, e, w) = a_plat_opt_slice(t);
%                 end
%             end
%         else
%             % If no non-zero value exists, keep the slice as is
%             a_plat_opt_lc(:, e, w) = a_plat_opt_slice;
%         end
%     end
% end
% 
% % Calculate mean platform actions over episodes
% a_plat_opt_lc_avg = mean(a_plat_opt_lc, 2);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% % Plot
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Plot platform optimal action learning curves
% figure;
% hold on;
% for w = 1:length(omega)
%     if omega(w) == 4/5
%         omega_label = "4/5";
%     else
%         omega_label = num2str(omega(w));
%     end
%     plot(a_plat_opt_lc_avg(:, w), "--", "Color", colors(w), "LineWidth", 2, ...
%         "DisplayName", strcat("\omega = ", omega_label));
% end
% hold off;
% ylabel("Proportion");
% xtickangle(45);
% xticks([1, 20, 40, 60, 80, 100]); 
% xticklabels({"10", "200", "400", "600", "800", "1,000"}); 
% xlabel("Time Period \times 1,000");
% legend("Orientation", "horizontal", "Location", "southoutside");
% grid on;
% 
% % Save platform action learning curve
% saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
%     version, "\", version, "_LC_Platform_Optimal_Actions.png"));


% %%%%%%%%%%%%%%%%%%%%%%
% % Platform Actions Bar
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Read in platform action selection data
% a_plat_store_Q_plat_0_data = load(plat_actions_file_name);
% a_plat_store_Q_plat_0 = a_plat_store_Q_plat_0_data.a_plat_store_Q_plat_0;
% 
% % Platform's actions
% a_plat_store = a_plat_store_Q_plat_0(1, :, :);
% 
% % Initial platform Q-matrix
% Q_plat_0 = a_plat_store_Q_plat_0(2, :, :);
% 
% % Highlight the bar that is the platform's profit maximizing action for each omega
% highlighted_bar_1 = find(Q_plat_0(1, :, 1) == max(Q_plat_0(1, :, 1)));
% highlighted_bar_2 = find(Q_plat_0(1, :, 2) == max(Q_plat_0(1, :, 2)));
% highlighted_bar_3 = find(Q_plat_0(1, :, 3) == max(Q_plat_0(1, :, 3)));
% highlighted_bar_4 = find(Q_plat_0(1, :, 4) == max(Q_plat_0(1, :, 4)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% % Plot
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Define tiled layout
% figure;
% t = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");
% 
% % Loop through omega values
% for w = 1:length(omega)
%     nexttile;
% 
%     % Extract data for the current subplot
%     action_counts = a_plat_store(1, :, w); 
%     profit_values = Q_plat_0(1, :, w);
% 
%     % Create bar plot
%     maxProfit = max(profit_values);
%     highlighted_bars = find(profit_values == maxProfit);
%     h1 = bar(1:length(action_counts), action_counts, "FaceColor", color_1, "EdgeColor", "black");
%     hold on;
%     for j = highlighted_bars
%         bar(j, action_counts(j), "FaceColor", color_2);
%     end
%     hold off;
%     ylabel("Proportion");
%     if omega(w) == 1/3
%         omega_label = "1/3";
%     elseif omega(w) == 2/3
%         omega_label = "2/3";
%     else
%         omega_label = num2str(omega(w));
%     end
%     title(strcat("\omega = ", omega_label), "FontWeight", "normal");
%     ylim([0, 1]);
%     xticks(1:length(a_plat_store));
%     xticklabels(arrayfun(@num2str, 1:length(a_plat_store), 'UniformOutput', false));
%     xlabel("Action");
%     grid on;
% end
% lgd = legend({"Non-Profit Maximizing Action", "Profit Maximizing Action"}, ...
%     "Orientation", "horizontal");
% lgd.Layout.Tile = "south";
% 
% % Save platform action selection charts
% saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
%     version, "\", version, "_Platform_Actions_Bar.png"));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Platform Profits LC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Read in learning curve data for platform profits
% pi_plat_lc_data = load(pi_plat_lc_file_name);
% pi_plat_lc = pi_plat_lc_data.pi_plat_lc;
% 
% % Replace zero values for in learning curve platform profits with last non-zero element for each episode
% for w = 1:length(omega)
%     for e = 1:E
%         % Get the time-series slice for seller i, episode e, omega value w
%         pi_plat_slice = pi_plat_lc(:, e, w);
% 
%         % Find the last non-zero value
%         last_non_zero_pi_plat = find(pi_plat_slice ~= 0, 1, "last");
% 
%         if ~isempty(last_non_zero_pi_plat)
%             % Replace all zero entries with the last non-zero value
%             for t = 1:length(pi_plat_slice)
%                 if pi_plat_slice(t) == 0
%                     pi_plat_lc(t, e, w) = pi_plat_slice(last_non_zero_pi_plat);
%                 else
%                     pi_plat_lc(t, e, w) = pi_plat_slice(t);
%                 end
%             end
%         else
%             % If no non-zero value exists, keep the slice as is
%             pi_plat_lc(:, e, w) = pi_plat_slice;
%         end
%     end
% end
% 
% % Calculate mean platform profits over episodes
% pi_plat_lc_avg = mean(pi_plat_lc, 2);
% 
% % Normalize platform profits
% for w = 1:length(omega)
%     slice = pi_plat_lc_avg(:, :, w); 
%     min_val = min(slice(:)); 
%     max_val = max(slice(:)); 
%     pi_plat_lc_avg(:, :, w) = (slice - min_val) / (max_val - min_val);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% % Plot
% %%%%%%%%%%%%%%%%%%%%%%
% 
% % Plot platform profits learning curves
% figure;
% hold on;
% for w = 1:length(omega)
%     if omega(w) == 4/5
%         omega_label = "4/5";
%     else
%         omega_label = num2str(omega(w));
%     end
%     plot(pi_plat_lc_avg(:, w), "--", "Color", colors(w), "LineWidth", 2, ...
%         "DisplayName", strcat("\omega = ", omega_label));
% end
% hold off;
% ylabel("Platform Profits");
% xtickangle(45);
% xticks([1, 20, 40, 60, 80, 100]); 
% xticklabels({"10", "200", "400", "600", "800", "1,000"}); 
% xlabel("Time Period \times 1,000");
% legend("Location", "Northeast"); 
% grid on;
% 
% % Save platform profits learning curves
% saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
%     version, "\", version, "_LC_Platform_Profits.png"));













