%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No RS Het Results
% William Brasic 
% The University of Arizona
% wbrasic@arizona.edu 
% williambrasic.com
% January 2024
%
% This script obtains results for No_RS_Het.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clear;  

% Colors for plots
color_1 = "blue";
color_2 = "red";
color_3 = "magenta";
colors = [color_1, color_2, color_3];

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

% Version of simulation
version = "No_RS_Het";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Equilibrium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path to logit equilibrium solver scripts
addpath(strcat("LE_", version));

% Solve for logit collusive equilibrium
run(strcat("LE_", version, "\LE_Coll_", version, ".m"));

% Solve for logit competitive equilibrium
run(strcat("LE_", version, "\LE_Comp_", version, ".m"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without Platform Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path files
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

% File name for learning curve results for consumer surplus
cs_lc_file_name = strcat(version, "\", ...
    version, "_LC_CS_", num2str(R), ".mat");

% File name for time steps until convergence
converge_file_name = strcat(version, "\", ...
    version, "_Converge_", num2str(R), ".mat");

% File name for trigger strategy results
rp_file_name = strcat(version, "\", ...
    version, "_RP_", num2str(R), ".mat");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without Platform Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in final results 
results_final = readtable(results_final_file_name);

% Mean convergence            
mean_converge = results_final.Mean_Convergence;

% Mean market and firm profits
mean_market_pi = results_final.Mean_Market_Profit;
mean_firm_pi = [results_final.Mean_Firm_1_Profit; results_final.Mean_Firm_2_Profit];

% Mean market and firm revenues
mean_market_rvn = results_final.Mean_Market_Revenue;
mean_firm_rvn = [results_final.Mean_Firm_1_Revenue; results_final.Mean_Firm_2_Revenue];

% Mean market and firm demand
mean_market_d = results_final.Mean_Market_Quantity;
mean_firm_d = [results_final.Mean_Firm_1_Quantity; results_final.Mean_Firm_2_Quantity];

% Mean market and firm prices
mean_market_p = results_final.Mean_Market_Price;
mean_firm_p = [results_final.Mean_Firm_1_Price; results_final.Mean_Firm_2_Price];

% Mean consumer surplus
mean_cs = results_final.Mean_CS;

% Average results across *E* episodes
fprintf(1,"\n*********************************************************\n");
fprintf(1,"* OVERAL RESULTS                  ***********************\n");
fprintf(1,"* (Averaged across %4.0f episodes) ***********************\n",E);
fprintf(1,"*********************************************************\n");
fprintf(1,"\nAverage number of time steps until convergence: %1.0f\n", mean_converge);
fprintf(1,"\n                                 Firms                 \n");
fprintf(1,"                       ----------------------------------\n");
fprintf(1,"             Tot/Avg");
fprintf(1,"            %1.0f", [1:n]')
fprintf(1,"\n---------------------------------------------------------\n");
fprintf(1,"\nProfits            %1.4f", mean_market_pi)
fprintf(1,"      %1.4f", mean_firm_pi)
fprintf(1,"\nRevnues            %1.4f", mean_market_rvn)
fprintf(1,"      %1.4f", mean_firm_rvn)
fprintf(1,"\nDemand             %1.4f", mean_market_d)
fprintf(1,"      %1.4f", mean_firm_d)
fprintf(1,"\nPrices             %1.4f", mean_market_p)
fprintf(1,"      %1.4f", mean_firm_p)
fprintf(1,"\nCS                 %1.4f", mean_cs)
fprintf(1,"\n---------------------------------------------------------\n");

% Average percentage change from competitive outcome across *E* episodes
fprintf(1,"\n*********************************************************\n");
fprintf(1,"* PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ************\n");
fprintf(1,"* (Averaged across %4.0f episodes)            ************\n",E);
fprintf(1,"*********************************************************\n");
fprintf(1,"\nAverage number of time steps until convergence: %1.0f\n", mean_converge);
fprintf(1,"\n                                 Firms                 \n");
fprintf(1,"                       ----------------------------------\n");
fprintf(1,"               Tot/Avg");
fprintf(1,"        %1.0f", [1:n]')
fprintf(1,"\n---------------------------------------------------------\n");
fprintf(1,"\nProfits         %1.2f%%", 100 * (mean_market_pi - sum(comp_pi))/sum(comp_pi))
fprintf(1,"    %1.2f%%", 100 * (mean_firm_pi - comp_pi')/comp_pi')
fprintf(1,"\nRevenues        %1.2f%%", 100 * (mean_market_rvn - sum(comp_rvn)) / sum(comp_rvn))
fprintf(1,"     %1.2f%%", 100 * (mean_firm_rvn - comp_rvn')/comp_rvn')
fprintf(1,"\nDemand          %1.2f%%", 100 * (mean_market_d - sum(comp_d))/sum(comp_d))
fprintf(1,"    %1.2f%%", 100 * (mean_firm_d - comp_d')/comp_d')
fprintf(1,"\nPrices          %1.2f%%", 100 * (mean_market_p - mean(comp_p))/mean(comp_p))
fprintf(1,"    %1.2f%%", 100 * (mean_firm_p - comp_p')/comp_p')
fprintf(1,"\nCS              %1.2f%%", 100 * (mean_cs - comp_cs)/comp_cs)
fprintf(1,"\n---------------------------------------------------------\n");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CS, Prices, and Demand LC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% CS
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for cs 
cs_lc_data = load(cs_lc_file_name);
cs_lc = cs_lc_data.cs_lc;

% Replace zero values for in learning curve cs with last non-zero element for each episode 
for e = 1:E
    % Get the time-series slice for each episode 
    cs_slice = cs_lc(:, e);

    % Find the last non-zero value
    last_non_zero_cs = find(cs_slice ~= 0, 1, "last");
    if ~isempty(last_non_zero_cs)
        % Replace all zero entries with the last non-zero value
        for t = 1:length(cs_slice)
            if cs_slice(t) == 0
                cs_lc(t, e) = cs_slice(last_non_zero_cs);
            else
                cs_lc(t, e) = cs_slice(t);
            end
        end
    else
        % If no non-zero value exists, keep the slice as is
        cs_lc(:, e) = cs_slice;
    end
end

% Calculate quantiles across episodes
Q1_cs = quantile(cs_lc, 0.25, 2)';
Q3_cs = quantile(cs_lc, 0.75, 2)';

% Calculate median cs over episodes
cs_lc_med = median(cs_lc, 2);


%%%%%%%%%%%%%%%%%%%%%%
% Prices
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for sellers' prices
p_sellers_lc_data = load(p_sellers_lc_file_name);
p_sellers_lc = p_sellers_lc_data.p_sellers_lc;

% Replace zero values for in learning curve price with last non-zero element for each episode
for e = 1:E
    for i = 1:n
        % Get the time-series slice for seller i, episode e, omega value w
        p_slice = p_sellers_lc(:, i, e);

        % Find the last non-zero value
        last_non_zero_p = find(p_slice ~= 0, 1, "last");
        if ~isempty(last_non_zero_p)
            % Replace all zero entries with the last non-zero value
            for t = 1:length(p_slice)
                if p_slice(t) == 0
                    p_sellers_lc(t, i, e) = p_slice(last_non_zero_p);
                else
                    p_sellers_lc(t, i, e) = p_slice(t);
                end
            end
        else
            % If no non-zero value exists, keep the slice as is
            p_sellers_lc(:, i, e) = p_slice;
        end
    end
end

% Calculate mean prices over sellers
p_sellers_lc_avg = mean(p_sellers_lc, 2);

% Calculate quantiles across episodes
Q1_p_sellers = quantile(p_sellers_lc_avg, 0.25, 3)';
Q3_p_sellers = quantile(p_sellers_lc_avg, 0.75, 3)';

% Calculate mean prices over sellers and median over episodes
p_sellers_lc_avg_med = median(p_sellers_lc_avg, 3);


%%%%%%%%%%%%%%%%%%%%%%
% Demand
%%%%%%%%%%%%%%%%%%%%%%

% Read in learning curve data for sellers' demand
d_sellers_lc_data = load(d_sellers_lc_file_name);
d_sellers_lc = d_sellers_lc_data.d_sellers_lc;

% Replace zero values for in learning curve demand with last non-zero element for each episode
for e = 1:E
    for i = 1:n
        % Get the time-series slice for seller and episode
        d_slice = d_sellers_lc(:, i, e);

        % Find the last non-zero value
        last_non_zero_d = find(d_slice ~= 0, 1, "last");
        if ~isempty(last_non_zero_d)
            % Replace all zero entries with the last non-zero value
            for t = 1:length(d_slice)
                if d_slice(t) == 0
                    d_sellers_lc(t, i, e) = d_slice(last_non_zero_d);
                else
                    d_sellers_lc(t, i, e) = d_slice(t);
                end
            end
        else
            % If no non-zero value exists, keep the slice as is
            d_sellers_lc(:, i, e) = d_slice;
        end
    end
end

% Calculate total demand over sellers
D_sellers_lc = sum(d_sellers_lc, 2);

% Calculate quantiles across episodes
Q1_D_sellers = quantile(D_sellers_lc, 0.25, 3)';
Q3_D_sellers = quantile(D_sellers_lc, 0.75, 3)';

% Calculate median demand over episodes
D_sellers_lc_med = median(D_sellers_lc, 3);


%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%

% Create tiled layout with 1 row and 3 columns
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Consumer Surplus plot
nexttile;
hold on;
h1 = yline(comp_cs, "-", "Color", colors(1), "LineWidth", 2);
h2 = yline(coll_cs, "-", "Color", colors(2), "LineWidth", 2);
x = 1:length(cs_lc_med);
x_fill = [x, fliplr(x)];
y_fill = [Q1_cs, fliplr(Q3_cs)];
fill(x_fill, y_fill, colors(3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(cs_lc_med, "--", "Color", colors(3), "LineWidth", 2);
hold off;
ylabel("Consumer Surplus");
xtickangle(45);
xticks([1, 30, 60, 90, 120]); 
xticklabels({"1", "300", "600", "900", "1,200"});
xlabel("Time Period \times 1,000");
grid on;
ylim([coll_cs - 0.1, comp_cs + 0.1]);

% Prices plot
nexttile;
hold on;
yline(mean(comp_p), "-", "Color", colors(1), "LineWidth", 2);
yline(mean(coll_p), "-", "Color", colors(2), "LineWidth", 2);
x = 1:length(p_sellers_lc_avg_med);
x_fill = [x, fliplr(x)];
y_fill = [Q1_p_sellers, fliplr(Q3_p_sellers)];
fill(x_fill, y_fill, colors(3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(p_sellers_lc_avg_med, "--", "Color", colors(3), "LineWidth", 2);
hold off;
ylabel("Prices");
xtickangle(45);
xticks([1, 30, 60, 90, 120]); 
xticklabels({"1", "300", "600", "900", "1,200"});
xlabel("Time Period \times 1,000");
grid on;
ylim([mean(comp_p) - 0.1, mean(coll_p) + 0.1]);

% Demand plot
nexttile;
hold on;
x = 1:length(D_sellers_lc_med);
x_fill = [x, fliplr(x)];
y_fill = [Q1_D_sellers, fliplr(Q3_D_sellers)];
fill(x_fill, y_fill, colors(3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
yline(sum(comp_d), "-", "Color", colors(1), "LineWidth", 2);
yline(sum(coll_d), "-", "Color", colors(2), "LineWidth", 2);
plot(D_sellers_lc_med, "--", "Color", colors(3), "LineWidth", 2);
hold off;
ylabel("Total Output");
xtickangle(45);
xticks([1, 30, 60, 90, 120]); 
xticklabels({"1", "300", "600", "900", "1,200"});
xlabel("Time Period \times 1,000");
grid on;
ylim([sum(coll_d) - 0.1, sum(comp_d) + 0.1]);

% Create a global legend at the bottom using only the Bertrand-Nash and Joint-Collusive handles
lgd = legend([h1, h2], "Bertrand-Nash", "Joint-Collusive", "Orientation", "horizontal");
lgd.Layout.Tile = "south";

% Save learning curve plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_LC.png"));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trigger Strategy Test Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in trigger strategy test data with platform
p_sellers_rp_data = load(rp_file_name);
p_sellers_rp = p_sellers_rp_data.p_sellers_rp;

% Calculate mean trigger strategy test prices over episodes
p_sellers_rp_avg = mean(p_sellers_rp, 3);

% Plot trigger strategy test
figure;
hold on;
plot(p_sellers_rp_avg(1, :, 1), "--o", "Color", "magenta", "DisplayName", "Deviating Seller");
plot(p_sellers_rp_avg(2, :, 1), "--o", "Color", "cyan", "DisplayName", "Non-Deviating Seller");
yline(mean(comp_p), "-", "Color", colors(1), "LineWidth", 2, "DisplayName", "Bertrand-Nash");
yline(mean(coll_p), "-", "Color", colors(2), "LineWidth", 2, "DisplayName", "Joint-Collusive");
hold off;
ylabel("Price");
ylim([mean(comp_p) - 0.1, mean(coll_p) + 0.1]);
xlabel("Time Period");
xtickangle(45);
grid on;

% Create a horizontal legend centered at the bottom of the figure
legend("show", "Orientation", "horizontal", "Location", "southoutside", "NumColumns", 2);

% Save trigger strategy test plot
saveas(gcf, strcat("C:\Users\wbras\OneDrive\Documents\Desktop\UA\3rd_Year_Paper\3rd_Year_Paper\3rd_Year_Paper_Pictures\", ...
    version, "\", version, "_Trigger_Strategy.png"));
