%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Competitive Equilibrium Solver for Varying (theta, gamma)
% William Brasic
% The University of Arizona
% wbrasic@arizona.edu
% williambrasic.com
% January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not show warnings
warning off all;   

% Numbers are rounded
format longG; 

% Specify version
version = "RS_Het_theta_gamma";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Primitives   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of products/firms
n = 2;

% Number of consumers on platform
k = 2;   

% Share of consumers of type 1
gamma = linspace(0, 1, 100);
num_gamma = length(gamma);

% Share of consumers shown only recommended seller
tau = 0;

% Platform fee
f = 0.2;

% Price coefficients for consumers
theta = [0.9, 1.0; 
        1.0, 1.0; 
        1.1, 1.0];
num_theta = size(theta, 1);

% Search cost for consumers
c = 1/4 .* ones(k, 1);

% Product preference matrix (columns refer to consumer types and rows refer to sellers)
a = [2, 1.9; 
    1.9, 2];

% Marginal costs for each firm
mc = ones(n, 1);

% Horizontal differentiation index (Scale parameter of T1EV distribution)
mu = 1/4;

% Inverse index of aggregate demand (Mean utility of outside option)
a0 = 0;

% Initial prices for each firm
p0 = 2.5 .* ones(n, 1);

% When max distance is smaller than this number we converged
tol = 1e-8;

% "Dampening" parameter
eta = 0.5;

% Generate unique (gamma, theta_1) pairs
[gamma_grid, theta_1_grid] = meshgrid(gamma, theta(:,1));
gamma_theta_pairs = [gamma_grid(:), theta_1_grid(:)];

% Initialize result arrays with first two columns for gamma and theta_1
comp_p_results = [gamma_theta_pairs, zeros(num_gamma * num_theta, n)];
comp_d_results = [gamma_theta_pairs, zeros(num_gamma * num_theta, n)];
comp_rvn_results = [gamma_theta_pairs, zeros(num_gamma * num_theta, n)];
comp_pi_results = [gamma_theta_pairs, zeros(num_gamma * num_theta, n)];
comp_cs_results = [gamma_theta_pairs, zeros(num_gamma * num_theta, 1)]; % 1 column for CS

% Counter for storing results
index = 1;

% Iterate over each gamma-theta pair
for g = 1:num_gamma
    for t = 1:num_theta

        % Current values of gamma and theta_1
        gamma_curr = gamma(g);
        theta_curr = theta(t, :);
        
        % Initialize price iteration variables
        pj0 = p0;
        norm = 1;
        avgnorm = 1;
        i = 0;

        % Fixed point iteration for equilibrium prices
        while norm > tol && avgnorm > tol
            pj = pj0;

            % Demand from recommendation for each seller
            comp_r_1_1 = exp((a(1, 1) - (theta_curr(1) .* pj(1))) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* pj(1))) ./ mu) + exp(a0 ./ mu)) ;
            comp_r_1_2 = 0;
            comp_r_2_1 = 0;
            comp_r_2_2 = exp((a(2, 2) - (theta_curr(2) .* pj(2))) ./ mu) ./ (exp((a(2, 2) - (theta_curr(2) .* pj(2))) ./ mu) + exp(a0 ./ mu)) ;
        
            % Demand from search for each seller from each consumer type
            comp_s_1_1 = exp((a(1, 1) - (theta_curr(1) .* pj(1)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* pj(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta_curr(1) .* pj(1)) - c(1)) ./ mu) + exp(a0 ./ mu));
            comp_s_1_2 = exp((a(1, 2) - (theta_curr(2) .* pj(1)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta_curr(2) .* pj(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta_curr(2) .* pj(1)) - c(2)) ./ mu) + exp(a0 ./ mu));
            comp_s_2_1 = exp((a(2, 1) - (theta_curr(1) .* pj(2)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* pj(2)) - c(1)) ./ mu) + exp((a(2, 1) - (theta_curr(1) .* pj(2)) - c(1)) ./ mu) + exp(a0 ./ mu));
            comp_s_2_2 = exp((a(2, 2) - (theta_curr(2) .* pj(2)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta_curr(2) .* pj(2)) - c(2)) ./ mu) + exp((a(2, 2) - (theta_curr(2) .* pj(2)) - c(2)) ./ mu) + exp(a0 ./ mu));
        
            % Total demand
            comp_d = [tau .* (gamma_curr .* comp_r_1_1 + (1 - gamma_curr) .* comp_r_1_2) + (1 - tau) .* (gamma_curr .* comp_s_1_1 + (1 - gamma_curr) .* comp_s_1_2); 
                tau .* (gamma_curr .* comp_r_2_1 + (1 - gamma_curr) .* comp_r_2_2) + (1 - tau) .* (gamma_curr .* comp_s_2_1 + (1 - gamma_curr) .* comp_s_2_2)];
        
            % Compute prices (unique symmetric eq. price) with eta parameter
            pj0_1 = eta * pj(1) + (1 - eta) * ((mc(1) ./ (1 - f)) + (mu .* comp_d(1)) ./ ( tau .* (gamma_curr .* theta_curr(1) .* comp_r_1_1 .* (1 - comp_r_1_1) + (1 - gamma_curr) .* theta_curr(2) .* comp_r_1_2 .* (1 - comp_r_1_2)) + ...
                (1 - tau) .* (gamma_curr .* theta_curr(1) .* comp_s_1_1 .* (1 - comp_s_1_1) + (1 - gamma_curr) .* theta_curr(2) .* comp_s_1_2 .* (1 - comp_s_1_2) ) ) );
            pj0_2 = eta * pj(2) + (1 - eta) * ((mc(2) ./ (1 - f)) + (mu .* comp_d(2)) ./ ( tau .* (gamma_curr .* theta_curr(1) .* comp_r_2_1 .* (1 - comp_r_2_1) + (1 - gamma_curr) .* theta_curr(2) .* comp_r_2_2 .* (1 - comp_r_2_2)) + ...
                (1 - tau) .* (gamma_curr .* theta_curr(1) .* comp_s_2_1 .* (1 - comp_s_2_1) + (1 - gamma_curr) .* theta_curr(2) .* comp_s_2_2 .* (1 - comp_s_2_2) ) ) );
            pj0 = [pj0_1; pj0_2];

            % Obtain the distance between old prices pj and new prices pj0
            t = abs(pj0 - pj);
        
            % Get the max distance out of the two pricing distances
            norm = max(t);
        
            % Get the average distance out of the two pricing distances
            avgnorm = mean(t);
        
            % Increase the counter
            i = i + 1;

            % If convergence fails
            if i == 1000
                fprintf("Not converging for gamma=%.3f, theta_1=%.2f\n", gamma_curr, theta_curr(1));
                norm = tol; avgnorm = tol;
            end
        end
        % Competitive prices across the firms 
        comp_p = pj0;
        
        % Competitive demand across the firms
        comp_d(1) = tau .* (gamma_curr .* ( exp((a(1, 1) - (theta_curr(1) .* comp_p(1))) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* comp_p(1))) ./ mu) + exp(a0 ./ mu)) ) ) ...
            + (1 - tau) .* (gamma_curr .* (exp((a(1, 1) - (theta_curr(1) .* comp_p(1)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta_curr(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) ...
                + (1 - gamma_curr) .* (exp((a(1, 2) - (theta_curr(2) .* comp_p(1)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta_curr(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta_curr(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu))) );
        comp_d(2) = tau .* ((1 - gamma_curr) .* ( exp((a(2, 2) - (theta_curr(2) .* comp_p(2))) ./ mu) ./ (exp((a(2, 2) - (theta_curr(2) .* comp_p(2))) ./ mu) + exp(a0 ./ mu)) ) ) ...
                + (1 - tau) .* (gamma_curr .* (exp((a(2, 1) - (theta_curr(1) .* comp_p(2)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta_curr(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta_curr(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) ...
                    + (1 - gamma_curr) .* (exp((a(2, 2) - (theta_curr(2) .* comp_p(2)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta_curr(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta_curr(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu))) ); 
        comp_d = [comp_d(1); comp_d(2)];
        
        % Competitive revenue across the firms
        comp_rvn = (1 - f) .* comp_p .* comp_d;
        
        % Competitive average profits across the firms
        comp_pi = ((1 - f) .* comp_p - mc) .* comp_d;
        
        % Competitive consumer surplus
        comp_cs = (tau .* ((mu ./ theta_curr(1)) .* gamma_curr .* (log(exp((a(1, 1) - (theta_curr(1) .* comp_p(1))) ./ mu) + exp(a0 ./ mu)) ) + ...
                    (mu ./ theta_curr(2)) .* (1 - gamma_curr) .* (log(exp((a(2, 2) - (theta_curr(2) .* comp_p(2))) ./ mu) + exp(a0 ./ mu))) ) ...
                    + (1 - tau) .* ((mu ./ theta_curr(1)) .* gamma_curr .* (log(exp((a(1, 1) - (theta_curr(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta_curr(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) + ...
                        (mu ./ theta_curr(2)) .* (1 - gamma_curr) .* (log(exp((a(1, 2) - (theta_curr(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta_curr(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu)) ) ) );

        % Store equilibrium results
        comp_p_results(index, 3:end) = comp_p;
        comp_d_results(index, 3:end) = comp_d;
        comp_rvn_results(index, 3:end) = comp_rvn;
        comp_pi_results(index, 3:end) = comp_pi;
        comp_cs_results(index, 3) = comp_cs;

        % Increment index
        index = index + 1;
    end
end

% Save results to a .mat file
save(version + "/" + version + "_Results/" + "No_" + version + "_Comp_Price.mat", "comp_p_results");
save(version + "/" + version + "_Results/" + "No_" + version + "_Comp_Demand.mat", "comp_d_results");
save(version + "/" + version + "_Results/" + "No_" + version + "_Comp_Revenue.mat", "comp_rvn_results");
save(version + "/" + version + "_Results/" + "No_" + version + "_Comp_Profit.mat", "comp_pi_results");
save(version + "/" + version + "_Results/" + "No_" + version + "_Comp_CS.mat", "comp_cs_results");



