%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Competitive Equilibrium Solver for Varying c
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
version = "No_RS_Het_c";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Primitives   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of products/firms
n = 2;

% Number of consumers on platform
k = 2;   

% Vary consumer disutility factor c.
% Define 100 equally spaced points from 0 to 1/2.
c_vector = linspace(0, 1/2, 100);
num_c = length(c_vector);

% Share of consumers shown only recommended seller
tau = 0;

% Consumer type shares (fixed)
gamma = 1/2 .* ones(k, 1);

% Platform fee
f = [0.2; 0.2];

% Price coefficients for consumers
theta = ones(k, 1);

% Search cost for consumers
c = 0 .* ones(k, 1);

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

% Initialize result arrays with first column for c values
comp_p_results = [c_vector', zeros(num_c, n)];
comp_d_results = [c_vector', zeros(num_c, n)];
comp_rvn_results = [c_vector', zeros(num_c, n)];
comp_pi_results = [c_vector', zeros(num_c, n)];
comp_cs_results = [c_vector', zeros(num_c, 1)]; 

% Counter for storing results
index = 1;

% Iterate over each c value
for i = 1:num_c

    % Update c
    c = [c_vector(i); c_vector(i)];
    
    % Initialize price iteration variables
    pj0 = p0;
    norm = 1;
    avgnorm = 1;
    i_iter = 0;

    % Fixed point iteration for equilibrium prices
    while norm > tol && avgnorm > tol
        pj = pj0;

        % Demand from recommendation for each seller
        comp_r_1_1 = exp((a(1, 1) - (theta(1) .* pj(1))) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* pj(1))) ./ mu) + exp(a0 ./ mu)) ;
        comp_r_1_2 = 0;
        comp_r_2_1 = 0;
        comp_r_2_2 = exp((a(2, 2) - (theta(2) .* pj(2))) ./ mu) ./ (exp((a(2, 2) - (theta(2) .* pj(2))) ./ mu) + exp(a0 ./ mu)) ;
        
        % Demand from search for each seller from each consumer type
        comp_s_1_1 = exp((a(1, 1) - (theta(1) .* pj(1)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* pj(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* pj(1)) - c(1)) ./ mu) + exp(a0 ./ mu));
        comp_s_1_2 = exp((a(1, 2) - (theta(2) .* pj(1)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* pj(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* pj(1)) - c(2)) ./ mu) + exp(a0 ./ mu));
        comp_s_2_1 = exp((a(2, 1) - (theta(1) .* pj(2)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* pj(2)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* pj(2)) - c(1)) ./ mu) + exp(a0 ./ mu));
        comp_s_2_2 = exp((a(2, 2) - (theta(2) .* pj(2)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* pj(2)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* pj(2)) - c(2)) ./ mu) + exp(a0 ./ mu));
    
        % Total demand
        comp_d = [tau .* (gamma(1) .* comp_r_1_1 + gamma(2) .* comp_r_1_2) + (1 - tau) .* (gamma(1) .* comp_s_1_1 + gamma(2) .* comp_s_1_2); 
            tau .* (gamma(1) .* comp_r_2_1 + gamma(2) .* comp_r_2_2) + (1 - tau) .* (gamma(1) .* comp_s_2_1 + gamma(2) .* comp_s_2_2)];
    
        % Compute prices (unique symmetric eq. price) with eta parameter
        pj0_1 = eta * pj(1) + (1 - eta) * ((mc(1) ./ (1 - f(1))) + (mu .* comp_d(1)) ./ ( tau .* (gamma(1) .* theta(1) .* comp_r_1_1 .* (1 - comp_r_1_1) + gamma(2) .* theta(2) .* comp_r_1_2 .* (1 - comp_r_1_2)) + ...
            (1 - tau) .* (gamma(1) .* theta(1) .* comp_s_1_1 .* (1 - comp_s_1_1) + gamma(2) .* theta(2) .* comp_s_1_2 .* (1 - comp_s_1_2) ) ) );
        pj0_2 = eta * pj(2) + (1 - eta) * ((mc(2) ./ (1 - f(2))) + (mu .* comp_d(2)) ./ ( tau .* (gamma(1) .* theta(1) .* comp_r_2_1 .* (1 - comp_r_2_1) + gamma(2) .* theta(2) .* comp_r_2_2 .* (1 - comp_r_2_2)) + ...
            (1 - tau) .* (gamma(1) .* theta(1) .* comp_s_2_1 .* (1 - comp_s_2_1) + gamma(2) .* theta(2) .* comp_s_2_2 .* (1 - comp_s_2_2) ) ) );
        pj0 = [pj0_1; pj0_2];

        % Obtain the distance between old prices pj and new prices pj0
        t = abs(pj0 - pj);
        
        % Get the max distance out of the two pricing distances
        norm = max(t);
        
        % Get the average distance out of the two pricing distances
        avgnorm = mean(t);
        
        % Increase the counter
        i_iter = i_iter + 1;

        % If convergence fails
        if i_iter == 1000
            fprintf("Not converging for a2=%.3f\n", a2_vector(i));
            norm = tol; avgnorm = tol;
        end
    end
    % Competitive prices across the firms 
    comp_p = pj0;
    
    % Competitive demand across the firms
    comp_d(1) = tau .* (gamma(1) .* ( exp((a(1, 1) - (theta(1) .* comp_p(1))) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* comp_p(1))) ./ mu) + exp(a0 ./ mu)) ) ) ...
        + (1 - tau) .* (gamma(1) .* (exp((a(1, 1) - (theta(1) .* comp_p(1)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) ...
            + gamma(2) .* (exp((a(1, 2) - (theta(2) .* comp_p(1)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu))) );
    comp_d(2) = tau .* (gamma(2) .* ( exp((a(2, 2) - (theta(2) .* comp_p(2))) ./ mu) ./ (exp((a(2, 2) - (theta(2) .* comp_p(2))) ./ mu) + exp(a0 ./ mu)) ) ) ...
            + (1 - tau) .* (gamma(1) .* (exp((a(2, 1) - (theta(1) .* comp_p(2)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) ...
                + gamma(2) .* (exp((a(2, 2) - (theta(2) .* comp_p(2)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu))) ); 
    comp_d = [comp_d(1); comp_d(2)];
    
    % Competitive revenue across the firms
    comp_rvn = (1 - f) .* comp_p .* comp_d;
    
    % Competitive average profits across the firms
    comp_pi = ((1 - f) .* comp_p - mc) .* comp_d;
    
    % Competitive consumer surplus
    comp_cs = mu .* (tau .* (gamma(1) .* (log(exp((a(1, 1) - (theta(1) .* comp_p(1))) ./ mu) + exp(a0 ./ mu)) ) + ...
                gamma(2) .* (log(exp((a(2, 2) - (theta(2) .* comp_p(2))) ./ mu) + exp(a0 ./ mu))) ) ...
                + (1 - tau) .* (gamma(1) .* (log(exp((a(1, 1) - (theta(1) .* comp_p(1)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* comp_p(2)) - c(1)) ./ mu) + exp(a0 ./ mu)) ) + ...
                    gamma(2) .* (log(exp((a(1, 2) - (theta(2) .* comp_p(1)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* comp_p(2)) - c(2)) ./ mu) + exp(a0 ./ mu)) ) ) );

    % Store equilibrium results
    comp_p_results(index, 2:end) = comp_p;
    comp_d_results(index, 2:end) = comp_d;
    comp_rvn_results(index, 2:end) = comp_rvn;
    comp_pi_results(index, 2:end) = comp_pi;
    comp_cs_results(index, 2) = comp_cs;

    % Increment index
    index = index + 1;
end

% Save results to a .mat file
save("LE_" + version + "/" + version + "_Comp_Price.mat", "comp_p_results");
save("LE_" + version + "/" + version + "_Comp_Demand.mat", "comp_d_results");
save("LE_" + version + "/" + version + "_Comp_Revenue.mat", "comp_rvn_results");
save("LE_" + version + "/" + version + "_Comp_Profit.mat", "comp_pi_results");
save("LE_" + version + "/" + version + "_Comp_CS.mat", "comp_cs_results");








