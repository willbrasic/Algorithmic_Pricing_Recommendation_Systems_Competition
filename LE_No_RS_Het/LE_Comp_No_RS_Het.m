%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Competitive Equilibrium Solver for Baseline Specification
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Primitives   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of products/firms
n = 2;

% Number of consumers on platform
k = 2;   

% Share of consumers of type 1
gamma = 1/2 .* ones(k, 1);

% Share of consumers shown only recommended seller
tau = 0;

% Platform fee
f = 0.2; 

% Price coefficients for consumers
theta = ones(k, 1);   

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
p0 = 2.5*ones(n, 1);

% Assign pj0 to the initialized prices
pj0 = p0;

% Initialize the max distance of the two prices
norm = 1;

% Initialize the mean distance of the two prices
avgnorm = 1;

% When max distance is smaller than this number we converged
tol = 1e-8;

% "Dampening" parameter
eta = 0.5;

% Counter for number of convergence iterations
i = 0;

% While both the max and average distances are greater than tol
while norm > tol && avgnorm > tol

    % Set pj equal to the initialized prices
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
    pj0_1 = eta * pj(1) + (1 - eta) * ((mc(1) ./ (1 - f)) + (mu .* comp_d(1)) ./ ( tau .* (gamma(1) .* theta(1) .* comp_r_1_1 .* (1 - comp_r_1_1) + gamma(2) .* theta(2) .* comp_r_1_2 .* (1 - comp_r_1_2)) + ...
        (1 - tau) .* (gamma(1) .* theta(1) .* comp_s_1_1 .* (1 - comp_s_1_1) + gamma(2) .* theta(2) .* comp_s_1_2 .* (1 - comp_s_1_2) ) ) );
    pj0_2 = eta * pj(2) + (1 - eta) * ((mc(2) ./ (1 - f)) + (mu .* comp_d(2)) ./ ( tau .* (gamma(1) .* theta(1) .* comp_r_2_1 .* (1 - comp_r_2_1) + gamma(2) .* theta(2) .* comp_r_2_2 .* (1 - comp_r_2_2)) + ...
        (1 - tau) .* (gamma(1) .* theta(1) .* comp_s_2_1 .* (1 - comp_s_2_1) + gamma(2) .* theta(2) .* comp_s_2_2 .* (1 - comp_s_2_2) ) ) );
    pj0 = [pj0_1; pj0_2];

    % Obtain the distance between old prices pj and new prices pj0
    t = abs(pj0 - pj);

    % Get the max distance out of the two pricing distances
    norm = max(t);

    % Get the average distance out of the two pricing distances
    avgnorm = mean(t);

    % Increase the counter
    i = i + 1;
    
    % Print the counter and the average distance out of the two prices
    %fprintf(1,'iteration %1.0f---avgnorm %1.12f.\n',i,avgnorm);

    % If we get to 1000 loops, we probably are not convering so quit the
    % program
    if i == 1000 
        fprintf(1,'Probably not converging---quiting fixed point routine (avgnorm = %1.6f).\n',avgnorm)
        norm = tol; avgnorm = tol; 
    end
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

fprintf(1,'\n****************************************************\n');
fprintf(1,'********** OVERALL RESULTS (Competitve) *************\n');
fprintf(1,'****************************************************\n');
fprintf(1,'\nNumber of repetitions until convergence: %1.0f\n',i);
fprintf(1,'\n                                 Firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'             Tot/Avg');
fprintf(1,'         %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nProfits         %1.4f', sum(comp_pi))
fprintf(1,'    %1.4f', comp_pi)
fprintf(1,'\nRevenues        %1.4f', sum(comp_rvn))
fprintf(1,'    %1.4f', comp_rvn)
fprintf(1,'\nDemand          %1.4f', sum(comp_d))
fprintf(1,'    %1.4f', comp_d)
fprintf(1,'\nPrices          %1.4f', mean(comp_p))
fprintf(1,'    %1.4f', comp_p)
fprintf(1,'\nCS              %1.4f', comp_cs)
fprintf(1,'\n---------------------------------------------------------\n');





