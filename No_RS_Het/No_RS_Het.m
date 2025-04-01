%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline No RS
% William Brasic 
% The University of Arizona
% wbrasic@arizona.edu 
% williambrasic.com
% January 2024
%
% This script simulates sellers using Q-learning pricing 
% algorithms engaging in price competition on a platform with no RS 
% with two consumer types. This is the baseline no RS specification 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Clear workspace
clear;  

% Do not show warnings
warning off all;    

% Numbers are rounded without scientific notation
format longG;      

% Reset random number generator
rng(0,"twister");

% Run (for saving results)
R = 1;  

% Number of episodes
E = 100;                    

% Version of simulation
version = "No_RS_Het";

% File name for storing final results averaged over episodes
results_final_file_name = strcat(version, "\", ...
    version, "_Results_Final_", num2str(R), ".csv");

% Delete old final results file storing results averaged over episodes
delete(results_final_file_name); 

% File name for storing learning curve results for sellers' prices
p_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Prices_Sellers_", num2str(R), ".mat");

% File name for storing learning curve results for sellers' revenues
rvn_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Revenues_Sellers_", num2str(R), ".mat");

% File name for storing learning curve results for sellers' demand
d_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Demand_Sellers_", num2str(R), ".mat");

% File name for storing learning curve results for sellers' profits
pi_sellers_lc_file_name = strcat(version, "\", ...
    version, "_LC_Profits_Sellers_", num2str(R), ".mat");

% File name for storing learning curve results for consumer surplus
cs_lc_file_name = strcat(version, "\", ...
    version, "_LC_CS_", num2str(R), ".mat");

% File name for storing time periods until convergence
converge_file_name = strcat(version, "\", ...
    version, "_Converge_", num2str(R), ".mat");

% File name for storing trigger strategy test results
rp_file_name = strcat(version, "\", ...
    version, "_RP_", num2str(R), ".mat");


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
% Simulation Primitives   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Platform's model parameters
k = 2;                              % Number of consumers types
f = 0.2;                           % Platform fee
gamma = 1/2;                        % Share of consumers of type 1
theta = ones(k, 1);                 % Price coefficients for consumers
c = 1/4 .* ones(k, 1);              % Search cost for consumers

% Firms' model parameters
n = 2;                              % Number of products/sellers
m = 15;                             % Number of equally spaced price points
a_0 = 0;                            % Inverse index of aggregate demand
mu = 1/4;                           % Horizontal differentiation index
mc = ones(n, 1);                    % Marginal costs

% Product preference matrix (columns refer to consumer types and rows refer to sellers)
a = [2, 1.9; 
    1.9, 2];    

% Learning parameters
q = 1;                              % Memory
delta = 0.95;                       % Discount factor
alpha = 0.15;                       % Learning rate
beta = 1e-5;                        % Experimentation rate

% Minimum and maximum prices for sellers' action spaces
A_seller_min = 1.0;
A_seller_max = 2.1; 

% Discretization of the sellers' action spaces (m equally spaced points)
A_sellers = linspace(A_seller_min, A_seller_max, m);

% Size of sellers' state spaces
S_seller_cardinality = m^(n * q);

% All possible combinations of actions for sellers 
A_seller_comb = A_sellers;
for i = 1:n-1
    A_seller_comb = combvec(A_seller_comb, A_sellers);
end

% sellers' prices in each state
p_sellers = A_seller_comb;

% Demand of each seller in each state
d_sellers(1, :) = gamma .* ( exp((a(1, 1) - (theta(1) .* p_sellers(1, :)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers(1, :)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, :)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(1, 2) - (theta(2) .* p_sellers(1, :)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers(1, :)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, :)) - c(2)) ./ mu) + exp(a_0 ./ mu)) );
d_sellers(2, :) = gamma .* ( exp((a(2, 1) - (theta(1) .* p_sellers(2, :)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers(1, :)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, :)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(2, 2) - (theta(2) .* p_sellers(2, :)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers(1, :)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, :)) - c(2)) ./ mu) + exp(a_0 ./ mu)) ); 

% sellers' profits in each state
pi_sellers = ((1 - f) .* p_sellers - mc) .* d_sellers;

% Consumer surplus in each state
cs = mu .* (gamma .*( log(exp((a(1, 1) - (theta(1) .* p_sellers(1, :)) - c(1)) / mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, :)) - c(1)) / mu) + exp(a_0 / mu))) + ...
    (1 - gamma) .* ( log(exp((a(1, 2) - (theta(2) .* p_sellers(1, :)) - c(2)) / mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, :)) - c(2)) / mu) + exp(a_0 / mu)) ) );

% Average profit of the states when firm *i* sets price *j*
for i = 1:n 
    for j = 1:m
        Q_sellers_0(i, j) = mean(pi_sellers(i, p_sellers(i, :) == A_sellers(j))) ./ (1 - delta);
    end
end

% Simulation header
fprintf(1,"\n**********************\n");
fprintf(1,"* SIMULATION RESULTS *\n");
fprintf(1,"**********************\n");

% Calculate argmax(Q_sellers) once every *convergence_check* repetitions to check for convergence
convergence_check = 100; 

% Number of time periods needed for argmax(Q_sellers) to be constant for convergence
norm = 100000 / convergence_check; 

% Number of time periods allowed per episode
maxt = 10000001;     

% Compute learning curve data every *lc_check* time periods
lc_check = 10000;

% Compute learning curve data until *lc_finish* time period
lc_finish = 1200000;

% Display results for episode *e* every *results_check* time periods
results_check = 100000;

% Initialize Q-matrix for sellers
Q_sellers = zeros(S_seller_cardinality, m, n, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning Curve Storage  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize matrix to store sellers' prices learning curve results
p_sellers_lc = zeros(lc_finish / lc_check, n, E);

% Initialize matrix to store sellers' revenues learning curve results
rvn_sellers_lc = zeros(lc_finish / lc_check, n, E);

% Initialize matrix to store sellers' demand learning curve results
d_sellers_lc = zeros(lc_finish / lc_check, n, E);

% Initialize matrix to store sellers' profits learning curve results
pi_sellers_lc = zeros(lc_finish / lc_check, n, E);

% Initialize matrix to store consumer surplus learning curve results
cs_lc = zeros(lc_finish / lc_check, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Storage  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize matrix to store time periods until convergence
converge_store = zeros(E, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trigger Stategy Storage  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time periods used for trigger strategy test (+ 1)
t_rp = 99;

% Price vector to store evolution before and after trigger strategy test
p_sellers_rp = zeros(n, t_rp, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Episode Loop  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start looping over episodes
for e = 1:E

    % Start timer
    tic                         

    % Clear variables at the start of each episode
    clear a_sellers p_sellers rvn_sellers pi_sellers d_sellers cs
    
    % Initialize each seller's Q-matrix
    for i = 1:n
        Q_sellers(:, :, i, e) = ones(S_seller_cardinality, 1) * Q_sellers_0(i, :); 
    end
    
    % Randomly determine intitial actions for sellers
    for i = 1:n
        a_sellers_t0(i, 1) = randperm(m, 1);
    end
 
    % Determine initial state for sellers
    temp = a_sellers_t0(1 , 1);
    for i = 2:n
        % Modify temp in some way that it lies in {1, 2, ..., S_seller_cardinality}
        temp = temp + (a_sellers_t0(i, 1) - 1) * (m^(i - 1));
    end

    % Determine inital state for sellers
    s_sellers(1) = temp;
    
    % Initalize convergence counter to 0 
    convergence_count = 0; 
    
    % Initialize time period to 1
    t = 1;

    % Initialize a counter for plotting training curves
    lc_count = 1;
    
    % Print episode *e* of E total
    fprintf("\nEpisode %1.0f of %1.0f total. Time Step: ", [e;E]);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin Time Step Loop  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % While max time periods not reached and convergence not achieved
    while ((t < maxt) && (convergence_count < norm))
        
        % Get actions for each seller
        for i = 1:n
            % Explore for sellers
            if rand(1, 1) < exp(-beta * t) 
                a_sellers(i, t) = randperm(m, 1);
            % Exploit for sellers
            else
                a_sellers(i, t) = find(Q_sellers(s_sellers(t), :, i, e) == max(Q_sellers(s_sellers(t), :, i, e)), 1); 
            end
        end
        
        % Determine subsequent seller state
        temp = a_sellers(1, t);
        for i = 2:n
            % Modify temp in some way that it lies in {1, 2, ..., S_seller_cardinality}
            temp = temp + (a_sellers(i, t) - 1) * (m^(i - 1));
        end
        s_sellers(t + 1) = temp;
             
        % Loop over sellers to get their price
        for i = 1:n
            p_sellers(i, t) = A_sellers(a_sellers(i, t)); 
        end
        
        % Demand for each seller at time period *t*
        d_sellers(1, t) = gamma .* ( exp((a(1, 1) - (theta(1) .* p_sellers(1, t)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers(1, t)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, t)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(1, 2) - (theta(2) .* p_sellers(1, t)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers(1, t)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, t)) - c(2)) ./ mu) + exp(a_0 ./ mu)) );
        d_sellers(2, t) = gamma .* ( exp((a(2, 1) - (theta(1) .* p_sellers(2, t)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers(1, t)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, t)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(2, 2) - (theta(2) .* p_sellers(2, t)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers(1, t)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, t)) - c(2)) ./ mu) + exp(a_0 ./ mu)) ); 
                  
        % Consumer surplus at time period *t*
        cs(t) = mu .* (gamma .*( log(exp((a(1, 1) - (theta(1) .* p_sellers(1, t)) - c(1)) / mu) + exp((a(2, 1) - (theta(1) .* p_sellers(2, t)) - c(1)) / mu) + exp(a_0 / mu))) + ...
            (1 - gamma) .* ( log(exp((a(1, 2) - (theta(2) .* p_sellers(1, t)) - c(2)) / mu) + exp((a(2, 2) - (theta(2) .* p_sellers(2, t)) - c(2)) / mu) + exp(a_0 / mu)) ) );

        % Revenues for each seller at time period *t*
        rvn_sellers(:, t) = (1 - f) .* p_sellers(:, t) .* d_sellers(:, t);
        
        % Profits for each seller at time period *t*
        pi_sellers(:, t) = ((1 - f) .* p_sellers(:, t) - mc) .* d_sellers(:, t);
        
        % Store data for learning curves every *lc_check* time periods
        if ((mod(t, lc_check) == 0) && (t <= lc_finish + 1))
            % Store sellers' prices averaged over last 10k time periods
            p_sellers_lc(lc_count, :, e) = mean(p_sellers(:, t-9999:t), 2);
            % Store sellers' revenues averaged over last 10k time periods
            rvn_sellers_lc(lc_count, :, e) = mean(rvn_sellers(:, t-9999:t), 2);
            % Store sellers' demand averaged over last 10k time periods
            d_sellers_lc(lc_count, :, e) = mean(d_sellers(:, t-9999:t), 2);
            % Store sellers' profits averaged over last 10k time periods
            pi_sellers_lc(lc_count, :, e) = mean(pi_sellers(:, t-9999:t), 2);
            % Store consumer surplus averaged over last 10k time periods
            cs_lc(lc_count, e) = mean(cs(t-9999:t)); 
            % Increment the learning curve counter
            lc_count = lc_count + 1;
        end
                 
        % Update Q-matrix for sellers
        for i = 1:n
            Q_sellers(s_sellers(t), a_sellers(i, t), i, e) = (1 - alpha) .* Q_sellers(s_sellers(t), a_sellers(i, t), i, e) + ...
                alpha .* (pi_sellers(i, t) + delta .* max(Q_sellers(s_sellers(t + 1), :, i, e)));
        end
        
        % Check for convergence every *convergence_check* time periods
        if mod(t, convergence_check) == 0
            if (t / convergence_check) > 1
                % Find only the indices of the maximal actions for each seller state for each seller
                [~, amax2] = max(Q_sellers(:, :, :, e), [], 2);
                % Check if prior period's optimal actions for each seller in each seller state is the same as this period
                if sum(sum(amax2 == amax1)) == n * S_seller_cardinality
                    convergence_count = convergence_count + 1;
                else
                    convergence_count = 0;
                end                
                % Set the old optimal seller action indices to the new ones
                amax1 = amax2;
            else
                % Find only the indices of the maximal actions for each seller state for each seller
                [~, amax1] = max(Q_sellers(:, :, :, e), [], 2);      
            end
        end
              
        % Display counter every *results_check* time periods
        if mod(t, results_check) == 0
            if t > results_check
                % Delete previous counter display
                for i = 0:log10(t - 1)
                    fprintf("\b"); 
                end
            end

            % Print the time
            fprintf("%d", t); 

            % Allows time for display to update
            pause(.05); 
        end
        
        % Update time period
        t = t + 1;

    % End while loop
    end
    
    % Print new line
    fprintf("\n")

    % Stop timer
    tt = toc; 

    % Average price over last *results_check* time periods of episode *e*
    p_sellers_e(:, e) = mean(p_sellers(:, t - results_check:t - 1), 2);

    % Average revenues over last *results_check* time periods of episode *e*
    rvn_sellers_e(:, e) = mean(rvn_sellers(:, t - results_check:t - 1), 2);

    % Average demand over last *results_check* time periods of episode *e*
    d_sellers_e(:, e) = mean(d_sellers(:, t - results_check:t - 1), 2);
        
    % Averge sellers' profits over last *results_check* time periods of episode *e*
    pi_sellers_e(:, e) = mean(pi_sellers(:, t - results_check:t - 1), 2);

    % Average consumer surplus over last *results_check* time periods of episode *e*
    cs_e(:, e) = mean(cs(:, t - results_check:t - 1), 2);

    % Track iterations until convergence for episode *e*
    converge(e) = t - 1;

    % Store time periods until convergence
    converge_store(e) = converge(e);

    % If t is the maximum time period (so the algorithm didn"t converge)
    if t == maxt
        fprintf(1, "Did not converge.\n");
    else
        fprintf(1, "# of time periods until convergance: %1.0f\n", t - 1);  
    end

    % Results for episode *e* of *E* averaged over last *results_check* time periods
    fprintf(1,"\nIt took about %1.0f minutes and %1.0f seconds until convergence.\n", ...
        floor(tt / 60), round(tt - floor(tt / 60) * 60))
    fprintf(1,"\nAveraged across last 100,000 time periods:\n");
    fprintf(1,"                                 Firms                 \n");
    fprintf(1,"                       ----------------------------------\n");
    fprintf(1,"                  Tot/Avg");
    fprintf(1,"           %1.0f", [1:n]')
    fprintf(1,"\n---------------------------------------------------------\n");
    fprintf(1,"\nProfits            %1.4f", sum(pi_sellers_e(:, e)))
    fprintf(1,"       %1.4f", pi_sellers_e(:, e))
    fprintf(1,"\nRevenues           %1.4f", sum(rvn_sellers_e(:, e)))
    fprintf(1,"       %1.4f", rvn_sellers_e(:, e))
    fprintf(1,"\nDemand             %1.4f", sum(d_sellers_e(:, e)))
    fprintf(1,"       %1.4f", d_sellers_e(:, e))
    fprintf(1,"\nPrices             %1.4f", mean(p_sellers_e(:, e)))
    fprintf(1,"       %1.4f", p_sellers_e(:, e))
    fprintf(1,"\nCS                 %1.4f", cs_e(:, e))
    fprintf(1,"\n---------------------------------------------------------\n");


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Trigger Strategy Test 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialization for trigger strategy test
    Q_sellers_rp = Q_sellers;
    s_sellers_rp = s_sellers;
    a_sellers_rp = a_sellers;
    pi_sellers_rp = pi_sellers;
    d_sellers_rp = d_sellers;
    cs_rp = cs;
    
    % Re-initialize *t* to convergence achieved for episode *e*
    t = converge(e);
    
    % Second index for *p_sellers_rp*
    p_rp_2 = 1 + t - converge(e);
    
    % While *t* is less than or equal to *t_rp*
    while t <= converge(e) + t_rp
        
        % Get the actions for each seller
        for i = 1:n
            % Explore for sellers
            if rand(1, 1) < exp(-beta * (t)) 
                a_sellers_rp(i, t) = randperm(m, 1); 
            % Exploit for sellers
            else 
                a_sellers_rp(i, t) = find(Q_sellers_rp(s_sellers_rp(t), :, i, e) ...
                    == max(Q_sellers_rp(s_sellers_rp(t), :, i, e)), 1);
            end

            % Force agent *i* to deviate at third time period
            if ((i == 1) && (t == converge(e) + 2))
                a_sellers_rp(i, t) = 9; % A_sellers(5) is roughly the Bertrand-Nash price
            end
        end

        % Determine subsequent seller state 
        temp = a_sellers_rp(1, t);
        for i = 2:n
            % Modify temp in some way that it lies in {1,..., S_seller_cardinality} 
            temp = temp + (a_sellers_rp(i, t) - 1) * (m^(i - 1));
        end
        s_sellers_rp(t + 1) = temp;

        % Loop over sellers to get their price
        for i = 1:n
            if ((i == 1) && (t == converge(e) + 2))
                p_sellers_rp(i, p_rp_2, e) = mean(comp_p); 
            else
                p_sellers_rp(i, p_rp_2, e) = A_sellers(a_sellers_rp(i, t)); 
            end
        end
        
        % Demand for each firm at time period *t*
        d_sellers_rp(1, t) = gamma .* ( exp((a(1, 1) - (theta(1) .* p_sellers_rp(1, p_rp_2, e)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers_rp(1, p_rp_2, e)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers_rp(2, p_rp_2, e)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(1, 2) - (theta(2) .* p_sellers_rp(1, p_rp_2, e)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers_rp(1, p_rp_2, e)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers_rp(2, p_rp_2, e)) - c(2)) ./ mu) + exp(a_0 ./ mu)) );
        d_sellers_rp(2, t) = gamma .* ( exp((a(2, 1) - (theta(1) .* p_sellers_rp(2, p_rp_2, e)) - c(1)) ./ mu) ./ (exp((a(1, 1) - (theta(1) .* p_sellers_rp(1, p_rp_2, e)) - c(1)) ./ mu) + exp((a(2, 1) - (theta(1) .* p_sellers_rp(2, p_rp_2, e)) - c(1)) ./ mu) + exp(a_0 ./ mu)) ) ...
            + (1 - gamma) .* ( exp((a(2, 2) - (theta(2) .* p_sellers_rp(2, p_rp_2, e)) - c(2)) ./ mu) ./ (exp((a(1, 2) - (theta(2) .* p_sellers_rp(1, p_rp_2, e)) - c(2)) ./ mu) + exp((a(2, 2) - (theta(2) .* p_sellers_rp(2, p_rp_2, e)) - c(2)) ./ mu) + exp(a_0 ./ mu)) ); 

        % Consumer surplus at time period *t*
        cs_rp(t) = mu .* (gamma .*( log(exp((a(1, 1) - (theta(1) .* p_sellers_rp(1, p_rp_2, e)) - c(1)) / mu) + exp((a(2, 1) - (theta(1) .* p_sellers_rp(2, p_rp_2, e)) - c(1)) / mu) + exp(a_0 / mu))) + ...
            (1 - gamma) .* ( log(exp((a(1, 2) - (theta(2) .* p_sellers_rp(1, p_rp_2, e)) - c(2)) / mu) + exp((a(2, 2) - (theta(2) .* p_sellers_rp(2, p_rp_2, e)) - c(2)) / mu) + exp(a_0 / mu)) ) );

        % Profits for each seller at time period *t*
        pi_sellers_rp(:, t) = ((1 - f) .* p_sellers_rp(:, p_rp_2, e) - mc) .* d_sellers_rp(:, t);
        
        % Update Q-matrix for sellers 
        for i = 1:n
            Q_sellers_rp(s_sellers_rp(t), a_sellers_rp(i, t), i, e) = (1 - alpha) ...
                .* Q_sellers_rp(s_sellers_rp(t), a_sellers_rp(i, t), i, e) ...
                + alpha.*(pi_sellers_rp(i, p_rp_2) + delta .* max(Q_sellers_rp(s_sellers_rp(t + 1), :, i, e)));
        end

        % Update time period
        t = t + 1;

        % Update second index of *p_sellers_rp*
        p_rp_2 = p_rp_2 + 1;

    % End RP Test
    end 

% End episode loop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average results across *E* episodes
fprintf(1,"\n*********************************************************\n");
fprintf(1,"* OVERAL RESULTS ****************************************\n");
fprintf(1,"* (Averaged across %4.0f episodes) ***********************\n",E);
fprintf(1,"*********************************************************\n");
fprintf(1,"\nAverage number of time periods until convergence: %1.0f\n", mean(converge));
fprintf(1,"\n                                 Firms                 \n");
fprintf(1,"                       ----------------------------------\n");
fprintf(1,"                  Tot/Avg");
fprintf(1,"           %1.0f", [1:n]')
fprintf(1,"\n---------------------------------------------------------\n");
fprintf(1,"\nProfits            %1.4f", mean(sum(pi_sellers_e), 2))
fprintf(1,"       %1.4f", mean(pi_sellers_e, 2))
fprintf(1,"\nRevenues           %1.4f", mean(sum(rvn_sellers_e), 2))
fprintf(1,"       %1.4f", mean(rvn_sellers_e, 2))
fprintf(1,"\nDemand             %1.4f", mean(sum(d_sellers_e), 2))
fprintf(1,"       %1.4f", mean(d_sellers_e, 2))
fprintf(1,"\nPrices             %1.4f", mean(mean(p_sellers_e, 2)))
fprintf(1,"       %1.4f", mean(p_sellers_e, 2))
fprintf(1,"\nCS                 %1.4f", mean(cs_e, 2))
fprintf(1,"\n---------------------------------------------------------\n");

% Average percentage change from competitive outcome across *E* episodes
fprintf(1,"\n*********************************************************\n");
fprintf(1,"* PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ************\n");
fprintf(1,"* (Averaged across %4.0f episodes)            ************\n",E);
fprintf(1,"*********************************************************\n");
fprintf(1,"\nAverage number of time periods until convergence: %1.0f\n", mean(converge));
fprintf(1,"\n                                 Firms                 \n");
fprintf(1,"                       ----------------------------------\n");
fprintf(1,"               Tot/Avg");
fprintf(1,"        %1.0f", [1:n]')
fprintf(1,"\n---------------------------------------------------------\n");
fprintf(1,"\nProfits         %1.2f%%", 100 * (mean(sum(pi_sellers_e), 2) - sum(comp_pi))/sum(comp_pi))
fprintf(1,"    %1.2f%%", 100 * (mean(pi_sellers_e, 2) - comp_pi) ./ comp_pi)
fprintf(1,"\nRevenues        %1.2f%%", 100 * (mean(sum(rvn_sellers_e), 2) - sum(comp_rvn)) / sum(comp_rvn))
fprintf(1,"     %1.2f%%", 100 * (mean(rvn_sellers_e, 2) - comp_rvn) ./ comp_rvn)
fprintf(1,"\nDemand          %1.2f%%", 100 * (mean(sum(d_sellers_e), 2) - sum(comp_d))/sum(comp_d))
fprintf(1,"    %1.2f%%", 100 * (mean(d_sellers_e, 2) - comp_d) ./ comp_d)
fprintf(1,"\nPrices           %1.2f%%", 100 * (mean(sum(p_sellers_e), 2) - sum(comp_p))/sum(comp_p))
fprintf(1,"     %1.2f%%", 100 * (mean(p_sellers_e, 2) - comp_p) ./ comp_p)

fprintf(1,"\nCS              %1.2f%%", 100 * (mean(cs_e, 2) - comp_cs)/comp_cs)
fprintf(1,"\n---------------------------------------------------------\n");  

% Start new final results file storing results averaged over *E* episodes
results_final = fopen(results_final_file_name, "at");  

% If file opened successfully, write file headers
if results_final ~= -1
    fprintf(results_final, "Mean_Convergence\t");
    fprintf(results_final, "Mean_Market_Profit\t");
    for i = 1:n
        fprintf(results_final, strcat("Mean_Firm_", num2str(i), "_Profit\t"));
    end
    fprintf(results_final, "Mean_Market_Revenue\t");
    for i = 1:n
        fprintf(results_final, strcat("Mean_Firm_", num2str(i), "_Revenue\t"));
    end
    fprintf(results_final, "Mean_Market_Quantity\t");
    for i = 1:n
        fprintf(results_final, strcat("Mean_Firm_", num2str(i), "_Quantity\t"));
    end
    fprintf(results_final, "Mean_Market_Price\t");
    for i = 1:n
        fprintf(results_final, strcat("Mean_Firm_", num2str(i), "_Price\t"));
    end
    fprintf(results_final, "Mean_CS\t");
    fprintf(results_final, "\n");
    fclose(results_final);
end

% Open file to write final results to
results_final = fopen(results_final_file_name, "at"); 

% If file opened successfully, write final results averaged over *E* episodes
if results_final ~= -1 
    fprintf(results_final, "%1.0f\t", mean(converge));
    fprintf(results_final, "%1.14f\t", mean(sum(pi_sellers_e), 2));
    fprintf(results_final, "%1.14f\t", mean(pi_sellers_e, 2));
    fprintf(results_final, "%1.14f\t", mean(sum(rvn_sellers_e), 2));
    fprintf(results_final, "%1.14f\t", mean(rvn_sellers_e, 2));
    fprintf(results_final, "%1.14f\t", mean(sum(d_sellers_e), 2));
    fprintf(results_final, "%1.14f\t", mean(d_sellers_e, 2)); 
    fprintf(results_final, "%1.14f\t", mean(mean(p_sellers_e, 2)));
    fprintf(results_final, "%1.14f\t", mean(p_sellers_e, 2));
    fprintf(results_final, "%1.14f\t", mean(cs_e, 2));
    fprintf(results_final, "\n");
    fclose(results_final); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning Curves 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save learning curve results
save(p_sellers_lc_file_name, "p_sellers_lc");
save(rvn_sellers_lc_file_name, "rvn_sellers_lc");
save(d_sellers_lc_file_name, "d_sellers_lc");
save(pi_sellers_lc_file_name, "pi_sellers_lc");
save(cs_lc_file_name, "cs_lc");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save convergence results
save(converge_file_name, "converge_store");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trigger Strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save RP results
save(rp_file_name, "p_sellers_rp");





