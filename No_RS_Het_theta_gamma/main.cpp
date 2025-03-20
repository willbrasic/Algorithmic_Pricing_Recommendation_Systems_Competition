///////////////////////////////////////////////////////////////////////////
// Sellers Q-learning Varying (gamma, theta_1)
//
// William Brasic
// The University of Arizona
// wbrasic@arizona.edu
// williambrasic.com
// January 2024
//
// This script simulates sellers using Q-learning pricing
// algorithms engaging in price competition with two
// consumer types. Price sensitivities and
// share of consumers allowed to vary.
///////////////////////////////////////////////////////////////////////////


//////////////////////////////
// Header Files
//////////////////////////////

#include "Functions.h"
#include "Global.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace Global;


//////////////////////////////
// Main Function
//////////////////////////////

int main()
{
    ////////////////////////////////////
    // Preliminaries
    ////////////////////////////////////

    // Set seed
    std::mt19937 gen(0);

    // Create discrete uniform distribution over {0, 1, ..., m - 1} to draw seller's action from
    std::uniform_int_distribution<> dist_1(0, m - 1);

    // Create continuous uniform distribution over the unit interval for epsilon-greedy exploration
    std::uniform_real_distribution<> dist_2(0, 1);

    // *Fixed* surpresses scientific notation and *setprecision* rounds to decimal places
    std::cout << std::fixed << std::setprecision(6);


    ////////////////////////////////////
    // main Initialization
    ////////////////////////////////////

    // Create sellers' action spaces
    std::vector<double> A_sellers = populate_A_sellers(m, A_sellers_min, A_sellers_step_size);

    // All possible combination of sellers' prices
    std::vector<std::vector<double>> p_sellers = populate_p_sellers(A_sellers, n, m);

    // Populate gamma vector
    std::vector<double> gamma_vector = populate_gamma(gamma_min, gamma_step_size, gamma_values);

    // Populate theta vector
    std::vector<std::vector<double>> theta_vector = populate_theta(theta_1_min, theta_1_step_size, theta_1_values, k);

    // Get exploration number for all possible time steps
    std::vector<double> minus_beta_times_t(maxt, 0.0);
    for (int t{ 0 }; t < maxt; ++t)
        minus_beta_times_t[t] = exp(-beta * t);

    // Values used to determine subsequent state
    for (int i{ 0 }; i < n; ++i)
        pow_vec[i] = static_cast<int>(pow(m, i));

    
    ////////////////////////////////////
    // Begin Simulation
    ////////////////////////////////////

    // Simulation header
    std::cout << "\n***********************\n";
    std::cout << "* Simulation Results *\n";
    std::cout << "***********************\n";


    ////////////////////////////////////
    // gamma Loop
    ////////////////////////////////////

    // Loop over values for gamma
    for (int g{ 0 }; g < gamma_vector.size(); ++g)
    {
        // Current value of gamma
        double gamma_current = gamma_vector[g];

        for(int theta_1{ 0 }; theta_1 < theta_vector.size(); ++theta_1)
        {

            // Current theta vector
            std::vector<double> theta(k, 0.0);
            for(int j{ 0 }; j < k; ++j)
                theta[j] = theta_vector[theta_1][j];

            // Populate initial values to obtain initial Q-matrices for sellers and the platform
            std::vector<std::vector<double>> d_sellers = populate_d_sellers(p_sellers, a, theta, c, gamma_current, mu, a_0, n, S_sellers_cardinality);
            std::vector<std::vector<double>> pi_sellers = populate_pi_sellers(d_sellers, p_sellers, mc, S_sellers_cardinality, n, f);
            std::vector<double> cs = populate_cs(p_sellers, a, theta, c, mu, gamma_current, a_0, S_sellers_cardinality, n, k);
            std::vector<std::vector<double>> Q_sellers_0 = populate_Q_sellers_0(pi_sellers, p_sellers, A_sellers, delta, S_sellers_cardinality, n, m);


            ////////////////////////////////////
            // Episode Loop
            ////////////////////////////////////

            // Start time measurement for episode
            auto start_time = std::chrono::steady_clock::now();

            // Loop over episodes
            for (int e{ 0 }; e < E; ++e)
            {
                // Start time measurement for episode
                auto start_time = std::chrono::steady_clock::now();

                // Print the current episode
                std::cout << "\nEpisode " << e + 1 << " of " << E << "\n";

                // Print the current value of gamma
                std::cout << "gamma " << g + 1 << " of " << gamma_vector.size() << " is: " << gamma_current << "\n";

                // Print the current value of theta_1
                std::cout << "theta_1 " << theta_1 + 1 << " of " << theta_vector.size() << " is: " << theta[0] << "\n";

                // Initialize convergence counter to 0
                int convergence_counter{ 0 };

                // Initialize time period to 0
                int t{ 0 };

                // Vectors to store each seller's Q-value maximizing action (used to test for convergence)
                std::vector<std::vector<int>> argmax_1(S_sellers_cardinality, std::vector<int>(n, 0));
                std::vector<std::vector<int>> argmax_2(S_sellers_cardinality, std::vector<int>(n, 0));

                // Populate sellers' Q-matrices
                std::vector<std::vector<std::vector<double>>> Q_sellers = populate_Q_sellers(Q_sellers_0, n, m, S_sellers_cardinality);

                // Vectors to store variables at time period *t*
                std::vector<int> s_sellers_t;                                           // Seller's state
                std::vector<std::vector<double>> pi_sellers_t;                          // Sellers' profits
                std::vector<std::vector<double>> p_sellers_t;                           // Sellers' prices
                std::vector<double> cs_t;                                               // Consumer surplus

                // Reserve memory in each vector up to *maxt* (increases speed)
                s_sellers_t.reserve(maxt);
                pi_sellers_t.reserve(maxt);
                p_sellers_t.reserve(maxt);
                cs_t.reserve(maxt);

                // Draw each seller's initial action to obtain the initial seller's state
                int a_sellers_init_t_0 = dist_1(gen);
                int a_sellers_init_t_1 = dist_1(gen);

                // Use preliminary sellers' actions to get the initial seller's state
                int s_sellers = a_sellers_init_t_0 * pow_vec[0] + a_sellers_init_t_1 * pow_vec[1];
                s_sellers_t.emplace_back(s_sellers);


                ////////////////////////////////////
                // While Loop
                ////////////////////////////////////

                // Start while loop
                while ((t < maxt) && (convergence_counter < norm))
                {
                    // Extract each seller's current state and the platform's current state at time period *t*
                    int s_sellers_current = s_sellers_t[t];

                    // Sellers' actions at time period *t*
                    std::vector<int> a_sellers(n, 0);
                    for (int i{ 0 }; i < n; ++i)
                    {
                        // Exploration condition for getting current action
                        bool explore_condition = dist_2(gen) < minus_beta_times_t[t];

                        // Store the action for seller *i* at time period *t*
                        int action_t_i;

                        // Explore at time period *t*
                        if (explore_condition)
                            action_t_i = dist_1(gen);
                        // Exploit at time period *t*
                        else
                        {
                            double max_Q_sellers = std::numeric_limits<double>::lowest();
                            for (int j{ 0 }; j < m; ++j)
                            {
                                if (Q_sellers[s_sellers_current][j][i] > max_Q_sellers)
                                {
                                    max_Q_sellers = Q_sellers[s_sellers_current][j][i];
                                    action_t_i = j;
                                }
                            }
                        }
                        a_sellers[i] = action_t_i;
                    }

                    // Seller's subsequent state
                    int a_sellers_t_0 = a_sellers[0];
                    int a_sellers_t_1 = a_sellers[1];
                    int s_sellers_next = a_sellers_t_0 * pow_vec[0] + a_sellers_t_1 * pow_vec[1];
                    s_sellers_t.emplace_back(s_sellers_next);

                    // Sellers' prices at time period *t*
                    double p_sellers_t_0 = A_sellers[a_sellers_t_0];
                    double p_sellers_t_1 = A_sellers[a_sellers_t_1];
                    p_sellers_t.emplace_back(std::vector<double>{p_sellers_t_0, p_sellers_t_1});

                    // Sellers' demand at time period *t*
                    double d_sellers_t_0, d_sellers_t_1;
                    d_sellers_t_0 =
                            gamma_current * (exp((a[0][0] - theta[0] * p_sellers_t_0 - c[0]) / mu) /
                                             (exp((a[0][0] - theta[0] * p_sellers_t_0 - c[0]) / mu) +
                                              exp((a[1][0] - theta[0] * p_sellers_t_1 - c[0]) / mu) + exp_a_0_mu)) +
                            (1 - gamma_current) * (exp((a[0][1] - theta[1] * p_sellers_t_0 - c[1]) / mu) /
                                                   (exp((a[0][1] - theta[1] * p_sellers_t_0 - c[1]) / mu) +
                                                    exp((a[1][1] - theta[1] * p_sellers_t_1 - c[1]) / mu) + exp_a_0_mu));
                    d_sellers_t_1 =
                            gamma_current * (exp((a[1][0] - theta[0] * p_sellers_t_1 - c[0]) / mu) /
                                             (exp((a[0][0] - theta[0] * p_sellers_t_0 - c[0]) / mu) +
                                              exp((a[1][0] - theta[0] * p_sellers_t_1 - c[0]) / mu) + exp_a_0_mu)) +
                            (1 - gamma_current) * (exp((a[1][1] - theta[1] * p_sellers_t_1 - c[1]) / mu) /
                                                   (exp((a[0][1] - theta[1] * p_sellers_t_0 - c[1]) / mu) +
                                                    exp((a[1][1] - theta[1] * p_sellers_t_1 - c[1]) / mu) + exp_a_0_mu));

                    // Consumer surplus at time period *t*
                    double exp_sum_gamma = exp((a[0][0] - theta[0] * p_sellers_t_0 - c[0]) / mu) + exp((a[1][0] - theta[0] * p_sellers_t_1 - c[0]) / mu) + exp_a_0_mu;
                    double exp_sum_one_minus_gamma = exp((a[0][1] - theta[1] * p_sellers_t_0 - c[1]) / mu) + exp((a[1][1] - theta[1] * p_sellers_t_1 - c[1]) / mu) + exp_a_0_mu;
                    double cs_current = mu * (gamma_current * log(exp_sum_gamma) + (1 - gamma_current) * log(exp_sum_one_minus_gamma));
                    cs_t.emplace_back(cs_current);

                    // Sellers' profits at time period *t*
                    double pi_sellers_t_0 = ( (1 - f) * p_sellers_t_0 - mc[0] ) * d_sellers_t_0;
                    double pi_sellers_t_1 = ( (1 - f) * p_sellers_t_1 - mc[1] ) * d_sellers_t_1;
                    pi_sellers_t.emplace_back(std::vector<double>{pi_sellers_t_0, pi_sellers_t_1});

                    // Each seller's Q-learning update (can do j += 2 since m - 1 is even)
                    for (int i{ 0 }; i < n; ++i)
                    {
                        double max_Q_sellers = std::numeric_limits<double>::lowest();
                        for (int j{ 0 }; j + 1 < m; j += 2)
                        {
                            double q_1 = Q_sellers[s_sellers_next][j][i];
                            double q_2 = Q_sellers[s_sellers_next][j + 1][i];
                            double local_max = std::max(q_1, q_2);
                            if (local_max > max_Q_sellers)
                                max_Q_sellers = local_max;
                        }
                        Q_sellers[s_sellers_current][a_sellers[i]][i] = (1 - alpha) * Q_sellers[s_sellers_current][a_sellers[i]][i]
                                                                        + alpha * (pi_sellers_t[t][i] + delta * max_Q_sellers);
                    }

                    // Check for convergence every *convergence_check* time steps
                    if (t % convergence_check == 0)
                    {
                        // If we have already checked for convergence once before, compare old argmax_a(Q(s,a)) to current argmax_a(Q(s,a))
                        if (t > convergence_check)
                        {
                            // Get actions that give maximum Q-values for each agent for each state
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_sellers_cardinality; ++s)
                                {
                                    int argmax{ 0 };
                                    double max_Q_sellers = std::numeric_limits<double>::lowest();
                                    for (int j{ 0 }; j + 1 < m; j+=2)
                                    {
                                        if (Q_sellers[s][j][i] > max_Q_sellers)
                                        {
                                            max_Q_sellers = Q_sellers[s][j][i];
                                            argmax = j;
                                        }
                                        if (Q_sellers[s][j + 1][i] > max_Q_sellers)
                                        {
                                            max_Q_sellers = Q_sellers[s][j + 1][i];
                                            argmax = j + 1;
                                        }
                                    }
                                    argmax_2[s][i] = argmax;
                                }
                            }
                            // Calculate how many times previous argmax calculation *convergence_check* time steps ago is the same as new argmax calculation
                            int sum_argmax{ 0 };
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_sellers_cardinality; ++s)
                                {
                                    if (argmax_2[s][i] == argmax_1[s][i])
                                        ++sum_argmax;
                                }
                            }
                            // If sellers' optimal actions have not changed, increase the convergence counter. Otherwise, reset it
                            if (sum_argmax == n * S_sellers_cardinality)
                                ++convergence_counter;
                            else
                                convergence_counter = 0;

                            // Set the old argmaxes to the new ones for next convergence check
                            argmax_1 = argmax_2;
                        }
                        // Calculate initial sellers' optimal actions at time period *convergence_check*
                        else
                        {
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_sellers_cardinality; ++s)
                                {
                                    int argmax{ 0 };
                                    double max_Q_sellers = std::numeric_limits<double>::lowest();
                                    for (int j{ 0 }; j < m; ++j)
                                    {
                                        if (Q_sellers[s][j][i] > max_Q_sellers)
                                        {
                                            max_Q_sellers = Q_sellers[s][j][i];
                                            argmax = j;
                                        }
                                    }
                                    argmax_1[s][i] = argmax;
                                }
                            }
                        }
                    }

                    // Increment t
                    ++t;

                // End while loop
                }

                // Average results over the last *results_check* time periods prior to convergence for episode *e*
                std::vector<double> pi_sellers_sum(n, 0.0);
                std::vector<double> p_sellers_sum(n, 0.0);
                double cs_sum = 0.0;
                for (int j = t - results_check; j < t; ++j)
                {
                    for (int i{ 0 }; i < n; ++i)
                    {
                        pi_sellers_sum[i] += pi_sellers_t[j][i];
                        p_sellers_sum[i] += p_sellers_t[j][i];

                    }
                    cs_sum += cs_t[j];
                }
                for (int i{ 0 }; i < n; ++i)
                {
                    pi_sellers_store[i][e][g][theta_1] = pi_sellers_sum[i] / results_check;
                    p_sellers_store[i][e][g][theta_1] = p_sellers_sum[i] / results_check;
                }
                cs_store[e][g][theta_1] = cs_sum / results_check;

                // End time measurement for episode loop
                auto end_time = std::chrono::steady_clock::now();

                // Calculate the duration of episode loop
                auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

                // Print the duration of the episode loop
                std::cout << "Seconds to completion of episode " << e + 1 << ": " << duration_seconds << std::endl;

            // End episode loop
            }

            // Average results over episodes
            std::vector<double> pi_sellers_sum(n, 0.0);
            std::vector<double> p_sellers_sum(n, 0.0);
            double cs_sum = 0.0;
            for (int e{ 0 }; e < E; ++e)
            {
                for (int i{ 0 }; i < n; ++i)
                {
                    pi_sellers_sum[i] += pi_sellers_store[i][e][g][theta_1];
                    p_sellers_sum[i] += p_sellers_store[i][e][g][theta_1];

                }
                cs_sum += cs_store[e][g][theta_1];
            }
            for (int i{ 0 }; i < n; ++i)
            {
                pi_sellers_store_avg_e[i][g][theta_1] = pi_sellers_sum[i] / E;
                p_sellers_store_avg_e[i][g][theta_1] = p_sellers_sum[i] / E;
            }
            cs_store_avg_e[g][theta_1] = cs_sum / E;

        // End theta_1 loop
        }

    // End gamma loop
    }

    // If running on HPC, set to true
    bool HPC{ true };

    // Declare file streams
    std::ofstream pi_sellers_file, p_sellers_file, cs_file;

    // Set correct file path prefix
    std::string path_prefix = HPC ? "./" : "../";
    std::string version = "No_RS_Het_theta_gamma";

    // Open files
    pi_sellers_file.open(path_prefix + version + "_Results/" + version + "_Profits_Sellers.csv");
    p_sellers_file.open(path_prefix + version + "_Results/" + version + "_Prices_Sellers.csv");
    cs_file.open(path_prefix + version + "_Results/" + version + "_CS.csv");

    // Writing headers
    pi_sellers_file << "gamma,theta_1,firm,profit\n";
    p_sellers_file << "gamma,theta_1,firm,price\n";
    cs_file << "gamma,theta_1,cs\n";

    for (int g{ 0 }; g < gamma_vector.size(); ++g)
    {
        for(int theta_1{ 0 }; theta_1 < theta_vector.size(); ++theta_1)
        {
            for (int i{ 0 }; i < n; ++i)
            {
                pi_sellers_file << gamma_vector[g] << "," << theta_vector[theta_1][0] << "," << i + 1 << "," << pi_sellers_store_avg_e[i][g][theta_1] << "\n";
                p_sellers_file << gamma_vector[g] << "," << theta_vector[theta_1][0] << "," << i + 1 << "," << p_sellers_store_avg_e[i][g][theta_1] << "\n";
            }
            cs_file << gamma_vector[g] << "," << theta_vector[theta_1][0] << "," << cs_store_avg_e[g][theta_1] << "\n";
        }
    }

    return 0;
}


