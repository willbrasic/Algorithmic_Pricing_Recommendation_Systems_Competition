#include "Global.h"


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Number of episodes
    const int E{ 100 };

    // Firms' model parameters
    const int n{ 2 };                                                                       // Number of products/sellers
    const int m{ 15 };                                                                      // Number of equally spaced price points
    const int a_0{ 0 };                                                                     // Inverse index of aggregate demand
    const double mu{ 1.0 / 4.0 };                                                           // Horizontal differentiation index
    const double exp_a_0_mu{ exp(a_0 / mu) };                                            // Exponentiate outside option term
    const std::vector<int> mc(n, 1);                                                  // Marginal costs
    const int k{ 2 };                                                                       // Number of consumer types
    const std::vector<double> f{ 0.2, 0.0 };                                                  // Platform fee
    const double gammma{ 1.0 / 2.0 };
    const std::vector<double> theta(k, 1.0);
    const std::vector<double> c(k, 1.0 / 4.0);

    // Learning parameters
    const int q{ 1 };                                                                       // Memory
    const double delta{ 0.95 };                                                             // Discount factor
    const double alpha{ 0.15 };                                                             // Learning rate
    const double beta{ 1e-5 };                                                              // Experimentation rate

    // Size of sellers' state spaces
    const int S_sellers_cardinality = static_cast<int>(pow(m, n * q));
    
    // Minimum and maximum values in sellers' action spaces
    const double A_sellers_min{ 1.0 };
    const double A_sellers_max{ 2.1 };

    // Step size for sellers' action spaces
    const double A_sellers_step_size = (A_sellers_max - A_sellers_min) / (m - 1);

    // How many values of a
    const int a_values { 100 };

    // Minimum and maximum values for c
    const double a_min { 1.0 };
    const double a_max { 2.0 };

    // Step size for a
    const double a_step_size = (a_max - a_min) / (a_values - 1);

    // Calculate argmax_a(Q(a,s)) for each s once every *convergence_check* time periods to check for convergence
    const int convergence_check{ 100 };

    // Calculate results over last *results_check* time steps prior to convergence
    const int results_check{ 100000 };

    // Number of time periods needed for argmax_a(Q(a,s)) to be constant for each state to achieve convergence
    const int norm{ 100000 / convergence_check };

    // Number of time periods allowed per episode
    const int maxt{ 10000000 };

    // Vector of values used to determine next sellers' states
    std::vector<int> pow_vec(n, 0);

    // Time periods until convergence
    std::vector<std::vector<int>> converge_store(
            E,
            std::vector<int>(a_values, 0)
    );

    // Sellers' profits (averaged over the last *results_check* time periods)
    std::vector<std::vector<std::vector<double>>> pi_sellers_store(
            n,
            std::vector<std::vector<double>>(
                    E,
                    std::vector<double>(a_values, 0.0)
            )
    );

    // Sellers' revenues (averaged)
    std::vector<std::vector<std::vector<double>>> rvn_sellers_store(
            n,
            std::vector<std::vector<double>>(
                    E,
                    std::vector<double>(a_values, 0.0)
            )
    );

    // Sellers' demand (averaged)
    std::vector<std::vector<std::vector<double>>> d_sellers_store(
            n,
            std::vector<std::vector<double>>(
                    E,
                    std::vector<double>(a_values, 0.0)
            )
    );

    // Sellers' prices (averaged)
    std::vector<std::vector<std::vector<double>>> p_sellers_store(
            n,
            std::vector<std::vector<double>>(
                    E,
                    std::vector<double>(a_values, 0.0)
            )
    );

    // Consumer surplus (averaged)
    std::vector<std::vector<double>> cs_store(
            E,
            std::vector<double>(a_values, 0)
    );

    // Time periods until convergence averaged over episodes
    std::vector<double> converge_store_avg_e(
            a_values, 0
    );

    // Sellers' profits averaged over episodes
    std::vector<std::vector<double>> pi_sellers_store_avg_e(
            n,
            std::vector<double>(a_values, 0.0)
    );

    // Sellers' revenues averaged over episodes
    std::vector<std::vector<double>> rvn_sellers_store_avg_e(
            n,
            std::vector<double>(a_values, 0.0)
    );

    // Sellers' demand averaged over episodes
    std::vector<std::vector<double>> d_sellers_store_avg_e(
            n,
            std::vector<double>(a_values, 0.0)
    );

    // Sellers' prices averaged over episodes
    std::vector<std::vector<double>> p_sellers_store_avg_e(
            n,
            std::vector<double>(a_values, 0.0)
    );

    // Consumer surplus averaged over episodes
    std::vector<double> cs_store_avg_e(
            a_values, 0
    );





} //Global


