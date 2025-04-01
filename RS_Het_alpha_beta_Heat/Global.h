
#include <vector>
#include <cmath>

#ifndef GLOBAL_H
#define GLOBAL_H


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Number of episodes
    extern const int E;

    // Firms' model parameters
    extern const int n;                                                                 // Number of products/sellers
    extern const int m;                                                                 // Number of equally spaced price points
    extern const std::vector<std::vector<double>> a;                                    // Product quality index
    extern const int a_0;                                                               // Inverse index of aggregate demand
    extern const double mu;                                                             // Horizontal differentiation index
    extern const double exp_a_0_mu;                                                     // Exponentiate outside option term
    extern const std::vector<int> mc;                                                   // Marginal costs
    extern const int p_sellers_comb;                                                    // Number of price combinations

    // Platform's model parameters
    extern const int k;                                                                 // Number of consumer types
    extern const int plat_actions;                                                      // Number of platform actions
    extern const double f;                                                              // Platform fee
    extern const double gammma;                                                          // Share of consumers of type one
    extern const std::vector<double> theta;                                             // Price coefficient for consumers
    extern const std::vector<double> c;                                                 // Search cost for consumers
    extern const double tau;                                                            // Share of consumers who only see recommended product
    extern const std::vector<double> omega;                                             // Platform profit weight

    // Learning parameters
    extern const int q;                                                                 // Memory
    extern const double delta;                                                          // Discount factor
    extern const double alpha;                                                          // Learning rate
    extern const double beta;                                                           // Experimentation rate

    // Size of platform's state space
    extern const int S_plat_cardinality;

    // Size of sellers' state spaces
    extern const int S_sellers_cardinality;

    // Minimum and maximum values in sellers' action spaces
    extern const double A_sellers_min;
    extern const double A_sellers_max;

    // Step size for sellers' action spaces
    extern const double A_sellers_step_size;

    // How many values of alpha and beta for the heat map
    extern const int alpha_values;
    extern const int beta_values;

    // Minimum and maximum values for alpha and beta space
    extern const double alpha_min;
    extern const double alpha_max;
    extern const double beta_min;
    extern const double beta_max;

    // Step size for alpha and beta space
    extern const double alpha_step_size;
    extern const double beta_step_size;

    // Calculate argmax_a(Q(a,s)) for each s once every *convergence_check* time periods to check for convergence
    extern const int convergence_check;

    // Calculate results over last *results_check* time periods prior to convergence
    extern const int results_check;

    // Number of time periods needed for argmax_a(Q(a,s)) to be constant for each state to achieve convergence
    extern const int norm;

    // Number of time periods allowed per episode
    extern const int maxt;

    // Vector of values used to determine next sellers' states
    extern std::vector<int> pow_vec;

    // Time periods until convergence
    extern std::vector<std::vector<std::vector<std::vector<int>>>> converge_store;

    // Store platform's action
    extern std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> a_plat_store;

    // Variables to store that are averaged over the last *results_check* time periods prior to convergence
    extern std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> pi_sellers_store;  // Sellers' profits
    extern std::vector<std::vector<std::vector<std::vector<double>>>> pi_plat_store;                  // Platform's profit
    extern std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> rvn_sellers_store; // Sellers' revenues
    extern std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> d_sellers_store;   // Sellers' demand
    extern std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> p_sellers_store;   // Sellers' prices
    extern std::vector<std::vector<std::vector<std::vector<double>>>> cs_store;                       // Consumer surplus

    // Variables to average over episodes
    extern std::vector<std::vector<std::vector<double>>> converge_store_avg_e;
    extern std::vector<std::vector<std::vector<std::vector<double>>>> a_plat_store_avg_e;
    extern std::vector<std::vector<std::vector<std::vector<double>>>> pi_sellers_store_avg_e;
    extern std::vector<std::vector<std::vector<double>>> pi_plat_store_avg_e;
    extern std::vector<std::vector<std::vector<std::vector<double>>>> rvn_sellers_store_avg_e;
    extern std::vector<std::vector<std::vector<std::vector<double>>>> d_sellers_store_avg_e;
    extern std::vector<std::vector<std::vector<std::vector<double>>>> p_sellers_store_avg_e;
    extern std::vector<std::vector<std::vector<double>>> cs_store_avg_e;





} // namespace Global

#endif //GLOBAL_H
