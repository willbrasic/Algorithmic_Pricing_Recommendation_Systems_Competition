#ifndef FUNCTIONS_H
#define FUNCTIONS_H


//////////////////////////////
// Header Files
//////////////////////////////


#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>


//////////////////////////////
// Function Prototypes
//////////////////////////////


// Function to populate sellers' action spaces
std::vector<double> populate_A_sellers(
    int m,
    double A_sellers_min,
    double A_sellers_step_size
);

// Function to determine prices for each state for each firm
std::vector<std::vector<double>> populate_p_sellers(
    const std::vector<double>& A_sellers,
    int n,
    int m,
    int plat_actions
);

// Function to determine demand for each state for each firm given the platform's action
std::vector<std::vector<double>> populate_d_sellers(
    const std::vector<std::vector<double>>& p_sellers,
    const std::vector<std::vector<double>>& a,
    const std::vector<double>& theta,
    const std::vector<double>& c,
    double gammma,
    double mu,
    double tau,
    int a_0,
    int n,
    int S_sellers_cardinality,
    int p_sellers_comb
);

// Function to determine profit for each state for each firm given the platform's action
std::vector<std::vector<double>> populate_pi_sellers(
    const std::vector<std::vector<double>>& d_sellers,
    const std::vector<std::vector<double>>& p_sellers,
    const std::vector<int>& mc,
    const std::vector<double>& f,
    int S_sellers_cardinality,
    int n
);

// Function to determine consumer surplus for each state
std::vector<double> populate_cs(
    const std::vector<std::vector<double>>& p_sellers,
    const std::vector<std::vector<double>>& a,
    const std::vector<double>& theta,
    const std::vector<double>& c,
    double mu,
    double gammma,
    double tau,
    int a_0,
    int S_sellers_cardinality,
    int n,
    int p_sellers_comb
);

// Function to determine platform profits for each state given the platform's action
std::vector<double> populate_pi_plat(
    const std::vector<std::vector<double>>& d_sellers,
    const std::vector<std::vector<double>>& p_sellers,
    const std::vector<double>& cs,
    double omega,
    const std::vector<double>& f,
    int S_plat_cardinality,
    int n
);

// Function to determine average profits for each firm when setting each of the *m* possible prices for each platform's action
std::vector<std::vector<double>> populate_Q_sellers_0(
    const std::vector<std::vector<double>>& pi_sellers,
    const std::vector<std::vector<double>>& p_sellers,
    const std::vector<double>& A_sellers,
    double delta,
    int S_sellers_cardinality,
    int n,
    int m
);

// Function to determine average profits for the platform given their action
std::vector<double> populate_Q_plat_0(
        const std::vector<double>& pi_plat,
        double delta,
        int plat_actions,
        int p_sellers_comb
);

// Function to populate Q-matrix for each firm
std::vector<std::vector<std::vector<double>>> populate_Q_sellers(
        const std::vector<std::vector<double>>& Q_sellers_0,
        int n,
        int m,
        int S_sellers_cardinality
);

// Function to populate Q-matrix for the platform
std::vector<std::vector<double>> populate_Q_plat(
        const std::vector<double>& Q_plat_0,
        int plat_actions,
        int S_plat_cardinality
);

// Function to generate all possible values of tau
std::vector<double> populate_tau(
        double tau_min,
        double tau_step_size,
        int tau_values
);





#endif //FUNCTIONS_H