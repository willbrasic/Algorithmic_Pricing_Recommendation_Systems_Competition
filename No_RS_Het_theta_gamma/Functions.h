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
        int m
);

// Function to determine demand for each state for each firm
std::vector<std::vector<double>> populate_d_sellers(
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<std::vector<double>>& a,
        const std::vector<double>& theta,
        const std::vector<double>& c,
        double gamma,
        double mu,
        int a_0,
        int n,
        int S_sellers_cardinality
);

// Function to determine profit for each state for each firm
std::vector<std::vector<double>> populate_pi_sellers(
        const std::vector<std::vector<double>>& d_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<int>& mc,
        int S_sellers_cardinality,
        int n,
        double f
);

// Function to determine consumer surplus for each state
std::vector<double> populate_cs(
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<std::vector<double>>& a,
        const std::vector<double>& theta,
        const std::vector<double>& c,
        double mu,
        double gamma,
        int a_0,
        int S_sellers_cardinality,
        int n,
        int k
);

// Function to determine average profits for each firm when setting each of the *m* possible prices
std::vector<std::vector<double>> populate_Q_sellers_0(
        const std::vector<std::vector<double>>& pi_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<double>& A_sellers,
        double delta,
        int S_sellers_cardinality,
        int n,
        int m
);

// Function to populate Q-matrix for each firm
std::vector<std::vector<std::vector<double>>> populate_Q_sellers(
        const std::vector<std::vector<double>>& Q_sellers_0,
        int n,
        int m,
        int S_sellers_cardinality
);

// Function to generate all possible values of gamma
std::vector<double> populate_gamma(
        double gamma_min,
        double gamma_step_size,
        int gamma_values
);

// Function to generate all possible values of theta
std::vector<std::vector<double>> populate_theta(
        double theta_1_min,
        double theta_1_step_size,
        int theta_1_values,
        int k
);




#endif //FUNCTIONS_H