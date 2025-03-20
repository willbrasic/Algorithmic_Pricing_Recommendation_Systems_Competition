#include "Functions.h"


//////////////////////////////
// Function Definitions
//////////////////////////////


// Function to populate sellers' action spaces
std::vector<double> populate_A_sellers(
        int m,
        double A_sellers_min,
        double A_sellers_step_size
)
{
    // Initialize vector to store action space
    std::vector<double> action_space(m, 0.0);

    // Populate action space
    for (int j{ 0 }; j < m; ++j)
        action_space[j] = A_sellers_min + (j * A_sellers_step_size);

    return action_space;
}

// Function to determine prices for each state for each firm
std::vector<std::vector<double>> populate_p_sellers(
        const std::vector<double>& A_sellers,
        int n,
        int m
)
{
    // All possible combinations of actions for firms
    std::vector<std::vector<double>> A_seller_comb;
    for (int i{ 0 }; i < n; ++i)
    {
        if (i == 0)
        {
            for (int j = 0; j < m; ++j)
            {
                A_seller_comb.push_back({ A_sellers[j] });
            }
        }
        else
        {
            std::vector<std::vector<double>> newAs;
            for (int j{ 0 }; j < m; ++j)
            {
                for (const auto& prev : A_seller_comb)
                {
                    std::vector<double> combination = prev;
                    combination.push_back(A_sellers[j]);
                    newAs.push_back(combination);
                }
            }
            A_seller_comb = newAs;
        }
    }

    return A_seller_comb;
}

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
)
{
    // Initialize vector to store demand for each state for each firm
    std::vector<std::vector<double>> d_sellers(
            S_sellers_cardinality,
            std::vector<double>(n, 0.0)
    );

    for (int s{0}; s < S_sellers_cardinality; ++s)
    {
        d_sellers[s][0] =
                gamma * (exp((a[0][0] - theta[0] * p_sellers[s][0] - c[0]) / mu) /
                         (exp((a[0][0] - theta[0] * p_sellers[s][0] - c[0]) / mu) +
                          exp((a[1][0] - theta[0] * p_sellers[s][1] - c[0]) / mu) + exp(a_0 / mu))) +
                (1 - gamma) * (exp((a[0][1] - theta[1] * p_sellers[s][0] - c[1]) / mu) /
                               (exp((a[0][1] - theta[1] * p_sellers[s][0] - c[1]) / mu) +
                                exp((a[1][1] - theta[1] * p_sellers[s][1] - c[1]) / mu) + exp(a_0 / mu)));
        d_sellers[s][1] =
                gamma * (exp((a[1][0] - theta[0] * p_sellers[s][1] - c[0]) / mu) /
                         (exp((a[0][0] - theta[0] * p_sellers[s][0] - c[0]) / mu) +
                          exp((a[1][0] - theta[0] * p_sellers[s][1] - c[0]) / mu) + exp(a_0 / mu))) +
                (1 - gamma) * (exp((a[1][1] - theta[1] * p_sellers[s][1] - c[1]) / mu) /
                               (exp((a[0][1] - theta[1] * p_sellers[s][0] - c[1]) / mu) +
                                exp((a[1][1] - theta[1] * p_sellers[s][1] - c[1]) / mu) + exp(a_0 / mu)));
    }


    return d_sellers;
}

// Function to determine average profits for each firm when setting each of the *m* possible prices
std::vector<std::vector<double>> populate_pi_sellers(
        const std::vector<std::vector<double>>& d_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<int>& mc,
        int S_sellers_cardinality,
        int n,
        double f
)
{
    // Initialize vector to store profits for each state for each firm
    std::vector<std::vector<double>> pi_sellers(
            S_sellers_cardinality,
            std::vector<double>(n, 0.0)
    );
    for (int j{ 0 }; j < S_sellers_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            pi_sellers[j][i] = ((1 - f) * p_sellers[j][i] - mc[i]) * d_sellers[j][i];
        }
    }

    return pi_sellers;
}

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
)
{
    // Initialize consumer surplus vector
    std::vector<double> cs(S_sellers_cardinality, 0.0);

    for (int s{ 0 }; s < S_sellers_cardinality; ++s) {
        // Compute terms for gamma-weighted and (1-gamma)-weighted log sums
        double exp_sum_gamma = 0.0;
        double exp_sum_one_minus_gamma = 0.0;

        for (int i{ 0 }; i < n; ++i)
        {
            for (int j{ 0 }; j < k; ++j)
            {
                // First consumer type
                if (j == 0)
                    exp_sum_gamma += exp((a[i][j] - theta[j] * p_sellers[s][i] - c[0]) / mu);
                // Second consumer type
                if (j == 1)
                    exp_sum_one_minus_gamma += exp((a[i][j] - theta[j] * p_sellers[s][i] - c[1]) / mu);
            }
        }

        // Add outside option
        exp_sum_gamma += exp(a_0 / mu);
        exp_sum_one_minus_gamma += exp(a_0 / mu);

        // Calculate consumer surplus for state s
        cs[s] = mu * (gamma * log(exp_sum_gamma) + (1 - gamma) * log(exp_sum_one_minus_gamma));
    }

    return cs;
}

// Function to determine average profits for each firm when setting each of the *m* possible prices
std::vector<std::vector<double>> populate_Q_sellers_0(
        const std::vector<std::vector<double>>& pi_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<double>& A_sellers,
        double delta,
        int S_sellers_cardinality,
        int n,
        int m
)
{
    // Initialize matrix for storing average profits for each seller for each possible price
    std::vector<std::vector<double>> Q_sellers_0(m, std::vector<double>(n, 0.0));

    // Loop over sellers
    for (int i{0}; i < n; ++i)
    {
        // Loop over possible prices for each seller
        for (int j{0}; j < m; ++j)
        {
            // Initialize vector to store profits for all states where seller *i* sets price A_sellers[j]
            std::vector<double> relevant_profits;

            // Loop over states
            for (int s{0}; s < S_sellers_cardinality; ++s)
            {   // Match price
                if (p_sellers[s][i] == A_sellers[j])
                    // Add profit for this state directly to *relevant_profits*
                    relevant_profits.push_back(pi_sellers[s][i]);
            }

            // Calculate the mean of relevant profits and normalize by (1 - *delta*)
            if (!relevant_profits.empty())
            {
                double mean_profit = std::accumulate(relevant_profits.begin(), relevant_profits.end(), 0.0) /
                                     relevant_profits.size();
                Q_sellers_0[j][i] = mean_profit / (1.0 - delta);
            }
            else
            {
                Q_sellers_0[j][i] = 0.0; // If no relevant profits, set to 0
            }
        }
    }

    return Q_sellers_0;
}

// Function to populate Q-matrix for each firm
std::vector<std::vector<std::vector<double>>> populate_Q_sellers(
        const std::vector<std::vector<double>>& Q_sellers_0,
        int n,
        int m,
        int S_sellers_cardinality
)
{
    // Initialize Q-matrix for storage
    std::vector<std::vector<std::vector<double>>> Q_sellers(
            S_sellers_cardinality,
            std::vector<std::vector<double>>(m,
                                             std::vector<double>(n, 0.0))
    );

    // Initialize each firm's Q-matrix with average profits when setting each of the m possible prices
    for (int i{ 0 }; i < n; ++i)
    {
        for (int s{ 0 }; s < S_sellers_cardinality; ++s)
        {
            for (int j{ 0 }; j < m; ++j)
                // Initialize Q-matrix in episode *e* for firm *i* when setting price *j* in state *s* as the
                // average profit when agent *i* sets price *j*
                Q_sellers[s][j][i] = Q_sellers_0[j][i];
        }
    }

    return Q_sellers;
}

// Function to generate all possible values of gamma
std::vector<double> populate_gamma(
        double gamma_min,
        double gamma_step_size,
        int gamma_values
)
{
    // Create vector to store values of gamma
    std::vector<double> gamma_vector(gamma_values, 0.0);

    // Populate vector with *gamma_values* equally spaced points
    for (int gamma{ 0 }; gamma < gamma_values; ++gamma)
        gamma_vector[gamma] = gamma_min + (gamma * gamma_step_size);

    return gamma_vector;
}

// Function to generate all possible values of theta
std::vector<std::vector<double>> populate_theta(
        double theta_1_min,
        double theta_1_step_size,
        int theta_1_values,
        int k
)
{
    // Create vector to store values of theta
    std::vector<std::vector<double>> theta_vector(theta_1_values, std::vector<double>(k, 0.0));

    // Populate vector with *theta_1_values* equally spaced points
    for (int theta_1{ 0 }; theta_1 < theta_1_values; ++theta_1)
    {
        theta_vector[theta_1][0] = theta_1_min + (theta_1 * theta_1_step_size);
        theta_vector[theta_1][1] = 1.0;
    }

    return theta_vector;
}






