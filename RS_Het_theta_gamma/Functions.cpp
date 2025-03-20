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

// Function to generate all possible combinations of sellers' actions
std::vector<std::vector<double>> populate_p_sellers(
        const std::vector<double>& A_sellers,
        int n,
        int m,
        int plat_actions)
{
    // Generate all possible combinations of actions for sellers
    std::vector<std::vector<double>> A_seller_comb;

    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            for (int j = 0; j < m; ++j) {
                A_seller_comb.push_back({A_sellers[j]});
            }
        } else {
            std::vector<std::vector<double>> newAs;
            for (int j = 0; j < m; ++j) {
                for (const auto& prev : A_seller_comb) {
                    std::vector<double> combination = prev;
                    combination.push_back(A_sellers[j]);
                    newAs.push_back(combination);
                }
            }
            A_seller_comb = newAs;
        }
    }

    // Replicate rows for plat_actions times
    std::vector<std::vector<double>> p_sellers;
    for (int i = 0; i < plat_actions; ++i) {
        p_sellers.insert(p_sellers.end(), A_seller_comb.begin(), A_seller_comb.end());
    }

    return p_sellers;
}


// Function to determine demand for each state for each seller
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
) {
    // Initialize demand vector
    std::vector<std::vector<double>> d_sellers(S_sellers_cardinality, std::vector<double>(n, 0.0));

    for (int s{ 0 }; s < S_sellers_cardinality; ++s) {
        // Get current platform action
        int plat_action = std::ceil(static_cast<double>(s + 1) / p_sellers_comb) - 1;

        if (plat_action == 0)
        {
            d_sellers[s][0] = tau * (
                    gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) /
                             (std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu))) +
                    (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) /
                                   (std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu)))
            ) +
                              (1 - tau) * (
                                      gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) /
                                               (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                                std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                      (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) /
                                                     (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                      std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
                              );

            d_sellers[s][1] = (1 - tau) * (
                    gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) /
                             (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                              std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                    (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) /
                                   (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                    std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
            );
        }
        else if (plat_action == 1)
        {
            d_sellers[s][0] = tau * (
                    gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) /
                             (std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu)))
            ) +
                              (1 - tau) * (
                                      gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) /
                                               (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                                std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                      (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) /
                                                     (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                      std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
                              );

            d_sellers[s][1] = tau * (
                    (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) /
                                   (std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu)))
            ) +
                              (1 - tau) * (
                                      gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) /
                                               (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                                std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                      (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) /
                                                     (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                      std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
                              );
        }
        else if (plat_action == 2)
        {
            d_sellers[s][0] = tau * (
                    (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) /
                                   (std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu)))
            ) +
                              (1 - tau) * (
                                      gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) /
                                               (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                                std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                      (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) /
                                                     (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                      std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
                              );

            d_sellers[s][1] = tau * (
                    gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) /
                             (std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu)))
            ) +
                              (1 - tau) * (
                                      gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) /
                                               (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                                std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                      (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) /
                                                     (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                      std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))
                              );
        }
        else if (plat_action == 3) {
            d_sellers[s][0] =
                    (1 - tau) * (gammma * (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) /
                                          (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                           std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                 (1 - gammma) * (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) /
                                                (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                 std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu))));

            d_sellers[s][1] =
                    tau * (gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) /
                                    (std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu))) +
                           (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) /
                                          (std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu))))
                    +
                    (1 - tau) * (gammma * (std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) /
                                          (std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) +
                                           std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu))) +
                                 (1 - gammma) * (std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) /
                                                (std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) +
                                                 std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu))));
        }
    }

    return d_sellers;
}

// Function to determine average profits for each firm when setting each of the *m* possible prices for each platform's action
std::vector<std::vector<double>> populate_pi_sellers(
        const std::vector<std::vector<double>>& d_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<int>& mc,
        double f,
        int S_sellers_cardinality,
        int n
)
{
    // Initialize vector to store profits for each state for each firm given the platform's action
    std::vector<std::vector<double>> pi_sellers(
            S_sellers_cardinality, std::vector<double>(n, 0.0)
    );
    for (int s{ 0 }; s < S_sellers_cardinality; ++s)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            pi_sellers[s][i] = ((1 - f) * p_sellers[s][i] - mc[i]) * d_sellers[s][i];
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
        double gammma,
        double tau,
        int a_0,
        int S_sellers_cardinality,
        int n,
        int p_sellers_comb
)
{
    // Initialize consumer surplus vector
    std::vector<double> cs(S_sellers_cardinality, 0.0);

    for (int s{ 0 }; s < S_sellers_cardinality; ++s)
    {
        // Get current platform action
        int plat_action = std::ceil(static_cast<double>(s + 1) / p_sellers_comb) - 1;

        if (plat_action == 0)
        {
            cs[s] = mu * (tau * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu)))
                                 + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu))))
                          + (1 - tau) * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) + std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu)))
                                         + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) + std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))));
        }
        else if (plat_action == 1)
        {
            cs[s] = mu * (tau * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu)))
                                 + (1 - gammma) * (std::log(std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu))))
                    + (1 - tau) * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) + std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu)))
                                   + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) + std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))));
        }
        else if (plat_action == 2)
        {
            cs[s] = mu * (tau * (gammma * (std::log(std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu)))
                                 + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0])) / mu) + std::exp(a_0 / mu))))
                          + (1 - tau) * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) + std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu)))
                                         + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) + std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))));
        }
        else if (plat_action == 3)
        {
            cs[s] = mu * (tau * (gammma * (std::log(std::exp((a[1][0] - (theta[0] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu)))
                                 + (1 - gammma) * (std::log(std::exp((a[1][1] - (theta[1] * p_sellers[s][1])) / mu) + std::exp(a_0 / mu))))
                          + (1 - tau) * (gammma * (std::log(std::exp((a[0][0] - (theta[0] * p_sellers[s][0]) - c[0]) / mu) + std::exp((a[1][0] - (theta[0] * p_sellers[s][1]) - c[0]) / mu) + std::exp(a_0 / mu)))
                                         + (1 - gammma) * (std::log(std::exp((a[0][1] - (theta[1] * p_sellers[s][0]) - c[1]) / mu) + std::exp((a[1][1] - (theta[1] * p_sellers[s][1]) - c[1]) / mu) + std::exp(a_0 / mu)))));
        }
    }

    return cs;
}

// Function to determine average profits for the platform given their action
std::vector<double> populate_pi_plat(
        const std::vector<std::vector<double>>& d_sellers,
        const std::vector<std::vector<double>>& p_sellers,
        const std::vector<double>& cs,
        double omega,
        double f,
        int S_plat_cardinality,
        int n
)
{
    // Initialize platform profits vector
    std::vector<double> pi_plat(S_plat_cardinality, 0.0);

    for (int s { 0 }; s < S_plat_cardinality; ++s)
    {
        // Calculate revenues
        double rvn = 0.0;
        for (int i { 0 }; i < n; ++i)
            rvn += p_sellers[s][i] * d_sellers[s][i];

        // Compute platform profit for this state and action
        pi_plat[s] = omega * (f * rvn) + (1 - omega) * cs[s];
    }

    return pi_plat;
}

// Function to determine average profits for each firm when setting each of the *m* possible prices for each platform's action
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
    for (int i { 0 }; i < n; ++i)
    {
        // Loop over prices
        for (int j { 0 }; j < m; ++j)
        {
            // Initialize vector to store indices of states where seller *i* sets the *j*'th price
            std::vector<int> indices;

            // Loop over the possible states
            for (int s{ 0 }; s < S_sellers_cardinality; ++s)
            {
                // If in state *s* seller *i* sets price *j*
                if (p_sellers[s][i] == A_sellers[j])
                {
                    // Add this state index to indices vector
                    indices.push_back(s);
                }
            }

            // Calculate profit summation for seller *i* when setting price *j*
            double sum{ 0.0 };
            for (int index : indices) {
                sum += pi_sellers.at(index)[i];
            }
            // Find number of times seller *i* sets price *j*
            int count = indices.size();

            // Calculate mean profit for seller *i* when setting price *j* (avoiding division by 0)
            Q_sellers_0[j][i] = (count > 0) ? (sum / count) / (1.0 - delta) : 0.0;
        }
    }

    return Q_sellers_0;
}

// Function to determine average profits for the platform given their action
std::vector<double> populate_Q_plat_0(
        const std::vector<double>& pi_plat,
        double delta,
        int plat_actions,
        int p_sellers_comb
)
{
    // Initialize platform's initial Q-matrix
    std::vector<double> Q_plat_0(plat_actions, 0.0);

    for (int plat_action{ 0 }; plat_action < plat_actions; ++plat_action) {
        int start_idx = plat_action * p_sellers_comb;
        int end_idx = (plat_action + 1) * p_sellers_comb;

        double sum = 0.0;
        int count = 0;

        for (int i {start_idx}; i < end_idx; ++i)
        {
            sum += pi_plat[i];
            ++count;
        }

        Q_plat_0[plat_action] = (sum / count) / (1 - delta);
    }

    return Q_plat_0;
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
    std::vector<std::vector<std::vector<double>>> Q_sellers(S_sellers_cardinality, std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0)));

    // Initialize each firm's Q-matrix with average profits when setting each of the m possible prices given the platform's action
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

// Function to populate Q-matrix for the platform
std::vector<std::vector<double>> populate_Q_plat(
        const std::vector<double>& Q_plat_0,
        int plat_actions,
        int S_plat_cardinality
)
{
    // Initialize Q-matrix for storage
    std::vector<std::vector<double>> Q_plat(S_plat_cardinality, std::vector<double>(plat_actions, 0.0));

    // Initialize the platform's Q-matrix using Q_plat_0
    for (int s{ 0 }; s < S_plat_cardinality; ++s)
    {
        for (int j{ 0 }; j < plat_actions; ++j)
            Q_plat[s][j] = Q_plat_0[j];
    }

    return Q_plat;
}

// Function to generate all possible values of gammma
std::vector<double> populate_gammma(
        double gammma_min,
        double gammma_step_size,
        int gammma_values
)
{
    // Create vector to store values of gammma
    std::vector<double> gammma_vector(gammma_values, 0.0);

    // Populate vector with *gammma_values* equally spaced points
    for (int gammma{ 0 }; gammma < gammma_values; ++gammma)
        gammma_vector[gammma] = gammma_min + (gammma * gammma_step_size);

    return gammma_vector;
}

// Function to generate all possible values of theta_1
std::vector<double> populate_theta_1(
        double theta_1_min,
        double theta_1_step_size,
        int theta_1_values
)
{
    // Create vector to store values of theta_1
    std::vector<double> theta_1_vector(theta_1_values, 0.0);

    // Populate vector with *theta_1_values* equally spaced points
    for (int theta_1{ 0 }; theta_1 < theta_1_values; ++theta_1)
        theta_1_vector[theta_1] = theta_1_min + (theta_1 * theta_1_step_size);

    return theta_1_vector;
}








