// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for the Knapsack constraint for Knapsack-MWM problem
 */
class Knapsack {
   private:
    int numDepartments_;                         ///< number of departments
    int numTeams_;                               ///< number of teams
    bool filteringFlag_;                         ///< filtering flag
    impalib_type alpha_;                         ///< filtering parameter
    vector<vector<impalib_type>> extrinsic_old;  ///< messages from departments to teams equality constraint before filtering

   public:
    vector<vector<impalib_type>> forward(int, int, const vector<int> &, const int *, const vector<int> &,
                                         const vector<vector<impalib_type>> &);  ///< forward pass of forward-backward algorithm

    vector<vector<impalib_type>> backward(int, int, const vector<int> &, const int *, const vector<int> &,
                                          const vector<vector<impalib_type>> &);  ///< backward pass of forward-backward algorithm

    vector<impalib_type> extrinsic_output_department_lhs(const vector<int> &, const vector<vector<impalib_type>> &, const vector<vector<impalib_type>> &, int, const vector<vector<impalib_type>> &,
                                                         int);  ///< extrinsic_ output of department constraint

    void team_to_knapsack_update(const vector<vector<int>> &, vector<vector<impalib_type>> &, const vector<impalib_type> &, const vector<vector<impalib_type>> &,
                                 const vector<impalib_type> &);                                        ///< calculate messages from teams to knapsack constraints
    vector<impalib_type> process_extrinsic_output_department(int, int, const vector<impalib_type> &);  ///< perform filtering (if needed) on messages from departments to teams

    Knapsack(int N_DEPARTMENTS, int N_TEAMS, bool FILT_FLAG, impalib_type ALPHA);  ///< constructor
};

/**
 * Construct Knapsack object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] FILT_FLAG: filtering flag
 * @param[in] ALPHA: filtering parameter (between 0 and 1)
 * @param[out] numDepartments_: N_DEPARTMENTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] filteringFlag_: FILT_FLAG
 * @param[out] alpha_: ALPHA
 *
 */

Knapsack::Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), filteringFlag_(FILT_FLAG), alpha_(ALPHA), extrinsic_old(numDepartments_, vector<impalib_type>(numTeams_, 0)){};

/**
 * Perform forward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] idx: index of department (knapsack constraint)
 * @param[in] capacity: capacity of department under
 * @param[in] idx_nonzero: indices where connections are present between teams and departments
 * @param[in] idx_nonzero_sizes: sizes of non-zero weight indices
 * @param[in] rTeamsWeightsPerDepartment: weights of teams for the department under investigation
 * @param[in] team2knapsack: messages from teams to knapsack constraints
 * @param[out] rStageForwardMessages: calculated forward messages on the trellis
 *
 */

vector<vector<impalib_type>> Knapsack::forward(const int idx, const int capacity, const vector<int> &idx_nonzero, const int *idx_nonzero_sizes, const vector<int> &weights,
                                               const vector<vector<impalib_type>> &team2knapsack) {
    vector<int>::const_iterator upper;
    vector<vector<impalib_type>> forward(numTeams_ + 1, vector<impalib_type>(capacity + 1, zero_value));

    // Initialize initial forward messages
    vector<impalib_type> forward0(capacity + 1, zero_value);
    fill(forward0.begin() + 1, forward0.end(), value_inf);

    // Assign initial forward messages to the first stage
    forward[0] = forward0;

    // If the first non-zero weight index is not zero, assign initial messages to subsequent stages
    if (idx_nonzero[0] != 0) {
        for (int j = 0; j < idx_nonzero[0]; j++) {
            forward[j + 1] = forward0;
        }
    }

    for (int l = 0; l < idx_nonzero_sizes[idx]; l++) {
        int t = idx_nonzero[l];
        for (int a = 0; a <= capacity; a++) {
            if (a - weights[t] >= 0 && (!(t == 0 && a != weights[t]))) {
                // Update forward messages
                forward[t + 1][a] = min(forward0[a], forward0[a - weights[t]] + team2knapsack[idx][t]);
            } else {
                forward[t + 1][a] = forward0[a];
            }
        }

        // Update initial forward messages
        forward0 = forward[t + 1];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (idx_nonzero.size() != numTeams_) {
            if (!binary_search(idx_nonzero.begin(), idx_nonzero.end(), t + 1)) {
                if (t + 1 >= idx_nonzero.back() && t + 1 < numTeams_) {
                    for (int j = t + 1; j < numTeams_; j++) {
                        forward[j + 1] = forward0;
                    }
                } else if (t + 1 < idx_nonzero.back()) {
                    upper = upper_bound(idx_nonzero.begin(), idx_nonzero.end(), t);
                    size_t next_index = upper - idx_nonzero.begin();
                    for (int j = t + 1; j < idx_nonzero[next_index]; j++) {
                        forward[j + 1] = forward0;
                    }
                }
            }
        }
    }
    return forward;
}

/**
 * Perform backward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] idx: index of department (knapsack constraint)
 * @param[in] capacity: capacity of department under
 * @param[in] idx_nonzero: indices where connections are present between teams and departments
 * @param[in] idx_nonzero_sizes: sizes of non-zero weight indices
 * @param[in] weights: weights of teams for the department under investigation
 * @param[in] team2knapsack: messages from teams to knapsack constraints
 * @param[out] rStageBackwardMessages: calculated backward messages on the trellis
 *
 */

vector<vector<impalib_type>> Knapsack::backward(const int idx, const int capacity, const vector<int> &idx_nonzero, const int *idx_nonzero_sizes, const vector<int> &weights,
                                                const vector<vector<impalib_type>> &team2knapsack) {
    vector<int>::const_iterator upper;
    vector<vector<impalib_type>> backward(numTeams_ + 1, vector<impalib_type>(capacity + 1, zero_value));

    // Initialize initial backward messages
    vector<impalib_type> backward0(capacity + 1, zero_value);

    // Assign initial backward messages to the last stage
    backward[numTeams_] = backward0;

    // If the last non-zero weight index is not numTeams_ - 1, assign initial messages to previous stages
    if (idx_nonzero[idx_nonzero_sizes[idx] - 1] != numTeams_ - 1) {
        for (int j = numTeams_ - 1; j > idx_nonzero[0]; j--) {
            backward[j] = backward0;
        }
    }

    for (int l = idx_nonzero_sizes[idx] - 1; l >= 0; l--) {
        int t = idx_nonzero[l];
        for (int a = 0; a <= capacity; a++) {
            if ((t > 0 && a + weights[t] <= capacity) || (a == 0 && t == 0)) {
                // Update backward messages
                backward[t][a] = min(backward0[a], backward0[a + weights[t]] + team2knapsack[idx][t]);
            } else {
                backward[t][a] = backward0[a];
            }
        }

        // Update initial backward messages
        backward0 = backward[t];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (idx_nonzero_sizes[idx] - 1 != numTeams_) {
            if (!binary_search(idx_nonzero.begin(), idx_nonzero.end(), t - 1)) {
                if ((t - 1 <= idx_nonzero[0]) && t - 1 >= 0) {
                    for (int j = t - 1; j >= 0; j--) {
                        backward[j] = backward0;
                    }
                } else if (t - 1 > idx_nonzero[0]) {
                    upper = upper_bound(idx_nonzero.begin(), idx_nonzero.end(), t - 1);
                    size_t next_index = upper - idx_nonzero.begin();
                    for (int j = t - 1; j > idx_nonzero[next_index - 1]; j--) {
                        backward[j] = backward0;
                    }
                }
            }
        }
    }
    return backward;
}

/**
 * Calculate messages from departments to teams after obtaining forward and backward messages on a knapsack trellis following Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] weights: weights of teams for the department under investigation
 * @param[in] forward: calculated forward messages on the trellis
 * @param[in] team2knapsack: messages from teams to knapsack constraints
 * @param[in] idx: index of department (knapsack constraint)
 * @param[in] backward: calculated backward messages on the trellis
 * @param[in] capacity: capacity of department under
 * @param[out] rExtrinsicOutputDepartment: extrinsic_ output messages from department to teams
 *
 */

vector<impalib_type> Knapsack::extrinsic_output_department_lhs(const vector<int> &weights, const vector<vector<impalib_type>> &forward, const vector<vector<impalib_type>> &team2knapsack,
                                                               const int idx, const vector<vector<impalib_type>> &backward, const int capacity) {
    vector<impalib_type> extrinsic(weights.size());
    vector<impalib_type> solid, dash;

    for (int i = 0; i < weights.size(); i++) {
        dash.clear();
        solid.clear();

        // Check if the team has non-zero weight with the department
        if (weights[i] != 0) {
            // Calculate metric paths for solid (activated) and dashed (deactivated) states
            if (i == 0) {
                solid.push_back(forward[i][i] + backward[i + 1][weights[i]] + team2knapsack[idx][i]);
                dash.push_back(forward[i][i] + backward[i + 1][0]);
            } else {
                for (int j = 0; j <= capacity; j++) {
                    if (j + weights[i] <= capacity) {
                        solid.push_back(forward[i][j] + backward[i + 1][j + weights[i]] + team2knapsack[idx][i]);
                    }
                    dash.push_back(forward[i][j] + backward[i + 1][j]);
                }
            }
            // Calculate extrinsic_ output for the department
            extrinsic[i] = *min_element(solid.begin(), solid.end()) - *min_element(dash.begin(), dash.end()) - team2knapsack[idx][i];
        } else {
            extrinsic[i] = zero_value;
        }
    }
    return extrinsic;
}

/**
 * Calculate messages from teams to departments for the Knapsack-MWM problem
 *
 * @param[in] idx_nonzero: indices where connections are present between teams and departments
 * @param[in] team2knapsack: messages from teams to knapsack constraints
 * @param[in] reward: rewards of teams
 * @param[in] extrinsic: extrinsic_ output messages from department to teams
 * @param[out] oric2team: messages from ORIC to teams
 *
 */

void Knapsack::team_to_knapsack_update(const vector<vector<int>> &idx_nonzero, vector<vector<impalib_type>> &team2knapsack, const vector<impalib_type> &reward,
                                       const vector<vector<impalib_type>> &extrinsic, const vector<impalib_type> &oric2team) {
    for (int i = 0; i < extrinsic.size(); i++) {
        // Create a list of remaining departments
        vector<int> remaining(numDepartments_);
        iota(remaining.begin(), remaining.end(), 0);

        // Get unique edge indices for the current department
        vector<int> unique_edges(idx_nonzero[i]);
        remaining.erase(remaining.begin() + i);

        for (auto t : remaining) {
            // Calculate intersection of non-zero weight indices between the current department and other departments
            vector<int> intersection(numTeams_);
            vector<int>::iterator it;
            it = set_intersection(idx_nonzero[i].begin(), idx_nonzero[i].end(), idx_nonzero[t].begin(), idx_nonzero[t].end(), intersection.begin());
            intersection.resize(it - intersection.begin());

            // Update team to knapsack constraint messages
            for (auto l : intersection) {
                team2knapsack[i][l] = reward[l] + extrinsic[t][l] + oric2team[l];

                // Remove the edge index from the unique list
                vector<int>::iterator position = std::find(unique_edges.begin(), unique_edges.end(), l);
                unique_edges.erase(position);
            }
        }

        // Update team to knapsack constraint messages for remaining unique edge department indices
        for (auto u : unique_edges) {
            team2knapsack[i][u] = reward[u] + oric2team[u];
        }
    }
}

/**
 * Calculate messages from departments to teams based on filtering conditions for the Knapsack-MWM problem
 *
 * @param[in] idx: index of department (knapsack constraint)
 * @param[in] iter: iteration index
 * @param[in] extrinsic_in: extrinsic_ output messages from department to teams before filtering
 * @param[out] extrinsic_out: extrinsic_ output messages from department to teams after filtering
 *
 */

vector<impalib_type> Knapsack::process_extrinsic_output_department(const int idx, const int iter, const vector<impalib_type> &extrinsic_in) {
    if (!filteringFlag_) {
        return extrinsic_in;
    }

    // Calculate weighted extrinsic outputs
    auto extrinsic_out = extrinsic_in;
    for (int j = 0; j < extrinsic_out.size(); ++j) {
        extrinsic_out[j] = (1 - alpha_) * extrinsic_out[j] + (alpha_)*extrinsic_old[idx][j];
    }

    extrinsic_old[idx] = extrinsic_out;
    return extrinsic_out;
}