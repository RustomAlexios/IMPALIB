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
    int nDept_;                                     ///< number of departments
    int nTeams_;                                    ///< number of teams
    bool doFilter_;                                 ///< filtering flag
    impalib_type alpha_;                            ///< filtering parameter
    vector<vector<impalib_type>> extrinsicOutOld_;  ///< messages from departments to teams equality constraint before filtering

   public:
    void forward(int, vector<vector<impalib_type>> &, int, vector<vector<int>> &, const int *, const vector<vector<int>> &,
                 const vector<vector<impalib_type>> &) const;  ///< forward pass of forward-backward algorithm

    void backward(int, vector<vector<impalib_type>> &, int, vector<vector<int>> &, const int *, const vector<vector<int>> &,
                  const vector<vector<impalib_type>> &) const;  ///< backward pass of forward-backward algorithm

    void extrinsic_output_department_lhs(const vector<vector<int>> &, const vector<vector<impalib_type>> &, const vector<vector<impalib_type>> &, int, const vector<vector<impalib_type>> &, int,
                                         vector<vector<impalib_type>> &);  ///< extrinsic output of department constraint

    void team_to_knapsack_update(vector<vector<int>> &, vector<vector<impalib_type>> &, const vector<impalib_type> &, const vector<vector<impalib_type>> &,
                                 const vector<impalib_type> &) const;                                                    ///< calculate messages from teams to knapsack constraints
    void process_extrinsic_output_department(int, int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &);  ///< perform filtering (if needed) on messages from departments to teams

    Knapsack(int N_DEPARTMENTS, int N_TEAMS, bool FILT_FLAG, impalib_type ALPHA);  ///< constructor
};

/**
 * Construct Knapsack object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] FILT_FLAG: filtering flag
 * @param[in] ALPHA: filtering parameter (between 0 and 1)
 *
 */

inline Knapsack::Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA) : nDept_(N_DEPARTMENTS), nTeams_(N_TEAMS), doFilter_(FILT_FLAG), alpha_(ALPHA) {
    // Initialize extrinsic output department
    extrinsicOutOld_.reserve(nDept_);
    for (int department = 0; department < nDept_; department++) {
        extrinsicOutOld_.push_back(vector<impalib_type>(nTeams_, zero_value));
    }
};

/**
 * Perform forward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] department: index of department (knapsack constraint)
 * @param[in] capacity: capacity of department under
 * @param[in] nonzeroWeights: indices where connections are present between teams and departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of non-zero weight indices
 * @param[in] weights: weights of teams for the department under investigation
 * @param[in] team2KnapsackM: messages from teams to knapsack constraints
 * @param[out] forward: calculated forward messages on the trellis
 *
 */

inline void Knapsack::forward(const int department, vector<vector<impalib_type>> &forward, const int capacity, vector<vector<int>> &nonzeroWeights, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                              const vector<vector<int>> &weights, const vector<vector<impalib_type>> &team2KnapsackM) const {
    vector<int>::iterator upper;

    // Initialize initial forward messages
    vector<impalib_type> forward0(capacity + 1, zero_value);
    fill(forward0.begin() + 1, forward0.end(), value_inf);

    // Assign initial forward messages to the first stage
    forward[0] = forward0;

    // If the first non-zero weight index is not zero, assign initial messages to subsequent stages
    if (nonzeroWeights[department][0] != 0) {
        for (int j = 0; j < nonzeroWeights[department][0]; j++) {
            forward[j + 1] = forward0;
        }
    }

    for (int l = 0; l < pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department]; l++) {
        int t = nonzeroWeights[department][l];
        for (int a = 0; a <= capacity; a++) {
            if (a - weights[department][t] >= 0 && (!(t == 0 && a != weights[department][t]))) {
                // Update forward messages
                forward[t + 1][a] = min(forward0[a], forward0[a - weights[department][t]] + team2KnapsackM[department][t]);
            } else {
                forward[t + 1][a] = forward0[a];
            }
        }

        // Update initial forward messages
        forward0 = forward[t + 1];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (nonzeroWeights[department].size() != nTeams_) {
            if (!binary_search(nonzeroWeights[department].begin(), nonzeroWeights[department].end(), t + 1)) {
                if (t + 1 >= nonzeroWeights[department].back() && t + 1 < nTeams_) {
                    for (int j = t + 1; j < nTeams_; j++) {
                        forward[j + 1] = forward0;
                    }
                } else if (t + 1 < nonzeroWeights[department].back()) {
                    upper = upper_bound(nonzeroWeights[department].begin(), nonzeroWeights[department].end(), t);
                    size_t next_index = upper - nonzeroWeights[department].begin();
                    for (int j = t + 1; j < nonzeroWeights[department][next_index]; j++) {
                        forward[j + 1] = forward0;
                    }
                }
            }
        }
    }
}

/**
 * Perform backward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] department: index of department (knapsack constraint)
 * @param[in] capacity: capacity of department under
 * @param[in] nonzeroWeights: indices where connections are present between teams and departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of non-zero weight indices
 * @param[in] weights: weights of teams for the department under investigation
 * @param[in] team2KnapsackM: messages from teams to knapsack constraints
 * @param[out] backward: calculated backward messages on the trellis
 *
 */

inline void Knapsack::backward(const int department, vector<vector<impalib_type>> &backward, const int capacity, vector<vector<int>> &nonzeroWeights, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                               const vector<vector<int>> &weights, const vector<vector<impalib_type>> &team2KnapsackM) const {
    vector<int>::iterator upper;

    // Initialize initial backward messages
    vector<impalib_type> backward0(capacity + 1, zero_value);

    // Assign initial backward messages to the last stage
    backward[nTeams_] = backward0;

    // If the last non-zero weight index is not numTeams_ - 1, assign initial messages to previous stages
    if (nonzeroWeights[department][pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department] - 1] != nTeams_ - 1) {
        for (int j = nTeams_ - 1; j > nonzeroWeights[department][0]; j--) {
            backward[j] = backward0;
        }
    }

    for (int l = pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department] - 1; l >= 0; l--) {
        int t = nonzeroWeights[department][l];
        for (int a = 0; a <= capacity; a++) {
            if ((t > 0 && a + weights[department][t] <= capacity) || (a == 0 && t == 0)) {
                // Update backward messages
                backward[t][a] = min(backward0[a], backward0[a + weights[department][t]] + team2KnapsackM[department][t]);
            } else {
                backward[t][a] = backward0[a];
            }
        }

        // Update initial backward messages
        backward0 = backward[t];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department] - 1 != nTeams_) {
            if (!binary_search(nonzeroWeights[department].begin(), nonzeroWeights[department].end(), t - 1)) {
                if ((t - 1 <= nonzeroWeights[department][0]) && t - 1 >= 0) {
                    for (int j = t - 1; j >= 0; j--) {
                        backward[j] = backward0;
                    }
                } else if (t - 1 > nonzeroWeights[department][0]) {
                    upper = upper_bound(nonzeroWeights[department].begin(), nonzeroWeights[department].end(), t - 1);
                    size_t next_index = upper - nonzeroWeights[department].begin();
                    for (int j = t - 1; j > nonzeroWeights[department][next_index - 1]; j--) {
                        backward[j] = backward0;
                    }
                }
            }
        }
    }
}

/**
 * Calculate messages from departments to teams after obtaining forward and backward messages on a knapsack trellis following Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] weights: weights of teams for the department under investigation
 * @param[in] forward: calculated forward messages on the trellis
 * @param[in] team2KnapsackM: messages from teams to knapsack constraints
 * @param[in] department: index of department (knapsack constraint)
 * @param[in] backward: calculated backward messages on the trellis
 * @param[in] capacity: capacity of department under
 * @param[out] extrinsicOut: extrinsic output messages from department to teams
 *
 */

inline void Knapsack::extrinsic_output_department_lhs(const vector<vector<int>> &weights, const vector<vector<impalib_type>> &forward, const vector<vector<impalib_type>> &team2KnapsackM,
                                                      const int department, const vector<vector<impalib_type>> &backward, const int capacity, vector<vector<impalib_type>> &extrinsicOut) {
    // Initialize metric paths
    vector<impalib_type> metric_path_solid, metric_path_dash;

    for (int team = 0; team < weights[department].size(); team++) {
        metric_path_dash.clear();
        metric_path_solid.clear();

        // Check if the team has non-zero weight with the department
        if (weights[department][team] != 0) {
            // Calculate metric paths for solid (activated) and dashed (deactivated) states
            if (team == 0) {
                metric_path_solid.push_back(forward[team][team] + backward[team + 1][weights[department][team]] + team2KnapsackM[department][team]);
                metric_path_dash.push_back(forward[team][team] + backward[team + 1][0]);
            } else {
                for (int i = 0; i <= capacity; i++) {
                    if (i + weights[department][team] <= capacity) {
                        metric_path_solid.push_back(forward[team][i] + backward[team + 1][i + weights[department][team]] + team2KnapsackM[department][team]);
                    }
                    metric_path_dash.push_back(forward[team][i] + backward[team + 1][i]);
                }
            }
            // Calculate extrinsic output for the department
            extrinsicOut[department][team] =
                *min_element(metric_path_solid.begin(), metric_path_solid.end()) - *min_element(metric_path_dash.begin(), metric_path_dash.end()) - team2KnapsackM[department][team];
        } else {
            extrinsicOut[department][team] = zero_value;
        }
    }
}

/**
 * Calculate messages from teams to departments for the Knapsack-MWM problem
 *
 * @param[in] nonzeroWeights: indices where connections are present between teams and departments
 * @param[in] team2KnapsackM: messages from teams to knapsack constraints
 * @param[in] rewards: rewards of teams
 * @param[in] extrisicOut: extrinsic output messages from department to teams
 * @param[out] oric2TeamM: messages from ORIC to teams
 *
 */

inline void Knapsack::team_to_knapsack_update(vector<vector<int>> &nonzeroWeights, vector<vector<impalib_type>> &team2KnapsackM, const vector<impalib_type> &rewards,
                                              const vector<vector<impalib_type>> &extrisicOut, const vector<impalib_type> &oric2TeamM) const {
    for (int department = 0; department < extrisicOut.size(); department++) {
        // Create a list of remaining departments
        vector<int> remaining_departments(nDept_);
        iota(remaining_departments.begin(), remaining_departments.end(), 0);

        // Get unique edge indices for the current department
        vector<int> unique_edge_department(nonzeroWeights[department]);
        remaining_departments.erase(remaining_departments.begin() + department);

        for (auto t : remaining_departments) {
            // Calculate intersection of non-zero weight indices between the current department and other departments
            vector<int> intersection(nTeams_);
            vector<int>::iterator it;
            it = set_intersection(nonzeroWeights[department].begin(), nonzeroWeights[department].end(), nonzeroWeights[t].begin(), nonzeroWeights[t].end(), intersection.begin());
            intersection.resize(it - intersection.begin());

            // Update team to knapsack constraint messages
            for (auto l : intersection) {
                team2KnapsackM[department][l] = rewards[l] + extrisicOut[t][l] + oric2TeamM[l];

                // Remove the edge index from the unique list
                vector<int>::iterator position = std::find(unique_edge_department.begin(), unique_edge_department.end(), l);
                unique_edge_department.erase(position);
            }
        }

        // Update team to knapsack constraint messages for remaining unique edge department indices
        for (auto u : unique_edge_department) {
            team2KnapsackM[department][u] = rewards[u] + oric2TeamM[u];
        }
    }
}

/**
 * Calculate messages from departments to teams based on filtering conditions for the Knapsack-MWM problem
 *
 * @param[in] department: index of department (knapsack constraint)
 * @param[in] iter: iteration index
 * @param[in] extrinsicOutPre: extrinsic output messages from department to teams before filtering
 * @param[out] extrinsicOut: extrinsic output messages from department to teams after filtering
 *
 */

inline void Knapsack::process_extrinsic_output_department(const int department, const int iter, vector<vector<impalib_type>> &extrinsicOutPre, vector<vector<impalib_type>> &extrinsicOut) {
    if ((doFilter_) and (alpha_ != zero_value)) {
        vector<impalib_type> intermediate_dummy(extrinsicOutPre[department]), intermediate_old(extrinsicOutOld_[department]), intermediate_extrinsic;

        impalib_type w_1 = alpha_, w_2 = 1 - alpha_;

        // Calculate weighted extrinsic outputs
        transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
        transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

        if (iter == 0) {
            copy(intermediate_dummy.begin(), intermediate_dummy.end(), extrinsicOut[department].begin());
        } else {
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), std::back_inserter(intermediate_extrinsic), std::plus<impalib_type>());
            copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), extrinsicOut[department].begin());
        }

        copy(extrinsicOut[department].begin(), extrinsicOut[department].end(), extrinsicOutOld_[department].begin());
    }

    else {
        copy(extrinsicOutPre[department].begin(), extrinsicOutPre[department].end(), extrinsicOut[department].begin());
    }
}