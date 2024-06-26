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
class Knapsack
{
private:
    int                          numDepartments_; ///< number of departments
    int                          numTeams_; ///< number of teams
    bool                         filteringFlag_; ///< filtering flag
    impalib_type                 alpha_; ///< filtering parameter
    vector<vector<impalib_type>> extrinsicOutputDepartmentOld_; ///< messages from departments to teams equality constraint before filtering

public:
    void forward(int, vector<vector<impalib_type>> &, int, vector<vector<int>> &, const int *, const vector<vector<int>> &,
                 const vector<vector<impalib_type>> &) const; ///< forward pass of forward-backward algorithm

    void backward(int, vector<vector<impalib_type>> &, int, vector<vector<int>> &, const int *, const vector<vector<int>> &,
                  const vector<vector<impalib_type>> &) const; ///< backward pass of forward-backward algorithm

    void extrinsic_output_department_lhs(const vector<vector<int>> &, const vector<vector<impalib_type>> &,
                                         const vector<vector<impalib_type>> &, int, const vector<vector<impalib_type>> &, int,
                                         vector<vector<impalib_type>> &); ///< extrinsic output of department constraint

    void team_to_knapsack_update(vector<vector<int>> &, vector<vector<impalib_type>> &, const vector<impalib_type> &,
                                 const vector<vector<impalib_type>> &, const vector<impalib_type> &) const; ///< calculate messages from teams to knapsack constraints
    void process_extrinsic_output_department(int, int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &); ///< perform filtering (if needed) on messages from departments to teams

    Knapsack(int N_DEPARTMENTS, int N_TEAMS, bool FILT_FLAG, impalib_type ALPHA); ///< constructor
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

inline Knapsack::Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), filteringFlag_(FILT_FLAG), alpha_(ALPHA)
{

    // Initialize extrinsic output department
    extrinsicOutputDepartmentOld_.reserve(numDepartments_);
    for (int department_index = 0; department_index < numDepartments_; department_index++)
    {
        extrinsicOutputDepartmentOld_.push_back(vector<impalib_type>(numTeams_, zero_value));
    }
};

/**
 * Perform forward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] department_index: index of department (knapsack constraint)
 * @param[in] max_state_department: capacity of department under 
 * @param[in] rNonZeroWeightIndices: indices where connections are present between teams and departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of non-zero weight indices
 * @param[in] rTeamsWeightsPerDepartment: weights of teams for the department under investigation
 * @param[in] rTeam2KnapsackM: messages from teams to knapsack constraints
 * @param[out] rStageForwardMessages: calculated forward messages on the trellis
 * 
 */

inline void Knapsack::forward(const int department_index, vector<vector<impalib_type>> &rStageForwardMessages,
                       const int max_state_department, vector<vector<int>> &rNonZeroWeightIndices,
                       const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const vector<vector<int>> &rTeamsWeightsPerDepartment,
                       const vector<vector<impalib_type>> &rTeam2KnapsackM) const
{
    vector<int>::iterator upper;

    // Initialize initial forward messages
    vector<impalib_type> initial_forward_messages(max_state_department + 1, zero_value);
    fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

    // Assign initial forward messages to the first stage
    rStageForwardMessages[0] = initial_forward_messages;

    // If the first non-zero weight index is not zero, assign initial messages to subsequent stages
    if (rNonZeroWeightIndices[department_index][0] != 0)
    {
        for (int j = 0; j < rNonZeroWeightIndices[department_index][0]; j++)
        {
            rStageForwardMessages[j + 1] = initial_forward_messages;
        }
    }

    for (int l = 0; l < pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]; l++)
    {
        int t = rNonZeroWeightIndices[department_index][l];
        for (int a = 0; a <= max_state_department; a++)
        {
            if (a - rTeamsWeightsPerDepartment[department_index][t] >= 0
                && (!(t == 0 && a != rTeamsWeightsPerDepartment[department_index][t])))
            {
                // Update forward messages
                rStageForwardMessages[t + 1][a] =
                    min(initial_forward_messages[a],
                        initial_forward_messages[a - rTeamsWeightsPerDepartment[department_index][t]]
                            + rTeam2KnapsackM[department_index][t]);
            }
            else
            {
                rStageForwardMessages[t + 1][a] = initial_forward_messages[a];
            }
        }

        // Update initial forward messages
        initial_forward_messages = rStageForwardMessages[t + 1];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (rNonZeroWeightIndices[department_index].size() != numTeams_)
        {
            if (!binary_search(rNonZeroWeightIndices[department_index].begin(),
                               rNonZeroWeightIndices[department_index].end(), t + 1))
            {
                if (t + 1 >= rNonZeroWeightIndices[department_index].back() && t + 1 < numTeams_)
                {
                    for (int j = t + 1; j < numTeams_; j++)
                    {
                        rStageForwardMessages[j + 1] = initial_forward_messages;
                    }
                }
                else if (t + 1 < rNonZeroWeightIndices[department_index].back())
                {
                    upper             = upper_bound(rNonZeroWeightIndices[department_index].begin(),
                                                    rNonZeroWeightIndices[department_index].end(), t);
                    size_t next_index = upper - rNonZeroWeightIndices[department_index].begin();
                    for (int j = t + 1; j < rNonZeroWeightIndices[department_index][next_index]; j++)
                    {
                        rStageForwardMessages[j + 1] = initial_forward_messages;
                    }
                }
            }
        }
    }
}

/**
 * Perform backward pass for a department (trellis of knapsack) in Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] department_index: index of department (knapsack constraint)
 * @param[in] max_state_department: capacity of department under 
 * @param[in] rNonZeroWeightIndices: indices where connections are present between teams and departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of non-zero weight indices
 * @param[in] rTeamsWeightsPerDepartment: weights of teams for the department under investigation
 * @param[in] rTeam2KnapsackM: messages from teams to knapsack constraints
 * @param[out] rStageBackwardMessages: calculated backward messages on the trellis
 * 
 */

inline void Knapsack::backward(const int department_index, vector<vector<impalib_type>> &rStageBackwardMessages,
                        const int max_state_department, vector<vector<int>> &rNonZeroWeightIndices,
                        const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const vector<vector<int>> &rTeamsWeightsPerDepartment,
                        const vector<vector<impalib_type>> &rTeam2KnapsackM) const
{

    vector<int>::iterator upper;

    // Initialize initial backward messages
    vector<impalib_type>  initial_backward_messages(max_state_department + 1, zero_value);

    // Assign initial backward messages to the last stage
    rStageBackwardMessages[numTeams_] = initial_backward_messages;

    // If the last non-zero weight index is not numTeams_ - 1, assign initial messages to previous stages
    if (rNonZeroWeightIndices[department_index][pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] - 1]
        != numTeams_ - 1)
    {
        for (int j = numTeams_ - 1; j > rNonZeroWeightIndices[department_index][0]; j--)
        {
            rStageBackwardMessages[j] = initial_backward_messages;
        }
    }

    for (int l = pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] - 1; l >= 0; l--)
    {
        int t = rNonZeroWeightIndices[department_index][l];
        for (int a = 0; a <= max_state_department; a++)
        {
            if ((t > 0 && a + rTeamsWeightsPerDepartment[department_index][t] <= max_state_department)
                || (a == 0 && t == 0))
            {
                // Update backward messages
                rStageBackwardMessages[t][a] =
                    min(initial_backward_messages[a],
                        initial_backward_messages[a + rTeamsWeightsPerDepartment[department_index][t]]
                            + rTeam2KnapsackM[department_index][t]);
            }
            else
            {
                rStageBackwardMessages[t][a] = initial_backward_messages[a];
            }
        }

        // Update initial backward messages
        initial_backward_messages = rStageBackwardMessages[t];

        // Handle the case when non-zero weight indices are not equal to the number of teams
        if (pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] - 1 != numTeams_)
        {
            if (!binary_search(rNonZeroWeightIndices[department_index].begin(),
                               rNonZeroWeightIndices[department_index].end(), t - 1))
            {
                if ((t - 1 <= rNonZeroWeightIndices[department_index][0]) && t - 1 >= 0)
                {
                    for (int j = t - 1; j >= 0; j--)
                    {
                        rStageBackwardMessages[j] = initial_backward_messages;
                    }
                }
                else if (t - 1 > rNonZeroWeightIndices[department_index][0])
                {
                    upper             = upper_bound(rNonZeroWeightIndices[department_index].begin(),
                                                    rNonZeroWeightIndices[department_index].end(), t - 1);
                    size_t next_index = upper - rNonZeroWeightIndices[department_index].begin();
                    for (int j = t - 1; j > rNonZeroWeightIndices[department_index][next_index - 1]; j--)
                    {
                        rStageBackwardMessages[j] = initial_backward_messages;
                    }
                }
            }
        }
    }
}

/**
 * Calculate messages from departments to teams after obtaining forward and backward messages on a knapsack trellis following Forward-Backward algorithm for the Knapsack-MWM problem
 *
 * @param[in] rTeamsWeightsPerDepartment: weights of teams for the department under investigation
 * @param[in] rStageForwardMessages: calculated forward messages on the trellis
 * @param[in] rTeam2KnapsackM: messages from teams to knapsack constraints
 * @param[in] department_index: index of department (knapsack constraint)
 * @param[in] rStageBackwardMessages: calculated backward messages on the trellis
 * @param[in] max_state_department: capacity of department under 
 * @param[out] rExtrinsicOutputDepartment: extrinsic output messages from department to teams
 * 
 */

inline void Knapsack::extrinsic_output_department_lhs(const vector<vector<int>>          &rTeamsWeightsPerDepartment,
                                               const vector<vector<impalib_type>> &rStageForwardMessages,
                                               const vector<vector<impalib_type>> &rTeam2KnapsackM, const int department_index,
                                               const vector<vector<impalib_type>> &rStageBackwardMessages,
                                               const int                           max_state_department,
                                               vector<vector<impalib_type>> &rExtrinsicOutputDepartment)
{

    // Initialize metric paths
    vector<impalib_type> metric_path_solid, metric_path_dash;

    for (int team_index = 0; team_index < rTeamsWeightsPerDepartment[department_index].size(); team_index++)
    {
        metric_path_dash.clear();
        metric_path_solid.clear();

        // Check if the team has non-zero weight with the department
        if (rTeamsWeightsPerDepartment[department_index][team_index] != 0)
        {
            // Calculate metric paths for solid (activated) and dashed (deactivated) states
            if (team_index == 0)
            {
                metric_path_solid.push_back(
                    rStageForwardMessages[team_index][team_index]
                    + rStageBackwardMessages[team_index + 1][rTeamsWeightsPerDepartment[department_index][team_index]]
                    + rTeam2KnapsackM[department_index][team_index]);
                metric_path_dash.push_back(rStageForwardMessages[team_index][team_index]
                                           + rStageBackwardMessages[team_index + 1][0]);
            }
            else
            {
                for (int i = 0; i <= max_state_department; i++)
                {
                    if (i + rTeamsWeightsPerDepartment[department_index][team_index] <= max_state_department)
                    {
                        metric_path_solid.push_back(
                            rStageForwardMessages[team_index][i]
                            + rStageBackwardMessages[team_index + 1]
                                                    [i + rTeamsWeightsPerDepartment[department_index][team_index]]
                            + rTeam2KnapsackM[department_index][team_index]);
                    }
                    metric_path_dash.push_back(rStageForwardMessages[team_index][i]
                                               + rStageBackwardMessages[team_index + 1][i]);
                }
            }
            // Calculate extrinsic output for the department
            rExtrinsicOutputDepartment[department_index][team_index] =
                *min_element(metric_path_solid.begin(), metric_path_solid.end())
                - *min_element(metric_path_dash.begin(), metric_path_dash.end())
                - rTeam2KnapsackM[department_index][team_index];
        }
        else
        {
            rExtrinsicOutputDepartment[department_index][team_index] = zero_value;
        }
    }
}

/**
 * Calculate messages from teams to departments for the Knapsack-MWM problem
 *
 * @param[in] rNonZeroWeightIndices: indices where connections are present between teams and departments
 * @param[in] rTeam2KnapsackM: messages from teams to knapsack constraints
 * @param[in] rRewardTeam: rewards of teams
 * @param[in] rExtrinsicOutputDepartment: extrinsic output messages from department to teams
 * @param[out] mORIC2Team: messages from ORIC to teams
 * 
 */

inline void Knapsack::team_to_knapsack_update(vector<vector<int>>          &rNonZeroWeightIndices,
                                       vector<vector<impalib_type>> &rTeam2KnapsackM, const vector<impalib_type> &rRewardTeam,
                                       const vector<vector<impalib_type>> &rExtrinsicOutputDepartment,
                                       const vector<impalib_type>         &mORIC2Team) const
{
    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++)
    {
        // Create a list of remaining departments
        vector<int> remaining_departments(numDepartments_);
        iota(remaining_departments.begin(), remaining_departments.end(), 0);

        // Get unique edge indices for the current department
        vector<int> unique_edge_department(rNonZeroWeightIndices[department_index]);
        remaining_departments.erase(remaining_departments.begin() + department_index);

        for (auto t : remaining_departments)
        {
            // Calculate intersection of non-zero weight indices between the current department and other departments
            vector<int>           intersection(numTeams_);
            vector<int>::iterator it;
            it = set_intersection(rNonZeroWeightIndices[department_index].begin(),
                                  rNonZeroWeightIndices[department_index].end(), rNonZeroWeightIndices[t].begin(),
                                  rNonZeroWeightIndices[t].end(), intersection.begin());
            intersection.resize(it - intersection.begin());

            // Update team to knapsack constraint messages
            for (auto l : intersection)
            {

                rTeam2KnapsackM[department_index][l] =
                    rRewardTeam[l] + rExtrinsicOutputDepartment[t][l] + mORIC2Team[l];
                
                // Remove the edge index from the unique list
                vector<int>::iterator position =
                    std::find(unique_edge_department.begin(), unique_edge_department.end(), l);
                unique_edge_department.erase(position);
            }
        }

        // Update team to knapsack constraint messages for remaining unique edge department indices
        for (auto u : unique_edge_department)
        {
            rTeam2KnapsackM[department_index][u] = rRewardTeam[u] + mORIC2Team[u];
        }
    }
}

/**
 * Calculate messages from departments to teams based on filtering conditions for the Knapsack-MWM problem
 *
 * @param[in] department_index: index of department (knapsack constraint)
 * @param[in] iter: iteration index
 * @param[in] rExtrinsicOutputDepartmentDummy: extrinsic output messages from department to teams before filtering
 * @param[out] rExtrinsicOutputDepartment: extrinsic output messages from department to teams after filtering
 * 
 */

inline void Knapsack::process_extrinsic_output_department(const int department_index, const int iter,
                                                   vector<vector<impalib_type>> &rExtrinsicOutputDepartmentDummy,
                                                   vector<vector<impalib_type>> &rExtrinsicOutputDepartment)
{

    if ((filteringFlag_) and (alpha_ != zero_value))
    {

        vector<impalib_type> intermediate_dummy(rExtrinsicOutputDepartmentDummy[department_index]),
            intermediate_old(extrinsicOutputDepartmentOld_[department_index]), intermediate_extrinsic;

        impalib_type w_1 = alpha_, w_2 = 1 - alpha_;

        // Calculate weighted extrinsic outputs
        transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                  [w_2](const impalib_type &c) { return c * w_2; });
        transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                  [w_1](const impalib_type &c) { return c * w_1; });

        if (iter == 0)
        {
            copy(intermediate_dummy.begin(), intermediate_dummy.end(),
                 rExtrinsicOutputDepartment[department_index].begin());
        }
        else
        {
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(),
                      std::back_inserter(intermediate_extrinsic), std::plus<impalib_type>());
            copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(),
                 rExtrinsicOutputDepartment[department_index].begin());
        }
        
        copy(rExtrinsicOutputDepartment[department_index].begin(), rExtrinsicOutputDepartment[department_index].end(),
             extrinsicOutputDepartmentOld_[department_index].begin());
    }

    else
    {
        copy(rExtrinsicOutputDepartmentDummy[department_index].begin(),
             rExtrinsicOutputDepartmentDummy[department_index].end(),
             rExtrinsicOutputDepartment[department_index].begin());
    }
}