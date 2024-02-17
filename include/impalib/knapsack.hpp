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
    void forward(int, vector<vector<impalib_type>> &, int, vector<int> &, const int *, vector<int> &,
                 vector<vector<impalib_type>> &); ///< forward pass of forward-backward algorithm

    void backward(int, vector<vector<impalib_type>> &, int, vector<int> &, const int *, vector<int> &,
                  vector<vector<impalib_type>> &); ///< backward pass of forward-backward algorithm

    vector<impalib_type> extrinsic_output_department_lhs(vector<int> &, vector<vector<impalib_type>> &,
                                         vector<vector<impalib_type>> &, int, vector<vector<impalib_type>> &, int); ///< extrinsic output of department constraint

    void team_to_knapsack_update(vector<vector<int>> &, vector<vector<impalib_type>> &, vector<impalib_type> &,
                                 vector<vector<impalib_type>> &, vector<impalib_type> &); ///< calculate messages from teams to knapsack constraints
    void process_extrinsic_output_department(int, int, vector<impalib_type> &, vector<vector<impalib_type>> &); ///< perform filtering (if needed) on messages from departments to teams

    Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA); ///< constructor
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
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), filteringFlag_(FILT_FLAG), alpha_(ALPHA),
      extrinsicOutputDepartmentOld_(numDepartments_, vector<impalib_type>(numTeams_, 0))
{
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

void Knapsack::forward(int department_index, vector<vector<impalib_type>> &rStageForwardMessages,
                       int max_state_department, vector<int> &rNonZeroWeightIndices,
                       const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, vector<int> &rTeamWeights,
                       vector<vector<impalib_type>> &rTeam2KnapsackM)
{
    vector<int>::iterator upper;

    // Initialize initial forward messages
    vector<impalib_type> initial_forward_messages(max_state_department + 1, zero_value);
    fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

    // Assign initial forward messages to the first stage
    rStageForwardMessages[0] = initial_forward_messages;

    // If the first non-zero weight index is not zero, assign initial messages to subsequent stages
    if (rNonZeroWeightIndices[0] != 0)
    {
        for (int j = 0; j < rNonZeroWeightIndices[0]; j++)
        {
            rStageForwardMessages[j + 1] = initial_forward_messages;
        }
    }

    for (int l = 0; l < pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]; l++)
    {
        int t = rNonZeroWeightIndices[l];
        for (int a = 0; a <= max_state_department; a++)
        {
            if (a - rTeamWeights[t] >= 0
                && (!(t == 0 && a != rTeamWeights[t])))
            {
                // Update forward messages
                rStageForwardMessages[t + 1][a] =
                    min(initial_forward_messages[a],
                        initial_forward_messages[a - rTeamWeights[t]]
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
        if (rNonZeroWeightIndices.size() != numTeams_)
        {
            if (!binary_search(rNonZeroWeightIndices.begin(),
                               rNonZeroWeightIndices.end(), t + 1))
            {
                if (t + 1 >= rNonZeroWeightIndices.back() && t + 1 < numTeams_)
                {
                    for (int j = t + 1; j < numTeams_; j++)
                    {
                        rStageForwardMessages[j + 1] = initial_forward_messages;
                    }
                }
                else if (t + 1 < rNonZeroWeightIndices.back())
                {
                    upper             = upper_bound(rNonZeroWeightIndices.begin(),
                                                    rNonZeroWeightIndices.end(), t);
                    size_t next_index = upper - rNonZeroWeightIndices.begin();
                    for (int j = t + 1; j < rNonZeroWeightIndices[next_index]; j++)
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

void Knapsack::backward(int department_index, vector<vector<impalib_type>> &rStageBackwardMessages,
                        int max_state_department, vector<int> &rNonZeroWeightIndices,
                        const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, vector<int> &rTeamsWeightsPerDepartment,
                        vector<vector<impalib_type>> &rTeam2KnapsackM)
{

    vector<int>::iterator upper;

    // Initialize initial backward messages
    vector<impalib_type>  initial_backward_messages(max_state_department + 1, zero_value);

    // Assign initial backward messages to the last stage
    rStageBackwardMessages[numTeams_] = initial_backward_messages;

    // If the last non-zero weight index is not numTeams_ - 1, assign initial messages to previous stages
    if (rNonZeroWeightIndices[pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] - 1]
        != numTeams_ - 1)
    {
        for (int j = numTeams_ - 1; j > rNonZeroWeightIndices[0]; j--)
        {
            rStageBackwardMessages[j] = initial_backward_messages;
        }
    }

    for (int l = pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] - 1; l >= 0; l--)
    {
        int t = rNonZeroWeightIndices[l];
        for (int a = 0; a <= max_state_department; a++)
        {
            if ((t > 0 && a + rTeamsWeightsPerDepartment[t] <= max_state_department)
                || (a == 0 && t == 0))
            {
                // Update backward messages
                rStageBackwardMessages[t][a] =
                    min(initial_backward_messages[a],
                        initial_backward_messages[a + rTeamsWeightsPerDepartment[t]]
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
            if (!binary_search(rNonZeroWeightIndices.begin(),
                               rNonZeroWeightIndices.end(), t - 1))
            {
                if ((t - 1 <= rNonZeroWeightIndices[0]) && t - 1 >= 0)
                {
                    for (int j = t - 1; j >= 0; j--)
                    {
                        rStageBackwardMessages[j] = initial_backward_messages;
                    }
                }
                else if (t - 1 > rNonZeroWeightIndices[0])
                {
                    upper             = upper_bound(rNonZeroWeightIndices.begin(),
                                                    rNonZeroWeightIndices.end(), t - 1);
                    size_t next_index = upper - rNonZeroWeightIndices.begin();
                    for (int j = t - 1; j > rNonZeroWeightIndices[next_index - 1]; j--)
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

vector<impalib_type> Knapsack::extrinsic_output_department_lhs(vector<int>          &rTeamsWeightsPerDepartment,
                                               vector<vector<impalib_type>> &rStageForwardMessages,
                                               vector<vector<impalib_type>> &rTeam2KnapsackM, int department_index,
                                               vector<vector<impalib_type>> &rStageBackwardMessages,
                                               int                           max_state_department)
{
    vector<impalib_type> rExtrinsicOutput(rTeamsWeightsPerDepartment.size());
    vector<impalib_type> metric_path_solid, metric_path_dash;

    for (int i = 0; i < rTeamsWeightsPerDepartment.size(); i++)
    {
        metric_path_dash.clear();
        metric_path_solid.clear();

        // Check if the team has non-zero weight with the department
        if (rTeamsWeightsPerDepartment[i] != 0)
        {
            // Calculate metric paths for solid (activated) and dashed (deactivated) states
            if (i == 0)
            {
                metric_path_solid.push_back(
                    rStageForwardMessages[i][i]
                    + rStageBackwardMessages[i + 1][rTeamsWeightsPerDepartment[i]]
                    + rTeam2KnapsackM[department_index][i]);
                metric_path_dash.push_back(rStageForwardMessages[i][i]
                                           + rStageBackwardMessages[i + 1][0]);
            }
            else
            {
                for (int j = 0; j <= max_state_department; j++)
                {
                    if (j + rTeamsWeightsPerDepartment[i] <= max_state_department)
                    {
                        metric_path_solid.push_back(
                            rStageForwardMessages[i][j]
                            + rStageBackwardMessages[i + 1]
                                                    [j + rTeamsWeightsPerDepartment[i]]
                            + rTeam2KnapsackM[department_index][i]);
                    }
                    metric_path_dash.push_back(rStageForwardMessages[i][j]
                                               + rStageBackwardMessages[i + 1][j]);
                }
            }
            // Calculate extrinsic output for the department
            rExtrinsicOutput[i] =
                *min_element(metric_path_solid.begin(), metric_path_solid.end())
                - *min_element(metric_path_dash.begin(), metric_path_dash.end())
                - rTeam2KnapsackM[department_index][i];
        }
        else
        {
            rExtrinsicOutput[i] = zero_value;
        }
    }
    return rExtrinsicOutput;
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

void Knapsack::team_to_knapsack_update(vector<vector<int>>          &rNonZeroWeightIndices,
                                       vector<vector<impalib_type>> &rTeam2KnapsackM, vector<impalib_type> &rRewardTeam,
                                       vector<vector<impalib_type>> &rExtrinsicOutputDepartment,
                                       vector<impalib_type>         &mORIC2Team)
{
    for (int i = 0; i < rExtrinsicOutputDepartment.size(); i++)
    {
        // Create a list of remaining departments
        vector<int> remaining_departments(numDepartments_);
        iota(remaining_departments.begin(), remaining_departments.end(), 0);

        // Get unique edge indices for the current department
        vector<int> unique_edge_department(rNonZeroWeightIndices[i]);
        remaining_departments.erase(remaining_departments.begin() + i);

        for (auto t : remaining_departments)
        {
            // Calculate intersection of non-zero weight indices between the current department and other departments
            vector<int>           intersection(numTeams_);
            vector<int>::iterator it;
            it = set_intersection(rNonZeroWeightIndices[i].begin(),
                                  rNonZeroWeightIndices[i].end(), rNonZeroWeightIndices[t].begin(),
                                  rNonZeroWeightIndices[t].end(), intersection.begin());
            intersection.resize(it - intersection.begin());

            // Update team to knapsack constraint messages
            for (auto l : intersection)
            {

                rTeam2KnapsackM[i][l] =
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
            rTeam2KnapsackM[i][u] = rRewardTeam[u] + mORIC2Team[u];
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

void Knapsack::process_extrinsic_output_department(int department_index, int iter,
                                                   vector<impalib_type> &rExtrinsicOutput,
                                                   vector<vector<impalib_type>> &rExtrinsicOutputDepartment)
{

    if ((filteringFlag_) and (alpha_ != zero_value))
    {
        auto intermediate_dummy = rExtrinsicOutput;
        auto intermediate_old = extrinsicOutputDepartmentOld_[department_index];
        vector<impalib_type> intermediate_extrinsic;

        impalib_type w_1 = alpha_, w_2 = 1 - alpha_;

        // Calculate weighted extrinsic outputs
        transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                  [w_2](impalib_type &c) { return c * w_2; });
        transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                  [w_1](impalib_type &c) { return c * w_1; });

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
        copy(rExtrinsicOutput.begin(),
             rExtrinsicOutput.end(),
             rExtrinsicOutputDepartment[department_index].begin());
    }
}