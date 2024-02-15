// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for the subtour elimination constraint for the TSP
 */

class SubtourEliminationConstraint
{
private:
    int          numNodes_; ///< number of nodes
    int          numEdgeVariables_; ///< number of edges
    bool         filteringFlag_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter
    impalib_type initial_forward_message_; ///< initial value of forward message using forward-backward algorithm
    impalib_type initial_backward_message_; ///< initial value of backward message using forward-backward algorithm

public:
    vector<vector<impalib_type>> subtourConstraints2EdgeEcOld_; ///< messages from subtour constraints to edge equality constraint before filtering
    void subtour_constraints_to_edge_ec_update(vector<vector<impalib_type>> &, vector<vector<int>> &,
                                               vector<vector<impalib_type>> &); ///< calculate messages from subtour to edge equality constraints
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &); ///< perform filtering on messages from subtour to edge equality constraints

    SubtourEliminationConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                                 const impalib_type ALPHA); ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of connections between nodes
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[out] filteringFlag_: FILTERING_FLAG
 * @param[out] alpha_: ALPHA
 * @param[out] numNodes_: NUM_NODES
 * @param[out] numEdgeVariables_: NUM_EDGE_VARIABLES
 * @param[out] initial_forward_message_: set to infinity
 * @param[out] initial_backward_message_: set to infinity
 */
SubtourEliminationConstraint::SubtourEliminationConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                                           const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES),
      initial_forward_message_(value_inf), initial_backward_message_(value_inf){
                                           };

/**
 * Calculate messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] rEdgeEc2SubtourConstraintsM: messages from edge equality constraints to subtour elimination constraints
 * @param[in] rDeltaSIndicesList: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 * @param[out] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints
 * 
 */

void SubtourEliminationConstraint::subtour_constraints_to_edge_ec_update(
    vector<vector<impalib_type>> &rEdgeEc2SubtourConstraintsM, vector<vector<int>> &rDeltaSIndicesList,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM)
{
    // Define containers for forward and backward messages
    vector<impalib_type> stage_forward_messages(numEdgeVariables_ + 1, zero_value);
    vector<impalib_type> stage_backward_messages(numEdgeVariables_ + 1, zero_value);

    // Iterate over subtour constraints
    for (size_t index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
         index_subtour_constraint++)
    {
        // Initialize and calculate forward messages
        stage_forward_messages[rDeltaSIndicesList[index_subtour_constraint][0]] = initial_forward_message_;
        for (int stage = 1; stage < rDeltaSIndicesList[index_subtour_constraint].size(); stage++)
        {
            stage_forward_messages[rDeltaSIndicesList[index_subtour_constraint][stage]] =
                min(stage_forward_messages[rDeltaSIndicesList[index_subtour_constraint][stage - 1]],
                    rEdgeEc2SubtourConstraintsM[index_subtour_constraint]
                                               [rDeltaSIndicesList[index_subtour_constraint][stage - 1]]);
        }

        // Initialize and calculate backward messages
        stage_backward_messages[rDeltaSIndicesList[index_subtour_constraint]
                                                  [rDeltaSIndicesList[index_subtour_constraint].size() - 1]
                                + 1] = initial_backward_message_;

        for (size_t stage = rDeltaSIndicesList[index_subtour_constraint].size() - 1; stage >= 1; stage--)
        {
            stage_backward_messages[rDeltaSIndicesList[index_subtour_constraint][stage - 1] + 1] =
                min(stage_backward_messages[rDeltaSIndicesList[index_subtour_constraint][stage] + 1],
                    rEdgeEc2SubtourConstraintsM[index_subtour_constraint]
                                               [rDeltaSIndicesList[index_subtour_constraint][stage]]);
        }

        // Update subtour constraints to edge EC messages
        for (int index_edge_variable = 0; index_edge_variable < rDeltaSIndicesList[index_subtour_constraint].size();
             index_edge_variable++)
        {
            impalib_type minimumValue = zero_value;
            minimumValue =
                min(stage_forward_messages[rDeltaSIndicesList[index_subtour_constraint][index_edge_variable]],
                    stage_backward_messages[rDeltaSIndicesList[index_subtour_constraint][index_edge_variable] + 1]);
            minimumValue = min(-minimumValue, zero_value);
            rSubtourConstraints2EdgeEcM[index_subtour_constraint]
                                       [rDeltaSIndicesList[index_subtour_constraint][index_edge_variable]] =
                                           minimumValue;
        }
    }
}

/**
 * Process filtering on messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: index of iteration of IMPA
 * @param[in] rSubtourConstraints2EdgeEcDummyM: messages from subtour elimination constraints to edge equality constraints before filtering
 * @param[out] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints after filtering
 * @param[in] rDeltaSIndicesList: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 * 
 */

void SubtourEliminationConstraint::process_filtering(int                           iter,
                                                     vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcDummyM,
                                                     vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM,
                                                     vector<vector<int>>          &rDeltaSIndicesList)
{

    // Iterate over subtour constraints
    for (int index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
         index_subtour_constraint++)
    {

        // Apply filtering if enabled and alpha is not zero
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            // Calculate intermediate values based on alpha
            vector<impalib_type> intermediate_dummy(rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint]),
                intermediate_old(subtourConstraints2EdgeEcOld_[index_subtour_constraint]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                      [w_2](impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                      [w_1](impalib_type &c) { return c * w_1; });
            
            // Update the messages based on iteration
            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(),
                     rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(),
                          back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(),
                     rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin());
            }

            // Update the old messages for the next iteration
            copy(rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin(),
                 rSubtourConstraints2EdgeEcM[index_subtour_constraint].end(),
                 subtourConstraints2EdgeEcOld_[index_subtour_constraint].begin());
        }

        // If filtering is not enabled or alpha is zero, simply copy the dummy values
        else
        {
            copy(rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint].begin(),
                 rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint].end(),
                 rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin());
        }
    }
}