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
    void subtour_constraints_to_edge_ec_update(const vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                               vector<vector<impalib_type>> &) const; ///< calculate messages from subtour to edge equality constraints
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, const vector<vector<int>> &); ///< perform filtering on messages from subtour to edge equality constraints

    SubtourEliminationConstraint(int NUM_NODES, int NUM_EDGE_VARIABLES, bool FILTERING_FLAG,
                                 impalib_type ALPHA); ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of connections between nodes
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
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
    const vector<vector<impalib_type>> &rEdgeEc2SubtourConstraintsM, const vector<vector<int>> &rDeltaSIndicesList,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM) const
{
    vector<impalib_type> stage_forward_messages(numEdgeVariables_ + 1, zero_value);
    vector<impalib_type> stage_backward_messages(numEdgeVariables_ + 1, zero_value);

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

void SubtourEliminationConstraint::process_filtering(const int                           iter,
                                                     vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcDummyM,
                                                     vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM,
                                                     const vector<vector<int>>          &rDeltaSIndicesList)
{

    for (int index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
         index_subtour_constraint++)
    {

        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            vector<impalib_type> intermediate_dummy(rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint]),
                intermediate_old(subtourConstraints2EdgeEcOld_[index_subtour_constraint]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                      [w_2](const impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                      [w_1](const impalib_type &c) { return c * w_1; });

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

            copy(rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin(),
                 rSubtourConstraints2EdgeEcM[index_subtour_constraint].end(),
                 subtourConstraints2EdgeEcOld_[index_subtour_constraint].begin());
        }

        else
        {
            copy(rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint].begin(),
                 rSubtourConstraints2EdgeEcDummyM[index_subtour_constraint].end(),
                 rSubtourConstraints2EdgeEcM[index_subtour_constraint].begin());
        }
    }
}