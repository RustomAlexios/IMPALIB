// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for the degree constraint for the TSP
 */
class DegreeConstraint
{
private:
    int                          numNodes_; ///< number of nodes of TSP
    int                          numEdgeVariables_; ///< number of edge connections
    bool                         filteringFlag_; ///< filtering flag
    impalib_type                 alpha_; ///< filtering parameter
    vector<vector<impalib_type>> degreeConstraint2EqConstraintOld; ///< messages from degree constraints to equality constraints before filtering
    impalib_type                 initial_forward_message_; ///< initial forward message of forward-backward algorithm
    impalib_type                 initial_backward_message_; ///< initial backward message of forward-backward algorithm
    vector<vector<impalib_type>> stage_forward_messages; ///< forward messages of trellis representation
    vector<vector<impalib_type>> stage_backward_messages; ///< backward messages of trellis representation

public:
    void degree_constraint_to_edge_ec_update(vector<vector<impalib_type>> &, vector<vector<int>> &,
                                             vector<vector<impalib_type>> &); ///< calculate messages from degree constraint to edge equality constraint
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &); ///< process filtering on messages from degree constraint to edge equality constraint

    DegreeConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                     const impalib_type ALPHA); ///< constructor
};

/**
 * Construct DegreeConstraint object for the Knapsack-MWM problem
 *
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of connections between nodes
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[out] filteringFlag_: FILTERING_FLAG
 * @param[out] alpha_: ALPHA
 * @param[out] numNodes_: NUM_NODES
 * @param[out] numEdgeVariables_: NUM_EDGE_VARIABLES
 */

DegreeConstraint::DegreeConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                                   const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES)
{

    degreeConstraint2EqConstraintOld.reserve(numEdgeVariables_);
    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++)
    {
        degreeConstraint2EqConstraintOld.push_back(vector<impalib_type>(numNodes_, zero_value));
    }

    // Set initial forward and backward messages to infinity
    initial_forward_message_  = value_inf;
    initial_backward_message_ = value_inf;
};

/**
 * Calculate messages from degree constraints to edge equality constraints for the TSP
 *
 * @param[in] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[out] rDegreeConstraint2EqConstraintDummyM: messages from degree constraints to edge equality constraints before filtering
 * 
 */

void DegreeConstraint::degree_constraint_to_edge_ec_update(
    vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM, vector<vector<int>> &rEdgeConnections,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintDummyM)
{
    vector<impalib_type> stage_forward_messages(numEdgeVariables_ + 1, zero_value);
    vector<impalib_type> stage_backward_messages(numEdgeVariables_ + 1, zero_value);

    for (int i = 0; i < numNodes_; i++)
    {
        vector<int> connections_first, connections_second;

        // Populate connections_first and connections_second with edge indices
        for (int j = 0; j < rEdgeConnections.size(); j++)
        {
            if (i == rEdgeConnections[j][0])
            {
                connections_first.push_back(j);
            }
            else if (i == rEdgeConnections[j][1])
            {
                connections_second.push_back(j);
            }
        }

        // Calculate forward messages for connections_first
        stage_forward_messages[connections_first[0]] = initial_forward_message_;

        for (int s = 1; s < connections_first.size(); s++)
        {

            stage_forward_messages[connections_first[s]] =
                min(stage_forward_messages[connections_first[s - 1]],
                    rEdgeEc2DegreeConstraintM[connections_first[s - 1]][i]);
        }

        // Calculate backward messages for connections_first
        stage_backward_messages[connections_first[connections_first.size() - 1] + 1] = initial_backward_message_;

        for (size_t s = connections_first.size() - 1; s >= 1; s--)
        {
            stage_backward_messages[connections_first[s - 1] + 1] =
                min(stage_backward_messages[connections_first[s] + 1],
                    rEdgeEc2DegreeConstraintM[connections_first[s]][i]);
        }

        for (int j = 0; j < connections_first.size(); j++)
        {
            impalib_type minimumValue = zero_value;
            minimumValue              = min(stage_forward_messages[connections_first[j]],
                                            stage_backward_messages[connections_first[j] + 1]);
            rDegreeConstraint2EqConstraintDummyM[connections_first[j]][i] = -minimumValue;
        }

        fill(stage_forward_messages.begin(), stage_forward_messages.end(), zero_value);
        fill(stage_backward_messages.begin(), stage_backward_messages.end(), zero_value);

        // Calculate forward messages
        stage_forward_messages[connections_second[0]] = initial_forward_message_;

        for (int s = 1; s < connections_second.size(); s++)
        {

            stage_forward_messages[connections_second[s]] =
                min(stage_forward_messages[connections_second[s - 1]],
                    rEdgeEc2DegreeConstraintM[connections_second[s - 1]][i]);
        }

        // Calculate backward messages
        stage_backward_messages[connections_second[connections_second.size() - 1] + 1] = initial_backward_message_;

        for (size_t s = connections_second.size() - 1; s >= 1; s--)
        {
            stage_backward_messages[connections_second[s - 1] + 1] =
                min(stage_backward_messages[connections_second[s] + 1],
                    rEdgeEc2DegreeConstraintM[connections_second[s]][i]);
        }

        for (int j = 0; j < connections_second.size(); j++)
        {
            impalib_type minimumValue = zero_value;
            minimumValue              = min(stage_forward_messages[connections_second[j]],
                                            stage_backward_messages[connections_second[j] + 1]);
            rDegreeConstraint2EqConstraintDummyM[connections_second[j]][i] = -minimumValue;
        }
    }
}

/**
 * Perform filtering on messages from degree constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: iteration index of IMPA
 * @param[in] rDegreeConstraint2EqConstraintDummyM: messages from degree constraints to edge equality constraints before filtering
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to edge equality constraints after filtering
 * 
 */

void DegreeConstraint::process_filtering(int iter, vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintDummyM,
                                         vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM)
{
    for (int i = 0; i < numEdgeVariables_; i++)
    {
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            // Calculate weighted values for current and old messages
            vector<impalib_type> intermediate_dummy(rDegreeConstraint2EqConstraintDummyM[i]),
                intermediate_old(degreeConstraint2EqConstraintOld[i]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                      [w_2](impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                      [w_1](impalib_type &c) { return c * w_1; });

            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(),
                     rDegreeConstraint2EqConstraintM[i].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(),
                          back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(),
                     rDegreeConstraint2EqConstraintM[i].begin());
            }
            copy(rDegreeConstraint2EqConstraintM[i].begin(),
                 rDegreeConstraint2EqConstraintM[i].end(),
                 degreeConstraint2EqConstraintOld[i].begin());
        }

        else
        {
            copy(rDegreeConstraint2EqConstraintDummyM[i].begin(),
                 rDegreeConstraint2EqConstraintDummyM[i].end(),
                 rDegreeConstraint2EqConstraintM[i].begin());
        }
    }
}
