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
    void degree_constraint_to_edge_ec_update(const vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                             vector<vector<impalib_type>> &) const; ///< calculate messages from degree constraint to edge equality constraint
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &); ///< process filtering on messages from degree constraint to edge equality constraint

    DegreeConstraint(int NUM_NODES, int NUM_EDGE_VARIABLES, bool FILTERING_FLAG,
                     impalib_type ALPHA); ///< constructor
};

/**
 * Construct DegreeConstraint object for the Knapsack-MWM problem
 *
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of connections between nodes
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 */

inline DegreeConstraint::DegreeConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                                   const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES)
{

    degreeConstraint2EqConstraintOld.reserve(numEdgeVariables_);
    for (int edge = 0; edge < numEdgeVariables_; edge++)
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
 * @param[in] eq2DegM: messages from edge equality constraints to degree constraints
 * @param[in] connections: list of connections for each edge equality constraint
 * @param[out] deg2EqPreM: messages from degree constraints to edge equality constraints before filtering
 * 
 */

inline void DegreeConstraint::degree_constraint_to_edge_ec_update(
    const vector<vector<impalib_type>> &eq2DegM, const vector<vector<int>> &connections,
    vector<vector<impalib_type>> &deg2EqPreM) const
{
    vector<impalib_type> forward(numEdgeVariables_ + 1, zero_value);
    vector<impalib_type> backward(numEdgeVariables_ + 1, zero_value);

    for (int node = 0; node < numNodes_; node++)
    {
        vector<int> first, second;

        // Populate connections_first and connections_second with edge indices
        for (int i = 0; i < connections.size(); i++)
        {
            if (node == connections[i][0])
            {
                first.push_back(i);
            }
            else if (node == connections[i][1])
            {
                second.push_back(i);
            }
        }

        // Calculate forward messages for connections_first
        forward[first[0]] = initial_forward_message_;

        for (int stage = 1; stage < first.size(); stage++)
        {

            forward[first[stage]] =
                min(forward[first[stage - 1]],
                    eq2DegM[first[stage - 1]][node]);
        }

        // Calculate backward messages for connections_first
        backward[first[first.size() - 1] + 1] = initial_backward_message_;

        for (size_t stage = first.size() - 1; stage >= 1; stage--)
        {
            backward[first[stage - 1] + 1] =
                min(backward[first[stage] + 1],
                    eq2DegM[first[stage]][node]);
        }

        for (int edge = 0; edge < first.size(); edge++)
        {
            impalib_type minimumValue = zero_value;
            minimumValue              = min(forward[first[edge]],
                                            backward[first[edge] + 1]);
            deg2EqPreM[first[edge]][node] = -minimumValue;
        }

        fill(forward.begin(), forward.end(), zero_value);
        fill(backward.begin(), backward.end(), zero_value);

        // Calculate forward messages
        forward[second[0]] = initial_forward_message_;

        for (int stage = 1; stage < second.size(); stage++)
        {

            forward[second[stage]] =
                min(forward[second[stage - 1]],
                    eq2DegM[second[stage - 1]][node]);
        }

        // Calculate backward messages
        backward[second[second.size() - 1] + 1] = initial_backward_message_;

        for (size_t stage = second.size() - 1; stage >= 1; stage--)
        {
            backward[second[stage - 1] + 1] =
                min(backward[second[stage] + 1],
                    eq2DegM[second[stage]][node]);
        }

        for (int edge = 0; edge < second.size(); edge++)
        {
            impalib_type minimumValue = zero_value;
            minimumValue              = min(forward[second[edge]],
                                            backward[second[edge] + 1]);
            deg2EqPreM[second[edge]][node] = -minimumValue;
        }
    }
}

/**
 * Perform filtering on messages from degree constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: iteration index of IMPA
 * @param[in] deg2EqPreM: messages from degree constraints to edge equality constraints before filtering
 * @param[in] deg2EqM: messages from degree constraints to edge equality constraints after filtering
 * 
 */

inline void DegreeConstraint::process_filtering(const int iter, vector<vector<impalib_type>> &deg2EqPreM,
                                         vector<vector<impalib_type>> &deg2EqM)
{
    for (int edge = 0; edge < numEdgeVariables_; edge++)
    {
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            // Calculate weighted values for current and old messages
            vector<impalib_type> intermediate_dummy(deg2EqPreM[edge]);
            vector<impalib_type> intermediate_old(degreeConstraint2EqConstraintOld[edge]);
            vector<impalib_type> intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                      [w_2](const impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                      [w_1](const impalib_type &c) { return c * w_1; });

            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(),
                     deg2EqM[edge].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(),
                          back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(),
                     deg2EqM[edge].begin());
            }
            copy(deg2EqM[edge].begin(),
                 deg2EqM[edge].end(),
                 degreeConstraint2EqConstraintOld[edge].begin());
        }

        else
        {
            copy(deg2EqPreM[edge].begin(),
                 deg2EqPreM[edge].end(),
                 deg2EqM[edge].begin());
        }
    }
}
