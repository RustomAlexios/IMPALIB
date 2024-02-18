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

class SubtourEliminationConstraint {
   private:
    int n_Nodes;             ///< number of nodes
    int n_Edges;             ///< number of edges
    bool doFilter;           ///< filtering flag
    impalib_type alpha_;     ///< filtering parameter
    impalib_type forward0;   ///< initial value of forward message using forward-backward algorithm
    impalib_type backward0;  ///< initial value of backward message using forward-backward algorithm

   public:
    vector<vector<impalib_type>> M_subtour2edge_old;  ///< messages from subtour constraints to edge equality constraint before filtering
    void messages_to_edge_ec(const vector<vector<impalib_type>> &edge2subtour, const vector<vector<int>> &deltaS,
                             vector<vector<impalib_type>> &subtour2edge);  ///< calculate messages from subtour to edge equality constraints
    void process_filtering(int, const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                           const vector<vector<int>> &);  ///< perform filtering on messages from subtour to edge equality constraints

    SubtourEliminationConstraint(int NUM_NODES, int NUM_EDGE_VARIABLES, bool FILTERING_FLAG,
                                 impalib_type ALPHA);  ///< constructor
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
SubtourEliminationConstraint::SubtourEliminationConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG, const impalib_type ALPHA)
    : doFilter(FILTERING_FLAG), alpha_(ALPHA), n_Nodes(NUM_NODES), n_Edges(NUM_EDGE_VARIABLES), forward0(value_inf), backward0(value_inf){};

/**
 * Calculate messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] edge2subtour: messages from edge equality constraints to subtour elimination constraints
 * @param[in] deltaS: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 * @param[out] subtour2edge: messages from subtour elimination constraints to edge equality constraints
 *
 */

void SubtourEliminationConstraint::messages_to_edge_ec(const vector<vector<impalib_type>> &edge2subtour, const vector<vector<int>> &deltaS, vector<vector<impalib_type>> &subtour2edge) {
    vector<impalib_type> forward(n_Edges + 1, zero_value);
    vector<impalib_type> backward(n_Edges + 1, zero_value);

    for (size_t i = 0; i < deltaS.size(); i++) {
        // Initialize and calculate forward messages
        forward[deltaS[i][0]] = forward0;
        for (int stage = 1; stage < deltaS[i].size(); stage++) {
            forward[deltaS[i][stage]] = min(forward[deltaS[i][stage - 1]], edge2subtour[i][deltaS[i][stage - 1]]);
        }

        // Initialize and calculate backward messages
        backward[deltaS[i][deltaS[i].size() - 1] + 1] = backward0;

        for (size_t stage = deltaS[i].size() - 1; stage >= 1; stage--) {
            backward[deltaS[i][stage - 1] + 1] = min(backward[deltaS[i][stage] + 1], edge2subtour[i][deltaS[i][stage]]);
        }

        // Update subtour constraints to edge EC messages
        for (int idx_edge = 0; idx_edge < deltaS[i].size(); idx_edge++) {
            impalib_type min_val = zero_value;
            min_val = min(forward[deltaS[i][idx_edge]], backward[deltaS[i][idx_edge] + 1]);
            min_val = min(-min_val, zero_value);
            subtour2edge[i][deltaS[i][idx_edge]] = min_val;
        }
    }
}

/**
 * Process filtering on messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: index of iteration of IMPA
 * @param[in] subtour2edge_in: messages from subtour elimination constraints to edge equality constraints before filtering
 * @param[out] subtour2edge_out: messages from subtour elimination constraints to edge equality constraints after filtering
 * @param[in] deltaS: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 *
 */

void SubtourEliminationConstraint::process_filtering(const int iter, const vector<vector<impalib_type>> &subtour2edge_in, vector<vector<impalib_type>> &subtour2edge_out,
                                                     const vector<vector<int>> &deltaS) {
    for (int i = 0; i < deltaS.size(); i++) {
        if ((doFilter) and (alpha_ != zero_value)) {
            vector<impalib_type> temp(subtour2edge_in[i]), temp_old(M_subtour2edge_old[i]), temp_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(temp.begin(), temp.end(), temp.begin(), [w_2](impalib_type &c) { return c * w_2; });
            transform(temp_old.begin(), temp_old.end(), temp_old.begin(), [w_1](impalib_type &c) { return c * w_1; });

            if (iter == 0) {
                copy(temp.begin(), temp.end(), subtour2edge_out[i].begin());
            } else {
                transform(temp.begin(), temp.end(), temp_old.begin(), back_inserter(temp_extrinsic), plus<impalib_type>());
                copy(temp_extrinsic.begin(), temp_extrinsic.end(), subtour2edge_out[i].begin());
            }

            copy(subtour2edge_out[i].begin(), subtour2edge_out[i].end(), M_subtour2edge_old[i].begin());
        }

        else {
            copy(subtour2edge_in[i].begin(), subtour2edge_in[i].end(), subtour2edge_out[i].begin());
        }
    }
}