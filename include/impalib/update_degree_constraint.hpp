// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include <cassert>

#include "impalib/impalib.hpp"

/**
 * Represents a class for the degree constraint for the TSP
 */
class DegreeConstraint {
   private:
    int numNodes_;                              ///< number of nodes of TSP
    int numEdgeVariables_;                      ///< number of edge connections
    bool filteringFlag_;                        ///< filtering flag
    impalib_type alpha_;                        ///< filtering parameter
    vector<vector<impalib_type>> M_deg2eq_old;  ///< messages from degree constraints to equality constraints before filtering
    impalib_type M_forward_init = value_inf;    ///< initial forward message of forward-backward algorithm
    impalib_type M_backward_init = value_inf;   ///< initial backward message of forward-backward algorithm

   public:
    void messages_to_edge_ec(const vector<vector<impalib_type>> &edge2degree, const vector<vector<int>> &edges,
                             vector<vector<impalib_type>> &degree2eq);                                        ///< calculate messages from degree constraint to edge equality constraint
    vector<vector<impalib_type>> process_filtering(int iter, const vector<vector<impalib_type>> &deg2eq_in);  ///< process filtering on messages from degree constraint to edge equality constraint

    DegreeConstraint(int NUM_NODES, int NUM_EDGE_VARIABLES, bool FILTERING_FLAG,
                     impalib_type ALPHA);  ///< constructor
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

DegreeConstraint::DegreeConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES), M_deg2eq_old(numEdgeVariables_, vector<impalib_type>(numNodes_, 0)){};

/**
 * Calculate messages from degree constraints to edge equality constraints for the TSP
 *
 * @param[in] edge2degree: messages from edge equality constraints to degree constraints
 * @param[in] edges: list of connections for each edge equality constraint
 * @param[out] degree2eq: messages from degree constraints to edge equality constraints before filtering
 *
 */

void DegreeConstraint::messages_to_edge_ec(const vector<vector<impalib_type>> &edge2degree, const vector<vector<int>> &edges, vector<vector<impalib_type>> &degree2eq) {
    vector<impalib_type> forward(numEdgeVariables_ + 1, zero_value);
    vector<impalib_type> backward(numEdgeVariables_ + 1, zero_value);

    for (int i = 0; i < numNodes_; i++) {
        vector<int> first, second;

        // Populate first and second with edge indices
        for (int j = 0; j < edges.size(); j++) {
            if (i == edges[j][0]) {
                first.push_back(j);
            } else if (i == edges[j][1]) {
                second.push_back(j);
            }
        }

        // Calculate forward messages for first
        forward[first[0]] = M_forward_init;

        for (int s = 1; s < first.size(); s++) {
            forward[first[s]] = min(forward[first[s - 1]], edge2degree[first[s - 1]][i]);
        }

        // Calculate backward messages for first
        backward[first[first.size() - 1] + 1] = M_backward_init;

        for (size_t s = first.size() - 1; s >= 1; s--) {
            backward[first[s - 1] + 1] = min(backward[first[s] + 1], edge2degree[first[s]][i]);
        }

        for (int j = 0; j < first.size(); j++) {
            impalib_type minimumValue = zero_value;
            minimumValue = min(forward[first[j]], backward[first[j] + 1]);
            degree2eq[first[j]][i] = -minimumValue;
        }

        fill(forward.begin(), forward.end(), zero_value);
        fill(backward.begin(), backward.end(), zero_value);

        // Calculate forward messages
        forward[second[0]] = M_forward_init;

        for (int s = 1; s < second.size(); s++) {
            forward[second[s]] = min(forward[second[s - 1]], edge2degree[second[s - 1]][i]);
        }

        // Calculate backward messages
        backward[second[second.size() - 1] + 1] = M_backward_init;

        for (size_t s = second.size() - 1; s >= 1; s--) {
            backward[second[s - 1] + 1] = min(backward[second[s] + 1], edge2degree[second[s]][i]);
        }

        for (int j = 0; j < second.size(); j++) {
            impalib_type minimumValue = zero_value;
            minimumValue = min(forward[second[j]], backward[second[j] + 1]);
            degree2eq[second[j]][i] = -minimumValue;
        }
    }
}

/**
 * Perform filtering on messages from degree constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: iteration index of IMPA
 * @param[in] deg2eq_in: messages from degree constraints to edge equality constraints before filtering
 * @param[in] deg2eq_out: messages from degree constraints to edge equality constraints after filtering
 *
 */

vector<vector<impalib_type>> DegreeConstraint::process_filtering(const int iter, const vector<vector<impalib_type>> &deg2eq_in) {
    assert(deg2eq_in.size() == numEdgeVariables_);
    if (!filteringFlag_) {
        return deg2eq_in;
    }

    // Calculate weighted values for current and old messages
    auto deg2eq_out = deg2eq_in;
    for (int i = 0; i < numEdgeVariables_; ++i) {
        for (int j = 0; j < deg2eq_out.size(); ++j) {
            deg2eq_out[i][j] = (1 - alpha_) * deg2eq_out[i][j] + (alpha_)*M_deg2eq_old[i][j];
        }
    }
    M_deg2eq_old = deg2eq_out;
    return deg2eq_out;
}
