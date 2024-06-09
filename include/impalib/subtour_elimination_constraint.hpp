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
    int          nNodes_; ///< number of nodes
    int          nEdgeVars_; ///< number of edges
    bool         doFilter_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter
    impalib_type forward0_; ///< initial value of forward message using forward-backward algorithm
    impalib_type backward0_; ///< initial value of backward message using forward-backward algorithm

public:
    vector<vector<impalib_type>> subtour2EqOldM_; ///< messages from subtour constraints to edge equality constraint before filtering
    void subtour_constraints_to_edge_ec_update(const vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                               vector<vector<impalib_type>> &) const; ///< calculate messages from subtour to edge equality constraints
    void process_filtering(int, const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, const vector<vector<int>> &); ///< perform filtering on messages from subtour to edge equality constraints

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
inline SubtourEliminationConstraint::SubtourEliminationConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                                           const bool FILTERING_FLAG, const impalib_type ALPHA)
    : doFilter_(FILTERING_FLAG), alpha_(ALPHA), nNodes_(NUM_NODES), nEdgeVars_(NUM_EDGE_VARIABLES),
      forward0_(value_inf), backward0_(value_inf){
                                           };

/**
 * Calculate messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] eq2SubtourM: messages from edge equality constraints to subtour elimination constraints
 * @param[in] deltaS: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 * @param[out] subtour2EqM: messages from subtour elimination constraints to edge equality constraints
 * 
 */

inline void SubtourEliminationConstraint::subtour_constraints_to_edge_ec_update(
    const vector<vector<impalib_type>> &eq2SubtourM, const vector<vector<int>> &deltaS,
    vector<vector<impalib_type>> &subtour2EqM) const
{
    vector<impalib_type> forward(nEdgeVars_ + 1, zero_value);
    vector<impalib_type> backward(nEdgeVars_ + 1, zero_value);

    for (size_t subtour = 0; subtour < deltaS.size();
         subtour++)
    {
        // Initialize and calculate forward messages
        forward[deltaS[subtour][0]] = forward0_;
        for (int stage = 1; stage < deltaS[subtour].size(); stage++)
        {
            forward[deltaS[subtour][stage]] =
                min(forward[deltaS[subtour][stage - 1]],
                    eq2SubtourM[subtour]
                                               [deltaS[subtour][stage - 1]]);
        }

        // Initialize and calculate backward messages
        backward[deltaS[subtour]
                                                  [deltaS[subtour].size() - 1]
                                + 1] = backward0_;

        for (size_t stage = deltaS[subtour].size() - 1; stage >= 1; stage--)
        {
            backward[deltaS[subtour][stage - 1] + 1] =
                min(backward[deltaS[subtour][stage] + 1],
                    eq2SubtourM[subtour]
                                               [deltaS[subtour][stage]]);
        }

        // Update subtour constraints to edge EC messages
        for (int edge = 0; edge < deltaS[subtour].size();
             edge++)
        {
            impalib_type minimumValue =
                min(forward[deltaS[subtour][edge]],
                    backward[deltaS[subtour][edge] + 1]);
            minimumValue = min(-minimumValue, zero_value);
            subtour2EqM[subtour]
                                       [deltaS[subtour][edge]] =
                                           minimumValue;
        }
    }
}

/**
 * Process filtering on messages from subtour elimination constraints to edge equality constraints for the TSP
 *
 * @param[in] iter: index of iteration of IMPA
 * @param[in] subtour2EqPreM: messages from subtour elimination constraints to edge equality constraints before filtering
 * @param[out] subtour2EqM: messages from subtour elimination constraints to edge equality constraints after filtering
 * @param[in] deltaS: list of sub-list of indices. Each sub-list is asscoiated with a subtour constraint, and contains indices of edges related to this constraint
 * 
 */

inline void SubtourEliminationConstraint::process_filtering(const int                           iter,
                                                     const vector<vector<impalib_type>> &subtour2EqPreM,
                                                     vector<vector<impalib_type>> &subtour2EqM,
                                                     const vector<vector<int>>          &deltaS)
{
    if (!doFilter_) {
        subtour2EqM = subtour2EqPreM;
        return;
    }

    // Calculate weighted values for current and old messages
    for (int subtour = 0; subtour < deltaS.size(); subtour++) {
        if (iter == 0) {
            subtour2EqM[subtour] = subtour2EqPreM[subtour];
        } else {
            vector<impalib_type> weighted(subtour2EqPreM[subtour].size());
            for (int i=0; i<weighted.size(); ++i) {
                weighted[i] = alpha_*subtour2EqOldM_[subtour][i] + (1-alpha_)*subtour2EqPreM[subtour][i];
            }
            subtour2EqM[subtour] = weighted;
        }
        subtour2EqOldM_[subtour] = subtour2EqM[subtour];
    }
}