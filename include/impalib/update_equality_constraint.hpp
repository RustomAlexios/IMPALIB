// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class EqualityConstraint
{
private:
    int nTeams_; ///< number of teams
    int nProj_; ///< number of projects
    int nDept_; ///< number of departments
    int          nNodes_; ///< number of nodes
    int          nEdgeVars_; ///< number of edge connections
    bool         doFilter_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter
    int nVars_;
    int nConstraints_;
    int k_;
public:

    EqualityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : nProj_(N_PROJECTS), nTeams_(N_TEAMS), nDept_(N_DEPARTMENTS),
        doFilter_(false), alpha_(zero_value), nNodes_(0), nEdgeVars_(0),
        nVars_(0), nConstraints_(0), k_(0){};

    EqualityConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : nProj_(0), nTeams_(0), nDept_(0),
      doFilter_(FILTERING_FLAG), alpha_(ALPHA), nNodes_(NUM_NODES),
      nEdgeVars_(NUM_EDGE_VARIABLES), nVars_(0), nConstraints_(0), k_(0){};

    EqualityConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : nProj_(0), nTeams_(0), nDept_(0),
      doFilter_(FILTERING_FLAG), alpha_(ALPHA), nNodes_(0),
      nEdgeVars_(0), nVars_(NUM_VARIABLES), nConstraints_(NUM_CONSTRAINTS), k_(K_VARIABLE){};

    void team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                           vector<impalib_type> &) const; ///< calculate messages from team equality constraint to ORIC

    static void project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &); ///< calculate messages from project equality constraint to ORIC

    void edge_ec_to_degree_constraint_relaxed_graph_update(const vector<vector<int>> &, vector<vector<impalib_type>> &,
                                                           const vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &) const; ///< calculate messages from edge to degree constraints for relaxed TSP

    static void flip_matrix(const vector<vector<impalib_type>> &, const vector<vector<int>> &, vector<vector<impalib_type>> &); ///< flip matrix
    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_update(const vector<vector<int>> &, const vector<impalib_type> &,
                                                                       const vector<vector<impalib_type>> &,
                                                                       const vector<vector<impalib_type>> &,
                                                                       const vector<vector<int>> &) const; ///< calculate message from edge to subtour constraint
    void                         edge_ec_to_degree_constraint_augmented_graph_update(const vector<vector<impalib_type>> &,
                                                                                     const vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                                                                     const vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &) const; ///< calculate messages from edge to degree constraints for augmented TSP

    void variable_ec_to_ksat_constraint_update(const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<int> &, const vector<impalib_type> &, const vector<vector<int>> &) const;
};

/**
 * Calculate messages from team equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] extrinsicOut: messages from departments to teams
 * @param[out] team2OricM: messages from team equality constraint to ORIC
 * @param[in] rewards: rewards of teams
 */

inline void EqualityConstraint::team_eq_constraint_to_oric_update(
    vector<vector<impalib_type>> &extrinsicOut, vector<impalib_type> &team2OricM,
    vector<impalib_type> &rewards) const
{
    vector<impalib_type> intermediate_team_to_oric_m(nTeams_, 0);

    for (int department = 0; department < extrinsicOut.size(); department++)
    {
        transform(extrinsicOut[department].begin(),
                  extrinsicOut[department].end(), intermediate_team_to_oric_m.begin(),
                  intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
    }
    transform(intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.end(), rewards.begin(),
              team2OricM.begin(), std::plus<impalib_type>());
}

/**
 * Calculate messages from project equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] project2EqM: messages from projects inequality constraints to project equality constraint
 * @param[out] eq2OricM: messages from project equality constraints to ORIC
 * @param[in] rewards: rewards of teams-projects combinations
 */

inline void EqualityConstraint::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &project2EqM,
                                                                   vector<vector<impalib_type>> &eq2OricM,
                                                                   vector<vector<impalib_type>> &rewards)
{
    for (int project = 0; project < project2EqM.size(); project++)
    {
        transform(project2EqM[project].begin(), project2EqM[project].end(),
                  rewards[project].begin(), eq2OricM[project].begin(),
                  std::plus<impalib_type>());
    }
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the relaxed TSP
 *
 * @param[in] connections: list of connections for each edge equality constraint
 * @param[in] cost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[in] deg2EqM: messages from degree constraints to equality constraints
 * @param[out] eq2DegM: messages from edge equality constraints to degree constraints
 *
 */

inline void EqualityConstraint::edge_ec_to_degree_constraint_relaxed_graph_update(
    const vector<vector<int>> &connections, vector<vector<impalib_type>> &cost,
    const vector<vector<impalib_type>> &deg2EqM,
    vector<vector<impalib_type>> &eq2DegM) const
{

    vector<vector<impalib_type>> deg2EqM_flipped = deg2EqM;
    flip_matrix(deg2EqM, connections, deg2EqM_flipped);

    for (int edge = 0; edge < nEdgeVars_; edge++)
    {
        transform(deg2EqM_flipped[edge].begin(),
                  deg2EqM_flipped[edge].end(),
                  cost[edge].begin(),
                  eq2DegM[edge].begin(), plus<impalib_type>());
    }
}


/**
 * Calculate messages from edge equality constraints to subtour elimination constraints for the TSP
 *
 * @param[in] deltaS: list of edge indices forming the subtour elimination constraints
 * @param[in] costs: costs of each edge variable
 * @param[in] deg2EqM: messages from degree constraint to equality constraint
 * @param[in] subtour2EqM: messages from subtour constraints to edge equality constraints
 * @param[in] connections: list of connections for each edge equality constraint
 * @return edge_ec_to_subtour_constraints_list: messages from edge equality constraint to subtour elimination constraints
 *
 */

inline vector<vector<impalib_type>> EqualityConstraint::edge_ec_to_subtour_constraints_update(
    const vector<vector<int>> &deltaS, const vector<impalib_type> &costs,
    const vector<vector<impalib_type>> &deg2EqM,
    const vector<vector<impalib_type>> &subtour2EqM, const vector<vector<int>> &connections) const
{

    vector<vector<impalib_type>> eq2SubtourM;

    if (deltaS.size() == 1)
    {
        vector<impalib_type> edge_ec_to_subtour_constraints_m(nEdgeVars_, zero_value);

        vector<impalib_type> deg2EqM_sum(nEdgeVars_, zero_value);
        for (size_t edge = 0; edge < nEdgeVars_; ++edge)
        {
            deg2EqM_sum[edge] =
                deg2EqM[edge][connections[edge][0]]
                + deg2EqM[edge][connections[edge][1]];
        }

        for (size_t i = 0; i < deltaS[0].size(); i++)
        {
            edge_ec_to_subtour_constraints_m[deltaS[0][i]] =
                deg2EqM_sum[deltaS[0][i]]
                + costs[deltaS[0][i]];
        }
        eq2SubtourM.push_back(edge_ec_to_subtour_constraints_m);
    }

    else
    {
        vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(nEdgeVars_, zero_value);
        for (const auto &row : subtour2EqM)
        {
            transform(combined_subtour_constraints_to_edge_ec_m.begin(),
                      combined_subtour_constraints_to_edge_ec_m.end(), row.begin(),
                      combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
        }

        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(nEdgeVars_, zero_value);
        for (size_t edge = 0; edge < nEdgeVars_; ++edge)
        {
            combined_degree_constraint_to_eq_constraint_m[edge] =
                deg2EqM[edge][connections[edge][0]]
                + deg2EqM[edge][connections[edge][1]];
        }

        for (size_t subtour = 0; subtour < deltaS.size();
             subtour++)
        {
            vector<impalib_type> edge_ec_to_subtour_constraints_m(nEdgeVars_, zero_value);

            for (size_t i = 0; i < deltaS[subtour].size(); i++)
            {
                edge_ec_to_subtour_constraints_m[deltaS[subtour][i]] =
                    combined_subtour_constraints_to_edge_ec_m[deltaS[subtour][i]]
                    + combined_degree_constraint_to_eq_constraint_m[deltaS[subtour][i]]
                    + costs[deltaS[subtour][i]]
                    - subtour2EqM[subtour]
                                                 [deltaS[subtour][i]];
            }
            eq2SubtourM.push_back(edge_ec_to_subtour_constraints_m);
        }
    }
    return eq2SubtourM;
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the augmented TSP
 *
 * @param[in] deg2EqM: messages from degree constraints to equality constraints
 * @param[in] subtour2EqM: messages from subtour elimination constraints to edge equality constraints
 * @param[in] connections: list of connections for each edge equality constraint
 * @param[in] cost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[out] eq2DegreeM: messages from edge equality constraints to degree constraints
 *
 */

inline void EqualityConstraint::edge_ec_to_degree_constraint_augmented_graph_update(
    const vector<vector<impalib_type>> &deg2EqM,
    const vector<vector<impalib_type>> &subtour2EqM, const vector<vector<int>> &connections,
    const vector<vector<impalib_type>> &cost, vector<vector<impalib_type>> &eq2DegreeM) const
{

    vector<impalib_type> subtour2EqM_sum(nEdgeVars_, zero_value);
    for (const auto &row : subtour2EqM)
    {
        transform(subtour2EqM_sum.begin(), subtour2EqM_sum.end(),
                  row.begin(), subtour2EqM_sum.begin(), std::plus<impalib_type>());
    }

    for (size_t i = 0; i < connections.size(); i++)
    {

        // Update the message for the first node of the edge
        eq2DegreeM[i][connections[i][0]] =
            subtour2EqM_sum[i] + deg2EqM[i][connections[i][1]]
            + cost[i][connections[i][0]];

        // Update the message for the second node of the edge
        eq2DegreeM[i][connections[i][1]] =
            subtour2EqM_sum[i] + deg2EqM[i][connections[i][0]]
            + cost[i][connections[i][1]];
    }
}

/**
 * Flip a matrix to facilitate message updates for the TSP. Matrices will be in the same format during IMPA
 *
 * @param[in] mat: matrix that requires flipping
 * @param[in] connections: list of connections for each edge equality constraint
 * @param[out] out: flipped matrix
 *
 */

inline void EqualityConstraint::flip_matrix(const vector<vector<impalib_type>> &mat, const vector<vector<int>> &connections,
                                        vector<vector<impalib_type>> &out)
{
    // Iterate over each edge connection
    for (size_t l = 0; l < connections.size(); ++l)
    {
        // Row index
        int row                = connections[l][0];
        // Column index
        int col                = connections[l][1];
        out[l][row] = mat[l][col];
        out[l][col] = mat[l][row];
    }
}

/**
 * Calculate messages from variable equality constraints to k-sat constraints for the K-SAT problem
 *
 * @param[in] ksat2EqM: messages from k-sat constraints to variable equality constraints
 * @param[out] var2KsatM: messages variable equality constraints to from k-sat constraints
 * @param[in] used: used variables in creating the constraints
 * @param[in] connections: constraints connections for each variable
 *
 */

inline void EqualityConstraint::variable_ec_to_ksat_constraint_update(const vector<vector<impalib_type>> &ksat2EqM, vector<vector<impalib_type>> &var2KsatM, vector<int> &used, const vector<impalib_type> &costs, const vector<vector<int>> &connections) const
{

    for(auto& row : var2KsatM) {
        row.assign(row.size(), zero_value);
    }

    vector<impalib_type> used_incoming_metrics_cost(nVars_, zero_value);

    for_each(used.begin(), used.end(), [&](int n) {
        used_incoming_metrics_cost[n] = costs[n];
    });

    vector<impalib_type> sum_messages(nVars_, zero_value);

    for (int i = 0; i < nVars_; ++i) {
        for (int j = 0; j < connections[i].size(); ++j) {
            sum_messages[i] += ksat2EqM[connections[i][j]][i];
        }
        sum_messages[i] +=used_incoming_metrics_cost[i];
    }

    for (int index_variable = 0; index_variable < used.size(); ++index_variable) {
        int variable = used[index_variable];
        for (int i = 0; i < connections[variable].size(); ++i) {
            int constraint = connections[variable][i];
            // This if statement check was added to account for the fact that a constraint can have the same variable more than once,
            // like in the benchmarks datasets. However, in practical cases, a variable cannot appear more than once in a constraint
            // and thus this if statement check can be dropped
            if (abs(var2KsatM[constraint][variable])<abs(sum_messages[variable] - ksat2EqM[constraint][variable])){
                var2KsatM[constraint][variable] = sum_messages[variable] - ksat2EqM[constraint][variable];
            }
        }
    }
}