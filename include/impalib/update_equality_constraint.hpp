// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"
/**
 * Represents a class for the equality constraint for the Knapsack-MWM problem
 */
class EqualityConstraintKcMwm {
   private:
    int numTeams_;        ///< number of teams
    int numProjects_;     ///< number of projects
    int numDepartments_;  ///< number of departments

   public:
    vector<impalib_type> team_messages_to_oric(const vector<vector<impalib_type>> &extrinsic,
                                               const vector<impalib_type> &reward_team);  ///< calculate messages from team equality constraint to ORIC

    void project_messages_to_oric(const vector<vector<impalib_type>> &proj2eq, vector<vector<impalib_type>> &eq2oric,
                                  const vector<vector<impalib_type>> &reward_proj);  ///< calculate messages from project equality constraint to ORIC

    EqualityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);  ///< constructor
};

/**
 * Construct Equality constraint object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numDepartments_: N_DEPARTMENTS
 */

EqualityConstraintKcMwm::EqualityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS) : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){};

/**
 * Calculate messages from team equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] extrinsic: messages from departments to teams
 * @param[out] rTeam2OricM: messages from team equality constraint to ORIC
 * @param[in] reward_team: rewards of teams
 */

vector<impalib_type> EqualityConstraintKcMwm::team_messages_to_oric(const vector<vector<impalib_type>> &extrinsic, const vector<impalib_type> &reward_team) {
    vector<impalib_type> team2oric(numTeams_, 0);

    for (int i = 0; i < extrinsic.size(); i++) {
        transform(extrinsic[i].begin(), extrinsic[i].end(), team2oric.begin(), team2oric.begin(), std::plus<impalib_type>());
    }
    transform(team2oric.begin(), team2oric.end(), reward_team.begin(), team2oric.begin(), std::plus<impalib_type>());
    return team2oric;
}

/**
 * Calculate messages from project equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] proj2eq: messages from projects inequality constraints to project equality constraint
 * @param[out] eq2oric: messages from project equality constraints to ORIC
 * @param[in] reward_proj: rewards of teams-projects combinations
 */

void EqualityConstraintKcMwm::project_messages_to_oric(const vector<vector<impalib_type>> &proj2eq, vector<vector<impalib_type>> &eq2oric, const vector<vector<impalib_type>> &reward_proj) {
    for (int i = 0; i < proj2eq.size(); i++) {
        transform(proj2eq[i].begin(), proj2eq[i].end(), reward_proj[i].begin(), eq2oric[i].begin(), std::plus<impalib_type>());
    }
}

/**
 * Represents a class for the equality constraint for the TSP
 */
class EqualityConstraintTsp {
   private:
    int numNodes_;          ///< number of nodes
    int numEdgeVariables_;  ///< number of edge connections
    bool filteringFlag_;    ///< filtering flag
    impalib_type alpha_;    ///< filtering parameter

   public:
    void messages_to_degree_relaxed(const vector<vector<int>> &edges, const vector<vector<impalib_type>> &cost_edges, const vector<vector<impalib_type>> &deg2eq,
                                    vector<vector<impalib_type>> &edge2deg);                          ///< calculate messages from edge to degree constraints for relaxed TSP
    vector<vector<impalib_type>> flip_matrix(const vector<vector<impalib_type>> &, const vector<vector<int>> &);  ///< flip matrix
    vector<vector<impalib_type>> messages_to_subtour(const vector<vector<int>> &deltaS, const vector<impalib_type> &edge_costs, const vector<vector<impalib_type>> &deg2edge, const vector<vector<impalib_type>> &subtour2edge,
                                                     const vector<vector<int>> &edges);  ///< calculate message from edge to subtour constraint
    void messages_to_degree_augmented(const vector<vector<impalib_type>> &deg2eq, const vector<vector<impalib_type>> &subtour2edge, const vector<vector<int>> &edges, const vector<vector<impalib_type>> &edge_costs,
                                      vector<vector<impalib_type>> &edge2deg);  ///< calculate messages from edge to degree constraints for augmented TSP

    EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                          const impalib_type ALPHA);  ///< constructor
};

/**
 * Construct equality constraint object for the TSP
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

EqualityConstraintTsp::EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES){};

/**
 * Calculate messages from edge equality constraints to degree constraints for the relaxed TSP
 *
 * @param[in] edges: list of connections for each edge equality constraint
 * @param[in] cost_edges: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[in] deg2eq: messages from degree constraints to equality constraints
 * @param[out] edge2deg: messages from edge equality constraints to degree constraints
 *
 */

void EqualityConstraintTsp::messages_to_degree_relaxed(const vector<vector<int>> &edges, const vector<vector<impalib_type>> &cost_edges, const vector<vector<impalib_type>> &deg2eq,
                                                       vector<vector<impalib_type>> &edge2deg) {
    auto reversed = flip_matrix(deg2eq, edges);

    for (int i = 0; i < numEdgeVariables_; i++) {
        transform(reversed[i].begin(), reversed[i].end(), cost_edges[i].begin(), edge2deg[i].begin(), plus<impalib_type>());
    }
}

/**
 * Calculate messages from edge equality constraints to subtour elimination constraints for the TSP
 *
 * @param[in] deltaS: list of edge indices forming the subtour elimination constraints
 * @param[in] edge_costs: costs of each edge variable
 * @param[in] deg2edge: messages from degree constraint to equality constraint
 * @param[in] subtour2edge: messages from subtour constraints to edge equality constraints
 * @param[in] edges: list of connections for each edge equality constraint
 * @return edge_ec_to_subtour_constraints_list: messages from edge equality constraint to subtour elimination constraints
 *
 */

vector<vector<impalib_type>> EqualityConstraintTsp::messages_to_subtour(const vector<vector<int>> &deltaS, const vector<impalib_type> &edge_costs, const vector<vector<impalib_type>> &deg2edge,
                                                                        const vector<vector<impalib_type>> &subtour2edge, const vector<vector<int>> &edges) {
    vector<vector<impalib_type>> edge2subtour_list;

    if (deltaS.size() == 1) {
        vector<impalib_type> edge2subtour(numEdgeVariables_, zero_value);

        vector<impalib_type> deg2eq_sum(numEdgeVariables_, zero_value);
        for (size_t i = 0; i < numEdgeVariables_; ++i) {
            deg2eq_sum[i] = deg2edge[i][edges[i][0]] + deg2edge[i][edges[i][1]];
        }

        for (size_t i = 0; i < deltaS[0].size(); i++) {
            edge2subtour[deltaS[0][i]] = deg2eq_sum[deltaS[0][i]] + edge_costs[deltaS[0][i]];
        }
        edge2subtour_list.push_back(edge2subtour);
    }

    else {
        vector<impalib_type> subtour2edge_sum(numEdgeVariables_, zero_value);
        for (const auto &row : subtour2edge) {
            transform(subtour2edge_sum.begin(), subtour2edge_sum.end(), row.begin(), subtour2edge_sum.begin(), std::plus<impalib_type>());
        }

        vector<impalib_type> deg2eq_sum(numEdgeVariables_, zero_value);
        for (size_t i = 0; i < numEdgeVariables_; ++i) {
            deg2eq_sum[i] = deg2edge[i][edges[i][0]] + deg2edge[i][edges[i][1]];
        }

        for (size_t i = 0; i < deltaS.size(); i++) {
            vector<impalib_type> edge2subtour(numEdgeVariables_, zero_value);

            for (size_t j = 0; j < deltaS[i].size(); j++) {
                edge2subtour[deltaS[i][j]] = subtour2edge_sum[deltaS[i][j]] + deg2eq_sum[deltaS[i][j]] + edge_costs[deltaS[i][j]] - subtour2edge[i][deltaS[i][j]];
            }
            edge2subtour_list.push_back(edge2subtour);
        }
    }
    return edge2subtour_list;
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the augmented TSP
 *
 * @param[in] deg2eq: messages from degree constraints to equality constraints
 * @param[in] subtour2edge: messages from subtour elimination constraints to edge equality constraints
 * @param[in] edges: list of connections for each edge equality constraint
 * @param[in] edge_costs: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[out] edge2deg: messages from edge equality constraints to degree constraints
 *
 */

void EqualityConstraintTsp::messages_to_degree_augmented(const vector<vector<impalib_type>> &deg2eq, const vector<vector<impalib_type>> &subtour2edge, const vector<vector<int>> &edges,
                                                         const vector<vector<impalib_type>> &edge_costs, vector<vector<impalib_type>> &edge2deg) {
    vector<impalib_type> subtour2edge_sum(numEdgeVariables_, zero_value);
    for (const auto &row : subtour2edge) {
        transform(subtour2edge_sum.begin(), subtour2edge_sum.end(), row.begin(), subtour2edge_sum.begin(), std::plus<impalib_type>());
    }

    for (size_t i = 0; i < edges.size(); i++) {
        // Update the message for the first node of the edge
        const auto first_node = edges[i][0];
        const auto second_node = edges[i][1];
        edge2deg[i][first_node] = subtour2edge_sum[i] + deg2eq[i][second_node] + edge_costs[i][first_node];

        // Update the message for the second node of the edge
        edge2deg[i][second_node] = subtour2edge_sum[i] + deg2eq[i][first_node] + edge_costs[i][second_node];
    }
}

/**
 * Flip a matrix to facilitate message updates for the TSP. Matrices will be in the same format during IMPA
 *
 * @param[in] M: matrix that requires flipping
 * @param[in] edges: list of connections for each edge equality constraint
 * @param[out] rFlippedMatrix: flipped matrix
 *
 */

vector<vector<impalib_type>> EqualityConstraintTsp::flip_matrix(const vector<vector<impalib_type>> &M, const vector<vector<int>> &edges) {
    auto flipped = M;
    for (size_t l = 0; l < edges.size(); ++l) {
        int row = edges[l][0];
        int col = edges[l][1];
        flipped[l][row] = M[l][col];
        flipped[l][col] = M[l][row];
    }
    return flipped;
}