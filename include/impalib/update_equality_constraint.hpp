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
class EqualityConstraintKcMwm
{
private:
    int numTeams_; ///< number of teams
    int numProjects_; ///< number of projects
    int numDepartments_; ///< number of departments

public:
    void team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                           vector<impalib_type> &); ///< calculate messages from team equality constraint to ORIC

    void project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &); ///< calculate messages from project equality constraint to ORIC

    EqualityConstraintKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS); ///< constructor
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

EqualityConstraintKcMwm::EqualityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                    };

/**
 * Calculate messages from team equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] rExtrinsicOutputDepartment: messages from departments to teams
 * @param[out] rTeam2OricM: messages from team equality constraint to ORIC
 * @param[in] rewardTeam: rewards of teams
 */

void EqualityConstraintKcMwm::team_eq_constraint_to_oric_update(
    vector<vector<impalib_type>> &rExtrinsicOutputDepartment, vector<impalib_type> &rTeam2OricM,
    vector<impalib_type> &rewardTeam)
{
    vector<impalib_type> intermediate_team_to_oric_m(numTeams_, 0);

    for (int i = 0; i < rExtrinsicOutputDepartment.size(); i++)
    {
        transform(rExtrinsicOutputDepartment[i].begin(),
                  rExtrinsicOutputDepartment[i].end(), intermediate_team_to_oric_m.begin(),
                  intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
    }
    transform(intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.end(), rewardTeam.begin(),
              rTeam2OricM.begin(), std::plus<impalib_type>());
}

/**
 * Calculate messages from project equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] rProject2EqConstraintM: messages from projects inequality constraints to project equality constraint
 * @param[out] rEqConstraint2OricM: messages from project equality constraints to ORIC
 * @param[in] rewardProject: rewards of teams-projects combinations
 */

void EqualityConstraintKcMwm::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rProject2EqConstraintM,
                                                                   vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                                   vector<vector<impalib_type>> &rewardProject)
{
    for (int i = 0; i < rProject2EqConstraintM.size(); i++)
    {
        transform(rProject2EqConstraintM[i].begin(), rProject2EqConstraintM[i].end(),
                  rewardProject[i].begin(), rEqConstraint2OricM[i].begin(),
                  std::plus<impalib_type>());
    }
}

/**
 * Represents a class for the equality constraint for the TSP
 */
class EqualityConstraintTsp
{
private:
    int          numNodes_; ///< number of nodes
    int          numEdgeVariables_; ///< number of edge connections
    bool         filteringFlag_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter

public:
    void edge_ec_to_degree_constraint_relaxed_graph_update(vector<vector<int>> &, vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &); ///< calculate messages from edge to degree constraints for relaxed TSP
    vector<vector<impalib_type>> flip_matrix(vector<vector<impalib_type>> &, vector<vector<int>> &); ///< flip matrix
    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_update(vector<vector<int>> &, vector<impalib_type> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<int>> &); ///< calculate message from edge to subtour constraint
    void                         edge_ec_to_degree_constraint_augmented_graph_update(vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &, vector<vector<int>> &,
                                                                                     vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &); ///< calculate messages from edge to degree constraints for augmented TSP

    EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                          const impalib_type ALPHA); ///< constructor
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

EqualityConstraintTsp::EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES),
      numEdgeVariables_(NUM_EDGE_VARIABLES){
      };

/**
 * Calculate messages from edge equality constraints to degree constraints for the relaxed TSP
 *
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[in] rEdgeDegreeConstraintCost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to equality constraints
 * @param[out] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * 
 */

void EqualityConstraintTsp::edge_ec_to_degree_constraint_relaxed_graph_update(
    vector<vector<int>> &rEdgeConnections, vector<vector<impalib_type>> &rEdgeDegreeConstraintCost,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM)
{

    auto reversed = flip_matrix(rDegreeConstraint2EqConstraintM, rEdgeConnections);

    for (int i = 0; i < numEdgeVariables_; i++)
    {
        transform(reversed[i].begin(),
                  reversed[i].end(),
                  rEdgeDegreeConstraintCost[i].begin(),
                  rEdgeEc2DegreeConstraintM[i].begin(), plus<impalib_type>());
    }
}


/**
 * Calculate messages from edge equality constraints to subtour elimination constraints for the TSP
 *
 * @param[in] rDeltaSIndicesList: list of edge indices forming the subtour elimination constraints
 * @param[in] rCostEdgeVaribale: costs of each edge variable
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraint to equality constraint
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour constraints to edge equality constraints
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @return edge_ec_to_subtour_constraints_list: messages from edge equality constraint to subtour elimination constraints
 * 
 */

vector<vector<impalib_type>> EqualityConstraintTsp::edge_ec_to_subtour_constraints_update(
    vector<vector<int>> &rDeltaSIndicesList, vector<impalib_type> &rCostEdgeVaribale,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, vector<vector<int>> &rEdgeConnections)
{

    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_list;

    if (rDeltaSIndicesList.size() == 1)
    {
        vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t i = 0; i < numEdgeVariables_; ++i)
        {
            combined_degree_constraint_to_eq_constraint_m[i] =
                rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][0]]
                + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][1]];
        }

        for (size_t i = 0; i < rDeltaSIndicesList[0].size(); i++)
        {
            edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[0][i]] =
                combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[0][i]]
                + rCostEdgeVaribale[rDeltaSIndicesList[0][i]];
        }
        edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
    }
    
    else
    {
        vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
        for (const auto &row : rSubtourConstraints2EdgeEcM)
        {
            transform(combined_subtour_constraints_to_edge_ec_m.begin(),
                      combined_subtour_constraints_to_edge_ec_m.end(), row.begin(),
                      combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
        }

        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t i = 0; i < numEdgeVariables_; ++i)
        {
            combined_degree_constraint_to_eq_constraint_m[i] =
                rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][0]]
                + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][1]];
        }

        for (size_t i = 0; i < rDeltaSIndicesList.size(); i++)
        {
            vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

            for (size_t j = 0; j < rDeltaSIndicesList[i].size(); j++)
            {
                edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[i][j]] =
                    combined_subtour_constraints_to_edge_ec_m[rDeltaSIndicesList[i][j]]
                    + combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[i][j]]
                    + rCostEdgeVaribale[rDeltaSIndicesList[i][j]]
                    - rSubtourConstraints2EdgeEcM[i]
                                                 [rDeltaSIndicesList[i][j]];
            }
            edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
        }
    }
    return edge_ec_to_subtour_constraints_list;
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the augmented TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to equality constraints
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[in] rEdgeDegreeConstraintCost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[out] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * 
 */

void EqualityConstraintTsp::edge_ec_to_degree_constraint_augmented_graph_update(
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, vector<vector<int>> &rEdgeConnections,
    vector<vector<impalib_type>> &rEdgeDegreeConstraintCost, vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM)
{

    vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
    for (const auto &row : rSubtourConstraints2EdgeEcM)
    {
        transform(combined_subtour_constraints_to_edge_ec_m.begin(), combined_subtour_constraints_to_edge_ec_m.end(),
                  row.begin(), combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
    }

    for (size_t i = 0; i < rEdgeConnections.size(); i++)
    {

        // Update the message for the first node of the edge
        const auto first_node = rEdgeConnections[i][0];
        const auto second_node = rEdgeConnections[i][1];
        rEdgeEc2DegreeConstraintM[i][first_node] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][second_node]
            + rEdgeDegreeConstraintCost[i][first_node];

        // Update the message for the second node of the edge
        rEdgeEc2DegreeConstraintM[i][second_node] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][first_node]
            + rEdgeDegreeConstraintCost[i][second_node];
    }
}

/**
 * Flip a matrix to facilitate message updates for the TSP. Matrices will be in the same format during IMPA
 *
 * @param[in] rMatrix: matrix that requires flipping
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[out] rFlippedMatrix: flipped matrix
 * 
 */

vector<vector<impalib_type>> EqualityConstraintTsp::flip_matrix(vector<vector<impalib_type>> &rMatrix, vector<vector<int>> &rEdgeConnections)
{
    auto flipped = rMatrix;
    for (size_t l = 0; l < rEdgeConnections.size(); ++l)
    {
        int row                = rEdgeConnections[l][0];
        int col                = rEdgeConnections[l][1];
        flipped[l][row] = rMatrix[l][col];
        flipped[l][col] = rMatrix[l][row];
    }
    return flipped;
}