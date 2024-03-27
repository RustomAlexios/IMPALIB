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
    int numTeams_; ///< number of teams
    int numProjects_; ///< number of projects
    int numDepartments_; ///< number of departments
    int          numNodes_; ///< number of nodes
    int          numEdgeVariables_; ///< number of edge connections
    bool         filteringFlag_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter

public:

    EqualityConstraint(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS), 
        filteringFlag_(false), alpha_(zero_value), numNodes_(0), numEdgeVariables_(0){};

    EqualityConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : numProjects_(0), numTeams_(0), numDepartments_(0),
      filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES),
      numEdgeVariables_(NUM_EDGE_VARIABLES){};

    void team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                           vector<impalib_type> &); ///< calculate messages from team equality constraint to ORIC

    void project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &); ///< calculate messages from project equality constraint to ORIC

    void edge_ec_to_degree_constraint_relaxed_graph_update(vector<vector<int>> &, vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &); ///< calculate messages from edge to degree constraints for relaxed TSP
    
    void flip_matrix(vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<impalib_type>> &); ///< flip matrix
    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_update(vector<vector<int>> &, vector<impalib_type> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<int>> &); ///< calculate message from edge to subtour constraint
    void                         edge_ec_to_degree_constraint_augmented_graph_update(vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &, vector<vector<int>> &,
                                                                                     vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &); ///< calculate messages from edge to degree constraints for augmented TSP                                             
};

/**
 * Calculate messages from team equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] rExtrinsicOutputDepartment: messages from departments to teams
 * @param[out] rTeam2OricM: messages from team equality constraint to ORIC
 * @param[in] rewardTeam: rewards of teams
 */

void EqualityConstraint::team_eq_constraint_to_oric_update(
    vector<vector<impalib_type>> &rExtrinsicOutputDepartment, vector<impalib_type> &rTeam2OricM,
    vector<impalib_type> &rewardTeam)
{
    // Vector to store intermediate values for team to ORIC messages
    vector<impalib_type> intermediate_team_to_oric_m(numTeams_, 0);

    // Iterate over each department
    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++)
    {
        // Sum rExtrinsicOutputDepartment
        transform(rExtrinsicOutputDepartment[department_index].begin(),
                  rExtrinsicOutputDepartment[department_index].end(), intermediate_team_to_oric_m.begin(),
                  intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
    }
    // Calculate rTeam2OricM
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

void EqualityConstraint::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rProject2EqConstraintM,
                                                                   vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                                   vector<vector<impalib_type>> &rewardProject)
{
    // Iterate over each project
    for (int project_index = 0; project_index < rProject2EqConstraintM.size(); project_index++)
    {
        // Transform and sum the messages from project to equality constraints and rewards for each project to get rEqConstraint2OricM
        transform(rProject2EqConstraintM[project_index].begin(), rProject2EqConstraintM[project_index].end(),
                  rewardProject[project_index].begin(), rEqConstraint2OricM[project_index].begin(),
                  std::plus<impalib_type>());
    }
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the relaxed TSP
 *
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[in] rEdgeDegreeConstraintCost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to equality constraints
 * @param[out] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * 
 */

void EqualityConstraint::edge_ec_to_degree_constraint_relaxed_graph_update(
    vector<vector<int>> &rEdgeConnections, vector<vector<impalib_type>> &rEdgeDegreeConstraintCost,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM)
{

    // Create a copy of the messages from degree constraints to equality constraints
    vector<vector<impalib_type>> flipped_degree_constraint_to_eq_constraint_m = rDegreeConstraint2EqConstraintM;

    // Flip the matrix representing the messages from degree constraints to equality constraints
    flip_matrix(rDegreeConstraint2EqConstraintM, rEdgeConnections, flipped_degree_constraint_to_eq_constraint_m);

    // Update the messages from edge equality constraints to degree constraints
    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++)
    {
        transform(flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].begin(),
                  flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].end(),
                  rEdgeDegreeConstraintCost[edge_variable_index].begin(),
                  rEdgeEc2DegreeConstraintM[edge_variable_index].begin(), plus<impalib_type>());
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

vector<vector<impalib_type>> EqualityConstraint::edge_ec_to_subtour_constraints_update(
    vector<vector<int>> &rDeltaSIndicesList, vector<impalib_type> &rCostEdgeVaribale,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, vector<vector<int>> &rEdgeConnections)
{

    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_list;

    // If there's only one subtour constraint
    if (rDeltaSIndicesList.size() == 1)
    {
        // Initialize a vector to store the messages from edge equality constraints to subtour constraints
        vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

        // Initialize a vector to store the combined messages from degree constraints to equality constraints
        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            // Combine the messages from degree constraints to equality constraints for each edge variable
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
        }

        // Update the messages from edge equality constraints to subtour constraints
        for (size_t i = 0; i < rDeltaSIndicesList[0].size(); i++)
        {
            edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[0][i]] =
                combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[0][i]]
                + rCostEdgeVaribale[rDeltaSIndicesList[0][i]];
        }
        // Add the updated messages to the list
        edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
    }
    
    else
    {
        // Initialize a vector to store the combined messages from subtour constraints to edge equality constraints
        vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
        for (const auto &row : rSubtourConstraints2EdgeEcM)
        {
            // Sum the messages from subtour constraints to edge equality constraints for each edge variable
            transform(combined_subtour_constraints_to_edge_ec_m.begin(),
                      combined_subtour_constraints_to_edge_ec_m.end(), row.begin(),
                      combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
        }

        // Initialize a vector to store the combined messages from degree constraints to equality constraints
        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            // Combine the messages from degree constraints to equality constraints for each edge variable
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
        }

        // Iterate over each subtour constraint
        for (size_t index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
             index_subtour_constraint++)
        {
            // Initialize a vector to store the messages from edge equality constraints to subtour constraints
            vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

            // Update the messages from edge equality constraints to subtour constraints
            for (size_t i = 0; i < rDeltaSIndicesList[index_subtour_constraint].size(); i++)
            {
                edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[index_subtour_constraint][i]] =
                    combined_subtour_constraints_to_edge_ec_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + rCostEdgeVaribale[rDeltaSIndicesList[index_subtour_constraint][i]]
                    - rSubtourConstraints2EdgeEcM[index_subtour_constraint]
                                                 [rDeltaSIndicesList[index_subtour_constraint][i]];
            }
            // Add the updated messages to the list
            edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
        }
    }
    // Return the list of udpated messages from edge equality constraints to subtour constraints
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

void EqualityConstraint::edge_ec_to_degree_constraint_augmented_graph_update(
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, vector<vector<int>> &rEdgeConnections,
    vector<vector<impalib_type>> &rEdgeDegreeConstraintCost, vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM)
{

    // Initialize a vector to store the combined messages from subtour constraints to edge equality constraints
    vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
    for (const auto &row : rSubtourConstraints2EdgeEcM)
    {
        // Sum the messages from subtour constraints to edge equality constraints for each edge variable
        transform(combined_subtour_constraints_to_edge_ec_m.begin(), combined_subtour_constraints_to_edge_ec_m.end(),
                  row.begin(), combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
    }

    // Update the messages from edge equality constraints to degree constraints
    for (size_t i = 0; i < rEdgeConnections.size(); i++)
    {

        // Update the message for the first node of the edge
        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][0]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][1]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][0]];

        // Update the message for the second node of the edge
        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][1]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][0]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][1]];
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

void EqualityConstraint::flip_matrix(vector<vector<impalib_type>> &rMatrix, vector<vector<int>> &rEdgeConnections,
                                        vector<vector<impalib_type>> &rFlippedMatrix)
{
    // Iterate over each edge connection
    for (size_t l = 0; l < rEdgeConnections.size(); ++l)
    {
        // Row index
        int row                = rEdgeConnections[l][0];
        // Column index
        int col                = rEdgeConnections[l][1];
        rFlippedMatrix[l][row] = rMatrix[l][col];
        rFlippedMatrix[l][col] = rMatrix[l][row];
    }
}