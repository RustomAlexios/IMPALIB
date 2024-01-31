// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class EqualityConstraintKcMwm
{
private:
    int numTeams_;
    int numProjects_;
    int numDepartments_;

public:
    void team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                           vector<impalib_type> &);

    void project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &);

    EqualityConstraintKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS);
};

EqualityConstraintKcMwm::EqualityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                        // numProjects_ = N_PROJECTS;
                                                        // numTeams_ = N_TEAMS;
                                                        // numDepartments_ = N_DEPARTMENTS;

                                                    };

void EqualityConstraintKcMwm::team_eq_constraint_to_oric_update(
    vector<vector<impalib_type>> &rExtrinsicOutputDepartment, vector<impalib_type> &rTeam2OricM,
    vector<impalib_type> &rewardTeam)
{

    vector<impalib_type> intermediate_team_to_oric_m(numTeams_, 0);

    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++)
    {
        transform(rExtrinsicOutputDepartment[department_index].begin(),
                  rExtrinsicOutputDepartment[department_index].end(), intermediate_team_to_oric_m.begin(),
                  intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
    }
    transform(intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.end(), rewardTeam.begin(),
              rTeam2OricM.begin(), std::plus<impalib_type>());
}

void EqualityConstraintKcMwm::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rProject2EqConstraintM,
                                                                   vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                                   vector<vector<impalib_type>> &rewardProject)
{
    for (int project_index = 0; project_index < rProject2EqConstraintM.size(); project_index++)
    {
        transform(rProject2EqConstraintM[project_index].begin(), rProject2EqConstraintM[project_index].end(),
                  rewardProject[project_index].begin(), rEqConstraint2OricM[project_index].begin(),
                  std::plus<impalib_type>());
    }
}

class EqualityConstraintTsp
{
private:
    int          numNodes_;
    int          numEdgeVariables_;
    bool         filteringFlag_;
    impalib_type alpha_;

public:
    void edge_ec_to_degree_constraint_relaxed_graph_update(vector<vector<int>> &, vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &);
    void flip_matrix(vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<impalib_type>> &);
    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_update(vector<vector<int>> &, vector<impalib_type> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<impalib_type>> &,
                                                                       vector<vector<int>> &);
    void                         edge_ec_to_degree_constraint_augmented_graph_update(vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &, vector<vector<int>> &,
                                                                                     vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &);

    EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool FILTERING_FLAG,
                          const impalib_type ALPHA);
};

EqualityConstraintTsp::EqualityConstraintTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES),
      numEdgeVariables_(NUM_EDGE_VARIABLES){
          // filteringFlag_ = FILTERING_FLAG;
          // alpha_ = ALPHA;
          // numNodes_ = NUM_NODES;
          // numEdgeVariables_ = NUM_EDGE_VARIABLES;
      };

void EqualityConstraintTsp::edge_ec_to_degree_constraint_relaxed_graph_update(
    vector<vector<int>> &rEdgeConnections, vector<vector<impalib_type>> &rEdgeDegreeConstraintCost,
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM)
{

    vector<vector<impalib_type>> flipped_degree_constraint_to_eq_constraint_m = rDegreeConstraint2EqConstraintM;

    flip_matrix(rDegreeConstraint2EqConstraintM, rEdgeConnections, flipped_degree_constraint_to_eq_constraint_m);

    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++)
    {
        transform(flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].begin(),
                  flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].end(),
                  rEdgeDegreeConstraintCost[edge_variable_index].begin(),
                  rEdgeEc2DegreeConstraintM[edge_variable_index].begin(), plus<impalib_type>());
    }
}

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
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
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
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
        }

        for (size_t index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
             index_subtour_constraint++)
        {

            vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);
            for (size_t i = 0; i < rDeltaSIndicesList[index_subtour_constraint].size(); i++)
            {
                edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[index_subtour_constraint][i]] =
                    combined_subtour_constraints_to_edge_ec_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + rCostEdgeVaribale[rDeltaSIndicesList[index_subtour_constraint][i]]
                    - rSubtourConstraints2EdgeEcM[index_subtour_constraint]
                                                 [rDeltaSIndicesList[index_subtour_constraint][i]];
            }
            edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
        }
    }
    return edge_ec_to_subtour_constraints_list;
}

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

        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][0]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][1]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][0]];

        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][1]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][0]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][1]];
    }
}

void EqualityConstraintTsp::flip_matrix(vector<vector<impalib_type>> &rMatrix, vector<vector<int>> &rEdgeConnections,
                                        vector<vector<impalib_type>> &rFlippedMatrix)
{

    for (size_t l = 0; l < rEdgeConnections.size(); ++l)
    {
        int row                = rEdgeConnections[l][0];
        int col                = rEdgeConnections[l][1];
        rFlippedMatrix[l][row] = rMatrix[l][col];
        rFlippedMatrix[l][col] = rMatrix[l][row];
    }
}