// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for inputs of Knapsack-MWMW problem
 */
class InputsKcMwm
{
private:
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int numProjects_; ///< number of projects
    int maxSizeNonzeroWeights_; ///< maximum # of non-zero weights over all departments

public:
    vector<impalib_type>         RewardTeam; ///< rewards of team equality constraints
    vector<vector<impalib_type>> Team2KnapsackM; ///< messages from teams to knapsack constraints
    vector<vector<int>>          TeamsWeightsPerDepartment; ///< weights of teams per each department
    vector<vector<impalib_type>> RewardProject; ///< reward of project equality constraint
    vector<int>                  MaxState; ///< vector of capacities of departments
    vector<vector<int>>          NonZeroWeightIndices; ///< indices of non-zero weights per each department

    void process_inputs(const impalib_type *, impalib_type *, const int *, const int *, const int *,
                        const impalib_type *, const int *); ///< process input of graphical model

    InputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS); ///< constructor
};

/**
 * Represents a class for outputs of Knapsack-MWMW problem
 */
class OutputsKcMwm
{
private:
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int numProjects_; ///< number of projects

public:
    vector<impalib_type> ExtrinsicOutputTeam; ///< extrinsic output of team equality constraints
    vector<impalib_type> IntrinsicOutMwm; ///< intrinsic outputs of project equality constraint 
    void                 intrinsic_out_mwm_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                                  vector<vector<impalib_type>> &); ///< calculate intrinsic outputs of project equality constraints
    void                 extrinsic_output_team_update(vector<vector<impalib_type>> &, vector<impalib_type> &); ///< calculate extrinsic output of team equality constraints
    OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS); ///< constructor
};

/**
 * Construct Input object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[in] MAX_SIZE_NON_ZERO_WEIGHTS: maximum number of connections between teams and departments
 * @param[out] numDepartments_: N_DEPARTMENTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] maxSizeNonzeroWeights_: MAX_SIZE_NON_ZERO_WEIGHTS
 * 
 */

InputsKcMwm::InputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS,
                         const int MAX_SIZE_NON_ZERO_WEIGHTS)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), numProjects_(N_PROJECTS),
      maxSizeNonzeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS)
{
    TeamsWeightsPerDepartment.reserve(numDepartments_);
    NonZeroWeightIndices.reserve(numDepartments_);
    RewardProject.reserve(numProjects_);
    Team2KnapsackM.reserve(numDepartments_);

    // Initialize Team2KnapsackM and TeamsWeightsPerDepartment vectors
    for (int department_index = 0; department_index < numDepartments_; department_index++)
    {
        Team2KnapsackM.push_back(vector<impalib_type>(numTeams_, zero_value));
        TeamsWeightsPerDepartment.push_back(vector<int>(numTeams_, 0));
    }

    // Initialize RewardProject vector
    for (int project_index = 0; project_index < numProjects_; project_index++)
    {
        RewardProject.push_back(vector<impalib_type>(numTeams_, zero_value));
    }
};

/**
 * Process inputs from python for the Knapsack-MWM problem
 *
 * @param[in] pREWARD_TEAM_PY: rewards of teams
 * @param[in] pTransition_model_py: teams to knapsack messages
 * @param[in] pTEAMS_WEIGHTS_PER_DEPARTMENT_PY: teams weights per department: weights of each team associated with all departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of connections between teams and departments
 * @param[in] p_NON_ZERO_WEIGHT_INDICES_PY: non-zero connections between teams and departments
 * @param[in] pREWARD_PROJECT_PY: rewards for project-team combination
 * @param[in] pMAX_STATE_PY: contains maximum capacity of departments
 * 
 */

void InputsKcMwm::process_inputs(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py,
                                 const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                                 const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY,
                                 const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY)
{

    // Copy reward values for each team and maximum states per department
    copy(pREWARD_TEAM_PY, pREWARD_TEAM_PY + numTeams_, back_inserter(RewardTeam));
    copy(pMAX_STATE_PY, pMAX_STATE_PY + numDepartments_, back_inserter(MaxState));

    for (int department_index = 0; department_index < numDepartments_; department_index++)
    {   
        // Populate NonZeroWeightIndices for the department
        NonZeroWeightIndices.push_back(vector<int>(pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index], 0));
        copy(pTransition_model_py + numTeams_ * department_index,
             pTransition_model_py + numTeams_ * (department_index + 1), Team2KnapsackM[department_index].begin());
        copy(pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * department_index,
             pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * (department_index + 1),
             TeamsWeightsPerDepartment[department_index].begin());
        copy(p_NON_ZERO_WEIGHT_INDICES_PY + maxSizeNonzeroWeights_ * department_index,
             p_NON_ZERO_WEIGHT_INDICES_PY + pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]
                 + maxSizeNonzeroWeights_ * department_index,
             NonZeroWeightIndices[department_index].begin());
    }
    for (int project_index = 0; project_index < numProjects_; project_index++)
    {
        // Copy reward values for each project-team combination
        copy(pREWARD_PROJECT_PY + numTeams_ * project_index, pREWARD_PROJECT_PY + numTeams_ * (project_index + 1),
             RewardProject[project_index].begin());
    }
}

/**
 * Construct outputs class for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[out] numDepartments_: N_DEPARTMENTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numProjects_: N_PROJECTS
 * 
 */

OutputsKcMwm::OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), numProjects_(N_PROJECTS)
{
    ExtrinsicOutputTeam.reserve(numTeams_);
    ExtrinsicOutputTeam.resize(numTeams_);
    fill(ExtrinsicOutputTeam.begin(), ExtrinsicOutputTeam.begin() + numTeams_, zero_value);

    IntrinsicOutMwm.reserve(numProjects_ * numTeams_);
    IntrinsicOutMwm.resize(numProjects_ * numTeams_);
    fill(IntrinsicOutMwm.begin(), IntrinsicOutMwm.begin() + numProjects_ * numTeams_, zero_value);
};

/**
 * Calculate instrinsic messages of project equality constraint for Knapsack-MWM problem
 *
 * @param[in] rOric2EqConstraintM: messages for ORIC to equality constraint of the projects
 * @param[in] rProject2EqConstraintM: messages from projects to equality constraints of the projects
 * @param[in] rRewardProject: rewards of project-team combinations
 * 
 */

void OutputsKcMwm::intrinsic_out_mwm_update(vector<vector<impalib_type>> &rOric2EqConstraintM,
                                            vector<vector<impalib_type>> &rProject2EqConstraintM,
                                            vector<vector<impalib_type>> &rRewardProject)
{
    for (int project_index = 0; project_index < rRewardProject.size(); project_index++)
    {
        for (int team_index = 0; team_index < rRewardProject[project_index].size(); team_index++)
        {
            // Calculate the intrinsic output for the project-team combination
            IntrinsicOutMwm[project_index + team_index + project_index * (numTeams_ - 1)] =
                rOric2EqConstraintM[project_index][team_index] + rProject2EqConstraintM[project_index][team_index]
                + rRewardProject[project_index][team_index];
        }
    }
}

/**
 * Calculate extrinsic messages of team equality constraint for Knapsack-MWM problem
 *
 * @param[in] rExtrinsicOutputDepartment: messages from knapsack constraints to team equality constraints
 * @param[in] rOric2TeamM: messages from ORIC to team equality constraints
 * 
 */

void OutputsKcMwm::extrinsic_output_team_update(vector<vector<impalib_type>> &rExtrinsicOutputDepartment,
                                                vector<impalib_type>         &rOric2TeamM)
{
    copy(rOric2TeamM.begin(), rOric2TeamM.end(), ExtrinsicOutputTeam.begin());

    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++)
    {
        // Calculate rExtrinsicOutputDepartment by adding rOric2TeamM
        transform(rExtrinsicOutputDepartment[department_index].begin(),
                  rExtrinsicOutputDepartment[department_index].end(), ExtrinsicOutputTeam.begin(),
                  ExtrinsicOutputTeam.begin(), std::plus<impalib_type>());
    }
}

/**
 * Represents a class for inputs of TSP
 */
class InputsTsp
{
private:
    int numNodes_; ///< number of nodes of TSP
    int numEdgeVariables_; ///< number of edges of TSP
    int numNodesPerEdge_ = 2; ///< number of nodes per edge

public:
    vector<vector<int>>          EdgeConnections; ///< constituent nodes per each edge
    vector<impalib_type>         CostEdgeVariable; ///< cost for each edge equality constraint
    vector<vector<impalib_type>> CostMatrix; ///< cost matrix (n_nodesxn_nodes)
    vector<vector<impalib_type>> EdgeDegreeConstraintCost; ///< cost matrix represented as num_edges x num_nodes to facilitate message updates
    vector<vector<impalib_type>> EdgeEc2DegreeConstraintM; ///< messages from edge equality constraint to degree constraint

    void process_inputs(const int *, const impalib_type *, const impalib_type *, impalib_type *, const impalib_type *); ///< process inputs of TSP graphical model

    InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES); ///< constructor
};

/**
 * Construct Input object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 * @param[out] numNodes_: NUM_NODES
 * @param[out] numEdgeVariables_: NUM_EDGE_VARIABLES
 * 
 */

InputsTsp::InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES)
    : numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES){
                            };

/**
 * Represents a class for outputs of TSP
 */
class OutputsTsp
{
private:
    int numNodes_; ///< number of nodes
    int numEdgeVariables_; ///< number of edges

public:
    vector<impalib_type> ExtrinsicOutputEdgeEc; ///< extrinsic output of edge equality constraint
    vector<impalib_type> IntrinsicOutputEdgeEc; ///< intrinsic output of edge equality constraint
    void                 extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>> &); ///< calculate extrinsic output of edge equality constraint for a relaxed TSP
    void                 extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>> &,
                                                                         vector<vector<impalib_type>> &); ///< calculate extrinsic output of edge equality constraint for augmented TSP
    void                 intrinsic_output_edge_ec_update(vector<impalib_type> &); ///< calculate intrinsic output of edge equality constraint for augmented TSP
    OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES); ///< constructor
};

/**
 * Construct Output object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 * @param[out] numNodes_: NUM_NODES
 * @param[out] numEdgeVariables_: NUM_EDGE_VARIABLES
 * 
 */

OutputsTsp::OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES)
    : numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES)
{

    // Reserve and initialize memory for extrinsic and intrinsic edge equality constraints outputs
    ExtrinsicOutputEdgeEc.reserve(numEdgeVariables_);
    ExtrinsicOutputEdgeEc.resize(numEdgeVariables_);
    IntrinsicOutputEdgeEc.reserve(numEdgeVariables_);
    IntrinsicOutputEdgeEc.resize(numEdgeVariables_);
    fill(ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.begin() + numEdgeVariables_, zero_value);
    fill(IntrinsicOutputEdgeEc.begin(), IntrinsicOutputEdgeEc.begin() + numEdgeVariables_, zero_value);
};

/**
 * Process inputs from python for TSP
 *
 * @param[in] pEDGE_CONNECTIONS_PY: contains all possible connections between nodes. Each edge has its constituent nodes
 * @param[in] pCOST_EDGE_VARIABLE_PY: cost for each possible connection between nodes. This has size of number of edges
 * @param[in] pCOST_MATRIX_PY: cost matrix of size number of nodes x number of nodes
 * @param[in] pEdge_ec_to_degree_constraint_m_py: messages from edges equality constraints to degree constraints
 * @param[in] pEDGE_DEGREE_CONSTRAINT_COST_PY: another matrix of costs that has size number of edges x number of nodes. 
 * Refer to example of TSP to understand how the various costs differ. The various various will facilitate computations in the IMPA
 * 
 */

void InputsTsp::process_inputs(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY,
                               const impalib_type *pCOST_MATRIX_PY, impalib_type *pEdge_ec_to_degree_constraint_m_py,
                               const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY)
{

    // Copy cost edge variable data
    copy(pCOST_EDGE_VARIABLE_PY, pCOST_EDGE_VARIABLE_PY + numEdgeVariables_, back_inserter(CostEdgeVariable));

    // Process input variables
    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++)
    {
        // Copy edge connections
        EdgeConnections.push_back(vector<int>(numNodesPerEdge_, 0));
        copy(pEDGE_CONNECTIONS_PY + numNodesPerEdge_ * edge_variable_index,
             pEDGE_CONNECTIONS_PY + numNodesPerEdge_ * (edge_variable_index + 1),
             EdgeConnections[edge_variable_index].begin());
        // Copy edge ec to degree constraint messages
        EdgeEc2DegreeConstraintM.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pEdge_ec_to_degree_constraint_m_py + numNodes_ * edge_variable_index,
             pEdge_ec_to_degree_constraint_m_py + numNodes_ * (edge_variable_index + 1),
             EdgeEc2DegreeConstraintM[edge_variable_index].begin());
        // Copy edge degree constraint cost
        EdgeDegreeConstraintCost.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pEDGE_DEGREE_CONSTRAINT_COST_PY + numNodes_ * edge_variable_index,
             pEDGE_DEGREE_CONSTRAINT_COST_PY + numNodes_ * (edge_variable_index + 1),
             EdgeDegreeConstraintCost[edge_variable_index].begin());
    }

    // Process cost matrix
    for (int node_index = 0; node_index < numNodes_; node_index++)
    {
        CostMatrix.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pCOST_MATRIX_PY + numNodes_ * node_index, pCOST_MATRIX_PY + numNodes_ * (node_index + 1),
             CostMatrix[node_index].begin());
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for relaxed TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to edge equality constraints
 * 
 */

void OutputsTsp::extrinsic_output_edge_ec_relaxed_graph_update(
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM)
{

    for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size();
         edge_variable_index++)
    {
        ExtrinsicOutputEdgeEc[edge_variable_index] =
            accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(),
                       rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for augmented TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to edge equality constraints
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints
 * 
 */

void OutputsTsp::extrinsic_output_edge_ec_augmented_graph_update(
    vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM)
{

    // Update extrinsic output for edge equality constraints based on messages from degree constraints to equality constraints
    for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size();
         edge_variable_index++)
    {
        ExtrinsicOutputEdgeEc[edge_variable_index] =
            accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(),
                       rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
    }

    // Update extrinsic output for edge equality constraints based on messages from subtour elimination constraints
    for (int subtour_constraint_index = 0; subtour_constraint_index < rSubtourConstraints2EdgeEcM.size();
         subtour_constraint_index++)
    {
        transform(rSubtourConstraints2EdgeEcM[subtour_constraint_index].begin(),
                  rSubtourConstraints2EdgeEcM[subtour_constraint_index].end(), ExtrinsicOutputEdgeEc.begin(),
                  ExtrinsicOutputEdgeEc.begin(), plus<impalib_type>());
    }
}

/**
 * Calculate output intrinsic messages of edge equality constraint for TSP
 *
 * @param[in] rCostEdgeVariable: cost of each edge equality constraint
 * 
 */

void OutputsTsp::intrinsic_output_edge_ec_update(vector<impalib_type> &rCostEdgeVariable)
{

    // Calculate intrinsic output for edge equality constraints by adding extrinsic output for edge equality constraints to cost edge variable
    transform(ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.end(), rCostEdgeVariable.begin(),
              IntrinsicOutputEdgeEc.begin(), plus<impalib_type>());
}