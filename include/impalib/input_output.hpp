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
class InputsKcMwm {
   private:
    int numTeams_;               ///< number of teams
    int numDepartments_;         ///< number of departments
    int numProjects_;            ///< number of projects
    int maxSizeNonzeroWeights_;  ///< maximum # of non-zero weights over all departments

   public:
    vector<impalib_type> RewardTeam;                ///< rewards of team equality constraints
    vector<vector<impalib_type>> Team2KnapsackM;    ///< messages from teams to knapsack constraints
    vector<vector<int>> TeamsWeightsPerDepartment;  ///< weights of teams per each department
    vector<vector<impalib_type>> RewardProject;     ///< reward of project equality constraint
    vector<int> MaxState;                           ///< vector of capacities of departments
    vector<vector<int>> NonZeroWeightIndices;       ///< indices of non-zero weights per each department

    void process_inputs(const impalib_type *, impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *);  ///< process input of graphical model

    InputsKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS, int MAX_SIZE_NON_ZERO_WEIGHTS);  ///< constructor
};

/**
 * Represents a class for outputs of Knapsack-MWMW problem
 */
class OutputsKcMwm {
   private:
    int numTeams_;        ///< number of teams
    int numDepartments_;  ///< number of departments
    int numProjects_;     ///< number of projects

   public:
    vector<impalib_type> ExtrinsicOutputTeam;  ///< extrinsic output of team equality constraints
    vector<impalib_type> IntrinsicOutMwm;      ///< intrinsic outputs of project equality constraint
    void intrinsic_out_mwm_update(const vector<vector<impalib_type>> &, const vector<vector<impalib_type>> &,
                                  const vector<vector<impalib_type>> &);                              ///< calculate intrinsic outputs of project equality constraints
    void extrinsic_output_team_update(vector<vector<impalib_type>> &, vector<impalib_type> &);  ///< calculate extrinsic output of team equality constraints
    OutputsKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);             ///< constructor
};

/**
 * Construct Input object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[in] MAX_SIZE_NON_ZERO_WEIGHTS: maximum number of connections between teams and departments
 *
 */

inline InputsKcMwm::InputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), numProjects_(N_PROJECTS), maxSizeNonzeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS) {
    TeamsWeightsPerDepartment.reserve(numDepartments_);
    NonZeroWeightIndices.reserve(numDepartments_);
    RewardProject.reserve(numProjects_);
    Team2KnapsackM.reserve(numDepartments_);

    for (int department_index = 0; department_index < numDepartments_; department_index++) {
        Team2KnapsackM.push_back(vector<impalib_type>(numTeams_, zero_value));
        TeamsWeightsPerDepartment.push_back(vector<int>(numTeams_, 0));
    }
    for (int project_index = 0; project_index < numProjects_; project_index++) {
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

inline void InputsKcMwm::process_inputs(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                                 const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
    copy(pREWARD_TEAM_PY, pREWARD_TEAM_PY + numTeams_, back_inserter(RewardTeam));
    copy(pMAX_STATE_PY, pMAX_STATE_PY + numDepartments_, back_inserter(MaxState));

    for (int department_index = 0; department_index < numDepartments_; department_index++) {
        NonZeroWeightIndices.push_back(vector<int>(pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index], 0));
        copy(pTransition_model_py + numTeams_ * department_index, pTransition_model_py + numTeams_ * (department_index + 1), Team2KnapsackM[department_index].begin());
        copy(pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * department_index, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * (department_index + 1),
             TeamsWeightsPerDepartment[department_index].begin());
        copy(p_NON_ZERO_WEIGHT_INDICES_PY + maxSizeNonzeroWeights_ * department_index,
             p_NON_ZERO_WEIGHT_INDICES_PY + pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index] + maxSizeNonzeroWeights_ * department_index, NonZeroWeightIndices[department_index].begin());
    }
    for (int project_index = 0; project_index < numProjects_; project_index++) {
        // Copy reward values for each project-team combination
        copy(pREWARD_PROJECT_PY + numTeams_ * project_index, pREWARD_PROJECT_PY + numTeams_ * (project_index + 1), RewardProject[project_index].begin());
    }
}

/**
 * Construct outputs class for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 *
 */

inline OutputsKcMwm::OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), numProjects_(N_PROJECTS) {
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

inline void OutputsKcMwm::intrinsic_out_mwm_update(const vector<vector<impalib_type>> &rOric2EqConstraintM, const vector<vector<impalib_type>> &rProject2EqConstraintM, const vector<vector<impalib_type>> &rRewardProject) {
    for (int project_index = 0; project_index < rRewardProject.size(); project_index++) {
        for (int team_index = 0; team_index < rRewardProject[project_index].size(); team_index++) {
            IntrinsicOutMwm[project_index + team_index + project_index * (numTeams_ - 1)] =
                rOric2EqConstraintM[project_index][team_index] + rProject2EqConstraintM[project_index][team_index] + rRewardProject[project_index][team_index];
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

inline void OutputsKcMwm::extrinsic_output_team_update(vector<vector<impalib_type>> &rExtrinsicOutputDepartment, vector<impalib_type> &rOric2TeamM) {
    copy(rOric2TeamM.begin(), rOric2TeamM.end(), ExtrinsicOutputTeam.begin());

    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++) {
        transform(rExtrinsicOutputDepartment[department_index].begin(), rExtrinsicOutputDepartment[department_index].end(), ExtrinsicOutputTeam.begin(), ExtrinsicOutputTeam.begin(),
                  std::plus<impalib_type>());
    }
}

/**
 * Represents a class for inputs of TSP
 */
class InputsTsp {
   private:
    int numNodes_;             ///< number of nodes of TSP
    int numEdgeVariables_;     ///< number of edges of TSP
    int numNodesPerEdge_ = 2;  ///< number of nodes per edge

   public:
    vector<vector<int>> EdgeConnections;                    ///< constituent nodes per each edge
    vector<impalib_type> CostEdgeVariable;                  ///< cost for each edge equality constraint
    vector<vector<impalib_type>> CostMatrix;                ///< cost matrix (n_nodesxn_nodes)
    vector<vector<impalib_type>> EdgeDegreeConstraintCost;  ///< cost matrix represented as num_edges x num_nodes to facilitate message updates
    vector<vector<impalib_type>> EdgeEc2DegreeConstraintM;  ///< messages from edge equality constraint to degree constraint

    void process_inputs(const int *, const impalib_type *, const impalib_type *, impalib_type *, const impalib_type *);  ///< process inputs of TSP graphical model

    InputsTsp(int NUM_NODES, int NUM_EDGE_VARIABLES);  ///< constructor
};

/**
 * Construct Input object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 *
 */

inline InputsTsp::InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES){};

/**
 * Represents a class for outputs of TSP
 */
class OutputsTsp {
   private:
    int numNodes_;          ///< number of nodes
    int numEdgeVariables_;  ///< number of edges

   public:
    vector<impalib_type> ExtrinsicOutputEdgeEc;                                          ///< extrinsic output of edge equality constraint
    vector<impalib_type> IntrinsicOutputEdgeEc;                                          ///< intrinsic output of edge equality constraint
    void extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>> &);  ///< calculate extrinsic output of edge equality constraint for a relaxed TSP
    void extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>> &,
                                                         vector<vector<impalib_type>> &);  ///< calculate extrinsic output of edge equality constraint for augmented TSP
    void intrinsic_output_edge_ec_update(vector<impalib_type> &);                          ///< calculate intrinsic output of edge equality constraint for augmented TSP
    OutputsTsp(int NUM_NODES, int NUM_EDGE_VARIABLES);                         ///< constructor
};

/**
 * Construct Output object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 *
 */

inline OutputsTsp::OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES) {
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

inline void InputsTsp::process_inputs(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY, impalib_type *pEdge_ec_to_degree_constraint_m_py,
                               const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
    copy(pCOST_EDGE_VARIABLE_PY, pCOST_EDGE_VARIABLE_PY + numEdgeVariables_, back_inserter(CostEdgeVariable));

    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++) {
        EdgeConnections.push_back(vector<int>(numNodesPerEdge_, 0));
        copy(pEDGE_CONNECTIONS_PY + numNodesPerEdge_ * edge_variable_index, pEDGE_CONNECTIONS_PY + numNodesPerEdge_ * (edge_variable_index + 1), EdgeConnections[edge_variable_index].begin());

        EdgeEc2DegreeConstraintM.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pEdge_ec_to_degree_constraint_m_py + numNodes_ * edge_variable_index, pEdge_ec_to_degree_constraint_m_py + numNodes_ * (edge_variable_index + 1),
             EdgeEc2DegreeConstraintM[edge_variable_index].begin());

        EdgeDegreeConstraintCost.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pEDGE_DEGREE_CONSTRAINT_COST_PY + numNodes_ * edge_variable_index, pEDGE_DEGREE_CONSTRAINT_COST_PY + numNodes_ * (edge_variable_index + 1),
             EdgeDegreeConstraintCost[edge_variable_index].begin());
    }

    for (int node_index = 0; node_index < numNodes_; node_index++) {
        CostMatrix.push_back(vector<impalib_type>(numNodes_, zero_value));
        copy(pCOST_MATRIX_PY + numNodes_ * node_index, pCOST_MATRIX_PY + numNodes_ * (node_index + 1), CostMatrix[node_index].begin());
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for relaxed TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to edge equality constraints
 *
 */

inline void OutputsTsp::extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM) {
    for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size(); edge_variable_index++) {
        ExtrinsicOutputEdgeEc[edge_variable_index] = accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(), rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for augmented TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to edge equality constraints
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints
 *
 */

inline void OutputsTsp::extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM, vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM) {
    for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size(); edge_variable_index++) {
        ExtrinsicOutputEdgeEc[edge_variable_index] = accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(), rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
    }

    for (int subtour_constraint_index = 0; subtour_constraint_index < rSubtourConstraints2EdgeEcM.size(); subtour_constraint_index++) {
        transform(rSubtourConstraints2EdgeEcM[subtour_constraint_index].begin(), rSubtourConstraints2EdgeEcM[subtour_constraint_index].end(), ExtrinsicOutputEdgeEc.begin(),
                  ExtrinsicOutputEdgeEc.begin(), plus<impalib_type>());
    }
}

/**
 * Calculate output intrinsic messages of edge equality constraint for TSP
 *
 * @param[in] rCostEdgeVariable: cost of each edge equality constraint
 *
 */

inline void OutputsTsp::intrinsic_output_edge_ec_update(vector<impalib_type> &rCostEdgeVariable) {
    transform(ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.end(), rCostEdgeVariable.begin(), IntrinsicOutputEdgeEc.begin(), plus<impalib_type>());
}

/**
 * Represents a class for inputs of K-SAT
 */
class InputsKsat {
   private:
    int numVariables_;      ///< total number of variables
    int numConstraints_;    ///< number of constraints
    int kVariable_;         ///< number of variables per constraint
    int numUsedVariables_;  ///< number of variables used to construct the formula

   public:
    vector<int> UsedVariables;                                ///< used variables to construct the formula
    vector<impalib_type> IncomingMetricsCost;                 ///< incoming metrics for varibales
    vector<vector<int>> ConstraintsConnections;               ///< connections to variables for each constraint
    vector<vector<int>> ConstraintsConnectionsType;           ///< types of connections to variables for each constraint
    vector<vector<int>> VariablesConnections;                 ///< connections to constraints for each variable
    vector<int> VariablesConnectionsSizes;                    ///< sizes of connections to constraints for each variable
    vector<vector<impalib_type>> VariableEc2KsatConstraintM;  ///< messages from variables equality constraints to k-sat constraints

    void process_inputs(const int *, const int *, const int *, const int *, const int *, const impalib_type *, impalib_type *);  ///< process inputs from python

    InputsKsat(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE, int NUM_USED_VARIABLES);  ///< constructor
};

/**
 * Construct Input object for k-sat problem
 *
 * @param[in] NUM_VARIABLES: total number of variables
 * @param[in] NUM_CONSTRAINTS: number of constraints
 * @param[in] K_VARIABLE: ///< number of variables per constraint
 * @param[in] NUM_USED_VARIABLES: number of variables used to construct the formula
 *
 */

inline InputsKsat::InputsKsat(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const int NUM_USED_VARIABLES)
    : numVariables_(NUM_VARIABLES),
      numConstraints_(NUM_CONSTRAINTS),
      kVariable_(K_VARIABLE),
      numUsedVariables_(NUM_USED_VARIABLES){

      };

/**
 * Represents a class for outputs of K-SAT
 */
class OutputsKsat {
   private:
    int numVariables_;    ///< total number of variables
    int numConstraints_;  ///< number of constraints
    int kVariable_;       ///< number of variables per constraint

   public:
    vector<impalib_type> ExtrinsicOutputVariableEc;                                         ///< extrinsic messages of variables equality constraints
    OutputsKsat(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE);  ///< constructor
    void update_extrinsic(const vector<vector<impalib_type>> &);                                  ///< calculate extrinsic messages of variables equality constraints
};

/**
 * Construct Output object for k-sat problem
 *
 * @param[in] NUM_VARIABLES: total number of variables
 * @param[in] NUM_CONSTRAINTS: number of constraints
 * @param[in] K_VARIABLE: ///< number of variables per constraint
 *
 */

inline OutputsKsat::OutputsKsat(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE)
    : numVariables_(NUM_VARIABLES), numConstraints_(NUM_CONSTRAINTS), kVariable_(K_VARIABLE), ExtrinsicOutputVariableEc(numVariables_, zero_value){};

/**
 * Process inputs from python for the K-SAT problem
 *
 * @param[in] pUSED_VARIABLES_PY: variables used in building the formula
 * @param[in] pVARIABLES_CONNECTIONS_PY: connections to constraints for each variable
 * @param[in] pVARIABLES_CONNECTIONS_SIZES: size of connections to constraints for each variable
 * @param[in] pCONSTRAINTS_CONNECTIONS: connections to variables for each constraint
 * @param[in] pCONSTRAINTS_CONNECTIONS_TYPE: types of connections to variables for each constraint
 * @param[in] pINCOMING_METRICS_COST: incoming metrics for each variable
 * @param[in] pVariable_ec_to_ksat_constraint_m_py: initial messages from variables equality constraints to k-sat constraints
 *
 */

inline void InputsKsat::process_inputs(const int *pUSED_VARIABLES_PY, const int *pVARIABLES_CONNECTIONS_PY, const int *pVARIABLES_CONNECTIONS_SIZES, const int *pCONSTRAINTS_CONNECTIONS,
                                const int *pCONSTRAINTS_CONNECTIONS_TYPE, const impalib_type *pINCOMING_METRICS_COST, impalib_type *pVariable_ec_to_ksat_constraint_m_py) {
    copy(pUSED_VARIABLES_PY, pUSED_VARIABLES_PY + numUsedVariables_, back_inserter(UsedVariables));

    copy(pINCOMING_METRICS_COST, pINCOMING_METRICS_COST + numVariables_, back_inserter(IncomingMetricsCost));

    copy(pVARIABLES_CONNECTIONS_SIZES, pVARIABLES_CONNECTIONS_SIZES + numVariables_, back_inserter(VariablesConnectionsSizes));

    for (int i = 0; i < numConstraints_; i++) {
        ConstraintsConnections.push_back(vector<int>(kVariable_, 0));

        copy(pCONSTRAINTS_CONNECTIONS + kVariable_ * i, pCONSTRAINTS_CONNECTIONS + kVariable_ * (i + 1), ConstraintsConnections[i].begin());

        ConstraintsConnectionsType.push_back(vector<int>(kVariable_, 0));

        copy(pCONSTRAINTS_CONNECTIONS_TYPE + kVariable_ * i, pCONSTRAINTS_CONNECTIONS_TYPE + kVariable_ * (i + 1), ConstraintsConnectionsType[i].begin());

        VariableEc2KsatConstraintM.push_back(vector<impalib_type>(numVariables_, zero_value));

        copy(pVariable_ec_to_ksat_constraint_m_py + numVariables_ * i, pVariable_ec_to_ksat_constraint_m_py + numVariables_ * (i + 1), VariableEc2KsatConstraintM[i].begin());
    }

    int conx_size_old = 0;

    for (int j = 0; j < numVariables_; j++) {
        if (find(UsedVariables.begin(), UsedVariables.end(), j) != UsedVariables.end()) {
            int conx_size = VariablesConnectionsSizes[j];

            VariablesConnections.push_back(vector<int>(conx_size, 0));

            copy(pVARIABLES_CONNECTIONS_PY + conx_size_old, pVARIABLES_CONNECTIONS_PY + conx_size_old + conx_size, VariablesConnections[j].begin());

            conx_size_old += conx_size;

        } else {
            VariablesConnections.push_back(vector<int>());
        }
    }
}

/**
 * Calculate output extrinsic messages of variable equality constraints for K-SAT
 *
 * @param[in] rKsatConstraint2EqConstraintM: messages from k-sat constraint to variable equality constraint
 *
 */

inline void OutputsKsat::update_extrinsic(const vector<vector<impalib_type>> &rKsatConstraint2EqConstraintM) {
    // Sum all messages coming into equality constraint except the
    // incoming message on the edge of interest
    for (int i = 0; i < numVariables_; i++) {
        for (int j = 0; j < numConstraints_; j++) {
            ExtrinsicOutputVariableEc[i] += rKsatConstraint2EqConstraintM[j][i];
        }
    }
}


/**
 * Represents a class for inputs of MOBARP
 */
class InputsMOBARP {
   private:
    int numFixedTx_;
    int numMobileTx_;
    int numBands_;
    int numTimeSteps_;
    int numRxLocs_;
    int numMobileTxLocs_;

   public:
    vector<vector<impalib_type>> RCosts;
    vector<vector<impalib_type>> ZCosts;
    vector<vector<vector<impalib_type>>> MobileTxCosts;
    vector<vector<vector<impalib_type>>> FixedTxCosts;
    vector<vector<vector<vector<int>>>> ConnectivityFixedTx;
    vector<vector<vector<vector<int>>>> ConnectivityMobileTx;
    vector<vector<impalib_type>> REqConst2AuxiliaryConstM;
    vector<vector<vector<impalib_type>>> MobileXEqConst2MobileCapacConstM;
    vector<vector<vector<impalib_type>>> FixedXEqConst2FixedCapacConstM;
    vector<int> FixedCapacConstraints; 
    vector<int> MobileCapacConstraints;
    vector<vector<int>> ConxMobTxPerNumMobTxLocs;
    vector<vector<int>> ConxFixedTxPerNumRXLocs;
    vector<vector<vector<vector<int>>>> ConxMobTxRx;
    vector<vector<vector<vector<int>>>> TransposedConxMobTxRx;
    vector<vector<int>> ConxMobTxR;
    vector<vector<vector<vector<int>>>> TempConxMobTxRx;
    vector<vector<int>> TempReshapedConxMobTxRx;
    vector<vector<int>> TempReshapedConnectivityFixedTx;


    void process_inputs(impalib_type *, impalib_type *, impalib_type *, const impalib_type *, const impalib_type *, const impalib_type *, const impalib_type *, const int *, const int *,
                            const int *, const int *, const int *, const int *, const int *);  ///< process inputs from python

    InputsMOBARP(int NUM_FIXED_TX, int NUM_MOBILE_TX, int NUM_BANDS, int NUM_TIME_STEPS, int NUM_RX_LOCS, int NUM_MOBILE_TX_LOCS);  ///< constructor
};


inline InputsMOBARP::InputsMOBARP(const int NUM_FIXED_TX, const int NUM_MOBILE_TX, const int NUM_BANDS, const int NUM_TIME_STEPS, const int NUM_RX_LOCS, const int NUM_MOBILE_TX_LOCS)
        :numFixedTx_(NUM_FIXED_TX),
        numMobileTx_(NUM_MOBILE_TX),
        numBands_(NUM_BANDS),
        numTimeSteps_(NUM_TIME_STEPS),
        numRxLocs_(NUM_RX_LOCS),
        numMobileTxLocs_(NUM_MOBILE_TX_LOCS),
        TransposedConxMobTxRx(vector<vector<vector<vector<int>>>>(NUM_BANDS*NUM_MOBILE_TX, vector<vector<vector<int>>>(NUM_RX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_MOBILE_TX_LOCS, 0))))),
        ConxMobTxR(vector<vector<int>>(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_TIME_STEPS, 0))),
        TempConxMobTxRx(vector<vector<vector<vector<int>>>>(NUM_RX_LOCS, vector<vector<vector<int>>>(NUM_MOBILE_TX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))))),
        TempReshapedConxMobTxRx(vector<vector<int>>(NUM_MOBILE_TX_LOCS*NUM_BANDS*NUM_MOBILE_TX, vector<int>(NUM_RX_LOCS*NUM_TIME_STEPS, 0))),
        TempReshapedConnectivityFixedTx(vector<vector<int>>(NUM_FIXED_TX*NUM_BANDS, vector<int>(NUM_TIME_STEPS*NUM_RX_LOCS, 0))){};

inline void InputsMOBARP::process_inputs(impalib_type *pFixed_x_eq_const_to_fixed_capac_const_m_py, impalib_type *pMobile_x_eq_const_to_mobile_capac_const_m_py,
                                            impalib_type *pR_eq_const_to_auxiliary_const_m_py, const impalib_type *pFIXED_X_COSTS_PY, const impalib_type *pMOBILE_X_COSTS_PY,
                                            const impalib_type *pZ_COSTS_PY, const impalib_type *pR_COSTS_PY, const int *pCONNECTIVITY_FIXED_TX_PY,
                                            const int *pCONNECTIVITY_MOBILE_TX_PY, const int *pFIXED_CAPACITY_CONSTRAINTS_PY, const int *pMOBILE_CAPACITY_CONSTRAINTS_PY,
                                            const int *pCONX_MOB_TX_PER_NUM_MOB_TX_LOCS_PY, const int *pCONX_FIXED_TX_PER_NUM_RX_LOCS_PY, const int *pCONX_MOB_TX_RX_PY) {
    
    for (int i=0; i< numMobileTx_; i++){
        RCosts.push_back(vector<impalib_type>(numMobileTxLocs_, 0));
        copy(pR_COSTS_PY + numMobileTxLocs_ * i, pR_COSTS_PY + numMobileTxLocs_ * (i + 1), RCosts[i].begin());
    }

    for (int i=0; i< numMobileTx_*numBands_*numTimeSteps_; i++){
        ZCosts.push_back(vector<impalib_type>(numMobileTxLocs_, 0));
        copy(pZ_COSTS_PY + numMobileTxLocs_ * i, pZ_COSTS_PY + numMobileTxLocs_ * (i + 1), ZCosts[i].begin());
    }

    for (int i=0; i< numMobileTx_; i++){
        MobileTxCosts.push_back(vector<vector<impalib_type>>(numBands_, vector<impalib_type>(numTimeSteps_, 0)));
        for (int j=0; j< numBands_; j++){
            copy(pMOBILE_X_COSTS_PY + numTimeSteps_ * j + numBands_*numTimeSteps_*i, pMOBILE_X_COSTS_PY + numTimeSteps_ * (j + 1) + numBands_*numTimeSteps_*i, MobileTxCosts[i][j].begin());
        }
    }

    for (int i=0; i< numFixedTx_; i++){
        FixedTxCosts.push_back(vector<vector<impalib_type>>(numBands_, vector<impalib_type>(numTimeSteps_, 0)));
        for (int j=0; j< numBands_; j++){
            copy(pFIXED_X_COSTS_PY + numTimeSteps_ * j + numBands_*numTimeSteps_*i, pFIXED_X_COSTS_PY + numTimeSteps_ * (j + 1) + numBands_*numTimeSteps_*i, FixedTxCosts[i][j].begin());
        }
    }

    for (int i=0; i< numFixedTx_; i++){
        ConnectivityFixedTx.push_back(vector<vector<vector<int>>>(numBands_, vector<vector<int>>(numTimeSteps_, vector<int>(numRxLocs_, 0))));
        for (int j=0; j<numBands_; j++){
            for (int k=0; k< numTimeSteps_; k++){
                copy(pCONNECTIVITY_FIXED_TX_PY + numRxLocs_*k + numRxLocs_*numTimeSteps_ * j + numRxLocs_*numBands_*numTimeSteps_*i, pCONNECTIVITY_FIXED_TX_PY + numRxLocs_*(k+1) + numRxLocs_*numTimeSteps_ * j + numRxLocs_*numBands_*numTimeSteps_*i, ConnectivityFixedTx[i][j][k].begin());
        }
        }
    }

    for (int n=0; n< numMobileTxLocs_; n++){
        ConnectivityMobileTx.push_back(vector<vector<vector<int>>>(numBands_, vector<vector<int>>(numTimeSteps_, vector<int>(numRxLocs_, 0))));
        for (int j=0; j<numBands_; j++){
            for (int k=0; k< numTimeSteps_; k++){
                copy(pCONNECTIVITY_MOBILE_TX_PY + numRxLocs_*k + numRxLocs_*numTimeSteps_ * j + numRxLocs_*numBands_*numTimeSteps_*n, pCONNECTIVITY_MOBILE_TX_PY + numRxLocs_*(k+1) + numRxLocs_*numTimeSteps_ * j + numRxLocs_*numBands_*numTimeSteps_*n, ConnectivityMobileTx[n][j][k].begin());
        }
        }
    }

    for (int i=0; i< numMobileTx_*numMobileTxLocs_; i++){
        REqConst2AuxiliaryConstM.push_back(vector<impalib_type>(numBands_*numTimeSteps_, 0));
        copy(pR_eq_const_to_auxiliary_const_m_py + numBands_*numTimeSteps_ * i, pR_eq_const_to_auxiliary_const_m_py + numBands_*numTimeSteps_ * (i + 1), REqConst2AuxiliaryConstM[i].begin());
    }

    for (int i=0; i< numMobileTx_; i++){
        MobileXEqConst2MobileCapacConstM.push_back(vector<vector<impalib_type>>(numBands_, vector<impalib_type>(numTimeSteps_, 0)));
        for (int j=0; j< numBands_; j++){
            copy(pMobile_x_eq_const_to_mobile_capac_const_m_py + numTimeSteps_ * j + numBands_*numTimeSteps_*i, pMobile_x_eq_const_to_mobile_capac_const_m_py + numTimeSteps_ * (j + 1) + numBands_*numTimeSteps_*i, MobileXEqConst2MobileCapacConstM[i][j].begin());
        }
    }

    for (int i=0; i< numFixedTx_; i++){
        FixedXEqConst2FixedCapacConstM.push_back(vector<vector<impalib_type>>(numBands_, vector<impalib_type>(numTimeSteps_, 0)));
        for (int j=0; j< numBands_; j++){
            copy(pFixed_x_eq_const_to_fixed_capac_const_m_py + numTimeSteps_ * j + numBands_*numTimeSteps_*i, pFixed_x_eq_const_to_fixed_capac_const_m_py + numTimeSteps_ * (j + 1) + numBands_*numTimeSteps_*i, FixedXEqConst2FixedCapacConstM[i][j].begin());
        }
    }

    copy(pFIXED_CAPACITY_CONSTRAINTS_PY, pFIXED_CAPACITY_CONSTRAINTS_PY + numFixedTx_, back_inserter(FixedCapacConstraints));

    copy(pMOBILE_CAPACITY_CONSTRAINTS_PY, pMOBILE_CAPACITY_CONSTRAINTS_PY + numMobileTx_, back_inserter(MobileCapacConstraints));

    for (int i=0; i< numMobileTx_*numBands_*numTimeSteps_; i++){
        ConxMobTxPerNumMobTxLocs.push_back(vector<int>(numMobileTxLocs_, 0));
        copy(pCONX_MOB_TX_PER_NUM_MOB_TX_LOCS_PY + numMobileTxLocs_ * i, pCONX_MOB_TX_PER_NUM_MOB_TX_LOCS_PY + numMobileTxLocs_ * (i + 1), ConxMobTxPerNumMobTxLocs[i].begin());
    }

    for (int i=0; i< numFixedTx_*numBands_*numTimeSteps_; i++){
        ConxFixedTxPerNumRXLocs.push_back(vector<int>(numRxLocs_, 0));
        copy(pCONX_FIXED_TX_PER_NUM_RX_LOCS_PY + numRxLocs_ * i, pCONX_FIXED_TX_PER_NUM_RX_LOCS_PY + numRxLocs_ * (i + 1), ConxFixedTxPerNumRXLocs[i].begin());
    }

    //self.conx_mob_tx_rx = np.concatenate([np.transpose(self.connectivity_mobile_tx, (2,3,0,1))]*self.num_mobile_tx, axis=3)

    for (int k=0; k< numTimeSteps_; k++){
        ConxMobTxRx.push_back(vector<vector<vector<int>>>(numRxLocs_, vector<vector<int>>(numMobileTxLocs_, vector<int>(numBands_*numMobileTx_, 0))));
        for (int l=0; l<numRxLocs_; l++){
            for (int n=0; n< numMobileTxLocs_; n++){
                copy(pCONX_MOB_TX_RX_PY + numBands_*numMobileTx_*n + numBands_*numMobileTx_*numMobileTxLocs_*l + numBands_*numMobileTx_*numMobileTxLocs_*numRxLocs_*k, pCONX_MOB_TX_RX_PY + numBands_*numMobileTx_*(n+1) + numBands_*numMobileTx_*numMobileTxLocs_*l + numBands_*numMobileTx_*numMobileTxLocs_*numRxLocs_*k, ConxMobTxRx[k][l][n].begin());
        }
        }
    }

    // temp_connectivity_mobile_tx_rx = self.conx_mob_tx_rx.transpose(3, 1, 0, 2)
    //need to add this expression here instead of update_equality_constraint.hpp: reshaped_temp_connectivity_mobile_tx_rx = np.array([temp_connectivity_mobile_tx_rx[:, i].flatten() for i in range(self.num_rx_locs)]).transpose()
    
    for (int k=0; k< numTimeSteps_; k++){
        for (int l=0; l<numRxLocs_; l++){
            for (int n=0; n< numMobileTxLocs_; n++){
                for (int j_i=0; j_i <numBands_*numMobileTx_; j_i++){
                    TransposedConxMobTxRx[j_i][l][k][n] = ConxMobTxRx[k][l][n][j_i];
                }
            }
        }
    }

    //modelEqConstraint_.r_eq_const_activation
    vector<vector<vector<int>>> reshaped_1(numMobileTx_, vector<vector<int>>(numBands_*numTimeSteps_, vector<int>(numMobileTxLocs_, 0)));

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
            for (size_t n = 0; n < numMobileTxLocs_; n++) {
                reshaped_1[i][j_k][n] = ConxMobTxPerNumMobTxLocs[i * numBands_ * numTimeSteps_ + j_k][n];
            }
        }
    }

    vector<vector<vector<int>>> reshaped_2(numMobileTx_, vector<vector<int>>(numMobileTxLocs_, vector<int>(numBands_*numTimeSteps_, 0)));


    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t n = 0; n < numMobileTxLocs_; n++) {
            for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
                reshaped_2[i][n][j_k] = reshaped_1[i][j_k][n];
            }
        }
    }

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t n = 0; n < numMobileTxLocs_; n++) {
            for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
                ConxMobTxR[i * numMobileTxLocs_ + n][j_k] = reshaped_2[i][n][j_k];
            }
        }
    }

    // temp_conx_mob_tx_rx = self.conx_mob_tx_rx.transpose(1,2,0,3)

    for (int k=0; k< numTimeSteps_; k++){
        for (int l=0; l<numRxLocs_; l++){
            for (int n=0; n< numMobileTxLocs_; n++){
                for (int j_i=0; j_i <numBands_*numMobileTx_; j_i++){
                    TempConxMobTxRx[l][n][k][j_i] = ConxMobTxRx[k][l][n][j_i];
                }
            }
        }
    }
    

    //TempReshapedConxMobTxRx(vector<vector<int>>(NUM_MOBILE_TX_LOCS*NUM_BANDS*NUM_MOBILE_TX, vector<int>(NUM_RX_LOCS*NUM_TIME_STEPS, 0)))
    // temp_reshaped_conx_mob_tx_rx = np.swapaxes(temp_conx_mob_tx_rx, 0, 3).reshape(-1,self.num_time_steps*self.num_rx_locs)

    vector<vector<vector<vector<int>>>> swaped_temp_conx_mob_tx_rx(numBands_*numMobileTx_, vector<vector<vector<int>>>(numMobileTxLocs_, vector<vector<int>>(numTimeSteps_, vector<int>(numRxLocs_, 0))));
    
    for (int k=0; k< numTimeSteps_; k++){
        for (int l=0; l<numRxLocs_; l++){
            for (int n=0; n< numMobileTxLocs_; n++){
                for (int j_i=0; j_i <numBands_*numMobileTx_; j_i++){
                    swaped_temp_conx_mob_tx_rx[j_i][n][k][l] = TempConxMobTxRx[l][n][k][j_i];
                }
            }
        }
    }

    for (int j_i = 0; j_i < numBands_*numMobileTx_; ++j_i) {
        for (int n = 0; n < numMobileTxLocs_; ++n) {
            int row_index = j_i * numMobileTxLocs_ + n;
            int col_index = 0;
            for (int k = 0; k < numTimeSteps_; ++k) {
                for (int l = 0; l < numRxLocs_; ++l) {
                    TempReshapedConxMobTxRx[row_index][col_index] = swaped_temp_conx_mob_tx_rx[j_i][n][k][l];
                    ++col_index;
                }
            }
        }
    }

    //temp_reshaped_connectivity_fixed_tx = self.conx_fixed_tx_per_num_rx_locs.reshape(self.num_fixed_tx*self.num_bands, -1)
    //use ConxFixedTxPerNumRXLocs(numFixedTx_*numBands_*numTimeSteps_, numRxLocs)
    //TempReshapedConnectivityFixedTx(vector<vector<int>>(NUM_FIXED_TX*NUM_BANDS, vector<int>(NUM_TIME_STEPS*NUM_RX_LOCS, 0)))

    vector<int> flatData;
    for (const auto& row : ConxFixedTxPerNumRXLocs) {
        for (int val : row) {
            flatData.push_back(val);
        }
    }

    int index = 0;
    for (int i = 0; i < numFixedTx_*numBands_; ++i) {
        for (int j = 0; j < numTimeSteps_*numRxLocs_; ++j) {
            TempReshapedConnectivityFixedTx[i][j] = flatData[index++];
        }
    }

}


class OutputsMOBARP {
   private:
    int numFixedTx_;
    int numMobileTx_;
    int numBands_;
    int numTimeSteps_;
    int numRxLocs_;
    int numMobileTxLocs_;
    bool excludeCapFlag_;

   public:
    vector<impalib_type> ExtrinsicFixedX; 
    vector<impalib_type> ExtrinsicMobileX; 
    vector<impalib_type> ExtrinsicR; 
    vector<vector<impalib_type>> ExtrinsicZ;
    void extrinsic_update(vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &, vector<vector<impalib_type>>&, vector<vector<vector<impalib_type>>>&,
                                            vector<vector<impalib_type>>&, vector<vector<impalib_type>>& ,vector<vector<impalib_type>>&, vector<vector<vector<vector<impalib_type>>>>&,
                                            vector<vector<int>>&, vector<vector<int>>&, vector<vector<vector<vector<int>>>>&, vector<vector<int>>&);

    vector<vector<impalib_type>> transpose_reshape(vector<vector<vector<impalib_type>>>) const;

    OutputsMOBARP(const int, const int, const int, const int, const int, const int, const bool);
};

inline OutputsMOBARP::OutputsMOBARP(const int NUM_FIXED_TX, const int NUM_MOBILE_TX, const int NUM_BANDS, const int NUM_TIME_STEPS, const int NUM_RX_LOCS, const int NUM_MOBILE_TX_LOCS, const bool EXCLUDE_CAP_FLAG)
        :numFixedTx_(NUM_FIXED_TX),
        numMobileTx_(NUM_MOBILE_TX),
        numBands_(NUM_BANDS),
        numTimeSteps_(NUM_TIME_STEPS),
        numRxLocs_(NUM_RX_LOCS),
        numMobileTxLocs_(NUM_MOBILE_TX_LOCS),
        excludeCapFlag_(EXCLUDE_CAP_FLAG),
        ExtrinsicFixedX(NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS, zero_value),
        ExtrinsicMobileX(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, zero_value),
        ExtrinsicR(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, zero_value),
        ExtrinsicZ(NUM_BANDS*NUM_MOBILE_TX*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, zero_value)){};

inline vector<vector<impalib_type>> OutputsMOBARP::transpose_reshape(vector<vector<vector<impalib_type>>> rSetCoverIneqConst2FixedXEqConstM) const
{

    vector<vector<vector<impalib_type>>> transposed_matrix(numFixedTx_*numBands_, vector<vector<impalib_type>>(numTimeSteps_, vector<impalib_type>(numRxLocs_, 0)));

    vector<vector<impalib_type>> reshaped_matrix;

    for (int i = 0; i < numTimeSteps_; ++i) {
        for (int j = 0; j < numRxLocs_; ++j) {
            for (int k = 0; k < numFixedTx_*numBands_; ++k) {
                transposed_matrix[k][i][j] = rSetCoverIneqConst2FixedXEqConstM[i][j][k];
            }
        }
    }

    for (int k = 0; k < numFixedTx_*numBands_; ++k) {
        for (int i = 0; i < numTimeSteps_; ++i) {
            vector<impalib_type> row;
            row.reserve(numRxLocs_);
            for (int j = 0; j < numRxLocs_; ++j) {
                row.push_back(transposed_matrix[k][i][j]);
            }
            reshaped_matrix.push_back(row);
        }
    }
    
    return reshaped_matrix;
}

inline void OutputsMOBARP::extrinsic_update(vector<vector<vector<impalib_type>>>& rFixedCapacConst2FixedXEqConstM, vector<vector<vector<impalib_type>>>& rMobileCapacConst2MobileXEqConstM,
                                        vector<vector<impalib_type>>& rAuxiliaryConst2MobileXEqConstM, vector<vector<vector<impalib_type>>>& rSetCoverIneqConst2FixedXEqConstM,
                                        vector<vector<impalib_type>>& rAuxiliaryConst2REqConstM, vector<vector<impalib_type>>& rMobileLocEqConst2REqConstM,
                                        vector<vector<impalib_type>>& rAuxiliaryConst2ZEqConstM, vector<vector<vector<vector<impalib_type>>>>& rSetCoverIneqConst2ZEqConstM,
                                        vector<vector<int>>& rConxMobTxPerNumMobTxLocs, vector<vector<int>>& rConxFixedTxPerNumRXLocs,
                                        vector<vector<vector<vector<int>>>>& rConxMobTxRx, vector<vector<int>>& rConxMobTxR){

    vector<int> conx_rows_mobile;

    for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++) {
        int sum = accumulate(rConxMobTxPerNumMobTxLocs[i].begin(), rConxMobTxPerNumMobTxLocs[i].end(), 0);
        if (sum !=0){
            conx_rows_mobile.push_back(i);
        }
    }

    vector<impalib_type> flattened_mobile_capac_const_to_mobile_x_eq_const_m;

    for (auto& matrix : rMobileCapacConst2MobileXEqConstM) {
        for (auto& row : matrix) {
            flattened_mobile_capac_const_to_mobile_x_eq_const_m.insert(flattened_mobile_capac_const_to_mobile_x_eq_const_m.end(), row.begin(), row.end());
        }
    }
    

    // vector<impalib_type> temp_extrinsic_mobile_x(numMobileTx_*numBands_*numTimeSteps_, 0); //extrinsic mobile x

    for (int i=0; i<conx_rows_mobile.size(); i++) {
        impalib_type row_sum = 0;
        for (int j = 0; j < numMobileTxLocs_; j++) {
            if (rConxMobTxPerNumMobTxLocs[conx_rows_mobile[i]][j] == 1) {
                row_sum += rAuxiliaryConst2MobileXEqConstM[conx_rows_mobile[i]][j];
            }
        }
        if (!(excludeCapFlag_)){
            ExtrinsicMobileX[conx_rows_mobile[i]] = row_sum + flattened_mobile_capac_const_to_mobile_x_eq_const_m[conx_rows_mobile[i]];
        }
        
        else{
            ExtrinsicMobileX[conx_rows_mobile[i]] = row_sum;
        }
    }


    vector<int> conx_rows_fixed;

    for (int i=0; i<numFixedTx_*numBands_*numTimeSteps_; i++) {
        int sum = std::accumulate(rConxFixedTxPerNumRXLocs[i].begin(), rConxFixedTxPerNumRXLocs[i].end(), 0);
        if (sum !=0){
            conx_rows_fixed.push_back(i);
        }
    }

    // vector<impalib_type> temp_extrinsic_fixed_x(numFixedTx_*numBands_*numTimeSteps_, 0);
    vector<impalib_type> flattened_fixed_capac_const_to_fixed_x_eq_const_m;

    auto reshaped_set_cover_ineq_const_to_fixed_x_eq_const = transpose_reshape(rSetCoverIneqConst2FixedXEqConstM);

    for (auto& matrix : rFixedCapacConst2FixedXEqConstM) {
        for (auto& row : matrix) {
            flattened_fixed_capac_const_to_fixed_x_eq_const_m.insert(flattened_fixed_capac_const_to_fixed_x_eq_const_m.end(), row.begin(), row.end());
        }
    }

    for (int i=0; i<conx_rows_fixed.size(); i++) {
        impalib_type row_sum = 0;
        for (int j = 0; j < numRxLocs_; j++) {
            if (rConxFixedTxPerNumRXLocs[conx_rows_fixed[i]][j] == 1) {
                row_sum += reshaped_set_cover_ineq_const_to_fixed_x_eq_const[conx_rows_fixed[i]][j];
            }
        }

        if (!(excludeCapFlag_)){
            ExtrinsicFixedX[conx_rows_fixed[i]] = row_sum + flattened_fixed_capac_const_to_fixed_x_eq_const_m[conx_rows_fixed[i]];
        }
        
        else{
            ExtrinsicFixedX[conx_rows_fixed[i]] = row_sum;
        }
    }


    vector<vector<vector<impalib_type>>> reshaped_1(numMobileTx_, vector<vector<impalib_type>>(numBands_*numTimeSteps_, vector<impalib_type>(numMobileTxLocs_, 0)));

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
            for (size_t n = 0; n < numMobileTxLocs_; n++) {
                reshaped_1[i][j_k][n] = rAuxiliaryConst2REqConstM[i * numBands_ * numTimeSteps_ + j_k][n];
            }
        }
    }

    vector<vector<vector<impalib_type>>> reshaped_2(numMobileTx_, vector<vector<impalib_type>>(numMobileTxLocs_, vector<impalib_type>(numBands_*numTimeSteps_, 0)));


    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t n = 0; n < numMobileTxLocs_; n++) {
            for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
                reshaped_2[i][n][j_k] = reshaped_1[i][j_k][n];
            }
        }
    }

    vector<vector<impalib_type>> reshaped_auxiliary_const_to_r_eq_const_m(numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numBands_*numTimeSteps_, 0));

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t n = 0; n < numMobileTxLocs_; n++) {
            for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
                reshaped_auxiliary_const_to_r_eq_const_m[i * numMobileTxLocs_ + n][j_k] = reshaped_2[i][n][j_k];
            }
        }
    }


    // vector<impalib_type> temp_extrinsic_r(numMobileTx_*numMobileTxLocs_, 0);

    vector<int> sums_conx_per_row(numMobileTx_*numMobileTxLocs_, 0);
    vector<int> conx_rows;

    for (int i=0; i<sums_conx_per_row.size(); i++){
        sums_conx_per_row[i] = accumulate(rConxMobTxR[i].begin(), rConxMobTxR[i].end(), 0);
        if (sums_conx_per_row[i] !=0) {conx_rows.push_back(i);}
    }

    vector<impalib_type> flattened_mobile_loc_eq_const_to_r_eq_const_m;

    for (const auto& row : rMobileLocEqConst2REqConstM) {
        flattened_mobile_loc_eq_const_to_r_eq_const_m.insert(flattened_mobile_loc_eq_const_to_r_eq_const_m.end(), row.begin(), row.end());
    }

    for (int i=0; i<conx_rows.size(); i++){
        
        impalib_type sum_elements = 0;
        
        for (size_t l = 0; l < numBands_*numTimeSteps_; l++) {
            if (rConxMobTxR[conx_rows[i]][l] == 1) {
                sum_elements += reshaped_auxiliary_const_to_r_eq_const_m[conx_rows[i]][l];
            }
        }
        
        ExtrinsicR[conx_rows[i]] = sum_elements + flattened_mobile_loc_eq_const_to_r_eq_const_m[conx_rows[i]];
    }


    // vector<vector<impalib_type>> reshaped_set_cover_ineq_const_to_z_eq_const_m(numBands_*numMobileTx_*numTimeSteps_, vector<impalib_type>(numMobileTxLocs_, 0));

    // vector<vector<impalib_type>> temp_extrinsic_z(reshaped_set_cover_ineq_const_to_z_eq_const_m);

    for (int j_i=0; j_i< numBands_*numMobileTx_; j_i++){
        for (int k=0; k< numTimeSteps_; k++){
            for (int n=0; n< numMobileTxLocs_; n++){
                impalib_type temp_sum = 0;
                for (int l=0; l< numRxLocs_; l++){
                    if (rConxMobTxRx[k][l][n][j_i] == 1){
                        temp_sum += rSetCoverIneqConst2ZEqConstM[k][l][n][j_i];
                    }
                }
            // reshaped_set_cover_ineq_const_to_z_eq_const_m[j_i*numTimeSteps_ + k][n] =  temp_sum;
            ExtrinsicZ[j_i*numTimeSteps_ + k][n] = temp_sum + rAuxiliaryConst2ZEqConstM[j_i*numTimeSteps_ + k][n];
            }
        }
    }


    // fstream file_output_1("./ut_results/extrinsic_fixed_x_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<ExtrinsicFixedX.size(); i++){
    //         file_output_1.write((char*)(&ExtrinsicFixedX[i]), sizeof(ExtrinsicFixedX[i]));}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

    // fstream file_output_2("./ut_results/extrinsic_mobile_x_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_2.is_open()) {
    //     for (int i=0; i<ExtrinsicMobileX.size(); i++){
    //         file_output_2.write((char*)(&ExtrinsicMobileX[i]), sizeof(ExtrinsicMobileX[i]));}
    //         file_output_2.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

    // fstream file_output_3("./ut_results/extrinsic_r_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_3.is_open()) {
    //     for (int i=0; i<ExtrinsicR.size(); i++){
    //         file_output_3.write((char*)(&ExtrinsicR[i]), sizeof(ExtrinsicR[i]));}
    //         file_output_3.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

    // fstream file_output_4("./ut_results/extrinsic_z_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_4.is_open()) {
    //     for (int i=0; i<ExtrinsicZ.size(); i++){
    //     for (int j=0; j<ExtrinsicZ[0].size(); j++){
    //         file_output_4.write((char*)(&ExtrinsicZ[i][j]), sizeof(ExtrinsicZ[i][j]));}}
    //         file_output_4.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

}

