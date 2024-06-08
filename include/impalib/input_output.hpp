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
    int nTeams_;    ///< number of teams
    int nDept_;     ///< number of departments
    int nProj_;     ///< number of projects
    int nNonzero_;  ///< maximum # of non-zero weights over all departments

   public:
    vector<impalib_type> rewardsTeam_;             ///< rewards of team equality constraints
    vector<vector<impalib_type>> team2KnapsackM_;  ///< messages from teams to knapsack constraints
    vector<vector<int>> weights_;                  ///< weights of teams per each department
    vector<vector<impalib_type>> rewardsProj_;     ///< reward of project equality constraint
    vector<int> capacities_;                       ///< vector of capacities of departments
    vector<vector<int>> nonzero_;                  ///< indices of non-zero weights per each department

    void process_inputs(const impalib_type *, impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *);  ///< process input of graphical model

    InputsKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS, int MAX_SIZE_NON_ZERO_WEIGHTS);  ///< constructor
};

/**
 * Represents a class for outputs of Knapsack-MWMW problem
 */
class OutputsKcMwm {
   private:
    int nTeams_;  ///< number of teams
    int nDept_;   ///< number of departments
    int nProj_;   ///< number of projects

   public:
    vector<impalib_type> extrinsicOut_;  ///< extrinsic output of team equality constraints
    vector<impalib_type> intrinsicOut_;  ///< intrinsic outputs of project equality constraint
    void intrinsic_out_mwm_update(const vector<vector<impalib_type>> &, const vector<vector<impalib_type>> &,
                                  const vector<vector<impalib_type>> &);                        ///< calculate intrinsic outputs of project equality constraints
    void extrinsic_output_team_update(vector<vector<impalib_type>> &, vector<impalib_type> &);  ///< calculate extrinsic output of team equality constraints
    OutputsKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);                               ///< constructor
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
    : nDept_(N_DEPARTMENTS), nTeams_(N_TEAMS), nProj_(N_PROJECTS), nNonzero_(MAX_SIZE_NON_ZERO_WEIGHTS) {
    weights_.reserve(nDept_);
    nonzero_.reserve(nDept_);
    rewardsProj_.reserve(nProj_);
    team2KnapsackM_.reserve(nDept_);

    for (int department = 0; department < nDept_; department++) {
        team2KnapsackM_.push_back(vector<impalib_type>(nTeams_, zero_value));
        weights_.push_back(vector<int>(nTeams_, 0));
    }
    for (int project = 0; project < nProj_; project++) {
        rewardsProj_.push_back(vector<impalib_type>(nTeams_, zero_value));
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

inline void InputsKcMwm::process_inputs(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                                        const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
    copy(pREWARD_TEAM_PY, pREWARD_TEAM_PY + nTeams_, back_inserter(rewardsTeam_));
    copy(pMAX_STATE_PY, pMAX_STATE_PY + nDept_, back_inserter(capacities_));

    for (int department = 0; department < nDept_; department++) {
        nonzero_.push_back(vector<int>(pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department], 0));
        copy(pTransition_model_py + nTeams_ * department, pTransition_model_py + nTeams_ * (department + 1), team2KnapsackM_[department].begin());
        copy(pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + nTeams_ * department, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + nTeams_ * (department + 1), weights_[department].begin());
        copy(p_NON_ZERO_WEIGHT_INDICES_PY + nNonzero_ * department, p_NON_ZERO_WEIGHT_INDICES_PY + pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department] + nNonzero_ * department,
             nonzero_[department].begin());
    }
    for (int project = 0; project < nProj_; project++) {
        // Copy reward values for each project-team combination
        copy(pREWARD_PROJECT_PY + nTeams_ * project, pREWARD_PROJECT_PY + nTeams_ * (project + 1), rewardsProj_[project].begin());
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

inline OutputsKcMwm::OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : nDept_(N_DEPARTMENTS), nTeams_(N_TEAMS), nProj_(N_PROJECTS) {
    extrinsicOut_.reserve(nTeams_);
    extrinsicOut_.resize(nTeams_);
    fill(extrinsicOut_.begin(), extrinsicOut_.begin() + nTeams_, zero_value);

    intrinsicOut_.reserve(nProj_ * nTeams_);
    intrinsicOut_.resize(nProj_ * nTeams_);
    fill(intrinsicOut_.begin(), intrinsicOut_.begin() + nProj_ * nTeams_, zero_value);
};

/**
 * Calculate instrinsic messages of project equality constraint for Knapsack-MWM problem
 *
 * @param[in] oric2EqM: messages for ORIC to equality constraint of the projects
 * @param[in] project2EqM: messages from projects to equality constraints of the projects
 * @param[in] rewards: rewards of project-team combinations
 *
 */

inline void OutputsKcMwm::intrinsic_out_mwm_update(const vector<vector<impalib_type>> &oric2EqM, const vector<vector<impalib_type>> &project2EqM, const vector<vector<impalib_type>> &rewards) {
    for (int project = 0; project < rewards.size(); project++) {
        for (int team = 0; team < rewards[project].size(); team++) {
            intrinsicOut_[project + team + project * (nTeams_ - 1)] = oric2EqM[project][team] + project2EqM[project][team] + rewards[project][team];
        }
    }
}

/**
 * Calculate extrinsic messages of team equality constraint for Knapsack-MWM problem
 *
 * @param[in] extrinsicOut: messages from knapsack constraints to team equality constraints
 * @param[in] oric2TeamM: messages from ORIC to team equality constraints
 *
 */

inline void OutputsKcMwm::extrinsic_output_team_update(vector<vector<impalib_type>> &extrinsicOut, vector<impalib_type> &oric2TeamM) {
    copy(oric2TeamM.begin(), oric2TeamM.end(), extrinsicOut_.begin());

    for (int department = 0; department < extrinsicOut.size(); department++) {
        transform(extrinsicOut[department].begin(), extrinsicOut[department].end(), extrinsicOut_.begin(), extrinsicOut_.begin(), std::plus<impalib_type>());
    }
}

/**
 * Represents a class for inputs of TSP
 */
class InputsTsp {
   private:
    int nNodes_;             ///< number of nodes of TSP
    int nEdges_;             ///< number of edges of TSP
    int nNodesPerEdge_ = 2;  ///< number of nodes per edge

   public:
    vector<vector<int>> connections_;             ///< constituent nodes per each edge
    vector<impalib_type> costs_;                  ///< cost for each edge equality constraint
    vector<vector<impalib_type>> costMat_;        ///< cost matrix (n_nodesxn_nodes)
    vector<vector<impalib_type>> costMatUpdate_;  ///< cost matrix represented as num_edges x num_nodes to facilitate message updates
    vector<vector<impalib_type>> eq2DegM_;        ///< messages from edge equality constraint to degree constraint

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

inline InputsTsp::InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : nNodes_(NUM_NODES), nEdges_(NUM_EDGE_VARIABLES){};

/**
 * Represents a class for outputs of TSP
 */
class OutputsTsp {
   private:
    int nNodes_;  ///< number of nodes
    int nEdges_;  ///< number of edges

   public:
    vector<impalib_type> extrinsicOut_;                                                  ///< extrinsic output of edge equality constraint
    vector<impalib_type> intrinsicOut_;                                                  ///< intrinsic output of edge equality constraint
    void extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>> &);  ///< calculate extrinsic output of edge equality constraint for a relaxed TSP
    void extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>> &,
                                                         vector<vector<impalib_type>> &);  ///< calculate extrinsic output of edge equality constraint for augmented TSP
    void intrinsic_output_edge_ec_update(vector<impalib_type> &);                          ///< calculate intrinsic output of edge equality constraint for augmented TSP
    OutputsTsp(int NUM_NODES, int NUM_EDGE_VARIABLES);                                     ///< constructor
};

/**
 * Construct Output object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 *
 */

inline OutputsTsp::OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : nNodes_(NUM_NODES), nEdges_(NUM_EDGE_VARIABLES) {
    extrinsicOut_.reserve(nEdges_);
    extrinsicOut_.resize(nEdges_);
    intrinsicOut_.reserve(nEdges_);
    intrinsicOut_.resize(nEdges_);
    fill(extrinsicOut_.begin(), extrinsicOut_.begin() + nEdges_, zero_value);
    fill(intrinsicOut_.begin(), intrinsicOut_.begin() + nEdges_, zero_value);
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

inline void InputsTsp::process_inputs(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY,
                                      impalib_type *pEdge_ec_to_degree_constraint_m_py, const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
    copy(pCOST_EDGE_VARIABLE_PY, pCOST_EDGE_VARIABLE_PY + nEdges_, back_inserter(costs_));

    for (int edge = 0; edge < nEdges_; edge++) {
        connections_.push_back(vector<int>(nNodesPerEdge_, 0));
        copy(pEDGE_CONNECTIONS_PY + nNodesPerEdge_ * edge, pEDGE_CONNECTIONS_PY + nNodesPerEdge_ * (edge + 1), connections_[edge].begin());

        eq2DegM_.push_back(vector<impalib_type>(nNodes_, zero_value));
        copy(pEdge_ec_to_degree_constraint_m_py + nNodes_ * edge, pEdge_ec_to_degree_constraint_m_py + nNodes_ * (edge + 1), eq2DegM_[edge].begin());

        costMatUpdate_.push_back(vector<impalib_type>(nNodes_, zero_value));
        copy(pEDGE_DEGREE_CONSTRAINT_COST_PY + nNodes_ * edge, pEDGE_DEGREE_CONSTRAINT_COST_PY + nNodes_ * (edge + 1), costMatUpdate_[edge].begin());
    }

    for (int node = 0; node < nNodes_; node++) {
        costMat_.push_back(vector<impalib_type>(nNodes_, zero_value));
        copy(pCOST_MATRIX_PY + nNodes_ * node, pCOST_MATRIX_PY + nNodes_ * (node + 1), costMat_[node].begin());
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for relaxed TSP
 *
 * @param[in] deg2EqM: messages from degree constraints to edge equality constraints
 *
 */

inline void OutputsTsp::extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>> &deg2EqM) {
    for (int edge = 0; edge < deg2EqM.size(); edge++) {
        extrinsicOut_[edge] = accumulate(deg2EqM[edge].begin(), deg2EqM[edge].end(), zero_value);
    }
}

/**
 * Calculate output extrinsic messages of edge equality constraint for augmented TSP
 *
 * @param[in] deg2EqM: messages from degree constraints to edge equality constraints
 * @param[in] subtour2EdgeM: messages from subtour elimination constraints to edge equality constraints
 *
 */

inline void OutputsTsp::extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>> &deg2EqM, vector<vector<impalib_type>> &subtour2EdgeM) {
    for (int edge = 0; edge < deg2EqM.size(); edge++) {
        extrinsicOut_[edge] = accumulate(deg2EqM[edge].begin(), deg2EqM[edge].end(), zero_value);
    }

    for (int subtour = 0; subtour < subtour2EdgeM.size(); subtour++) {
        transform(subtour2EdgeM[subtour].begin(), subtour2EdgeM[subtour].end(), extrinsicOut_.begin(), extrinsicOut_.begin(), plus<impalib_type>());
    }
}

/**
 * Calculate output intrinsic messages of edge equality constraint for TSP
 *
 * @param[in] costs: cost of each edge equality constraint
 *
 */

inline void OutputsTsp::intrinsic_output_edge_ec_update(vector<impalib_type> &costs) {
    transform(extrinsicOut_.begin(), extrinsicOut_.end(), costs.begin(), intrinsicOut_.begin(), plus<impalib_type>());
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
    int nVars_;         ///< total number of variables
    int nConstraints_;  ///< number of constraints
    int k_;             ///< number of variables per constraint

   public:
    vector<impalib_type> extrinsicOut_;                                   ///< extrinsic messages of variables equality constraints
    OutputsKsat(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE);  ///< constructor
    void update_extrinsic(const vector<vector<impalib_type>> &);          ///< calculate extrinsic messages of variables equality constraints
};

/**
 * Construct Output object for k-sat problem
 *
 * @param[in] NUM_VARIABLES: total number of variables
 * @param[in] NUM_CONSTRAINTS: number of constraints
 * @param[in] K_VARIABLE: ///< number of variables per constraint
 *
 */

inline OutputsKsat::OutputsKsat(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE) : nVars_(NUM_VARIABLES), nConstraints_(NUM_CONSTRAINTS), k_(K_VARIABLE), extrinsicOut_(nVars_, zero_value){};

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
 * @param[in] ksat2EqM: messages from k-sat constraint to variable equality constraint
 *
 */

inline void OutputsKsat::update_extrinsic(const vector<vector<impalib_type>> &ksat2EqM) {
    // Sum all messages coming into equality constraint except the
    // incoming message on the edge of interest
    for (int i = 0; i < nVars_; i++) {
        for (int j = 0; j < nConstraints_; j++) {
            extrinsicOut_[i] += ksat2EqM[j][i];
        }
    }
}