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
    vector<impalib_type> RewardTeam;               ///< rewards of team equality constraints
    vector<vector<impalib_type>> M_team2knapsack;  ///< messages from teams to knapsack constraints
    vector<vector<int>> weights;                   ///< weights of teams per each department
    vector<vector<impalib_type>> RewardProject;    ///< reward of project equality constraint
    vector<int> capacity;                          ///< vector of capacities of departments
    vector<vector<int>> NonZeroWeightIndices;      ///< indices of non-zero weights per each department

    void process_inputs(const impalib_type *, const impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *);  ///< process input of graphical model

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
    vector<impalib_type> extrinsic;  ///< extrinsic_ output of team equality constraints
    vector<impalib_type> intrinsic;  ///< intrinsic outputs of project equality constraint
    void update_intrinsic(const vector<vector<impalib_type>> &oric2eq, const vector<vector<impalib_type>> &project2eq,
                          const vector<vector<impalib_type>> &reward);                                                ///< calculate intrinsic outputs of project equality constraints
    void update_extrinsic(const vector<vector<impalib_type>> &extrinsic_out, const vector<impalib_type> &oric2team);  ///< calculate extrinsic_ output of team equality constraints
    OutputsKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);                                                     ///< constructor
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
    weights.reserve(numDepartments_);
    NonZeroWeightIndices.reserve(numDepartments_);
    RewardProject.reserve(numProjects_);
    M_team2knapsack.reserve(numDepartments_);

    // Initialize M_team2knapsack and weights vectors
    for (int i = 0; i < numDepartments_; i++) {
        M_team2knapsack.push_back(vector<impalib_type>(numTeams_, zero_value));
        weights.push_back(vector<int>(numTeams_, 0));
    }

    // Initialize RewardProject vector
    for (int i = 0; i < numProjects_; i++) {
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
inline void InputsKcMwm::process_inputs(const impalib_type *pREWARD_TEAM_PY, const impalib_type *pTransition_model_py, const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                                 const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
    // Copy reward values for each team and maximum states per department
    copy(pREWARD_TEAM_PY, pREWARD_TEAM_PY + numTeams_, back_inserter(RewardTeam));
    copy(pMAX_STATE_PY, pMAX_STATE_PY + numDepartments_, back_inserter(capacity));

    for (int i = 0; i < numDepartments_; i++) {
        // Populate NonZeroWeightIndices for the department
        NonZeroWeightIndices.push_back(vector<int>(pNON_ZERO_WEIGHT_INDICES_SIZES_PY[i], 0));
        copy(pTransition_model_py + numTeams_ * i, pTransition_model_py + numTeams_ * (i + 1), M_team2knapsack[i].begin());
        copy(pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * i, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_ * (i + 1), weights[i].begin());
        copy(p_NON_ZERO_WEIGHT_INDICES_PY + maxSizeNonzeroWeights_ * i, p_NON_ZERO_WEIGHT_INDICES_PY + pNON_ZERO_WEIGHT_INDICES_SIZES_PY[i] + maxSizeNonzeroWeights_ * i,
             NonZeroWeightIndices[i].begin());
    }
    for (int i = 0; i < numProjects_; i++) {
        // Copy reward values for each project-team combination
        copy(pREWARD_PROJECT_PY + numTeams_ * i, pREWARD_PROJECT_PY + numTeams_ * (i + 1), RewardProject[i].begin());
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
inline OutputsKcMwm::OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numDepartments_(N_DEPARTMENTS), numTeams_(N_TEAMS), numProjects_(N_PROJECTS), extrinsic(numTeams_, 0), intrinsic(numProjects_ * numTeams_, 0){};

/**
 * Calculate instrinsic messages of project equality constraint for Knapsack-MWM problem
 *
 * @param[in] oric2eq: messages for ORIC to equality constraint of the projects
 * @param[in] project2eq: messages from projects to equality constraints of the projects
 * @param[in] reward: rewards of project-team combinations
 *
 */
inline void OutputsKcMwm::update_intrinsic(const vector<vector<impalib_type>> &oric2eq, const vector<vector<impalib_type>> &project2eq, const vector<vector<impalib_type>> &reward) {
    for (int i = 0; i < reward.size(); i++)  // over projects
    {
        for (int j = 0; j < reward[i].size(); j++)  // over teams
        {
            // Calculate the intrinsic output for the project-team combination
            intrinsic[i + j + i * (numTeams_ - 1)] = oric2eq[i][j] + project2eq[i][j] + reward[i][j];
        }
    }
}

/**
 * Calculate extrinsic_ messages of team equality constraint for Knapsack-MWM problem
 *
 * @param[in] extrinsic_out: messages from knapsack constraints to team equality constraints
 * @param[in] oric2team: messages from ORIC to team equality constraints
 *
 */
inline void OutputsKcMwm::update_extrinsic(const vector<vector<impalib_type>> &extrinsic_out, const vector<impalib_type> &oric2team) {
    extrinsic = oric2team;

    for (int i = 0; i < extrinsic_out.size(); i++) {
        transform(extrinsic_out[i].begin(), extrinsic_out[i].end(), extrinsic.begin(), extrinsic.begin(), std::plus<impalib_type>());
    }
}

/**
 * Represents a class for inputs of TSP
 */
class InputsTsp {
   private:
    int n_Nodes;             ///< number of nodes of TSP
    int n_Edges;             ///< number of edges of TSP
    int n_NodesPerEdge = 2;  ///< number of nodes per edge

   public:
    vector<vector<int>> edges;                   ///< constituent nodes per each edge
    vector<impalib_type> cost_edge;              ///< cost for each edge equality constraint
    vector<vector<impalib_type>> cost;           ///< cost matrix (n_nodesxn_nodes)
    vector<vector<impalib_type>> cost_edge_mat;  ///< cost matrix represented as num_edges x num_nodes to facilitate message updates
    vector<vector<impalib_type>> M_edge2degree;  ///< messages from edge equality constraint to degree constraint

    void process_inputs(const int *, const impalib_type *, const impalib_type *, const impalib_type *, const impalib_type *);  ///< process inputs of TSP graphical model

    InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES);  ///< constructor
};

/**
 * Construct Input object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 *
 */
inline InputsTsp::InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : n_Nodes(NUM_NODES), n_Edges(NUM_EDGE_VARIABLES){};

/**
 * Represents a class for outputs of TSP
 */
class OutputsTsp {
   private:
    int numNodes_;          ///< number of nodes
    int numEdgeVariables_;  ///< number of edges

   public:
    vector<impalib_type> extrinsic;                                       ///< extrinsic_ output of edge equality constraint
    vector<impalib_type> intrinsic;                                       ///< intrinsic output of edge equality constraint
    void update_extrinsic_relaxed(vector<vector<impalib_type>> &deg2eq);  ///< calculate extrinsic_ output of edge equality constraint for a relaxed TSP
    void update_extrinsic_augmented(vector<vector<impalib_type>> &deg2eq,
                                    vector<vector<impalib_type>> &subtour2edge);  ///< calculate extrinsic_ output of edge equality constraint for augmented TSP
    void update_intrinsic(vector<impalib_type> &costs);                           ///< calculate intrinsic output of edge equality constraint for augmented TSP
    OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES);                ///< constructor
};

/**
 * Construct Output object for the TSP
 *
 * @param[in] NUM_NODES: number of nodes
 * @param[in] NUM_EDGE_VARIABLES: number of edge variables (edge connections)
 *
 */
inline OutputsTsp::OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES) : numNodes_(NUM_NODES), numEdgeVariables_(NUM_EDGE_VARIABLES) {
    // Reserve and initialize memory for extrinsic_ and intrinsic edge equality constraints outputs
    extrinsic.reserve(numEdgeVariables_);
    extrinsic.resize(numEdgeVariables_);
    intrinsic.reserve(numEdgeVariables_);
    intrinsic.resize(numEdgeVariables_);
    fill(extrinsic.begin(), extrinsic.begin() + numEdgeVariables_, zero_value);
    fill(intrinsic.begin(), intrinsic.begin() + numEdgeVariables_, zero_value);
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
inline void InputsTsp::process_inputs(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY, const impalib_type *pEdge_ec_to_degree_constraint_m_py,
                               const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
    // Copy cost edge variable data
    copy(pCOST_EDGE_VARIABLE_PY, pCOST_EDGE_VARIABLE_PY + n_Edges, back_inserter(cost_edge));

    // Process input variables
    for (int i = 0; i < n_Edges; i++) {
        // Copy edge connections
        edges.push_back(vector<int>(n_NodesPerEdge, 0));
        copy(pEDGE_CONNECTIONS_PY + n_NodesPerEdge * i, pEDGE_CONNECTIONS_PY + n_NodesPerEdge * (i + 1), edges[i].begin());
        // Copy edge ec to degree constraint messages
        M_edge2degree.push_back(vector<impalib_type>(n_Nodes, zero_value));
        copy(pEdge_ec_to_degree_constraint_m_py + n_Nodes * i, pEdge_ec_to_degree_constraint_m_py + n_Nodes * (i + 1), M_edge2degree[i].begin());
        // Copy edge degree constraint cost
        cost_edge_mat.push_back(vector<impalib_type>(n_Nodes, zero_value));
        copy(pEDGE_DEGREE_CONSTRAINT_COST_PY + n_Nodes * i, pEDGE_DEGREE_CONSTRAINT_COST_PY + n_Nodes * (i + 1), cost_edge_mat[i].begin());
    }

    // Process cost matrix
    for (int i = 0; i < n_Nodes; i++) {
        cost.push_back(vector<impalib_type>(n_Nodes, zero_value));
        copy(pCOST_MATRIX_PY + n_Nodes * i, pCOST_MATRIX_PY + n_Nodes * (i + 1), cost[i].begin());
    }
}

/**
 * Calculate output extrinsic_ messages of edge equality constraint for relaxed TSP
 *
 * @param[in] deg2eq: messages from degree constraints to edge equality constraints
 *
 */
inline void OutputsTsp::update_extrinsic_relaxed(vector<vector<impalib_type>> &deg2eq) {
    for (int i = 0; i < deg2eq.size(); i++) {
        extrinsic[i] = accumulate(deg2eq[i].begin(), deg2eq[i].end(), zero_value);
    }
}

/**
 * Calculate output extrinsic_ messages of edge equality constraint for augmented TSP
 *
 * @param[in] deg2eq: messages from degree constraints to edge equality constraints
 * @param[in] subtour2edge: messages from subtour elimination constraints to edge equality constraints
 *
 */
inline void OutputsTsp::update_extrinsic_augmented(vector<vector<impalib_type>> &deg2eq, vector<vector<impalib_type>> &subtour2edge) {
    // Update extrinsic_ output for edge equality constraints based on messages from degree constraints to equality constraints
    for (int i = 0; i < deg2eq.size(); i++) {
        extrinsic[i] = accumulate(deg2eq[i].begin(), deg2eq[i].end(), zero_value);
    }

    // Update extrinsic_ output for edge equality constraints based on messages from subtour elimination constraints
    for (int i = 0; i < subtour2edge.size(); i++) {
        transform(subtour2edge[i].begin(), subtour2edge[i].end(), extrinsic.begin(), extrinsic.begin(), plus<impalib_type>());
    }
}

/**
 * Calculate output intrinsic messages of edge equality constraint for TSP
 *
 * @param[in] costs: cost of each edge equality constraint
 *
 */
inline void OutputsTsp::update_intrinsic(vector<impalib_type> &costs) {
    // Calculate intrinsic output for edge equality constraints by adding extrinsic_ output for edge equality constraints to cost edge variable
    transform(extrinsic.begin(), extrinsic.end(), costs.begin(), intrinsic.begin(), plus<impalib_type>());
}