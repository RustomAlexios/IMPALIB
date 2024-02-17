// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a graphical model for the Knapsack-assignment problem for
 * computing the optimum configuration of assigning teams from departments to
 * projects.
 */
class GraphicalModelKcMwm {
   private:
    int numProjects_;  ///< number of projects
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int maxSizeNonZeroWeights_; ///< maximum size of non-zero weights per department
    int numIterations_; ///< number of iterations of IMPA
    bool filteringFlag_; ///< filtering flag of knapsack constraints
    impalib_type alpha_; ///< filtering parameter
    vector<vector<impalib_type>> extrinsicOutputDepartmentDummy_; ///< messages from departments to teams before filtering
    vector<vector<impalib_type>> extrinsicOutputDepartment_; ///< messages from departments to teams after filtering
    vector<impalib_type> oric2PackageM_; ///< messages from ORIC to packages
    vector<vector<impalib_type>> eqConstraint2OricM_; ///< messages from team equality constraint to ORIC
    vector<vector<impalib_type>> oric2EqConstraintM_; ///< messages from ORIC to team equality constraint
    vector<vector<impalib_type>> eqConstraint2ProjectM_; ///< messages from project equality constraint to project inequality constraint
    vector<vector<impalib_type>> project2EqConstraintM_; ///< messages from project inequality constraint to project equality constraint
    vector<impalib_type> team2OricM_; ///< messages from team equality constraint to ORIC
    Knapsack modelKnapsacks_; ///< Knapsack object
    InequalityConstraint projectIneqConstraint_; ///< Project Inequality constraint object
    EqualityConstraintKcMwm modelEqConstraint_; ///< Equality constraint object
    OrInequalityConstraint modelOric_; ///< ORIC object

   public:
    OutputsKcMwm outputs; ///< Graphical model outputs objects
    InputsKcMwm modelInputs_; ///< Graphical model inputs objects
    void initialize(const impalib_type *, impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *); ///< initialize graphical model
    void iterate(const int *); ///< iterate over graphical model
    GraphicalModelKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG, const impalib_type ALPHA); ///< constructor
};

/**
 * Construct Graphical Model for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[in] MAX_SIZE_NON_ZERO_WEIGHTS: maximum number of connections between
 * teams and departments
 * @param[in] N_ITERATIONS: number of iterations for running IMPA
 * @param[in] FILT_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[out] modelInputs_: object of type InputsKcMwm: inputs of graphical
 * model
 * @param[out] outputs: object of type OutputsKcMwm: outputs of graphical model
 * @param[out] projectIneqConstraint_: object of type InequalityConstraint:
 * object associated with inequality constraint operations
 * @param[out] modelEqConstraint_: object of type EqualityConstraintKcMwm:
 * object associated with equality constraint operations
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numDepartments_: N_DEPARTMENTS
 * @param[out] maxSizeNonZeroWeights_: MAX_SIZE_NON_ZERO_WEIGHTS
 * @param[out] numIterations_: N_ITERATIONS
 * @param[out] filteringFlag_: FILT_FLAG
 * @param[out] alpha_: ALPHA
 *
 */

GraphicalModelKcMwm::GraphicalModelKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG,
                                         const impalib_type ALPHA)
    : modelInputs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS),
      outputs(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      modelKnapsacks_(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA),
      projectIneqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      modelEqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      modelOric_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      numProjects_(N_PROJECTS),
      numTeams_(N_TEAMS),
      numDepartments_(N_DEPARTMENTS),
      maxSizeNonZeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS),
      numIterations_(N_ITERATIONS),
      filteringFlag_(FILT_FLAG),
      alpha_(ALPHA),
      oric2PackageM_(numTeams_, 0),
      team2OricM_(numTeams_, 0),
      eqConstraint2OricM_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      oric2EqConstraintM_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      eqConstraint2ProjectM_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      project2EqConstraintM_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      extrinsicOutputDepartmentDummy_(numDepartments_, vector<impalib_type>(numTeams_, 0)),
      extrinsicOutputDepartment_(numDepartments_, vector<impalib_type>(numTeams_, 0))
{
};

/**
 * Initialize Graphical Model inputs for the Knapsack-MWM problem. The inputs
 * are processed by reading inputs from python
 *
 * @param[in] pREWARD_TEAM_PY: rewards of teams
 * @param[in] pTransition_model_py: teams to knapsack messages
 * @param[in] pITEMS_WEIGHTS_PER_DEPARTMENT_PY: teams weights per department:
 * weights of each team associated with all departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: sizes of connections between
 * teams and departments
 * @param[in] p_NON_ZERO_WEIGHT_INDICES_PY: non-zero connections between teams
 * and departments
 * @param[in] pREWARD_PROJECT_PY: rewards for project-team combination
 * @param[in] pMAX_STATE_PY: contains maximum capacity of departments
 *
 */

void GraphicalModelKcMwm::initialize(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pITEMS_WEIGHTS_PER_DEPARTMENT_PY, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                                     const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
    /// calls a method process_inputs() on an object modelInputs_, passing
    /// several pointers to Python objects as arguments
    modelInputs_.process_inputs(pREWARD_TEAM_PY, pTransition_model_py, pITEMS_WEIGHTS_PER_DEPARTMENT_PY, pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY,
                                pMAX_STATE_PY);
}

/**
 * Run IMPA on Knapsack-MWM graphical model. This will propagate messages for a
 * certain number of iterations across the whole graphical model
 *
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: indices of non-zero connections
 * between teams and departments
 *
 */

void GraphicalModelKcMwm::iterate(const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY) {
    for (int i = 0; i < numIterations_; i++) {
        for (int j = 0; j < numDepartments_; j++) {
            int max_state_department = modelInputs_.MaxState[j];
            auto& idx_nonzero_dept = modelInputs_.NonZeroWeightIndices[j];
            auto& team_weights = modelInputs_.TeamsWeightsPerDepartment[j];

            // Initialize forward and backward message vectors
            vector<vector<impalib_type>> stage_forward_messages(numTeams_ + 1, vector<impalib_type>(max_state_department + 1, zero_value));
            vector<vector<impalib_type>> stage_backward_messages(numTeams_ + 1, vector<impalib_type>(max_state_department + 1, zero_value));

            modelKnapsacks_.forward(j, stage_forward_messages, max_state_department, idx_nonzero_dept, pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                                    team_weights, modelInputs_.Team2KnapsackM);
            modelKnapsacks_.backward(j, stage_backward_messages, max_state_department, idx_nonzero_dept, pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                                     team_weights, modelInputs_.Team2KnapsackM);

            modelKnapsacks_.extrinsic_output_department_lhs(modelInputs_.TeamsWeightsPerDepartment, stage_forward_messages, modelInputs_.Team2KnapsackM, j, stage_backward_messages,
                                                            max_state_department, extrinsicOutputDepartmentDummy_);
            modelKnapsacks_.process_extrinsic_output_department(j, i, extrinsicOutputDepartmentDummy_, extrinsicOutputDepartment_);
        }

        modelEqConstraint_.team_eq_constraint_to_oric_update(extrinsicOutputDepartment_, team2OricM_, modelInputs_.RewardTeam);
        modelOric_.oric_to_project_eq_constraint_update(eqConstraint2OricM_, team2OricM_, oric2EqConstraintM_, eqConstraint2ProjectM_, modelInputs_.RewardProject);
        projectIneqConstraint_.project_inequality_constraint_update(eqConstraint2ProjectM_, project2EqConstraintM_);
        modelEqConstraint_.project_eq_constraint_to_oric_update(project2EqConstraintM_, eqConstraint2OricM_, modelInputs_.RewardProject);

        modelOric_.oric_to_team_update(eqConstraint2OricM_, oric2PackageM_);
        outputs.update_intrinsic(oric2EqConstraintM_, project2EqConstraintM_, modelInputs_.RewardProject);
        modelKnapsacks_.team_to_knapsack_update(modelInputs_.NonZeroWeightIndices, modelInputs_.Team2KnapsackM, modelInputs_.RewardTeam, extrinsicOutputDepartment_, oric2PackageM_);
    }
    outputs.update_extrinsic(extrinsicOutputDepartment_, oric2PackageM_);
}

/**
 * Represents a graphical model for the TSP for computing the shortest route
 * between nodes (while visiting each node exactly once).
 */
class GraphicalModelTsp {
   private:
    int numNodes_; ///< number of nodes of TSP
    int numIterations_; ///< number of iterations of IMPA
    int numEdgeVariables_; ///< number of edges in TSP
    bool resetFlag_; ///< flag for resetting messages after each augmentation step
    bool filteringFlag_; ///< flag for filtering messages from degree/subtour constraints to edge equality constraints
    bool augmentationFlag_; ///< whether to activate augmentation
    vector<int> hard_decision; ///< hard decision vector of IMPA solution
    impalib_type alpha_; ///< filtering parameter
    impalib_type threshold_; ///< threshold on hard decision
    EqualityConstraintTsp modelEqConstraint_; ///< Equality Constraint object for TSP
    DegreeConstraint modelDegreeConstraint_; ///< Degree Constraint object for TSP
    SubtourEliminationConstraint modelSubtourEliminationConstraint_; ///< Subtour constraint for TSP
    vector<vector<impalib_type>> DegreeConstraint2EqConstraintDummyM_; ///< messages from degree constraint to team equality constraint before filtering
    vector<vector<impalib_type>> DegreeConstraint2EqConstraintM_; ///< messages from degree constraint to team equality constraint afters filtering
    vector<vector<int>> selected_edges_old_; ///< old list of selected edges (used for failure investigation)
    vector<vector<int>> selected_edges_old_old_; ///< another old list of selected edges (used for failure investigation)
    int maxCount_; ///< maximum count of failures
    bool tourImpaFlag_ = false; ///< set flag for detected tour to be false (will be updated later)
    vector<vector<int>> delta_S_indices_list; ///< will contain indices of edges that contribute to each subtour constraint
    vector<vector<impalib_type>> subtourConstraints2EdgeEcDummyM_; ///< messages from subtour constraints to edge equality constraint before filtering
    vector<vector<impalib_type>> edgeEc2SubtourConstraintsM_; ///< messages from edge equality constraint to subtour constraints
    void iterate_augmented_graph(); ///< function of IMPA on augmented graph
    void subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &, vector<vector<int>> &); ///< analysis of subtour constraints
    void hard_decision_analysis(vector<vector<int>> &); ///< function for hard decision solution on IMPA solution
    bool isSubsequence(const vector<int> &, const vector<int> &, int); ///< function for post-processing loops
    vector<vector<int>> get_closed_loops(unordered_map<int, vector<int>> &, vector<vector<int>> &); ///< function for getting loops
    vector<int> find_closed_loop(unordered_map<int, vector<int>> &, int, int, unordered_set<int>, vector<int>, vector<int> &); ///< function for finding loops
    InputsTsp modelInputs_; ///< Graphical Model Input object
    vector<vector<int>> selectedEdges_; ///< activated edges of IMPA
    int numAugmentations_ = 0; ///< number of performed augmentations in IMPA
    int noConsClosedLoopsCount_ = 0; ///< count for failure case (no consecutive loop detection and no tour)
    int solOscCount_ = 0; ///< count for failure case (oscillation in the solution)
    bool noConsClosedLoopsCountExcFlag_ = false; ///< flag for failure case (no consecutive loop detection and no tour)
    int noImprovSolCount_ = 0; ///< count for failure case (no solution improvement)
    bool noImprovSolCountExcFlag_ = false; ///< flag for failure case (no solution improvement)
    bool solOscCountExcFlag_ = false; ///< flag for failure case (oscillation in the solution)
    vector<int> tourImpa_; ///< list of nodes of tour (if detected)
    vector<vector<int>> subtourPaths_; ///< list of detected subtours (if detected)
    vector<int> closedPathsSize_; ///< list of sizes of loops
    impalib_type costImpa_ = zero_value; ///< cost of IMPA solution
    vector<vector<impalib_type>> subtourConstraints2EdgeEcM_; ///< messages from subtour constraints to edge equality constraint

   public:
    bool subtourConstraintsSatisfiedFlag = false; ///< initially set this flag for satisfied subtour constraints to false
    OutputsTsp outputs; ///< TSP graphical model outputs object
    void initialize(const int *, const impalib_type *, const impalib_type *, impalib_type *, const impalib_type *); ///< initialize graphical model
    void iterate_relaxed_graph(); ///< iterate over relaxed graphical model
    void perform_augmentation(const int); ///< perform augmentation on graphical model
    void process_ouputs(impalib_type *, int *, int *, int *, impalib_type *, bool *, bool *, bool *, int *, int *, int *, int *); ///< process outputs of Graphical Model
    GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG,
                      const impalib_type ALPHA, const impalib_type THRESHOLD, const int MAX_COUNT); ///< Constructor
};

/**
 * Construct Graphical Model for the TSP
 *
 * @param[in] NUM_ITERATIONS: number of iterations of IMPA
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of connections between nodes
 * @param[in] AUGMENTATION_FLAG: augmentation flag activated or not
 * @param[in] RESET_FLAG: reset flag for resetting messages after each
 * augmentation step
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[in] THRESHOLD: threshold for hard decision on intrinsic messages
 * @param[in] MAX_COUNT: maximum count of failure cases before exiting the code
 * @param[out] modelDegreeConstraint_: object of type DegreeConstraint: object
 * associated with degree constraint operations
 * @param[out] modelEqConstraint_: object of type EqualityConstraintTsp: object
 * associated with equality constraint operations
 * @param[out] modelSubtourEliminationConstraint_: object of type
 * SubtourEliminationConstraint: object associated with subtour elimination
 * constraint operations
 * @param[out] modelInputs_: object of type InputsTsp: object associated with
 * inputs of tsp
 * @param[out] outputs: object of type OutputsTsp: object associated with
 * outputs of tsp
 * @param[out] numIterations_: NUM_ITERATIONS
 * @param[out] filteringFlag_: FILTERING_FLAG
 * @param[out] alpha_: ALPHA
 * @param[out] numNodes_: NUM_NODES
 * @param[out] augmentationFlag_: AUGMENTATION_FLAG
 * @param[out] resetFlag_: RESET_FLAG
 * @param[out] numEdgeVariables_: NUM_EDGE_VARIABLES
 * @param[out] threshold_: THRESHOLD
 * @param[out] maxCount_: MAX_COUNT
 *
 */

GraphicalModelTsp::GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG,
                                     const impalib_type ALPHA, const impalib_type THRESHOLD, const int MAX_COUNT)
    : modelDegreeConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      modelEqConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      modelSubtourEliminationConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      modelInputs_(NUM_NODES, NUM_EDGE_VARIABLES),
      outputs(NUM_NODES, NUM_EDGE_VARIABLES),
      numIterations_(NUM_ITERATIONS),
      filteringFlag_(FILTERING_FLAG),
      alpha_(ALPHA),
      numNodes_(NUM_NODES),
      augmentationFlag_(AUGMENTATION_FLAG),
      resetFlag_(RESET_FLAG),
      numEdgeVariables_(NUM_EDGE_VARIABLES),
      threshold_(THRESHOLD),
      maxCount_(MAX_COUNT),
      DegreeConstraint2EqConstraintDummyM_(numEdgeVariables_, vector<impalib_type>(numNodes_, 0)),
      DegreeConstraint2EqConstraintM_(numEdgeVariables_, vector<impalib_type>(numNodes_, 0))
{
};

/**
 * Initialize Graphical Model inputs for the TSP. The inputs are processed by
 * reading inputs from python
 *
 * @param[in] pEDGE_CONNECTIONS_PY: contains all possible connections between
 * nodes. Each edge has its constituent nodes
 * @param[in] pCOST_EDGE_VARIABLE_PY: cost for each possible connection between
 * nodes. This has size of number of edges
 * @param[in] pCOST_MATRIX_PY: cost matrix of size number of nodes x number of
 * nodes
 * @param[in] pEdge_ec_to_degree_constraint_m_py: messages from edges equality
 * constraints to degree constraints
 * @param[in] pEDGE_DEGREE_CONSTRAINT_COST_PY: another matrix of costs that has
 * size number of edges x number of nodes. Refer to example of TSP to understand
 * how the various costs differ. The various various will facilitate
 * computations in the IMPA
 *
 */

void GraphicalModelTsp::initialize(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY, impalib_type *pEdge_ec_to_degree_constraint_m_py,
                                   const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
    // Populate model data
    modelInputs_.process_inputs(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY, pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);
}

/**
 * Run IMPA on relaxed Grpahical Model. Only degree constraints are included. No
 * subtour elimination constraints This will propagate messages for a certain
 * number of iterations across the whole graphical model between degree
 * constraints and edge equality constraints
 *
 */
void GraphicalModelTsp::iterate_relaxed_graph() {
    for (int i = 0; i < numIterations_; i++) {
        modelDegreeConstraint_.degree_constraint_to_edge_ec_update(modelInputs_.EdgeEc2DegreeConstraintM, modelInputs_.EdgeConnections, DegreeConstraint2EqConstraintDummyM_);
        modelDegreeConstraint_.process_filtering(i, DegreeConstraint2EqConstraintDummyM_, DegreeConstraint2EqConstraintM_);
        modelEqConstraint_.edge_ec_to_degree_constraint_relaxed_graph_update(modelInputs_.EdgeConnections, modelInputs_.EdgeDegreeConstraintCost, DegreeConstraint2EqConstraintM_,
                                                                             modelInputs_.EdgeEc2DegreeConstraintM);
    }
    outputs.update_extrinsic_relaxed(DegreeConstraint2EqConstraintM_);
    outputs.update_intrinsic(modelInputs_.CostEdgeVariable);

    // Perform hard decision analysis to select edges
    vector<vector<int>> selected_edges;
    hard_decision_analysis(selected_edges);

    // Initialize graph for subtour elimination constraints analysis
    unordered_map<int, vector<int>> graph;
    subtour_elimination_constraints_analysis(graph, selected_edges);

    selectedEdges_ = selected_edges;
}

/**
 * Run IMPA on Augmentated Grpahical Model. Degree constraints are included.
 * Subtour elimination constraints are added if detected This will propagate
 * messages for a certain number of iterations across the whole graphical model
 * between degree constraints, edge equality constraints, and subtour
 * elimination constraints
 * @param[in] MAX_AUGM_COUNT: maximum number of augmentation steps
 *
 */

void GraphicalModelTsp::perform_augmentation(const int MAX_AUGM_COUNT) {
    // Add empty vectors to subtour constraints if needed
    if (delta_S_indices_list.size() > 0) {
        vector<vector<impalib_type>> temp(delta_S_indices_list.size(), vector<impalib_type>(numEdgeVariables_, zero_value));
        subtourConstraints2EdgeEcM_.insert(subtourConstraints2EdgeEcM_.end(), temp.begin(), temp.end());
        subtourConstraints2EdgeEcDummyM_ = subtourConstraints2EdgeEcM_;
    }

    // Continue augmentation
    while (!subtourConstraintsSatisfiedFlag && augmentationFlag_ && !tourImpaFlag_) {
        if (noConsClosedLoopsCountExcFlag_ || noImprovSolCountExcFlag_ || solOscCountExcFlag_) {
            cout << "noConsClosedLoopsCountExcFlag_: " << noConsClosedLoopsCountExcFlag_ << '\n';
            cout << "noImprovSolCountExcFlag_: " << noImprovSolCountExcFlag_ << '\n';
            cout << "solOscCountExcFlag_: " << solOscCountExcFlag_ << '\n';
            break;
        }
        if (costImpa_ == zero_value) {
            cout << "Cost is zero" << '\n';
            cout << "Possibly Nan" << '\n';
            exit(0);
        }

        iterate_augmented_graph();

        if (numAugmentations_ == MAX_AUGM_COUNT) {
            cout << "MAX_AUGM_COUNT reached" << '\n';
            break;
        }

        // Add empty vectors to subtour constraints if needed
        if (subtourConstraints2EdgeEcM_.size() != delta_S_indices_list.size()) {
            size_t numLists2Add = delta_S_indices_list.size() - subtourConstraints2EdgeEcM_.size();
            vector<vector<impalib_type>> temp(numLists2Add, vector<impalib_type>(numEdgeVariables_, zero_value));
            subtourConstraints2EdgeEcM_.insert(subtourConstraints2EdgeEcM_.end(), temp.begin(), temp.end());
            subtourConstraints2EdgeEcDummyM_ = subtourConstraints2EdgeEcM_;
        }
    }
}

/**
 * After IMPA is completed, this will process the outputs of the TSP which will
 * be fed to IMPA
 * @param[out] pExtrinsic_output_edge_ec: output extrinsic messages of edge
 * equality constraints
 * @param[out] pNum_augmentations: number of performed augmentations when IMPA
 * is completed
 * @param[out] pNum_added_constraints: number of added subtour elimination
 * constraints when IMPA is completed
 * @param[out] pTour_impa: number of nodes in tour if found
 * @param[out] pCost_impa: cost obtained from the activated edges (tour or not)
 * @param[out] pNo_improv_sol_count_exc_flag: failure flag when solution does
 * not improve when max_count is reached
 * @param[out] pNo_cons_loops_count_exc_flag: failure flag when IMPA does not
 * detect subtours but fails to have a tour
 * @param[out] pSol_osc_count_exc_flag: failure flag when solution oscillates
 * between two states
 * @param[out] pSelected_edges: activated edges after IMPA is completed
 * @param[out] pSelected_edges_size: number of activated edges
 * @param[out] pSubtour_paths: subtour paths, if available, after IMPA is
 * completed. If not empty, this will be used for post-processing
 * @param[out] pSubtour_paths_size: number of subtour paths
 *
 */

void GraphicalModelTsp::process_ouputs(impalib_type *pExtrinsic_output_edge_ec, int *pNum_augmentations, int *pNum_added_constraints, int *pTour_impa, impalib_type *pCost_impa,
                                       bool *pNo_improv_sol_count_exc_flag, bool *pNo_cons_loops_count_exc_flag, bool *pSol_osc_count_exc_flag, int *pSelected_edges, int *pSelected_edges_size,
                                       int *pSubtour_paths, int *pSubtour_paths_size) {
    // Copy extrinsic output for edge equality constraints
    copy(outputs.ExtrinsicOutputEdgeEc.begin(), outputs.ExtrinsicOutputEdgeEc.begin() + numEdgeVariables_, pExtrinsic_output_edge_ec);
    *pNum_augmentations = numAugmentations_;
    *pNum_added_constraints = static_cast<int>(subtourConstraints2EdgeEcM_.size());

    copy(tourImpa_.begin(), tourImpa_.begin() + static_cast<int>(tourImpa_.size()), pTour_impa);
    *pCost_impa = costImpa_;

    *pNo_improv_sol_count_exc_flag = noImprovSolCountExcFlag_;
    *pNo_cons_loops_count_exc_flag = noConsClosedLoopsCountExcFlag_;
    *pSol_osc_count_exc_flag = solOscCountExcFlag_;

    // Copy selected edges
    vector<int> flattened_selected_edges = accumulate(selectedEdges_.begin(), selectedEdges_.end(), vector<int>{}, [](vector<int> &acc, const vector<int> &inner) {
        acc.insert(acc.end(), inner.begin(), inner.end());
        return acc;
    });

    copy(flattened_selected_edges.begin(), flattened_selected_edges.begin() + static_cast<int>(flattened_selected_edges.size()), pSelected_edges);
    *pSelected_edges_size = static_cast<int>(flattened_selected_edges.size());

    // Copy subtour edges
    vector<int> flattened_closed_paths = accumulate(subtourPaths_.begin(), subtourPaths_.end(), vector<int>{}, [](vector<int> &acc, const vector<int> &inner) {
        acc.insert(acc.end(), inner.begin(), inner.end());
        return acc;
    });
    copy(flattened_closed_paths.begin(), flattened_closed_paths.begin() + static_cast<int>(flattened_closed_paths.size()), pSubtour_paths);
    copy(closedPathsSize_.begin(), closedPathsSize_.begin() + static_cast<int>(closedPathsSize_.size()), pSubtour_paths_size);
}

/**
 * This function will run IMPA on augmented grahical model
 *
 */
void GraphicalModelTsp::iterate_augmented_graph() {
    // Reset messages if flag is set
    if (resetFlag_) {
        for (size_t i = 0; i < modelInputs_.EdgeConnections.size(); i++) {
            auto connection = modelInputs_.EdgeConnections[i];
            impalib_type cost = modelInputs_.CostMatrix[connection[0]][connection[1]];
            modelInputs_.EdgeEc2DegreeConstraintM[i][connection[0]] = cost;
            modelInputs_.EdgeEc2DegreeConstraintM[i][connection[1]] = cost;
        }
    }

    cout << "-----------------------" << '\n';
    cout << "iterate_augmented_graph" << '\n';
    cout << "-----------------------" << '\n';

    // Increment augmentation count
    numAugmentations_ += 1;
    cout << "Augmentation count: " << numAugmentations_ << '\n';
    cout << "delta_S_indices_list.size(): " << delta_S_indices_list.size() << '\n';

    // Reserve memory for old subtour constraints
    modelSubtourEliminationConstraint_.subtourConstraints2EdgeEcOld_.reserve(delta_S_indices_list.size());

    // Initialize old subtour constraints
    for (int i = 0; i < delta_S_indices_list.size(); i++) {
        modelSubtourEliminationConstraint_.subtourConstraints2EdgeEcOld_.push_back(vector<impalib_type>(numEdgeVariables_, zero_value));
    }

    for (int i = 0; i < numIterations_; i++) {
        modelDegreeConstraint_.degree_constraint_to_edge_ec_update(modelInputs_.EdgeEc2DegreeConstraintM, modelInputs_.EdgeConnections, DegreeConstraint2EqConstraintDummyM_);
        modelDegreeConstraint_.process_filtering(i, DegreeConstraint2EqConstraintDummyM_, DegreeConstraint2EqConstraintM_);

        edgeEc2SubtourConstraintsM_ = modelEqConstraint_.edge_ec_to_subtour_constraints_update(delta_S_indices_list, modelInputs_.CostEdgeVariable, DegreeConstraint2EqConstraintM_,
                                                                                               subtourConstraints2EdgeEcM_, modelInputs_.EdgeConnections);
        modelSubtourEliminationConstraint_.subtour_constraints_to_edge_ec_update(edgeEc2SubtourConstraintsM_, delta_S_indices_list, subtourConstraints2EdgeEcDummyM_);
        modelSubtourEliminationConstraint_.process_filtering(i, subtourConstraints2EdgeEcDummyM_, subtourConstraints2EdgeEcM_, delta_S_indices_list);
        modelEqConstraint_.edge_ec_to_degree_constraint_augmented_graph_update(DegreeConstraint2EqConstraintM_, subtourConstraints2EdgeEcM_, modelInputs_.EdgeConnections,
                                                                               modelInputs_.EdgeDegreeConstraintCost, modelInputs_.EdgeEc2DegreeConstraintM);
    }

    // Check for NaN values
    bool flag_nan = false;
    for (int i = 0; i < modelInputs_.EdgeEc2DegreeConstraintM.size(); i++) {
        for (int j = 0; j < modelInputs_.EdgeEc2DegreeConstraintM[i].size(); j++) {
            if (isnan(modelInputs_.EdgeEc2DegreeConstraintM[i][j])) {
                flag_nan = true;
            }
        }
    }

    if (flag_nan) {
        cout << "NaN detected" << '\n';
        exit(0);
    }

    outputs.update_extrinsic_augmented(DegreeConstraint2EqConstraintM_, subtourConstraints2EdgeEcM_);
    outputs.update_intrinsic(modelInputs_.CostEdgeVariable);

    selectedEdges_.clear();
    vector<vector<int>> selected_edges;
    hard_decision_analysis(selected_edges);

    // Check for solution improvement
    if (!selected_edges_old_.empty()) {
        bool areEqual = (selected_edges == selected_edges_old_);
        if (areEqual) {
            cout << "Exited: No Improvement of IMPA Solution" << '\n';
            noImprovSolCount_ += 1;
            if (noImprovSolCount_ > maxCount_) {
                cout << "Exited: noImprovSolCount_ Exceeded maximum allowable "
                        "number "
                     << maxCount_ << '\n';
                noImprovSolCountExcFlag_ = true;
            }
        } else {
            noImprovSolCount_ = 0;
        }
    }

    // Check for solution oscillation
    if (!selected_edges_old_old_.empty()) {
        bool oscillation_flag = ((selected_edges == selected_edges_old_old_) && (selected_edges_old_old_ != selected_edges_old_));
        if (oscillation_flag) {
            cout << "Exited: Solution is oscillating" << '\n';
            solOscCount_ += 1;
            if (solOscCount_ > maxCount_) {
                cout << "Exited: solOscCount_ Exceeded maximum allowable number " << maxCount_ << '\n';
                solOscCountExcFlag_ = true;
            }
        } else {
            solOscCount_ = 0;
        }
    }

    if (!selected_edges_old_.empty()) {
        selected_edges_old_old_ = selected_edges_old_;
    }
    selected_edges_old_ = selected_edges;
    selectedEdges_ = selected_edges;
    unordered_map<int, vector<int>> graph;
    subtour_elimination_constraints_analysis(graph, selected_edges);
}

/**
 * This will get the hard decision on edges by investigating the sign
 * IntrinsicOutputEdgeEc
 * @param[out] rSelectedEdges: activated edges after running IMPA
 *
 */
void GraphicalModelTsp::hard_decision_analysis(vector<vector<int>> &rSelectedEdges) {
    hard_decision.resize(numEdgeVariables_);
    fill(hard_decision.begin(), hard_decision.begin() + numEdgeVariables_, numeric_limits<int>::max());

    // Apply threshold to intrinsic output to determine hard decision
    transform(outputs.IntrinsicOutputEdgeEc.begin(), outputs.IntrinsicOutputEdgeEc.end(), hard_decision.begin(), [&](impalib_type value) { return value > threshold_ ? zero_value : 1; });

    // Select edges based on hard decision
    for (int i = 0; i < numEdgeVariables_; i++) {
        if (hard_decision[i] == 1) {
            rSelectedEdges.push_back(modelInputs_.EdgeConnections[i]);
        }
    }

    // Collect unique nodes from selected edges
    set<int> uniqueNodes;

    cout << "selected_edges: [";
    for (const vector<int> &edge : rSelectedEdges) {
        for (int node : edge) {
            if (uniqueNodes.find(node) == uniqueNodes.end()) {
                uniqueNodes.insert(node);
            }
        }
        cout << "[";
        for (size_t i = 0; i < edge.size(); i++) {
            cout << edge[i];
            if (i != edge.size() - 1) {
                cout << ", ";
            }
        }

        cout << "]";

        if (&edge != &rSelectedEdges.back()) {
            cout << ", ";
        }
    }
    cout << "]" << '\n';
    // cout << "Number of activated nodes: " << uniqueNodes.size()  << '\n';

    // Calculate cost based on selected edges
    impalib_type cost_impa = zero_value;
    for (size_t i = 0; i < hard_decision.size(); ++i) {
        if (hard_decision[i] == 1) {
            cost_impa += modelInputs_.CostEdgeVariable[i];
        }
    }

    cout << "C++ cost_impa: " << cost_impa << '\n';
    costImpa_ = cost_impa;
}

/**
 * Analyze the activated edges. If tour found, return tour. If tour not found,
 * detect subtours. Count failure cases if present
 * @param[in] rSelectedEdges: activated edges in the graphical model
 * @param[out] rGraph: graphical model of activated egdes. Defines connections
 * between nodes
 *
 */

void GraphicalModelTsp::subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &rGraph, vector<vector<int>> &rSelectedEdges) {
    closedPathsSize_.clear();  // just store it at the end if applicable

    // Get closed loops from the graph
    vector<vector<int>> loops_list = get_closed_loops(rGraph, rSelectedEdges);
    subtourPaths_.clear();
    subtourPaths_ = loops_list;

    if (loops_list.empty()) {
        cout << "Exited: Cannot find tours using get_closed_loops()" << '\n';
        noConsClosedLoopsCount_ += 1;
        if (noConsClosedLoopsCount_ > maxCount_) {
            cout << "Exited: noConsClosedLoopsCount_ Exceeded maximum allowable "
                    "number: "
                 << maxCount_ << '\n';
            noConsClosedLoopsCountExcFlag_ = true;
        }
    }

    else if (loops_list.size() == 1 && loops_list[0].size() == numNodes_) {
        tourImpaFlag_ = true;
        subtourConstraintsSatisfiedFlag = true;
        tourImpa_ = loops_list[0];
        tourImpa_.push_back(loops_list[0][0]);
        cout << "Tour found" << '\n';

        cout << "tour_impa: [";
        for (size_t i = 0; i < tourImpa_.size(); ++i) {
            cout << tourImpa_[i];
            if (i != tourImpa_.size() - 1) {
                cout << ", ";
            }
        }
        cout << "]" << '\n';
    }
    // If multiple subtours are found
    else {
        noConsClosedLoopsCount_ = 0;
        for (const auto &loop : loops_list) {
            // Store the size of each subtour
            closedPathsSize_.push_back(static_cast<int>(loop.size()));
            cout << "Subtour of size " << loop.size() << " detected @: ";
            cout << "[";
            for (size_t i = 0; i < loop.size(); ++i) {
                cout << loop[i];
                if (i != loop.size() - 1) {
                    cout << ", ";
                }
            }
            cout << "]" << '\n';

            // Find delta S indices for each subtour
            vector<int> delta_S_indices;
            for (size_t i = 0; i < modelInputs_.EdgeConnections.size(); ++i) {
                const auto &connection = modelInputs_.EdgeConnections[i];
                if (find(loop.begin(), loop.end(), connection[0]) != loop.end() && find(loop.begin(), loop.end(), connection[1]) == loop.end()) {
                    delta_S_indices.push_back(static_cast<int>(i));
                }
            }
            // Add delta S indices to delta_S_indices_list
            if (delta_S_indices_list.empty() && !delta_S_indices.empty()) {  // added !delta_S_indices.empty() to make
                                                                             // sure if a full tour was detected so it
                                                                             // is not be added with subtours
                delta_S_indices_list.push_back(delta_S_indices);
            } else {
                // Add delta S indices if not already present
                bool found = false;
                for (const auto &existing_indices : delta_S_indices_list) {
                    if (existing_indices == delta_S_indices) {
                        found = true;
                        break;
                    }
                }
                if (!found && !delta_S_indices.empty()) {  // added !delta_S_indices.empty() to
                                                           // make sure if a full tour was detected
                                                           // so it is not be added with subtours
                    delta_S_indices_list.push_back(delta_S_indices);
                }
            }
        }
    }
}

/**
 * Get closed loops. Function for building the graphical model of activated
 * edges and obtain loops
 * @param[in] rSelectedEdges: activated edges in the graphical model
 * @param[out] rGraph: graphical model of activated egdes. Defines connections
 * between nodes. Will be used for detecting of subtours
 * @return new_loops_list: list of detected loops in rGraph
 *
 */
vector<vector<int>> GraphicalModelTsp::get_closed_loops(unordered_map<int, vector<int>> &rGraph, vector<vector<int>> &rSelectedEgdes) {
    // Update the graph based on selected edges
    for (const auto &connection : rSelectedEgdes) {
        if (rGraph.find(connection[0]) != rGraph.end()) {
            rGraph[connection[0]].push_back(connection[1]);
        } else {
            rGraph[connection[0]] = {connection[1]};
        }
    }

    vector<vector<int>> loops_list;
    vector<int> visited_nodes;

    for (const auto &connection : rGraph) {
        int start_node = connection.first;
        const vector<int> &end_nodes = connection.second;
        // Skip if node is already visited
        if (find(visited_nodes.begin(), visited_nodes.end(), start_node) != visited_nodes.end()) {
            continue;
        } else {
            for (int j = 0; j < end_nodes.size(); j++) {
                unordered_set<int> visited_set;
                vector<int> path;

                vector<int> closed_loop = find_closed_loop(rGraph, start_node, end_nodes[j], visited_set, path, visited_nodes);
                if (!closed_loop.empty()) {
                    loops_list.push_back(closed_loop);
                }
            }
        }
    }

    // Remove duplicate/subset loops
    vector<int> list_indices_to_remove;
    vector<vector<int>> new_loops_list;

    for (int i = 0; i < static_cast<int>(loops_list.size()); i++) {
        const auto &list_1 = loops_list[i];
        if (list_1.empty()) {
            list_indices_to_remove.push_back(i);
            continue;
        } else if (find(list_indices_to_remove.begin(), list_indices_to_remove.end(), i) != list_indices_to_remove.end()) {
            continue;
        }
        vector<int> double_list_1 = list_1;
        double_list_1.insert(double_list_1.end(), list_1.begin(), list_1.end());
        // Compare current loop with others to check for duplicates/subsets
        for (int j = i + 1; j < static_cast<int>(loops_list.size()); j++) {
            const auto &list_2 = loops_list[j];

            if (list_2.empty()) {
                list_indices_to_remove.push_back(j);
            }

            if (double_list_1.size() < list_2.size()) {  // since if seq.size() - subseq.size()<0, this
                                                         // would give an arbitrary large number and
                                                         // core dumped
                continue;
            }

            else if (isSubsequence(double_list_1, list_2, j)) {
                list_indices_to_remove.push_back(j);
            }
        }
    }

    // Collect non-duplicate/non-subset loops
    for (size_t i = 0; i < loops_list.size(); ++i) {
        if (find(list_indices_to_remove.begin(), list_indices_to_remove.end(), i) == list_indices_to_remove.end()) {
            new_loops_list.push_back(loops_list[i]);
        }
    }

    return new_loops_list;
}

/**
 * Function used during the detection of loops. This function efficiently checks
 * for the presence of a subsequence within a larger sequence by iterating over
 * all possible starting positions for the subsequence and comparing elements at
 * corresponding positions
 * @param[in] seq: subsequence
 * @param[in] subseq: subsequence
 * @param[in] j: index in loops_list. For printing purposes, if needed. It is
 * not used in the code.
 * @return false or true if subseq is a isSubsequence of seq
 *
 */
bool GraphicalModelTsp::isSubsequence(const vector<int> &seq, const vector<int> &subseq, int j) {
    for (int i = 0; i < static_cast<int>(seq.size() - subseq.size()); ++i) {
        if (equal(seq.begin() + i, seq.begin() + i + static_cast<int>(subseq.size()), subseq.begin())) {
            return true;
        }
    }
    return false;
}

/**
 * This function will be called to find closed loops. This is called multiple
 * times as shown in get_closed_loops function. It will find path between nodes.
 * Will return any path (which could be a tour)
 * @param[in] rGraph: graph of activated edges. Mapping between nodes and their
 * neighboring nodes
 * @param[in] start_node: node the path is starting from
 * @param[in] current: current node under investigation in the path
 * @param[in] visited_nodes: This includes all visited nodes. This can help in
 * skipping the investigation of nodes that are already in the path to avoid
 * processing of loops list in get_closed_loops. This was deactivated in this
 * code, and can be used in the future to reduce processing of detected loops.
 * @param[out] visited: visited nodes while constructing the path
 * @param[out] path: current detected path, which will be augmented
 * @return new_path or empty vector if no path is found. new_path is a path
 * between the nodes (if a tour is found, )
 *
 */
vector<int> GraphicalModelTsp::find_closed_loop(unordered_map<int, vector<int>> &rGraph, int start_node, int current, unordered_set<int> visited, vector<int> path, vector<int> &visited_nodes) {
    visited.insert(current);
    path.push_back(current);

    // Check if loop is found
    if (current == start_node) {
        return path;
    }

    // If current node has no outgoing connections
    if (rGraph.find(current) == rGraph.end()) {
        vector<int> new_path = {start_node};
        new_path.insert(new_path.end(), path.begin(), path.end());
        return vector<int>();
    }

    for (int node : rGraph.at(current)) {
        // If node is not visited, continue search
        if (visited.find(node) == visited.end()) {
            vector<int> new_path = find_closed_loop(rGraph, start_node, node, visited, path, visited_nodes);
            // Return new path if loop is found
            if (!new_path.empty()) {
                return new_path;
            }
        }
    }

    return vector<int>();
}