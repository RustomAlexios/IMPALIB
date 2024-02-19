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
    int numProjects_;            ///< number of projects
    int numTeams_;               ///< number of teams
    int numDepartments_;         ///< number of departments
    int maxSizeNonZeroWeights_;  ///< maximum size of non-zero weights per department
    int numIterations_;          ///< number of iterations of IMPA
    bool doFilter_;              ///< filtering flag of knapsack constraints
    impalib_type alpha_;         ///< filtering parameter

    vector<vector<impalib_type>> extrinsic_;  ///< messages from departments to teams after filtering
    vector<vector<impalib_type>> M_eq2oric_;  ///< messages from team equality constraint to ORIC
    vector<vector<impalib_type>> M_oric2eq_;  ///< messages from ORIC to team equality constraint
    vector<vector<impalib_type>> M_eq2ineq_;  ///< messages from project equality constraint to project inequality constraint
    vector<vector<impalib_type>> M_ineq2eq_;  ///< messages from project inequality constraint to project equality constraint

    Knapsack knapsack_;                 ///< Knapsack object
    InequalityConstraint projectIneq_;  ///< Project Inequality constraint object
    EqualityConstraintKcMwm EqKcMwm_;   ///< Equality constraint object
    OrInequalityConstraint Oric_;       ///< ORIC object

   public:
    OutputsKcMwm outputs;                                                                                                                                      ///< Graphical model outputs objects
    InputsKcMwm modelInputs_;                                                                                                                                  ///< Graphical model inputs objects
    void initialize(const impalib_type *, const impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *);                     ///< initialize graphical model
    void iterate(const int *);                                                                                                                                 ///< iterate over graphical model
    GraphicalModelKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS, int MAX_SIZE_NON_ZERO_WEIGHTS, int N_ITERATIONS, bool FILT_FLAG, impalib_type ALPHA);  ///< constructor
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
      knapsack_(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA),
      projectIneq_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      EqKcMwm_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      Oric_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      numProjects_(N_PROJECTS),
      numTeams_(N_TEAMS),
      numDepartments_(N_DEPARTMENTS),
      maxSizeNonZeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS),
      numIterations_(N_ITERATIONS),
      doFilter_(FILT_FLAG),
      alpha_(ALPHA),
      M_eq2oric_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      M_oric2eq_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      M_eq2ineq_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      M_ineq2eq_(numProjects_, vector<impalib_type>(numTeams_, 0)),
      extrinsic_(numDepartments_, vector<impalib_type>(numTeams_, 0)){};

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

void GraphicalModelKcMwm::initialize(const impalib_type *pREWARD_TEAM_PY, const impalib_type *pTransition_model_py, const int *pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                                     const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
    /// calls a method process_inputs() on an object inputs_, passing
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
    vector<impalib_type> oric2PackageM_(numTeams_, 0);
    for (int i = 0; i < numIterations_; i++) {
        for (int j = 0; j < numDepartments_; j++) {
            int max_state_department = modelInputs_.capacity[j];
            auto &idx_nonzero_dept = modelInputs_.NonZeroWeightIndices[j];
            auto &team_weights = modelInputs_.weights[j];

            auto stage_forward_messages = knapsack_.forward(j, max_state_department, idx_nonzero_dept, pNON_ZERO_WEIGHT_INDICES_SIZES_PY, team_weights, modelInputs_.M_team2knapsack);
            auto stage_backward_messages = knapsack_.backward(j, max_state_department, idx_nonzero_dept, pNON_ZERO_WEIGHT_INDICES_SIZES_PY, team_weights, modelInputs_.M_team2knapsack);

            auto extrinsic_out = knapsack_.extrinsic_output_department_lhs(team_weights, stage_forward_messages, modelInputs_.M_team2knapsack, j, stage_backward_messages, max_state_department);
            knapsack_.process_extrinsic_output_department(j, i, extrinsic_out, extrinsic_[j]);
        }

        auto team2OricM_ = EqKcMwm_.team_messages_to_oric(extrinsic_, modelInputs_.RewardTeam);
        Oric_.messages_to_project_eq(M_eq2oric_, team2OricM_, M_oric2eq_, M_eq2ineq_, modelInputs_.RewardProject);
        M_ineq2eq_ = projectIneq_.messages_to_equality(M_eq2ineq_);
        EqKcMwm_.project_messages_to_oric(M_ineq2eq_, M_eq2oric_, modelInputs_.RewardProject);

        oric2PackageM_ = Oric_.messages_to_team_eq(M_eq2oric_);
        outputs.update_intrinsic(M_oric2eq_, M_ineq2eq_, modelInputs_.RewardProject);
        knapsack_.team_to_knapsack_update(modelInputs_.NonZeroWeightIndices, modelInputs_.M_team2knapsack, modelInputs_.RewardTeam, extrinsic_, oric2PackageM_);
    }
    outputs.update_extrinsic(extrinsic_, oric2PackageM_);
}

/**
 * Represents a graphical model for the TSP for computing the shortest route
 * between nodes (while visiting each node exactly once).
 */
class GraphicalModelTsp {
   private:
    int numNodes_;                                       ///< number of nodes of TSP
    int numIterations_;                                  ///< number of iterations of IMPA
    int numEdgeVariables_;                               ///< number of edges in TSP
    bool doReset_;                                       ///< flag for resetting messages after each augmentation step
    bool doFilter_;                                      ///< flag for filtering messages from degree/subtour constraints to edge equality constraints
    bool doAugment_;                                     ///< whether to activate augmentation
    vector<int> hard_decision;                           ///< hard decision vector of IMPA solution
    impalib_type alpha_;                                 ///< filtering parameter
    impalib_type threshold_;                             ///< threshold on hard decision
    EqualityConstraintTsp EqTsp_;                        ///< Equality Constraint object for TSP
    DegreeConstraint degree;                             ///< Degree Constraint object for TSP
    SubtourEliminationConstraint subtour;                ///< Subtour constraint for TSP
    vector<vector<impalib_type>> M_degree2eq_before;     ///< messages from degree constraint to team equality constraint before filtering
    vector<vector<impalib_type>> M_degree2eq;            ///< messages from degree constraint to team equality constraint afters filtering
    vector<vector<int>> selected_edges_old_;             ///< old list of selected edges (used for failure investigation)
    vector<vector<int>> selected_edges_old_old_;         ///< another old list of selected edges (used for failure investigation)
    int maxCount_;                                       ///< maximum count of failures
    bool tourImpaFlag_ = false;                          ///< set flag for detected tour to be false (will be updated later)
    vector<vector<int>> idx_delta_S;                     ///< will contain indices of edges that contribute to each subtour constraint
    vector<vector<impalib_type>> M_subtour2edge_before;  ///< messages from subtour constraints to edge equality constraint before filtering
    vector<vector<impalib_type>> M_edge2subtour;         ///< messages from edge equality constraint to subtour constraints
    void iterate_augmented_graph();                      ///< function of IMPA on augmented graph
    void subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &, const vector<vector<int>> &);                    ///< analysis of subtour constraints
    vector<vector<int>> hard_decision_analysis();                                                                                     ///< function for hard decision solution on IMPA solution
    bool isSubsequence(const vector<int> &, const vector<int> &, int);                                                                ///< function for post-processing loops
    vector<vector<int>> get_closed_loops(unordered_map<int, vector<int>> &, const vector<vector<int>> &);                             ///< function for getting loops
    vector<int> find_closed_loop(const unordered_map<int, vector<int>> &, int, int, unordered_set<int>, vector<int>, vector<int> &);  ///< function for finding loops
    InputsTsp inputs_;                                                                                                                ///< Graphical Model Input object
    vector<vector<int>> selectedEdges_;                                                                                               ///< activated edges of IMPA
    int numAugmentations_ = 0;                                                                                                        ///< number of performed augmentations in IMPA
    int noConsClosedLoopsCount_ = 0;              ///< count for failure case (no consecutive loop detection and no tour)
    int n_osc = 0;                                ///< count for failure case (oscillation in the solution)
    bool n_closedLoops = false;                   ///< flag for failure case (no consecutive loop detection and no tour)
    int n_noimprove = 0;                          ///< count for failure case (no solution improvement)
    bool noimprove = false;                       ///< flag for failure case (no solution improvement)
    bool osc = false;                             ///< flag for failure case (oscillation in the solution)
    vector<int> tour;                             ///< list of nodes of tour (if detected)
    vector<vector<int>> subtours;                 ///< list of detected subtours (if detected)
    vector<int> loops;                            ///< list of sizes of loops
    impalib_type cost = zero_value;               ///< cost of IMPA solution
    vector<vector<impalib_type>> M_subtour2edge;  ///< messages from subtour constraints to edge equality constraint

   public:
    bool subtourConstraintsSatisfiedFlag = false;                                                                          ///< initially set this flag for satisfied subtour constraints to false
    OutputsTsp outputs;                                                                                                    ///< TSP graphical model outputs object
    void initialize(const int *, const impalib_type *, const impalib_type *, const impalib_type *, const impalib_type *);  ///< initialize graphical model
    void iterate_relaxed_graph();                                                                                          ///< iterate over relaxed graphical model
    void perform_augmentation(int);                                                                                        ///< perform augmentation on graphical model
    void process_ouputs(impalib_type *, int *, int *, int *, impalib_type *, bool *, bool *, bool *, int *, int *, int *, int *);  ///< process outputs of Graphical Model
    GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG,
                      const impalib_type ALPHA, const impalib_type THRESHOLD, const int MAX_COUNT);  ///< Constructor
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
    : degree(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      EqTsp_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      subtour(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      inputs_(NUM_NODES, NUM_EDGE_VARIABLES),
      outputs(NUM_NODES, NUM_EDGE_VARIABLES),
      numIterations_(NUM_ITERATIONS),
      doFilter_(FILTERING_FLAG),
      alpha_(ALPHA),
      numNodes_(NUM_NODES),
      doAugment_(AUGMENTATION_FLAG),
      doReset_(RESET_FLAG),
      numEdgeVariables_(NUM_EDGE_VARIABLES),
      threshold_(THRESHOLD),
      maxCount_(MAX_COUNT),
      M_degree2eq_before(numEdgeVariables_, vector<impalib_type>(numNodes_, 0)),
      M_degree2eq(numEdgeVariables_, vector<impalib_type>(numNodes_, 0)){};

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

void GraphicalModelTsp::initialize(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY,
                                   const impalib_type *pEdge_ec_to_degree_constraint_m_py, const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
    // Populate model data
    inputs_.process_inputs(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY, pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);
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
        degree.messages_to_edge_ec(inputs_.M_edge2degree, inputs_.edges, M_degree2eq_before);
        degree.process_filtering(i, M_degree2eq_before, M_degree2eq);
        EqTsp_.messages_to_degree_relaxed(inputs_.edges, inputs_.cost_edge_mat, M_degree2eq, inputs_.M_edge2degree);
    }
    outputs.update_extrinsic_relaxed(M_degree2eq);
    outputs.update_intrinsic(inputs_.cost_edge);

    auto selected_edges = hard_decision_analysis();

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
    if (idx_delta_S.size() > 0) {
        vector<vector<impalib_type>> temp(idx_delta_S.size(), vector<impalib_type>(numEdgeVariables_, zero_value));
        M_subtour2edge.insert(M_subtour2edge.end(), temp.begin(), temp.end());
        M_subtour2edge_before = M_subtour2edge;
    }

    // Continue augmentation
    while (!subtourConstraintsSatisfiedFlag && doAugment_ && !tourImpaFlag_) {
        if (n_closedLoops || noimprove || osc) {
            cout << "n_closedLoops: " << n_closedLoops << '\n';
            cout << "noimprove: " << noimprove << '\n';
            cout << "osc: " << osc << '\n';
            break;
        }
        if (cost == zero_value) {
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
        if (M_subtour2edge.size() != idx_delta_S.size()) {
            size_t numLists2Add = idx_delta_S.size() - M_subtour2edge.size();
            vector<vector<impalib_type>> temp(numLists2Add, vector<impalib_type>(numEdgeVariables_, zero_value));
            M_subtour2edge.insert(M_subtour2edge.end(), temp.begin(), temp.end());
            M_subtour2edge_before = M_subtour2edge;
        }
    }
}

/**
 * After IMPA is completed, this will process the outputs of the TSP which will
 * be fed to IMPA
 * @param[out] pExtrinsic_output_edge_ec: output extrinsic_ messages of edge
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
    // Copy extrinsic_ output for edge equality constraints
    copy(outputs.extrinsic.begin(), outputs.extrinsic.begin() + numEdgeVariables_, pExtrinsic_output_edge_ec);
    *pNum_augmentations = numAugmentations_;
    *pNum_added_constraints = static_cast<int>(M_subtour2edge.size());

    copy(tour.begin(), tour.begin() + static_cast<int>(tour.size()), pTour_impa);
    *pCost_impa = cost;

    *pNo_improv_sol_count_exc_flag = noimprove;
    *pNo_cons_loops_count_exc_flag = n_closedLoops;
    *pSol_osc_count_exc_flag = osc;

    // Copy selected edges
    vector<int> flattened_selected_edges = accumulate(selectedEdges_.begin(), selectedEdges_.end(), vector<int>{}, [](vector<int> &acc, const vector<int> &inner) {
        acc.insert(acc.end(), inner.begin(), inner.end());
        return acc;
    });

    copy(flattened_selected_edges.begin(), flattened_selected_edges.begin() + static_cast<int>(flattened_selected_edges.size()), pSelected_edges);
    *pSelected_edges_size = static_cast<int>(flattened_selected_edges.size());

    // Copy subtour edges
    vector<int> flattened_closed_paths = accumulate(subtours.begin(), subtours.end(), vector<int>{}, [](vector<int> &acc, const vector<int> &inner) {
        acc.insert(acc.end(), inner.begin(), inner.end());
        return acc;
    });
    copy(flattened_closed_paths.begin(), flattened_closed_paths.begin() + static_cast<int>(flattened_closed_paths.size()), pSubtour_paths);
    copy(loops.begin(), loops.begin() + static_cast<int>(loops.size()), pSubtour_paths_size);
}

/**
 * This function will run IMPA on augmented grahical model
 *
 */
void GraphicalModelTsp::iterate_augmented_graph() {
    // Reset messages if flag is set
    if (doReset_) {
        for (size_t i = 0; i < inputs_.edges.size(); i++) {
            auto connection = inputs_.edges[i];
            impalib_type cost = inputs_.cost[connection[0]][connection[1]];
            inputs_.M_edge2degree[i][connection[0]] = cost;
            inputs_.M_edge2degree[i][connection[1]] = cost;
        }
    }

    cout << "-----------------------" << '\n';
    cout << "iterate_augmented_graph" << '\n';
    cout << "-----------------------" << '\n';

    // Increment augmentation count
    numAugmentations_ += 1;
    cout << "Augmentation count: " << numAugmentations_ << '\n';
    cout << "idx_delta_S.size(): " << idx_delta_S.size() << '\n';

    // Reserve memory for old subtour constraints
    subtour.M_subtour2edge_old.reserve(idx_delta_S.size());

    // Initialize old subtour constraints
    for (int i = 0; i < idx_delta_S.size(); i++) {
        subtour.M_subtour2edge_old.push_back(vector<impalib_type>(numEdgeVariables_, zero_value));
    }

    for (int i = 0; i < numIterations_; i++) {
        degree.messages_to_edge_ec(inputs_.M_edge2degree, inputs_.edges, M_degree2eq_before);
        degree.process_filtering(i, M_degree2eq_before, M_degree2eq);

        M_edge2subtour = EqTsp_.messages_to_subtour(idx_delta_S, inputs_.cost_edge, M_degree2eq, M_subtour2edge, inputs_.edges);
        subtour.messages_to_edge_ec(M_edge2subtour, idx_delta_S, M_subtour2edge_before);
        subtour.process_filtering(i, M_subtour2edge_before, M_subtour2edge, idx_delta_S);
        EqTsp_.messages_to_degree_augmented(M_degree2eq, M_subtour2edge, inputs_.edges, inputs_.cost_edge_mat, inputs_.M_edge2degree);
    }

    // Check for NaN values
    bool flag_nan = false;
    for (int i = 0; i < inputs_.M_edge2degree.size(); i++) {
        for (int j = 0; j < inputs_.M_edge2degree[i].size(); j++) {
            if (isnan(inputs_.M_edge2degree[i][j])) {
                flag_nan = true;
            }
        }
    }

    if (flag_nan) {
        cout << "NaN detected" << '\n';
        exit(0);
    }

    outputs.update_extrinsic_augmented(M_degree2eq, M_subtour2edge);
    outputs.update_intrinsic(inputs_.cost_edge);

    selectedEdges_.clear();
    auto selected_edges = hard_decision_analysis();

    // Check for solution improvement
    if (!selected_edges_old_.empty()) {
        bool areEqual = (selected_edges == selected_edges_old_);
        if (areEqual) {
            cout << "Exited: No Improvement of IMPA Solution" << '\n';
            n_noimprove += 1;
            if (n_noimprove > maxCount_) {
                cout << "Exited: n_noimprove Exceeded maximum allowable "
                        "number "
                     << maxCount_ << '\n';
                noimprove = true;
            }
        } else {
            n_noimprove = 0;
        }
    }

    // Check for solution oscillation
    if (!selected_edges_old_old_.empty()) {
        bool oscillation_flag = ((selected_edges == selected_edges_old_old_) && (selected_edges_old_old_ != selected_edges_old_));
        if (oscillation_flag) {
            cout << "Exited: Solution is oscillating" << '\n';
            n_osc += 1;
            if (n_osc > maxCount_) {
                cout << "Exited: n_osc Exceeded maximum allowable number " << maxCount_ << '\n';
                osc = true;
            }
        } else {
            n_osc = 0;
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
 * intrinsic
 * @param[out] rSelectedEdges: activated edges after running IMPA
 *
 */
vector<vector<int>> GraphicalModelTsp::hard_decision_analysis() {
    vector<vector<int>> selected;
    hard_decision.resize(numEdgeVariables_);
    fill(hard_decision.begin(), hard_decision.begin() + numEdgeVariables_, numeric_limits<int>::max());

    // Apply threshold to intrinsic output to determine hard decision
    transform(outputs.intrinsic.begin(), outputs.intrinsic.end(), hard_decision.begin(), [&](impalib_type value) { return value > threshold_ ? zero_value : 1; });

    // Select edges based on hard decision
    for (int i = 0; i < numEdgeVariables_; i++) {
        if (hard_decision[i] == 1) {
            selected.push_back(inputs_.edges[i]);
        }
    }

    // Collect unique nodes from selected edges
    set<int> uniqueNodes;

    cout << "selected_edges: [";
    for (const vector<int> &edge : selected) {
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

        if (&edge != &selected.back()) {
            cout << ", ";
        }
    }
    cout << "]" << '\n';
    // cout << "Number of activated nodes: " << uniqueNodes.size()  << '\n';

    // Calculate cost based on selected edges
    impalib_type cost_impa = zero_value;
    for (size_t i = 0; i < hard_decision.size(); ++i) {
        if (hard_decision[i] == 1) {
            cost_impa += inputs_.cost_edge[i];
        }
    }

    cout << "C++ cost_impa: " << cost_impa << '\n';
    cost = cost_impa;

    return selected;
}

/**
 * Analyze the activated edges. If tour found, return tour. If tour not found,
 * detect subtours. Count failure cases if present
 * @param[in] rSelectedEdges: activated edges in the graphical model
 * @param[out] rGraph: graphical model of activated egdes. Defines connections
 * between nodes
 *
 */

void GraphicalModelTsp::subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &rGraph, const vector<vector<int>> &rSelectedEdges) {
    loops.clear();  // just store it at the end if applicable

    // Get closed loops from the graph
    vector<vector<int>> loops_list = get_closed_loops(rGraph, rSelectedEdges);
    subtours.clear();
    subtours = loops_list;

    if (loops_list.empty()) {
        cout << "Exited: Cannot find tours using get_closed_loops()" << '\n';
        noConsClosedLoopsCount_ += 1;
        if (noConsClosedLoopsCount_ > maxCount_) {
            cout << "Exited: noConsClosedLoopsCount_ Exceeded maximum allowable "
                    "number: "
                 << maxCount_ << '\n';
            n_closedLoops = true;
        }
    }

    else if (loops_list.size() == 1 && loops_list[0].size() == numNodes_) {
        tourImpaFlag_ = true;
        subtourConstraintsSatisfiedFlag = true;
        tour = loops_list[0];
        tour.push_back(loops_list[0][0]);
        cout << "Tour found" << '\n';

        cout << "tour_impa: [";
        for (size_t i = 0; i < tour.size(); ++i) {
            cout << tour[i];
            if (i != tour.size() - 1) {
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
            loops.push_back(static_cast<int>(loop.size()));
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
            for (size_t i = 0; i < inputs_.edges.size(); ++i) {
                const auto &connection = inputs_.edges[i];
                if (find(loop.begin(), loop.end(), connection[0]) != loop.end() && find(loop.begin(), loop.end(), connection[1]) == loop.end()) {
                    delta_S_indices.push_back(static_cast<int>(i));
                }
            }
            // Add delta S indices to idx_delta_S
            if (idx_delta_S.empty() && !delta_S_indices.empty()) {  // added !idx_delta_S.empty() to make
                                                                    // sure if a full tour was detected so it
                                                                    // is not be added with subtours
                idx_delta_S.push_back(delta_S_indices);
            } else {
                // Add delta S indices if not already present
                bool found = false;
                for (const auto &existing_indices : idx_delta_S) {
                    if (existing_indices == delta_S_indices) {
                        found = true;
                        break;
                    }
                }
                if (!found && !delta_S_indices.empty()) {  // added !idx_delta_S.empty() to
                                                           // make sure if a full tour was detected
                                                           // so it is not be added with subtours
                    idx_delta_S.push_back(delta_S_indices);
                }
            }
        }
    }
}

/**
 * Get closed loops. Function for building the graphical model of activated
 * edges and obtain loops
 * @param[in] rSelectedEdges: activated edges in the graphical model
 * @param[out] G: graphical model of activated egdes. Defines connections
 * between nodes. Will be used for detecting of subtours
 * @return new_loops_list: list of detected loops in G
 *
 */
vector<vector<int>> GraphicalModelTsp::get_closed_loops(unordered_map<int, vector<int>> &G, const vector<vector<int>> &edges) {
    // Update the graph based on selected edges
    for (const auto &connection : edges) {
        if (G.find(connection[0]) != G.end()) {
            G[connection[0]].push_back(connection[1]);
        } else {
            G[connection[0]] = {connection[1]};
        }
    }

    vector<vector<int>> loops;
    vector<int> visited;

    for (const auto &connection : G) {
        int start_node = connection.first;
        const vector<int> &end_nodes = connection.second;
        // Skip if node is already visited
        if (find(visited.begin(), visited.end(), start_node) != visited.end()) {
            continue;
        } else {
            for (int j = 0; j < end_nodes.size(); j++) {
                unordered_set<int> visited_set;
                vector<int> path;

                vector<int> closed_loop = find_closed_loop(G, start_node, end_nodes[j], visited_set, path, visited);
                if (!closed_loop.empty()) {
                    loops.push_back(closed_loop);
                }
            }
        }
    }

    // Remove duplicate/subset loops
    vector<int> to_remove;
    vector<vector<int>> loops_new;

    for (int i = 0; i < static_cast<int>(loops.size()); i++) {
        const auto &list_1 = loops[i];
        if (list_1.empty()) {
            to_remove.push_back(i);
            continue;
        } else if (find(to_remove.begin(), to_remove.end(), i) != to_remove.end()) {
            continue;
        }
        vector<int> double_list_1 = list_1;
        double_list_1.insert(double_list_1.end(), list_1.begin(), list_1.end());
        // Compare current loop with others to check for duplicates/subsets
        for (int j = i + 1; j < static_cast<int>(loops.size()); j++) {
            const auto &list_2 = loops[j];

            if (list_2.empty()) {
                to_remove.push_back(j);
            }

            if (double_list_1.size() < list_2.size()) {  // since if seq.size() - subseq.size()<0, this
                                                         // would give an arbitrary large number and
                                                         // core dumped
                continue;
            }

            else if (isSubsequence(double_list_1, list_2, j)) {
                to_remove.push_back(j);
            }
        }
    }

    // Collect non-duplicate/non-subset loops
    for (size_t i = 0; i < loops.size(); ++i) {
        if (find(to_remove.begin(), to_remove.end(), i) == to_remove.end()) {
            loops_new.push_back(loops[i]);
        }
    }

    return loops_new;
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
 * @param[in] G: graph of activated edges. Mapping between nodes and their
 * neighboring nodes
 * @param[in] start_node: node the path is starting from
 * @param[in] current: current node under investigation in the path
 * @param[in] visited: This includes all visited_all nodes. This can help in
 * skipping the investigation of nodes that are already in the path to avoid
 * processing of loops list in get_closed_loops. This was deactivated in this
 * code, and can be used in the future to reduce processing of detected loops.
 * @param[out] visited_all: visited_all nodes while constructing the path
 * @param[out] path: current detected path, which will be augmented
 * @return new_path or empty vector if no path is found. new_path is a path
 * between the nodes (if a tour is found, )
 *
 */
vector<int> GraphicalModelTsp::find_closed_loop(const unordered_map<int, vector<int>> &G, int start_node, int current, unordered_set<int> visited_all, vector<int> path, vector<int> &visited) {
    visited_all.insert(current);
    path.push_back(current);

    // Check if loop is found
    if (current == start_node) {
        return path;
    }

    // If current node has no outgoing connections
    if (G.find(current) == G.end()) {
        vector<int> new_path = {start_node};
        new_path.insert(new_path.end(), path.begin(), path.end());
        return vector<int>();
    }

    for (int node : G.at(current)) {
        // If node is not visited_all, continue search
        if (visited_all.find(node) == visited_all.end()) {
            vector<int> new_path = find_closed_loop(G, start_node, node, visited_all, path, visited);
            // Return new path if loop is found
            if (!new_path.empty()) {
                return new_path;
            }
        }
    }

    return vector<int>();
}