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
    int nProjects_;                                 ///< number of projects
    int nTeams_;                                    ///< number of teams
    int nDepartments_;                              ///< number of departments
    int maxSizeNonZeroWeights_;                     ///< maximum size of non-zero weights per department
    int nIter_;                                     ///< number of iterations of IMPA
    bool doFilter_;                                 ///< filtering flag of knapsack constraints
    impalib_type alpha_;                            ///< filtering parameter
    vector<vector<impalib_type>> extrinsicOutPre_;  ///< messages from departments to teams before filtering
    vector<vector<impalib_type>> extrinsicOut_;     ///< messages from departments to teams after filtering
    vector<impalib_type> oric2PackageM_;            ///< messages from ORIC to packages
    vector<vector<impalib_type>> eq2OricM_;         ///< messages from team equality constraint to ORIC
    vector<vector<impalib_type>> oric2EqM_;         ///< messages from ORIC to team equality constraint
    vector<vector<impalib_type>> eq2ProjectM_;      ///< messages from project equality constraint to project inequality constraint
    vector<vector<impalib_type>> project2EqM_;      ///< messages from project inequality constraint to project equality constraint
    vector<impalib_type> team2OricM_;               ///< messages from team equality constraint to ORIC
    Knapsack knapsacks_;                            ///< Knapsack object
    InequalityConstraint ineqs_;                    ///< Project Inequality constraint object
    // EqualityConstraintKcMwm modelEqConstraint_; ///< Equality constraint object
    EqualityConstraint eqs_;        ///< Equality constraint object
    OrInequalityConstraint orics_;  ///< ORIC object

   public:
    OutputsKcMwm outputs;                                                                                                             ///< Graphical model outputs objects
    InputsKcMwm modelInputs_;                                                                                                         ///< Graphical model inputs objects
    void initialize(const impalib_type *, impalib_type *, const int *, const int *, const int *, const impalib_type *, const int *);  ///< initialize graphical model
    void iterate(const int *);                                                                                                        ///< iterate over graphical model
    GraphicalModelKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS, int MAX_SIZE_NON_ZERO_WEIGHTS, int N_ITERATIONS, bool FILT_FLAG,
                        impalib_type ALPHA);  ///< constructor
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
 *
 */

inline GraphicalModelKcMwm::GraphicalModelKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG,
                                                const impalib_type ALPHA)
    : modelInputs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS),
      outputs(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      knapsacks_(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA),
      ineqs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      eqs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      orics_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS),
      nProjects_(N_PROJECTS),
      nTeams_(N_TEAMS),
      nDepartments_(N_DEPARTMENTS),
      maxSizeNonZeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS),
      nIter_(N_ITERATIONS),
      doFilter_(FILT_FLAG),
      alpha_(ALPHA) {
    /// initializes and prepares data structures for mapping relationships
    /// between teams, projects, and constraints
    oric2PackageM_.resize(nTeams_);
    team2OricM_.resize(nTeams_);
    fill(oric2PackageM_.begin(), oric2PackageM_.begin() + nTeams_, zero_value);
    fill(team2OricM_.begin(), team2OricM_.begin() + nTeams_, zero_value);

    eq2OricM_.reserve(nProjects_);
    oric2EqM_.reserve(nProjects_);
    eq2ProjectM_.reserve(nProjects_);
    project2EqM_.reserve(nProjects_);

    for (int project = 0; project < nProjects_; project++) {
        eq2OricM_.push_back(vector<impalib_type>(nTeams_, zero_value));
        oric2EqM_.push_back(vector<impalib_type>(nTeams_, zero_value));
        eq2ProjectM_.push_back(vector<impalib_type>(nTeams_, zero_value));
        project2EqM_.push_back(vector<impalib_type>(nTeams_, zero_value));
    }

    /// initializes and prepares data structures for messages from departments
    extrinsicOut_.reserve(nDepartments_);
    extrinsicOutPre_.reserve(nDepartments_);

    /// initializes and populates nested vectors for mapping relationships
    /// between departments and teams
    for (int department = 0; department < nDepartments_; department++) {
        extrinsicOut_.push_back(vector<impalib_type>(nTeams_, zero_value));
        extrinsicOutPre_.push_back(vector<impalib_type>(nTeams_, zero_value));
    }
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

inline void GraphicalModelKcMwm::initialize(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                                            const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, const int *pMAX_STATE_PY) {
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

inline void GraphicalModelKcMwm::iterate(const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY) {
    for (int iter = 0; iter < nIter_; iter++) {
        for (int department = 0; department < nDepartments_; department++) {
            int capacity = modelInputs_.capacities_[department];

            // Initialize forward and backward message vectors
            vector<vector<impalib_type>> stage_forward_messages(nTeams_ + 1, vector<impalib_type>(capacity + 1, zero_value));
            vector<vector<impalib_type>> stage_backward_messages(nTeams_ + 1, vector<impalib_type>(capacity + 1, zero_value));

            knapsacks_.forward(department, stage_forward_messages, capacity, modelInputs_.nonzero_, pNON_ZERO_WEIGHT_INDICES_SIZES_PY, modelInputs_.weights_, modelInputs_.team2KnapsackM_);
            knapsacks_.backward(department, stage_backward_messages, capacity, modelInputs_.nonzero_, pNON_ZERO_WEIGHT_INDICES_SIZES_PY, modelInputs_.weights_, modelInputs_.team2KnapsackM_);

            knapsacks_.extrinsic_output_department_lhs(modelInputs_.weights_, stage_forward_messages, modelInputs_.team2KnapsackM_, department, stage_backward_messages, capacity, extrinsicOutPre_);
            knapsacks_.process_extrinsic_output_department(department, iter, extrinsicOutPre_, extrinsicOut_);
        }

        eqs_.team_eq_constraint_to_oric_update(extrinsicOut_, team2OricM_, modelInputs_.rewardsTeam_);
        orics_.oric_to_project_eq_constraint_update(eq2OricM_, team2OricM_, oric2EqM_, eq2ProjectM_, modelInputs_.rewardsProj_);
        ineqs_.project_inequality_constraint_update(eq2ProjectM_, project2EqM_);
        eqs_.project_eq_constraint_to_oric_update(project2EqM_, eq2OricM_, modelInputs_.rewardsProj_);

        orics_.oric_to_team_update(eq2OricM_, oric2PackageM_);
        outputs.intrinsic_out_mwm_update(oric2EqM_, project2EqM_, modelInputs_.rewardsProj_);
        knapsacks_.team_to_knapsack_update(modelInputs_.nonzero_, modelInputs_.team2KnapsackM_, modelInputs_.rewardsTeam_, extrinsicOut_, oric2PackageM_);
    }
    outputs.extrinsic_output_team_update(extrinsicOut_, oric2PackageM_);
}

/**
 * Represents a graphical model for the TSP for computing the shortest route
 * between nodes (while visiting each node exactly once).
 */
class GraphicalModelTsp {
   private:
    int nNodes_;                                   ///< number of nodes of TSP
    int nIter_;                                    ///< number of iterations of IMPA
    int nEdgeVars_;                                ///< number of edges in TSP
    bool doReset_;                                 ///< flag for resetting messages after each augmentation step
    bool doFilter_;                                ///< flag for filtering messages from degree/subtour constraints to edge equality constraints
    bool doAugment_;                               ///< whether to activate augmentation
    vector<int> hard_decision;                     ///< hard decision vector of IMPA solution
    impalib_type alpha_;                           ///< filtering parameter
    impalib_type threshold_;                       ///< threshold on hard decision
    EqualityConstraint eqs_;                       ///< Equality Constraint object for TSP
    DegreeConstraint degrees_;                     ///< Degree Constraint object for TSP
    SubtourEliminationConstraint subtours_;        ///< Subtour constraint for TSP
    vector<vector<impalib_type>> deg2EqPreM_;      ///< messages from degree constraint to team equality constraint before filtering
    vector<vector<impalib_type>> deg2EqM_;         ///< messages from degree constraint to team equality constraint afters filtering
    vector<vector<int>> selectedOld_;              ///< old list of selected edges (used for failure investigation)
    vector<vector<int>> selectedOldOld_;           ///< another old list of selected edges (used for failure investigation)
    int maxFailures_;                              ///< maximum count of failures
    bool isTour_ = false;                          ///< set flag for detected tour to be false (will be updated later)
    vector<vector<int>> deltaS_;                   ///< will contain indices of edges that contribute to each subtour constraint
    vector<vector<impalib_type>> subtour2EqPreM_;  ///< messages from subtour constraints to edge equality constraint before filtering
    vector<vector<impalib_type>> eq2SubtourM_;     ///< messages from edge equality constraint to subtour constraints
    void iterate_augmented_graph();                ///< function of IMPA on augmented graph
    void subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &, const vector<vector<int>> &);              ///< analysis of subtour constraints
    void hard_decision_analysis(vector<vector<int>> &);                                                                         ///< function for hard decision solution on IMPA solution
    static bool isSubsequence(const vector<int> &, const vector<int> &, int);                                                   ///< function for post-processing loops
    vector<vector<int>> get_closed_loops(unordered_map<int, vector<int>> &, const vector<vector<int>> &);                       ///< function for getting loops
    vector<int> find_closed_loop(unordered_map<int, vector<int>> &, int, int, unordered_set<int>, vector<int>, vector<int> &);  ///< function for finding loops
    InputsTsp inputs_;                                                                                                          ///< Graphical Model Input object
    vector<vector<int>> selected_;                                                                                              ///< activated edges of IMPA
    int nAugment_;                                                                                                              ///< number of performed augmentations in IMPA
    int nLoop_;                                                                                                                 ///< count for failure case (no consecutive loop detection and no tour)
    int nOsc_;                                                                                                                  ///< count for failure case (oscillation in the solution)
    bool isFailLoop_;                                                                                                           ///< flag for failure case (no consecutive loop detection and no tour)
    int nImprove_;                                                                                                              ///< count for failure case (no solution improvement)
    bool isFailImprove_;                                                                                                        ///< flag for failure case (no solution improvement)
    bool isFailOsc_;                                                                                                            ///< flag for failure case (oscillation in the solution)
    vector<int> tour_;                                                                                                          ///< list of nodes of tour (if detected)
    vector<vector<int>> subtourPaths_;                                                                                          ///< list of detected subtours (if detected)
    vector<int> loopsSize_;                                                                                                     ///< list of sizes of loops
    impalib_type tourCost_;                                                                                                     ///< cost of IMPA solution
    vector<vector<impalib_type>> subtour2EqM_;                                                                                  ///< messages from subtour constraints to edge equality constraint

   public:
    bool subtourConstraintsSatisfiedFlag = false;                                                                    ///< initially set this flag for satisfied subtour constraints to false
    OutputsTsp outputs;                                                                                              ///< TSP graphical model outputs object
    void initialize(const int *, const impalib_type *, const impalib_type *, impalib_type *, const impalib_type *);  ///< initialize graphical model
    void iterate_relaxed_graph();                                                                                    ///< iterate over relaxed graphical model
    void perform_augmentation(int);                                                                                  ///< perform augmentation on graphical model
    void process_ouputs(impalib_type *, int *, int *, int *, impalib_type *, bool *, bool *, bool *, int *, int *, int *, int *);  ///< process outputs of Graphical Model
    GraphicalModelTsp(int NUM_ITERATIONS, int NUM_NODES, int NUM_EDGE_VARIABLES, bool AUGMENTATION_FLAG, bool RESET_FLAG, bool FILTERING_FLAG, impalib_type ALPHA, impalib_type THRESHOLD,
                      int MAX_COUNT);  ///< Constructor
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
 *
 */

inline GraphicalModelTsp::GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG,
                                            const impalib_type ALPHA, const impalib_type THRESHOLD, const int MAX_COUNT)
    : degrees_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      eqs_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      subtours_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
      inputs_(NUM_NODES, NUM_EDGE_VARIABLES),
      outputs(NUM_NODES, NUM_EDGE_VARIABLES),
      nIter_(NUM_ITERATIONS),
      doFilter_(FILTERING_FLAG),
      alpha_(ALPHA),
      nNodes_(NUM_NODES),
      doAugment_(AUGMENTATION_FLAG),
      doReset_(RESET_FLAG),
      nEdgeVars_(NUM_EDGE_VARIABLES),
      threshold_(THRESHOLD),
      maxFailures_(MAX_COUNT) {
    deg2EqPreM_.reserve(nEdgeVars_);
    deg2EqM_.reserve(nEdgeVars_);

    // Populate degree constraint to equality constraint messages
    for (int edge = 0; edge < nEdgeVars_; edge++) {
        deg2EqPreM_.push_back(vector<impalib_type>(nNodes_, zero_value));
        deg2EqM_.push_back(vector<impalib_type>(nNodes_, zero_value));
    }

    nAugment_ = 0;
    tourCost_ = zero_value;
    nLoop_ = 0;
    nImprove_ = 0;
    nOsc_ = 0;

    isFailLoop_ = false;
    isFailImprove_ = false;
    isFailOsc_ = false;
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

inline void GraphicalModelTsp::initialize(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY,
                                          impalib_type *pEdge_ec_to_degree_constraint_m_py, const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY) {
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
inline void GraphicalModelTsp::iterate_relaxed_graph() {
    for (int iter = 0; iter < nIter_; iter++) {
        degrees_.degree_constraint_to_edge_ec_update(inputs_.eq2DegM_, inputs_.connections_, deg2EqPreM_);
        degrees_.process_filtering(iter, deg2EqPreM_, deg2EqM_);
        eqs_.edge_ec_to_degree_constraint_relaxed_graph_update(inputs_.connections_, inputs_.costMatUpdate_, deg2EqM_, inputs_.eq2DegM_);
    }
    outputs.extrinsic_output_edge_ec_relaxed_graph_update(deg2EqM_);
    outputs.intrinsic_output_edge_ec_update(inputs_.costs_);

    // Perform hard decision analysis to select edges
    vector<vector<int>> selected_edges;
    hard_decision_analysis(selected_edges);

    // Initialize graph for subtour elimination constraints analysis
    unordered_map<int, vector<int>> graph;
    subtour_elimination_constraints_analysis(graph, selected_edges);

    selected_ = selected_edges;
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

inline void GraphicalModelTsp::perform_augmentation(const int MAX_AUGM_COUNT) {
    // Add empty vectors to subtour constraints if needed
    if (deltaS_.size() > 0) {
        vector<vector<impalib_type>> temp(deltaS_.size(), vector<impalib_type>(nEdgeVars_, zero_value));
        subtour2EqM_.insert(subtour2EqM_.end(), temp.begin(), temp.end());
        subtour2EqPreM_ = subtour2EqM_;
    }

    // Continue augmentation
    while (!subtourConstraintsSatisfiedFlag && doAugment_ && !isTour_) {
        if (isFailLoop_ || isFailImprove_ || isFailOsc_) {
            cout << "noConsClosedLoopsCountExcFlag_: " << isFailLoop_ << '\n';
            cout << "noImprovSolCountExcFlag_: " << isFailImprove_ << '\n';
            cout << "solOscCountExcFlag_: " << isFailOsc_ << '\n';
            break;
        }
        if (tourCost_ == zero_value) {
            cout << "Cost is zero" << '\n';
            cout << "Possibly Nan" << '\n';
            exit(0);
        }

        iterate_augmented_graph();

        if (nAugment_ == MAX_AUGM_COUNT) {
            cout << "MAX_AUGM_COUNT reached" << '\n';
            break;
        }

        // Add empty vectors to subtour constraints if needed
        if (subtour2EqM_.size() != deltaS_.size()) {
            size_t numLists2Add = deltaS_.size() - subtour2EqM_.size();
            vector<vector<impalib_type>> temp(numLists2Add, vector<impalib_type>(nEdgeVars_, zero_value));
            subtour2EqM_.insert(subtour2EqM_.end(), temp.begin(), temp.end());
            subtour2EqPreM_ = subtour2EqM_;
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

inline void GraphicalModelTsp::process_ouputs(impalib_type *pExtrinsic_output_edge_ec, int *pNum_augmentations, int *pNum_added_constraints, int *pTour_impa, impalib_type *pCost_impa,
                                              bool *pNo_improv_sol_count_exc_flag, bool *pNo_cons_loops_count_exc_flag, bool *pSol_osc_count_exc_flag, int *pSelected_edges, int *pSelected_edges_size,
                                              int *pSubtour_paths, int *pSubtour_paths_size) {
    // Copy extrinsic output for edge equality constraints
    copy(outputs.extrinsicOut_.begin(), outputs.extrinsicOut_.begin() + nEdgeVars_, pExtrinsic_output_edge_ec);
    *pNum_augmentations = nAugment_;
    *pNum_added_constraints = static_cast<int>(subtour2EqM_.size());

    copy(tour_.begin(), tour_.begin() + static_cast<int>(tour_.size()), pTour_impa);
    *pCost_impa = tourCost_;

    *pNo_improv_sol_count_exc_flag = isFailImprove_;
    *pNo_cons_loops_count_exc_flag = isFailLoop_;
    *pSol_osc_count_exc_flag = isFailOsc_;

    // Copy selected edges
    vector<int> flattened_selected_edges = accumulate(selected_.begin(), selected_.end(), vector<int>{}, [](vector<int> &acc, const vector<int> &inner) {
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
    copy(loopsSize_.begin(), loopsSize_.begin() + static_cast<int>(loopsSize_.size()), pSubtour_paths_size);
}

/**
 * This function will run IMPA on augmented grahical model
 *
 */
inline void GraphicalModelTsp::iterate_augmented_graph() {
    // Reset messages if flag is set
    if (doReset_) {
        for (size_t i = 0; i < inputs_.connections_.size(); i++) {
            auto connection = inputs_.connections_[i];
            impalib_type cost = inputs_.costMat_[connection[0]][connection[1]];
            inputs_.eq2DegM_[i][connection[0]] = cost;
            inputs_.eq2DegM_[i][connection[1]] = cost;
        }
    }

    cout << "-----------------------" << '\n';
    cout << "iterate_augmented_graph" << '\n';
    cout << "-----------------------" << '\n';

    // Increment augmentation count
    nAugment_ += 1;
    cout << "Augmentation count: " << nAugment_ << '\n';
    cout << "delta_S_indices_list.size(): " << deltaS_.size() << '\n';

    // Reserve memory for old subtour constraints
    subtours_.subtour2EqOldM_.reserve(deltaS_.size());

    // Initialize old subtour constraints
    for (int subtour = 0; subtour < deltaS_.size(); subtour++) {
        subtours_.subtour2EqOldM_.push_back(vector<impalib_type>(nEdgeVars_, zero_value));
    }

    for (int iter = 0; iter < nIter_; iter++) {
        degrees_.degree_constraint_to_edge_ec_update(inputs_.eq2DegM_, inputs_.connections_, deg2EqPreM_);
        degrees_.process_filtering(iter, deg2EqPreM_, deg2EqM_);

        eq2SubtourM_ = eqs_.edge_ec_to_subtour_constraints_update(deltaS_, inputs_.costs_, deg2EqM_, subtour2EqM_, inputs_.connections_);
        subtours_.subtour_constraints_to_edge_ec_update(eq2SubtourM_, deltaS_, subtour2EqPreM_);
        subtours_.process_filtering(iter, subtour2EqPreM_, subtour2EqM_, deltaS_);
        eqs_.edge_ec_to_degree_constraint_augmented_graph_update(deg2EqM_, subtour2EqM_, inputs_.connections_, inputs_.costMatUpdate_, inputs_.eq2DegM_);
    }

    // Check for NaN values
    bool flag_nan = false;
    for (int i = 0; i < inputs_.eq2DegM_.size(); i++) {
        for (int j = 0; j < inputs_.eq2DegM_[i].size(); j++) {
            if (isnan(inputs_.eq2DegM_[i][j])) {
                flag_nan = true;
            }
        }
    }

    if (flag_nan) {
        cout << "NaN detected" << '\n';
        exit(0);
    }

    outputs.extrinsic_output_edge_ec_augmented_graph_update(deg2EqM_, subtour2EqM_);
    outputs.intrinsic_output_edge_ec_update(inputs_.costs_);

    selected_.clear();
    vector<vector<int>> selected_edges;
    hard_decision_analysis(selected_edges);

    // Check for solution improvement
    if (!selectedOld_.empty()) {
        bool areEqual = (selected_edges == selectedOld_);
        if (areEqual) {
            cout << "Exited: No Improvement of IMPA Solution" << '\n';
            nImprove_ += 1;
            if (nImprove_ > maxFailures_) {
                cout << "Exited: noImprovSolCount_ Exceeded maximum allowable "
                        "number "
                     << maxFailures_ << '\n';
                isFailImprove_ = true;
            }
        } else {
            nImprove_ = 0;
        }
    }

    // Check for solution oscillation
    if (!selectedOldOld_.empty()) {
        bool oscillation_flag = ((selected_edges == selectedOldOld_) && (selectedOldOld_ != selectedOld_));
        if (oscillation_flag) {
            cout << "Exited: Solution is oscillating" << '\n';
            nOsc_ += 1;
            if (nOsc_ > maxFailures_) {
                cout << "Exited: solOscCount_ Exceeded maximum allowable number " << maxFailures_ << '\n';
                isFailOsc_ = true;
            }
        } else {
            nOsc_ = 0;
        }
    }

    if (!selectedOld_.empty()) {
        selectedOldOld_ = selectedOld_;
    }
    selectedOld_ = selected_edges;
    selected_ = selected_edges;
    unordered_map<int, vector<int>> graph;
    subtour_elimination_constraints_analysis(graph, selected_edges);
}

/**
 * This will get the hard decision on edges by investigating the sign
 * IntrinsicOutputEdgeEc
 * @param[out] selected: activated edges after running IMPA
 *
 */
inline void GraphicalModelTsp::hard_decision_analysis(vector<vector<int>> &selected) {
    hard_decision.resize(nEdgeVars_);
    fill(hard_decision.begin(), hard_decision.begin() + nEdgeVars_, numeric_limits<int>::max());

    // Apply threshold to intrinsic output to determine hard decision
    transform(outputs.intrinsicOut_.begin(), outputs.intrinsicOut_.end(), hard_decision.begin(), [&](impalib_type value) { return value > threshold_ ? zero_value : 1; });

    // Select edges based on hard decision
    for (int edge = 0; edge < nEdgeVars_; edge++) {
        if (hard_decision[edge] == 1) {
            selected.push_back(inputs_.connections_[edge]);
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
            cost_impa += inputs_.costs_[i];
        }
    }

    cout << "C++ cost_impa: " << cost_impa << '\n';
    tourCost_ = cost_impa;
}

/**
 * Analyze the activated edges. If tour found, return tour. If tour not found,
 * detect subtours. Count failure cases if present
 * @param[in] selected: activated edges in the graphical model
 * @param[out] G: graphical model of activated egdes. Defines connections
 * between nodes
 *
 */

inline void GraphicalModelTsp::subtour_elimination_constraints_analysis(unordered_map<int, vector<int>> &G, const vector<vector<int>> &selected) {
    loopsSize_.clear();  // just store it at the end if applicable

    // Get closed loops from the graph
    vector<vector<int>> loops_list = get_closed_loops(G, selected);
    subtourPaths_.clear();
    subtourPaths_ = loops_list;

    if (loops_list.empty()) {
        cout << "Exited: Cannot find tours using get_closed_loops()" << '\n';
        nLoop_ += 1;
        if (nLoop_ > maxFailures_) {
            cout << "Exited: noConsClosedLoopsCount_ Exceeded maximum allowable "
                    "number: "
                 << maxFailures_ << '\n';
            isFailLoop_ = true;
        }
    }

    else if (loops_list.size() == 1 && loops_list[0].size() == nNodes_) {
        isTour_ = true;
        subtourConstraintsSatisfiedFlag = true;
        tour_ = loops_list[0];
        tour_.push_back(loops_list[0][0]);
        cout << "Tour found" << '\n';

        cout << "tour_impa: [";
        for (size_t i = 0; i < tour_.size(); ++i) {
            cout << tour_[i];
            if (i != tour_.size() - 1) {
                cout << ", ";
            }
        }
        cout << "]" << '\n';
    }
    // If multiple subtours are found
    else {
        nLoop_ = 0;
        for (const auto &loop : loops_list) {
            // Store the size of each subtour
            loopsSize_.push_back(static_cast<int>(loop.size()));
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
            for (size_t i = 0; i < inputs_.connections_.size(); ++i) {
                const auto &connection = inputs_.connections_[i];
                if (find(loop.begin(), loop.end(), connection[0]) != loop.end() && find(loop.begin(), loop.end(), connection[1]) == loop.end()) {
                    delta_S_indices.push_back(static_cast<int>(i));
                }
            }
            // Add delta S indices to delta_S_indices_list
            if (deltaS_.empty() && !delta_S_indices.empty()) {  // added !delta_S_indices.empty() to make
                                                                // sure if a full tour was detected so it
                                                                // is not be added with subtours
                deltaS_.push_back(delta_S_indices);
            } else {
                // Add delta S indices if not already present
                bool found = false;
                for (const auto &existing_indices : deltaS_) {
                    if (existing_indices == delta_S_indices) {
                        found = true;
                        break;
                    }
                }
                if (!found && !delta_S_indices.empty()) {  // added !delta_S_indices.empty() to
                                                           // make sure if a full tour was detected
                                                           // so it is not be added with subtours
                    deltaS_.push_back(delta_S_indices);
                }
            }
        }
    }
}

/**
 * Get closed loops. Function for building the graphical model of activated
 * edges and obtain loops
 * @param[in] selected: activated edges in the graphical model
 * @param[out] G: graphical model of activated egdes. Defines connections
 * between nodes. Will be used for detecting of subtours
 * @return new_loops_list: list of detected loops in rGraph
 *
 */
inline vector<vector<int>> GraphicalModelTsp::get_closed_loops(unordered_map<int, vector<int>> &G, const vector<vector<int>> &selected) {
    // Update the graph based on selected edges
    for (const auto &connection : selected) {
        if (G.find(connection[0]) != G.end()) {
            G[connection[0]].push_back(connection[1]);
        } else {
            G[connection[0]] = {connection[1]};
        }
    }

    vector<vector<int>> loops_list;
    vector<int> visited_nodes;

    for (const auto &connection : G) {
        int start_node = connection.first;
        const vector<int> &end_nodes = connection.second;
        // Skip if node is already visited
        if (find(visited_nodes.begin(), visited_nodes.end(), start_node) != visited_nodes.end()) {
            continue;
        } else {
            for (int j = 0; j < end_nodes.size(); j++) {
                unordered_set<int> visited_set;
                vector<int> path;

                vector<int> closed_loop = find_closed_loop(G, start_node, end_nodes[j], visited_set, path, visited_nodes);
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
inline bool GraphicalModelTsp::isSubsequence(const vector<int> &seq, const vector<int> &subseq, int j) {
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
 * @param[in] start: node the path is starting from
 * @param[in] cur: current node under investigation in the path
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
inline vector<int> GraphicalModelTsp::find_closed_loop(unordered_map<int, vector<int>> &G, const int start, const int cur, unordered_set<int> visited, vector<int> path, vector<int> &visited_nodes) {
    visited.insert(cur);
    path.push_back(cur);

    // Check if loop is found
    if (cur == start) {
        return path;
    }

    // If current node has no outgoing connections
    if (G.find(cur) == G.end()) {
        vector<int> new_path = {start};
        new_path.insert(new_path.end(), path.begin(), path.end());
        return vector<int>();
    }

    for (int node : G.at(cur)) {
        // If node is not visited, continue search
        if (visited.find(node) == visited.end()) {
            vector<int> new_path = find_closed_loop(G, start, node, visited, path, visited_nodes);
            // Return new path if loop is found
            if (!new_path.empty()) {
                return new_path;
            }
        }
    }

    return vector<int>();
}

/**
 * Represents a graphical model for the k-sat problem for
 * computing the assignment of variables to satisfy constraints.
 */
class GraphicalModelKsat {
   private:
    int nIter_;               ///< number of iterations
    int nVar_;                ///< total number of variables
    int nConstraints_;        ///< number of constraints
    int k_;                   ///< number of variables per constraint
    bool filteringFlag_;      ///< whether to apply filtering or not
    impalib_type alpha_;      ///< filtering parameter
    int nUsed_;               ///< number of used variables to build a formula
    EqualityConstraint eqs_;  ///< Equality constraint object
    KsatConstraint ksats_;    ///< ksat constraint object

   public:
    InputsKsat inputs_;                         ///< input object
    vector<vector<impalib_type>> ksat2EqPreM_;  ///< messages from ksat constraints to equality constraints before filtering
    vector<vector<impalib_type>> ksat2EqM_;     ///< messages from ksat constraints to equality constraints after filtering
    OutputsKsat outputs_;                       ///< output object
    void initialize(const int *, const int *, const int *, const int *, const int *, const impalib_type *, impalib_type *);  ///< initialize graphical model
    void iterate();                                                                                                          ///< iterate function
    GraphicalModelKsat(int NUM_ITERATIONS, int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE, bool FILTERING_FLAG, impalib_type ALPHA,
                       int NUM_USED_VARIABLES);  ///< constructor for k-sat graphical model
};

/**
 * Construct Graphical Model for the K-SAT problem
 *
 * @param[in] NUM_ITERATIONS: number of iterations
 * @param[in] NUM_VARIABLES: total number of variables
 * @param[in] NUM_CONSTRAINTS: number of constraints
 * @param[in] K_VARIABLE: number of variables per constraint
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[in] NUM_USED_VARIABLES: number of variables used to construct the formula
 *
 */

inline GraphicalModelKsat::GraphicalModelKsat(const int NUM_ITERATIONS, const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG, const impalib_type ALPHA,
                                              const int NUM_USED_VARIABLES)
    : inputs_(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, nUsed_),
      eqs_(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILTERING_FLAG, ALPHA),
      ksats_(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILTERING_FLAG, ALPHA),
      outputs_(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE),
      nIter_(NUM_ITERATIONS),
      nVar_(NUM_VARIABLES),
      nConstraints_(NUM_CONSTRAINTS),
      k_(K_VARIABLE),
      filteringFlag_(FILTERING_FLAG),
      alpha_(ALPHA),
      nUsed_(NUM_USED_VARIABLES),
      ksat2EqPreM_(nConstraints_, vector<impalib_type>(nVar_, zero_value)),
      ksat2EqM_(nConstraints_, vector<impalib_type>(nVar_, zero_value)){};

/**
 * Initialize Graphical Model inputs for the k-sat problem. The inputs
 * are processed by reading inputs from python
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

inline void GraphicalModelKsat::initialize(const int *pUSED_VARIABLES_PY, const int *pVARIABLES_CONNECTIONS_PY, const int *pVARIABLES_CONNECTIONS_SIZES, const int *pCONSTRAINTS_CONNECTIONS,
                                           const int *pCONSTRAINTS_CONNECTIONS_TYPE, const impalib_type *pINCOMING_METRICS_COST, impalib_type *pVariable_ec_to_ksat_constraint_m_py) {
    /// calls a method process_inputs() on an object modelInputs_, passing
    /// several pointers to Python objects as arguments
    inputs_.process_inputs(pUSED_VARIABLES_PY, pVARIABLES_CONNECTIONS_PY, pVARIABLES_CONNECTIONS_SIZES, pCONSTRAINTS_CONNECTIONS, pCONSTRAINTS_CONNECTIONS_TYPE, pINCOMING_METRICS_COST,
                           pVariable_ec_to_ksat_constraint_m_py);
}

/**
 * Run IMPA on k-sat-MWM graphical model. This will propagate messages for a
 * certain number of iterations across the whole graphical model
 *
 */

inline void GraphicalModelKsat::iterate() {
    for (int iter = 0; iter < nIter_; iter++) {
        /// update messages from k-sat constraints to variable equality constraints
        ksats_.ksat_constraint_to_variable_ec_update(inputs_.varEq2KsatM_, ksat2EqPreM_, inputs_.connectionsToVars_, inputs_.types_);
        /// perform filtering on messages from k-sat constraints to variable equality constraints
        ksats_.process_filtering(iter, ksat2EqPreM_, ksat2EqM_);
        /// update messages from variable equality constraints to k-sat constraints
        eqs_.variable_ec_to_ksat_constraint_update(ksat2EqM_, inputs_.varEq2KsatM_, inputs_.used_, inputs_.metrics_, inputs_.connectionsFromVars_);
    }
    /// process outputs of k-sat problem
    outputs_.update_extrinsic(ksat2EqM_);
}