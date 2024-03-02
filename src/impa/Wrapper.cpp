// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"

/**
 * This function will run IMPA on Knapsack-MWM graphical model in C++
 *
 * @param[in] NUM_ITERATIONS: number of iterations
 * @param[in] NUM_DEPARTMENTS: number of departments
 * @param[in] NUM_TEAMS: number of teams
 * @param[in] pTransition_model_py: messages from teams to departments
 * @param[in] pMAX_STATE_PY: contains capacities of departments
 * @param[in] pREWARD_TEAM_PY: contains rewards of teams
 * @param[in] pREWARD_PROJECT_PY: contains rewards of projects-teams combinations
 * @param[in] pTEAMS_WEIGHTS_PER_DEPARTMENT_PY: contains weights of teams to each department
 * @param[in] p_NON_ZERO_WEIGHT_INDICES_PY: contains indices of connections of teams to departments
 * @param[in] pNON_ZERO_WEIGHT_INDICES_SIZES_PY: contains sizes of non-zero weight indices of teams to departments
 * @param[in] MAX_SIZE_NON_ZERO_WEIGHTS: contains maximum size of non-zero weights or connections between teams and departments
 * @param[out] pExtrinsic_output_team: pointer that will store the extrinsic output messages of teams
 * @param[out] pIntrinsic_out_mwm: pointer that will store the extrinsic output messages of projects-teams combinations
 * @param[in] FILTERING_FLAG: filtering flag
 * @param[in] ALPHA: filtering parameter (between 0 and 1)
 * 
 */

extern "C" void WrapperKcMwm(const int NUM_ITERATIONS, const int NUM_DEPARTMENTS, const int NUM_TEAMS,
                             const int NUM_PROJECTS, impalib_type *pTransition_model_py, const int *pMAX_STATE_PY,
                             const impalib_type *pREWARD_TEAM_PY, const impalib_type *pREWARD_PROJECT_PY,
                             const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY,
                             const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int MAX_SIZE_NON_ZERO_WEIGHTS,
                             impalib_type *pExtrinsic_output_team, impalib_type *pIntrinsic_out_mwm,
                             const bool FILTERING_FLAG, const impalib_type ALPHA)
{

    // Instantiate a GraphicalModelKcMwm object
    GraphicalModelKcMwm model_graph(NUM_DEPARTMENTS, NUM_TEAMS, NUM_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS, NUM_ITERATIONS,
                                    FILTERING_FLAG, ALPHA);

    // Initialize the graphical model with provided parameters
    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                           pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY,
                           pMAX_STATE_PY);

    // Iterate over the graphical model
    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    // Copy the extrinsic output for each team to the provided variable
    copy(model_graph.outputs.extrinsic.begin(), model_graph.outputs.extrinsic.begin() + NUM_TEAMS,
         pExtrinsic_output_team);
    
    // Copy the intrinsic output for each team-project pair to the provided variable
    copy(model_graph.outputs.intrinsic.begin(),
         model_graph.outputs.intrinsic.begin() + NUM_TEAMS * NUM_PROJECTS, pIntrinsic_out_mwm);
}

/**
 * This function will run IMPA on TSP graphical model in C++
 *
 * @param[in] NUM_ITERATIONS: number of iterations
 * @param[in] NUM_NODES: number of nodes of TSP
 * @param[in] NUM_EDGE_VARIABLES: number of edge connections in the graph
 * @param[in] AUGMENTATION_FLAG: whether the user specifies augmentation or not
 * @param[in] RESET_FLAG: whether messages after each augmentation step are reset to zero or not
 * @param[in] FILTERING_FLAG: whether filtering on constraints is activated or not
 * @param[in] ALPHA: filtering parameter (between 0 and 1)
 * @param[in] pEDGE_CONNECTIONS_PY: edge connections (each element corresponds to an edge and has its constitudent nodes)
 * @param[in] pCOST_EDGE_VARIABLE_PY: cost of each edge connection
 * @param[in] pCOST_MATRIX_PY: cost matrix (originally has size number of nodes x number of nodes)
 * @param[in] pEdge_ec_to_degree_constraint_m_py: initial messages from edge equality constraints to degree constraints
 * @param[in] pEDGE_DEGREE_CONSTRAINT_COST_PY: cost matrix stored in a way function of number of nodes and number of edges to facilitate message update computations
 * @param[out] pExtrinsic_output_edge_ec: extrinsic output messages from edge equality constraints
 * @param[in] THRESHOLD: threshold for hard decision on intrinsic messages (usually set to be negative and close to zero)
 * @param[out] pNum_augmentations: number of augmentation steps performed by IMPA
 * @param[out] pTour_impa: elements of tour (if found) after running IMPA
 * @param[out] pCost_impa: sum of costs of the activated edges after running IMPA
 * @param[out] pNo_improv_sol_count_exc_flag: will be set to true if failure occurs due to exceeding count of no improvement of the solution after running IMPA
 * @param[out] pNo_cons_loops_count_exc_flag: will be set to true if failure occurs due to exceeding count of no consecutive detection of loops (not tour) after running IMPA
 * @param[out] pSelected_edges: activated edges after running IMPA
 * @param[out] pSelected_edges_size: number of activated edges after running IMPA
 * @param[out] pSol_osc_count_exc_flag: will be set to true if failure occurs due to exceeding count of oscillations in the solution after running IMPA
 * @param[out] pSubtour_paths: detected subtours after running IMPA
 * @param[out] pSubtour_paths_size: number of detected subtours after running IMPA
 * @param[in] MAX_COUNT: maximum counts of allowed failure
 * @param[out] pNum_added_constraints: number of added subtour elimination constraints after running IMPA
 * @param[in] MAX_AUGM_COUNT: maximum number of augmentation before stopping IMPA
 * 
 */

extern "C" void WrapperTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                           const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG,
                           const impalib_type ALPHA, const int *pEDGE_CONNECTIONS_PY,
                           const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY,
                           impalib_type       *pEdge_ec_to_degree_constraint_m_py,
                           const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY, impalib_type *pExtrinsic_output_edge_ec,
                           const impalib_type THRESHOLD, int *pNum_augmentations, int *pTour_impa,
                           impalib_type *pCost_impa, bool *pNo_improv_sol_count_exc_flag,
                           bool *pNo_cons_loops_count_exc_flag, int *pSelected_edges, int *pSelected_edges_size,
                           bool *pSol_osc_count_exc_flag, int *pSubtour_paths, int *pSubtour_paths_size,
                           const int MAX_COUNT, int *pNum_added_constraints, const int MAX_AUGM_COUNT)
{
    // Instantiate a GraphicalModelTsp object
    GraphicalModelTsp model_graph(NUM_ITERATIONS, NUM_NODES, NUM_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG,
                                  FILTERING_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

    // Initialize the graphical model with provided data
    model_graph.initialize(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY,
                           pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);

    // Perform IMPA on the relaxed graph
    auto out = model_graph.iterate_relaxed_graph();

    // Check and perform augmentation if needed
    if (!out.subtourConstraintsSatisfiedFlag && AUGMENTATION_FLAG)
    {
        model_graph.perform_augmentation(MAX_AUGM_COUNT);

    }

    // Process the outputs of the model
    model_graph.process_outputs(pExtrinsic_output_edge_ec, pNum_augmentations, pNum_added_constraints, pTour_impa, pCost_impa, pNo_improv_sol_count_exc_flag, \
                              pNo_cons_loops_count_exc_flag, pSol_osc_count_exc_flag, pSelected_edges, pSelected_edges_size, pSubtour_paths, pSubtour_paths_size);
}