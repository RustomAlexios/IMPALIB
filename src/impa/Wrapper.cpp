// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"

extern "C" void WrapperKcMwm(const int NUM_ITERATIONS, const int NUM_DEPARTMENTS, const int NUM_TEAMS,
                             const int NUM_PROJECTS, impalib_type *pTransition_model_py, const int *pMAX_STATE_PY,
                             const impalib_type *pREWARD_TEAM_PY, const impalib_type *pREWARD_PROJECT_PY,
                             const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY,
                             const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int MAX_SIZE_NON_ZERO_WEIGHTS,
                             impalib_type *pExtrinsic_output_team, impalib_type *pIntrinsic_out_mwm,
                             const bool FILTERING_FLAG, const impalib_type ALPHA)
{

    GraphicalModelKcMwm model_graph(NUM_DEPARTMENTS, NUM_TEAMS, NUM_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS, NUM_ITERATIONS,
                                    FILTERING_FLAG, ALPHA);

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                           pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY,
                           pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    copy(model_graph.outputs.ExtrinsicOutputTeam.begin(), model_graph.outputs.ExtrinsicOutputTeam.begin() + NUM_TEAMS,
         pExtrinsic_output_team);
    copy(model_graph.outputs.IntrinsicOutMwm.begin(),
         model_graph.outputs.IntrinsicOutMwm.begin() + NUM_TEAMS * NUM_PROJECTS, pIntrinsic_out_mwm);
}

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

    GraphicalModelTsp model_graph(NUM_ITERATIONS, NUM_NODES, NUM_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG,
                                  FILTERING_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

    model_graph.initialize(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY,
                           pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);

    model_graph.iterate_relaxed_graph();

    if (!model_graph.subtourConstraintsSatisfiedFlag && AUGMENTATION_FLAG)
    {
        model_graph.perform_augmentation(MAX_AUGM_COUNT);

    }

    model_graph.process_ouputs(pExtrinsic_output_edge_ec, pNum_augmentations, pNum_added_constraints, pTour_impa, pCost_impa, pNo_improv_sol_count_exc_flag, \
                              pNo_cons_loops_count_exc_flag, pSol_osc_count_exc_flag, pSelected_edges, pSelected_edges_size, pSubtour_paths, pSubtour_paths_size);
}