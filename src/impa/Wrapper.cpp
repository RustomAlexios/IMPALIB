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
        if (model_graph.delta_S_indices_list.size() > 0)
        {
            vector<vector<impalib_type>> temp(model_graph.delta_S_indices_list.size(),
                                              vector<impalib_type>(NUM_EDGE_VARIABLES, zero_value));
            model_graph.subtourConstraints2EdgeEcM.insert(model_graph.subtourConstraints2EdgeEcM.end(), temp.begin(),
                                                          temp.end());
            model_graph.subtourConstraints2EdgeEcDummyM = model_graph.subtourConstraints2EdgeEcM;
        }
    }

    while (!model_graph.subtourConstraintsSatisfiedFlag && AUGMENTATION_FLAG && !model_graph.tourImpaFlag)
    {

        if (model_graph.noConsClosedLoopsCountExcFlag_ || model_graph.noImprovSolCountExcFlag_
            || model_graph.solOscCountExcFlag_)
        {
            cout << "model_graph.noConsClosedLoopsCountExcFlag_: " << model_graph.noConsClosedLoopsCountExcFlag_
                 << '\n';
            cout << "model_graph.noImprovSolCountExcFlag_: " << model_graph.noImprovSolCountExcFlag_ << '\n';
            cout << "model_graph.solOscCountExcFlag_: " << model_graph.solOscCountExcFlag_ << '\n';
            break;
        }
        if (model_graph.costImpa_ == zero_value)
        {
            cout << "Cost is zero" << '\n';
            cout << "Possibly Nan" << '\n';
            exit(0);
        }
        model_graph.iterate_augmented_graph();

        if (model_graph.numAugmentations_ == MAX_AUGM_COUNT)
        {
            cout << "MAX_AUGM_COUNT reached" << '\n';
            break;
        }

        if (model_graph.subtourConstraints2EdgeEcM.size() != model_graph.delta_S_indices_list.size())
        {
            size_t numLists2Add =
                model_graph.delta_S_indices_list.size() - model_graph.subtourConstraints2EdgeEcM.size();
            vector<vector<impalib_type>> temp(numLists2Add, vector<impalib_type>(NUM_EDGE_VARIABLES, zero_value));
            model_graph.subtourConstraints2EdgeEcM.insert(model_graph.subtourConstraints2EdgeEcM.end(), temp.begin(),
                                                          temp.end());
            model_graph.subtourConstraints2EdgeEcDummyM = model_graph.subtourConstraints2EdgeEcM;
        }
    }

    copy(model_graph.outputs.ExtrinsicOutputEdgeEc.begin(),
         model_graph.outputs.ExtrinsicOutputEdgeEc.begin() + NUM_EDGE_VARIABLES, pExtrinsic_output_edge_ec);
    *pNum_augmentations     = model_graph.numAugmentations_;
    *pNum_added_constraints = static_cast<int>(model_graph.subtourConstraints2EdgeEcM.size());
    copy(model_graph.tourImpa_.begin(), model_graph.tourImpa_.begin() + static_cast<int>(model_graph.tourImpa_.size()),
         pTour_impa);
    *pCost_impa                    = model_graph.costImpa_;
    *pNo_improv_sol_count_exc_flag = model_graph.noImprovSolCountExcFlag_;
    *pNo_cons_loops_count_exc_flag = model_graph.noConsClosedLoopsCountExcFlag_;
    *pSol_osc_count_exc_flag       = model_graph.solOscCountExcFlag_;

    vector<int> flattened_selected_edges =
        accumulate(model_graph.selectedEdges_.begin(), model_graph.selectedEdges_.end(), vector<int>{},
                   [](vector<int> &acc, const vector<int> &inner)
                   {
                       acc.insert(acc.end(), inner.begin(), inner.end());
                       return acc;
                   });

    copy(flattened_selected_edges.begin(),
         flattened_selected_edges.begin() + static_cast<int>(flattened_selected_edges.size()), pSelected_edges);
    *pSelected_edges_size = static_cast<int>(flattened_selected_edges.size());

    vector<int> flattened_closed_paths =
        accumulate(model_graph.subtourPaths_.begin(), model_graph.subtourPaths_.end(), vector<int>{},
                   [](vector<int> &acc, const vector<int> &inner)
                   {
                       acc.insert(acc.end(), inner.begin(), inner.end());
                       return acc;
                   });

    copy(flattened_closed_paths.begin(),
         flattened_closed_paths.begin() + static_cast<int>(flattened_closed_paths.size()), pSubtour_paths);
    copy(model_graph.closedPathsSize_.begin(),
         model_graph.closedPathsSize_.begin() + static_cast<int>(model_graph.closedPathsSize_.size()),
         pSubtour_paths_size);
}