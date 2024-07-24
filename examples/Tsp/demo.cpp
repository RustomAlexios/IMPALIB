// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"

int main() {
    const int N_NODES = 5;                                     ///< number of nodes
    const bool FILT_FLAG = true;                               ///< whether filtering is activated or not
    const int N_ITER = 200;                                    ///< number of iterations of IMPA
    const int N_EDGE_VARIABLES = N_NODES * N_NODES - N_NODES;  ///< number of edge variables
    const impalib_type ALPHA = 0.5;                            ///< filtering parameter
    const bool SYM_FLAG = true;                                ///< symmetric flag
    const bool RESET_FLAG = false;                             ///< reset flag
    const impalib_type THRESHOLD = -0.0001;                    ///< threshold parameter
    const bool AUGM_FLAG = true;                               ///< augmentation flag
    const int MAX_AUGM_COUNT = 50;                             ///< maximum augmentation count
    const int MAX_FAILURE_COUNT = 50;                          ///< maximum count for failure

    // connections between nodes for each edge variable
    array<array<int, 2>, N_EDGE_VARIABLES> edge_connections = {
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 0}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 3}, {2, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 4}, {4, 0}, {4, 1}, {4, 2}, {4, 3}}};

    // cost matrix of the tsp problem
    array<array<impalib_type, N_NODES>, N_NODES> cost_matrix = {{{0.0, 151.75578773, 56.18610887, 718.31915651, 293.02503715},
                                                                 {151.75578773, 0.0, 568.4231286, 83.25740946, 545.45357536},
                                                                 {56.18610887, 568.4231286, 0.0, 445.55107005, 888.09445172},
                                                                 {718.31915651, 83.25740946, 445.55107005, 0.0, 719.05730714},
                                                                 {293.02503715, 545.45357536, 888.09445172, 719.05730714, 0.0}}};

    // cost_edge_variable constructed from edge_connections and cost_matrix
    array<impalib_type, N_EDGE_VARIABLES> cost_edge_variable = {151.75578773, 56.18610887,  718.31915651, 293.02503715, 151.75578773, 568.4231286,  83.25740946,
                                                                545.45357536, 56.18610887,  568.4231286,  445.55107005, 888.09445172, 718.31915651, 83.25740946,
                                                                445.55107005, 719.05730714, 293.02503715, 545.45357536, 888.09445172, 719.05730714};

    // edge_degree_constraint_cost constructed from edge_connections and cost_matrix
    array<array<impalib_type, N_NODES>, N_EDGE_VARIABLES> edge_degree_constraint_cost = {
        {{151.75578773, 151.75578773, 0.0, 0.0, 0.0}, {56.18610887, 0.0, 56.18610887, 0.0, 0.0},   {718.31915651, 0.0, 0.0, 718.31915651, 0.0}, {293.02503715, 0.0, 0.0, 0.0, 293.02503715},
         {151.75578773, 151.75578773, 0.0, 0.0, 0.0}, {0.0, 568.4231286, 568.4231286, 0.0, 0.0},   {0.0, 83.25740946, 0.0, 83.25740946, 0.0},   {0.0, 545.45357536, 0.0, 0.0, 545.45357536},
         {56.18610887, 0.0, 56.18610887, 0.0, 0.0},   {0.0, 568.4231286, 568.4231286, 0.0, 0.0},   {0.0, 0.0, 445.55107005, 445.55107005, 0.0}, {0.0, 0.0, 888.09445172, 0.0, 888.09445172},
         {718.31915651, 0.0, 0.0, 718.31915651, 0.0}, {0.0, 83.25740946, 0.0, 83.25740946, 0.0},   {0.0, 0.0, 445.55107005, 445.55107005, 0.0}, {0.0, 0.0, 0.0, 719.05730714, 719.05730714},
         {293.02503715, 0.0, 0.0, 0.0, 293.02503715}, {0.0, 545.45357536, 0.0, 0.0, 545.45357536}, {0.0, 0.0, 888.09445172, 0.0, 888.09445172}, {0.0, 0.0, 0.0, 719.05730714, 719.05730714}}};

    // Extract data pointers from various structures
    const int *pEDGE_CONNECTIONS_PY = addressof(get<0>(edge_connections[0]));
    const impalib_type *pCOST_MATRIX_PY = addressof(get<0>(cost_matrix[0]));
    const impalib_type *pCOST_EDGE_VARIABLE_PY = cost_edge_variable.data();
    impalib_type *pEdge_ec_to_degree_constraint_m_py = addressof(get<0>(edge_degree_constraint_cost[0]));
    const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY = addressof(get<0>(edge_degree_constraint_cost[0]));

    GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGM_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD, MAX_FAILURE_COUNT);

    // Initialize the model with provided data pointers
    model_graph.initialize(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY, pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);

    // Iterate through the relaxed graph
    model_graph.iterate_relaxed_graph();

    // Check if subtour constraints are not satisfied and augmentation is enabled
    if (!model_graph.subtourConstraintsSatisfiedFlag && AUGM_FLAG) {
        // Perform augmentation if conditions are met
        model_graph.perform_augmentation(MAX_AUGM_COUNT);
    }
}