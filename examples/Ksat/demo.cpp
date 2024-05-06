// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"

int main()
{
    const int  NUM_CONSTRAINTS          = 2;    ///< number of constraints
    const bool FILTERING_FLAG        = true; ///< whether filtering is activated or not
    const int  NUM_ITERATIONS           = 200;  ///< number of iterations of IMPA
    const int  NUM_VARIABLES = 5;  ///< total number of variables
    const int K_VARIABLE = 3; ///< number of variables per constraint
    const impalib_type ALPHA             = 0.5;     ///< filtering parameter
    const int         NUM_USED_VARIABLES          = 4;    ///< number of used variables to build the formula
    const impalib_type THRESHOLD         = -0.0001; ///< threshold parameter

    array<int, NUM_USED_VARIABLES>             used_variables    = {0, 2, 3, 4};
    vector<vector<int>> variables_connections = {{1}, {}, {0, 1}, {0, 1}, {0}}; ///< each list represents connections to constraints for each variable

    vector<int> flattened_connections;
    for (const auto& variable : variables_connections) {
        flattened_connections.insert(flattened_connections.end(), variable.begin(), variable.end());
    }

    array<int, NUM_VARIABLES> variables_connections_sizes = {1, 0, 2, 2, 1}; ///< each element represents size of connections to constraints for each variable

    array<array<int,K_VARIABLE>, NUM_CONSTRAINTS> constraints_connections = {{{4, 2, 3}, {2,0,3}}}; ///< each list represents connections to variables for each constraint

    array<array<int,K_VARIABLE>, NUM_CONSTRAINTS> constraints_connections_type = {{{-1, -1, -1}, {1,1,-1}}}; ///< each list represents type of connections to variables for each constraint

    array<impalib_type, NUM_VARIABLES>  incoming_metrics_cost    = {-2.66808398, -0.11790238, 4.37984813, -7.13066003, -2.15600797}; ///< incoming metric for each variable

    array<array<impalib_type, NUM_VARIABLES>, NUM_CONSTRAINTS> variable_ec_to_ksat_constraint_m = {
        {{0., 0., 4.37984813, -7.13066003, -2.15600797},
         {-2.66808398, 0., 4.37984813, -7.13066003, 0.}}};

    // Extract data pointers from various structures
    const int          *pUSED_VARIABLES_PY               = used_variables.data();
    const int *pVARIABLES_CONNECTIONS_PY = flattened_connections.data();
    const int *pVARIABLES_CONNECTIONS_SIZES                    = variables_connections_sizes.data();
    const int *pCONSTRAINTS_CONNECTIONS                    = addressof(get<0>(constraints_connections[0]));
    const int *pCONSTRAINTS_CONNECTIONS_TYPE                    = addressof(get<0>(constraints_connections_type[0]));
    const impalib_type *pINCOMING_METRICS_COST = incoming_metrics_cost.data();
    impalib_type *pVariable_ec_to_ksat_constraint_m_py = addressof(get<0>(variable_ec_to_ksat_constraint_m[0]));


    GraphicalModelKsat model_graph(NUM_ITERATIONS, NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILTERING_FLAG, ALPHA, NUM_USED_VARIABLES);

    model_graph.initialize(pUSED_VARIABLES_PY, pVARIABLES_CONNECTIONS_PY, pVARIABLES_CONNECTIONS_SIZES, pCONSTRAINTS_CONNECTIONS, pCONSTRAINTS_CONNECTIONS_TYPE, pINCOMING_METRICS_COST,
                           pVariable_ec_to_ksat_constraint_m_py);

    model_graph.iterate();

    for (int i = 0; i < NUM_VARIABLES; i++)
    {
        cout << model_graph.outputs.ExtrinsicOutputVariableEc[i] << '\n';
    }

}