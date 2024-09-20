// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_subtour_constraint(string&);

void ut_subtour_constraint(string& ut_name) {
    const char* n_nodes_bash = getenv("N_NODES");
    if (n_nodes_bash == NULL) {
        cout << "n_nodes_bash not available\n";
    }

    const char* n_subtours_bash = getenv("N_SUBTOURS");
    if (n_subtours_bash == NULL) {
        cout << "n_subtours_bash not available\n";
    }

    const char* filt_flag_bash = getenv("FILT_FLAG");
    if (filt_flag_bash == NULL) {
        cout << "filt_flag_bash not available\n";
    }

    const int N_NODES = atoi(n_nodes_bash);
    const int N_SUBTOURS = atoi(n_subtours_bash);
    const int N_EDGE_VARIABLES = N_NODES * N_NODES - N_NODES;
    const bool FILT_FLAG(filt_flag_bash);

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    SubtourEliminationConstraint modelSubtourConstraint(N_NODES, N_EDGE_VARIABLES, FILT_FLAG, ALPHA);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/edge_ec_to_subtour_constraints_m_pure.npy");
    impalib_type* edge_ec_to_subtour_constraints_m_pure = input1.data<impalib_type>();

    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES, zero_value));

    for (int subtour_index = 0; subtour_index < N_SUBTOURS; subtour_index++) {
        copy(edge_ec_to_subtour_constraints_m_pure + N_EDGE_VARIABLES * subtour_index, edge_ec_to_subtour_constraints_m_pure + N_EDGE_VARIABLES * (subtour_index + 1),
             edge_ec_to_subtour_constraints_m[subtour_index].begin());
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/delta_S_indices_list_sizes_pure.npy");
    const int* delta_S_indices_list_sizes_pure = input2.data<int>();

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/delta_S_indices_list_pure.npy");
    const int* delta_S_indices_list_pure = input3.data<int>();

    vector<vector<int>> delta_S_indices_list;

    int idx = 0;
    for (int subtour_index = 0; subtour_index < N_SUBTOURS; subtour_index++) {
        int size = delta_S_indices_list_sizes_pure[subtour_index];
        vector<int> temp_vector(delta_S_indices_list_pure + idx, delta_S_indices_list_pure + idx + size);
        delta_S_indices_list.push_back(temp_vector);
        idx += size;
    }

    if (ut_name == "SubtourConstraints2EdgeEcUpdate") {
        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES, zero_value));

        modelSubtourConstraint.subtour_constraints_to_edge_ec_update(edge_ec_to_subtour_constraints_m, delta_S_indices_list, subtour_constraints_to_edge_ec_m);

        fstream file_output("../ut_results/subtour_constraints_to_edge_ec_m_wrapper", ios::out | ios::binary | ios::trunc);
        if (file_output.is_open()) {
            for (int i = 0; i < N_SUBTOURS; i++) {
                for (int j = 0; j < N_EDGE_VARIABLES; j++) {
                    file_output.write((char*)(&subtour_constraints_to_edge_ec_m[i][j]), sizeof(subtour_constraints_to_edge_ec_m[i][j]));
                }
            }
            file_output.close();
        } else {
            cout << "Error! File cannot be opened!"
                 << "\n";
        }
    }
}
