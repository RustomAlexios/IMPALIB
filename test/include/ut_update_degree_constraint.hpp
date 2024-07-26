// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
// https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_degree_constraint(string&);

void ut_degree_constraint(string& ut_name) {
    const char* n_nodes_bash = getenv("N_NODES");
    if (n_nodes_bash == NULL) {
        cout << "n_nodes_bash not available\n";
    }

    const char* filt_flag_bash = getenv("FILT_FLAG");
    if (filt_flag_bash == NULL) {
        cout << "filt_flag_bash not available\n";
    }

    const int N_NODES = atoi(n_nodes_bash);
    const int N_EDGE_VARIABLES = N_NODES * N_NODES - N_NODES;
    const bool FILT_FLAG(filt_flag_bash);

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/edge_ec_to_degree_constraint_m_pure.npy");
    impalib_type* edge_ec_to_degree_constraint_m_pure = input1.data<impalib_type>();

    vector<vector<impalib_type>> edge_ec_to_degree_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

    for (int edge_variable_index = 0; edge_variable_index < N_EDGE_VARIABLES; edge_variable_index++) {
        copy(edge_ec_to_degree_constraint_m_pure + N_NODES * edge_variable_index, edge_ec_to_degree_constraint_m_pure + N_NODES * (edge_variable_index + 1),
             edge_ec_to_degree_constraint_m[edge_variable_index].begin());
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/edge_connections_pure.npy");
    int* edge_connections_pure = input2.data<int>();

    int num_connections = 2;

    vector<vector<int>> edge_connections(N_EDGE_VARIABLES, vector<int>(num_connections, 0));

    for (int edge_variable_index = 0; edge_variable_index < N_EDGE_VARIABLES; edge_variable_index++) {
        copy(edge_connections_pure + num_connections * edge_variable_index, edge_connections_pure + num_connections * (edge_variable_index + 1), edge_connections[edge_variable_index].begin());
    }

    DegreeConstraint modelDegreeConstraint(N_NODES, N_EDGE_VARIABLES, FILT_FLAG, ALPHA);

    if (ut_name == "DegreeConstraint2EdgeEcUpdate") {
        vector<vector<impalib_type>> degree_constraint_to_eq_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

        modelDegreeConstraint.degree_constraint_to_edge_ec_update(edge_ec_to_degree_constraint_m, edge_connections, degree_constraint_to_eq_constraint_m);

        fstream file_output("../ut_results/degree_constraint_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios::trunc);
        if (file_output.is_open()) {
            for (int i = 0; i < N_EDGE_VARIABLES; i++) {
                for (int j = 0; j < N_NODES; j++) {
                    file_output.write((char*)(&degree_constraint_to_eq_constraint_m[i][j]), sizeof(degree_constraint_to_eq_constraint_m[i][j]));
                }
            }
            file_output.close();
        } else {
            cout << "Error! File cannot be opened!"
                 << "\n";
        }
    }
}
