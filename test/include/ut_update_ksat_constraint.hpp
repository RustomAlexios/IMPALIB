// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_ksat_constraint(string&);

void ut_ksat_constraint(string& ut_name) {
    const char* n_variables_bash = getenv("NUM_VARIABLES");
    if (n_variables_bash == NULL) {
        cout << "n_variables_bash not available\n";
    }

    const char* k_variable_bash = getenv("K_VARIABLE");
    if (k_variable_bash == NULL) {
        cout << "k_variable_bash not available\n";
    }

    const char* filt_flag_bash = getenv("FILT_FLAG");
    if (filt_flag_bash == NULL) {
        cout << "filt_flag_bash not available\n";
    }

    const int NUM_VARIABLES = atoi(n_variables_bash);
    const int K_VARIABLE = atoi(k_variable_bash);
    const bool FILT_FLAG(filt_flag_bash);

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input_num_constraints = cnpy::npy_load("../ut_inputs/num_constraints.npy");
    int* num_constraints_pure = input_num_constraints.data<int>();
    const int NUM_CONSTRAINTS = *num_constraints_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/constraints_connections_pure.npy");
    int* constraints_connections_pure = input1.data<int>();
    vector<vector<int>> constraints_connections(NUM_CONSTRAINTS, vector<int>(K_VARIABLE, 0));

    for (int constraint_index = 0; constraint_index < NUM_CONSTRAINTS; constraint_index++) {
        copy(constraints_connections_pure + K_VARIABLE * constraint_index, constraints_connections_pure + K_VARIABLE * (constraint_index + 1), constraints_connections[constraint_index].begin());
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/constraints_connections_type_pure.npy");
    int* constraints_connections_type_pure = input2.data<int>();
    vector<vector<int>> constraints_connections_type(NUM_CONSTRAINTS, vector<int>(K_VARIABLE, 0));

    for (int constraint_index = 0; constraint_index < NUM_CONSTRAINTS; constraint_index++) {
        copy(constraints_connections_type_pure + K_VARIABLE * constraint_index, constraints_connections_type_pure + K_VARIABLE * (constraint_index + 1),
             constraints_connections_type[constraint_index].begin());
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/variable_ec_to_ksat_constraint_m_pure.npy");
    impalib_type* variable_ec_to_ksat_constraint_m_pure = input3.data<impalib_type>();
    vector<vector<impalib_type>> variable_ec_to_ksat_constraint_m(NUM_CONSTRAINTS, vector<impalib_type>(NUM_VARIABLES, zero_value));

    for (int constraint_index = 0; constraint_index < NUM_CONSTRAINTS; constraint_index++) {
        copy(variable_ec_to_ksat_constraint_m_pure + NUM_VARIABLES * constraint_index, variable_ec_to_ksat_constraint_m_pure + NUM_VARIABLES * (constraint_index + 1),
             variable_ec_to_ksat_constraint_m[constraint_index].begin());
    }

    KsatConstraint modelKsatConstraint(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILT_FLAG, ALPHA);

    if (ut_name == "KsatConstraint2VariableEcUpdate") {
        vector<vector<impalib_type>> ksat_constraint_to_eq_constraint_m(NUM_CONSTRAINTS, vector<impalib_type>(NUM_VARIABLES, zero_value));

        modelKsatConstraint.ksat_constraint_to_variable_ec_update(variable_ec_to_ksat_constraint_m, ksat_constraint_to_eq_constraint_m, constraints_connections, constraints_connections_type);

        fstream file_output("../ut_results/ksat_constraint_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios::trunc);
        if (file_output.is_open()) {
            for (int i = 0; i < NUM_CONSTRAINTS; i++) {
                for (int j = 0; j < NUM_VARIABLES; j++) {
                    file_output.write((char*)(&ksat_constraint_to_eq_constraint_m[i][j]), sizeof(ksat_constraint_to_eq_constraint_m[i][j]));
                }
            }
            file_output.close();
        } else {
            cout << "Error! File cannot be opened!"
                 << "\n";
        }
    }
}
