# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")

from ut_utils import *
import update_ksat_constraint as ksat_constraint
            
def ut_ksat_constraint(ut_name, n_variables, n_constraints, filt_flag, alpha, constraints_connections, constraint_connections_type, k_variable):
    
    NUM_VARIABLES = n_variables
    NUM_CONSTRAINTS = n_constraints
    FILT_FLAG = filt_flag
    ALPHA = alpha
    CONSTRAINTS_CONNECTIONS = constraints_connections
    CONSTRAINTS_CONNECTIONS_TYPE = constraint_connections_type
    
    model_ksat_constraint = ksat_constraint.KSatConstraint(NUM_VARIABLES, NUM_CONSTRAINTS, FILT_FLAG, ALPHA, CONSTRAINTS_CONNECTIONS, CONSTRAINTS_CONNECTIONS_TYPE)

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)
    
    f_input_n_constraints = os.getcwd() + "/../ut_inputs/num_constraints.npy"
    np.save(f_input_n_constraints, NUM_CONSTRAINTS)
    
    flattened_constraints_connections = [item for sublist in CONSTRAINTS_CONNECTIONS for item in sublist]
    f_input1 = os.getcwd() + "/../ut_inputs/constraints_connections_pure.npy"
    np.save(f_input1, np.array(flattened_constraints_connections, dtype=np.int32).flatten())

    flattened_constraints_connections_type = [item for sublist in CONSTRAINTS_CONNECTIONS_TYPE for item in sublist]
    f_input1 = os.getcwd() + "/../ut_inputs/constraints_connections_type_pure.npy"
    np.save(f_input1, np.array(flattened_constraints_connections_type, dtype=np.int32).flatten())

    f_input2 = os.getcwd() + "/../ut_inputs/variable_ec_to_ksat_constraint_m_pure.npy"
    variable_ec_to_ksat_constraint_m_pure = np.random.uniform(-500, 500, size=(NUM_CONSTRAINTS, NUM_VARIABLES))
    for row_idx, connections in enumerate(CONSTRAINTS_CONNECTIONS):
        for col_idx in range(k_variable):
            if col_idx not in connections:
                variable_ec_to_ksat_constraint_m_pure[row_idx, col_idx] = 0
    variable_ec_to_ksat_constraint_m_pure = variable_ec_to_ksat_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input2, variable_ec_to_ksat_constraint_m_pure.flatten())

    if ut_name == "KsatConstraint2VariableEcUpdate":
        model_ksat_constraint.ksat_constraint_to_variable_ec_update(variable_ec_to_ksat_constraint_m_pure)
        ksat_constraint_to_eq_constraint_m_pure = model_ksat_constraint.ksat_constraint_to_eq_constraint_m_dummy
        f_output_path = os.getcwd() + "/../ut_results/ksat_constraint_to_eq_constraint_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, ksat_constraint_to_eq_constraint_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
