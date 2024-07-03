# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
import update_equality_constraint as equality_constraint
from ut_utils import *

sys.path.append(sys.path[0] + "/../src")

def ut_equality_constraint(
    ut_name,
    n_variables,
    n_constraints,
    filt_flag,
    alpha,
    used_variables,
    variables_connections,
    constraints_connections,
    incoming_metrics_cost,
):
    NUM_VARIABLES = n_variables
    NUM_CONSTRAINTS = n_constraints
    FILT_FLAG = filt_flag
    USED_VARIABLES = used_variables
    ALPHA = alpha
    VARIABLES_CONNECTIONS = variables_connections
    CONSTRAINT_CONNECTIONS = constraints_connections
    INCOMING_METRICS_COST = incoming_metrics_cost

    model_equality_constraint = equality_constraint.EqualityConstraint(NUM_VARIABLES, NUM_CONSTRAINTS, FILT_FLAG, ALPHA, USED_VARIABLES, VARIABLES_CONNECTIONS, CONSTRAINT_CONNECTIONS, INCOMING_METRICS_COST)
    
    f_input1 = os.getcwd() + "/../ut_inputs/used_variables_pure.npy"
    np.save(f_input1, np.array(USED_VARIABLES).astype(np.int32).flatten())
    
    f_input2 = os.getcwd() + "/../ut_inputs/size_used_variables_pure.npy"
    np.save(f_input2, len(USED_VARIABLES))

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)
    
    f_input_n_constraints = os.getcwd() + "/../ut_inputs/num_constraints.npy"
    np.save(f_input_n_constraints, NUM_CONSTRAINTS)

    flattened_variables_connections = [item for sublist in VARIABLES_CONNECTIONS for item in sublist]
    f_input3 = os.getcwd() + "/../ut_inputs/variables_connections_pure.npy"
    np.save(f_input3, np.array(flattened_variables_connections).astype(np.int32).flatten())
    
    sizes_of_variables_connections = [len(sublist) for sublist in VARIABLES_CONNECTIONS]
    f_input4 = os.getcwd() + "/../ut_inputs/sizes_of_variables_connections_pure.npy"
    np.save(f_input4, np.array(sizes_of_variables_connections).astype(np.int32).flatten())
    
    flattened_constraints_connections = [item for sublist in CONSTRAINT_CONNECTIONS for item in sublist]
    f_input5 = os.getcwd() + "/../ut_inputs/constraints_connections_pure.npy"
    np.save(f_input5, np.array(flattened_constraints_connections, dtype=np.int32).flatten())

    f_input6 = os.getcwd() + "/../ut_inputs/incoming_metrics_cost_pure.npy"
    np.save(f_input6, np.array(INCOMING_METRICS_COST, dtype=np_impa_lib).flatten())

    f_input7 = os.getcwd() + "/../ut_inputs/ksat_constraint_to_eq_constraint_m_pure.npy"
    ksat_constraint_to_eq_constraint_m_pure = np.random.uniform(-500, 500, size=(NUM_CONSTRAINTS, NUM_VARIABLES))
    for row_idx, connections in enumerate(CONSTRAINT_CONNECTIONS):
        for col_idx in range(NUM_VARIABLES):
            if col_idx not in connections:
                ksat_constraint_to_eq_constraint_m_pure[row_idx, col_idx] = 0
    ksat_constraint_to_eq_constraint_m_pure = ksat_constraint_to_eq_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input7, ksat_constraint_to_eq_constraint_m_pure.flatten())
    
    #print(f"used_variables: {USED_VARIABLES}")
    #print(f"incoming_metrics_cost: {INCOMING_METRICS_COST}")
    #print(f"VARIABLES_CONNECTIONS: {VARIABLES_CONNECTIONS}")
    
    if ut_name == "VariableEc2KsatConstraintUpdate":
        variable_ec_to_ksat_constraint_m_pure = model_equality_constraint.variable_ec_to_ksat_constraint_update(ksat_constraint_to_eq_constraint_m_pure)
        variable_ec_to_ksat_constraint_m_pure = np.array(variable_ec_to_ksat_constraint_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/variable_ec_to_ksat_constraint_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, variable_ec_to_ksat_constraint_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()