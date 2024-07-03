# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
import graphical_model as model_graph
from ut_utils import *

sys.path.append(sys.path[0] + "/../src")

def ut_model_graph(ut_name, n_variables, n_constraints, filt_flag, alpha, threshold, n_iter, k_variable):
    
    print("-------")
    print("Python Pure")
    
    NUM_VARIABLES = n_variables
    NUM_CONSTRAINTS = n_constraints
    FILTERING_FLAG = filt_flag
    ALPHA = alpha
    NUM_ITERATIONS = n_iter
    THRESHOLD = threshold
    K_VARIABLE = k_variable
    
    RANDOM_TEST_FLAG = True

    EXACT_SOLVER_FLAG = False
    
    POST_PROCESS_FLAG = True

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    f_input_threshold = os.getcwd() + "/../ut_inputs/threshold.npy"
    np.save(f_input_threshold, THRESHOLD)

    if ut_name == "Iterate":

        ModelIMPA = model_graph.GraphicalModel(
            NUM_ITERATIONS,
            NUM_VARIABLES,
            NUM_CONSTRAINTS,
            K_VARIABLE,
            EXACT_SOLVER_FLAG,
            FILTERING_FLAG,
            ALPHA,
            THRESHOLD,
            RANDOM_TEST_FLAG,
            POST_PROCESS_FLAG
        )

        ModelIMPA.initialize()
        
        NUM_USED_VARIABLES = len(ModelIMPA.used_variables)
    
        USED_VARIABLES = ModelIMPA.used_variables
        f_input1 = os.getcwd() + "/../ut_inputs/used_variables_pure.npy"
        np.save(f_input1, np.array(USED_VARIABLES).astype(np.int32).flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/size_used_variables_pure.npy"
        np.save(f_input2, len(USED_VARIABLES))
        
        f_input_n_constraints = os.getcwd() + "/../ut_inputs/num_constraints.npy"
        np.save(f_input_n_constraints, ModelIMPA.num_constraints)
        
        VARIABLES_CONNECTIONS = ModelIMPA.variables_connections
        flattened_variables_connections = [item for sublist in VARIABLES_CONNECTIONS for item in sublist]
        f_input3 = os.getcwd() + "/../ut_inputs/variables_connections_pure.npy"
        np.save(f_input3, np.array(flattened_variables_connections).astype(np.int32).flatten())
        
        sizes_of_variables_connections = [len(sublist) for sublist in VARIABLES_CONNECTIONS]
        f_input4 = os.getcwd() + "/../ut_inputs/sizes_of_variables_connections_pure.npy"
        np.save(f_input4, np.array(sizes_of_variables_connections).astype(np.int32).flatten())
    
        CONSTRAINT_CONNECTIONS = ModelIMPA.constraints_connections
        flattened_constraints_connections = [item for sublist in CONSTRAINT_CONNECTIONS for item in sublist]
        f_input5 = os.getcwd() + "/../ut_inputs/constraints_connections_pure.npy"
        np.save(f_input5, np.array(flattened_constraints_connections, dtype=np.int32).flatten())

        CONSTRAINTS_CONNECTIONS_TYPE = ModelIMPA.constraints_connections_type
        flattened_constraints_connections_type = [item for sublist in CONSTRAINTS_CONNECTIONS_TYPE for item in sublist]
        f_input6 = os.getcwd() + "/../ut_inputs/constraints_connections_type_pure.npy"
        np.save(f_input6, np.array(flattened_constraints_connections_type, dtype=np.int32).flatten())
        
        INCOMING_METRICS_COST = ModelIMPA.incoming_metrics_cost
        f_input7 = os.getcwd() + "/../ut_inputs/incoming_metrics_cost_pure.npy"
        np.save(f_input7, np.array(INCOMING_METRICS_COST, dtype=np_impa_lib).flatten())

        f_input8 = os.getcwd() + "/../ut_inputs/variable_ec_to_ksat_constraint_m_pure.npy"
        np.save(f_input8, np.array(ModelIMPA.variable_ec_to_ksat_constraint_m, dtype=np_impa_lib).flatten())
        
        ModelIMPA.run_impa()
        
        extrinsic_out_variable_ec_pure = ModelIMPA.extrinsic_output_variable_ec
        f_output_path = os.getcwd() + "/../ut_results/extrinsic_out_variable_ec_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, extrinsic_out_variable_ec_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
