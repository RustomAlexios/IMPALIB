# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")

import graphical_model as model_graph

from ut_utils import *


def ut_model_graph(ut_name, n_nodes, filt_flag, alpha, threshold, random_test_flag, n_iter, sym_flag):
    N_NODES = n_nodes
    FILT_FLAG = filt_flag
    ALPHA = alpha
    THRESHOLD = threshold
    RANDOM_TEST_FLAG = random_test_flag
    N_ITER = n_iter
    SYM_FLAG = sym_flag

    EXACT_SOLVER_FLAG = False
    LKH_SOLVER_FLAG = False
    SIM_ANNEALING_FLAG = False
    NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE = 500
    RESET_FLAG = False

    f_input_alpha = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    f_input_threshold = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/threshold.npy"
    np.save(f_input_threshold, THRESHOLD)

    if ut_name == "IterateRelaxedGraph":
        print("-------")
        print("Python")
        AUGM_FLAG = False

        ModelIMPA = model_graph.GraphicalModel(
            N_ITER,
            N_NODES,
            SYM_FLAG,
            AUGM_FLAG,
            EXACT_SOLVER_FLAG,
            LKH_SOLVER_FLAG,
            SIM_ANNEALING_FLAG,
            NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE,
            RESET_FLAG,
            FILT_FLAG,
            ALPHA,
            THRESHOLD,
            RANDOM_TEST_FLAG,
        )

        ModelIMPA.initialize()

        f_input1 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_connections_pure.npy"
        np.save(f_input1, np.array(ModelIMPA.edge_connections, dtype=np.int32).flatten())

        f_input2 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/cost_edge_variable_pure.npy"
        np.save(f_input2, np.array(ModelIMPA.cost_edge_variable, dtype=np_impa_lib).flatten())

        f_input3 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/cost_matrix_pure.npy"
        np.save(f_input3, np.array(ModelIMPA.cost_matrix, dtype=np_impa_lib).flatten())

        f_input4 = (
            os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_ec_to_degree_constraint_m_pure.npy"
        )
        np.save(f_input4, np.array(ModelIMPA.edge_ec_to_degree_constraint_m, dtype=np_impa_lib).flatten())

        f_input5 = (
            os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_degree_constraint_cost_pure.npy"
        )
        np.save(f_input5, np.array(ModelIMPA.edge_degree_constraint_cost, dtype=np_impa_lib).flatten())

        ModelIMPA.iterate_relaxed_graph()

        ModelIMPA.hard_decision_analysis()

        ModelIMPA.subtour_elimination_constraints_analysis()

        intrinsic_out_edge_ec_pure = ModelIMPA.intrinsic_out_edge_ec
        f_output_path = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_edge_ec_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, intrinsic_out_edge_ec_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()

    elif ut_name == "IterateAugmentedGraph":
        print("-------")
        print("Python")

        AUGM_FLAG = True

        ModelIMPA = model_graph.GraphicalModel(
            N_ITER,
            N_NODES,
            SYM_FLAG,
            AUGM_FLAG,
            EXACT_SOLVER_FLAG,
            LKH_SOLVER_FLAG,
            SIM_ANNEALING_FLAG,
            NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE,
            RESET_FLAG,
            FILT_FLAG,
            ALPHA,
            THRESHOLD,
            RANDOM_TEST_FLAG,
        )

        ModelIMPA.initialize()

        ModelIMPA.save_flag = False

        f_input1 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_connections_pure.npy"
        np.save(f_input1, np.array(ModelIMPA.edge_connections, dtype=np.int32).flatten())

        f_input2 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/cost_edge_variable_pure.npy"
        np.save(f_input2, np.array(ModelIMPA.cost_edge_variable, dtype=np_impa_lib).flatten())

        f_input3 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/cost_matrix_pure.npy"
        np.save(f_input3, np.array(ModelIMPA.cost_matrix, dtype=np_impa_lib).flatten())

        f_input4 = (
            os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_ec_to_degree_constraint_m_pure.npy"
        )
        np.save(f_input4, np.array(ModelIMPA.edge_ec_to_degree_constraint_m, dtype=np_impa_lib).flatten())

        f_input5 = (
            os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/edge_degree_constraint_cost_pure.npy"
        )
        np.save(f_input5, np.array(ModelIMPA.edge_degree_constraint_cost, dtype=np_impa_lib).flatten())

        ModelIMPA.iterate_relaxed_graph()

        ModelIMPA.pre_analysis()

        intrinsic_out_edge_ec_pure = ModelIMPA.intrinsic_out_edge_ec
        f_output_path = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_edge_ec_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, intrinsic_out_edge_ec_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
