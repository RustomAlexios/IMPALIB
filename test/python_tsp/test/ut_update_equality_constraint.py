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
    n_edge_variables,
    n_nodes,
    filt_flag,
    alpha,
    n_subtours,
    edge_connections,
    augm_flag,
    edge_degree_constraint_cost,
    cost_edge_variable,
):
    N_EDGE_VARIABLES = n_edge_variables
    N_NODES = n_nodes
    FILT_FLAG = filt_flag
    EDGE_CONNECTIONS = edge_connections
    ALPHA = alpha
    N_SUBTOURS = n_subtours
    AUGM_FLAG = augm_flag

    model_equality_constraint = equality_constraint.EqualityConstraint(N_NODES, N_EDGE_VARIABLES, 0, EDGE_CONNECTIONS, edge_degree_constraint_cost, AUGM_FLAG, FILT_FLAG, ALPHA)

    delta_S_indices_list_length = N_SUBTOURS
    delta_S_indices_list = [random.sample(range(N_EDGE_VARIABLES), k=random.randint(2, N_NODES)) for _ in range(delta_S_indices_list_length)]

    delta_S_indices_list_sizes = [len(lst) for lst in delta_S_indices_list]
    delta_S_indices_list_flattened = sum(delta_S_indices_list, [])

    f_input1 = os.getcwd() + "/../ut_inputs/degree_constraint_to_eq_constraint_m_pure.npy"
    degree_constraint_to_eq_constraint_m_pure = np.random.uniform(10, 500, size=(N_EDGE_VARIABLES, N_NODES))
    for row_idx, edge_indices in enumerate(EDGE_CONNECTIONS):
        for col_idx in range(N_NODES):
            if col_idx not in edge_indices:
                degree_constraint_to_eq_constraint_m_pure[row_idx, col_idx] = 0
    degree_constraint_to_eq_constraint_m_pure = degree_constraint_to_eq_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input1, degree_constraint_to_eq_constraint_m_pure.flatten())

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    f_input2 = os.getcwd() + "/../ut_inputs/edge_degree_constraint_cost_pure.npy"
    edge_degree_constraint_cost_pure = edge_degree_constraint_cost.astype(np_impa_lib)
    np.save(f_input2, edge_degree_constraint_cost_pure.flatten())

    f_input3 = os.getcwd() + "/../ut_inputs/edge_connections_pure.npy"
    np.save(f_input3, np.array(EDGE_CONNECTIONS, dtype=np.int32).flatten())

    if ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate":
        edge_ec_to_degree_constraint_m_pure = model_equality_constraint.edge_ec_to_degree_constraint_relaxed_graph_update(degree_constraint_to_eq_constraint_m_pure)
        edge_ec_to_degree_constraint_m_pure = np.array(edge_ec_to_degree_constraint_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/edge_ec_to_degree_constraint_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, edge_ec_to_degree_constraint_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
    elif ut_name == "EdgeEc2SubtourConstraintsUpdate":
        f_input4 = os.getcwd() + "/../ut_inputs/delta_S_indices_list_pure.npy"
        np.save(f_input4, np.array(delta_S_indices_list_flattened, dtype=np.int32).flatten())

        f_input5 = os.getcwd() + "/../ut_inputs/delta_S_indices_list_sizes_pure.npy"
        np.save(f_input5, np.array(delta_S_indices_list_sizes, dtype=np.int32).flatten())

        f_input6 = os.getcwd() + "/../ut_inputs/subtour_constraints_to_edge_ec_m_pure.npy"
        subtour_constraints_to_edge_ec_m_pure = np.random.uniform(10, 500, size=(N_SUBTOURS, N_EDGE_VARIABLES))
        for row_idx, subtour_indices in enumerate(delta_S_indices_list):
            for col_idx in range(N_EDGE_VARIABLES):
                if col_idx not in subtour_indices:
                    subtour_constraints_to_edge_ec_m_pure[row_idx, col_idx] = 0
        subtour_constraints_to_edge_ec_m_pure = subtour_constraints_to_edge_ec_m_pure.astype(np_impa_lib)
        np.save(f_input6, subtour_constraints_to_edge_ec_m_pure.flatten())

        f_input7 = os.getcwd() + "/../ut_inputs/cost_edge_variable_pure.npy"
        cost_edge_variable_pure = cost_edge_variable.astype(np_impa_lib)
        np.save(f_input7, cost_edge_variable_pure.flatten())

        edge_ec_to_subtour_constraints_m_pure = model_equality_constraint.edge_ec_to_subtour_constraints_update(
            degree_constraint_to_eq_constraint_m_pure,
            subtour_constraints_to_edge_ec_m_pure,
            delta_S_indices_list,
            cost_edge_variable_pure,
        )
        edge_ec_to_subtour_constraints_m_pure = np.array(edge_ec_to_subtour_constraints_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/edge_ec_to_subtour_constraints_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, edge_ec_to_subtour_constraints_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
    elif ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate":
        f_input8 = os.getcwd() + "/../ut_inputs/subtour_constraints_to_edge_ec_m_pure.npy"
        subtour_constraints_to_edge_ec_m_pure = np.random.uniform(10, 500, size=(N_SUBTOURS, N_EDGE_VARIABLES))
        for row_idx, subtour_indices in enumerate(delta_S_indices_list):
            for col_idx in range(N_EDGE_VARIABLES):
                if col_idx not in subtour_indices:
                    subtour_constraints_to_edge_ec_m_pure[row_idx, col_idx] = 0
        subtour_constraints_to_edge_ec_m_pure = subtour_constraints_to_edge_ec_m_pure.astype(np_impa_lib)
        np.save(f_input8, subtour_constraints_to_edge_ec_m_pure.flatten())

        edge_ec_to_degree_constraint_m_pure = model_equality_constraint.edge_ec_to_degree_constraint_augmented_graph_update(
            degree_constraint_to_eq_constraint_m_pure, subtour_constraints_to_edge_ec_m_pure
        )
        edge_ec_to_degree_constraint_m_pure = np.array(edge_ec_to_degree_constraint_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/edge_ec_to_degree_constraint_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, edge_ec_to_degree_constraint_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
