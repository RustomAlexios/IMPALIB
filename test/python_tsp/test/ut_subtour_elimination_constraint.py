# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
import subtour_elimination_constraint as subtour_constraint
from ut_utils import *

sys.path.append(sys.path[0] + "/../src")

def ut_subtour_constraint(ut_name, n_edge_variables, n_nodes, filt_flag, alpha, n_subtours):
    N_EDGE_VARIABLES = n_edge_variables
    N_NODES = n_nodes
    FILT_FLAG = filt_flag
    ALPHA = alpha
    N_SUBTOURS = n_subtours

    model_subtour_constraint = subtour_constraint.SubtourEliminationConstraint(N_NODES, N_EDGE_VARIABLES, FILT_FLAG, ALPHA)

    delta_S_indices_list_length = N_SUBTOURS
    delta_S_indices_list = [random.sample(range(N_EDGE_VARIABLES), k=random.randint(2, N_NODES)) for _ in range(delta_S_indices_list_length)]

    f_input1 = os.getcwd() + "/../ut_inputs/edge_ec_to_subtour_constraints_m_pure.npy"
    edge_ec_to_subtour_constraints_m_pure = np.random.uniform(10, 500, size=(delta_S_indices_list_length, N_EDGE_VARIABLES))
    edge_ec_to_subtour_constraints_m_pure = edge_ec_to_subtour_constraints_m_pure.astype(np_impa_lib)
    np.save(f_input1, edge_ec_to_subtour_constraints_m_pure.flatten())

    delta_S_indices_list_sizes = [len(lst) for lst in delta_S_indices_list]
    delta_S_indices_list_flattened = sum(delta_S_indices_list, [])

    f_input2 = os.getcwd() + "/../ut_inputs/delta_S_indices_list_pure.npy"
    np.save(f_input2, np.array(delta_S_indices_list_flattened, dtype=np.int32).flatten())

    f_input3 = os.getcwd() + "/../ut_inputs/delta_S_indices_list_sizes_pure.npy"
    np.save(f_input3, np.array(delta_S_indices_list_sizes, dtype=np.int32).flatten())

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    if ut_name == "SubtourConstraints2EdgeEcUpdate":
        subtour_constraints_to_edge_ec_m_pure = model_subtour_constraint.subtour_constraints_to_edge_ec_update(edge_ec_to_subtour_constraints_m_pure, delta_S_indices_list)
        subtour_constraints_to_edge_ec_m_pure = np.array(subtour_constraints_to_edge_ec_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/subtour_constraints_to_edge_ec_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, subtour_constraints_to_edge_ec_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
