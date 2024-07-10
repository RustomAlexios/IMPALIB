# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
import update_degree_constraint as degree_constraint
from ut_utils import *


def ut_degree_constraint(ut_name, n_edge_variables, n_nodes, edge_connections, filt_flag, alpha):
    N_EDGE_VARIABLES = n_edge_variables
    N_NODES = n_nodes
    EDGE_CONNECTIONS = edge_connections
    FILT_FLAG = filt_flag
    ALPHA = alpha

    f_input1 = os.getcwd() + "/../ut_inputs/edge_ec_to_degree_constraint_m_pure.npy"
    edge_ec_to_degree_constraint_m_pure = np.random.uniform(10, 500, size=(N_EDGE_VARIABLES, N_NODES))
    edge_ec_to_degree_constraint_m_pure = edge_ec_to_degree_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input1, edge_ec_to_degree_constraint_m_pure.flatten())

    f_input2 = os.getcwd() + "/../ut_inputs/edge_connections_pure.npy"
    np.save(f_input2, np.array(EDGE_CONNECTIONS, dtype=np.int32).flatten())

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    model_degree_constraint = degree_constraint.DegreeConstraint(N_NODES, N_EDGE_VARIABLES, EDGE_CONNECTIONS, FILT_FLAG, ALPHA)

    if ut_name == "DegreeConstraint2EdgeEcUpdate":
        model_degree_constraint.degree_constraint_to_edge_ec_update(edge_ec_to_degree_constraint_m_pure)
        degree_constraint_to_eq_constraint_m_pure = model_degree_constraint.degree_constraint_to_eq_constraint_m_dummy
        f_output_path = os.getcwd() + "/../ut_results/degree_constraint_to_eq_constraint_m_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, degree_constraint_to_eq_constraint_m_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
