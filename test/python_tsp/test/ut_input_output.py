# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + '/../src')

import input_output
from ut_utils import *

    
def ut_io(ut_name, n_edge_variables, n_nodes, n_subtours, edge_connections):
    
    N_EDGE_VARIABLES = n_edge_variables
    N_NODES = n_nodes
    N_SUBTOURS = n_subtours
    EDGE_CONNECTIONS = edge_connections 

    f_input1  = os.getcwd() + '/../ut_inputs/ut_InputOutput/ut_'+ ut_name+'/degree_constraint_to_eq_constraint_m_pure.npy'
    degree_constraint_to_eq_constraint_m_pure = np.random.uniform(10,500, size = (N_EDGE_VARIABLES, N_NODES))
    for row_idx, edge_indices in enumerate(EDGE_CONNECTIONS):
        for col_idx in range(N_NODES):
            if col_idx not in edge_indices:
                degree_constraint_to_eq_constraint_m_pure[row_idx, col_idx] = 0 
    degree_constraint_to_eq_constraint_m_pure = degree_constraint_to_eq_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input1, degree_constraint_to_eq_constraint_m_pure.flatten())
    
    outputs = input_output.OutputsImpa(N_NODES, N_EDGE_VARIABLES, 0, 0)
    
    if (ut_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate"):        
        extrinsic_output_edge_ec_relaxed_graph_pure = outputs.extrinsic_output_edge_ec_relaxed_graph_update(degree_constraint_to_eq_constraint_m_pure)
        f_output_path  = os.getcwd() + '/../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_relaxed_graph_pure'
        output_file_python_pure = open(f_output_path, 'wb')
        np.save(output_file_python_pure, extrinsic_output_edge_ec_relaxed_graph_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()

        
    elif (ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"):
        delta_S_indices_list_length = N_SUBTOURS
        delta_S_indices_list = [random.sample(range(N_EDGE_VARIABLES), k=random.randint(2, N_NODES)) for _ in range(delta_S_indices_list_length)]
        subtour_constraints_to_edge_ec_m_pure = np.random.uniform(10,500, size = (N_SUBTOURS, N_EDGE_VARIABLES))
        for row_idx, subtour_indices in enumerate(delta_S_indices_list):
            for col_idx in range(N_EDGE_VARIABLES):
                if col_idx not in subtour_indices:
                    subtour_constraints_to_edge_ec_m_pure[row_idx, col_idx] = 0
        subtour_constraints_to_edge_ec_m_pure = subtour_constraints_to_edge_ec_m_pure.astype(np_impa_lib)
        f_input2 = os.getcwd() + '/../ut_inputs/ut_InputOutput/ut_'+ ut_name+'/subtour_constraints_to_edge_ec_m_pure.npy'
        np.save(f_input2, subtour_constraints_to_edge_ec_m_pure.flatten())
        
        extrinsic_output_edge_ec_augmented_graph_pure = outputs.extrinsic_output_edge_ec_augmented_graph_update(degree_constraint_to_eq_constraint_m_pure, subtour_constraints_to_edge_ec_m_pure)
        f_output_path  = os.getcwd() + '/../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_augmented_graph_pure'
        output_file_python_pure = open(f_output_path, 'wb')
        np.save(output_file_python_pure, extrinsic_output_edge_ec_augmented_graph_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
        
        