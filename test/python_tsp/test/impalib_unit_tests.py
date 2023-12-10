# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ut_input_output
import ut_update_degree_constraint
import ut_subtour_elimination_constraint
import ut_update_equality_constraint
import ut_graphical_model
from environmentModule import *


parser = argparse.ArgumentParser()
parser.add_argument('--sub_test_num', type=int,default=1,help="Sub-Test of Unit Test")
parser.add_argument('--sub_tests_total', type=int,default=1,help="Total Number of Sub-Tests of Unit Test")
parser.add_argument('--ut_name', type=str,help="Unit Test Name")
parser.add_argument('--nNodes', type=int, help="Number of Nodes")
parser.add_argument('--nSubtours', type=int, help="Number of Subtours")
parser.add_argument('--symFlag', help='Symmetric FLAG', type=bool)
parser.add_argument('--fFlag', help='FILTERING FLAG', type=bool)
parser.add_argument('--alpha', help='ALPHA', type=np_impa_lib)
parser.add_argument('--augmFlag', help='AUGMENTATION FLAG', type=bool)
parser.add_argument('--threshold', help='THRESHOLD', type=np_impa_lib)
parser.add_argument('--nIter', help='Number of Iterations', type=int)

if __name__ == '__main__':
    args = parser.parse_args() 
    sub_test_num = args.sub_test_num 
    total_sub_tests = args.sub_tests_total 
    ut_name = args.ut_name
    N_NODES = args.nNodes
    N_SUBTOURS = args.nSubtours
    SYM_FLAG = args.symFlag
    FILT_FLAG = args.fFlag
    ALPHA = args.alpha
    AUGM_FLAG = args.augmFlag
    THRESHOLD = args.threshold
    N_ITER = args.nIter
    
    N_EDGE_VARIABLES = N_NODES**2 - N_NODES
    
    upper_triangle = np.random.uniform(1, 1000, size=(N_NODES, N_NODES))
    if (SYM_FLAG):
        matrix = np.triu(upper_triangle) + np.triu(upper_triangle, k=1).T
        np.fill_diagonal(matrix, zero_value)
    else:
        matrix = upper_triangle.copy()
        np.fill_diagonal(matrix, zero_value)   
    
    EDGE_CONNECTIONS = np.transpose(np.nonzero(matrix))

    edge_degree_constraint_cost = np.zeros((N_EDGE_VARIABLES, N_NODES), dtype = np_impa_lib)
    cost_edge_variable = np.zeros(N_EDGE_VARIABLES, dtype = np_impa_lib)
    
    for i in range(len(EDGE_CONNECTIONS)):
        connection = EDGE_CONNECTIONS[i]
        cost = matrix[connection[0], connection[1]]
        cost_edge_variable[i] = cost
        edge_degree_constraint_cost[i][connection[0]] = cost
        edge_degree_constraint_cost[i][connection[1]] = cost
    

    if (ut_name ==  "ExtrinsicOutputEdgeEcRelaxedGraphUpdate" or ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"):
        ut_input_output.ut_io(ut_name = ut_name, n_edge_variables = N_EDGE_VARIABLES, n_nodes = N_NODES, n_subtours = N_SUBTOURS, edge_connections = EDGE_CONNECTIONS)
    elif (ut_name == "DegreeConstraint2EdgeEcUpdate"):
        ut_update_degree_constraint.ut_degree_constraint(ut_name = ut_name, n_edge_variables = N_EDGE_VARIABLES, n_nodes = N_NODES, edge_connections = EDGE_CONNECTIONS, filt_flag = FILT_FLAG, alpha= ALPHA)
    elif (ut_name == "SubtourConstraints2EdgeEcUpdate"):
        ut_subtour_elimination_constraint.ut_subtour_constraint(ut_name = ut_name, n_edge_variables = N_EDGE_VARIABLES, n_nodes = N_NODES, filt_flag = FILT_FLAG, alpha= ALPHA, n_subtours = N_SUBTOURS)
    elif (ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate" or ut_name == "EdgeEc2SubtourConstraintsUpdate" or ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate"):
        ut_update_equality_constraint.ut_equality_constraint(ut_name = ut_name, n_edge_variables = N_EDGE_VARIABLES, n_nodes = N_NODES, filt_flag = FILT_FLAG, alpha= ALPHA, n_subtours = N_SUBTOURS, edge_connections = EDGE_CONNECTIONS, augm_flag = AUGM_FLAG, \
                                                            edge_degree_constraint_cost = edge_degree_constraint_cost, cost_edge_variable = cost_edge_variable)
    elif (ut_name == "IterateRelaxedGraph" or ut_name == "IterateAugmentedGraph"):
        ut_graphical_model.ut_model_graph(ut_name = ut_name, n_nodes = N_NODES, filt_flag = FILT_FLAG, alpha= ALPHA, threshold = THRESHOLD, random_test_flag = True, n_iter = N_ITER, sym_flag = SYM_FLAG)