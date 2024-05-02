# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ut_input_output
import ut_update_ksat_constraint
import ut_update_equality_constraint
import ut_graphical_model
from environmentModule import *

def create_incoming_metrics_cost(num_variables, num_constraints, k_variable, var, type_metrics = 1):
    
    constraints_connections = [[] for _ in range(num_constraints)]
    constraints_connections_type = [[] for _ in range(num_constraints)]
    incoming_metrics_cost = np.zeros(num_variables)
    variables_connections = [[] for _ in range(num_variables)] 
    variables_connections_type = [[] for _ in range(num_variables)]
    used_variables = []
    
    valid_sol = np.random.randint(2, size=num_variables, dtype=int)
    #exit()
    normal_variance = var
    for i, truth_value in enumerate(valid_sol):
        if (type_metrics == 1): # biased
            mean = -normal_variance/2 if truth_value == 1 else normal_variance/2
            incoming_metrics_cost[i] = np.random.normal(loc=mean, scale=np.sqrt(normal_variance))
        elif (type_metrics == 0): # unbiased
            random_sign = np.random.choice([-1, 1])
            mean = random_sign*normal_variance/2
            incoming_metrics_cost[i] = np.random.normal(mean, np.sqrt(normal_variance))

    index_constraint = 0
    for constraint_idx in range(num_constraints):
        connections = []
        connections_types = []
        
        satisfying_variable_index = random.randint(0, num_variables - 1)
        satisfying_variable_value = valid_sol[satisfying_variable_index]
        connection_type = -1 if satisfying_variable_value == 0 else 1
        connections.append(satisfying_variable_index)
        connections_types.append(connection_type)
        used_variables.append(satisfying_variable_index)
        variables_connections[satisfying_variable_index].append(index_constraint)
        variables_connections_type[satisfying_variable_index].append(connection_type)
        
        for k_variable_idx in range(k_variable-1):
            choice_variable = random.choice([idx for idx in range(num_variables) if idx not in connections])
            connections.append(choice_variable)
            choice_connection_type = random.choice([-1, 1])
            connections_types.append(choice_connection_type)
            
            variables_connections[choice_variable].append(index_constraint)
            used_variables.append(choice_variable)
            variables_connections_type[choice_variable].append(choice_connection_type)

        found = False
        for i, sublist_i in enumerate(constraints_connections):
            if(set(sublist_i) == set(connections)):
                index_found = i
                found = True
                #exit()
                
        if (not found):
            constraints_connections[index_constraint] = connections
            constraints_connections_type[index_constraint] = connections_types 
            index_constraint+=1
        else:
            for variable in connections:
                index = variables_connections[variable].index(index_constraint)
                variables_connections[variable].pop(index)
                variables_connections_type[variable].pop(index)
            
    used_variables = list(set(used_variables))
    constraints_connections = [sublist for sublist in constraints_connections if sublist]
    constraints_connections_type = [sublist for sublist in constraints_connections_type if sublist]
    return constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, valid_sol
    
parser = argparse.ArgumentParser()
parser.add_argument("--sub_test_num", type=int, default=1, help="Sub-Test of Unit Test")
parser.add_argument("--sub_tests_total", type=int, default=1, help="Total Number of Sub-Tests of Unit Test")
parser.add_argument("--ut_name", type=str, help="Unit Test Name")
parser.add_argument("--nITER", type=int, default=200, help="Number of Iterations of IMPA")
parser.add_argument("--nVariables", type=int, default=10, help="Number of variables in k-SAT")
parser.add_argument("--nConstraints", type=int, default=3, help="Number of constraints in k-SAT")
parser.add_argument("--kVariable", type=int, default=3, help="Number of variables per constraint")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument("--PPElements", type=int, default=2, help="Number of elements in a PP combination")
parser.add_argument("--filteringFlag", type=bool, default=True, help="Activate Filtering or not")
parser.add_argument("--alpha", type=np_impa_lib, default=0.5, help="Filtering Rate [0,1]")
parser.add_argument("--typeMetrics", type=bool, default=False, help="True: biased, False: unbiased")
parser.add_argument("--var", type=np_impa_lib, default=3, help="Variance of incoming metrics")

if __name__ == "__main__":
    args = parser.parse_args()
    sub_test_num = args.sub_test_num
    total_sub_tests = args.sub_tests_total
    ut_name = args.ut_name
    NUM_ITERATIONS = args.nITER
    NUM_VARIABLES = args.nVariables
    NUM_CONSTRAINTS = args.nConstraints
    K_VARIABLE = args.kVariable
    THRESHOLD = args.threshold
    FILTERING_FLAG = args.filteringFlag
    ALPHA = args.alpha
    POST_PROCESS_FLAG = True
    TYPE_METRICS = args.typeMetrics
    PP_ELEMENTS = args.PPElements
    IM_VARIANCE = args.var
    
    if (ut_name == "VariableEc2KsatConstraintUpdate" or ut_name == "KsatConstraint2VariableEcUpdate" or ut_name=="Iterate"):
        constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, valid_sol = create_incoming_metrics_cost(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, IM_VARIANCE, type_metrics = TYPE_METRICS)
        NUM_CONSTRAINTS = len(constraints_connections)
    
    if ut_name == "ExtrinsicOutputVariableEcUpdate":
        ut_input_output.ut_io(
            ut_name=ut_name,
            n_variables=NUM_VARIABLES,
            n_constraints=NUM_CONSTRAINTS,
        )
    elif ut_name == "VariableEc2KsatConstraintUpdate":
        ut_update_equality_constraint.ut_equality_constraint(
            ut_name=ut_name,
            n_variables=NUM_VARIABLES,
            n_constraints=NUM_CONSTRAINTS,
            filt_flag=FILTERING_FLAG,
            alpha=ALPHA,
            used_variables = used_variables,
            variables_connections = variables_connections,
            constraints_connections = constraints_connections,
            incoming_metrics_cost = incoming_metrics_cost,
        )
    elif ut_name == "KsatConstraint2VariableEcUpdate":
        ut_update_ksat_constraint.ut_ksat_constraint(
            ut_name=ut_name,
            n_variables=NUM_VARIABLES,
            n_constraints=NUM_CONSTRAINTS,
            filt_flag=FILTERING_FLAG,
            alpha=ALPHA,
            constraints_connections = constraints_connections,
            constraint_connections_type = constraints_connections_type,
            k_variable = K_VARIABLE
        )
    elif ut_name == "Iterate":
        ut_graphical_model.ut_model_graph(
            ut_name=ut_name,
            n_variables=NUM_VARIABLES,
            n_constraints=NUM_CONSTRAINTS,
            filt_flag=FILTERING_FLAG,
            alpha=ALPHA,
            threshold=THRESHOLD,
            n_iter=NUM_ITERATIONS,
            k_variable=K_VARIABLE,
        )
