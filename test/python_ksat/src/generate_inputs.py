# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)
          
import numpy as np
import random
import os
import pickle as pkl

def create_incoming_metrics_cost(num_variables, num_constraints, k_variable, type_metrics = 0, normal_variance=3):
    
    constraints_connections = [[] for _ in range(num_constraints)]
    constraints_connections_type = [[] for _ in range(num_constraints)]
    incoming_metrics_cost = np.zeros(num_variables)
    variables_connections = [[] for _ in range(num_variables)] 
    variables_connections_type = [[] for _ in range(num_variables)]
    used_variables = []
    
    valid_sol = np.random.randint(2, size=num_variables, dtype=int)
    
    #normal_variance = 3
    
    for i, truth_value in enumerate(valid_sol):
        if (type_metrics == 1): # correctly biased
            mean = -normal_variance/2 if truth_value == 1 else normal_variance/2
            incoming_metrics_cost[i] = np.random.normal(loc=mean, scale=np.sqrt(normal_variance))
        elif (type_metrics == 0): # random 
            random_sign = np.random.choice([-1, 1])
            mean = random_sign*normal_variance/2
            incoming_metrics_cost[i] = np.random.normal(mean, np.sqrt(normal_variance))
        elif (type_metrics == 2): #normal(0,sigma^2)
            incoming_metrics_cost[i] = np.random.normal(0, np.sqrt(normal_variance))
    
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

if __name__ == "__main__":
    
    np_impa_lib = np.float32
    zero_value = np_impa_lib(0)
    samples = 500
    #type_metrics = 1 #biased
    #type_metrics = 0 #unbiased
    type_metrics = 2 #N(0,sigma^2)
    
    N_VARIABLES_MAX = 200 #200 for unbiased_500 and biased_500
    N_CONSTRAINTS_MAX = 500 #200 for unbiased_500 and biased_500
    K_VARIABLE_MIN = 2 #2 for unbiased_500 and biased_500
    
    if (type_metrics==1):
        type = "biased"
    elif (type_metrics==0):
        type = "unbiased"
    elif (type_metrics==2):
        type = "zeromean"
    
    for index_sample in range(0, samples):
        
        normal_variance = 10
        num_variables = np.random.randint(5, N_VARIABLES_MAX+1)
        num_constraints = np.random.randint(2, N_CONSTRAINTS_MAX+1)
        k_variable = np.random.randint(K_VARIABLE_MIN, num_variables)
        constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, valid_sol = create_incoming_metrics_cost(num_variables, num_constraints, k_variable, type_metrics, normal_variance)
        num_constraints = len(constraints_connections)
        input = [
            num_variables,
            num_constraints,
            k_variable,
            constraints_connections,
            constraints_connections_type, 
            incoming_metrics_cost,
            used_variables,
            variables_connections, 
            variables_connections_type,
            type_metrics,
            valid_sol
        ]
        similar_elements = []
        
        for i, sublist_i in enumerate(constraints_connections):
            for j, sublist_j in enumerate(constraints_connections):
                if i<j and set(sublist_i) == set(sublist_j):
                    similar_elements.append((i, j))
                    print("ERROR")
                    exit()
                    
        if not (os.path.isdir(f"../../../data/inputs_ksat_{type}_var_{normal_variance}_samples_{samples}")):
            os.makedirs(f"../../../data/inputs_ksat_{type}_var_{normal_variance}_samples_{samples}")

        with open(
            f"../../../data/inputs_ksat_{type}_var_{normal_variance}_samples_{samples}/inputs_set{str(index_sample)}.pkl",
            "wb",
        ) as f:
            print(f"Test file: {index_sample} & num_variables: {num_variables} & num_constraints: {num_constraints} & k_variable: {k_variable}")
            pkl.dump(input, f)