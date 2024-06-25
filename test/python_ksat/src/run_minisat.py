# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)
          
import numpy as np
import os
import pickle as pkl
from pysat.solvers import Minisat22
import time

def read_cnf_file(file_path, solver):
    constraints_connections = []
    constraints_connections_type = []
    used_variables = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()

        for index, line in enumerate(lines):
            if line.startswith('p cnf'):
                index_start = index +1
                _, _, num_variables, num_constraints = line.strip().split()
                num_variables = int(num_variables)
                num_constraints = int(num_constraints)
                break
        
        variables_connections = [[] for _ in range(num_variables)] 
        variables_connections_type = [[] for _ in range(num_variables)]
        
        for index, line in enumerate(lines):
            if line.startswith('c') or line.startswith('p') or line.startswith('%') or line.startswith('\n') or line.startswith('0'):
                continue
            else:
                constraint = line.strip().split()[:-1]
                constraint_number_list = [int(x) for x in constraint]
                solver.add_clause(constraint_number_list)
                connections = []
                connections_type = []
                constraint_index = index - index_start
                for variable in constraint:
                    if variable.startswith('-'):
                        variable_index = int(variable[1:])-1
                        variable_type = -1
                        used_variables.append(variable_index)
                        variables_connections[variable_index].append(constraint_index)
                        connections.append(variable_index)
                        connections_type.append(variable_type)
                        variables_connections_type[variable_index].append(variable_type)
                    else:
                        variable_index = int(variable)-1
                        variable_type = 1
                        used_variables.append(variable_index)
                        variables_connections[variable_index].append(constraint_index)
                        connections.append(variable_index)
                        connections_type.append(variable_type)
                        variables_connections_type[variable_index].append(variable_type)
                constraints_connections.append(connections)
                constraints_connections_type.append(connections_type)
                if constraint_index == 0:
                    k_variable = len(connections)
                elif k_variable != len(connections):
                    raise ValueError("All constraints must have the same number of variables.")
                
                
    used_variables = list(set(used_variables))
    
    return num_variables, num_constraints, constraints_connections, constraints_connections_type, used_variables, variables_connections, k_variable, variables_connections_type, solver


if __name__ == "__main__":
    np_impa_lib = np.float32
    zero_value = np_impa_lib(0)
    normal_variance = 10
    #type_metrics = 1 #biased
    type_metrics = 0 #unbiased
    #type_metrics = 2 #N(0,sigma^2)
    save_flag = False
    n_samples = 48
    
    if (type_metrics==1):
        type = "biased"
    elif (type_metrics==0):
        type = "unbiased"
    elif (type_metrics==2):
        type = "zeromean"
    
    file_name = "aim" #uf.50.218.1000 #uf125.538.100 #aim
    folder_path = f"../../../data/inputs_{type}_var{normal_variance}_{file_name}"
    folder_outputs = f"../../../data/outputs_{type}_var{normal_variance}_{file_name}_minisat"
    #folder_path = f"../../../data/inputs_ksat_{type}_var_{normal_variance}_samples_{n_samples}"
    #folder_outputs = f"../../../data/outputs_ksat_{type}_var_{normal_variance}_samples_{n_samples}_minisat"
    all_files = os.listdir(folder_path)

    files = [os.path.join(folder_path, file) for file in all_files if os.path.isfile(os.path.join(folder_path, file)) if file.endswith('.cnf')]

    if not (os.path.exists(f"{folder_outputs}")) and save_flag:
        os.makedirs(f"{folder_outputs}")
            
    for index_sample, file_path in enumerate(files):

        # do not use number detection when dealing with aim dataset 
        #file_name = os.path.basename(file_path)
        #number = file_name.replace("inputs_set", "").replace(".cnf", "")
        #print(number)
        results_composed = []
        
        solver = Minisat22()

        num_variables, num_constraints, constraints_connections, constraints_connections_type, used_variables, \
            variables_connections, k_variable, variables_connections_type, solver = read_cnf_file(file_path, solver)
        
        start_time = time.time()
        solver.solve()
        end_time = time.time()
        
        variables = solver.get_model()
        solver.delete()
        
        formula_satisfied = False
        
        hard_decision = [1 if x > 0 else 0 for x in variables]
        
        active_variables = [i for i, value in enumerate(hard_decision) if value]
        inactive_variables = [i for i, value in enumerate(hard_decision) if not value]
        
        indices_satisfied_constraints_solid = [(i, j) for i, sublist in enumerate(variables_connections_type) for j, element in enumerate(sublist) if element == 1 and i in active_variables]
        satisfied_constraints_solid = list(set([variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_solid]))
        
        indices_satisfied_constraints_dashed = [(i, j) for i, sublist in enumerate(variables_connections_type) for j, element in enumerate(sublist) if element == -1 and i in inactive_variables]
        satisfied_constraints_dashed = list(set([variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_dashed])) 
        
        satisfied_constraints = list(set(satisfied_constraints_solid + satisfied_constraints_dashed))
        
        unsatisfied_constraints = [constraint for constraint in range(num_constraints) if constraint not in satisfied_constraints]
        
        if (len(unsatisfied_constraints)):
            formula_satisfied = False
            print(file_path)
        else:
            formula_satisfied = True
        
        results_composed.append(
                (
                    type_metrics,
                    num_variables,
                    num_constraints,
                    k_variable,
                    constraints_connections,
                    constraints_connections_type,
                    used_variables,
                    variables_connections,
                    variables_connections_type,
                    active_variables,
                    inactive_variables,
                    satisfied_constraints,
                    unsatisfied_constraints,
                    formula_satisfied,
                    end_time - start_time
                )
            )
        
        if (save_flag):
            with open(f"{folder_outputs}/outputs_set{index_sample}.pkl", "wb") as f: # for aim dataset
            #with open(f"{folder_outputs}/outputs_set{number}.pkl", "wb") as f:
                pkl.dump(results_composed, f)