import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import re

def create_cnf_file(constraints_connections, constraints_connections_type, num_variables, num_constraints, cnf_file):
    
    with open(cnf_file, 'w') as f:
        f.write(f"p cnf {num_variables} {num_constraints}\n")
        
        for i, clause_type in enumerate(constraints_connections_type):
            clause = constraints_connections[i]
            for j, var_type in enumerate(clause_type):
                if var_type == -1:
                    f.write('-')
                # add +1 since Mini-SAT solver starts from 1 for indexing variables
                f.write(str(clause[j]+1) + ' ')
            f.write('0\n')

if __name__ == "__main__":
    
    type = "biased"
    n_samples = 500
    var = 10

    folder_inputs = f"../../../data/inputs_ksat_{type}_var_{var}_samples_{n_samples}"
    
    for test_file in range(n_samples):
        
        file_input = f"{folder_inputs}/inputs_set{test_file}.pkl"

        with open(file_input, "rb") as f:
            loaded_data = pkl.load(f)

        num_variables, num_constraints, k_variable, constraints_connections, constraints_connections_type, incoming_metrics_cost, \
                used_variables, variables_connections, variables_connections_type, type_metrics, valid_sol = loaded_data
        
        cnf_file = f"{folder_inputs}/inputs_set{test_file}.cnf"
        create_cnf_file(constraints_connections, constraints_connections_type, num_variables, num_constraints, cnf_file)
        