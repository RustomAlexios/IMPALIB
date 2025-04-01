from packageGenerator import *
import numpy as np
import math
import itertools
from itertools import combinations, permutations
import pandas as pd 
from ortools.sat.python import cp_model 
import pickle as pkl 
import pdb
import warnings
import time
warnings.simplefilter(action='ignore', category=FutureWarning)


def parser(N_u, package_types, number_of_targets, package_costs=None, target_costs=None, alpha=1, beta=1, rho=1):

    #Number of Aircrafts in air units
    #N_u = np.array((20, 10, 10, 4, 30, 10, 4, 20, 20, 20, 20), dtype = int)

    #Package Types
    #package_types = [1,2,3,4,5]

    units = range(1,len(N_u)+1)

    #Parameters in the reward function
    F_u = np.array((0, 0, 0, 4, 0, 0, 4, 0, 0, 1, 1), dtype = int)
    W_u = np.array((2, 2, 2, 0, 2, 2, 0, 1, 1, 0, 0), dtype = int)

    distance_metric = rho*np.array(([0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                                [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                                [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                                [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                                [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                                [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                                [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                                [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2],
                                [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2]))



    distance_metric_rendezvouz = np.array((['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E','E', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E','E', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E','E', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E','E', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'W', 'W', 'W', 'W'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'W', 'W', 'W', 'W'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'W', 'W', 'W', 'W'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'W', 'W', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'W', 'W', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E', 'E', 'E', 'E'],
                                        ['E', 'E', 'E', 'E', 'W', 'W', 'W', 'E', 'E', 'E', 'E']))

    permutations_units = [p for p in itertools.product(units, repeat=2)]
    combinations_units_prePruning = list(combinations(units,2))

    #PRUNING STEP
    list_pruning = []
    for i in range(0, len(combinations_units_prePruning)):
        combination = combinations_units_prePruning[i]
        r_1 = math.exp(-F_u[combination[0]-1])*math.exp(-W_u[combination[1]-1])
        r_2 = math.exp(-F_u[combination[1]-1])*math.exp(-W_u[combination[0]-1])
        
        if (r_1 > r_2 and distance_metric[combination[1]-1][combination[0]-1] == distance_metric[combination[0]-1][combination[1]-1] and \
            distance_metric_rendezvouz[combination[1]-1][combination[0]-1] == distance_metric_rendezvouz[combination[0]-1][combination[1]-1]):
            list_pruning.append(combination[::-1])
        elif (r_2 > r_1 and distance_metric[combination[1]-1][combination[0]-1] == distance_metric[combination[0]-1][combination[1]-1]and \
            distance_metric_rendezvouz[combination[1]-1][combination[0]-1] == distance_metric_rendezvouz[combination[0]-1][combination[1]-1]):
            list_pruning.append(combination)
        elif (r_2 == r_1 and \
            distance_metric_rendezvouz[combination[1]-1][combination[0]-1] == distance_metric_rendezvouz[combination[0]-1][combination[1]-1]):
            list_pruning.append(combination[::-1])

    available_combinations = [x for x in permutations_units if x not in list_pruning]

    packages_dict, reward_package, items_weights_per_bin, items_types_per_bin, combination_units, \
        combination_package_size, combination_package, sum_packages, index_map  = \
            package_reward_generation(available_combinations, N_u, package_types, rho)

    sum_packages = np.array(sum_packages)
    last_index = []
    for i in range(0,len(sum_packages)):
        if (i==0):
            last_index.append(sum_packages[i]-1)
        else:
            last_index.append(sum(sum_packages[0:i]) + sum_packages[i]-1)

    items_types = items_types_per_bin[0]

    #number_of_targets  = 2
    rendezvouz_points = ['CV', 'RR']

    reward_target, l_t, k_t, s_t = package_to_target_rewards(last_index, combination_units, items_types, number_of_targets, rendezvouz_points, alpha, beta)

    targets_dict = {}
    i = 0 
    for package_assignment in packages_dict.keys():
        for target_idx, reward in enumerate(reward_target[i]):
            temp = list(package_assignment)
            temp.append(target_idx+1)
            targets_dict[tuple(temp)] = reward 

            
    if not(package_costs is None):
        temp_dict = {}
        if package_costs == 'Random':
            temp_list = [] 
            for key, value in packages_dict.items():
                #temp_dict[key] = (value[0], value[1],  -np.random.uniform(0, 5))
                cost = 8.5*(np.random.uniform(key[2]-1, key[2]) - 5)
                temp_dict[key] = (value[0], value[1],  cost) 
                temp_list.append(cost)
            
            reward_package = temp_list
      
        else: 
            i = 0 
            for key, value in packages_dict.items():
                temp_dict[key] = (value[0], value[1],  package_costs[i])
                i += 1 
        
        packages_dict = temp_dict

    if not(target_costs is None):
        temp_dict = {}
        if target_costs == 'Random':
            temp_array = np.zeros((len(packages_dict), number_of_targets))
            for key, value in targets_dict.items():
                #temp_dict[key] = np.random.uniform(0,5)
                map = [5,4,3,2,1]
                cost = 10*np.random.uniform(map[key[2]-1]-1, map[key[2]-1])
                temp_dict[key] = cost
                temp_array[key[3]-1, key[4]-1] = cost 

            reward_target = temp_array
      
        else: 
            i = 0 
            target_costs = target_costs.flatten()
            for key, value in targets_dict.items():
                temp_dict[key] = target_costs[i]
                i += 1

        targets_dict = temp_dict 

    packages_targets_dict = {}

    for target_assignment, target_cost in targets_dict.items():
        package_assignment = target_assignment[0:4]
        target_package_params = list(packages_dict[package_assignment])
        target_package_params.append(target_cost)
        packages_targets_dict[target_assignment] = tuple(target_package_params)


    return packages_dict, targets_dict, packages_targets_dict, reward_package, reward_target, l_t, k_t, s_t

def extract_results(solver, vars):
    results_list = []
    for assignment in vars.keys():
        if solver.BooleanValue(vars[assignment]):
            results_list.append(assignment)
    return results_list

def LHS(packages_dict, N_u):
    packages_df = pd.DataFrame(packages_dict.keys())
    model = cp_model.CpModel() 

    packages_assignment_vars = {} 
    for assignment in packages_dict.keys():
        packages_assignment_vars[tuple(assignment)] = model.NewBoolVar(str(assignment[0]) + " " + str(assignment[1]) + " "
                                                                     + str(assignment[2]) + " " + str(assignment[3]))

    # Add capacity constraint for each unit. 
    for unit in range(len(N_u)):
        u_list = list(packages_df[packages_df[0] == unit + 1].values)
        v_list = list(packages_df[packages_df[1] == unit + 1].values)

        u_prod = [packages_dict[tuple(package)][0] * packages_assignment_vars[tuple(package)] for package in u_list]
        v_prod = [packages_dict[tuple(package)][1] * packages_assignment_vars[tuple(package)] for package in v_list]
        model.Add(sum(u_prod + v_prod) <= N_u[unit])

    # Add our objective. 
    #model.Minimize(sum(int(1e6*packages_dict[tuple(package)][2]) * packages_assignment_vars[tuple(package)] for package in packages_dict.keys()))
    model.Minimize(sum(int(100*packages_dict[tuple(package)][2]) * packages_assignment_vars[tuple(package)] for package in packages_dict.keys()))
    # Solve our model.
    solver = cp_model.CpSolver()
    #print('LHS num_search_workers: ', solver.parameters.num_search_workers)
    #solver.parameters.linearization_level = 0  
    solver.Solve(model)
    print(f'Solve status: {solver.StatusName()}')
    return extract_results(solver, packages_assignment_vars)

def RHS(targets_dict):
    targets_df = pd.DataFrame(targets_dict.keys())
    model = cp_model.CpModel() 

    targets_assignment_vars = {} 
    for assignment in targets_dict.keys():
        targets_assignment_vars[tuple(assignment)] = model.NewBoolVar(str(assignment[0]) + " " + str(assignment[1]) + " "
                                                                        + str(assignment[2]) + " " + str(assignment[3]) + " "
                                                                        + str(assignment[4]))


    for target in set(targets_df[4]):
        packages = list(targets_df[targets_df[4] == target].values)
        model.Add(sum(targets_assignment_vars[tuple(assignment)] for assignment in packages) <= 1)

    for package in np.unique(targets_df.iloc[:,[0,1,2,3]].values, axis=0):
        targets = list(targets_df[(targets_df[0]==package[0]) &
                            (targets_df[1]==package[1]) &
                            (targets_df[2]==package[2]) &
                            (targets_df[3]==package[3])].values)
        model.Add(sum(targets_assignment_vars[tuple(assignment)]for assignment in targets) <= 1)

    model.Minimize(sum(int(-100*targets_dict[tuple(assignment)]) * targets_assignment_vars[tuple(assignment)] for assignment in targets_dict.keys())) # NEGATIVE COEF FOR TESTING ONLY
    solver = cp_model.CpSolver()
    #print('RHS num_search_workers: ', solver.parameters.num_search_workers)
    #solver.parameters.linearization_level = 0  
    solver.Solve(model)
    print(f'Solve status: {solver.StatusName()}')
    return extract_results(solver, targets_assignment_vars)

def composed(pt_dict, N_u):
    pt_df = pd.DataFrame(pt_dict.keys())
    model = cp_model.CpModel() 
    #print('Composed num_search_workers: ', solver.parameters.num_search_workers)
    # Define variables.
    pt_assignment_vars = {} 
    for assignment in pt_dict.keys():
        pt_assignment_vars[tuple(assignment)] = model.NewBoolVar(str(assignment[0]) + " " + str(assignment[1]) + " "
                                                                        + str(assignment[2]) + " " + str(assignment[3]) + " "
                                                                        + str(assignment[4]))

    # Constraint 3. 
    for unit in range(len(N_u)):
        u_list = list(pt_df[pt_df[0] == unit + 1].values)
        v_list = list(pt_df[pt_df[1] == unit + 1].values)

        u_prod = [pt_dict[tuple(package)][0] * pt_assignment_vars[tuple(package)] for package in u_list]
        v_prod = [pt_dict[tuple(package)][1] * pt_assignment_vars[tuple(package)] for package in v_list]

        model.Add(sum(u_prod + v_prod) <= N_u[unit])

    # Constraint 8a.
    for package in np.unique(pt_df.iloc[:,[0,1,2,3]].values, axis=0):
        targets = list(pt_df[(pt_df[0]==package[0]) &
                            (pt_df[1]==package[1]) &
                            (pt_df[2]==package[2]) &
                            (pt_df[3]==package[3])].values)
        model.Add(sum(pt_assignment_vars[tuple(assignment)]for assignment in targets) <= 1)

    # Constraint 8b.
    for target in set(pt_df[4]):
        packages = list(pt_df[pt_df[4] == target].values)
        model.Add(sum(pt_assignment_vars[tuple(assignment)] for assignment in packages) <= 1)

    lhs_obj = sum(int(100*pt_dict[tuple(assignment)][2]) * pt_assignment_vars[tuple(assignment)] for assignment in pt_dict.keys())
    rhs_obj = sum(int(100*pt_dict[tuple(assignment)][3]) * pt_assignment_vars[tuple(assignment)] for assignment in pt_dict.keys()) 

    model.Minimize(lhs_obj + rhs_obj) 
    solver = cp_model.CpSolver()
    #print(solver.num_search_parameters)
    solver.parameters.num_search_workers = 1
    # solver.parameters.max_time_in_seconds = 100
    #print('Composed num_search_workers: ', solver.parameters.num_search_workers)
    solver.parameters.linearization_level = 2
    start_time = time.time()
    solver.Solve(model)
    runtime = time.time() - start_time
    print(f'Solve status: {solver.StatusName()}')
    return runtime, extract_results(solver, pt_assignment_vars)

def human(pt_dict, N_u, l, k, s):
    df = pd.DataFrame(pt_dict.keys())
    targets = list(set(df[4]))
    package_type_map = np.array([[1,0], [1,1], [2,0], [2,1], [2,2]])
    counts = np.copy(N_u)
    results = [] 
    for target in targets:
        k_t = k[target-1]
        s_t = s[target-1]
        l_t = l[target-1]

        package = [k_t, s_t]
        
        if package == [1,2]:
            package = [2,2]
            k_t = 2 
            s_t = 2 

        package_type = np.where(np.all(package_type_map == package, axis=1))[0][0] + 1 

        # Select west package. 
        # Select fighters
        
        if counts[4] >= k_t:
            #counts[4] -= k_t 
            w_k_unit = 5 

        elif counts[5] >= k_t:
            #counts[5] -= k_t
            w_k_unit = 6 

        elif counts[7] >= k_t:
            #counts[7] -= k_t
            w_k_unit = 8

        elif counts[8] >= k_t:
            #counts[8] -= k_t
            w_k_unit = 9

        else:
            w_k_unit = False

        # Select weasles.
        if s_t != 0: 
            if counts[6] >= s_t:
                #counts[6] -= s_t
                w_s_unit = 7

            elif counts[9] >= s_t:
                #counts[9] -= s_t
                w_s_unit = 10

            elif counts[10] >= s_t:
                #counts[10] -= s_t
                w_s_unit = 11

            else: 
                w_s_unit = False

        else:
            # If we arent selecting weasles, we set the w_s_unit to the w_k_unit per the
            # convention in pt_dict. 
            w_s_unit = w_k_unit

      # Select east package. 
        # Select fighters
        if counts[0] >= k_t:
            #counts[0] -= k_t 
            e_k_unit = 1 

        elif counts[1] >= k_t:
            #counts[1] -= k_t
            e_k_unit = 2

        elif counts[2] >= k_t:
            #counts[2] -= k_t
            e_k_unit = 3

        elif counts[7] >= k_t:
            #counts[7] -= k_t
            e_k_unit = 8

        elif counts[8] >= k_t:
            #counts[8] -= k_t
            e_k_unit = 9

        else:
            e_k_unit = False

        # Select weasles.
        if s_t != 0: 
            if counts[3] >= s_t:
                #counts[3] -= s_t
                e_s_unit = 4

            elif counts[9] >= s_t:
                #counts[9] -= s_t
                e_s_unit = 10

            elif counts[10] >= s_t:
                #counts[10] -= s_t
                e_s_unit = 11

            else: 
                e_s_unit = False

        else:
            # If we arent selecting weasles, we set the w_s_unit to the w_k_unit per the
            # convention in pt_dict. 
            e_s_unit = e_k_unit


        if np.all([w_k_unit, w_s_unit, e_k_unit, e_s_unit]):
            # Since both west and east packages are feasible, break the tie using distance.
            if l_t < .5: 
                assignment = (w_k_unit, w_s_unit, package_type, target)
            else:
                assignment = (e_k_unit, e_s_unit, package_type, target) 
            
        elif np.all([w_k_unit, w_s_unit]):
            assignment = (w_k_unit, w_s_unit, package_type, target)

        elif np.all([e_k_unit, e_s_unit]): 
            assignment = (e_k_unit, e_s_unit, package_type, target) 

        else: 
            # No assignment is supported, so skip target. 
            continue
        
        # This is just to deal with the indexing of pt_dict.
        assignment = tuple(df[(df[0]==assignment[0])&(df[1]==assignment[1])
                        &(df[2]==assignment[2])&(df[4]==assignment[3])].iloc[0])


        counts[assignment[0]-1] -= package[0]
        counts[assignment[1]-1] -= package[1]
        if np.any(counts<0):
            pdb.set_trace()
        results.append(assignment)

    return results


if __name__ == '__main__':
    N_u = np.array((20, 10, 10, 4, 30, 10, 4, 20, 20, 20, 20), dtype = int)
    print(N_u.shape)
    package_types = [1,2,3,4,5]
    packages_dict, targets_dict, pt_dict, reward_package, reward_target = parser(N_u, package_types, number_of_targets = 4, package_costs='Random', target_costs='Random')
    results_lhs = LHS(packages_dict, N_u)
    print(pd.DataFrame(results_lhs))
    results_rhs = RHS(targets_dict)
    print(pd.DataFrame(results_rhs))
    results_composed = composed(pt_dict, N_u)
    print(pd.DataFrame(results_composed))
    inputs  = [N_u, package_types, reward_package, reward_target]

    with open('inputs.pkl', 'wb') as f:
        pkl.dump(inputs, f)

    with open('outputs.pkl', 'wb') as f:
        pkl.dump(results_composed, f)

    with open('inputs.pkl', 'rb') as f:
        input_load = pkl.load(f)

    with open('outputs.pkl', 'rb') as f:
        output_load = pkl.load(f)

    pdb.set_trace()
    