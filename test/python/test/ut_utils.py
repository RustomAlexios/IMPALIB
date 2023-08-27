# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *
from initializationModule import *

def check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=1e-05, atol = 1e-08):
    ut_failed = True
    y_pure = deepcopy(y_pure.flatten())
    y_wrapper = deepcopy(y_wrapper.flatten())
    assert y_pure.shape == y_wrapper.shape, f'Shape mismatch: Python {y_pure.shape} and C++ {y_wrapper.shape}'
    
    max_absolute_error = np.max(abs(y_pure-y_wrapper))
    max_relative_error = np.max(abs(y_pure-y_wrapper)/abs(np.max(y_pure)+1e-30));
    #index = np.argmax(abs(y_pure-y_wrapper)/abs(np.max(y_pure)+1e-30))
    #print('y_pure: ', y_pure[index])
    #print('y_wrapper: ', y_wrapper[index])
    #print('np.max(y_pure)', abs(np.max(y_pure)+1e-30))
    
    if (max_absolute_error> atol or max_relative_error>rtol):
        ut_failed = True
    else:
        ut_failed = False
    
    if (ut_failed):
        print(f'FAILED SUB-TEST {sub_test_num} out of {total_sub_tests}:: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}')
        #print('y_p: ', y_pure)
        #print('y_w: ', y_wrapper)
    else:
        print(f'PASSED SUB-TEST {sub_test_num} out of {total_sub_tests}:: Test Name: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}')
        #print('y_p: ', y_pure)
        #print('y_w: ', y_wrapper)
        
def check_forward_backward_results(sub_test_num, total_sub_tests, ut_name, rtol=1e-05, atol = 1e-08):
    regex = re.compile(r'\d+')
    if (ut_name == "KnapsackForward"):
        file_names = [file_name for file_name in os.listdir('../ut_results/ut_Knapsack/ut_'+ ut_name) if file_name.startswith("forward_pure")]
    elif (ut_name == "KnapsackBackward"):
        file_names = [file_name for file_name in os.listdir('../ut_results/ut_Knapsack/ut_'+ ut_name) if file_name.startswith("backward_pure")]
    
    file_numbers = []
    for i in range(0, len(file_names)):
        file_numbers.append([int(x) for x in regex.findall(file_names[i])][0])
        
    max_absolute_error = []; ut_failed = []; max_relative_error=[]
    for number in file_numbers:
        flag_failed = True
        if (ut_name == "KnapsackForward"):
            f_path_pure  = '../ut_results/ut_Knapsack/ut_'+ ut_name+'/forward_pure'+str(number)
            f_path_wrapper  = '../ut_results/ut_Knapsack/ut_'+ ut_name+'/forward_wrapper'+str(number)
        elif (ut_name == "KnapsackBackward"):
            f_path_pure  = '../ut_results/ut_Knapsack/ut_'+ ut_name+'/backward_pure'+str(number)
            f_path_wrapper  = '../ut_results/ut_Knapsack/ut_'+ ut_name+'/backward_wrapper'+str(number)
        
        file_array_pure = open(f_path_pure, 'rb'); y_pure = np.load(file_array_pure); y_pure = deepcopy(y_pure.flatten())
        y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib); y_wrapper = deepcopy(y_wrapper.flatten())
        error = abs(y_pure[y_pure!= np_impa_lib(inf)]-y_wrapper[y_wrapper!= np_impa_lib(inf)])
        max_absolute_error.append(np.nanmax(error[error != np_impa_lib(inf)]))
        max_relative_error.append(np.nanmax(error[error != np_impa_lib(inf)])/abs(np.max(y_pure[y_pure!= np_impa_lib(inf)])+1e-30))
        assert y_pure.shape == y_wrapper.shape, f'Shape mismatch: Python {y_pure.shape} and C++ {y_wrapper.shape}'
        if (any(max_absolute_error)> atol or any(max_relative_error)>rtol):
            flag_failed = True
        else:
            flag_failed = False
        ut_failed.append(flag_failed)
        #if (np.allclose(y_pure, y_wrapper, rtol=rtol, atol = atol)):
        #    flag_failed = False
        
        
    if (any(ut_failed)):
        print(f'FAILED SUB-TEST {sub_test_num} out of {total_sub_tests}:: {ut_name}, Max. Abs. Error: {np.max(max_absolute_error):.4e}, Max. Rel. Error: {np.max(max_relative_error):.4e}')
    else:
        print(f'PASSED SUB-TEST {sub_test_num} out of {total_sub_tests}:: Test Name: {ut_name}, Max. Abs. Error: {np.max(max_absolute_error):.4e}, Max. Rel. Error: {np.max(max_relative_error):.4e}')

def prune_teams(N_u):
    
    units = range(1,len(N_u)+1)
    
    permutations_units = [p for p in itertools.product(units, repeat=2)]
    combinations_units_prePruning = list(combinations(units,2))

    list_pruning = []
    '''
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
    '''
    available_combinations = [x for x in permutations_units if x not in list_pruning]

    return available_combinations

def team_reward_generation(available_combinations, N_u, team_types):
    
    teams_weights_per_department = []
    teams_types_per_department = []


    for l in range(0,len(N_u)):
        teams_weights_per_department.append([])
        teams_types_per_department.append([])

    indices_departments = np.array(range(1,len(N_u)+1))

    for combination in available_combinations:
        u = combination[0]; v = combination[1]
        for type in team_types:
            if (type == 1):
                if (u == v):
                    team_size = N_u[u-1]
                    teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                    remaining_departments = np.setdiff1d(indices_departments, v)
                else: 
                    team_size = 0    
            elif (type == 2):
                if (u == v):
                    team_size = math.floor(N_u[u-1]/2)
                    teams_weights_per_department[u-1].append(fighters_count[type-1] + weasels_count[type-1]); \
                        teams_weights_per_department[u-1].extend([fighters_count[type-1] + weasels_count[type-1] for i in range(0, team_size-1)])
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                    remaining_departments = np.setdiff1d(indices_departments, v)
                else:
                    team_size = min(N_u[u-1], N_u[v-1])
                    teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                    teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)])
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                    teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)])
                    remaining_departments = np.setdiff1d(indices_departments, [u, v])                         
            elif (type == 3):
                if (u == v):
                    team_size = math.floor(N_u[u-1]/2)
                    teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                    remaining_departments = np.setdiff1d(indices_departments, v)
                else:
                    team_size = 0
            elif (type == 4):
                if (u == v):
                    team_size = math.floor(N_u[u-1]/3)
                    teams_weights_per_department[u-1].append(fighters_count[type-1]+ weasels_count[type-1]);\
                        teams_weights_per_department[u-1].extend([fighters_count[type-1]+ weasels_count[type-1] for i in range(0, team_size-1)])
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                    remaining_departments = np.setdiff1d(indices_departments, v)   
                else:
                    team_size = math.floor(min(N_u[u-1]/2, N_u[v-1]))
                    teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)]) 
                    teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)]) 
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                    teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)]) 
                    remaining_departments = np.setdiff1d(indices_departments, [u, v])
            elif (type == 5):
                if (u == v):
                    team_size = math.floor(N_u[u-1]/4)
                    teams_weights_per_department[u-1].append(fighters_count[type-1]+ weasels_count[type-1]); \
                        teams_weights_per_department[u-1].extend([fighters_count[type-1]+ weasels_count[type-1] for i in range(0, team_size-1)]) 
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                    remaining_departments = np.setdiff1d(indices_departments, v)
                else:
                    team_size = math.floor(min(N_u[u-1]/2, N_u[v-1]/2))
                    teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)]) 
                    teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)]) 
                    teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                    teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)]) 
                    remaining_departments = np.setdiff1d(indices_departments, [u, v])
                    
            if (team_size != 0):
                for l in remaining_departments:
                    teams_weights_per_department[l-1].append(0); teams_weights_per_department[l-1].extend([0 for i in range(0, team_size-1)]) 
                    teams_types_per_department[l-1].append(type); teams_types_per_department[l-1].extend([type for i in range(0, team_size-1)])

    return teams_weights_per_department, teams_types_per_department