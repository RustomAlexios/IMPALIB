# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impalib import *

class GraphicalModelKcMwm:
    
    def __init__(self,NUM_ITERATIONS, FILTERING_FLAG, POST_PROCESS_FLAG, ALPHA, THRESHOLD):
        
        self.num_iterations = NUM_ITERATIONS
        self.filtering_flag = FILTERING_FLAG 
        self.post_process_flag = POST_PROCESS_FLAG
        self.alpha = ALPHA
        self.threshold = THRESHOLD
    
    def initialize(self, input_load, test_flag = False):
        
        max_state = np.array(input_load[0])
        team_types = input_load[1]#.tolist()
        num_projects = len(input_load[3][0])
        
        self.team_types = team_types
        self.num_projects = num_projects
        
        #print('N_u: ', max_state)
        #print('team_types: ', team_types)
        #print('number_of_targets: ', num_projects)
        
        #print('Number of Iterations: ', self.num_iterations)
        
        #if (self.filtering_flag):
        #    print('Alpha: ', self.alpha)

        #if (self.post_process_flag):
        #    print('post_processing_Flag Enabled')
            
        self.max_state = max_state
        num_departments = max_state.size
        self.num_departments = num_departments
        
        l_target = np.random.uniform(0, 1, num_projects)
        k_target = np.random.randint(min(fighters_count), max(fighters_count)+1, num_projects)
        s_target = np.random.randint(min(weasels_count),max(weasels_count)+1 , num_projects)
        
        self.prune_teams(test_flag=test_flag)
        
        reward_team, teams_weights_per_department, teams_types_per_department, last_index = \
                        self.team_reward_generation()

        self.teams_weights_per_department = teams_weights_per_department
        self.last_index = last_index
        #print('Total Number of Targets: ', num_projects)
        
        self.teams_types = teams_types_per_department[0]

        reward_team = np.array(input_load[2], dtype = np_impa_lib)
        self.reward_team = reward_team
        reward_project = np.array(input_load[3], dtype=np_impa_lib).T
        self.reward_project = reward_project
        
        num_teams = len(reward_team)
        self.num_teams = num_teams
        
        #print('-------')

        self.intrinsic_out_mwm = np.zeros((self.num_projects, self.num_teams), dtype=np_impa_lib)
        
        if (self.num_projects != self.num_teams):
            self.unbalanced_flag = True
        else:
            self.unbalanced_flag = False
            
        assert self.intrinsic_out_mwm.shape[0] == self.reward_project.shape[0], f'Row Shape mismatch between MWM and Rewards'
        assert self.intrinsic_out_mwm.shape[1] == self.reward_project.shape[1], f'Column Shape mismatch between MWM and Rewards'
        
        all_projects = np.array(range(0,num_projects))
        self.all_projects = all_projects
        
        self.team_to_oric_m = np.zeros(num_teams, dtype = np_impa_lib)
        self.oric_to_team_m = np.zeros(num_teams, dtype = np_impa_lib)
        self.eq_constraint_to_oric_m = np.zeros((num_projects, num_teams), dtype = np_impa_lib)
        self.project_to_eq_constraint_m = np.zeros((num_projects, num_teams), dtype = np_impa_lib)
        self.eq_constraint_to_project_m = np.zeros((num_projects, num_teams), dtype = np_impa_lib)
        self.oric_to_eq_constraint_m  = np.zeros((num_projects, num_teams), dtype = np_impa_lib)
        self.intrinsic_out_mwm = np.zeros((num_projects, num_teams), dtype = np_impa_lib)
        
        team_to_knapsack_m = np.zeros((num_departments,num_teams), dtype = np_impa_lib)
        
        for i in range(0,num_departments):
            for j in range(0,num_teams):
                if (teams_weights_per_department[i][j] != 0):
                    team_to_knapsack_m[i][j] = reward_team[j]

        self.team_to_knapsack_m = team_to_knapsack_m
        
        self.model_knapsacks = Knapsack(self.num_departments, self.num_teams, self.filtering_flag, self.alpha, self.reward_team)
        self.project_ineq_constraint = InequalityConstraint(self.num_departments, self.num_teams, self.num_projects, self.unbalanced_flag)
        self.model_oric = OrInequalityConstraint(self.num_departments, self.num_teams, self.num_projects, self.unbalanced_flag, self.reward_project)
        self.model_eq_constraint = EqualityConstraint(self.num_departments, self.num_teams, self.num_projects, self.reward_team, self.reward_project)
        self.outputs = OutputsImpa(self.num_departments, self.num_teams, self.num_projects, self.reward_project)
    
    def prune_teams(self, test_flag = False):
    
        units = range(1,len(self.max_state)+1)
        permutations_units = [p for p in itertools.product(units, repeat=2)]
        combinations_units_pre_pruning = list(combinations(units,2))

        list_pruning = []
        if (not test_flag):
            for i in range(0, len(combinations_units_pre_pruning)):
                combination = combinations_units_pre_pruning[i]
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
            
        self.available_combinations = [x for x in permutations_units if x not in list_pruning]
    
    def team_reward_generation(self, test_flag = False):
        
        teams_weights_per_department = []
        teams_types_per_department = []
        reward_team = []
        total_team_array = []

        for l in range(0,len(self.max_state)):
            teams_weights_per_department.append([])
            teams_types_per_department.append([])

        sum_teams = []
        
        max_state = self.max_state

        indices_departments = np.array(range(1,len(max_state)+1))

        for combination in self.available_combinations:
            u = combination[0]; v = combination[1]
            team_array = []
            for type in self.team_types:
                if (type == 1):
                    if (u == v):
                        team_size = max_state[u-1]
                        teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1)
                    else: 
                        team_size = 0    
                elif (type == 2):
                    if (u == v):
                        team_size = math.floor(max_state[u-1]/2)
                        teams_weights_per_department[u-1].append(fighters_count[type-1] + weasels_count[type-1]); \
                            teams_weights_per_department[u-1].extend([fighters_count[type-1] + weasels_count[type-1] for i in range(0, team_size-1)])
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])
                    else:
                        team_size = min(max_state[u-1], max_state[v-1])
                        teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                        teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)])
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                        teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])                         
                elif (type == 3):
                    if (u == v):
                        team_size = math.floor(max_state[u-1]/2)
                        teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)])
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1)
                    else:
                        team_size = 0
                elif (type == 4):
                    if (u == v):
                        team_size = math.floor(max_state[u-1]/3)
                        teams_weights_per_department[u-1].append(fighters_count[type-1]+ weasels_count[type-1]);\
                            teams_weights_per_department[u-1].extend([fighters_count[type-1]+ weasels_count[type-1] for i in range(0, team_size-1)])
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])  
                    else:
                        team_size = math.floor(min(max_state[u-1]/2, max_state[v-1]))
                        teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)]) 
                        teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)]) 
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                        teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)]) 
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])
                elif (type == 5):
                    if (u == v):
                        team_size = math.floor(max_state[u-1]/4)
                        teams_weights_per_department[u-1].append(fighters_count[type-1]+ weasels_count[type-1]); \
                            teams_weights_per_department[u-1].extend([fighters_count[type-1]+ weasels_count[type-1] for i in range(0, team_size-1)]) 
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])
                    else:
                        team_size = math.floor(min(max_state[u-1]/2, max_state[v-1]/2))
                        teams_weights_per_department[u-1].append(fighters_count[type-1]); teams_weights_per_department[u-1].extend([fighters_count[type-1] for i in range(0, team_size-1)]) 
                        teams_weights_per_department[v-1].append(weasels_count[type-1]); teams_weights_per_department[v-1].extend([weasels_count[type-1] for i in range(0, team_size-1)]) 
                        teams_types_per_department[u-1].append(type); teams_types_per_department[u-1].extend([type for i in range(0, team_size-1)]) 
                        teams_types_per_department[v-1].append(type); teams_types_per_department[v-1].extend([type for i in range(0, team_size-1)]) 
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        if (not test_flag):
                            reward = np_impa_lib(distance_metric[v-1][u-1] +  F_u[u-1] - K_1 + W_u[v-1])
                        
                if (team_size != 0):
                    if (not test_flag):
                        reward_team.append(np_impa_lib(reward))
                        reward_team.extend([np_impa_lib(reward) for i in range(0, team_size-1)])
                    for l in remaining_departments:
                        teams_weights_per_department[l-1].append(0); teams_weights_per_department[l-1].extend([0 for i in range(0, team_size-1)]) 
                        teams_types_per_department[l-1].append(type); teams_types_per_department[l-1].extend([type for i in range(0, team_size-1)]) 
                team_array.append(team_size)
            
            sum_teams.append(sum(team_array))
            total_team_array.append(team_array)
            #print(' u: ', u-1, '  v: ', v-1, '  Packages Size:  ', team_array, 'with sum: ', sum(team_array), '& last index: ', np.sum(sum_teams)-1)  

        #print('Total Number of Packages: ', np.sum(sum_teams))
        
        sum_teams = np.array(sum_teams, dtype = np_impa_lib)
        last_index = []
        for i in range(0,len(sum_teams)):
            if (i==0):
                last_index.append(sum_teams[i]-1)
            else:
                last_index.append(sum(sum_teams[0:i]) + sum_teams[i]-1)

        return reward_team, teams_weights_per_department, teams_types_per_department, last_index
    
    def iterate(self):
        
        for iter in range(0, self.num_iterations):
            
            for bin_index in range(0, self.num_departments):

                self.model_knapsacks.forward_backward(self.max_state[bin_index], self.team_to_knapsack_m[bin_index], np.array(self.teams_weights_per_department[bin_index]))
                
                self.model_knapsacks.extrinsic_output_department_lhs(bin_index, self.max_state[bin_index], self.team_to_knapsack_m[bin_index], np.array(self.teams_weights_per_department[bin_index]))
                
            self.model_knapsacks.process_extrinsic_output_department(iter)
            
            self.team_to_oric_m = self.model_eq_constraint.team_ec_to_oric_update(self.model_knapsacks.extrinsic_output_department)
            
            self.oric_to_eq_constraint_m, self.eq_constraint_to_project_m  = self.model_oric.oric_to_project_ec_update(self.eq_constraint_to_oric_m, self.team_to_oric_m)
            
            self.project_to_eq_constraint_m  = self.project_ineq_constraint.project_inequality_constraint_update(self.eq_constraint_to_project_m)
            
            self.eq_constraint_to_oric_m = self.model_eq_constraint.project_eq_const_to_oric_update(self.project_to_eq_constraint_m)

            self.oric_to_team_m = self.model_oric.oric_to_team_update(self.eq_constraint_to_oric_m)
            
            self.intrinsic_out_mwm = self.outputs.intrinsic_out_mwm_update(self.oric_to_eq_constraint_m, self.project_to_eq_constraint_m)
            
            self.team_to_knapsack_m = self.model_knapsacks.team_to_department_update(self.teams_weights_per_department, self.oric_to_team_m)
            
            self.extrinsic_output_team = self.outputs.extrinsic_output_team_update(self.model_knapsacks.extrinsic_output_department, self.oric_to_team_m)

        self.intrinsic_output = self.extrinsic_output_team + self.reward_team
    
    def pre_analysis(self):
    
        intrinsic_output = self.intrinsic_output
        threshold = self.threshold
        
        valid_match, impa_metric, p_mwm = self.check_match()
        self.matching_array = np.where(p_mwm==1)
        selected_teams_indices = np.where(p_mwm==1)[1]
        self.selected_teams_indices = np.unique(selected_teams_indices)
        
        hard_decision = np.array(deepcopy(intrinsic_output))
        hard_decision[intrinsic_output>threshold] = 0
        
        self.pre_processing_analysis()
        
    def check_match(self):
        
        A = self.intrinsic_out_mwm
        
        unbalanced_flag = self.unbalanced_flag
        
        match_metric  = sum(x for x in A.flatten() if x<0)
        P = deepcopy(A)
        threshold = self.threshold
        P[P>threshold] = 0
        P[P<threshold] = 1
        valid_match = False

        min_N = min(A.shape[0],A.shape[1])

        if (unbalanced_flag):
            P_interm = self.check_unbalanced_flag(P)
        else:
            P_interm = P
        row_sums = np.sum(P_interm, axis = 1)
        col_sums = np.sum(P_interm, axis = 0)

        if np.array_equal(col_sums, np.ones(min_N)):
            if np.array_equal(row_sums, np.ones(min_N)):
                valid_match = True   
        return valid_match, match_metric, P
    
    def check_unbalanced_flag(self, P):
        idx = 0; idy = 0
        P_interm = deepcopy(P)
        idy = np.argwhere(np.all(P_interm[...,:] == 0, axis=0)) #find which columns there are zeros
        P_interm = np.delete(P_interm, idy, axis=1)
        idx = np.argwhere(np.all(P_interm[:,...] == 0, axis=1)) #find which rows there are zeross
        P_interm = np.delete(P_interm, idx, axis=0)
        return P_interm
    
    def pre_processing_analysis(self):
        
        num_departments = self.num_departments
        num_projects = self.num_projects
        
        selected_teams_indices = self.selected_teams_indices
        matching_array = self.matching_array
        teams_types = self.teams_types
        teams_weights_per_department = self.teams_weights_per_department
        reward_team = self.reward_team
        last_index = self.last_index
        reward_project = self.reward_project
        available_combinations = self.available_combinations
        max_state = self.max_state
        
        capacity_exceeded_flag = False
        project_matched_flag = False
        ic_violated_flag = False
        
        all_projects  = np.array(range(0,num_projects))
        
        team_departments_weights = np.zeros((num_departments), dtype = int)

        department_weights =[]

        for j in range(0, num_departments):
            department_weights.append([])
        
        print('Matching Array: ', matching_array)
        
        print('-------------')
        print('Air Units Assignment')
        print('-------------')
        
        team_type = []
        team_locations = []
        total_value = 0
        total_weight = 0
        for i in range(0,len(selected_teams_indices)):
            team_departments = []
            index = selected_teams_indices[i]
            project_index = np.where(matching_array[1]==index)[0]
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            for j in range(0, num_departments):
                weight = teams_weights_per_department[j][index]
                if (weight !=0):
                    team_departments.append(j)
                    department_weights[j].append(weight)
                    team_departments_weights[j] = team_departments_weights[j] + weight
            location_team  = next(x[0] for x in enumerate(last_index) if x[1] >= index)
            team_locations.append(location_team)
            print('Team ', index, ' has type ', teams_types[index], ' & is assigned to project ',matching_array[0][project_index] ,' & is in units ', np.array(available_combinations[location_team])-1)

        print('-------')
        print('BEFORE POST-PROCESSING')
        for i in range(0, num_departments):
            print('Unit ', i, ' has a used capacity ', team_departments_weights[i])
            print('Unit ', i, ' has a MAX capacity ', max_state[i] )
            if (team_departments_weights[i] > max_state[i]):
                print('Unit ', i, ' needs post-processing')
                capacity_exceeded_flag = True
        print('Exceeded Capacity Flag: ', capacity_exceeded_flag)   
        print('-------')
        print('Total Value LHS: ', -total_value)
        print('Total number of selected teams: ', len(selected_teams_indices), 'out of ', len(reward_team))
        print('Total number of Matched Projects: ', len(matching_array[0]), 'out of ', num_projects)
        print('Total number of Unmatched Projects: ', len(np.setdiff1d(all_projects, matching_array[0])))

        print('-------')
        print('Types of Teams: ', np.unique(team_type, return_counts=True)[0], 'Counts: ', np.unique(team_type, return_counts=True)[1])
        print('-------------')

        total_weight = 0
        visited_teams = []
        reused_teams = []
        combination_reused_team = []
        visited_combinations = []
        for i in range(0, len(matching_array[0])):
            team_index = matching_array[1][i]
            project_index = matching_array[0][i]
            if (team_index not in visited_teams):
                visited_teams.append(team_index)
                visited_combinations.append((project_index, team_index))
            else:
                reused_teams.append(team_index)
                index_used = visited_teams.index(team_index)
                used_combination = visited_combinations[index_used]
                combination_reused_team.append([(project_index, team_index), used_combination])
                ic_violated_flag = True
            total_weight = total_weight + reward_project[project_index, team_index]

        print('MWM')
        print('-------------')

        if (len(matching_array[0])==num_projects):
            project_matched_flag = True

        print(f'Matching All Projects to Teams : {project_matched_flag}')
        print('Team assigned to multiple projects: ', ic_violated_flag, 'with index: ', np.unique(reused_teams))
        print('Total Weight RHS: ', total_weight)
        
        self.capacity_exceeded_flag = capacity_exceeded_flag
        self.ic_violated_flag = ic_violated_flag
        self.combination_reused_team = combination_reused_team
        self.team_departments_weights = team_departments_weights
        self.team_locations = team_locations
        self.reused_teams = reused_teams
    
    def post_analysis(self):
        
        post_process_flag = self.post_process_flag
        
        if (post_process_flag):
                self.apply_post_processing()
        else:
            self.intermediate_selected_teams_indices = self.selected_teams_indices
            self.post_processing_ic_combination = []
            
        self.results_analysis()
            
    def apply_post_processing(self):

        max_state = self.max_state
        ic_violated_flag = self.ic_violated_flag
        combination_reused_team = self.combination_reused_team
        intrinsic_out_mwm = self.intrinsic_out_mwm
        reward_project = self.reward_project
        team_departments_weights = self.team_departments_weights
        selected_teams_indices = self.selected_teams_indices
        intrinsic_output = self.intrinsic_output
        available_combinations = self.available_combinations
        team_locations = self.team_locations
        reward_team = self.reward_team
        teams_types = self.teams_types
        capacity_exceeded_flag = self.capacity_exceeded_flag
        
        print('-------')
        post_processing_ic_combination = []
        if (ic_violated_flag):
            post_processing_ic_combination = [0]*len(combination_reused_team)
            for i in range(0, len(combination_reused_team)):
                min_intrinsic = 1000000000
                for j in range(0, len(combination_reused_team[i])):
                    combination = combination_reused_team[i][j]
                    if (intrinsic_out_mwm[combination[0], combination[1]]<min_intrinsic):
                        min_intrinsic = intrinsic_out_mwm[combination[0], combination[1]]
                        post_processing_ic_combination[i] = combination
                    print('Cost of ', combination, 'is: ', reward_project[combination[0], combination[1]])
                    print('Intrinsic Output RHS', intrinsic_out_mwm[combination[0], combination[1]])
            print('post_processing_IC_combination: ', post_processing_ic_combination)
        
        intermediate_team_department_weights = deepcopy(team_departments_weights)
        intermediate_selected_teams_indices = deepcopy(selected_teams_indices)
        
        if (capacity_exceeded_flag):
            print('-------')
            print('POST-PROCESSING STARTED')
            print('-------')
            while (any(intermediate_team_department_weights>max_state)):
                post_processing_departments = np.where(intermediate_team_department_weights>max_state)[0]
                print('Post-processing BINS: ', post_processing_departments)
                print('-------')
                if (len(post_processing_departments)):
                    intrinsic_output_selected = intrinsic_output[intermediate_selected_teams_indices]
                    processing_department_index = post_processing_departments[0]
                    print('Post-processing bin: ', processing_department_index)
                    packages_processing_bin_index = [i for i in intermediate_selected_teams_indices \
                        if (processing_department_index in np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices==i)[0][0]]])-1)]
                    print('max(reward_team[packages_processing_bin_index]: ', max(reward_team[packages_processing_bin_index]))
                    removed_team = intermediate_selected_teams_indices[np.where(intrinsic_output_selected==max(intrinsic_output_selected[intermediate_selected_teams_indices.searchsorted(packages_processing_bin_index)]))[0][0]]
                    print('reward[removed_team]: ', reward_team[removed_team])
                    removed_team_location = team_locations[np.where(intermediate_selected_teams_indices==removed_team)[0][0]]
                    removed_team_combination = np.array(available_combinations[removed_team_location])-1
                    removed_team_type = teams_types[removed_team]
                    print('Removing Package ', removed_team, 'of type ' ,removed_team_type , \
                    ' from bins ', removed_team_combination)
                    intermediate_selected_teams_indices = np.setdiff1d(intermediate_selected_teams_indices, removed_team)
                    team_locations.remove(removed_team_location)
                    u = removed_team_combination[0]; v = removed_team_combination[1]
                    if (removed_team_type==1):
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                    elif(removed_team_type==2):
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif(removed_team_type==3):
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                    elif(removed_team_type==4):
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif(removed_team_type==5):
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 2
        
            print('-------')
            print('POST-PROCESSING ENDED')
            print('-------')

        self.intermediate_selected_teams_indices = intermediate_selected_teams_indices
        self.post_processing_ic_combination = post_processing_ic_combination
        
    def results_analysis(self):

        intermediate_selected_teams_indices = self.intermediate_selected_teams_indices
        reused_teams = self.reused_teams
        post_process_flag = self.post_process_flag
        matching_array = self.matching_array
        post_processing_ic_combination = self.post_processing_ic_combination
        reward_team = self.reward_team
        teams_types = self.teams_types
        teams_weights_per_department = self.teams_weights_per_department
        team_locations = self.team_locations
        available_combinations = self.available_combinations
        num_projects = self.num_projects
        reward_project = self.reward_project
        max_state = self.max_state
        num_departments = self.num_departments
        
        results_composed = []
        project_matched_flag = False
        post_processing_capacity_exceeded_flag = False
        
        post_processing_package_bins_weights = np.zeros((num_departments), dtype = int)
        team_type = []
        total_value = 0
        project_array = []
        for i in range(0,len(intermediate_selected_teams_indices)):
            index = intermediate_selected_teams_indices[i]
            if (index in reused_teams and post_process_flag):
                project_index = np.where(matching_array[0]==[i[0] for i in post_processing_ic_combination if i[1]==index][0])[0]
            else:
                project_index = np.where(matching_array[1]==index)[0]
            project_array.append(project_index)
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            for j in range(0, num_departments):
                weight = teams_weights_per_department[j][index]
                if (weight !=0):
                    post_processing_package_bins_weights[j] = post_processing_package_bins_weights[j] + weight
            location_team = team_locations[i]
            print('Team ', index, ' has type ', teams_types[index], ' & is assigned to project ',matching_array[0][project_index] ,' & is in units ', np.array(available_combinations[location_team])-1)
            results_composed.append((available_combinations[location_team][0], available_combinations[location_team][1], teams_types[index], (matching_array[1][project_index]+1).tolist(), (matching_array[0][project_index]+1).tolist()))
        print('-------')
        for i in range(0, num_departments):
            print('Unit ', i, ' has a used capacity ', post_processing_package_bins_weights[i])
            print('Unit ', i, ' has a MAX capacity ', max_state[i] )
            if (post_processing_package_bins_weights[i] > max_state[i]):
                print('Unit ', i, ' needs post-processing')
                post_processing_capacity_exceeded_flag = True
        print('Post-Processing Exceeded Capacity Flag: ', post_processing_capacity_exceeded_flag) 

        print('-------')
        print('Types of Teams: ', np.unique(team_type, return_counts=True)[0], 'Counts: ', np.unique(team_type, return_counts=True)[1])
        print('-------------')

        print('-------')
        print('Total Value LHS: ', -total_value)
        print('Total number of selected teams: ', len(intermediate_selected_teams_indices), 'out of ', len(reward_team))
        print('Total number of Matched Projects: ', len(project_array), 'out of ', num_projects)

        total_weight = 0
        for i in range(0, len(matching_array[1])):
            team_index = matching_array[1][i]
            if (team_index in reused_teams and post_process_flag):
                project_index = [i[0] for i in post_processing_ic_combination if i[1]==team_index][0]
            else:
                project_index = matching_array[0][i]
            if (team_index in intermediate_selected_teams_indices):
                total_weight = total_weight + reward_project[project_index, team_index]

        print('MWM')
        print('-------------')

        if (len(project_array)==num_projects):
            project_matched_flag = True
            

        print(f'Matching All Targets to Packages : {project_matched_flag}')
        print('Total Weight RHS: ', total_weight)