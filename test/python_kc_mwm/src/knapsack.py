# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *
from environmentModule import zero_value

class Knapsack:
    
    def __init__(self, NUM_DEPARTMENTS, NUM_TEAMS, FILTERING_FLAG, ALPHA, reward_team):
        self.num_departments = NUM_DEPARTMENTS
        self.num_teams = NUM_TEAMS
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        
        self.reward_team = reward_team
        
        self.extrinsic_output_department_dummy = np.zeros((self.num_departments, self.num_teams), dtype = np_impa_lib)
        self.extrinsic_output_department = np.zeros((self.num_departments, self.num_teams), dtype = np_impa_lib)
        all_departments  = np.array(range(0,self.num_departments))
        self.all_departments = all_departments
        
    def forward(self, max_state_department, transition_model, weight_vector):

        forward_messages_initial_stage = inf*np.ones(max_state_department, dtype = np_impa_lib)

        forward_messages_initial_stage = np.insert(forward_messages_initial_stage, 0, zero_value)

        initial_forward_messages = deepcopy(forward_messages_initial_stage)
        
        non_zero_weight_department = [i for i, e in enumerate(weight_vector) if e != 0]
        
        number_of_stages = len(transition_model)

        stage_forward_messages = []
        
        if (non_zero_weight_department[0] != 0):
            stage_forward_messages.extend([initial_forward_messages.tolist() for i in range(non_zero_weight_department[0])])
        
        for t in non_zero_weight_department:
            forward_messages = []
            for a in range(0,max_state_department+1):
                if (a - weight_vector[t] >= 0 and (not (t== 0 and a != weight_vector[t]))):
                    min_forward = min(initial_forward_messages[a - weight_vector[t]]+transition_model[t], \
                        initial_forward_messages[a])
                    forward_messages.append(min_forward)
                else:
                    forward_messages.append(initial_forward_messages[a])                   
            initial_forward_messages = forward_messages
            stage_forward_messages.append(initial_forward_messages)
            if (len(non_zero_weight_department) != number_of_stages):
                if (t+1 not in non_zero_weight_department):
                    if (t+1 >= non_zero_weight_department[-1] and t+1 <number_of_stages):
                        stage_forward_messages.extend([initial_forward_messages for i in range(number_of_stages-(t+1))])
                    elif (t+1 < non_zero_weight_department[-1]):
                        next_index = next(x for x, val in enumerate(non_zero_weight_department) if val > t)
                        stage_forward_messages.extend([initial_forward_messages for i in range(non_zero_weight_department[next_index]-t-1)])
        returned_list = list(np_impa_lib([forward_messages_initial_stage.tolist()] + stage_forward_messages))
        self.stage_forward_messages = returned_list
        
    def backward(self, max_state_department, transition_model, weight_vector):
        
        backward_messages_initial_stage = np.zeros(max_state_department+1, dtype = np_impa_lib)
        
        initial_backward_messages = deepcopy(backward_messages_initial_stage)
        
        non_zero_weight_department = [i for i, e in enumerate(weight_vector) if e != 0]
        
        reversed_non_zero_weight_bin = non_zero_weight_department[::-1]

        number_of_stages = len(transition_model)

        stage_backward_messages = []

        if (reversed_non_zero_weight_bin[0] != number_of_stages-1):
            stage_backward_messages.extend([initial_backward_messages.tolist() for i in range(number_of_stages -1 - reversed_non_zero_weight_bin[0])])

        for t in reversed_non_zero_weight_bin:
        
            backward_messages = []

            for a in range(0,max_state_department+1):
                if ((t>0 and a + weight_vector[t] <= max_state_department) or (a==0 and t==0)):
                    min_backward = min(initial_backward_messages[a], initial_backward_messages[a + weight_vector[t]]+transition_model[t])
                    backward_messages.append(min_backward)
                else:
                    backward_messages.append(initial_backward_messages[a])
                
            initial_backward_messages = backward_messages
            stage_backward_messages.append(initial_backward_messages)

            if (len(non_zero_weight_department) != number_of_stages):
                if (t-1 not in non_zero_weight_department):
                    if (t-1 <= reversed_non_zero_weight_bin[-1] and t-1 >=0):
                        stage_backward_messages.extend([initial_backward_messages for i in range(reversed_non_zero_weight_bin[-1])])
                    elif (t-1 > reversed_non_zero_weight_bin[-1]):
                        old_index = next(x for x, val in enumerate(non_zero_weight_department) if val > t-1)
                        stage_backward_messages.extend([initial_backward_messages for i in range(t-non_zero_weight_department[old_index -1]-1)]) 
        stage_backward_messages = stage_backward_messages[::-1]
        
        returned_list = list(np_impa_lib(stage_backward_messages + [backward_messages_initial_stage.tolist()]))
        self.stage_backward_messages = returned_list

    def forward_backward(self, max_state_department, transition_model, weight_vector):
        
        self.forward(max_state_department, transition_model, weight_vector)
        
        self.backward(max_state_department, transition_model, weight_vector)

    
    def extrinsic_output_department_lhs(self, department_index, max_state_department, transition_model, weight_vector):
        
        fv = self.stage_forward_messages
        bv = self.stage_backward_messages
        
        num_teams = self.num_teams
        
        non_zero_weight_department = [i for i, e in enumerate(weight_vector) if e != 0]

        extrinsic_output = np.zeros(num_teams, dtype = np_impa_lib)
        
        for team_index in non_zero_weight_department:
                metric_path_solid = []
                metric_path_dash = []
                if (team_index == 0):
                    metric_path_solid.append(fv[team_index][team_index]+bv[team_index+1][weight_vector[team_index]]+transition_model[team_index])
                    metric_path_dash.append(fv[team_index][team_index]+bv[team_index+1][0])
                else:
                    for i in range(0,max_state_department+1):
                        if (i + weight_vector[team_index]<=max_state_department):
                            metric_path_solid.append(fv[team_index][i]+bv[team_index+1][i+weight_vector[team_index]]+transition_model[team_index])
                        metric_path_dash.append(fv[team_index][i]+bv[team_index+1][i])
                extrinsic_output[team_index] = (min(metric_path_solid) - min(metric_path_dash) - transition_model[team_index])
        
        self.extrinsic_output_department_dummy[department_index] = extrinsic_output
        
    def process_extrinsic_output_department(self, iter):
        
        alpha = self.alpha
        filtering_flag = self.filtering_flag
        
        if (iter==0 and filtering_flag and alpha !=zero_value):
            extrinsic_output_department = (1-alpha)*self.extrinsic_output_department_dummy
            self.extrinsic_output_department_old = deepcopy(extrinsic_output_department)
        elif (iter>0 and filtering_flag and alpha !=0):
            extrinsic_output_department = alpha*self.extrinsic_output_department_old + (1-alpha)*self.extrinsic_output_department_dummy
            self.extrinsic_output_department_old = deepcopy(extrinsic_output_department)
        elif (not filtering_flag):
            extrinsic_output_department = deepcopy(self.extrinsic_output_department_dummy)
            
        self.extrinsic_output_department = extrinsic_output_department
        
        
    def team_to_department_update(self, items_weights_per_bin, oric_to_package_m):
        
        num_departments = self.num_departments
        num_teams = self.num_teams
        all_departments = self.all_departments
        
        team_to_knapsack_m = np.zeros((num_departments,num_teams), dtype = np_impa_lib)
        
        for department_index in range(0, num_departments):
            remaining_inputs = np.setdiff1d(all_departments, department_index)
            total_intersection = []
            non_zero_weight_department = [i for i, e in enumerate(items_weights_per_bin[department_index]) if e != 0]
            for u in remaining_inputs:
                intersection = np.intersect1d(np.nonzero(items_weights_per_bin[department_index])[0], np.nonzero(items_weights_per_bin[u])[0])
                total_intersection.extend(intersection)
                team_to_knapsack_m[department_index][intersection] = self.reward_team[intersection] + self.extrinsic_output_department[u][intersection] + oric_to_package_m[intersection]
            unique_edge_bin = np.setdiff1d(non_zero_weight_department, total_intersection)
            team_to_knapsack_m[department_index][unique_edge_bin] = self.reward_team[unique_edge_bin] + oric_to_package_m[unique_edge_bin]
        
        return team_to_knapsack_m