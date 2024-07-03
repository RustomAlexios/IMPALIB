# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impalib import *

class GraphicalModel:
    def __init__(
        self,
        NUM_ITERATIONS,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        K_VARIABLE,
        EXACT_SOLVER_FLAG,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
        POST_PROCESS_FLAG
    ):
        self.num_iterations = NUM_ITERATIONS
        self.num_variables = NUM_VARIABLES
        self.num_constraints = NUM_CONSTRAINTS
        self.k_variable = K_VARIABLE
        self.threshold = THRESHOLD
        self.exact_solver_flag = EXACT_SOLVER_FLAG
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.random_test_flag = RANDOM_TEST_FLAG

    def initialize(self):
        
        self.results_composed = []

        if not self.random_test_flag:
            self.num_variables = self.input_load[0]
            self.num_constraints = self.input_load[1]
            self.k_variable = self.input_load[2]

        num_variables = self.num_variables
        num_constraints = self.num_constraints
        k_variable = self.k_variable

        constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, truth_vector = self.create_incoming_metrics_cost(num_variables, num_constraints, k_variable)

        if not self.random_test_flag:
            self.constraints_connections = self.input_load[3]
            self.constraints_connections_type = self.input_load[4]
            self.incoming_metrics_cost = self.input_load[5]
            self.used_variables = self.input_load[6]
            self.variables_connections = self.input_load[7]
            self.variables_connections_type = self.input_load[8]
            constraints_connections = self.constraints_connections
            constraints_connections_type = self.constraints_connections_type
            incoming_metrics_cost = self.incoming_metrics_cost
            used_variables = self.used_variables
            variables_connections = self.variables_connections
            variables_connections_type = self.variables_connections_type

        self.constraints_connections = constraints_connections
        self.constraints_connections_type = constraints_connections_type
        self.incoming_metrics_cost = incoming_metrics_cost
        self.variables_connections = variables_connections
        self.variables_connections_type = variables_connections_type
        self.truth_vector = truth_vector
        
        self.used_variables = used_variables

        self.intrinsic_out_variable_ec = np.zeros(
            num_variables,
            dtype=np_impa_lib,
        )

        variable_ec_to_ksat_constraint_m = np.zeros(
            (
                num_constraints,
                num_variables,
            ),
            dtype=np_impa_lib,
        )

        for i in range(len(constraints_connections)):
            connection = constraints_connections[i]
            for variable in connection:
                cost = incoming_metrics_cost[variable]
                variable_ec_to_ksat_constraint_m[i][variable] = cost
                variable_ec_to_ksat_constraint_m[i][variable] = cost
        
        self.variable_ec_to_ksat_constraint_m = variable_ec_to_ksat_constraint_m
        
        self.model_eq_constraint = EqualityConstraint(
            num_variables,
            num_constraints,
            self.filtering_flag,
            self.alpha,
            self.used_variables,
            self.variables_connections,
            self.constraints_connections,
            self.incoming_metrics_cost
        )
        
        self.model_ksat_constraint = KSatConstraint(
            num_variables,
            num_constraints,
            self.filtering_flag,
            self.alpha,
            self.constraints_connections,
            self.constraints_connections_type, 
        )
        self.outputs = OutputsImpa(
            num_variables,
            num_constraints,
            self.incoming_metrics_cost
        )

    def create_incoming_metrics_cost(self, num_variables, num_constraints, k_variable):
        
        constraints_connections = [[] for _ in range(num_constraints)]
        constraints_connections_type = [[] for _ in range(num_constraints)]
        incoming_metrics_cost = np.zeros(num_variables)
        variables_connections = [[] for _ in range(num_variables)] 
        variables_connections_type = [[] for _ in range(num_variables)]
        used_variables = []
        
        truth_vector = np.random.randint(2, size=num_variables, dtype=int)
        
        normal_variance = 3
       
        for i, truth_value in enumerate(truth_vector):
            mean = -normal_variance/2 if truth_value == 1 else normal_variance/2
            incoming_metrics_cost[i] = np.random.normal(loc=mean, scale=np.sqrt(normal_variance))
        
        for constraint_idx in range(num_constraints):
            connections = []
            connections_types = []
            
            satisfying_variable_index = random.randint(0, num_variables - 1)
            satisfying_variable_value = truth_vector[satisfying_variable_index]
            connection_type = -1 if satisfying_variable_value == 0 else 1
            connections.append(satisfying_variable_index)
            connections_types.append(connection_type)
            used_variables.append(satisfying_variable_index)
            variables_connections[satisfying_variable_index].append(constraint_idx)
            variables_connections_type[satisfying_variable_index].append(connection_type)
            
            for k_variable_idx in range(k_variable-1):
                choice_variable = random.choice([idx for idx in range(num_variables) if idx not in connections])
                connections.append(choice_variable)
                choice_connection_type = random.choice([-1, 1])
                connections_types.append(choice_connection_type)
                
                variables_connections[choice_variable].append(constraint_idx)
                used_variables.append(choice_variable)
                variables_connections_type[choice_variable].append(choice_connection_type)

            constraints_connections[constraint_idx] = connections
            constraints_connections_type[constraint_idx] = connections_types 

        used_variables = list(set(used_variables))
        
        return constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, truth_vector
    
    def run_impa(self):
        
        extrinsic_output_variable_ec = np.zeros(self.num_variables, dtype = np_impa_lib)
        
        is_within_tolerance_count = 0
        
        for iter in range(0, self.num_iterations,):
            
            extrinsic_output_variable_ec_old = deepcopy(extrinsic_output_variable_ec)
            
            self.model_ksat_constraint.ksat_constraint_to_variable_ec_update(self.variable_ec_to_ksat_constraint_m)
            
            self.model_ksat_constraint.process_filtering(iter)
            
            self.variable_ec_to_ksat_constraint_m = self.model_eq_constraint.variable_ec_to_ksat_constraint_update(self.model_ksat_constraint.ksat_constraint_to_eq_constraint_m)
            
            extrinsic_output_variable_ec = self.outputs.extrinsic_output_variable_ec_update(self.model_ksat_constraint.ksat_constraint_to_eq_constraint_m)
       
        used_incoming_metrics_cost = np.zeros(self.num_variables)
        
        used_incoming_metrics_cost[self.used_variables] = self.incoming_metrics_cost[self.used_variables]
        
        self.extrinsic_output_variable_ec = extrinsic_output_variable_ec
        
        intrinsic_out_variable_ec = extrinsic_output_variable_ec + used_incoming_metrics_cost
        
        self.intrinsic_out_variable_ec = intrinsic_out_variable_ec
        
        self.hard_decision_analysis()
        
    def hard_decision_analysis(self,):
        
        intrinsic_out_variable_ec = self.intrinsic_out_variable_ec

        hard_decision = np.array(deepcopy(intrinsic_out_variable_ec))
        hard_decision[intrinsic_out_variable_ec > self.threshold] = 0
        hard_decision[intrinsic_out_variable_ec <= self.threshold] = 1
        
        self.hard_decision = hard_decision
        
        self.active_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if hard_decision[self.used_variables][i]]
        self.inactive_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if not hard_decision[self.used_variables][i]]
        
        print('--------')
        print("Active variables: ",self.active_variables,)
        print("Inactive variables: ", self.inactive_variables,)
        print('--------')
        
        indices_satisfied_constraints_solid = [(i, j) for i, sublist in enumerate(self.variables_connections_type) for j, element in enumerate(sublist) if element == 1 and i in self.active_variables]
        satisfied_constraints_solid = list(set([self.variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_solid]))
        
        indices_satisfied_constraints_dashed = [(i, j) for i, sublist in enumerate(self.variables_connections_type) for j, element in enumerate(sublist) if element == -1 and i in self.inactive_variables]
        satisfied_constraints_dashed = list(set([self.variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_dashed])) 
        
        satisfied_constraints = list(set(satisfied_constraints_solid + satisfied_constraints_dashed))
        
        unsatisfied_constraints = [constraint for constraint in range(self.num_constraints) if constraint not in satisfied_constraints]
        
        if (len(unsatisfied_constraints)):
            print(f"{len(satisfied_constraints)} satisfied constraints: {satisfied_constraints}")
            print(f"{len(unsatisfied_constraints)} unsatisfied constraints: {unsatisfied_constraints}")
        else:
            print(f"All {self.num_constraints} constraints are satisfied")