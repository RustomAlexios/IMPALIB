# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *


class EqualityConstraint:
    def __init__(
        self,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        FILTERING_FLAG,
        ALPHA,
        USED_VARIABLES,
        VARIABLES_CONNECTIONS,
        CONSTRAINT_CONNECTIONS,
        INCOMING_METRICS_COST
    ):
        self.num_variables = NUM_VARIABLES
        self.num_constraints = NUM_CONSTRAINTS
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.used_variables = USED_VARIABLES
        self.variables_connections = VARIABLES_CONNECTIONS
        self.constraint_connections = CONSTRAINT_CONNECTIONS
        self.incoming_metrics_cost = INCOMING_METRICS_COST
    
    def variable_ec_to_ksat_constraint_update(self, ksat_constraint_to_eq_constraint_m):
        
        used_variables = self.used_variables
        self.variable_ec_to_ksat_constraint_m = np.zeros((self.num_constraints, self.num_variables), dtype=np_impa_lib)
        
        used_incoming_metrics_cost = np.zeros(self.num_variables)
        used_incoming_metrics_cost[used_variables] = self.incoming_metrics_cost[used_variables]
        
        sum_messages = np.sum(ksat_constraint_to_eq_constraint_m, axis=0) + used_incoming_metrics_cost
        
        for index_variable in range(len(used_variables)):
            variable = used_variables[index_variable]
            for constraint in self.variables_connections[variable]:
                self.variable_ec_to_ksat_constraint_m[constraint][variable] = sum_messages[variable] - ksat_constraint_to_eq_constraint_m[constraint][variable]
                
        return self.variable_ec_to_ksat_constraint_m