# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *


class KSatConstraint:
    def __init__(
        self,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        FILTERING_FLAG,
        ALPHA,
        CONSTRAINT_CONNECTIONS,
        CONSTRAINT_CONNECTIONS_TYPE, 
    ):
        self.num_variables = NUM_VARIABLES
        self.num_constraints = NUM_CONSTRAINTS
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.constraint_connections = CONSTRAINT_CONNECTIONS
        self.constraint_connections_type = CONSTRAINT_CONNECTIONS_TYPE
    
    def ksat_constraint_to_variable_ec_update(self, variable_ec_to_ksat_constraint_m):
        
        ksat_constraint_to_eq_constraint_m = np.zeros((self.num_constraints, self.num_variables,), dtype=np_impa_lib,)
        
        for index_constraint in range(self.num_constraints):
            connections = self.constraint_connections[index_constraint]
            connections_type = self.constraint_connections_type[index_constraint]
            for index_variable, variable in enumerate(connections):
                connections_solid = [node for i, node in enumerate(connections) if connections_type[i] == 1 and node !=variable]
                connections_dashed = [node for i, node in enumerate(connections) if connections_type[i] == -1 and node !=variable] 
                if (not len(connections_solid)):
                    optimum_solid_messages = np.where(connections_type[index_variable] == 1,
                                                    -np_impa_lib(np.inf),
                                                    np_impa_lib(np.inf))
                else:
                    optimum_solid_messages = np.where(connections_type[index_variable] == 1,
                                                    np.max(-variable_ec_to_ksat_constraint_m[index_constraint, connections_solid]),
                                                    np.min(variable_ec_to_ksat_constraint_m[index_constraint, connections_solid]))
                if (not len(connections_dashed)):
                    optimum_dashed_messages = np.where(connections_type[index_variable] == 1,
                                                    -np_impa_lib(np.inf),
                                                    np_impa_lib(np.inf))
                else:
                    optimum_dashed_messages = np.where(connections_type[index_variable] == 1,
                                                    np.max(variable_ec_to_ksat_constraint_m[index_constraint, connections_dashed]),
                                                    np.min(-variable_ec_to_ksat_constraint_m[index_constraint, connections_dashed]))
                
                if (connections_type[index_variable] == 1):
                    ksat_constraint_to_eq_constraint_m[index_constraint, variable] = min(zero_value, max(optimum_solid_messages, optimum_dashed_messages))   
                else:
                    ksat_constraint_to_eq_constraint_m[index_constraint, variable] = max(zero_value, min(optimum_solid_messages, optimum_dashed_messages))

        self.ksat_constraint_to_eq_constraint_m_dummy = deepcopy(ksat_constraint_to_eq_constraint_m)
        
    def process_filtering(self, iter):
        alpha = self.alpha
        filtering_flag = self.filtering_flag

        if iter == 0 and filtering_flag and alpha != zero_value:
            ksat_constraint_to_eq_constraint_m = (1 - alpha) * self.ksat_constraint_to_eq_constraint_m_dummy
            self.ksat_constraint_to_eq_constraint_m_old = deepcopy(ksat_constraint_to_eq_constraint_m)
        elif iter > 0 and filtering_flag and alpha != 0:
            ksat_constraint_to_eq_constraint_m = alpha * self.ksat_constraint_to_eq_constraint_m_old + (1 - alpha) * self.ksat_constraint_to_eq_constraint_m_dummy
            self.ksat_constraint_to_eq_constraint_m_old = deepcopy(ksat_constraint_to_eq_constraint_m)
        elif not filtering_flag:
            ksat_constraint_to_eq_constraint_m = deepcopy(self.ksat_constraint_to_eq_constraint_m_dummy)

        self.ksat_constraint_to_eq_constraint_m = ksat_constraint_to_eq_constraint_m
        