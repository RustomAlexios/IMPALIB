# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *

class OutputsImpa:
    def __init__(
        self,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        INCOMING_METRICS_COST
    ):
        self.num_variables = NUM_VARIABLES
        self.num_constraints = NUM_CONSTRAINTS
        self.incoming_metrics_cost = INCOMING_METRICS_COST
        
    def extrinsic_output_variable_ec_update(self, ksat_constraint_to_eq_constraint_m):
        return np.sum(ksat_constraint_to_eq_constraint_m, axis=0,)