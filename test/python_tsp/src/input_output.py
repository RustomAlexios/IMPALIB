# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *

class OutputsImpa:
    def __init__(self, NUM_NODES, num_edge_variables, cost_matrix, augmentation_flag):
        self.num_nodes = NUM_NODES
        self.num_edge_variables = num_edge_variables
        self.cost_matrix = cost_matrix
        self.augmentation_flag = augmentation_flag
        
    def extrinsic_output_edge_ec_relaxed_graph_update(self, degree_constraint_to_eq_constraint_m):
        return np.sum(degree_constraint_to_eq_constraint_m,axis=1)
    
    def extrinsic_output_edge_ec_augmented_graph_update(self, degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m):
        return np.sum(degree_constraint_to_eq_constraint_m,axis=1) + np.sum(np.array(subtour_constraints_to_edge_ec_m).T, axis=1)