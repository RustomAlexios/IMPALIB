# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *


class EqualityConstraint:
    def __init__(
        self,
        NUM_NODES,
        num_edge_variables,
        cost_matrix,
        edge_connections,
        edge_degree_constraint_cost,
        AUGMENTATION_FLAG,
        FILTERING_FLAG,
        ALPHA,
    ):
        self.num_nodes = NUM_NODES
        self.num_edge_variables = num_edge_variables
        self.cost_matrix = cost_matrix
        self.edge_degree_constraint_cost = edge_degree_constraint_cost
        self.augmentation_flag = AUGMENTATION_FLAG
        self.edge_connections = edge_connections
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA

    def edge_ec_to_degree_constraint_relaxed_graph_update(
        self,
        degree_constraint_to_eq_constraint_m,
    ):
        flipped_degree_constraint_to_eq_constraint_m = self.flip_matrix(degree_constraint_to_eq_constraint_m)
        self.edge_ec_to_degree_constraint_m = flipped_degree_constraint_to_eq_constraint_m + self.edge_degree_constraint_cost
        return self.edge_ec_to_degree_constraint_m

    def edge_ec_to_subtour_constraints_update(
        self,
        degree_constraint_to_eq_constraint_m,
        subtour_constraints_to_edge_ec_m,
        delta_S_indices_list,
        cost_edge_variable,
    ):
        all_subtour_constraints = np.array(range(len(delta_S_indices_list)))
        
        if len(delta_S_indices_list) == 1:
            edge_ec_to_subtour_constraints_m = np.zeros(
                (self.num_edge_variables),
                dtype=np_impa_lib,
            )
            edge_ec_to_subtour_constraints_m[delta_S_indices_list[0]] = (
                np.sum(
                    degree_constraint_to_eq_constraint_m,
                    axis=1,
                )[delta_S_indices_list[0]]
                + cost_edge_variable[delta_S_indices_list[0]]
            )
            self.edge_ec_to_subtour_constraints_m = [edge_ec_to_subtour_constraints_m.tolist()]
        else:
            self.edge_ec_to_subtour_constraints_m = []
            for index_subtour_constraint in range(len(delta_S_indices_list)):
                remaining_subtour_constraints_indices = np.setdiff1d(
                    all_subtour_constraints,
                    index_subtour_constraint,
                )
                sum_degree_constraint_to_eq_constraint_m = np.zeros(self.num_edge_variables)
                sum_degree_constraint_to_eq_constraint_m[delta_S_indices_list[index_subtour_constraint]] = np.sum(
                    degree_constraint_to_eq_constraint_m[
                        delta_S_indices_list[index_subtour_constraint],
                        :,
                    ],
                    axis=1,
                )
                subtour_cost_edge_variable = np.zeros(self.num_edge_variables)
                subtour_cost_edge_variable[delta_S_indices_list[index_subtour_constraint]] = cost_edge_variable[delta_S_indices_list[index_subtour_constraint]]
                sum_subtours = np.zeros(self.num_edge_variables)
                sum_subtours[delta_S_indices_list[index_subtour_constraint]] = np.sum(
                    np.array(subtour_constraints_to_edge_ec_m).T[
                        :,
                        remaining_subtour_constraints_indices,
                    ],
                    axis=1,
                )[delta_S_indices_list[index_subtour_constraint]]
                self.edge_ec_to_subtour_constraints_m.append((sum_degree_constraint_to_eq_constraint_m + sum_subtours + subtour_cost_edge_variable).tolist())
                # self.edge_ec_to_subtour_constraints_m.append((sum_degree_constraint_to_eq_constraint_m + np.sum(np.array(subtour_constraints_to_edge_ec_m).T[:, remaining_subtour_constraints_indices], axis=1) + subtour_cost_edge_variable).tolist())
                # self.edge_ec_to_subtour_constraints_m.append((np.sum(degree_constraint_to_eq_constraint_m, axis=1) + np.sum(np.array(subtour_constraints_to_edge_ec_m).T[:, remaining_subtour_constraints_indices], axis=1) + cost_edge_variable).tolist())
        
        return self.edge_ec_to_subtour_constraints_m

    def edge_ec_to_degree_constraint_augmented_graph_update(
        self,
        degree_constraint_to_eq_constraint_m,
        subtour_constraints_to_edge_ec_m,
    ):
        combined_subtour_constraints_to_edge_ec_m = np.sum(
            np.array(subtour_constraints_to_edge_ec_m).T,
            axis=1,
        )
        flipped_degree_constraint_to_eq_constraint_m = self.flip_matrix(degree_constraint_to_eq_constraint_m)
        reshaped_combined_subtour_constraints_to_edge_ec_m = np.zeros(flipped_degree_constraint_to_eq_constraint_m.shape)
        edge_connections = self.edge_connections
        
        for i in range(len(edge_connections)):
            connection = edge_connections[i]
            reshaped_combined_subtour_constraints_to_edge_ec_m[i][connection[0]] = combined_subtour_constraints_to_edge_ec_m[i]
            reshaped_combined_subtour_constraints_to_edge_ec_m[i][connection[1]] = combined_subtour_constraints_to_edge_ec_m[i]
        
        self.edge_ec_to_degree_constraint_m = reshaped_combined_subtour_constraints_to_edge_ec_m + flipped_degree_constraint_to_eq_constraint_m + self.edge_degree_constraint_cost
        
        return self.edge_ec_to_degree_constraint_m

    def flip_matrix(self, matrix):
        flipped_matrix = np.copy(matrix)
        
        for ele in range(len(self.edge_connections)):
            (row,col,) = self.edge_connections[ele]
            flipped_matrix[ele,row,] = matrix[ele,col,]
            flipped_matrix[ele,col,] = matrix[ele,row,]
        
        return flipped_matrix
