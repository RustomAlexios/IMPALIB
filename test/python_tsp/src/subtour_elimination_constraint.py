# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *


class SubtourEliminationConstraint:
    def __init__(
        self,
        NUM_NODES,
        num_edge_variables,
        FILTERING_FLAG,
        ALPHA,
    ):
        self.num_nodes = NUM_NODES
        self.num_edge_variables = num_edge_variables
        # self.symmetric_flag = SYMMETRIC_FLAG
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA

    def subtour_constraints_to_edge_ec_update(
        self,
        edge_ec_to_subtour_constraints_m,
        delta_S_indices_list,
    ):
        # if (self.symmetric_flag):
        #    subtour_constraints_to_edge_ec_m = np.zeros((self.num_edge_variables, len(delta_S_indices_list)), dtype = np_impa_lib)
        #    for index_subtour_constraint in range(len(delta_S_indices_list)):
        #        delta_S_indices = delta_S_indices_list[index_subtour_constraint]
        #        for index_edge_variable in range(len(delta_S_indices)):
        #            remaining_connections = delta_S_indices[:index_edge_variable]+delta_S_indices[index_edge_variable+1:]
        #            subtour_constraints_to_edge_ec_m[delta_S_indices[index_edge_variable], index_subtour_constraint] = -np.maximum(np.sort(np.array(edge_ec_to_subtour_constraints_m).T[remaining_connections,index_subtour_constraint])[1], zero_value)
        #    self.subtour_constraints_to_edge_ec_m = subtour_constraints_to_edge_ec_m.T.tolist()
        # else:
        subtour_constraints_to_edge_ec_m = np.zeros(
            (
                self.num_edge_variables,
                len(delta_S_indices_list),
            ),
            dtype=np_impa_lib,
        )
        for index_subtour_constraint in range(len(delta_S_indices_list)):
            delta_S_indices = delta_S_indices_list[index_subtour_constraint]
            for index_edge_variable in range(len(delta_S_indices)):
                remaining_connections = delta_S_indices[:index_edge_variable] + delta_S_indices[index_edge_variable + 1 :]
                subtour_constraints_to_edge_ec_m[
                    delta_S_indices[index_edge_variable],
                    index_subtour_constraint,
                ] = -np.maximum(
                    min(
                        np.array(edge_ec_to_subtour_constraints_m).T[
                            remaining_connections,
                            index_subtour_constraint,
                        ]
                    ),
                    zero_value,
                )
        self.subtour_constraints_to_edge_ec_m = subtour_constraints_to_edge_ec_m.T.tolist()

        self.subtour_constraints_to_edge_ec_m_dummy = np.array(self.subtour_constraints_to_edge_ec_m)

        return self.subtour_constraints_to_edge_ec_m

    ### I NEED TO DO MORE TEST FOR THIS FILTERING FUNCTION
    def process_filtering(self, iter):
        alpha = self.alpha
        filtering_flag = self.filtering_flag

        if iter == 0 and filtering_flag and alpha != zero_value:
            subtour_constraints_to_edge_ec_m = (1 - alpha) * self.subtour_constraints_to_edge_ec_m_dummy
            self.subtour_constraints_to_edge_ec_m_old = deepcopy(subtour_constraints_to_edge_ec_m)
        elif iter > 0 and filtering_flag and alpha != 0:
            subtour_constraints_to_edge_ec_m = alpha * self.subtour_constraints_to_edge_ec_m_old + (1 - alpha) * self.subtour_constraints_to_edge_ec_m_dummy
            self.subtour_constraints_to_edge_ec_m_old = deepcopy(subtour_constraints_to_edge_ec_m)
        elif not filtering_flag:
            subtour_constraints_to_edge_ec_m = deepcopy(self.subtour_constraints_to_edge_ec_m_dummy)

        return subtour_constraints_to_edge_ec_m.tolist()
