# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *


class DegreeConstraint:
    # def __init__(self, NUM_NODES, num_edge_variables, edge_connections, SYMMETRIC_FLAG, FILTERING_FLAG, ALPHA):
    def __init__(
        self,
        NUM_NODES,
        num_edge_variables,
        edge_connections,
        FILTERING_FLAG,
        ALPHA,
    ):
        self.num_nodes = NUM_NODES
        self.num_edge_variables = num_edge_variables
        self.edge_connections = edge_connections
        # self.symmetric_flag = SYMMETRIC_FLAG
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA

    def degree_constraint_to_edge_ec_update(
        self,
        edge_ec_to_degree_constraint_m,
    ):
        degree_constraint_to_eq_constraint_m = np.zeros(
            (
                self.num_edge_variables,
                self.num_nodes,
            ),
            dtype=np_impa_lib,
        )
        for index_node in range(self.num_nodes):
            # if (self.symmetric_flag):
            #    connections = [i for i, pair in enumerate(self.edge_connections) if index_node in pair]
            #    for index_edge_variable in range(len(connections)):
            #        remaining_connections = connections[:index_edge_variable]+connections[index_edge_variable+1:]
            #        degree_constraint_to_eq_constraint_m[connections[index_edge_variable], index_node] = -np.sort(edge_ec_to_degree_constraint_m[remaining_connections,index_node])[1]
            # else:
            connections = [i for i, pair in enumerate(self.edge_connections) if index_node == pair[0]]
            for index_edge_variable in range(len(connections)):
                remaining_connections = connections[:index_edge_variable] + connections[index_edge_variable + 1 :]
                degree_constraint_to_eq_constraint_m[
                    connections[index_edge_variable],
                    index_node,
                ] = -np.min(
                    edge_ec_to_degree_constraint_m[
                        remaining_connections,
                        index_node,
                    ]
                )
            connections = [i for i, pair in enumerate(self.edge_connections) if index_node == pair[1]]
            for index_edge_variable in range(len(connections)):
                remaining_connections = connections[:index_edge_variable] + connections[index_edge_variable + 1 :]
                degree_constraint_to_eq_constraint_m[
                    connections[index_edge_variable],
                    index_node,
                ] = -np.min(
                    edge_ec_to_degree_constraint_m[
                        remaining_connections,
                        index_node,
                    ]
                )

        self.degree_constraint_to_eq_constraint_m_dummy = deepcopy(degree_constraint_to_eq_constraint_m)

        # return degree_constraint_to_eq_constraint_m

    def process_filtering(self, iter):
        alpha = self.alpha
        filtering_flag = self.filtering_flag

        if iter == 0 and filtering_flag and alpha != zero_value:
            degree_constraint_to_eq_constraint_m = (1 - alpha) * self.degree_constraint_to_eq_constraint_m_dummy
            self.degree_constraint_to_eq_constraint_m_old = deepcopy(degree_constraint_to_eq_constraint_m)
        elif iter > 0 and filtering_flag and alpha != 0:
            degree_constraint_to_eq_constraint_m = alpha * self.degree_constraint_to_eq_constraint_m_old + (1 - alpha) * self.degree_constraint_to_eq_constraint_m_dummy
            self.degree_constraint_to_eq_constraint_m_old = deepcopy(degree_constraint_to_eq_constraint_m)
        elif not filtering_flag:
            degree_constraint_to_eq_constraint_m = deepcopy(self.degree_constraint_to_eq_constraint_m_dummy)

        self.degree_constraint_to_eq_constraint_m = degree_constraint_to_eq_constraint_m
