# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impa.environmentModule import (
    np,
    math,
    np_impa_lib,
    zero_value,
    deepcopy,
    combinations,
    itertools,
    elkai,
    solve_tsp_brute_force,
    solve_tsp_simulated_annealing,
    islice,
    copy,
    time,
    permutations,
    pkl,
    chain,
)
from impa.initializationModule import (
    fighters_count,
    weasels_count,
    F_u,
    W_u,
    distance_metric,
    distance_metric_rendezvouz,
    K_1,
)
from impa.cFunctionAPI import (
    c_int_p,
    c_impa_lib_type_p,
    c_bool_p,
    WrapperTsp,
    WrapperKcMwm,
)


class GraphicalModelTsp:
    def __init__(
        self,
        NUM_ITERATIONS,
        NUM_NODES,
        SYMMETRIC_FLAG,
        AUGMENTATION_FLAG,
        EXACT_SOLVER_FLAG,
        LKH_SOLVER_FLAG,
        SIM_ANNEALING_FLAG,
        RESET_FLAG,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
        MAX_COUNT,
        POST_PROCESS_FLAG,
        K_OPT_FLAG,
        MAX_AUGM_COUNT,
    ):
        self.num_iterations = NUM_ITERATIONS
        self.num_nodes = NUM_NODES
        self.threshold = THRESHOLD
        self.symmetric_flag = SYMMETRIC_FLAG
        self.augmentation_flag = AUGMENTATION_FLAG
        self.exact_solver_flag = EXACT_SOLVER_FLAG
        self.lkh_solver_flag = LKH_SOLVER_FLAG
        self.sim_annealing_flag = SIM_ANNEALING_FLAG
        self.reset_flag = RESET_FLAG
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.random_test_flag = RANDOM_TEST_FLAG
        self.max_count = MAX_COUNT
        self.post_process_flag = POST_PROCESS_FLAG
        self.k_opt_flag = K_OPT_FLAG
        self.max_augm_count = MAX_AUGM_COUNT

    def initialize(self):
        self.results_composed = []

        self.num_augmentations = 0
        self.num_added_constraints = 0
        self.no_improvement_sol_count_exceeded_flag = False
        self.no_consecutive_loops_count_exceeded_flag = False
        self.sol_oscillation_count_exceeded_flag = False

        if not self.random_test_flag:
            self.num_nodes = self.input_load[0]
            self.symmetric_flag = self.input_load[1]

        num_nodes = self.num_nodes
        symmetric_flag = self.symmetric_flag

        print(f"num_nodes: {self.num_nodes}")
        print(f"symmetric_flag: {self.symmetric_flag}")

        cost_matrix, edge_connections = self.create_cost_matrix(symmetric_flag)
        self.tour_impa_flag = False
        self.tour_impa = np.zeros(num_nodes + 1, dtype=np.int32)  # +1 to account for returning edge
        self.cost_impa = np.zeros(1, dtype=np_impa_lib)
        self.selected_edges = np.zeros(math.comb(self.num_nodes, 2), dtype=np.int32)
        self.subtour_paths = np.zeros(self.num_nodes**2, dtype=np.int32)
        self.open_paths = np.zeros(self.num_nodes**2, dtype=np.int32)

        if not self.random_test_flag:
            self.cost_matrix = self.input_load[2]
            self.edge_connections = self.input_load[3]
            cost_matrix = self.cost_matrix
            edge_connections = self.edge_connections

        self.cost_matrix = cost_matrix
        self.edge_connections = edge_connections

        num_edge_variables = num_nodes**2 - num_nodes

        self.num_edge_variables = num_edge_variables

        self.intrinsic_out_edge_ec = np.zeros(num_edge_variables, dtype=np_impa_lib)

        cost_edge_variable = np.zeros(num_edge_variables, dtype=np_impa_lib)

        edge_ec_to_degree_constraint_m = np.zeros((num_edge_variables, num_nodes), dtype=np_impa_lib)

        for i in range(len(edge_connections)):
            connection = edge_connections[i]
            cost = cost_matrix[connection[0], connection[1]]
            cost_edge_variable[i] = cost
            edge_ec_to_degree_constraint_m[i][connection[0]] = cost
            edge_ec_to_degree_constraint_m[i][connection[1]] = cost

        self.cost_edge_variable = cost_edge_variable
        self.edge_degree_constraint_cost = edge_ec_to_degree_constraint_m

        self.edge_ec_to_degree_constraint_m = edge_ec_to_degree_constraint_m

    def create_cost_matrix(self, symmetric_flag):
        n = self.num_nodes
        upper_triangle = np.random.uniform(1, 1000, size=(n, n))
        if symmetric_flag:
            matrix = np.triu(upper_triangle) + np.triu(upper_triangle, k=1).T
            np.fill_diagonal(matrix, zero_value)
        else:
            matrix = upper_triangle.copy()
            np.fill_diagonal(matrix, zero_value)
        edge_connections = np.transpose(np.nonzero(matrix))
        return matrix, edge_connections

    def run_impa(self):
        (
            edge_connection_flatten_p,
            cost_edge_variable_flatten_p,
            cost_matrix_flatten_p,
            edge_ec_to_degree_constraint_m_flatten_p,
            edge_degree_constraint_cost_flatten_p,
            extrinsic_output_edge_ec_p,
            num_augmentations_p,
            tour_impa_flatten_p,
            cost_impa_p,
            selected_edges_p,
            selected_edges_size_p,
            no_improvement_sol_count_exceeded_flag_p,
            no_consecutive_loops_count_exceeded_flag_p,
            sol_oscillation_count_exceeded_flag_p,
            subtour_paths_flatten_p,
            subtour_paths_size_flatten_p,
            num_added_constraints_p,
        ) = self.process_inputs_ctypes()

        WrapperTsp(
            np.int32(self.num_iterations),
            np.int32(self.num_nodes),
            np.int32(self.num_edge_variables),
            self.augmentation_flag,
            self.reset_flag,
            np_impa_lib(self.alpha),
            self.filtering_flag,
            edge_connection_flatten_p,
            cost_edge_variable_flatten_p,
            cost_matrix_flatten_p,
            edge_ec_to_degree_constraint_m_flatten_p,
            edge_degree_constraint_cost_flatten_p,
            extrinsic_output_edge_ec_p,
            np_impa_lib(self.threshold),
            num_augmentations_p,
            tour_impa_flatten_p,
            cost_impa_p,
            no_improvement_sol_count_exceeded_flag_p,
            no_consecutive_loops_count_exceeded_flag_p,
            selected_edges_p,
            selected_edges_size_p,
            sol_oscillation_count_exceeded_flag_p,
            subtour_paths_flatten_p,
            subtour_paths_size_flatten_p,
            np.int32(self.max_count),
            num_added_constraints_p,
            np.int32(self.max_augm_count),
        )

        self.process_outputs_ctypes(
            sol_oscillation_count_exceeded_flag_p,
            selected_edges_size_p,
            selected_edges_p,
            subtour_paths_size_flatten_p,
            subtour_paths_flatten_p,
            no_improvement_sol_count_exceeded_flag_p,
            no_consecutive_loops_count_exceeded_flag_p,
            tour_impa_flatten_p,
            num_augmentations_p,
            cost_impa_p,
            extrinsic_output_edge_ec_p,
            num_added_constraints_p,
        )

    def process_inputs_ctypes(self):
        edge_connection_flatten = self.edge_connections.flatten().astype(np.int32)
        edge_connection_flatten_p = edge_connection_flatten.ctypes.data_as(c_int_p)
        cost_edge_variable_flatten = self.cost_edge_variable.flatten().astype(np_impa_lib)
        cost_edge_variable_flatten_p = cost_edge_variable_flatten.ctypes.data_as(c_impa_lib_type_p)
        cost_matrix_flatten = self.cost_matrix.flatten().astype(np_impa_lib)
        cost_matrix_flatten_p = cost_matrix_flatten.ctypes.data_as(c_impa_lib_type_p)
        edge_ec_to_degree_constraint_m_flatten = self.edge_ec_to_degree_constraint_m.flatten().astype(np_impa_lib)
        edge_ec_to_degree_constraint_m_flatten_p = edge_ec_to_degree_constraint_m_flatten.ctypes.data_as(
            c_impa_lib_type_p
        )
        edge_degree_constraint_cost_flatten = self.edge_degree_constraint_cost.flatten().astype(np_impa_lib)
        edge_degree_constraint_cost_flatten_p = edge_degree_constraint_cost_flatten.ctypes.data_as(c_impa_lib_type_p)
        extrinsic_output_edge_ec = np.zeros((self.num_edge_variables))
        extrinsic_output_edge_ec = np.array(extrinsic_output_edge_ec).flatten().astype(np_impa_lib)
        extrinsic_output_edge_ec_p = extrinsic_output_edge_ec.ctypes.data_as(c_impa_lib_type_p)

        num_augmentations_p = np.array(self.num_augmentations).ctypes.data_as(c_int_p)
        num_added_constraints_p = np.array(self.num_added_constraints).ctypes.data_as(c_int_p)
        tour_impa_flatten = self.tour_impa.flatten().astype(np.int32)
        tour_impa_flatten_p = tour_impa_flatten.ctypes.data_as(c_int_p)
        cost_impa_p = np.array(self.cost_impa).ctypes.data_as(c_impa_lib_type_p)
        selected_edges_flatten = np.array(self.selected_edges).flatten().astype(np.int32)
        selected_edges_p = selected_edges_flatten.ctypes.data_as(c_int_p)
        selected_edges_size_p = np.array(0).ctypes.data_as(c_int_p)

        no_improvement_sol_count_exceeded_flag_p = np.array(self.no_improvement_sol_count_exceeded_flag).ctypes.data_as(
            c_bool_p
        )
        no_consecutive_loops_count_exceeded_flag_p = np.array(
            self.no_consecutive_loops_count_exceeded_flag
        ).ctypes.data_as(c_bool_p)
        sol_oscillation_count_exceeded_flag_p = np.array(self.sol_oscillation_count_exceeded_flag).ctypes.data_as(
            c_bool_p
        )

        subtour_paths_flatten = self.subtour_paths.flatten().astype(np.int32)
        subtour_paths_flatten_p = subtour_paths_flatten.ctypes.data_as(c_int_p)
        subtour_paths_size_flatten = np.zeros(self.num_nodes**2).flatten().astype(np.int32)
        subtour_paths_size_flatten_p = subtour_paths_size_flatten.ctypes.data_as(c_int_p)

        return (
            edge_connection_flatten_p,
            cost_edge_variable_flatten_p,
            cost_matrix_flatten_p,
            edge_ec_to_degree_constraint_m_flatten_p,
            edge_degree_constraint_cost_flatten_p,
            extrinsic_output_edge_ec_p,
            num_augmentations_p,
            tour_impa_flatten_p,
            cost_impa_p,
            selected_edges_p,
            selected_edges_size_p,
            no_improvement_sol_count_exceeded_flag_p,
            no_consecutive_loops_count_exceeded_flag_p,
            sol_oscillation_count_exceeded_flag_p,
            subtour_paths_flatten_p,
            subtour_paths_size_flatten_p,
            num_added_constraints_p,
        )

    def process_outputs_ctypes(
        self,
        sol_oscillation_count_exceeded_flag_p,
        selected_edges_size_p,
        selected_edges_p,
        subtour_paths_size_flatten_p,
        subtour_paths_flatten_p,
        no_improvement_sol_count_exceeded_flag_p,
        no_consecutive_loops_count_exceeded_flag_p,
        tour_impa_flatten_p,
        num_augmentations_p,
        cost_impa_p,
        extrinsic_output_edge_ec_p,
        num_added_constraints_p,
    ):
        selected_edges_size = list(selected_edges_size_p.__dict__.values())[0]
        selected_edges_flatten = list(selected_edges_p.__dict__.values())[0][:selected_edges_size]
        self.selected_edges = [
            [selected_edges_flatten[i], selected_edges_flatten[i + 1]] for i in range(0, len(selected_edges_flatten), 2)
        ]

        subtour_paths_size = list(subtour_paths_size_flatten_p.__dict__.values())[0]
        subtour_paths_size = [x for x in subtour_paths_size if x != 0]
        if len(subtour_paths_size) != 0:
            subtour_paths_flatten = list(subtour_paths_flatten_p.__dict__.values())[0][: sum(subtour_paths_size)]
            self.subtour_paths = [
                list(
                    islice(
                        subtour_paths_flatten,
                        sum(subtour_paths_size[:i]),
                        sum(subtour_paths_size[: i + 1]),
                    )
                )
                for i in range(len(subtour_paths_size))
            ]
        else:
            self.subtour_paths = None

        self.sol_oscillation_count_exceeded_flag = list(sol_oscillation_count_exceeded_flag_p.__dict__.values())[0]
        self.no_improvement_sol_count_exceeded_flag = list(no_improvement_sol_count_exceeded_flag_p.__dict__.values())[
            0
        ]
        self.no_consecutive_loops_count_exceeded_flag = list(
            no_consecutive_loops_count_exceeded_flag_p.__dict__.values()
        )[0]
        tour_impa = list(tour_impa_flatten_p.__dict__.values())[0]
        if all(element == 0 for element in tour_impa):
            self.tour_impa = None
        else:
            self.tour_impa = [x for x in tour_impa]

        self.num_augmentations = int(list(num_augmentations_p.__dict__.values())[0])
        if self.augmentation_flag:
            self.num_added_constraints = int(list(num_added_constraints_p.__dict__.values())[0])
        else:
            self.num_added_constraints = 0

        self.cost_impa = list(cost_impa_p.__dict__.values())[0][0]
        extrinsic_output_edge_ec = list(extrinsic_output_edge_ec_p.__dict__.values())[0]
        intrinsic_out_edge_ec = extrinsic_output_edge_ec + self.cost_edge_variable
        self.intrinsic_out_edge_ec = intrinsic_out_edge_ec
        if self.subtour_paths is not None:
            self.subtour_paths = [sublist + [sublist[0]] for sublist in self.subtour_paths if sublist]
        if self.tour_impa is None and self.post_process_flag:
            self.run_post_processing_improved()
            self.pp_performed = True
        else:
            self.pp_performed = False

    def undirected_graph(self, selected_edges):
        graph = {}
        for connection in selected_edges:
            if connection[0] in graph:
                graph[connection[0]].append(connection[1])
            else:
                graph[connection[0]] = [connection[1]]

            if connection[1] in graph:
                graph[connection[1]].append(connection[0])
            else:
                graph[connection[1]] = [connection[0]]
        return graph

    def paths_cleaning(self, removed_edges, investigated_paths):
        if investigated_paths is None:
            return None, None
        removed_edges = [tuple(edge) for edge in removed_edges]
        new_investigated_paths = []
        removed_edge_flag_list = []
        for path in investigated_paths:
            removed_edge_flag = False
            new_path = []
            i = 0
            while i < len(path) - 1:
                if (path[i], path[i + 1]) in removed_edges:
                    removed_edge_flag = True
                    if new_path:
                        new_path.append(path[i])
                        new_investigated_paths.append(new_path)
                    new_path = []
                    i += 1
                else:
                    new_path.append(path[i])
                    i += 1
            if i == len(path) - 1:  # Add the last node if there's one left
                new_path.append(path[-1])
            if new_path:  # Add the last path if it's not empty
                new_investigated_paths.append(new_path)
            removed_edge_flag_list.append(removed_edge_flag)
        final_investigated_paths = [
            path
            for i, path in enumerate(new_investigated_paths)
            if not any(set(path).issubset(p) for p in new_investigated_paths[:i] + new_investigated_paths[i + 1 :])
        ]
        return final_investigated_paths, removed_edge_flag_list

    def construct_graph(self, selected_edges):
        graph = {}
        for connection in selected_edges:
            if connection[0] in graph:
                graph[connection[0]][1].append(connection[1])
            else:
                graph[connection[0]] = ([], [connection[1]])

            if connection[1] in graph:
                graph[connection[1]][0].append(connection[0])
            else:
                graph[connection[1]] = ([connection[0]], [])
        return graph

    def run_post_processing_improved(self):
        print("--------")
        print("Post Processing Started:")
        print("--------")
        num_nodes = self.num_nodes
        # edge_connections = self.edge_connections
        # cost_edge_variable = self.cost_edge_variable
        subtour_paths = self.subtour_paths

        if subtour_paths is not None:
            for i, subtour_path in enumerate(subtour_paths):
                print(f"Subtour path of size {len(subtour_path)-1}: {subtour_path}")
        print("--------")
        self.selected_edges_pp = copy.deepcopy(self.selected_edges)
        selected_edges_pp = self.selected_edges_pp[:]
        degree_constraints_satisfied = False

        graph = self.construct_graph(selected_edges_pp)

        degree_constraints_satisfied = True
        keys = []
        for node, (inward_edges_nodes, outward_edges_nodes) in graph.items():
            if len(inward_edges_nodes) > 1 or len(outward_edges_nodes) > 1:
                keys.append(node)
            if len(inward_edges_nodes) != 1 or len(outward_edges_nodes) != 1:
                degree_constraints_satisfied = False

        print(f"Degree Constraints Satisfied? {degree_constraints_satisfied}")
        removed_edges_violated_dc = []
        if len(keys) != 0:
            print("Investigating Violated Degree Constraints with degree >1")
            for key in keys:
                inward_edges_nodes = graph[key][0]
                outward_edges_nodes = graph[key][1]
                print(
                    f"Node: {key}, Inward Edges Nodes: {inward_edges_nodes}, Outward Edges Nodes: {outward_edges_nodes}"
                )
                for i, connections in enumerate(graph[key]):
                    removed_edges_per_node = []  # ; indices_used = [];
                    while len(connections) - len(removed_edges_per_node) > 1:
                        indices_in_edge_connections = []
                        for node in connections:
                            if i == 0:
                                investigated_edge = [node, key]
                            else:
                                investigated_edge = [key, node]
                            if investigated_edge not in removed_edges_violated_dc:  # not in removed_edges_node):
                                index_edge_connections = next(
                                    (
                                        i
                                        for i, edge in enumerate(self.edge_connections)
                                        if np.array_equal(edge, investigated_edge)
                                    ),
                                    None,
                                )
                                indices_in_edge_connections.append(index_edge_connections)
                        # python3 main_wrapper_tsp.py --testFile=36 --alpha=0.5 --nITER=200 --inputPath=inputs_random_1000 --outputPath=outputs_random_1000_maxAugmCount50 --saveFlag=True --augmFlag=True --maxAugmCount=50 --filteringFlag=True
                        if (
                            indices_in_edge_connections
                        ):  # added since experienced a failure scenario where all outward edges where already removed
                            # flag_emptied_connection = False
                            # valid_indices = [
                            #    index for index in indices_in_edge_connections
                            # ]
                            min_index = indices_in_edge_connections[
                                np.argmin(np.abs(self.intrinsic_out_edge_ec[indices_in_edge_connections]))
                            ]
                            removed_edge = list(self.edge_connections[min_index])
                            removed_edge_cost = self.cost_edge_variable[min_index]
                            if removed_edge not in removed_edges_per_node:
                                removed_edges_per_node.append(removed_edge)
                            if removed_edge not in removed_edges_violated_dc:
                                print(
                                    f"Remove: {removed_edge} with cost: {removed_edge_cost}"
                                )  # or {self.cost_matrix[removed_edge[0]][removed_edge[1]]}
                                removed_edges_violated_dc.append(removed_edge)
                                indices_in_edge_connections.remove(min_index)
                            else:
                                print(f"Edge {removed_edge} with cost: {removed_edge_cost} already considered")
                        else:
                            print("index_edge_connections empty. All edges of violation were already removed")
                            break
        else:
            print("No Violated Degree Constraints with degree >1")

        subtour_paths, removed_edge_flag_list = self.paths_cleaning(removed_edges_violated_dc, self.subtour_paths)
        print(f"subtour_paths after removing removed_edges_violated_dc {removed_edges_violated_dc}: \n {subtour_paths}")
        print("------")

        if self.subtour_paths is not None:
            subtour_paths_disconnected_flag_list = removed_edge_flag_list[:]
            removed_subtour_paths_edges = []
            for i, is_true in enumerate(subtour_paths_disconnected_flag_list):
                if is_true:
                    print(f"subtour paths of index: {i}  was already disconnected")
                else:
                    print(f"Investigating connected subtour path of index: {i}")
                    subtour_path = self.subtour_paths[i]
                    pair_occurrences = self.find_indices_edges_from_tour(subtour_path)
                    min_index = pair_occurrences[np.argmin(np.abs(self.intrinsic_out_edge_ec[pair_occurrences]))]
                    removed_edge = list(self.edge_connections[min_index])
                    removed_edge_cost = self.cost_edge_variable[min_index]
                    print(f"Remove: {removed_edge} with cost: {removed_edge_cost}")
                    removed_subtour_paths_edges.append(removed_edge)

            subtour_paths, removed_edge_flag_list = self.paths_cleaning(removed_subtour_paths_edges, subtour_paths)
            print(
                f"subtour_paths after removing removed_subtour_paths_edges {removed_subtour_paths_edges}: \n {subtour_paths}"
            )

            selected_edges_pp = [edge for edge in selected_edges_pp if edge not in removed_subtour_paths_edges]

        selected_edges_pp = [edge for edge in selected_edges_pp if edge not in removed_edges_violated_dc]

        print("------")

        graph = self.construct_graph(selected_edges_pp)
        final_paths = self.get_paths(graph)

        degree_constraints_satisfied = True
        for node, (inward_edges_nodes, outward_edges_nodes) in graph.items():
            if len(inward_edges_nodes) != 1 or len(outward_edges_nodes) != 1:
                degree_constraints_satisfied = False

        print("Degree Constraints Satisfied?", degree_constraints_satisfied)

        for i, path in enumerate(final_paths):
            print(f"path {i} of size {len(path)}: {path}")

        combined_final_paths = list(chain.from_iterable(final_paths))
        missing_nodes = [num for num in range(num_nodes) if num not in combined_final_paths]
        print(f"missing_nodes: {missing_nodes}")
        unique_nodes, counts = np.unique(np.array(combined_final_paths), return_counts=True)
        duplicate_nodes_counts = list(zip(unique_nodes[counts > 1], counts[counts > 1]))
        if duplicate_nodes_counts:
            for node, count in duplicate_nodes_counts:
                print(f"Node {node} occurs {count} times.")

        if len(final_paths) > 1:
            print("------")
            print("Connecting Paths to form a tour")
            lst = range(0, len(final_paths))
            perms = permutations(lst)
            results = [perm for perm in perms if perm[0] == 0]

            possible_tours = []
            cost_tours = []

            for i, result in enumerate(results):
                print(f"Investigating connection {i} out of {len(results)} connections")
                tour = []
                count_used_paths = 0
                while count_used_paths < len(final_paths):
                    tour.extend(final_paths[result[count_used_paths]])
                    count_used_paths += 1

                if len(missing_nodes) != 0:
                    tour = self.add_missing_nodes(tour, missing_nodes)

                pair_occurrences = self.find_indices_edges_from_tour(tour)
                possible_tours.append(tour)
                cost_tours.append(sum(self.cost_edge_variable[pair_occurrences]))

            index_best_tour = np.argmin(cost_tours)
            tour_impa_pp = possible_tours[index_best_tour]
            cost_impa_pp = cost_tours[index_best_tour]
            min_cost = cost_impa_pp
            best_tour = tour_impa_pp[:]

            if self.k_opt_flag:
                print(f"tour_impa_pp: {tour_impa_pp+ [possible_tours[index_best_tour][0]]}")
                print(f"cost_impa_pp: {cost_impa_pp}")
                best_tour, min_cost = self.perform_k_opt(best_tour, min_cost)

            tour_impa_pp = best_tour[:]
            tour_impa_pp = tour_impa_pp + [tour_impa_pp[0]]
            cost_impa_pp = min_cost
            print(f"tour_impa_pp of size {len(tour_impa_pp)-1}: {tour_impa_pp}")
            print(f"cost_impa_pp: {cost_impa_pp}")

        else:
            print("------")
            print("Forming a tour from one path")
            if len(missing_nodes) != 0:
                print("Adding Missing Nodes")
                tour = self.add_missing_nodes(final_paths[0], missing_nodes)
            else:
                tour = final_paths[0]

            pair_occurrences = self.find_indices_edges_from_tour(tour)
            cost_impa_pp = sum(self.cost_edge_variable[pair_occurrences])
            tour_impa_pp = tour[:]

            if self.k_opt_flag:
                print(f"tour_impa_pp: {tour_impa_pp + [tour_impa_pp[0]]}")
                print(f"cost_impa_pp: {cost_impa_pp}")
                best_tour, min_cost = self.perform_k_opt(tour_impa_pp, cost_impa_pp)
            else:
                best_tour = tour_impa_pp[:]
                min_cost = cost_impa_pp

            tour_impa_pp = best_tour[:]
            tour_impa_pp = tour_impa_pp + [tour_impa_pp[0]]
            cost_impa_pp = min_cost
            print(f"tour_impa_pp of size {len(tour_impa_pp)-1}: {tour_impa_pp}")
            print(f"cost_impa_pp: {cost_impa_pp}")

        self.cost_impa = cost_impa_pp
        self.tour_impa = tour_impa_pp

    def perform_k_opt(self, best_tour, min_cost):
        print("------")
        print("Performing K-OPT")
        for k in [3, 2]:
            print(f"{k}-opt")
            candidate_paths = self.generate_k_opt_paths(best_tour, k)
            candidate_path = candidate_paths[0]
            for candidate_path in candidate_paths:
                pair_occurrences = self.find_indices_edges_from_tour(candidate_path)
                cost = sum(self.cost_edge_variable[element] for element in pair_occurrences)
                if cost < min_cost:
                    min_cost = cost
                    best_tour = candidate_path[:]
        return best_tour, min_cost

    def add_missing_nodes(self, tour, missing_nodes):
        added_missing_nodes = []
        while len(added_missing_nodes) != len(missing_nodes):
            missing_combinations_left = [[node, tour[0]] for node in missing_nodes if node not in added_missing_nodes]
            missing_combinations_right = [[tour[-1], node] for node in missing_nodes if node not in added_missing_nodes]
            min_cost_left = np_impa_lib("inf")
            min_cost_right = np_impa_lib("inf")
            for combination in missing_combinations_left:
                cost = self.cost_matrix[combination[0]][combination[1]]
                if cost < min_cost_left:
                    min_cost_left = cost
                    combination_left = combination
            for combination in missing_combinations_right:
                cost = self.cost_matrix[combination[0]][combination[1]]
                if cost < min_cost_right:
                    min_cost_right = cost
                    combination_right = combination
            costs_left_right = [min_cost_left, min_cost_right]
            # added_cost = np.min(costs_left_right)
            added_cost_index = np.argmin(costs_left_right)
            if added_cost_index == 0:  # left
                added_missing_nodes.append(combination_left[0])
                tour.insert(0, combination_left[0])
            else:
                added_missing_nodes.append(combination_right[1])
                tour.append(combination_right[1])
        return tour

    def find_indices_edges_from_tour(self, tour):
        pairs = np.column_stack((tour, np.roll(tour, -1)))
        edge_connections_arr = np.array(self.edge_connections)
        matching_pairs = (edge_connections_arr[:, None, :] == pairs).all(axis=-1)
        pair_occurrences = np.where(matching_pairs)[0]
        return pair_occurrences

    def find_paths(self, graph, node, visited, path=[]):
        visited.add(node)

        self.visited_nodes.append(node)

        path = path + [node]

        if len(graph[node][1]) == 0:
            self.open_paths_dummy.append(path)
            return path

        for neighbor in graph[node][1]:
            if neighbor not in visited:
                new_path = self.find_paths(graph, neighbor, visited.copy(), path.copy())
                if new_path:
                    return new_path

    def get_paths(self, graph):
        self.open_paths_dummy = []
        path_lists = []
        self.visited_nodes = []
        for start_node in graph.keys():
            if start_node in self.visited_nodes:
                continue
            else:
                path = self.find_paths(graph, start_node, set(), [])
                if path:
                    path_lists.append(path)
        final_paths = [
            path
            for i, path in enumerate(path_lists)
            if not any(set(path).issubset(p) for p in path_lists[:i] + path_lists[i + 1 :])
        ]
        return final_paths

    def generate_k_opt_paths(self, path, k):
        n = len(path)
        if k < 2 or k > n - 2:
            raise ValueError("k must be between 2 and n-2 for a valid k-opt operation.")

        swap_combinations = np.array(list(combinations(range(n), k)), dtype=np.int32)
        new_paths = []

        for swap_indices in swap_combinations:
            new_path = np.array(path, dtype=np.int32)
            swap_pairs = np.array(list(combinations(swap_indices, 2)), dtype=np.int32)
            new_path[swap_pairs[:, 0]], new_path[swap_pairs[:, 1]] = (
                new_path[swap_pairs[:, 1]],
                new_path[swap_pairs[:, 0]],
            )
            new_paths.append(new_path.tolist())

        return new_paths

    def run_pre_analysis(self):
        self.runtime_impa = time.time() - self.start_time
        print(f"IMPA Time: {self.runtime_impa}")
        print("--------")

        if self.exact_solver_flag:
            self.solve_tsp_exact()
        if self.lkh_solver_flag:
            start_time = time.time()
            self.tour_lkh, self.cost_lkh = self.solve_tsp_lkh()
            self.runtime_lkh = time.time() - start_time
            print(f"LKH Time: {self.runtime_lkh}")
        if self.sim_annealing_flag:
            self.solve_tsp_simulated_annealing()
        if self.save_flag:
            self.results_composed.append(
                (
                    self.reset_flag,
                    self.filtering_flag,
                    self.alpha,
                    self.cost_matrix,
                    self.tour_impa,
                    self.cost_impa,
                    self.tour_lkh,
                    self.cost_lkh,
                    self.pp_performed,
                    self.num_augmentations,
                    self.selected_edges,
                    self.sol_oscillation_count_exceeded_flag,
                    self.no_improvement_sol_count_exceeded_flag,
                    self.no_consecutive_loops_count_exceeded_flag,
                    self.intrinsic_out_edge_ec,
                    self.runtime_impa,
                    self.runtime_lkh,
                    self.num_added_constraints,
                )
            )

    def find_loops(self, graph, start, current, visited, path=[]):
        visited.add(current)

        self.visited_nodes_loop.append(current)

        path = path + [current]

        if current == start:
            visited.add(current)
            return path

        if current not in graph:
            return None

        for node in graph[current]:
            if node not in visited:
                new_path = self.find_loops(graph, start, node, visited.copy(), path)
                if new_path:
                    return new_path

    def get_loops(self, selected_edges, path_length):
        graph = {}
        for connection in selected_edges:
            if connection[0] in graph:
                graph[connection[0]].append(connection[1])
            else:
                graph[connection[0]] = [connection[1]]

        self.graph = graph
        self.visited_nodes_loop = []
        for start_node, end_node in graph.items():
            if start_node in self.visited_nodes_loop:
                continue
            else:
                for j in range(len(end_node)):
                    closed_loop = self.find_loops(graph, start_node, end_node[j], set())
                    if closed_loop is not None and len(closed_loop) == path_length:
                        return closed_loop

    def solve_tsp_exact(self):
        print("--------")
        print("Exact Solver: ")
        print("--------")
        tour_exact, cost_exact = solve_tsp_brute_force(self.cost_matrix)
        print("tour_exact: ", tour_exact)
        print("cost_exact: ", cost_exact)
        return tour_exact, cost_exact

    def solve_tsp_lkh(self):
        print("--------")
        print("LKH Solver: ")
        print("--------")
        cost_matrix_h = elkai.DistanceMatrix(self.cost_matrix)
        tour_h = cost_matrix_h.solve_tsp()
        cost_h = sum(self.cost_matrix[tour_h[i], tour_h[i + 1]] for i in range(len(tour_h) - 1))
        print("tour_h: ", tour_h)
        print("cost_h: ", cost_h)
        return tour_h, cost_h

    def solve_tsp_simulated_annealing(self):
        print("--------")
        print("Simulated Annealing Solver:")
        print("--------")
        tour_h, cost_h = solve_tsp_simulated_annealing(self.cost_matrix)
        print("tour_h: ", tour_h)
        print("cost_h: {:.4f}".format(cost_h))

    def save_outputs(self):
        with open(f"{self.folder_outputs}/outputs_set{self.test_file}.pkl", "wb") as f:
            pkl.dump(self.results_composed, f)


class ImpaKcMwm:
    def __init__(
        self,
        NUM_ITERATIONS,
        FILTERING_FLAG,
        POST_PROCESS_FLAG,
        ALPHA,
        THRESHOLD,
        PP_OPTION,
    ):
        self.num_iterations = NUM_ITERATIONS
        self.filtering_flag = FILTERING_FLAG
        self.post_process_flag = POST_PROCESS_FLAG
        self.alpha = ALPHA
        self.threshold = THRESHOLD
        self.post_process_option = PP_OPTION

    def initialize(self, input_load):
        max_state = np.array(input_load[0])
        teams_types = input_load[1]  # .tolist()
        num_projects = len(input_load[3][0])

        # max_state = np.array([3,3])
        # teams_types = np.array([2])
        # num_projects = 2

        self.teams_types = teams_types
        self.num_projects = num_projects

        print("N_u: ", max_state)
        print("teams_types: ", teams_types)
        print("number_of_projects: ", num_projects)

        print("Number of Iterations: ", self.num_iterations)

        if self.filtering_flag:
            print("Alpha: ", self.alpha)

        if self.post_process_flag:
            print("post_processing_Flag Enabled")
            if self.post_process_option == 1:  # brute force post-processing on departments
                print("Brute Force Post-Processing on Departments")
            elif self.post_process_option == 2:  # brute force post-processing on teams
                print("Brute Force Post-Processing on Teams")
        self.max_state = max_state
        self.num_departments = self.max_state.size

        # l_target = np.random.uniform(0, 1, num_projects)
        # k_target = np.random.randint(
        #    min(fighters_count), max(fighters_count) + 1, num_projects
        # )
        # s_target = np.random.randint(
        #    min(weasels_count), max(weasels_count) + 1, num_projects
        # )

        self.prune_teams()

        (
            reward_team,
            teams_weights_per_department,
            teams_types_per_department,
            last_index,
        ) = self.team_reward_generation()

        self.teams_weights_per_department = teams_weights_per_department
        self.last_index = last_index
        print("Total Number of Projects: ", num_projects)

        self.teams_types = teams_types_per_department[0]

        self.num_teams = len(reward_team)

        self.reward_team = np.array(input_load[2], dtype=np_impa_lib)
        self.reward_project = np.array(input_load[3], dtype=np_impa_lib).T
        # self.reward_team = -np.random.randint(0,100,size=self.num_teams)
        # reward_project = np.random.randint(0,100,size=(self.num_teams, self.num_projects))
        # self.reward_project = reward_project.T
        # print('self.reward_project: ', self.reward_project)
        # print('self.reward_team: ', self.reward_team)

        print("-------")

        self.intrinsic_out_mwm = np.zeros((self.num_projects, self.num_teams))

        if self.num_projects != self.num_teams:
            self.unbalanced_flag = True
        else:
            self.unbalanced_flag = False

        assert (
            self.intrinsic_out_mwm.shape[0] == self.reward_project.shape[0]
        ), "Row Shape mismatch between MWM and Rewards"
        assert (
            self.intrinsic_out_mwm.shape[1] == self.reward_project.shape[1]
        ), "Column Shape mismatch between MWM and Rewards"

    def prune_teams(self):
        units = range(1, len(self.max_state) + 1)

        permutations_units = [p for p in itertools.product(units, repeat=2)]
        combinations_units_pre_pruning = list(combinations(units, 2))

        list_pruning = []

        for i in range(0, len(combinations_units_pre_pruning)):
            combination = combinations_units_pre_pruning[i]
            r_1 = math.exp(-F_u[combination[0] - 1]) * math.exp(-W_u[combination[1] - 1])
            r_2 = math.exp(-F_u[combination[1] - 1]) * math.exp(-W_u[combination[0] - 1])

            if (
                r_1 > r_2
                and distance_metric[combination[1] - 1][combination[0] - 1]
                == distance_metric[combination[0] - 1][combination[1] - 1]
                and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1]
                == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]
            ):
                list_pruning.append(combination[::-1])
            elif (
                r_2 > r_1
                and distance_metric[combination[1] - 1][combination[0] - 1]
                == distance_metric[combination[0] - 1][combination[1] - 1]
                and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1]
                == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]
            ):
                list_pruning.append(combination)
            elif (
                r_2 == r_1
                and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1]
                == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]
            ):
                list_pruning.append(combination[::-1])

        self.available_combinations = [x for x in permutations_units if x not in list_pruning]

    def team_reward_generation(self):
        teams_weights_per_department = []
        teams_types_per_department = []
        reward_team = []
        total_team_array = []

        for state_index in range(0, len(self.max_state)):
            teams_weights_per_department.append([])
            teams_types_per_department.append([])

        sum_teams = []

        max_state = self.max_state

        indices_departments = np.array(range(1, len(max_state) + 1))

        for combination in self.available_combinations:
            u = combination[0]
            v = combination[1]
            team_array = []
            for type in self.teams_types:
                if type == 1:
                    if u == v:
                        team_size = max_state[u - 1]
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1
                    else:
                        team_size = 0
                elif type == 2:
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 2)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = min(max_state[u - 1], max_state[v - 1])
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend(
                            [weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        teams_types_per_department[v - 1].append(type)
                        teams_types_per_department[v - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                elif type == 3:
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 2)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1
                    else:
                        team_size = 0
                elif type == 4:
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 3)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = math.floor(min(max_state[u - 1] / 2, max_state[v - 1]))
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend(
                            [weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        teams_types_per_department[v - 1].append(type)
                        teams_types_per_department[v - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                elif type == 5:
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 4)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = math.floor(min(max_state[u - 1] / 2, max_state[v - 1] / 2))
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend(
                            [fighters_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend(
                            [weasels_count[type - 1] for i in range(0, team_size - 1)]
                        )
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        teams_types_per_department[v - 1].append(type)
                        teams_types_per_department[v - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]

                if team_size != 0:
                    reward_team.append(reward)
                    reward_team.extend([reward for i in range(0, team_size - 1)])
                    for departmemt_index in remaining_departments:
                        teams_weights_per_department[departmemt_index - 1].append(0)
                        teams_weights_per_department[departmemt_index - 1].extend([0 for i in range(0, team_size - 1)])
                        teams_types_per_department[departmemt_index - 1].append(type)
                        teams_types_per_department[departmemt_index - 1].extend([type for i in range(0, team_size - 1)])
                team_array.append(team_size)

            sum_teams.append(sum(team_array))
            total_team_array.append(team_array)
            # print(' u: ', u-1, '  v: ', v-1, '  Packages Size:  ', team_array, 'with sum: ', sum(team_array), '& last index: ', np.sum(sum_teams)-1)

        print("Total Number of Packages: ", np.sum(sum_teams))

        sum_teams = np.array(sum_teams)
        last_index = []
        for i in range(0, len(sum_teams)):
            if i == 0:
                last_index.append(sum_teams[i] - 1)
            else:
                last_index.append(sum(sum_teams[0:i]) + sum_teams[i] - 1)

        return (
            reward_team,
            teams_weights_per_department,
            teams_types_per_department,
            last_index,
        )

    def process_inputs_ctypes(self):
        teams_weights_per_department = self.teams_weights_per_department
        reward_project = self.reward_project
        reward_team = self.reward_team
        max_state = self.max_state
        num_departments = self.num_departments
        num_teams = self.num_teams
        num_projects = self.num_projects

        message_team_to_department = np.zeros((num_departments, num_teams), dtype=np_impa_lib)

        intrinsic_out_mwm = np.zeros((num_projects, num_teams), dtype=np_impa_lib)

        for i in range(0, num_departments):
            for j in range(0, num_teams):
                if self.teams_weights_per_department[i][j] != 0:
                    message_team_to_department[i][j] = reward_team[j]

        non_zero_weight_indices = []

        for department_index in range(0, num_departments):
            indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
            non_zero_weight_indices = non_zero_weight_indices + [indices]

        extrinsic_output_team = np.zeros((num_teams))
        transitional_model = message_team_to_department.flatten().astype(np_impa_lib)
        transitional_model_p = transitional_model.ctypes.data_as(c_impa_lib_type_p)
        max_state = max_state.flatten().astype(np.int32)
        max_state_p = max_state.ctypes.data_as(c_int_p)
        reward_team = reward_team.flatten().astype(np_impa_lib)
        reward_team_p = reward_team.ctypes.data_as(c_impa_lib_type_p)
        reward_project_flatten = reward_project.flatten().astype(np_impa_lib)
        reward_project_flatten_p = reward_project_flatten.ctypes.data_as(c_impa_lib_type_p)
        teams_weights_per_department_flatten = np.array(teams_weights_per_department).flatten().astype(np.int32)
        teams_weights_per_department_flatten_p = teams_weights_per_department_flatten.ctypes.data_as(c_int_p)
        extrinsic_output_team = np.array(extrinsic_output_team).flatten().astype(np_impa_lib)
        extrinsic_output_team_p = extrinsic_output_team.ctypes.data_as(c_impa_lib_type_p)
        intrinsic_out_mwm = np.array(intrinsic_out_mwm).flatten().astype(np_impa_lib)
        intrisic_out_mwm_p = intrinsic_out_mwm.ctypes.data_as(c_impa_lib_type_p)

        non_zero_weight_indices_sizes = [len(indices) for indices in non_zero_weight_indices]

        max_size_nonzero_weights = max(non_zero_weight_indices_sizes)
        non_zero_weight_indices_arr = np.zeros(
            (num_departments, max_size_nonzero_weights), dtype=np.int32
        )  # added dtype=np.int32

        for i in range(len(non_zero_weight_indices)):
            for j in range(len(non_zero_weight_indices[i])):
                non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

        non_zero_weight_indices = np.array(non_zero_weight_indices_arr).flatten().astype(np.int32)
        non_zero_weight_indices_p = non_zero_weight_indices.ctypes.data_as(c_int_p)
        non_zero_weight_indices_sizes = np.array(non_zero_weight_indices_sizes).flatten().astype(np.int32)
        non_zero_weight_indices_sizes_p = non_zero_weight_indices_sizes.ctypes.data_as(c_int_p)

        return (
            transitional_model_p,
            max_state_p,
            reward_team_p,
            reward_project_flatten_p,
            teams_weights_per_department_flatten_p,
            extrinsic_output_team_p,
            intrisic_out_mwm_p,
            non_zero_weight_indices_p,
            non_zero_weight_indices_sizes_p,
            max_size_nonzero_weights,
        )

    def iterate(self):
        (
            transitional_model_p,
            max_state_p,
            reward_team_p,
            reward_project_flatten_p,
            teams_weights_per_department_flatten_p,
            extrinsic_output_team_p,
            intrisic_out_mwm_p,
            non_zero_weight_indices_p,
            non_zero_weight_indices_sizes_p,
            max_size_nonzero_weights,
        ) = self.process_inputs_ctypes()

        WrapperKcMwm(
            np.int32(self.num_iterations),
            np.int32(self.num_departments),
            np.int32(self.num_teams),
            np.int32(self.num_projects),
            transitional_model_p,
            max_state_p,
            reward_team_p,
            reward_project_flatten_p,
            teams_weights_per_department_flatten_p,
            non_zero_weight_indices_p,
            non_zero_weight_indices_sizes_p,
            np.int32(max_size_nonzero_weights),
            extrinsic_output_team_p,
            intrisic_out_mwm_p,
            self.filtering_flag,
            np_impa_lib(self.alpha),
        )

        extrinsic_output_team = list(extrinsic_output_team_p.__dict__.values())[0]
        intrinsic_out_mwm = np.reshape(
            list(intrisic_out_mwm_p.__dict__.values())[0],
            (self.num_projects, self.num_teams),
        )
        self.intrinsic_out_mwm = intrinsic_out_mwm

        self.intrinsic_output = extrinsic_output_team + self.reward_team

    def pre_analysis(self):
        intrinsic_output = self.intrinsic_output
        threshold = self.threshold

        valid_match, impa_metric, p_mwm = self.check_match()
        self.matching_array = np.where(p_mwm == 1)
        selected_teams_indices = np.where(p_mwm == 1)[1]
        self.selected_teams_indices = np.unique(selected_teams_indices)
        hard_decision = np.array(deepcopy(intrinsic_output))
        hard_decision[intrinsic_output > threshold] = 0

        self.pre_processing_analysis()

    def check_match(self):
        A = self.intrinsic_out_mwm

        unbalanced_flag = self.unbalanced_flag

        match_metric = sum(x for x in A.flatten() if x < 0)
        P = deepcopy(A)
        threshold = self.threshold
        P[P > threshold] = 0
        P[P < threshold] = 1
        valid_match = False

        min_N = min(A.shape[0], A.shape[1])

        if unbalanced_flag:
            P_interm = self.check_unbalanced_flag(P)
        else:
            P_interm = P
        row_sums = np.sum(P_interm, axis=1)
        col_sums = np.sum(P_interm, axis=0)

        if np.array_equal(col_sums, np.ones(min_N)):
            if np.array_equal(row_sums, np.ones(min_N)):
                valid_match = True
        return valid_match, match_metric, P

    def check_unbalanced_flag(self, P):
        idx = 0
        idy = 0
        P_interm = deepcopy(P)
        idy = np.argwhere(np.all(P_interm[..., :] == 0, axis=0))  # find which columns there are zeros
        P_interm = np.delete(P_interm, idy, axis=1)
        idx = np.argwhere(np.all(P_interm[:, ...] == 0, axis=1))  # find which rows there are zeross
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

        all_projects = np.array(range(0, num_projects))

        team_departments_weights = np.zeros((num_departments), dtype=int)

        department_weights = []

        for j in range(0, num_departments):
            department_weights.append([])

        print("Matching Array: ", matching_array)

        print("-------------")
        print("Air Units Assignment")
        print("-------------")

        team_type = []
        team_locations = []
        total_value = 0
        total_weight = 0
        for i in range(0, len(selected_teams_indices)):
            team_departments = []
            index = selected_teams_indices[i]
            project_index = np.where(matching_array[1] == index)[0]
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            for j in range(0, num_departments):
                weight = teams_weights_per_department[j][index]
                if weight != 0:
                    team_departments.append(j)
                    department_weights[j].append(weight)
                    team_departments_weights[j] = team_departments_weights[j] + weight
            location_package = next(x[0] for x in enumerate(last_index) if x[1] >= index)
            team_locations.append(location_package)
            print(
                "Package ",
                index,
                " has type ",
                teams_types[index],
                " & is assigned to target ",
                matching_array[0][project_index],
                " & is in units ",
                np.array(available_combinations[location_package]) - 1,
            )

        print("-------")
        print("BEFORE POST-PROCESSING")
        for i in range(0, num_departments):
            print("Unit ", i, " has a used capacity ", team_departments_weights[i])
            print("Unit ", i, " has a MAX capacity ", max_state[i])
            if team_departments_weights[i] > max_state[i]:
                print("Unit ", i, " needs post-processing")
                capacity_exceeded_flag = True
        print("Exceeded Capacity Flag: ", capacity_exceeded_flag)
        print("-------")
        print("Total Value LHS: ", -total_value)
        print(
            "Total number of selected teams: ",
            len(selected_teams_indices),
            "out of ",
            len(reward_team),
        )
        print(
            "Total number of Matched Projects: ",
            len(matching_array[0]),
            "out of ",
            num_projects,
        )
        print(
            "Total number of Unmatched Projects: ",
            len(np.setdiff1d(all_projects, matching_array[0])),
        )

        print("-------")
        print(
            "Types of Teams: ",
            np.unique(team_type, return_counts=True)[0],
            "Counts: ",
            np.unique(team_type, return_counts=True)[1],
        )
        print("-------------")

        total_weight = 0
        visited_teams = []
        reused_teams = []
        combination_reused_team = []
        visited_combinations = []
        for i in range(0, len(matching_array[0])):
            team_index = matching_array[1][i]
            project_index = matching_array[0][i]
            if team_index not in visited_teams:
                visited_teams.append(team_index)
                visited_combinations.append((project_index, team_index))
            else:
                reused_teams.append(team_index)
                index_used = visited_teams.index(team_index)
                used_combination = visited_combinations[index_used]
                combination_reused_team.append([(project_index, team_index), used_combination])
                ic_violated_flag = True
            total_weight = total_weight + reward_project[project_index, team_index]

        print("MWM")
        print("-------------")

        if len(matching_array[0]) == num_projects:
            project_matched_flag = True

        print(f"Matching All Projects to Teams : {project_matched_flag}")
        print(
            "Team assigned to multiple Projects: ",
            ic_violated_flag,
            "with index: ",
            np.unique(reused_teams),
        )
        print("Total Weight RHS: ", total_weight)

        self.capacity_exceeded_flag = capacity_exceeded_flag
        self.ic_violated_flag = ic_violated_flag
        self.combination_reused_team = combination_reused_team
        self.team_departments_weights = team_departments_weights
        self.team_locations = team_locations
        self.reused_teams = reused_teams

    def post_analysis(self):
        post_process_flag = self.post_process_flag

        if post_process_flag:
            if self.post_process_option == 1:  # brute force post-processing on bins
                self.apply_post_processing_departments_brute()
            elif self.post_process_option == 2:  # brute force post-processing on packages
                self.apply_post_processing_teams_brute()
                # exit()
        else:
            self.intermediate_selected_teams_indices = self.selected_teams_indices
            self.post_processing_ic_combination = []

        self.results_analysis()

    def apply_post_processing_departments_brute(self):
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

        print("-------")
        post_processing_ic_combination = []
        if ic_violated_flag:
            post_processing_ic_combination = [0] * len(combination_reused_team)
            for i in range(0, len(combination_reused_team)):
                min_intrinsic = 1000000000
                for j in range(0, len(combination_reused_team[i])):
                    combination = combination_reused_team[i][j]
                    if intrinsic_out_mwm[combination[0], combination[1]] < min_intrinsic:
                        min_intrinsic = intrinsic_out_mwm[combination[0], combination[1]]
                        post_processing_ic_combination[i] = combination
                    print(
                        "Cost of ",
                        combination,
                        "is: ",
                        reward_project[combination[0], combination[1]],
                    )
                    print(
                        "Intrinsic Output RHS",
                        intrinsic_out_mwm[combination[0], combination[1]],
                    )
            print("post_processing_IC_combination: ", post_processing_ic_combination)

        intermediate_team_department_weights = deepcopy(team_departments_weights)
        intermediate_selected_teams_indices = deepcopy(selected_teams_indices)

        if capacity_exceeded_flag:
            print("-------")
            print("POST-PROCESSING STARTED")
            print("-------")
            while any(intermediate_team_department_weights > max_state):
                post_processing_departments = np.where(intermediate_team_department_weights > max_state)[0]
                print("Post-processing Departments: ", post_processing_departments)
                print("-------")
                if len(post_processing_departments):
                    intrinsic_output_selected = intrinsic_output[intermediate_selected_teams_indices]
                    processing_department_index = post_processing_departments[0]
                    print("Post-processing department: ", processing_department_index)
                    packages_processing_department_index = [
                        i
                        for i in intermediate_selected_teams_indices
                        if (
                            processing_department_index
                            in np.array(
                                available_combinations[
                                    team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                ]
                            )
                            - 1
                        )
                    ]
                    print(
                        "max(reward_team[packages_processing_department_index]: ",
                        max(reward_team[packages_processing_department_index]),
                    )
                    removed_team = intermediate_selected_teams_indices[
                        np.where(
                            intrinsic_output_selected
                            == max(
                                intrinsic_output_selected[
                                    intermediate_selected_teams_indices.searchsorted(
                                        packages_processing_department_index
                                    )
                                ]
                            )
                        )[0][0]
                    ]
                    print("reward[removed_team]: ", reward_team[removed_team])
                    removed_team_location = team_locations[
                        np.where(intermediate_selected_teams_indices == removed_team)[0][0]
                    ]
                    removed_team_combination = np.array(available_combinations[removed_team_location]) - 1
                    removed_team_type = teams_types[removed_team]
                    print(
                        "Removing Team ",
                        removed_team,
                        "of type ",
                        removed_team_type,
                        " from bins ",
                        removed_team_combination,
                    )
                    intermediate_selected_teams_indices = np.setdiff1d(
                        intermediate_selected_teams_indices, removed_team
                    )
                    team_locations.remove(removed_team_location)
                    u = removed_team_combination[0]
                    v = removed_team_combination[1]
                    if removed_team_type == 1:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                    elif removed_team_type == 2:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif removed_team_type == 3:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                    elif removed_team_type == 4:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif removed_team_type == 5:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 2

            print("-------")
            print("POST-PROCESSING ENDED")
            print("-------")

        self.intermediate_selected_teams_indices = intermediate_selected_teams_indices
        self.post_processing_ic_combination = post_processing_ic_combination

    def apply_post_processing_teams_brute(self):
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
        # reward_team = self.reward_team
        teams_types = self.teams_types
        capacity_exceeded_flag = self.capacity_exceeded_flag

        print("-------")
        post_processing_ic_combination = []
        if ic_violated_flag:
            post_processing_ic_combination = [0] * len(combination_reused_team)
            for i in range(0, len(combination_reused_team)):
                min_intrinsic = 1000000000
                for j in range(0, len(combination_reused_team[i])):
                    combination = combination_reused_team[i][j]
                    if intrinsic_out_mwm[combination[0], combination[1]] < min_intrinsic:
                        min_intrinsic = intrinsic_out_mwm[combination[0], combination[1]]
                        post_processing_ic_combination[i] = combination
                    print(
                        "Cost of ",
                        combination,
                        "is: ",
                        reward_project[combination[0], combination[1]],
                    )
                    print(
                        "Intrinsic Output RHS",
                        intrinsic_out_mwm[combination[0], combination[1]],
                    )
            print("post_processing_IC_combination: ", post_processing_ic_combination)

        intermediate_team_department_weights = deepcopy(team_departments_weights)
        intermediate_selected_teams_indices = deepcopy(selected_teams_indices)

        if capacity_exceeded_flag:
            print("-------")
            print("POST-PROCESSING STARTED")
            print("-------")
            while any(intermediate_team_department_weights > max_state):
                post_processing_departments = np.where(intermediate_team_department_weights > max_state)[0]
                print("Post-processing Departments: ", post_processing_departments)
                print("-------")
                list_teams_processing_department_indices_combinations = []
                if len(post_processing_departments):
                    intrinsic_output_selected = intrinsic_output[intermediate_selected_teams_indices]
                    for department in post_processing_departments:
                        teams_processing_department_indices_combinations = [
                            [
                                i,
                                (
                                    np.array(
                                        available_combinations[
                                            team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                        ]
                                    )
                                    - 1
                                ).tolist(),
                            ]
                            for i in intermediate_selected_teams_indices
                            if (
                                (
                                    department
                                    in np.array(
                                        available_combinations[
                                            team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                        ]
                                    )
                                    - 1
                                )
                                and available_combinations[
                                    team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                ][0]
                                - 1
                                in post_processing_departments
                                and available_combinations[
                                    team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                ][1]
                                - 1
                                in post_processing_departments
                            )
                        ]  # and
                        # len(set(np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices==i)[0][0]]])-1)) != 1)]
                        if not len(teams_processing_department_indices_combinations):
                            teams_processing_department_indices_combinations = [
                                [
                                    i,
                                    (
                                        np.array(
                                            available_combinations[
                                                team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                            ]
                                        )
                                        - 1
                                    ).tolist(),
                                ]
                                for i in intermediate_selected_teams_indices
                                if (
                                    department
                                    in np.array(
                                        available_combinations[
                                            team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]
                                        ]
                                    )
                                    - 1
                                )
                            ]
                        # print('department: ',department, '\n teams_processing_department_indices_combinations: \n', teams_processing_department_indices_combinations)
                        # print(len(teams_processing_department_indices_combinations))
                        list_teams_processing_department_indices_combinations.extend(
                            y
                            for y in teams_processing_department_indices_combinations
                            if y not in list_teams_processing_department_indices_combinations
                        )
                    print(
                        "\nlist_packages_processing_bin_indices_combinations:\n ",
                        list_teams_processing_department_indices_combinations,
                    )
                    team_indices_processing = [i[0] for i in list_teams_processing_department_indices_combinations]
                    print("team_indices_processing: ", team_indices_processing)
                    team_combinations_processing = [i[1] for i in list_teams_processing_department_indices_combinations]
                    print("team_combinations_processing: ", team_combinations_processing)
                    flag_combinations = [0 if i[0] == i[1] else 1 for i in team_combinations_processing]
                    print("flag_combinations: ", flag_combinations)
                    if np.all(np.array(flag_combinations) == 0):
                        index_processing_combinations = np.where(np.array(flag_combinations, dtype=np.int) == 0)[0]
                        # print('index_processing_combinations: ', index_processing_combinations)
                    else:
                        index_processing_combinations = np.where(np.array(flag_combinations, dtype=np.int) != 0)[0]
                    print(
                        "Set of teams to remove from: ",
                        np.array(team_indices_processing)[index_processing_combinations],
                    )
                    print(
                        "intrinsic_output_selected :",
                        intrinsic_output_selected[
                            intermediate_selected_teams_indices.searchsorted(
                                np.array(team_indices_processing)[index_processing_combinations]
                            )
                        ],
                    )
                    removed_team = intermediate_selected_teams_indices[
                        np.where(
                            intrinsic_output_selected
                            == max(
                                intrinsic_output_selected[
                                    intermediate_selected_teams_indices.searchsorted(
                                        np.array(team_indices_processing)[index_processing_combinations]
                                    )
                                ]
                            )
                        )[0][0]
                    ]
                    print("removed_team:", removed_team)
                    # print('reward[removed_team]: ', reward_team[removed_team])
                    removed_team_location = team_locations[
                        np.where(intermediate_selected_teams_indices == removed_team)[0][0]
                    ]
                    removed_team_combination = np.array(available_combinations[removed_team_location]) - 1
                    removed_team_type = teams_types[removed_team]
                    print(
                        "Removing Package ",
                        removed_team,
                        "of type ",
                        removed_team_type,
                        " from bins ",
                        removed_team_combination,
                    )
                    intermediate_selected_teams_indices = np.setdiff1d(
                        intermediate_selected_teams_indices, removed_team
                    )
                    team_locations.remove(removed_team_location)
                    u = removed_team_combination[0]
                    v = removed_team_combination[1]
                    if removed_team_type == 1:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                    elif removed_team_type == 2:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 1
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif removed_team_type == 3:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                    elif removed_team_type == 4:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 1
                    elif removed_team_type == 5:
                        intermediate_team_department_weights[u] = intermediate_team_department_weights[u] - 2
                        intermediate_team_department_weights[v] = intermediate_team_department_weights[v] - 2
            print("-------")
            print("POST-PROCESSING ENDED")
            print("-------")

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

        post_processing_team_departments_weights = np.zeros((num_departments), dtype=int)
        team_type = []
        total_value = 0
        project_array = []
        # test=[]
        for i in range(0, len(intermediate_selected_teams_indices)):
            index = intermediate_selected_teams_indices[i]
            if index in reused_teams and post_process_flag:
                project_index = np.where(
                    matching_array[0] == [i[0] for i in post_processing_ic_combination if i[1] == index][0]
                )[0]
            else:
                project_index = np.where(matching_array[1] == index)[0]
            project_array.append(project_index)
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            for j in range(0, num_departments):
                weight = teams_weights_per_department[j][index]
                if weight != 0:
                    post_processing_team_departments_weights[j] = post_processing_team_departments_weights[j] + weight
            location_package = team_locations[i]
            print(
                "Team ",
                index,
                " has type ",
                teams_types[index],
                " & is assigned to project ",
                matching_array[0][project_index],
                " & is in units ",
                np.array(available_combinations[location_package]) - 1,
            )
            # print('reward_team[index]: ', reward_team[index])
            # test.append(reward_team[index])
            results_composed.append(
                (
                    available_combinations[location_package][0],
                    available_combinations[location_package][1],
                    teams_types[index],
                    (matching_array[1][project_index] + 1).tolist(),
                    (matching_array[0][project_index] + 1).tolist(),
                )
            )
        # print(test)
        print("-------")
        for i in range(0, num_departments):
            print(
                "Unit ",
                i,
                " has a used capacity ",
                post_processing_team_departments_weights[i],
            )
            print("Unit ", i, " has a MAX capacity ", max_state[i])
            if post_processing_team_departments_weights[i] > max_state[i]:
                print("Unit ", i, " needs post-processing")
                post_processing_capacity_exceeded_flag = True
        print(
            "Post-Processing Exceeded Capacity Flag: ",
            post_processing_capacity_exceeded_flag,
        )

        # Added to compare OR-tools & IMPA
        self.results_composed = results_composed

        print("-------")
        print(
            "Types of Teams: ",
            np.unique(team_type, return_counts=True)[0],
            "Counts: ",
            np.unique(team_type, return_counts=True)[1],
        )
        print("-------------")

        print("-------")
        print("Total Value LHS: ", -total_value)
        print(
            "Total number of selected teams: ",
            len(intermediate_selected_teams_indices),
            "out of ",
            len(reward_team),
        )
        print(
            "Total number of Matched Projects: ",
            len(project_array),
            "out of ",
            num_projects,
        )

        total_weight = 0
        total_weight_used_teams = []
        for i in range(0, len(matching_array[1])):
            team_index = matching_array[1][i]
            if team_index in reused_teams and post_process_flag:
                project_index = [i[0] for i in post_processing_ic_combination if i[1] == team_index][0]
            else:
                project_index = matching_array[0][i]
            if team_index in intermediate_selected_teams_indices and team_index not in total_weight_used_teams:
                total_weight = total_weight + reward_project[project_index, team_index]
                total_weight_used_teams.append(team_index)

        print("MWM")
        print("-------------")

        if len(project_array) == num_projects:
            project_matched_flag = True

        print(f"Matching All Projects to Teams : {project_matched_flag}")
        print("Total Weight RHS: ", total_weight)
