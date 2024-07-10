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
        NUM_NODES,
        SYMMETRIC_FLAG,
        AUGMENTATION_FLAG,
        EXACT_SOLVER_FLAG,
        LKH_SOLVER_FLAG,
        SIM_ANNEALING_FLAG,
        NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE,
        RESET_FLAG,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
    ):
        self.num_iterations = NUM_ITERATIONS
        self.num_nodes = NUM_NODES
        self.threshold = THRESHOLD
        self.symmetric_flag = SYMMETRIC_FLAG
        self.augmentation_flag = AUGMENTATION_FLAG
        self.exact_solver_flag = EXACT_SOLVER_FLAG
        self.lkh_solver_flag = LKH_SOLVER_FLAG
        self.sim_annealing_flag = SIM_ANNEALING_FLAG
        self.num_augmented_samples_per_subtour_size = NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE
        self.reset_flag = RESET_FLAG
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.random_test_flag = RANDOM_TEST_FLAG

    def initialize(
        self,
    ):
        self.results_composed = []

        if not self.random_test_flag:
            self.num_nodes = self.input_load[0]
            self.symmetric_flag = self.input_load[1]

        num_nodes = self.num_nodes
        symmetric_flag = self.symmetric_flag

        self.subtour_constraints_satisfied_flag = False

        (
            cost_matrix,
            edge_connections,
        ) = self.create_cost_matrix(symmetric_flag)
        self.tour_impa_flag = False

        if not self.random_test_flag:
            self.cost_matrix = self.input_load[2]
            self.edge_connections = self.input_load[3]
            cost_matrix = self.cost_matrix
            edge_connections = self.edge_connections

        self.cost_matrix = cost_matrix
        self.edge_connections = edge_connections

        num_edge_variables = num_nodes**2 - num_nodes

        self.num_edge_variables = num_edge_variables

        self.intrinsic_out_edge_ec = np.zeros(
            num_edge_variables,
            dtype=np_impa_lib,
        )

        cost_edge_variable = np.zeros(
            num_edge_variables,
            dtype=np_impa_lib,
        )

        edge_ec_to_degree_constraint_m = np.zeros(
            (
                num_edge_variables,
                num_nodes,
            ),
            dtype=np_impa_lib,
        )

        for i in range(len(edge_connections)):
            connection = edge_connections[i]
            cost = cost_matrix[
                connection[0],
                connection[1],
            ]
            cost_edge_variable[i] = cost
            edge_ec_to_degree_constraint_m[i][connection[0]] = cost
            edge_ec_to_degree_constraint_m[i][connection[1]] = cost

        self.cost_edge_variable = cost_edge_variable
        self.edge_degree_constraint_cost = edge_ec_to_degree_constraint_m

        self.edge_ec_to_degree_constraint_m = edge_ec_to_degree_constraint_m
        self.model_eq_constraint = EqualityConstraint(
            num_nodes,
            num_edge_variables,
            cost_matrix,
            self.edge_connections,
            self.edge_degree_constraint_cost,
            self.augmentation_flag,
            self.filtering_flag,
            self.alpha,
        )
        self.model_degree_constraint = DegreeConstraint(
            num_nodes,
            num_edge_variables,
            edge_connections,
            self.filtering_flag,
            self.alpha,
        )
        self.outputs = OutputsImpa(
            num_nodes,
            num_edge_variables,
            cost_matrix,
            self.augmentation_flag,
        )
        self.delta_S_indices_list = []
        if self.augmentation_flag:
            self.model_subtour_constraint = SubtourEliminationConstraint(
                num_nodes,
                num_edge_variables,
                self.filtering_flag,
                self.alpha,
            )

    def create_cost_matrix(
        self,
        symmetric_flag,
    ):
        n = self.num_nodes
        upper_triangle = np.random.uniform(
            1,
            1000,
            size=(
                n,
                n,
            ),
        )
        if symmetric_flag:
            matrix = (
                np.triu(upper_triangle)
                + np.triu(
                    upper_triangle,
                    k=1,
                ).T
            )
            np.fill_diagonal(
                matrix,
                zero_value,
            )
        else:
            matrix = upper_triangle.copy()
            np.fill_diagonal(
                matrix,
                zero_value,
            )
        edge_connections = np.transpose(np.nonzero(matrix))
        return (
            matrix,
            edge_connections,
        )

    def iterate_relaxed_graph(
        self,
    ):
        for iter in range(
            0,
            self.num_iterations,
        ):
            self.model_degree_constraint.degree_constraint_to_edge_ec_update(self.edge_ec_to_degree_constraint_m)

            self.model_degree_constraint.process_filtering(iter)

            self.edge_ec_to_degree_constraint_m = self.model_eq_constraint.edge_ec_to_degree_constraint_relaxed_graph_update(self.model_degree_constraint.degree_constraint_to_eq_constraint_m)

        extrinsic_output_edge_ec = self.outputs.extrinsic_output_edge_ec_relaxed_graph_update(self.model_degree_constraint.degree_constraint_to_eq_constraint_m)
        intrinsic_out_edge_ec = extrinsic_output_edge_ec + self.cost_edge_variable
        self.intrinsic_out_edge_ec = intrinsic_out_edge_ec

    def pre_analysis(
        self,
    ):
        self.hard_decision_analysis()

        self.subtour_elimination_constraints_analysis()

        if (not self.subtour_constraints_satisfied_flag) and self.augmentation_flag:
            self.subtour_constraints_to_edge_ec_m = [[zero_value] * self.num_edge_variables for _ in range(len(self.delta_S_indices_list))]

        while (not self.subtour_constraints_satisfied_flag) and self.augmentation_flag and not self.tour_impa_flag:
            self.selected_edges_old = deepcopy(self.selected_edges)
            self.iterate_augmented_graph()
            self.hard_decision_analysis()
            if self.selected_edges.shape == self.selected_edges_old.shape and {frozenset(sub) for sub in self.selected_edges} == {frozenset(sub) for sub in self.selected_edges_old}:
                print("Exited: No improvement in IMPA solution")
                if np.any(np.isnan(self.edge_ec_to_degree_constraint_m)):
                    print("Exited: No improvement in IMPA solution")
                    print("NAN detected")
                    exit("Exited: No improvement in IMPA solution")

            self.subtour_elimination_constraints_analysis()
            if len(self.subtour_constraints_to_edge_ec_m) != len(self.delta_S_indices_list):
                self.subtour_constraints_to_edge_ec_m.extend([[0] * self.num_edge_variables for _ in range(len(self.delta_S_indices_list) - len(self.subtour_constraints_to_edge_ec_m))])

        if self.augmentation_flag:
            print(f"tour_impa: {self.tour_impa}")
        print(f"cost_impa: {self.cost_impa}")

        if self.exact_solver_flag:
            self.solve_tsp_exact()
        if self.lkh_solver_flag:
            start_time = time.time()
            (
                self.tour_lkh,
                self.cost_lkh,
            ) = self.solve_tsp_lkh()
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
                    self.open_paths,
                    self.runtime_impa,
                    self.runtime_lkh,
                )
            )

    def iterate_augmented_graph(
        self,
    ):
        if self.reset_flag:
            for i in range(len(self.edge_connections)):
                connection = self.edge_connections[i]
                cost = self.cost_matrix[
                    connection[0],
                    connection[1],
                ]
                self.edge_ec_to_degree_constraint_m[i][connection[0]] = cost
                self.edge_ec_to_degree_constraint_m[i][connection[1]] = cost

        print("-----------------------")
        print("iterate_augmented_graph")
        print("-----------------------")
        for iter in range(
            0,
            self.num_iterations,
        ):
            self.model_degree_constraint.degree_constraint_to_edge_ec_update(self.edge_ec_to_degree_constraint_m)

            self.model_degree_constraint.process_filtering(iter)

            self.edge_ec_to_subtour_constraints_m = self.model_eq_constraint.edge_ec_to_subtour_constraints_update(
                self.model_degree_constraint.degree_constraint_to_eq_constraint_m,
                self.subtour_constraints_to_edge_ec_m,
                self.delta_S_indices_list,
                self.cost_edge_variable,
            )

            self.subtour_constraints_to_edge_ec_m_dummy = self.model_subtour_constraint.subtour_constraints_to_edge_ec_update(
                self.edge_ec_to_subtour_constraints_m,
                self.delta_S_indices_list,
            )

            self.subtour_constraints_to_edge_ec_m = self.model_subtour_constraint.process_filtering(iter)

            self.edge_ec_to_degree_constraint_m = self.model_eq_constraint.edge_ec_to_degree_constraint_augmented_graph_update(
                self.model_degree_constraint.degree_constraint_to_eq_constraint_m,
                self.subtour_constraints_to_edge_ec_m,
            )

        extrinsic_output_edge_ec = self.outputs.extrinsic_output_edge_ec_augmented_graph_update(
            self.model_degree_constraint.degree_constraint_to_eq_constraint_m,
            self.subtour_constraints_to_edge_ec_m,
        )
        intrinsic_out_edge_ec = extrinsic_output_edge_ec + self.cost_edge_variable
        self.intrinsic_out_edge_ec = intrinsic_out_edge_ec

    def find_closed_loop_new(
        self,
        graph,
        start,
        current,
        visited,
        path=[],
    ):
        visited.add(current)

        path = path + [current]

        if current == start:
            visited.add(current)
            return path

        if current not in graph:
            self.open_paths_dummy.append([start] + path)
            return None

        for node in graph[current]:
            if node not in visited:
                new_path = self.find_closed_loop_new(
                    graph,
                    start,
                    node,
                    visited.copy(),
                    path,
                )
                if new_path:
                    return new_path

    def get_closed_loops(
        self,
    ):
        self.open_paths_dummy = []
        graph = {}
        for connection in self.selected_edges:
            if connection[0] in graph:
                graph[connection[0]].append(connection[1])
            else:
                graph[connection[0]] = [connection[1]]

        self.graph = graph
        loops_list = []
        self.visited_nodes = []
        for (
            start_node,
            end_node,
        ) in graph.items():
            if start_node in self.visited_nodes:
                continue
            else:
                for j in range(len(end_node)):
                    closed_loop = self.find_closed_loop_new(
                        graph,
                        start_node,
                        end_node[j],
                        set(),
                    )
                    loops_list.append(closed_loop)

        new_loops_list = []

        list_indices_to_remove = []
        for (
            i,
            list_1,
        ) in enumerate(loops_list):
            if list_1 is None:
                list_indices_to_remove.append(i)
                continue
            elif i in list_indices_to_remove:
                continue
            double_list_1 = list_1 + list_1
            for (
                j,
                list_2,
            ) in enumerate(
                loops_list[i + 1 :],
                start=i + 1,
            ):
                if list_2 is None:
                    list_indices_to_remove.append(j)
                elif any(double_list_1[ele : ele + len(list_2)] == list_2 for ele in range(len(list_1))):
                    list_indices_to_remove.append(j)

        new_loops_list = [element for i, element in enumerate(loops_list) if i not in list_indices_to_remove]

        for sublist in new_loops_list:
            sublist.insert(
                0,
                sublist.pop(-1),
            )

        self.loops_list = new_loops_list

        self.open_paths = []
        if len(self.open_paths_dummy):
            open_paths_dummy = self.open_paths_dummy
            paths_indices_to_remove = []
            for (
                index,
                element,
            ) in enumerate(open_paths_dummy):
                size = len(element)
                for (
                    ele,
                    other_element,
                ) in enumerate(open_paths_dummy[:index] + open_paths_dummy[index + 1 :]):
                    if len(other_element) > size and element in [other_element[i : i + size] for i in range(len(other_element) - size + 1)]:
                        paths_indices_to_remove.append(index)
                        break

            self.open_paths = [element for i, element in enumerate(open_paths_dummy) if i not in paths_indices_to_remove]

            print("--------")
            print(f"Number of nodes: {self.num_nodes}")
            for path in self.open_paths:
                print(f"Open Path: {path} with size {len(path)}")

    def subtour_elimination_constraints_analysis(
        self,
    ):
        self.get_closed_loops()

        loops_list = self.loops_list

        if not loops_list:
            print("Exited: Cannot find tours using self.get_closed_loops()")
        elif len(loops_list) == 1 and len(loops_list[0]) == self.num_nodes:
            self.tour_impa_flag = True
            self.subtour_constraints_satisfied_flag = True
            self.tour_impa = loops_list[0] + [loops_list[0][0]]
            print("Tour found")
        else:
            self.subtour_constraints_satisfied_flag = False
            self.tour_impa_flag = False
            for sublist in loops_list:
                print(f"Subtour of size {len(sublist)} detected @: {sublist}")
                pairs = [
                    [
                        sublist[i],
                        sublist[i + 1],
                    ]
                    for i in range(len(sublist) - 1)
                ]
                pairs.append(
                    [
                        sublist[-1],
                        sublist[0],
                    ]
                )
                indices = [self.edge_connections.tolist().index(pair) for pair in pairs if pair in self.edge_connections]
                if np.sum(self.hard_decision[indices]) <= len(sublist) - 1:
                    print("Exited: DETECTION OF SUBTOURS IS WRONG")
                delta_S_indices = [
                    i
                    for i, connection in enumerate(self.edge_connections)
                    if tuple(connection) not in list(
                        itertools.permutations(
                            sublist,
                            2,
                        )
                    )
                    and connection[0] in sublist
                ]
                if delta_S_indices not in self.delta_S_indices_list:
                    self.delta_S_indices_list.append(delta_S_indices)

    def hard_decision_analysis(
        self,
    ):
        intrinsic_out_edge_ec = self.intrinsic_out_edge_ec

        hard_decision = np.array(deepcopy(intrinsic_out_edge_ec))
        hard_decision[intrinsic_out_edge_ec > self.threshold] = 0
        hard_decision[intrinsic_out_edge_ec <= self.threshold] = 1
        self.hard_decision = hard_decision
        self.selected_edges = self.edge_connections[hard_decision == 1]
        self.cost_impa = np.sum(self.cost_edge_variable[hard_decision == 1])
        print(
            "selected_edges: ",
            self.selected_edges.tolist(),
        )
        print(
            "cost_impa: ",
            self.cost_impa,
        )

    def solve_tsp_exact(
        self,
    ):
        print("--------")
        print("Exact Solver: ")
        print("--------")
        (
            tour_exact,
            cost_exact,
        ) = solve_tsp_brute_force(self.cost_matrix)
        print(
            "tour_exact: ",
            tour_exact,
        )
        print(
            "cost_exact: ",
            cost_exact,
        )
        return (
            tour_exact,
            cost_exact,
        )

    def solve_tsp_lkh(
        self,
    ):
        print("--------")
        print("LKH Solver: ")
        print("--------")
        cost_matrix_h = elkai.DistanceMatrix(self.cost_matrix)
        tour_h = cost_matrix_h.solve_tsp()
        cost_h = sum(
            self.cost_matrix[
                tour_h[i],
                tour_h[i + 1],
            ]
            for i in range(len(tour_h) - 1)
        )
        print(
            "tour_h: ",
            tour_h,
        )
        print(
            "cost_h: ",
            cost_h,
        )
        return (
            tour_h,
            cost_h,
        )

    def solve_tsp_simulated_annealing(
        self,
    ):
        print("--------")
        print("Simulated Annealing Solver:")
        print("--------")
        (
            tour_h,
            cost_h,
        ) = solve_tsp_simulated_annealing(self.cost_matrix)
        print(
            "tour_h: ",
            tour_h,
        )
        print("cost_h: {:.4f}".format(cost_h))

    def save_outputs(
        self,
    ):
        with open(
            f"{self.folder_outputs}/outputs_set{self.test_file}.pkl",
            "wb",
        ) as f:
            pkl.dump(
                self.results_composed,
                f,
            )
