# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

# Import necessary modules and functions from the environmentModule
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
    random
)

# Import necessary variables and functions from the initializationModule
from impa.initializationModule import (
    fighters_count,
    weasels_count,
    F_u,
    W_u,
    distance_metric,
    distance_metric_rendezvouz,
    K_1,
)

# Import necessary types and function prototypes from the cFunctionAPI
from impa.cFunctionAPI import (
    c_int_p,
    c_impa_lib_type_p,
    c_bool_p,
    WrapperTsp,
    WrapperKcMwm,
    WrapperKsat,
)

"""
Graphical Model class for the TSP

Attributes:
    NUM_ITERATIONS: number of iterations of IMPA
    NUM_NODES: number of nodes of TSP
    SYMMETRIC_FLAG: symmetric or asymmetric TSP
    AUGMENTATION_FLAG: perform augmentation or not
    EXACT_SOLVER_FLAG: run exact solver on TSP
    LKH_SOLVER_FLAG: run LKH solver on TSP
    SIM_ANNEALING_FLAG: run simulated annealing on TSP
    RESET_FLAG: reset messages after each augmentation step
    FILTERING_FLAG: filtering messages of degree and subtour elimination constraint
    ALPHA: filtering parameter (between 0 and 1)
    THRESHOLD: threshold for hard-decision of IMPA (usually a negative value close to zero)
    RANDOM_TEST_FLAG: run IMPA on random TSP
    MAX_COUNT: maximum count of failure cases
    POST_PROCESS_FLAG: run post-processing
    K_OPT_FLAG: run k-opt on post-processed solution
    MAX_AUGM_COUNT: maximum number of augmentations
"""
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
        """
        Initialize Graphical Model for the TSP
        """

        # Initialize variables and flags
        self.results_composed = []

        self.num_augmentations = 0
        self.num_added_constraints = 0
        self.no_improvement_sol_count_exceeded_flag = False
        self.no_consecutive_loops_count_exceeded_flag = False
        self.sol_oscillation_count_exceeded_flag = False

        # If not performing random testing, load data from input_load
        if not self.random_test_flag:
            self.num_nodes = self.input_load[0]
            self.symmetric_flag = self.input_load[1]

        # Extract necessary parameters
        num_nodes = self.num_nodes
        symmetric_flag = self.symmetric_flag

        print(f"num_nodes: {self.num_nodes}")
        print(f"symmetric_flag: {self.symmetric_flag}")

        # Create cost matrix and edge connections
        cost_matrix, edge_connections = self.create_cost_matrix(symmetric_flag)
        self.tour_impa_flag = False
        self.tour_impa = np.zeros(num_nodes + 1, dtype=np.int32)  # +1 to account for returning edge
        self.cost_impa = np.zeros(1, dtype=np_impa_lib)
        self.selected_edges = np.zeros(math.comb(self.num_nodes, 2), dtype=np.int32)
        self.subtour_paths = np.zeros(self.num_nodes**2, dtype=np.int32)
        self.open_paths = np.zeros(self.num_nodes**2, dtype=np.int32)

        # If not performing random testing, use loaded data for cost_matrix and edge_connections
        if not self.random_test_flag:
            self.cost_matrix = self.input_load[2]
            self.edge_connections = self.input_load[3]
            cost_matrix = self.cost_matrix
            edge_connections = self.edge_connections

        # Assign cost_matrix and edge_connections
        self.cost_matrix = cost_matrix
        self.edge_connections = edge_connections

        # Calculate number of edge variables
        num_edge_variables = num_nodes**2 - num_nodes

        self.num_edge_variables = num_edge_variables

        # Initialize arrays for edge variables
        self.intrinsic_out_edge_ec = np.zeros(num_edge_variables, dtype=np_impa_lib)

        cost_edge_variable = np.zeros(num_edge_variables, dtype=np_impa_lib)

        edge_ec_to_degree_constraint_m = np.zeros((num_edge_variables, num_nodes), dtype=np_impa_lib)

        # Populate arrays based on edge connections and cost matrix
        for i in range(len(edge_connections)):
            connection = edge_connections[i]
            cost = cost_matrix[connection[0], connection[1]]
            cost_edge_variable[i] = cost
            edge_ec_to_degree_constraint_m[i][connection[0]] = cost
            edge_ec_to_degree_constraint_m[i][connection[1]] = cost

        # Assign calculated arrays
        self.cost_edge_variable = cost_edge_variable
        self.edge_degree_constraint_cost = edge_ec_to_degree_constraint_m

        self.edge_ec_to_degree_constraint_m = edge_ec_to_degree_constraint_m

    def create_cost_matrix(self, symmetric_flag):
        """
        Create cost matrix for the TSP

        Args:
            symmetric_flag: symmetric TSP or asymmetric

        Returns:
           matrix: cost matrix
           edge_connections: list of edges. Each element consist of the nodes of the edge
        """
        # Generate a random matrix for the upper triangle
        n = self.num_nodes
        upper_triangle = np.random.uniform(1, 1000, size=(n, n))

        # Make the matrix symmetric if symmetric_flag is True
        if symmetric_flag:
            matrix = np.triu(upper_triangle) + np.triu(upper_triangle, k=1).T
            np.fill_diagonal(matrix, zero_value)
        else:
            # Otherwise, keep it asymmetric
            matrix = upper_triangle.copy()
            np.fill_diagonal(matrix, zero_value)
        # Find edge connections from the cost matrix
        edge_connections = np.transpose(np.nonzero(matrix))
        return matrix, edge_connections

    def run_impa(self):
        """
        Run impa in the C++ code
        """
        # Process inputs and get ctypes pointers
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

        # Call the IMPA wrapper function with the necessary arguments
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

        # Process outputs from ctypes pointers
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
        """
        This function will process the inputs to the C++ code

        Returns:
            edge_connection_flatten_p: a pointer for flattened edge connections (edge connections contains lists, each list has the nodes of an edge)
            cost_edge_variable_flatten_p: a pointer for flattened cost of each edge equality constraint
            cost_matrix_flatten_p: a pointer for a flattened cost matrix
            edge_ec_to_degree_constraint_m_flatten_p: a pointer for flattened messages from edge equality constraints to degree constraints
            edge_degree_constraint_cost_flatten_p: a pointer for flattened edge degree constraint cost (this cost before flattening has size function of number of nodes and number of edges)
            extrinsic_output_edge_ec_p: a pointer that will store extrinsic output messages from edge equality constraints which are calculated in C++
            num_augmentations_p: a pointer that will store the number of augmentations of IMPA obtained in C++
            tour_impa_flatten_p: a pointer that will store the elements of a tour (if found) obtained in C++
            cost_impa_p: a pointer that will store the cost of the activated edges (could be a tour) obtained in C++
            selected_edges_p: a pointer that will store the activated edges obtained after running IMPA in C++
            selected_edges_size_p: a pointer that will store the number of activated edges after running IMPA in C++
            no_improvement_sol_count_exceeded_flag_p: a pointer that will store a failure in IMPA due to exceeding the maximum count of no improvement in the IMPA solution
            no_consecutive_loops_count_exceeded_flag_p: a pointer that will store a failure in IMPA due to exceeding the maximum count of no consecutive detection of loops in IMPA (when no tour is found)
            sol_oscillation_count_exceeded_flag_p: a pointer that will store a failure in IMPA due to exceeding the maximum count of detection of consecutive oscillation in the IMPA solution
            subtour_paths_flatten_p: a pointer that will store detected subtours after running the IMPA in C++
            subtour_paths_size_flatten_p: a pointer that will store the number of detected subtours after running the IMPA in C++
            num_added_constraints_p: a pointer that will stire the number of added constraints after running the IMPA in C++
        """
        # Flatten and get pointers
        edge_connection_flatten = self.edge_connections.flatten().astype(np.int32)
        edge_connection_flatten_p = edge_connection_flatten.ctypes.data_as(c_int_p)
        cost_edge_variable_flatten = self.cost_edge_variable.flatten().astype(np_impa_lib)
        cost_edge_variable_flatten_p = cost_edge_variable_flatten.ctypes.data_as(c_impa_lib_type_p)
        cost_matrix_flatten = self.cost_matrix.flatten().astype(np_impa_lib)
        cost_matrix_flatten_p = cost_matrix_flatten.ctypes.data_as(c_impa_lib_type_p)
        edge_ec_to_degree_constraint_m_flatten = self.edge_ec_to_degree_constraint_m.flatten().astype(np_impa_lib)
        edge_ec_to_degree_constraint_m_flatten_p = edge_ec_to_degree_constraint_m_flatten.ctypes.data_as(c_impa_lib_type_p)
        edge_degree_constraint_cost_flatten = self.edge_degree_constraint_cost.flatten().astype(np_impa_lib)
        edge_degree_constraint_cost_flatten_p = edge_degree_constraint_cost_flatten.ctypes.data_as(c_impa_lib_type_p)
        extrinsic_output_edge_ec = np.zeros((self.num_edge_variables))
        extrinsic_output_edge_ec = np.array(extrinsic_output_edge_ec).flatten().astype(np_impa_lib)
        extrinsic_output_edge_ec_p = extrinsic_output_edge_ec.ctypes.data_as(c_impa_lib_type_p)

        # Get pointers for other variables
        num_augmentations_p = np.array(self.num_augmentations).ctypes.data_as(c_int_p)
        num_added_constraints_p = np.array(self.num_added_constraints).ctypes.data_as(c_int_p)
        tour_impa_flatten = self.tour_impa.flatten().astype(np.int32)
        tour_impa_flatten_p = tour_impa_flatten.ctypes.data_as(c_int_p)
        cost_impa_p = np.array(self.cost_impa).ctypes.data_as(c_impa_lib_type_p)
        selected_edges_flatten = np.array(self.selected_edges).flatten().astype(np.int32)
        selected_edges_p = selected_edges_flatten.ctypes.data_as(c_int_p)
        selected_edges_size_p = np.array(0).ctypes.data_as(c_int_p)

        no_improvement_sol_count_exceeded_flag_p = np.array(self.no_improvement_sol_count_exceeded_flag).ctypes.data_as(c_bool_p)
        no_consecutive_loops_count_exceeded_flag_p = np.array(self.no_consecutive_loops_count_exceeded_flag).ctypes.data_as(c_bool_p)
        sol_oscillation_count_exceeded_flag_p = np.array(self.sol_oscillation_count_exceeded_flag).ctypes.data_as(c_bool_p)

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
        """
        This function will process the outputs of the C++ code

        Args:
            sol_oscillation_count_exceeded_flag_p: a pointer that stores a failure in IMPA due to exceeding the maximum count of detection of consecutive oscillation in the IMPA solution
            selected_edges_size_p: a pointer that stores the number of activated edges after running IMPA in C++
            selected_edges_p: a pointer that stores the activated edges obtained after running IMPA in C++
            subtour_paths_size_flatten_p: a pointer that stores the number of detected subtours after running the IMPA in C++
            subtour_paths_flatten_p: a pointer that stores detected subtours after running the IMPA in C++
            no_improvement_sol_count_exceeded_flag_p: a pointer that stores a failure in IMPA due to exceeding the maximum count of no improvement in the IMPA solution
            no_consecutive_loops_count_exceeded_flag_p: a pointer that stores failure in IMPA due to exceeding the maximum count of no consecutive detection of loops in IMPA (when no tour is found)
            tour_impa_flatten_p: a pointer that stores the elements of a tour (if found) obtained in C++
            num_augmentations_p: a pointer that stores the number of augmentations of IMPA obtained in C++
            cost_impa_p: a pointer that stores the cost of the activated edges (could be a tour) obtained in C++
            extrinsic_output_edge_ec_p: a pointer that stores extrinsic output messages from edge equality constraints which are calculated in C++
            num_added_constraints_p: a pointer that stores the number of added constraints after running the IMPA in C++
        """
        # Extract values from pointers
        selected_edges_size = list(selected_edges_size_p.__dict__.values())[0]
        selected_edges_flatten = list(selected_edges_p.__dict__.values())[0][:selected_edges_size]
        self.selected_edges = [[selected_edges_flatten[i], selected_edges_flatten[i + 1]] for i in range(0, len(selected_edges_flatten), 2)]

        # Extract subtour paths
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

        # Update flags and tour
        self.sol_oscillation_count_exceeded_flag = list(sol_oscillation_count_exceeded_flag_p.__dict__.values())[0]
        self.no_improvement_sol_count_exceeded_flag = list(no_improvement_sol_count_exceeded_flag_p.__dict__.values())[0]
        self.no_consecutive_loops_count_exceeded_flag = list(no_consecutive_loops_count_exceeded_flag_p.__dict__.values())[0]
        tour_impa = list(tour_impa_flatten_p.__dict__.values())[0]
        if all(element == 0 for element in tour_impa):
            self.tour_impa = None
        else:
            self.tour_impa = [x for x in tour_impa]

        # Update augmentation and added constraints information
        self.num_augmentations = int(list(num_augmentations_p.__dict__.values())[0])
        if self.augmentation_flag:
            self.num_added_constraints = int(list(num_added_constraints_p.__dict__.values())[0])
        else:
            self.num_added_constraints = 0
        # Update cost information
        self.cost_impa = list(cost_impa_p.__dict__.values())[0][0]
        extrinsic_output_edge_ec = list(extrinsic_output_edge_ec_p.__dict__.values())[0]
        intrinsic_out_edge_ec = extrinsic_output_edge_ec + self.cost_edge_variable
        self.intrinsic_out_edge_ec = intrinsic_out_edge_ec
        if self.subtour_paths is not None:
            self.subtour_paths = [sublist + [sublist[0]] for sublist in self.subtour_paths if sublist]
        # If subtours detected and post-processing flag is set, run post-processing
        if self.tour_impa is None and self.post_process_flag:
            self.run_post_processing_improved()
            self.pp_performed = True
        else:
            self.pp_performed = False

    def undirected_graph(self, selected_edges):
        """
        This function will create an undirected graph based on the activated edges. This function is used in the pure python code

        Args:
            selected_edges: activated edges after running IMPA
        Returns:
            graph: undirected graph that contains the connection between the nodes
        """
        # Initialize an empty dictionary to represent the graph
        graph = {}

        # Iterate over each activated edge
        for connection in selected_edges:
            # Add the first node of the edge to the graph
            if connection[0] in graph:
                graph[connection[0]].append(connection[1])
            else:
                graph[connection[0]] = [connection[1]]
            # Add the second node of the edge to the graph
            if connection[1] in graph:
                graph[connection[1]].append(connection[0])
            else:
                graph[connection[1]] = [connection[0]]
        # Return the constructed undirected graph
        return graph

    def paths_cleaning(self, removed_edges, investigated_paths):
        """
        This function will process the subtour paths. This is used during post-processing where some edges will be removed, and new paths will be obtained after removing edges

        Args:
            removed_edges: edges to be removed from investigated_paths
            investigated_paths: paths to investigate and remove edges from
        Returns:
            final_investigated_paths: paths obtained after removing desired edges
            removed_edge_flag_list: list that indicates where an investigated path inside investigated_paths has been disconnected
        """
        # If investigated_paths is None, return None for both outputs
        if investigated_paths is None:
            return None, None
        # Convert removed_edges to a list of tuples for easy comparison
        removed_edges = [tuple(edge) for edge in removed_edges]
        # Initialize lists to store modified paths and flags indicating removed edges
        new_investigated_paths = []
        removed_edge_flag_list = []

        # Iterate over each path in investigated_paths
        for path in investigated_paths:
            removed_edge_flag = False  # Flag to indicate if any edge in the path is removed
            new_path = []  # List to store the modified path
            i = 0
            # Iterate over length of the path
            while i < len(path) - 1:
                # Check if the edge is in removed_edges
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
            # Remove any duplicate paths from new_investigated_paths
        final_investigated_paths = [path for i, path in enumerate(new_investigated_paths) if not any(set(path).issubset(p) for p in new_investigated_paths[:i] + new_investigated_paths[i + 1 :])]
        # Return the final lists
        return final_investigated_paths, removed_edge_flag_list

    def construct_graph(self, selected_edges):
        """
        This function will create a graph of the avtivated edges. This will be used for the detection of paths that will be connected to form a tour

        Args:
            selected_edges: activated edges forming the solution of IMPA
        Returns:
            graph: graph that will be used during post-processing for forming a tour of the TSP
        """
        # Initialize an empty dictionary to represent the graph
        graph = {}
        # Iterate over each activated edge
        for connection in selected_edges:
            if connection[0] in graph:
                graph[connection[0]][1].append(connection[1])
            else:
                graph[connection[0]] = ([], [connection[1]])

            if connection[1] in graph:
                graph[connection[1]][0].append(connection[0])
            else:
                graph[connection[1]] = ([connection[0]], [])
        # Return the constructed graph
        return graph

    def run_post_processing_improved(self):
        """
        This function will run post-processing on the solution obtained from running IMPA in C++
        """
        print("--------")
        print("Post Processing Started:")
        print("--------")
        # Retrieve necessary data
        num_nodes = self.num_nodes
        # edge_connections = self.edge_connections
        # cost_edge_variable = self.cost_edge_variable
        subtour_paths = self.subtour_paths
        # Print subtour paths if available
        if subtour_paths is not None:
            for i, subtour_path in enumerate(subtour_paths):
                print(f"Subtour path of size {len(subtour_path)-1}: {subtour_path}")
        print("--------")

        # Create a copy of selected edges
        self.selected_edges_pp = copy.deepcopy(self.selected_edges)
        selected_edges_pp = self.selected_edges_pp[:]
        degree_constraints_satisfied = False

        # Construct the graph
        graph = self.construct_graph(selected_edges_pp)

        # Check degree constraints
        degree_constraints_satisfied = True
        keys = []
        for node, (inward_edges_nodes, outward_edges_nodes) in graph.items():
            if len(inward_edges_nodes) > 1 or len(outward_edges_nodes) > 1:
                keys.append(node)
            if len(inward_edges_nodes) != 1 or len(outward_edges_nodes) != 1:
                degree_constraints_satisfied = False

        # Investigate violated degree constraints with degree > 1
        print(f"Degree Constraints Satisfied? {degree_constraints_satisfied}")
        removed_edges_violated_dc = []
        if len(keys) != 0:
            print("Investigating Violated Degree Constraints with degree >1")
            for key in keys:
                inward_edges_nodes = graph[key][0]
                outward_edges_nodes = graph[key][1]
                print(f"Node: {key}, Inward Edges Nodes: {inward_edges_nodes}, Outward Edges Nodes: {outward_edges_nodes}")
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
                                    (i for i, edge in enumerate(self.edge_connections) if np.array_equal(edge, investigated_edge)),
                                    None,
                                )
                                indices_in_edge_connections.append(index_edge_connections)
                        # python3 main_wrapper_tsp.py --testFile=36 --alpha=0.5 --nITER=200 --inputPath=inputs_random_1000 --outputPath=outputs_random_1000_maxAugmCount50 --saveFlag=True --augmFlag=True --maxAugmCount=50 --filteringFlag=True
                        if indices_in_edge_connections:  # added since experienced a failure scenario where all outward edges where already removed
                            # flag_emptied_connection = False
                            # valid_indices = [
                            #    index for index in indices_in_edge_connections
                            # ]
                            min_index = indices_in_edge_connections[np.argmin(np.abs(self.intrinsic_out_edge_ec[indices_in_edge_connections]))]
                            removed_edge = list(self.edge_connections[min_index])
                            removed_edge_cost = self.cost_edge_variable[min_index]
                            if removed_edge not in removed_edges_per_node:
                                removed_edges_per_node.append(removed_edge)
                            if removed_edge not in removed_edges_violated_dc:
                                print(f"Remove: {removed_edge} with cost: {removed_edge_cost}")  # or {self.cost_matrix[removed_edge[0]][removed_edge[1]]}
                                removed_edges_violated_dc.append(removed_edge)
                                indices_in_edge_connections.remove(min_index)
                            else:
                                print(f"Edge {removed_edge} with cost: {removed_edge_cost} already considered")
                        else:
                            print("index_edge_connections empty. All edges of violation were already removed")
                            break
        else:
            print("No Violated Degree Constraints with degree >1")

        # Clean subtour paths
        subtour_paths, removed_edge_flag_list = self.paths_cleaning(removed_edges_violated_dc, self.subtour_paths)
        print(f"subtour_paths after removing removed_edges_violated_dc {removed_edges_violated_dc}: \n {subtour_paths}")
        print("------")

        # Investigate connected subtour paths
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
            print(f"subtour_paths after removing removed_subtour_paths_edges {removed_subtour_paths_edges}: \n {subtour_paths}")

            selected_edges_pp = [edge for edge in selected_edges_pp if edge not in removed_subtour_paths_edges]

        selected_edges_pp = [edge for edge in selected_edges_pp if edge not in removed_edges_violated_dc]

        print("------")

        # Construct a graph from the updated selected edges
        graph = self.construct_graph(selected_edges_pp)

        # Get final paths from the updated graph
        final_paths = self.get_paths(graph)

        # Check degree constraints satisfaction
        degree_constraints_satisfied = True
        for node, (inward_edges_nodes, outward_edges_nodes) in graph.items():
            if len(inward_edges_nodes) != 1 or len(outward_edges_nodes) != 1:
                degree_constraints_satisfied = False

        print("Degree Constraints Satisfied?", degree_constraints_satisfied)

        # Print final paths and check for missing nodes and duplicate nodes
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

        # If there are multiple final paths, connect them to form a tour
        if len(final_paths) > 1:
            print("------")
            print("Connecting Paths to form a tour")
            lst = range(0, len(final_paths))
            perms = permutations(lst)
            results = [perm for perm in perms if perm[0] == 0]

            possible_tours = []
            cost_tours = []

            # Investigate different connections to find the optimal tour
            for i, result in enumerate(results):
                print(f"Investigating connection {i} out of {len(results)} connections")
                tour = []
                count_used_paths = 0
                while count_used_paths < len(final_paths):
                    tour.extend(final_paths[result[count_used_paths]])
                    count_used_paths += 1

                # Add missing nodes if any
                if len(missing_nodes) != 0:
                    tour = self.add_missing_nodes(tour, missing_nodes)

                # Calculate the cost of the tour
                pair_occurrences = self.find_indices_edges_from_tour(tour)
                possible_tours.append(tour)
                cost_tours.append(sum(self.cost_edge_variable[pair_occurrences]))

            # Choose the tour with the minimum cost
            index_best_tour = np.argmin(cost_tours)
            tour_impa_pp = possible_tours[index_best_tour]
            cost_impa_pp = cost_tours[index_best_tour]
            min_cost = cost_impa_pp
            best_tour = tour_impa_pp[:]

            # Perform k-opt optimization if enabled
            if self.k_opt_flag:
                print(f"tour_impa_pp: {tour_impa_pp+ [possible_tours[index_best_tour][0]]}")
                print(f"cost_impa_pp: {cost_impa_pp}")
                best_tour, min_cost = self.perform_k_opt(best_tour, min_cost)

            # Finalize the tour
            tour_impa_pp = best_tour[:]
            tour_impa_pp = tour_impa_pp + [tour_impa_pp[0]]
            cost_impa_pp = min_cost
            print(f"tour_impa_pp of size {len(tour_impa_pp)-1}: {tour_impa_pp}")
            print(f"cost_impa_pp: {cost_impa_pp}")
        # If there is only one final path, form a tour from it
        else:
            print("------")
            print("Forming a tour from one path")
            if len(missing_nodes) != 0:
                print("Adding Missing Nodes")
                tour = self.add_missing_nodes(final_paths[0], missing_nodes)
            else:
                tour = final_paths[0]

            # Calculate the cost of the tour
            pair_occurrences = self.find_indices_edges_from_tour(tour)
            cost_impa_pp = sum(self.cost_edge_variable[pair_occurrences])
            tour_impa_pp = tour[:]

            # Perform k-opt optimization if enabled
            if self.k_opt_flag:
                print(f"tour_impa_pp: {tour_impa_pp + [tour_impa_pp[0]]}")
                print(f"cost_impa_pp: {cost_impa_pp}")
                best_tour, min_cost = self.perform_k_opt(tour_impa_pp, cost_impa_pp)
            else:
                best_tour = tour_impa_pp[:]
                min_cost = cost_impa_pp

            # Finalize the tour
            tour_impa_pp = best_tour[:]
            tour_impa_pp = tour_impa_pp + [tour_impa_pp[0]]
            cost_impa_pp = min_cost
            print(f"tour_impa_pp of size {len(tour_impa_pp)-1}: {tour_impa_pp}")
            print(f"cost_impa_pp: {cost_impa_pp}")

        # Update the cost and tour attributes
        self.cost_impa = cost_impa_pp
        self.tour_impa = tour_impa_pp

    def perform_k_opt(self, best_tour, min_cost):
        """
        This function will perform k-opt on the solution obtained from post-processing

        Args:
            best_tour: best tour found after running post-processing on the IMPA solution
            min_cost: minimum cost of the tour found after running post-processing on the IMPA solution
        Returns:
            best_tour: best tour found after running the k-opt algorithm
            min_cost: best cost found after running the k-opt algorithm
        """
        print("------")
        print("Performing K-OPT")
        # Iterate over each value of k in [3, 2]
        for k in [3, 2]:
            print(f"{k}-opt")
            # Generate candidate paths using k-opt
            candidate_paths = self.generate_k_opt_paths(best_tour, k)
            candidate_path = candidate_paths[0]
            # Evaluate each candidate path
            for candidate_path in candidate_paths:
                # Calculate the cost of the candidate path
                pair_occurrences = self.find_indices_edges_from_tour(candidate_path)
                cost = sum(self.cost_edge_variable[element] for element in pair_occurrences)
                # Update minimum cost and best tour if the cost is lower
                if cost < min_cost:
                    min_cost = cost
                    best_tour = candidate_path[:]
        return best_tour, min_cost

    def add_missing_nodes(self, tour, missing_nodes):
        """
        This function will add missing nodes on paths obtained after doing post-processing since post-processing could eliminate a node from the graph

        Args:
            tour: connection of path that will form a tour
            missing_nodes: missing nodes from the graph
        Returns:
            tour: tour obtained after adding the missing nodes to the connected paths
        """
        added_missing_nodes = []
        # Loop until all missing nodes are added to the tour
        while len(added_missing_nodes) != len(missing_nodes):
            # Generate combinations of missing nodes with the start and end nodes of the tour
            missing_combinations_left = [[node, tour[0]] for node in missing_nodes if node not in added_missing_nodes]
            missing_combinations_right = [[tour[-1], node] for node in missing_nodes if node not in added_missing_nodes]

            # Initialize minimum costs for left and right combinations
            min_cost_left = np_impa_lib("inf")
            min_cost_right = np_impa_lib("inf")

            # Find the combination with the minimum cost for the left and right connections
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
            # Determine which combination has the minimum cost
            costs_left_right = [min_cost_left, min_cost_right]
            # added_cost = np.min(costs_left_right)
            added_cost_index = np.argmin(costs_left_right)

            # Add the corresponding missing node to the tour
            if added_cost_index == 0:  # left
                added_missing_nodes.append(combination_left[0])
                tour.insert(0, combination_left[0])
            else:  # right
                added_missing_nodes.append(combination_right[1])
                tour.append(combination_right[1])
        return tour

    def find_indices_edges_from_tour(self, tour):
        """
        This function will find the indices of edges forming the tour. This will be used in the calculation of the cost of the tour after running post-processing

        Args:
            tour: tour which consists of connected edges
        Returns:
            pair_occurrences: edges (or pair of nodes) forming the tour
        """
        # Create pairs of nodes representing edges in the tour
        pairs = np.column_stack((tour, np.roll(tour, -1)))
        edge_connections_arr = np.array(self.edge_connections)
        # Check for matching pairs (edges) in the list of edge connections
        matching_pairs = (edge_connections_arr[:, None, :] == pairs).all(axis=-1)
        # Get the indices of matching pairs
        pair_occurrences = np.where(matching_pairs)[0]
        return pair_occurrences

    def find_paths(self, graph, node, visited, path=[]):
        """
        This function will find paths in a graph

        Args:
            graph: graph of the activated edges
            node: node from which to start finding paths
            visited: list of nodes that are already visited
        Returns:
            path: detected path
        """
        # Mark the current node as visited
        visited.add(node)

        # Append the current node to the visited nodes list
        self.visited_nodes.append(node)

        # Add the current node to the current path
        path = path + [node]

        # If the current node has no outgoing edges, add the path to the list of open paths and return
        if len(graph[node][1]) == 0:
            self.open_paths_dummy.append(path)
            return path

        # Recursively explore neighbors of the current node
        for neighbor in graph[node][1]:
            if neighbor not in visited:
                # Recursively find paths starting from the neighbor node
                new_path = self.find_paths(graph, neighbor, visited.copy(), path.copy())
                # If a new path is found, return it
                if new_path:
                    return new_path

    def get_paths(self, graph):
        """
        This function will call find_paths recursively to find all paths in a graph

        Args:
            graph: graph of the activated edges
        Returns:
            final_paths: detected paths in the graph. These paths will be used in the post-processing
        """
        # Initialize variables
        self.open_paths_dummy = []
        path_lists = []
        self.visited_nodes = []
        # Iterate through all nodes in the graph
        for start_node in graph.keys():
            # If the start node has already been visited, skip it
            if start_node in self.visited_nodes:
                continue
            else:
                # Call the find_paths function to find paths starting from the current node
                path = self.find_paths(graph, start_node, set(), [])
                # If a path is found, add it to the list of path_lists
                if path:
                    path_lists.append(path)
        # Remove duplicate paths from path_lists
        final_paths = [path for i, path in enumerate(path_lists) if not any(set(path).issubset(p) for p in path_lists[:i] + path_lists[i + 1 :])]
        return final_paths

    def generate_k_opt_paths(self, path, k):
        """
        This function will generate the paths after running the k-opt algorithm

        Args:
            path: path on which to perform k-opt algorithm
            k: parameter of the k-opt algorithm (usually 2 or 3)
        Returns:
            new_paths: paths after running the k-opt algorithm
        """
        n = len(path)
        # Check if the value of k is within valid range
        if k < 2 or k > n - 2:
            raise ValueError("k must be between 2 and n-2 for a valid k-opt operation.")

        # Generate all combinations of k indices from the path
        swap_combinations = np.array(list(combinations(range(n), k)), dtype=np.int32)
        new_paths = []

        # Iterate over each combination of k indices
        for swap_indices in swap_combinations:
            new_path = np.array(path, dtype=np.int32)
            # Generate all possible pairs of indices to swap
            swap_pairs = np.array(list(combinations(swap_indices, 2)), dtype=np.int32)
            # Swap the elements at the selected indices
            new_path[swap_pairs[:, 0]], new_path[swap_pairs[:, 1]] = (
                new_path[swap_pairs[:, 1]],
                new_path[swap_pairs[:, 0]],
            )
            new_paths.append(new_path.tolist())

        return new_paths

    def run_pre_analysis(self):
        """
        This function will run an analysis for the TSP. Other known algorithms could be evaluated on the TSP, and the IMPA solution could be stored
        """
        # Calculate runtime for IMPA
        self.runtime_impa = time.time() - self.start_time
        print(f"IMPA Time: {self.runtime_impa}")
        print("--------")

        # Solve TSP using exact solver if enabled
        if self.exact_solver_flag:
            self.solve_tsp_exact()
        # Solve TSP using LKH solver if enabled
        if self.lkh_solver_flag:
            start_time = time.time()
            self.tour_lkh, self.cost_lkh = self.solve_tsp_lkh()
            self.runtime_lkh = time.time() - start_time
            print(f"LKH Time: {self.runtime_lkh}")
        # Solve TSP using simulated annealing if enabled
        if self.sim_annealing_flag:
            self.solve_tsp_simulated_annealing()
        # Store results if save flag is enabled
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
        """
        This function is used in finding loops in a graph in the pure python code

        Args:
            graph: graph of the activated edges
            start: node from which to start finding paths
            current: current node under investigation
            visited: list of visited nodes
            path: detected path
        Returns:
            path: closed loop if found
        """

        # Add the current node to the visited set
        visited.add(current)

        # Append the current node to the path
        self.visited_nodes_loop.append(current)

        path = path + [current]

        # If the current node is the start node, a loop is found
        if current == start:
            visited.add(current)
            return path

        # If the current node is not in the graph, return None
        if current not in graph:
            return None

        # Explore neighbors of the current node
        for node in graph[current]:
            # If the neighbor node has not been visited, explore it
            if node not in visited:
                # Recursively call find_loops with the neighbor node
                new_path = self.find_loops(graph, start, node, visited.copy(), path)
                # If a loop is found in the recursion, return the path
                if new_path:
                    return new_path

    def get_loops(self, selected_edges, path_length):
        """
        This function will call find_loops recursively to find all closed loops in a path. This function is used in the pure python code.

        Args:
            selected_edges: activated edges of the IMPA solution
            path_length: length of the path to compare it to the detected closed loop length
        Returns:
            closed_loop: detected closed loop
        """
        # Initialize an empty graph dictionary
        graph = {}
        # Populate the graph dictionary with the selected edges
        for connection in selected_edges:
            if connection[0] in graph:
                graph[connection[0]].append(connection[1])
            else:
                graph[connection[0]] = [connection[1]]

        self.graph = graph
        # Initialize an empty list to store visited nodes during loop detection
        self.visited_nodes_loop = []
        # Iterate through each start node in the graph
        for start_node, end_node in graph.items():
            # If the start node has already been visited, skip it
            if start_node in self.visited_nodes_loop:
                continue
            else:
                # Iterate through each end node connected to the start node
                for j in range(len(end_node)):
                    # Find closed loops starting from the current start node and end node
                    closed_loop = self.find_loops(graph, start_node, end_node[j], set())
                    # If a closed loop of the specified path length is found, return it
                    if closed_loop is not None and len(closed_loop) == path_length:
                        return closed_loop

    def solve_tsp_exact(self):
        """
        This function will run the exact solver on the TSP
        """
        print("--------")
        print("Exact Solver: ")
        print("--------")
        # Calling the brute force solver
        tour_exact, cost_exact = solve_tsp_brute_force(self.cost_matrix)
        print("tour_exact: ", tour_exact)
        print("cost_exact: ", cost_exact)
        return tour_exact, cost_exact

    def solve_tsp_lkh(self):
        """
        This function will run LKH solver on the TSP
        """
        print("--------")
        print("LKH Solver: ")
        print("--------")
        # Creating a distance matrix for LKH solver
        cost_matrix_h = elkai.DistanceMatrix(self.cost_matrix)
        # Solving TSP using LKH solver
        tour_h = cost_matrix_h.solve_tsp()
        # Calculating the cost of the tour
        cost_h = sum(self.cost_matrix[tour_h[i], tour_h[i + 1]] for i in range(len(tour_h) - 1))
        print("tour_h: ", tour_h)
        print("cost_h: ", cost_h)
        return tour_h, cost_h

    def solve_tsp_simulated_annealing(self):
        """
        This function will run simulated annealing solver on the TSP
        """
        print("--------")
        print("Simulated Annealing Solver:")
        print("--------")
        # Solving TSP using simulated annealing
        tour_h, cost_h = solve_tsp_simulated_annealing(self.cost_matrix)
        print("tour_h: ", tour_h)
        print("cost_h: {:.4f}".format(cost_h))

    def save_outputs(self):
        """
        This function will save results_composed (a dictionary containing variables of interest)
        """
        # Saving the results to a file
        with open(f"{self.folder_outputs}/outputs_set{self.test_file}.pkl", "wb") as f:
            pkl.dump(self.results_composed, f)


"""
Graphical Model class for the Knapsack-MWM problem

Attributes:
    NUM_ITERATIONS: number of iterations of IMPA
    FILTERING_FLAG: filtering messages of degree and subtour elimination constraint
    POST_PROCESS_FLAG: run post-processing
    ALPHA: filtering parameter (between 0 and 1)
    THRESHOLD: threshold for hard-decision of IMPA (usually a negative value close to zero)
    PP_OPTION: perform post-processing on teams or on departments
"""
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
        """
        Initialize Graphical Model for the Knapsack-MWM problem

        Args:
            input_load: inputs of the graphical model (input_load[0] contains capacities of departments, input_load[1] contains the types of teams, input_load[3][0] contains the number of projects)
        """
        # Extracting inputs from input_load
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

        # Pruning teams and generating team rewards
        self.prune_teams()

        (
            reward_team,
            teams_weights_per_department,
            teams_types_per_department,
            last_index,
        ) = self.team_reward_generation()

        # Setting additional instance variables
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

        # Initializing intrinsic_out_mwm
        self.intrinsic_out_mwm = np.zeros((self.num_projects, self.num_teams))

        # Checking if the problem is unbalanced
        if self.num_projects != self.num_teams:
            self.unbalanced_flag = True
        else:
            self.unbalanced_flag = False

        # Assertions for data consistency
        assert self.intrinsic_out_mwm.shape[0] == self.reward_project.shape[0], "Row Shape mismatch between MWM and Rewards"
        assert self.intrinsic_out_mwm.shape[1] == self.reward_project.shape[1], "Column Shape mismatch between MWM and Rewards"

    def prune_teams(self):
        """
        This function will prune teams in the graphical model based on constraints obtained from Caliola
        """
        units = range(1, len(self.max_state) + 1)

        permutations_units = [p for p in itertools.product(units, repeat=2)]
        combinations_units_pre_pruning = list(combinations(units, 2))

        list_pruning = []

        # Evaluate combinations and prune based on specified conditions
        for i in range(0, len(combinations_units_pre_pruning)):
            combination = combinations_units_pre_pruning[i]
            r_1 = math.exp(-F_u[combination[0] - 1]) * math.exp(-W_u[combination[1] - 1])
            r_2 = math.exp(-F_u[combination[1] - 1]) * math.exp(-W_u[combination[0] - 1])

            # Prune combinations based on ratios and distance metrics
            if (
                r_1 > r_2
                and distance_metric[combination[1] - 1][combination[0] - 1] == distance_metric[combination[0] - 1][combination[1] - 1]
                and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1] == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]
            ):
                list_pruning.append(combination[::-1])
            elif (
                r_2 > r_1
                and distance_metric[combination[1] - 1][combination[0] - 1] == distance_metric[combination[0] - 1][combination[1] - 1]
                and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1] == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]
            ):
                list_pruning.append(combination)
            elif r_2 == r_1 and distance_metric_rendezvouz[combination[1] - 1][combination[0] - 1] == distance_metric_rendezvouz[combination[0] - 1][combination[1] - 1]:
                list_pruning.append(combination[::-1])

        # Update available combinations after pruning
        self.available_combinations = [x for x in permutations_units if x not in list_pruning]

    def team_reward_generation(self):
        """
        This function will generate the rewards of the team of the graphical model
        """

        # Initialize lists to store team weights and types per department, rewards for each team, and other variables
        teams_weights_per_department = []
        teams_types_per_department = []
        reward_team = []
        total_team_array = []

        # Loop through each department to initialize team weight and type lists
        for state_index in range(0, len(self.max_state)):
            teams_weights_per_department.append([])
            teams_types_per_department.append([])

        sum_teams = []  # Initialize a list to store sum of teams per combination

        max_state = self.max_state  # Get maximum states for each department

        indices_departments = np.array(range(1, len(max_state) + 1))  # Get indices of departments

        # Iterate over available combinations of departments
        for combination in self.available_combinations:
            u = combination[0]  # Get first department index in the combination
            v = combination[1]  # Get second department index in the combination
            team_array = []
            # Iterate over team types
            for type in self.teams_types:
                # Perform calculations based on team type
                if type == 1:
                    # Calculations for type 1 teams
                    if u == v:
                        team_size = max_state[u - 1]
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1
                    else:
                        team_size = 0
                elif type == 2:
                    # Calculations for type 2 teams
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 2)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = min(max_state[u - 1], max_state[v - 1])
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] for i in range(0, team_size - 1)])
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend([weasels_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        teams_types_per_department[v - 1].append(type)
                        teams_types_per_department[v - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                elif type == 3:
                    # Calculations for type 3 teams
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 2)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1
                    else:
                        team_size = 0
                elif type == 4:
                    # Calculations for type 4 teams
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 3)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = math.floor(min(max_state[u - 1] / 2, max_state[v - 1]))
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] for i in range(0, team_size - 1)])
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend([weasels_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        teams_types_per_department[v - 1].append(type)
                        teams_types_per_department[v - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, [u, v])
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                elif type == 5:
                    # Calculations for type 5 teams
                    if u == v:
                        team_size = math.floor(max_state[u - 1] / 4)
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1] + weasels_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] + weasels_count[type - 1] for i in range(0, team_size - 1)])
                        teams_types_per_department[u - 1].append(type)
                        teams_types_per_department[u - 1].extend([type for i in range(0, team_size - 1)])
                        remaining_departments = np.setdiff1d(indices_departments, v)
                        reward = distance_metric[v - 1][u - 1] + F_u[u - 1] - K_1 + W_u[v - 1]
                    else:
                        team_size = math.floor(min(max_state[u - 1] / 2, max_state[v - 1] / 2))
                        teams_weights_per_department[u - 1].append(fighters_count[type - 1])
                        teams_weights_per_department[u - 1].extend([fighters_count[type - 1] for i in range(0, team_size - 1)])
                        teams_weights_per_department[v - 1].append(weasels_count[type - 1])
                        teams_weights_per_department[v - 1].extend([weasels_count[type - 1] for i in range(0, team_size - 1)])
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
        # Iterate over sum_teams to calculate last indices for each combination
        for i in range(0, len(sum_teams)):
            if i == 0:
                last_index.append(sum_teams[i] - 1)
            else:
                last_index.append(sum(sum_teams[0:i]) + sum_teams[i] - 1)
        # Return tuple containing calculated values
        return (
            reward_team,
            teams_weights_per_department,
            teams_types_per_department,
            last_index,
        )

    def process_inputs_ctypes(self):
        """
        This function will process the inputs that the C++ code will take

        Returns:
            transitional_model_p: a pointer of messages from teams to departments that will be udpated in C++
            max_state_p: a pointer that stores the capacities of departments
            reward_team_p: a pointer that stores the rewards of the teams
            reward_project_flatten_p: a pointer that stores the flattened rewards of the projects-teams
            teams_weights_per_department_flatten_p: a pointer that stores the flattened teams weights per departments (how much a team takes from departments)
            extrinsic_output_team_p: a pointer that will store the extrinsic output messages from teams after running IMPA in C++
            intrisic_out_mwm_p: a pointer that will store the intrinsic output messages from projects after running IMPA in C++
            non_zero_weight_indices_p: a pointer that stores the indices of non-zero weights or connections between teams and departments
            non_zero_weight_indices_sizes_p: a pointer that stores the sizes of non-zero weight indices for each department
            max_size_nonzero_weights: a pointer that stores the maximum size of non-zero weight indices across all departments
        """
        # Extract required data from object attributes
        teams_weights_per_department = self.teams_weights_per_department
        reward_project = self.reward_project
        reward_team = self.reward_team
        max_state = self.max_state
        num_departments = self.num_departments
        num_teams = self.num_teams
        num_projects = self.num_projects

        # Initialize arrays to store messages, rewards, etc.
        message_team_to_department = np.zeros((num_departments, num_teams), dtype=np_impa_lib)

        intrinsic_out_mwm = np.zeros((num_projects, num_teams), dtype=np_impa_lib)

        # Compute messages from teams to departments
        for i in range(0, num_departments):
            for j in range(0, num_teams):
                if self.teams_weights_per_department[i][j] != 0:
                    message_team_to_department[i][j] = reward_team[j]

        # Compute non-zero weight indices
        non_zero_weight_indices = []

        for department_index in range(0, num_departments):
            indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
            non_zero_weight_indices = non_zero_weight_indices + [indices]

        # Flatten and convert arrays to ctypes pointers
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
        non_zero_weight_indices_arr = np.zeros((num_departments, max_size_nonzero_weights), dtype=np.int32)  # added dtype=np.int32

        for i in range(len(non_zero_weight_indices)):
            for j in range(len(non_zero_weight_indices[i])):
                non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

        non_zero_weight_indices = np.array(non_zero_weight_indices_arr).flatten().astype(np.int32)
        non_zero_weight_indices_p = non_zero_weight_indices.ctypes.data_as(c_int_p)
        non_zero_weight_indices_sizes = np.array(non_zero_weight_indices_sizes).flatten().astype(np.int32)
        non_zero_weight_indices_sizes_p = non_zero_weight_indices_sizes.ctypes.data_as(c_int_p)

        # Return the ctypes pointers as a tuple
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
        """
        This function will run IMPA on the graphical model of the Knapsack-MWM problem
        """
        # Prepare inputs for C++ code
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

        # Call the C++ function to execute the IMPA algorithm
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

        # Retrieve and process outputs
        extrinsic_output_team = list(extrinsic_output_team_p.__dict__.values())[0]
        intrinsic_out_mwm = np.reshape(
            list(intrisic_out_mwm_p.__dict__.values())[0],
            (self.num_projects, self.num_teams),
        )
        # Update intrinsic output
        self.intrinsic_out_mwm = intrinsic_out_mwm

        self.intrinsic_output = extrinsic_output_team + self.reward_team

    def pre_analysis(self):
        """
        This function will run a pre-analysis on the C++ IMPA solution
        """
        # Retrieve necessary data
        intrinsic_output = self.intrinsic_output
        threshold = self.threshold

        # Check for valid match and calculate IMPA metric
        valid_match, impa_metric, p_mwm = self.check_match()

        # Identify selected teams indices
        self.matching_array = np.where(p_mwm == 1)
        selected_teams_indices = np.where(p_mwm == 1)[1]
        self.selected_teams_indices = np.unique(selected_teams_indices)

        # Apply threshold to obtain hard decisions
        hard_decision = np.array(deepcopy(intrinsic_output))
        hard_decision[intrinsic_output > threshold] = 0

        # Perform pre-processing analysis
        self.pre_processing_analysis()

    def check_match(self):
        """
        This function will check the assignment of the MWM problem

        Returns:
            valid_match: if constraints of MWM are true, then valid_match is true. Otherwise, false
            match_metric: this is the sum of intrinsic messages of the project equality constraint
            P: hard-decision on the MWM problem
        """

        # Retrieve intrinsic output matrix
        A = self.intrinsic_out_mwm

        # Retrieve unbalanced flag
        unbalanced_flag = self.unbalanced_flag

        # Calculate match metric
        match_metric = sum(x for x in A.flatten() if x < 0)
        P = deepcopy(A)
        threshold = self.threshold
        # Apply threshold to obtain hard decisions
        P[P > threshold] = 0
        P[P < threshold] = 1
        valid_match = False

        min_N = min(A.shape[0], A.shape[1])
        # Check for unbalanced flag
        if unbalanced_flag:
            P_interm = self.check_unbalanced_flag(P)
        else:
            P_interm = P
        # Calculate row and column sums
        row_sums = np.sum(P_interm, axis=1)
        col_sums = np.sum(P_interm, axis=0)

        # Check if constraints of the MWM problem are true
        if np.array_equal(col_sums, np.ones(min_N)):
            if np.array_equal(row_sums, np.ones(min_N)):
                valid_match = True
        return valid_match, match_metric, P

    def check_unbalanced_flag(self, P):
        """
        This function will find hard-decision if MWM imbalanced (number of teams is not equal to number of projects)

        Args:
            P: hard-decision on MWM
        Returns:
            P_interm: hard-decision on MWM when MWM is imbalanced (number of teams is not equal to number of projects)
        """
        idx = 0
        idy = 0
        P_interm = deepcopy(P)
        # Find columns with all zeros and delete them
        idy = np.argwhere(np.all(P_interm[..., :] == 0, axis=0))  # find which columns there are zeros
        P_interm = np.delete(P_interm, idy, axis=1)
        idx = np.argwhere(np.all(P_interm[:, ...] == 0, axis=1))  # find which rows there are zeros
        # Find rows with all zeros and delete them
        P_interm = np.delete(P_interm, idx, axis=0)
        return P_interm

    def pre_processing_analysis(self):
        """
        This function will run an analysis on the IMPA solution before post-processing
        """

        # Extract necessary attributes
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

        # Initialize flags
        capacity_exceeded_flag = False
        project_matched_flag = False
        ic_violated_flag = False

        # Create an array of all project indices
        all_projects = np.array(range(0, num_projects))

        # Initialize arrays
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
        # Loop through selected teams to analyze assignment
        for i in range(0, len(selected_teams_indices)):
            team_departments = []
            index = selected_teams_indices[i]
            project_index = np.where(matching_array[1] == index)[0]
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            # Store assigned departments and their weights
            for j in range(0, num_departments):
                weight = teams_weights_per_department[j][index]
                if weight != 0:
                    team_departments.append(j)
                    department_weights[j].append(weight)
                    team_departments_weights[j] = team_departments_weights[j] + weight
            # Compute package location
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
        # Check department capacities
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
        # Analyze project-team assignments
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

        # Check if all projects are matched to teams
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

        # Update class attributes with analysis results
        self.capacity_exceeded_flag = capacity_exceeded_flag
        self.ic_violated_flag = ic_violated_flag
        self.combination_reused_team = combination_reused_team
        self.team_departments_weights = team_departments_weights
        self.team_locations = team_locations
        self.reused_teams = reused_teams

    def post_analysis(self):
        """
        This function will run an analysis on the IMPA solution after post-processing
        """
        # Retrieve post-processing flags and options
        post_process_flag = self.post_process_flag

        # Execute post-processing based on the selected option
        if post_process_flag:
            if self.post_process_option == 1:  # brute force post-processing on bins
                self.apply_post_processing_departments_brute()
            elif self.post_process_option == 2:  # brute force post-processing on packages
                self.apply_post_processing_teams_brute()
                # exit()
        else:
            # If post-processing is not required, retain the intermediate results
            self.intermediate_selected_teams_indices = self.selected_teams_indices
            self.post_processing_ic_combination = []

        # Perform analysis on the post-processed results
        self.results_analysis()

    def apply_post_processing_departments_brute(self):
        """
        This function will run post-processing by investigating departments and removing least confident teams from these departments
        """
        # Retrieve necessary data
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
        # Handle cases where internal constraints are violated
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

        # Initialize intermediate variables for post-processing
        intermediate_team_department_weights = deepcopy(team_departments_weights)
        intermediate_selected_teams_indices = deepcopy(selected_teams_indices)

        # Perform post-processing if capacity constraints are violated
        if capacity_exceeded_flag:
            print("-------")
            print("POST-PROCESSING STARTED")
            print("-------")
            # Iterate until all departments satisfy capacity constraints
            while any(intermediate_team_department_weights > max_state):
                post_processing_departments = np.where(intermediate_team_department_weights > max_state)[0]
                print("Post-processing Departments: ", post_processing_departments)
                print("-------")
                if len(post_processing_departments):
                    # Find the least confident team in the department
                    intrinsic_output_selected = intrinsic_output[intermediate_selected_teams_indices]
                    processing_department_index = post_processing_departments[0]
                    print("Post-processing department: ", processing_department_index)
                    packages_processing_department_index = [
                        i
                        for i in intermediate_selected_teams_indices
                        if (processing_department_index in np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]]) - 1)
                    ]
                    print(
                        "max(reward_team[packages_processing_department_index]: ",
                        max(reward_team[packages_processing_department_index]),
                    )
                    removed_team = intermediate_selected_teams_indices[
                        np.where(intrinsic_output_selected == max(intrinsic_output_selected[intermediate_selected_teams_indices.searchsorted(packages_processing_department_index)]))[0][0]
                    ]
                    print("reward[removed_team]: ", reward_team[removed_team])
                    removed_team_location = team_locations[np.where(intermediate_selected_teams_indices == removed_team)[0][0]]
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
                    # Update weights and indices after removing the team
                    intermediate_selected_teams_indices = np.setdiff1d(intermediate_selected_teams_indices, removed_team)
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
        # Update class attributes with post-processing results
        self.intermediate_selected_teams_indices = intermediate_selected_teams_indices
        self.post_processing_ic_combination = post_processing_ic_combination

    def apply_post_processing_teams_brute(self):
        """
        This function will run post-processing by investigating all teams associated with a department with a violated capacity. Least confident ones belonging to multiple departments will be removed first
        """
        # Retrieve necessary data
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
        # Handle cases where inequality constraints are violated
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

        # Initialize intermediate variables for post-processing
        intermediate_team_department_weights = deepcopy(team_departments_weights)
        intermediate_selected_teams_indices = deepcopy(selected_teams_indices)

        # Perform post-processing if capacity constraints are violated
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
                        # Find teams associated with the current department and other departments
                        teams_processing_department_indices_combinations = [
                            [
                                i,
                                (np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]]) - 1).tolist(),
                            ]
                            for i in intermediate_selected_teams_indices
                            if (
                                (department in np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]]) - 1)
                                and available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]][0] - 1 in post_processing_departments
                                and available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]][1] - 1 in post_processing_departments
                            )
                        ]  # and
                        # len(set(np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices==i)[0][0]]])-1)) != 1)]
                        if not len(teams_processing_department_indices_combinations):
                            # If there are no teams associated with multiple departments, consider teams associated only with the current department
                            teams_processing_department_indices_combinations = [
                                [
                                    i,
                                    (np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]]) - 1).tolist(),
                                ]
                                for i in intermediate_selected_teams_indices
                                if (department in np.array(available_combinations[team_locations[np.where(intermediate_selected_teams_indices == i)[0][0]]]) - 1)
                            ]
                        # print('department: ',department, '\n teams_processing_department_indices_combinations: \n', teams_processing_department_indices_combinations)
                        # print(len(teams_processing_department_indices_combinations))
                        list_teams_processing_department_indices_combinations.extend(
                            y for y in teams_processing_department_indices_combinations if y not in list_teams_processing_department_indices_combinations
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
                        intrinsic_output_selected[intermediate_selected_teams_indices.searchsorted(np.array(team_indices_processing)[index_processing_combinations])],
                    )
                    # Remove the least confident team
                    removed_team = intermediate_selected_teams_indices[
                        np.where(
                            intrinsic_output_selected
                            == max(intrinsic_output_selected[intermediate_selected_teams_indices.searchsorted(np.array(team_indices_processing)[index_processing_combinations])])
                        )[0][0]
                    ]
                    print("removed_team:", removed_team)
                    # print('reward[removed_team]: ', reward_team[removed_team])
                    removed_team_location = team_locations[np.where(intermediate_selected_teams_indices == removed_team)[0][0]]
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
                    intermediate_selected_teams_indices = np.setdiff1d(intermediate_selected_teams_indices, removed_team)
                    team_locations.remove(removed_team_location)
                    u = removed_team_combination[0]
                    v = removed_team_combination[1]
                    # Update department weights
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

        # Store updated values
        self.intermediate_selected_teams_indices = intermediate_selected_teams_indices
        self.post_processing_ic_combination = post_processing_ic_combination

    def results_analysis(self):
        """
        This function will analyze the results of the IMPA solution
        """
        # Retrieve necessary data
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
        # Iterate through selected teams and analyze assignments
        for i in range(0, len(intermediate_selected_teams_indices)):
            index = intermediate_selected_teams_indices[i]
            # Determine the associated project index
            if index in reused_teams and post_process_flag:
                project_index = np.where(matching_array[0] == [i[0] for i in post_processing_ic_combination if i[1] == index][0])[0]
            else:
                project_index = np.where(matching_array[1] == index)[0]
            project_array.append(project_index)
            total_value = total_value + reward_team[index]
            team_type.append(teams_types[index])
            # Update department weights based on team assignments
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
        # Store results for comparison and print summary statistics
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

        # Calculate total weight of matched projects-teams
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

        # Display MWM matching results
        print("MWM")
        print("-------------")

        if len(project_array) == num_projects:
            project_matched_flag = True

        print(f"Matching All Projects to Teams : {project_matched_flag}")
        print("Total Weight RHS: ", total_weight)
        
        
class GraphicalModelKsat:
    def __init__(
        self,
        NUM_ITERATIONS,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        K_VARIABLE,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
        POST_PROCESS_FLAG,
        TYPE_METRICS,
        PP_ELEMENTS,
        IM_VARIANCE,
        OVERWRITE,
    ):
        self.num_iterations = NUM_ITERATIONS
        self.num_variables = NUM_VARIABLES
        self.num_constraints = NUM_CONSTRAINTS
        self.k_variable = K_VARIABLE
        self.threshold = THRESHOLD
        self.filtering_flag = FILTERING_FLAG
        self.alpha = ALPHA
        self.random_test_flag = RANDOM_TEST_FLAG
        self.type_metrics = TYPE_METRICS
        self.pp_elements = PP_ELEMENTS
        self.num_decisions = 2
        self.var = IM_VARIANCE
        self.overwrite = OVERWRITE
        self.post_process_flag = POST_PROCESS_FLAG
        
    def initialize(self):
        
        self.results_composed = []

        if not self.random_test_flag:
            self.num_variables = self.input_load[0]
            self.num_constraints = self.input_load[1]
            self.k_variable = self.input_load[2] 
            type_metrics = self.input_load[9]
            self.valid_sol = self.input_load[10]
        

        num_variables = self.num_variables
        num_constraints = self.num_constraints
        k_variable = self.k_variable

        print(f"num_variables: {num_variables}")
        print(f"num_constraints: {num_constraints}")
        print(f"k_variable: {k_variable}")
        if (not self.overwrite and not self.random_test_flag):
            print(f"type_metrics: {type_metrics}")
            self.type_metrics = type_metrics
        elif (not self.overwrite and self.random_test_flag):
            print(f"type_metrics: {self.type_metrics}")
        elif (self.overwrite and not self.random_test_flag):
            print(f"type_metrics changed from {type_metrics} to {self.type_metrics}")
        elif (self.overwrite and self.random_test_flag):
            print("Cannot Overwrite Metrics")
            print(f"type_metrics: {self.type_metrics}")
        
        type_metrics = self.type_metrics
        print(f"var: {self.var}")
        if (self.random_test_flag):
            constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, valid_sol = self.create_model_structure(num_variables, num_constraints, k_variable, type_metrics = type_metrics)
            self.valid_sol = valid_sol
            num_constraints = len(constraints_connections)
            if (num_constraints != self.num_constraints):
                print("Duplicate constraints were found. Pre-processing completed.")
                print(f"num_constraints: {num_constraints}")
            self.num_constraints  = num_constraints
        
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
            
            if (self.overwrite):
                incoming_metrics_cost = self.generate_incoming_metrics(self.valid_sol)
                self.incoming_metrics_cost = incoming_metrics_cost

        self.constraints_connections = constraints_connections
        self.constraints_connections_type = constraints_connections_type
        self.incoming_metrics_cost = incoming_metrics_cost
        self.variables_connections = variables_connections
        self.variables_connections_type = variables_connections_type
        valid_sol = self.valid_sol
        
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
        
        self.num_used_variables = len(self.used_variables)
        
        self.variables_connections_sizes = np.array([len(connection) for connection in self.variables_connections])    

    def generate_incoming_metrics(self, valid_sol):
        incoming_metrics_cost = np.zeros(self.num_variables)
        type_metrics = self.type_metrics
        normal_variance = self.var
        if (len(valid_sol)):
            for i, valid_sol_value in enumerate(valid_sol):
                if (type_metrics == 1): # correctly biased
                    mean = -normal_variance/2 if valid_sol_value == 1 else normal_variance/2
                    incoming_metrics_cost[i] = np.random.normal(loc=mean, scale=np.sqrt(normal_variance))
                elif (type_metrics == 0): # random IM
                    random_sign = np.random.choice([-1, 1])
                    mean = random_sign*normal_variance/2
                    incoming_metrics_cost[i] = np.random.normal(mean, np.sqrt(normal_variance))
                elif (type_metrics == 2): #normal(0,sigma^2)
                    incoming_metrics_cost[i] = np.random.normal(0, np.sqrt(normal_variance))
        else:
            if type_metrics == 1:
                raise ValueError("type_metrics cannot be equal to 1.")
            else:
                for i in range(self.num_variables):
                    if (type_metrics == 0): # random IM
                        random_sign = np.random.choice([-1, 1])
                        mean = random_sign*normal_variance/2
                        incoming_metrics_cost[i] = np.random.normal(mean, np.sqrt(normal_variance))
                    elif (type_metrics == 2): #normal(0,sigma^2)
                        incoming_metrics_cost[i] = np.random.normal(0, np.sqrt(normal_variance))   
        return incoming_metrics_cost
    
    def create_model_structure(self, num_variables, num_constraints, k_variable, type_metrics = 1):
        
        constraints_connections = [[] for _ in range(num_constraints)]
        constraints_connections_type = [[] for _ in range(num_constraints)]
        variables_connections = [[] for _ in range(num_variables)] 
        variables_connections_type = [[] for _ in range(num_variables)]
        used_variables = []
        
        valid_sol = np.random.randint(2, size=num_variables, dtype=int)
        incoming_metrics_cost = self.generate_incoming_metrics(valid_sol)     
        index_constraint = 0
        for constraint_idx in range(num_constraints):
            connections = []
            connections_types = []
            
            satisfying_variable_index = random.randint(0, num_variables - 1)
            satisfying_variable_value = valid_sol[satisfying_variable_index]
            connection_type = -1 if satisfying_variable_value == 0 else 1
            connections.append(satisfying_variable_index)
            connections_types.append(connection_type)
            used_variables.append(satisfying_variable_index)
            variables_connections[satisfying_variable_index].append(index_constraint)
            variables_connections_type[satisfying_variable_index].append(connection_type)
            
            for k_variable_idx in range(k_variable-1):
                choice_variable = random.choice([idx for idx in range(num_variables) if idx not in connections])
                connections.append(choice_variable)
                choice_connection_type = random.choice([-1, 1])
                connections_types.append(choice_connection_type)
                
                variables_connections[choice_variable].append(index_constraint)
                used_variables.append(choice_variable)
                variables_connections_type[choice_variable].append(choice_connection_type)
 
            found = False
            for i, sublist_i in enumerate(constraints_connections):
                if(set(sublist_i) == set(connections)):
                    found = True
                    
            if (not found):
                constraints_connections[index_constraint] = connections
                constraints_connections_type[index_constraint] = connections_types 
                index_constraint+=1
            else:
                for variable in connections:
                    index = variables_connections[variable].index(index_constraint)
                    variables_connections[variable].pop(index)
                    variables_connections_type[variable].pop(index)
                
        used_variables = list(set(used_variables))
        constraints_connections = [sublist for sublist in constraints_connections if sublist]
        constraints_connections_type = [sublist for sublist in constraints_connections_type if sublist]
        return constraints_connections, constraints_connections_type, incoming_metrics_cost, used_variables, variables_connections, variables_connections_type, valid_sol
    
    def run_impa(self):
        
        print("--------------")
        print("Running IMPA")
        print("--------------")
        
        (
            used_variables_flatten_p,
            variables_connections_flatten_p,
            variables_connections_sizes_flatten_p,
            constraints_connections_flatten_p,
            constraints_connections_type_flatten_p,
            incoming_metrics_cost_flatten_p,
            variable_ec_to_ksat_constraint_m_flatten_p,
            extrinsic_output_variable_ec_p,
        ) = self.process_inputs_ctypes()
        
        WrapperKsat(
            np.int32(self.num_iterations),
            np.int32(self.num_variables),
            np.int32(self.num_constraints),
            np.int32(self.num_used_variables),
            np_impa_lib(self.alpha),
            self.filtering_flag,
            used_variables_flatten_p,
            variables_connections_flatten_p,
            variables_connections_sizes_flatten_p,
            constraints_connections_flatten_p,
            constraints_connections_type_flatten_p,
            incoming_metrics_cost_flatten_p,
            variable_ec_to_ksat_constraint_m_flatten_p,
            extrinsic_output_variable_ec_p,
            np.int32(self.k_variable)
        )
        
        extrinsic_output_variable_ec = list(extrinsic_output_variable_ec_p.__dict__.values())[0]
    
        used_incoming_metrics_cost = np.zeros(self.num_variables)
        used_incoming_metrics_cost[self.used_variables] = self.incoming_metrics_cost[self.used_variables]
        intrinsic_out_variable_ec = extrinsic_output_variable_ec + used_incoming_metrics_cost
        self.intrinsic_out_variable_ec = intrinsic_out_variable_ec
        self.hard_decision_analysis()
        
        self.runtime_impa = time.time() - self.start_time
        print(f"IMPA Time: {self.runtime_impa}")
        print("--------")
        
        self.required_PP = False
        self.start_time_pp = time.time()
        self.runtime_pp = 0.0
        if (len(self.unsatisfied_constraints) and self.post_process_flag):
            self.required_PP = True
            self.run_post_processing()
            self.runtime_pp = time.time() - self.start_time_pp 
            print(f"PP Time: {self.runtime_pp}")
        
        print('--------')
        print("IMPA results:")
        
        if (not len(self.unsatisfied_constraints)):
            print("Problem satisfied")  
        else:
            print("Problem not satisfied") 
        
        
        #print("Active variables: ",self.active_variables,)
        #print("Inactive variables: ", self.inactive_variables,)
        
        if (len(self.valid_sol) and not len(self.unsatisfied_constraints)):
            self.valid_sol_active_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if self.valid_sol[self.used_variables][i]]
            self.valid_sol_inactive_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if not self.valid_sol[self.used_variables][i]]
            
            similarity_active = self.calculate_similarity(set(self.active_variables), set(self.valid_sol_active_variables))
            similarity_inactive = self.calculate_similarity(set(self.inactive_variables), set(self.valid_sol_inactive_variables))
            
            self.average_similarity = 100*(similarity_active+similarity_inactive)/2
            print(f"average_similarity: {self.average_similarity:.2f}")
        else:
            self.average_similarity = []
            
        if self.save_flag:
            self.results_composed.append(
                (
                    self.filtering_flag,
                    self.alpha,
                    self.type_metrics,
                    self.num_variables,
                    self.num_constraints,
                    self.k_variable,
                    self.constraints_connections,
                    self.constraints_connections_type,
                    self.incoming_metrics_cost,
                    self.used_variables,
                    self.variables_connections,
                    self.variables_connections_type,
                    self.valid_sol,
                    self.intrinsic_out_variable_ec,
                    self.active_variables,
                    self.inactive_variables,
                    self.satisfied_constraints,
                    self.unsatisfied_constraints,
                    self.formula_satisfied,
                    self.average_similarity,
                    self.required_PP,
                    self.runtime_impa,
                    self.runtime_pp
                )
            )

    def calculate_similarity(self,set_1, set_2):
        cardinality_intersection = len(set_1.intersection(set_2))
        cardinality_union_set = len(set_1.union(set_2))
        if (cardinality_union_set == 0):
            similarity = 1
        else:
            similarity = cardinality_intersection/cardinality_union_set
        return similarity
        
    def process_inputs_ctypes(self):
        
        used_variables_flatten = np.array(self.used_variables).flatten().astype(np.int32)
        used_variables_flatten_p = used_variables_flatten.ctypes.data_as(c_int_p)
        
        variables_connections_flatten = np.array([item for sublist in self.variables_connections for item in sublist]).flatten().astype(np.int32)
        variables_connections_flatten_p = variables_connections_flatten.ctypes.data_as(c_int_p)
        
        variables_connections_sizes_flatten = self.variables_connections_sizes.flatten().astype(np.int32)
        variables_connections_sizes_flatten_p = variables_connections_sizes_flatten.ctypes.data_as(c_int_p)
        
        constraints_connections_flatten = np.array(self.constraints_connections).flatten().astype(np.int32)
        constraints_connections_flatten_p = constraints_connections_flatten.ctypes.data_as(c_int_p)
        
        constraints_connections_type_flatten = np.array(self.constraints_connections_type).flatten().astype(np.int32)
        constraints_connections_type_flatten_p = constraints_connections_type_flatten.ctypes.data_as(c_int_p)
        
        incoming_metrics_cost_flatten = self.incoming_metrics_cost.flatten().astype(np_impa_lib)
        incoming_metrics_cost_flatten_p = incoming_metrics_cost_flatten.ctypes.data_as(c_impa_lib_type_p)
        
        variable_ec_to_ksat_constraint_m_flatten = self.variable_ec_to_ksat_constraint_m.flatten().astype(np_impa_lib)
        variable_ec_to_ksat_constraint_m_flatten_p = variable_ec_to_ksat_constraint_m_flatten.ctypes.data_as(c_impa_lib_type_p)

        extrinsic_output_variable_ec = np.zeros((self.num_variables))
        extrinsic_output_variable_ec = np.array(extrinsic_output_variable_ec).flatten().astype(np_impa_lib)
        extrinsic_output_variable_ec_p = extrinsic_output_variable_ec.ctypes.data_as(c_impa_lib_type_p)
        
        return (
            used_variables_flatten_p,
            variables_connections_flatten_p,
            variables_connections_sizes_flatten_p,
            constraints_connections_flatten_p,
            constraints_connections_type_flatten_p,
            incoming_metrics_cost_flatten_p,
            variable_ec_to_ksat_constraint_m_flatten_p,
            extrinsic_output_variable_ec_p,
        )

    def get_summary(self, hard_decision):
        
        self.active_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if hard_decision[self.used_variables][i]]
        self.inactive_variables = [self.used_variables[i] for i in range(len(self.used_variables)) if not hard_decision[self.used_variables][i]]
        
        indices_satisfied_constraints_solid = [(i, j) for i, sublist in enumerate(self.variables_connections_type) for j, element in enumerate(sublist) if element == 1 and i in self.active_variables]
        satisfied_constraints_solid = list(set([self.variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_solid]))
        
        indices_satisfied_constraints_dashed = [(i, j) for i, sublist in enumerate(self.variables_connections_type) for j, element in enumerate(sublist) if element == -1 and i in self.inactive_variables]
        satisfied_constraints_dashed = list(set([self.variables_connections[index[0]][index[1]] for index in indices_satisfied_constraints_dashed])) 
        
        satisfied_constraints = list(set(satisfied_constraints_solid + satisfied_constraints_dashed))
        self.satisfied_constraints = satisfied_constraints
        
        unsatisfied_constraints = [constraint for constraint in range(self.num_constraints) if constraint not in satisfied_constraints]
        self.unsatisfied_constraints = unsatisfied_constraints
                  
    def print_summary(self, unsatisfied_constraints):
        
        if (len(unsatisfied_constraints)):
            self.formula_satisfied = False
            print(f"{len(unsatisfied_constraints)} unsatisfied constraints: {unsatisfied_constraints}")
        else:
            self.formula_satisfied = True
            print(f"All {self.num_constraints} constraints are satisfied")
            
    def hard_decision_analysis(self,):
        
        intrinsic_out_variable_ec = self.intrinsic_out_variable_ec

        hard_decision = np.array(deepcopy(intrinsic_out_variable_ec))
        hard_decision[intrinsic_out_variable_ec > self.threshold] = 0
        hard_decision[intrinsic_out_variable_ec <= self.threshold] = 1
        
        self.hard_decision = hard_decision
        
        self.get_summary(self.hard_decision)
        self.print_summary(self.unsatisfied_constraints)
        
        similar_elements = []
        for i, sublist_i in enumerate(self.constraints_connections):
            for j, sublist_j in enumerate(self.constraints_connections):
                if i<j and set(sublist_i) == set(sublist_j):
                    similar_elements.append((i, j))
        if (len(similar_elements)):
            print("Duplicate constraints found")  
        
    def run_post_processing(self):
    
        print('--------')
        print("Before Post-Processing")
        #print("Active variables: ",self.active_variables,)
        #print("Inactive variables: ", self.inactive_variables,)
        
        print("--------------")
        print("Post-Processing started")
        print("--------------")
        
        common_elements_list = []
        
        print("Checking common variables across violated constraints")
        for violated_constraint in self.unsatisfied_constraints:
            remaining_violated_constraints = [self.constraints_connections[i] for i in self.unsatisfied_constraints if i != violated_constraint]
            remaining_violated_constraints_type = [self.constraints_connections_type[i] for i in self.unsatisfied_constraints if i != violated_constraint]
            for index_remaining_constraint, remaining_constraint_connections in enumerate(remaining_violated_constraints):
                common_elements = [(element_1, self.constraints_connections_type[violated_constraint][index_1]) for index_1, element_1 in enumerate(self.constraints_connections[violated_constraint]) for index_2, element_2 in enumerate(remaining_constraint_connections) if (element_1 == element_2 and self.constraints_connections_type[violated_constraint][index_1] == remaining_violated_constraints_type[index_remaining_constraint][index_2] \
                                    and ((element_1 in self.active_variables and self.constraints_connections_type[violated_constraint][index_1] == -1) or (element_1 in self.inactive_variables and self.constraints_connections_type[violated_constraint][index_1] == 1)))]
                if (common_elements):
                    common_elements_with_count = [(element[0], element[1]) for element in common_elements]
                    common_elements_list.extend([element for element in common_elements_with_count if element not in common_elements_list])
        
        old_length_unsatisfied_constraints = len(self.unsatisfied_constraints)
        for variable, type in common_elements_list:
            hard_decision_new = copy.deepcopy(self.hard_decision)
            new_variable_hard_decision = 0 if type == -1 else 1
            hard_decision_new[variable] = new_variable_hard_decision
            self.get_summary(hard_decision_new)
            if (len(self.unsatisfied_constraints)) < old_length_unsatisfied_constraints:
                print(f"Variable {variable} changed to {new_variable_hard_decision}")
                self.hard_decision[variable] = new_variable_hard_decision
                old_length_unsatisfied_constraints = len(self.unsatisfied_constraints)
                self.print_summary(self.unsatisfied_constraints)
                if (not len(self.unsatisfied_constraints)):
                    print("Post-Processing Ended")
                    break
        self.get_summary(self.hard_decision)
        combinations_list = []
        print('--------')
        print("Checking common variables across each violated constraint and its neighboring constraints (only satisfied by nodes in a violated constraint)")
        
        for violated_constraint in self.unsatisfied_constraints:
            combinations_list.append(self.constraints_connections[violated_constraint])
            constraints_connections_investigated = [[connections, self.constraints_connections_type[index], index] for index, connections in enumerate(self.constraints_connections) if any(element in connections for element in self.constraints_connections[violated_constraint]) and not all(element in connections for element in self.constraints_connections[violated_constraint])]
            count_satisfied = np.zeros(len(constraints_connections_investigated))
            satisfying_variables_constraints = [[] for _ in constraints_connections_investigated]
            for i in range(len(constraints_connections_investigated)):
                for j in range(len(constraints_connections_investigated[i][0])):
                    if (constraints_connections_investigated[i][0][j] not in self.constraints_connections[violated_constraint] and constraints_connections_investigated[i][1][j] ==1 and constraints_connections_investigated[i][0][j] in self.active_variables):
                        count_satisfied[i] +=1
                        satisfying_variables_constraints[i].append(constraints_connections_investigated[i][0][j])
                    elif (constraints_connections_investigated[i][0][j] not in self.constraints_connections[violated_constraint] and constraints_connections_investigated[i][1][j] == -1 and constraints_connections_investigated[i][0][j] in self.inactive_variables):
                        count_satisfied[i] +=1
                        satisfying_variables_constraints[i].append(constraints_connections_investigated[i][0][j])
            
            index_constraints_of_interest = np.where(count_satisfied==0)[0]
            overlap_common_type_list = []
            overlap_forbidden_list = []
            
            num_decisions = self.num_decisions
            for i in range(len(index_constraints_of_interest)):
                index = index_constraints_of_interest[i]
                index_constraint = constraints_connections_investigated[index][2]
                overlap = [x for x in self.constraints_connections[index_constraint] if x in self.constraints_connections[violated_constraint]]
                overlap_common_type = [x for x in overlap if self.constraints_connections_type[index_constraint][self.constraints_connections[index_constraint].index(x)] == self.constraints_connections_type[violated_constraint][self.constraints_connections[violated_constraint].index(x)]]
                if (overlap_common_type):
                    overlap_common_type_list.append(overlap_common_type)
                if (len(overlap)==1):
                    overlap_forbidden_list.append(overlap)
                combinations_list.append(self.constraints_connections[index_constraint])
            overlap_common_type_set = {item for sublist in overlap_common_type_list for item in sublist}
            overlap_forbidden_set = {item for sublist in overlap_forbidden_list for item in sublist}
            pool_variables = overlap_common_type_set.difference(overlap_forbidden_set)      
              
            unsatisfied_constraints_old = copy.deepcopy(self.unsatisfied_constraints)
            old_length_unsatisfied_constraints = len(unsatisfied_constraints_old)
            if (len(pool_variables)):
                for variable_in_pool in pool_variables:
                    hard_decision_new = copy.deepcopy(self.hard_decision)
                    new_variable_hard_decision = 1 if self.constraints_connections_type[violated_constraint][self.constraints_connections[violated_constraint].index(variable_in_pool)] ==1 else 0
                    hard_decision_new[variable_in_pool] = new_variable_hard_decision
                    self.get_summary(hard_decision_new)
                    if (len(self.unsatisfied_constraints)) < old_length_unsatisfied_constraints:
                        print(f"Variable {variable_in_pool} changed to {new_variable_hard_decision}")
                        self.hard_decision[variable_in_pool] = new_variable_hard_decision
                        old_length_unsatisfied_constraints = len(self.unsatisfied_constraints)
                        self.print_summary(self.unsatisfied_constraints)
                        if (not len(self.unsatisfied_constraints)):
                            print("Post-Processing Ended")
                            break
        
        print('--------')
        self.get_summary(self.hard_decision)
        print("Brute force search across variables in violated constraints and their neighboring constraints (only satisfied by nodes in a violated constraint)")
        if (len(self.unsatisfied_constraints)):
            unsatisfied_constraints_old = copy.deepcopy(self.unsatisfied_constraints)
            old_length_unsatisfied_constraints = len(unsatisfied_constraints_old)
            for num_elements in range(self.pp_elements, 1, -1):
                print(f"Investigating combinations of size {num_elements}")
                new_combinations_list = list(set(element for sublist in combinations_list for element in sublist))
                new_combinations_list = list(itertools.combinations(new_combinations_list, num_elements)) 
                print(f"Number of combinations: {len(new_combinations_list)}")
                for index_element, element in enumerate(new_combinations_list):
                    possible_element_decisions = list(itertools.product(range(num_decisions), repeat=num_elements))
                    for element_decision in possible_element_decisions:
                        hard_decision_new = copy.deepcopy(self.hard_decision)
                        for i, sub_element_decision in enumerate(element_decision):
                            hard_decision_new[element[i]] = sub_element_decision   
                        self.get_summary(hard_decision_new)
                        if (len(self.unsatisfied_constraints)) < old_length_unsatisfied_constraints:
                            for i, sub_element_decision in enumerate(element_decision):
                                if (self.hard_decision[element[i]] != hard_decision_new[element[i]]):
                                    print(f"Variable {element[i]} changed to {int(hard_decision_new[element[i]])}")
                            old_length_unsatisfied_constraints = len(self.unsatisfied_constraints)
                            self.hard_decision = hard_decision_new
                            self.print_summary(self.unsatisfied_constraints)
                            if (not len(self.unsatisfied_constraints)):
                                print("Post-Processing Ended")
                            break  
                    if (not len(self.unsatisfied_constraints)):
                        break          
        self.get_summary(self.hard_decision)
             
    def save_outputs(self):
        """
        This function will save results_composed (a dictionary containing variables of interest)
        """
        # Saving the results to a file
        with open(f"{self.folder_outputs}/outputs_set{self.test_file}.pkl", "wb") as f:
            pkl.dump(self.results_composed, f)