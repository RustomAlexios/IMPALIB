# **IMPALIB**

=============

Table of contents
-----------------

* [Introduction](#introduction)
* [Application](#application)
* [Supported Constraints](#supported-constraints)
* [Code Parameters](#code-parameters)
* [Usage](#usage)
* [Requirements and Installation](#requirements-and-installation)
* [License](#license)

### **Introduction**

---

Would you like to solve optimization problems using message-passing algorithms? And would you like perform Belief Propagation in C++ and python? Then IMPALIB (Iterative Message Passing Algorithm Library) is exactly for you.

IMPALIB:

- can be used in three different ways
  1. A small header-only library written in modern and pure C++
  2. A pure python code
  3. A C++ code with a python wrapper
- is very easy to integrate and use
- is self contained and depends only [cnpy](https://github.com/rogersce/cnpy) library for unit testing - also header-only library which is used for reading and writing numpy files
- supports inference for optimization problems with various possible constraints

### **Application**

---

- Application 1: We considered a first application of team formation and assignment problem based on a matrix management structure. We assumed that there are $N_d$ departments, with department $d$ comprising $S_d$ employees of a given specialty (e.g., hardware, systems, software, etc.). A total of $N_t$ teams are considered for formation by drawing from these departments. For example, a given team may require $2$ software engineers, $1$ hardware engineer and no systems engineers. The binary variable $y_t \in {0, 1}$ indicates whether team $t$ is formed. Each of the formed teams can be assigned to one of $N_p$ projects. The cost of forming team $t$ is $c_T(t)$ and the cost of matching team $t$ to project $p$ is $c_P(t, p)$. These costs will vary based on the team composition and how well a team matches the project needs. The problem of optimal team formation is a multiple knapsack problem and team-project assignment is a weighted matching problem. We focused on the joint optimization of these two problems.
- Application 2: We also considered solving yhe Traveling Salesperson Application (TSP) to determine the shortest possible route passing through a set of $N$ cities exactly once and returning to the starting point. In this case, cities are represented as nodes and their connections as edges within a graph. The challenge lies in finding an optimal Hamiltonian cycle in this graph while considering the distances between cities.

### **Supported Constraints**

---

Various constraints are implemented:

1. Equality Constraints which allow multiples copies of a variable
2. Inequality Constraints ($\le 1$):
   - The all-zero configuration has a zero configuration metric
   - The other configurations have exactly one edge variable taking the $1$ value
3. XNIC Constraints:
   - These constraints do not allow matching to a project if a team is not formed
   - These constraints only allow matching to exactly one project if a team is formed
4. Knapsack Contsraints:
   - These constraints enforce that the team assignments do not violate the capacity of each department
5. Degree Constraints which enforce that each city must be visited exactly once
6. Subtour Elimination constraints prevent the existence of smaller loops or subtours within potential solutions.

#### **Code Parameters**

---

1. Common variables to Problems $1$ and $2$:
   - nITER: Number of iterations of the IMPA
   - filteringFlag: whether filtering is activated or not
   - alpha: filtering coefficient ($0\le \alpha \lt 1$)
   - PPFlag: whether post-processing is activated or not
   - threshold: threshold value to make hard-decisions
2. Variables related to Application $1$:
   - PPOption: type of post-processing ($1$-perform post-processing on departments, $2$-perform post-processing on teams)
3. Variables related to Application $2$:
   - symFlag: whether to use symmetrical graphical model or not
   - augmFlag: whether to solve the augmented TSP or relaxed one
   - lkhSolFlag: Run LKH algorithm to find a tour
   - resetFlag: whether to reset messages at each augmentation step or not
   - randomTestFlag: whether to run IMPA on a graphical model with random cost values or from a pre-generated input files
   - inputPath: If randomTestFlag is not set to True, specify the inputs path
   - maxCount: Maximum count of invalid configurations before stopping the IMPA. This will allow IMPA to exit when singular solutions are encountered consecutively
   - maxAugmCount: Maximum number of augmentation steps. This will allow IMPA to exit when no tour was found up to a certain number of augmentation steps
   - KOPTFlag: whether to perform K-opt algorithm on the obtained tour (only $2$-opt and $3$-opt are investigated to avoid the curse of dimensionality)

### **Usage**

---

There are three different ways for implementing IMPALIB:

1. using a header-only C++ library
2. using a pure python code which is relatively slow
3. using a C++ code with a python wrapper which is relatively fast

##### **1. Header-only C++ library**

The headers in the `include` directory can be directly copied to your project
We assume in the code samples below you've copied them to an `impalib` subdirectory of one of your project's include directories

- Include header files of the library:

```cpp
    #include "impalib/impalib.hpp"
```

- Demo code of Application $1$:

  ```cpp
  // Copyright 2023, Alexios Rustom.
  // https://github.com/RustomAlexios/IMPALIB
  // Distributed under the MIT License.
  // (See accompanying LICENSE file or at
  //  https://opensource.org/licenses/MIT)

  #include "impalib/impalib.hpp"

  int main()
  {
      // Problem parameters
      const int                 N_PROJECTS                    = 2;      ///< number of projects
      const bool                FILT_FLAG                     = true;   ///< whether filtering is activated or not
      const int                 N_ITER                        = 400;    ///< number of iterations of IMPA
      const int                 N_DEPARTMENTS                 = 2;      ///< number of departments
      array<int, N_DEPARTMENTS> max_state                     = {3, 3}; ///< size of each department
      const impalib_type        ALPHA                         = 0.9;    ///< filtering parameter
      const int                 N_TEAMS                       = 5;      ///< number of teams
      array<int, N_DEPARTMENTS> non_zero_weight_indices_sizes = {4, 4}; ///< sizes of non-zero connections per department

      // Extracting data pointers
      const int *pNON_ZERO_WEIGHT_INDICES_SIZES = non_zero_weight_indices_sizes.data();

      // Calculate max size of non-zero weight
      const auto max_size_non_zero_weight_iter =
          max_element(non_zero_weight_indices_sizes.begin(), non_zero_weight_indices_sizes.end());
      int max_size_non_zero_weight = *max_size_non_zero_weight_iter;

      // Initialize the graphical model
      GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG,
                                      ALPHA);

      // Define rewards for teams and projects
      array<impalib_type, N_TEAMS>             reward_team    = {-22, -31, -68, -39, -84};
      array<impalib_type, N_PROJECTS *N_TEAMS> reward_project = {
          44, 1,  41, 10, 3,
          7,  17, 56, 98, 63};

      // Extracting data pointers
      const impalib_type *pREWARD_PROJECT = reward_project.data();
      const impalib_type *pREWARD_TEAM    = reward_team.data();

      // Defining and extracting variables
      array<impalib_type, N_DEPARTMENTS *N_TEAMS> transition_model  = {-22, -31, -68, -39, 0, 0, -31, -68, -39, -84};
      impalib_type                               *pTransition_model = transition_model.data();
      
      array<int, N_DEPARTMENTS *N_TEAMS> teams_weights_per_department = {2, 1, 1, 1, 0, 0, 1, 1, 1, 2};

      const int *pTEAMS_WEIGHTS_PER_DEPARTMENT = teams_weights_per_department.data();

      array<int, N_DEPARTMENTS *N_TEAMS> non_zero_weight_indices   = {0, 1, 2, 3, 1, 2, 3, 4};
      const int                         *p_NON_ZERO_WEIGHT_INDICES = non_zero_weight_indices.data();

      const int *pMAX_STATE = max_state.data();

      // Initialize and iterate through the model
      model_graph.initialize(pREWARD_TEAM, pTransition_model, pTEAMS_WEIGHTS_PER_DEPARTMENT,
                            pNON_ZERO_WEIGHT_INDICES_SIZES, p_NON_ZERO_WEIGHT_INDICES, pREWARD_PROJECT, pMAX_STATE);
      model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES);
      
      // Print results
      for (int i = 0; i < N_TEAMS; i++)
      {
          cout << model_graph.outputs.ExtrinsicOutputTeam[i] + model_graph.modelInputs_.RewardTeam[i] << '\n';
      }
  }
  ```

<!--Graphical Model of Application $1$:

- reward_team and reward_project are represented by arrows on the left and right equality constraints, respectively.
- teams_weights_per_department are represented by red (weight $=2$) or blue (weight $=1$) edges.

![graphicalModel](./img/demoGraphicalModel.png)-->

- Demo code of Application $2$:

  ```cpp
  // Copyright 2023, Alexios Rustom.
  // https://github.com/RustomAlexios/IMPALIB
  // Distributed under the MIT License.
  // (See accompanying LICENSE file or at
  //  https://opensource.org/licenses/MIT)

  #include "impalib/impalib.hpp"

  int main()
  {
      const int  N_NODES          = 5;    ///< number of nodes
      const bool FILT_FLAG        = true; ///< whether filtering is activated or not
      const int  N_ITER           = 200;  ///< number of iterations of IMPA
      const int  N_EDGE_VARIABLES = N_NODES * N_NODES - N_NODES;  ///< number of edge variables
      const impalib_type ALPHA             = 0.5;     ///< filtering parameter
      const bool         SYM_FLAG          = true;    ///< symmetric flag
      const bool         RESET_FLAG        = false;   ///< reset flag
      const impalib_type THRESHOLD         = -0.0001; ///< threshold parameter
      const bool         AUGM_FLAG         = true;    ///< augmentation flag
      const int          MAX_AUGM_COUNT    = 50;      ///< maximum augmentation count
      const int          MAX_FAILURE_COUNT = 50;      ///< maximum count for failure

      // connections between nodes for each edge variable
      array<array<int, 2>, N_EDGE_VARIABLES> edge_connections = {{{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 0}, {1, 2}, {1, 3},
                                                                  {1, 4}, {2, 0}, {2, 1}, {2, 3}, {2, 4}, {3, 0}, {3, 1},
                                                                  {3, 2}, {3, 4}, {4, 0}, {4, 1}, {4, 2}, {4, 3}}};

      // cost matrix of the tsp problem
      array<array<impalib_type, N_NODES>, N_NODES> cost_matrix = {
          {{0.0, 151.75578773, 56.18610887, 718.31915651, 293.02503715},
          {151.75578773, 0.0, 568.4231286, 83.25740946, 545.45357536},
          {56.18610887, 568.4231286, 0.0, 445.55107005, 888.09445172},
          {718.31915651, 83.25740946, 445.55107005, 0.0, 719.05730714},
          {293.02503715, 545.45357536, 888.09445172, 719.05730714, 0.0}}};

      // cost_edge_variable constructed from edge_connections and cost_matrix
      array<impalib_type, N_EDGE_VARIABLES> cost_edge_variable = {
          151.75578773, 56.18610887,  718.31915651, 293.02503715, 151.75578773, 568.4231286,  83.25740946,
          545.45357536, 56.18610887,  568.4231286,  445.55107005, 888.09445172, 718.31915651, 83.25740946,
          445.55107005, 719.05730714, 293.02503715, 545.45357536, 888.09445172, 719.05730714};

      // edge_degree_constraint_cost constructed from edge_connections and cost_matrix
      array<array<impalib_type, N_NODES>, N_EDGE_VARIABLES> edge_degree_constraint_cost = {
          {{151.75578773, 151.75578773, 0.0, 0.0, 0.0}, {56.18610887, 0.0, 56.18610887, 0.0, 0.0},
          {718.31915651, 0.0, 0.0, 718.31915651, 0.0}, {293.02503715, 0.0, 0.0, 0.0, 293.02503715},
          {151.75578773, 151.75578773, 0.0, 0.0, 0.0}, {0.0, 568.4231286, 568.4231286, 0.0, 0.0},
          {0.0, 83.25740946, 0.0, 83.25740946, 0.0},   {0.0, 545.45357536, 0.0, 0.0, 545.45357536},
          {56.18610887, 0.0, 56.18610887, 0.0, 0.0},   {0.0, 568.4231286, 568.4231286, 0.0, 0.0},
          {0.0, 0.0, 445.55107005, 445.55107005, 0.0}, {0.0, 0.0, 888.09445172, 0.0, 888.09445172},
          {718.31915651, 0.0, 0.0, 718.31915651, 0.0}, {0.0, 83.25740946, 0.0, 83.25740946, 0.0},
          {0.0, 0.0, 445.55107005, 445.55107005, 0.0}, {0.0, 0.0, 0.0, 719.05730714, 719.05730714},
          {293.02503715, 0.0, 0.0, 0.0, 293.02503715}, {0.0, 545.45357536, 0.0, 0.0, 545.45357536},
          {0.0, 0.0, 888.09445172, 0.0, 888.09445172}, {0.0, 0.0, 0.0, 719.05730714, 719.05730714}}};

      // Extract data pointers from various structures
      const int          *pEDGE_CONNECTIONS_PY               = addressof(get<0>(edge_connections[0]));
      const impalib_type *pCOST_MATRIX_PY                    = addressof(get<0>(cost_matrix[0]));
      const impalib_type *pCOST_EDGE_VARIABLE_PY             = cost_edge_variable.data();
      impalib_type       *pEdge_ec_to_degree_constraint_m_py = addressof(get<0>(edge_degree_constraint_cost[0]));
      const impalib_type *pEDGE_DEGREE_CONSTRAINT_COST_PY    = addressof(get<0>(edge_degree_constraint_cost[0]));

      GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGM_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD,
                                    MAX_FAILURE_COUNT);

      // Initialize the model with provided data pointers
      model_graph.initialize(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY,
                            pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);

      // Iterate through the relaxed graph
      model_graph.iterate_relaxed_graph();

      // Check if subtour constraints are not satisfied and augmentation is enabled
      if (!model_graph.subtourConstraintsSatisfiedFlag && AUGM_FLAG)
      {   
          // Perform augmentation if conditions are met
          model_graph.perform_augmentation(MAX_AUGM_COUNT);

      }
  }
  ```
- To run any of the above codes:

  - navigate to: ``IMPALIB/examples/KcMwm`` or ``IMPALIB/examples/Tsp``
  - Run: ``cmake -B build ``
  - Run: ``cmake --build build ``
  - Run: ``cd build``
  - Run: ``./demo ``

##### **2. Pure Python code**

To run pure code using sample datasets:

- Application $1$:

  - Navigate to ``IMPALIB/test/python_kc_mwm/src``
  - Run:
    ```bash
        python3 main_pure_optimized.py --nITER=400 --filteringFlag=True --alpha=0.9 --PPFlag=True --threshold=-0.0001
    ```
- Application $2$:

  - Navigate to ``IMPALIB/test/python_tsp/src``
  - Run:
    ```bash
        python3 main_tsp.py --nNodes=10 --filteringFlag=True --alpha=0.5 --augmFlag=True --threshold=-0.0001 --nITER=200 --randomTestFlag=True
    ```

##### **3. C++ code with a python wrapper**

To compile the C++ library and install the Python wrapper, navigate to the project root and use:

```bash
    python3 -m pip install . -v
```

- Navigate to ``IMPALIB/src/impa``
- For Application $1$:

  - To run wrapper code using sample datasets, Run:

    ```bash
        python3 main_kc_mwm.py --nITER=400 --filteringFlag=True --alpha=0.9 --PPFlag=True --PPOption=1 --threshold=-0.0001 
    ```
- For Application $2$:

  - To run wrapper code using sample datasets:
    ```bash
        python3 main_tsp.py --filteringFlag=True --alpha=0.5 --augmFlag=True --threshold=-0.0001 --nITER=200 --inputPath=inputs_1000_nNodes_random --testFile=0 --lkhSolFlag=True --maxAugmCount=10
    ```
  - To run wrapper code using randomly generated cost values:
    ```bash
        python3 main_tsp.py --nNodes=10 --filteringFlag=True --alpha=0.5 --augmFlag=True --threshold=-0.0001 --nITER=200 --randomTestFlag=True --lkhSolFlag=True --maxAugmCount=20
    ```

**Note**: currently this option looks for a relevant sample dataset in the `data` directory, one directory up from the current working directory.
This will be fixed in a future version.

### **Requirements and Installation**

- A C++ $11$ -compatible compiler
- Python $3.9.7$
- To perform unit testing: randomized simulations using a pure python code and a python wrapper around a C++ code are carried out for both applications. A checking routine on the stored numpy files is executed to compare results. An external library called [cnpy](https://github.com/rogersce/cnpy) is used to save and load numpy arrays in C++

### **Unit Testing**

- Refer to this README file for [Unit Testing](test/README.md) framework.

### **License**

Distributed under the MIT License.
See accompanying file [`LICENSE`](https://github.com/RustomAlexios/IMPALIB/blob/main/LICENSE) or at
([https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT))
