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
------------

Would you like to solve optimization problems using message-passing algorithms? And would you like perform Belief Propagation in C++ and python? Then IMPALIB (Iterative Message Passing Algorithm Library) is exactly for you.

IMPALIB:
- can be used in three different ways
    1. A small header-only library written in modern and pure C++
    2. A pure python code
    3. A C++ code with a python wrapper
- is very easy to integrate and use
- is self contained and depends only [cnpy](https://github.com/rogersce/cnpy) library for unit testing - also header-only library which is used for reading and writing numpy files
- supports inference for optimization problems with various possible constraints
- results in a much faster solution than Google OR-Tools

### **Application**
------------

We considered a canonical example team formation and assignment problem based on a matrix management structure. We assumed that there are $N_d$ departments, with department $d$ comprising $S_d$ employees of a given specialty (e.g., hardware, systems, software, etc.). A total of $N_t$ teams are considered for formation by drawing from these departments. For example, a given team may require $2$ software engineers, $1$ hardware engineer and no systems engineers. The binary variable $y_t \in {0, 1}$ indicates whether team $t$ is formed. Each of the formed teams can be assigned to one of $N_p$ projects. The cost of forming team $t$ is $c_T(t)$ and the cost of matching team $t$ to project $p$ is $c_P(t, p)$. These costs will vary based on the team composition and how well a team matches the project needs. The problem of optimal team formation is a multiple knapsack problem and team-project assignment is a weighted matching problem. We focused on the joint optimization of these two problems.

### **Supported Constraints**
------------

Various constraints are implemented:
1.  Equality Constraints which allow multiples copies of a variable
2. Inequality Constraints ($\le 1$):
    - The all-zero configuration has a zero configuration metric
    - The other configurations have exactly one edge variable taking the $1$ value
3. XNIC Constraints:
    - These constraints do not allow matching to a project if a team is not formed
    - These constraints only allow matching to exactly one project if a team is formed
4. Knapsack Contsraints:
    - These constraints enforce that the team assignments do not violate the capacity of each department

#### **Code Parameters**
------------
- nITER: Number of iterations of the IMPA
- filterFlag: whether filtering is activated or not
- alpha: filtering coefficient ($0\le \alpha \lt 1$)
- PPFlag: whether post-processing is activated or not
- PPOption: type of post-processing ($1$-perform post-processing on departments, $2$-perform post-processing on teams)
- threshold: threshold value to make hard-decisions

### **Usage**
------------

There are three different ways for implementing IMPALIB:
1. using a header-only C++ library 
2. using a pure python code which is relatively slow
3. using a C++ code with a python wrapper which is relatively fast

##### **1. Header-only C++ library**
The headers in the `include` directory can be directly copied to your project.
We assume in the code samples below you've copied them to an `impalib` subdirectory of one of your project's include directories.

- Include header files of the library:
```cpp
    #include "impalib/impalib.hpp"
```

- Demo code:
```cpp
    // Copyright 2023, Alexios Rustom.
    // https://github.com/RustomAlexios/IMPALIB
    // Distributed under the MIT License.
    // (See accompanying LICENSE file or at
    //  https://opensource.org/licenses/MIT)

    #include "impalib/impalib.hpp"

    int main(){
        const int N_PROJECTS = 2; //number of projects
        const bool FILT_FLAG = true; //whether filtering is activated or not
        const int N_ITER = 400; //number of iterations of IMPA
        const int N_DEPARTMENTS = 2; //number of departments
        int max_state[N_DEPARTMENTS] = {3,3}; //size of each department
        const impalib_type ALPHA = 0.9; //filtering parameter
        const int N_TEAMS = 5; //number of teams
        int non_zero_weight_indices_sizes[N_DEPARTMENTS] = {4, 4}; 

        const int* pNON_ZERO_WEIGHT_INDICES_SIZES = non_zero_weight_indices_sizes;
        int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES , pNON_ZERO_WEIGHT_INDICES_SIZES + N_DEPARTMENTS);

        GraphicalModel model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

        impalib_type reward_team[N_TEAMS] = {-22, -31, -68, -39, -84}; //cost of activating a team
        impalib_type reward_project[N_PROJECTS*N_TEAMS] = {44, 1, 41, 10, 3,
                                                        7, 17, 56, 98, 63}; //cost of assigning a team to a project
        
        const impalib_type* pREWARD_PROJECT = reward_project;
        const impalib_type* pREWARD_TEAM = reward_team;

        impalib_type transition_model[N_DEPARTMENTS*N_TEAMS]= {-22, -31, -68, -39, 0,        
                                                        0, -31, -68, -39, -84};
        impalib_type* pTransition_model = transition_model;

        int teams_weights_per_department[N_DEPARTMENTS*N_TEAMS]= {2,1,1,1,0,
                                                                0,1,1,1,2}; //red edges correspond to weight=2
                                                                            //blue edges correspond to weight=1
        const int* pTEAMS_WEIGHTS_PER_DEPARTMENT = teams_weights_per_department;

        int non_zero_weight_indices[N_DEPARTMENTS*N_TEAMS] = {0,1,2,3,1,2,3,4}; 
        const int* p_NON_ZERO_WEIGHT_INDICES = non_zero_weight_indices;


        const int* pMAX_STATE = max_state;
        
        model_graph.initialize(pREWARD_TEAM, pTransition_model, pTEAMS_WEIGHTS_PER_DEPARTMENT,
                            pNON_ZERO_WEIGHT_INDICES_SIZES, p_NON_ZERO_WEIGHT_INDICES, pREWARD_PROJECT, 
                            pMAX_STATE);
        model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES);
        for (int i=0; i<N_TEAMS; i++) {
            cout << model_graph.outputs.ExtrinsicOutputTeam[i]+model_graph.modelInputs_.RewardTeam[i]<<endl;
        }
        }
```

Graphical Model of the above code:
- reward_team and reward_project are represented by arrows on the left and right equality constraints, respectively.
- teams_weights_per_department are represented by red (weight $=2$) or blue (weight $=1$) edges.

![graphicalModel](./img/demoGraphicalModel.png)

##### **2. Pure Python code**

To run pure code using sample datasets:
- Navigate to ``IMPALIB/test/python/src``
- Run: 
```bash
    python3 main_pure_optimized.py --nITER=400 --filterFlag=True --alpha=0.9 --PPFlag=True --threshold=-0.0001
```

##### **3. C++ code with a python wrapper**
To compile the C++ library and install the Python wrapper, navigate to the project root and use: 
```bash
    python3 -m pip install . -v
```

- Navigate to ``IMPALIB/src/impa``

To run wrapper code using sample datasets:
- Run:  
```bash 
    python3 main_wrapper_optimized.py --nITER=400 --filterFlag=True --alpha=0.9 --PPFlag=True --PPOption=1 --threshold=-0.0001 
```

**Note**: currently this option looks for a relevant sample dataset in the `data` directory, one directory up from the current working directory.
This will be fixed in a future version.

### **Requirements and Installation**
- A C++ $11$ -compatible compiler
- Python $3.9.7$
- To perform unit testing: randomized simulations using a pure python code and a python wrapper around a C++ code are carried out. A checking routine on the stored numpy files is executed to compare results. An external library called [cnpy](https://github.com/rogersce/cnpy) is used to save and load numpy arrays in C++

### **License**
Distributed under the MIT License.
See accompanying file [`LICENSE`](https://github.com/RustomAlexios/IMPALIB/blob/main/LICENSE) or at
[https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT))


