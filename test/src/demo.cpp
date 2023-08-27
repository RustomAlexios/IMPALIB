// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "../include/impalib_unit_tests.hpp"

//demonstration of the use of the library in c++
//
//export LD_LIBRARY_PATH=../../build/cnpylib
//g++ -o demo demo.cpp -L../../build/cnpylib -lcnpy -lz --std=c++11

//to run wrapper code using sample datasets
//location: "pythonWrapperOptimized/src"
//python3 main_wrapper_optimized.py --nITER=400 --filterFlag=True --alpha=0.9 --PPFlag=True --PPOption=1 --threshold=-0.0001

//to run pure code using sample datasets
//location: "pythonWrapperOptimized/test/python/src"
//python3 main_pure_optimized.py --nITER=400 --filterFlag=True --alpha=0.9 --PPFlag=True --PPOption=2 --threshold=-0.0001

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
    impalib_type reward_project[N_PROJECTS*N_TEAMS] = {44, 1, 41, 10, 3, //cost of assigning a team to a project
                                                    7, 17, 56, 98, 63};
    
    const impalib_type* pREWARD_PROJECT = reward_project;
    const impalib_type* pREWARD_TEAM = reward_team;

    impalib_type transition_model[N_DEPARTMENTS*N_TEAMS]= {-22, -31, -68, -39, 0,        
                                                    0, -31, -68, -39, -84};
    impalib_type* pTransition_model = transition_model;

    int teams_weights_per_department[N_DEPARTMENTS*N_TEAMS]= {2,1,1,1,0,
                                                            0,1,1,1,2};
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT = teams_weights_per_department;

    int non_zero_weight_indices[N_DEPARTMENTS*N_TEAMS] = {0,1,2,3,1,2,3,4}; 
    const int* p_NON_ZERO_WEIGHT_INDICES = non_zero_weight_indices;


    const int* pMAX_STATE = max_state;
    
    model_graph.initialize(pREWARD_TEAM, pTransition_model, pTEAMS_WEIGHTS_PER_DEPARTMENT,
                        pNON_ZERO_WEIGHT_INDICES_SIZES, p_NON_ZERO_WEIGHT_INDICES, pREWARD_PROJECT, 
                        pMAX_STATE);
    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES);
    //for (int i=0; i<N_TEAMS; i++){
    //    cout<<model_graph.outputs.ExtrinsicOutputTeam[i]+model_graph.modelInputs_.RewardTeam[i]<<endl;
    //}
}