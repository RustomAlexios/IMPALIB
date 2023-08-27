// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"
#include "ut_utils.hpp"

void ut_iterate(string);
void ut_iterate_sample_graph(string);

void ut_iterate(string ut_name){

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const char *n_iter_bash=getenv("N_ITER");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";} 

    const int N_PROJECTS = atoi(n_projects_bash);
    const bool FILT_FLAG(filt_flag_bash);
    const int N_ITER = atoi(n_iter_bash);

    const char *Nu_bash=getenv("Nu");
    if(Nu_bash == NULL){cout << "Nu_bash not available\n";}
    string Nu_string = Nu_bash;
    vector<int> Nu = take_int(Nu_string);

    const int N_DEPARTMENTS = Nu.size();

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModel model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    fstream file_output1("../ut_results/ut_GraphicalModel/ut_"+ ut_name+"/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << endl;}

    fstream file_output2("../ut_results/ut_GraphicalModel/ut_"+ ut_name+"/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << endl;}
}

void ut_iterate_sample_graph(string ut_name){

    cnpy::NpyArray input_projects = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/N_PROJECTS_pure.npy");
    int* n_projects_pure = input_projects.data<int>();
    const int N_PROJECTS = *n_projects_pure;

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const char *n_iter_bash=getenv("N_ITER");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";} 

    const bool FILT_FLAG(filt_flag_bash);
    const int N_ITER = atoi(n_iter_bash);

    cnpy::NpyArray input_departments = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/N_DEPARTMENTS_pure.npy");
    int* n_departments_pure = input_departments.data<int>();
    const int N_DEPARTMENTS = *n_departments_pure;

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModel model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/ut_"+ ut_name+"/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    fstream file_output1("../ut_results/ut_GraphicalModel/ut_"+ ut_name+"/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << endl;}

    fstream file_output2("../ut_results/ut_GraphicalModel/ut_"+ ut_name+"/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << endl;}
}