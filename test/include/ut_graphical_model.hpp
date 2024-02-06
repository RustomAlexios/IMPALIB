// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"
#include "ut_utils.hpp"

void ut_iterate_kc_mwm(string);
void ut_iterate_sample_graph_kc_mwm(string);

void ut_iterate_kc_mwm(string ut_name){

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

    //cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/alpha.npy");
    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    //cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/N_TEAMS_pure.npy");
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    //cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/non_zero_weight_indices_sizes_pure.npy");
    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    //cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/reward_team_pure.npy");
    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    //cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/team_to_knapsack_m_pure.npy");
    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    //cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/teams_weights_per_department_pure.npy");
    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    //cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/non_zero_weight_indices_arr_pure.npy");
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    //cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/reward_project_pure.npy");
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    //cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/max_state_py.npy");
    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    //fstream file_output1("../ut_results/ut_GraphicalModel/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    fstream file_output1("../ut_results/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << endl;}

    //fstream file_output2("../ut_results/ut_GraphicalModel/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    fstream file_output2("../ut_results/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << endl;}
}

void ut_iterate_sample_graph_kc_mwm(string ut_name){

    //cnpy::NpyArray input_projects = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/N_PROJECTS_pure.npy");
    cnpy::NpyArray input_projects = cnpy::npy_load("../ut_inputs/N_PROJECTS_pure.npy");
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

    //cnpy::NpyArray input_departments = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/N_DEPARTMENTS_pure.npy");
    cnpy::NpyArray input_departments = cnpy::npy_load("../ut_inputs/N_DEPARTMENTS_pure.npy");
    int* n_departments_pure = input_departments.data<int>();
    const int N_DEPARTMENTS = *n_departments_pure;

    //cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/alpha.npy");
    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    //cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/N_TEAMS_pure.npy");
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    //cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/non_zero_weight_indices_sizes_pure.npy");
    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    //cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/reward_team_pure.npy");
    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    //cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/team_to_knapsack_m_pure.npy");
    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    //cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/teams_weights_per_department_pure.npy");
    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    //cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/non_zero_weight_indices_arr_pure.npy");
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    //cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/reward_project_pure.npy");
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    //cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/max_state_py.npy");
    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    //fstream file_output1("../ut_results/ut_GraphicalModel/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    fstream file_output1("../ut_results/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << endl;}

    //fstream file_output2("../ut_results/ut_GraphicalModel/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    fstream file_output2("../ut_results/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << endl;}
}

void ut_model_graph_tsp(string);

void ut_model_graph_tsp(string ut_name){

    const char *n_nodes_bash=getenv("N_NODES");
    if(n_nodes_bash == NULL)
    {cout << "n_nodes_bash not available\n";}

    const char *n_iterations_bash=getenv("N_ITER");
    if(n_iterations_bash == NULL)
    {cout << "n_iterations_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}

    const char *sym_flag_bash=getenv("SYM_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "sym_flag_bash not available\n";}
    
    const int N_NODES = atoi(n_nodes_bash);  
    const int N_ITER = atoi(n_iterations_bash);
    const int N_EDGE_VARIABLES = N_NODES*N_NODES-N_NODES;
    const bool FILT_FLAG(filt_flag_bash);
    const bool SYM_FLAG(sym_flag_bash);

    const bool RESET_FLAG = false;
    const int MAX_COUNT = 50;


    //cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/alpha.npy");
    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;
    
    //cnpy::NpyArray input_threshold = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/threshold.npy");
    cnpy::NpyArray input_threshold = cnpy::npy_load("../ut_inputs/threshold.npy");
    impalib_type* threshold_pure = input_threshold.data<impalib_type>();
    const impalib_type THRESHOLD = *threshold_pure;

    //cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/edge_connections_pure.npy");
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/edge_connections_pure.npy");
    int* edge_connections_pure = input1.data<int>();

    //cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/cost_edge_variable_pure.npy");
    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/cost_edge_variable_pure.npy");
    const impalib_type* cost_edge_variable_pure = input2.data<impalib_type>();

    //cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/cost_matrix_pure.npy");
    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/cost_matrix_pure.npy");
    const impalib_type* cost_matrix_pure = input3.data<impalib_type>();

    //cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/edge_ec_to_degree_constraint_m_pure.npy");
    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/edge_ec_to_degree_constraint_m_pure.npy");
    impalib_type* edge_ec_to_degree_constraint_m_pure = input4.data<impalib_type>();

    //cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_GraphicalModel/edge_degree_constraint_cost_pure.npy");
    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/edge_degree_constraint_cost_pure.npy");
    const impalib_type* edge_degree_constraint_cost_pure = input5.data<impalib_type>();


    if (ut_name == "IterateRelaxedGraph"){

        cout<<"-------"<<endl;
        cout<<"C++"<<endl;

        const bool AUGMENTATION_FLAG = false;

        GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

        model_graph.initialize(edge_connections_pure, cost_edge_variable_pure, cost_matrix_pure, edge_ec_to_degree_constraint_m_pure, edge_degree_constraint_cost_pure);

        model_graph.iterate_relaxed_graph();

        //fstream file_output("../ut_results/ut_GraphicalModel/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
        fstream file_output("../ut_results/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                    file_output.write((char*)(&model_graph.outputs.IntrinsicOutputEdgeEc[i]), sizeof(model_graph.outputs.IntrinsicOutputEdgeEc[i]));}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}

    }

    else if(ut_name == "IterateAugmentedGraph"){

        cout<<"-------"<<endl;
        cout<<"C++"<<endl;

        const bool AUGMENTATION_FLAG = true;

        const int MAX_AUGM_COUNT = 50;

        GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

        model_graph.initialize(edge_connections_pure, cost_edge_variable_pure, cost_matrix_pure, edge_ec_to_degree_constraint_m_pure, edge_degree_constraint_cost_pure);

        model_graph.iterate_relaxed_graph();

        if (!model_graph.subtourConstraintsSatisfiedFlag && AUGMENTATION_FLAG)
        {
            model_graph.perform_augmentation(MAX_AUGM_COUNT);

        }

        //fstream file_output("../ut_results/ut_GraphicalModel/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
        fstream file_output("../ut_results/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                    file_output.write((char*)(&model_graph.outputs.IntrinsicOutputEdgeEc[i]), sizeof(model_graph.outputs.IntrinsicOutputEdgeEc[i]));}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}


    }

}
