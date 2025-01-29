// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"
#include "ut_utils.hpp"

void ut_iterate_kc_mwm(string&);
void ut_iterate_sample_graph_kc_mwm(string&);

void ut_iterate_kc_mwm(string& ut_name){

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
    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");
    // cout<<"FILT_FLAG: "<<FILT_FLAG<<"\n";
    const int N_ITER = atoi(n_iter_bash);

    const char *Nu_bash=getenv("Nu");
    if(Nu_bash == NULL){cout << "Nu_bash not available\n";}
    string Nu_string = Nu_bash;
    vector<int> Nu = take_int(Nu_string);

    const int N_DEPARTMENTS = static_cast<int>(Nu.size());

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    fstream file_output1("../ut_results/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << "\n";}

    fstream file_output2("../ut_results/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << "\n";}
}

void ut_iterate_sample_graph_kc_mwm(string& ut_name){

    cnpy::NpyArray input_projects = cnpy::npy_load("../ut_inputs/N_PROJECTS_pure.npy");
    int* n_projects_pure = input_projects.data<int>();
    const int N_PROJECTS = *n_projects_pure;

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const char *n_iter_bash=getenv("N_ITER");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";} 

    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");
    // cout<<"FILT_FLAG: "<<FILT_FLAG<<"\n";
    const int N_ITER = atoi(n_iter_bash);

    cnpy::NpyArray input_departments = cnpy::npy_load("../ut_inputs/N_DEPARTMENTS_pure.npy");
    int* n_departments_pure = input_departments.data<int>();
    const int N_DEPARTMENTS = *n_departments_pure;

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_sizes_pure.npy");
    const int* pNON_ZERO_WEIGHT_INDICES_SIZES_PY = input2.data<int>();

    int max_size_non_zero_weight = *max_element(pNON_ZERO_WEIGHT_INDICES_SIZES_PY , pNON_ZERO_WEIGHT_INDICES_SIZES_PY + N_DEPARTMENTS);
    
    GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    const impalib_type* pREWARD_TEAM_PY = input3.data<impalib_type>();

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* pTransition_model_py = input4.data<impalib_type>();

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    const int* pTEAMS_WEIGHTS_PER_DEPARTMENT_PY = input5.data<int>();
    
    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    const int* p_NON_ZERO_WEIGHT_INDICES_PY = input6.data<int>();
    
    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    const impalib_type* pREWARD_PROJECT_PY = input7.data<impalib_type>();

    cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/max_state_py.npy");
    const int* pMAX_STATE_PY = input8.data<int>();

    model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

    fstream file_output1("../ut_results/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output1.is_open()) {
        for (int i=0; i<N_TEAMS; i++){
            file_output1.write((char*)(&model_graph.outputs.ExtrinsicOutputTeam[i]), sizeof(model_graph.outputs.ExtrinsicOutputTeam[i]));}
            file_output1.close();}
    else {cout << "Error! File cannot be opened!" << "\n";}

    fstream file_output2("../ut_results/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
    if (file_output2.is_open()) {
        for (int i=0; i<N_TEAMS*N_PROJECTS; i++){
            file_output2.write((char*)(&model_graph.outputs.IntrinsicOutMwm[i]), sizeof(model_graph.outputs.IntrinsicOutMwm[i]));}
            file_output2.close();}
    else {cout << "Error! File cannot be opened!" << "\n";}
}

void ut_model_graph_tsp(string&);

void ut_model_graph_tsp(string& ut_name){

    const char *n_nodes_bash=getenv("N_NODES");
    if(n_nodes_bash == NULL)
    {cout << "n_nodes_bash not available\n";}

    const char *n_iterations_bash=getenv("N_ITER");
    if(n_iterations_bash == NULL)
    {cout << "n_iterations_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const int N_NODES = atoi(n_nodes_bash);  
    const int N_ITER = atoi(n_iterations_bash);
    const int N_EDGE_VARIABLES = N_NODES*N_NODES-N_NODES;
    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");
    // cout<<"FILT_FLAG: "<<FILT_FLAG<<"\n";
    const bool RESET_FLAG = false;
    const int MAX_COUNT = 50;


    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;
    
    cnpy::NpyArray input_threshold = cnpy::npy_load("../ut_inputs/threshold.npy");
    impalib_type* threshold_pure = input_threshold.data<impalib_type>();
    const impalib_type THRESHOLD = *threshold_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/edge_connections_pure.npy");
    int* edge_connections_pure = input1.data<int>();

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/cost_edge_variable_pure.npy");
    const impalib_type* cost_edge_variable_pure = input2.data<impalib_type>();

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/cost_matrix_pure.npy");
    const impalib_type* cost_matrix_pure = input3.data<impalib_type>();

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/edge_ec_to_degree_constraint_m_pure.npy");
    impalib_type* edge_ec_to_degree_constraint_m_pure = input4.data<impalib_type>();

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/edge_degree_constraint_cost_pure.npy");
    const impalib_type* edge_degree_constraint_cost_pure = input5.data<impalib_type>();


    if (ut_name == "IterateRelaxedGraph"){

        cout<<"-------"<<"\n";
        cout<<"C++"<<"\n";

        const bool AUGMENTATION_FLAG = false;

        GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

        model_graph.initialize(edge_connections_pure, cost_edge_variable_pure, cost_matrix_pure, edge_ec_to_degree_constraint_m_pure, edge_degree_constraint_cost_pure);

        model_graph.iterate_relaxed_graph();

        fstream file_output("../ut_results/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                    file_output.write((char*)(&model_graph.outputs.IntrinsicOutputEdgeEc[i]), sizeof(model_graph.outputs.IntrinsicOutputEdgeEc[i]));}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}

    }

    else if(ut_name == "IterateAugmentedGraph"){

        cout<<"-------"<<"\n";
        cout<<"C++"<<"\n";

        const bool AUGMENTATION_FLAG = true;

        const int MAX_AUGM_COUNT = 50;

        GraphicalModelTsp model_graph(N_ITER, N_NODES, N_EDGE_VARIABLES, AUGMENTATION_FLAG, RESET_FLAG, FILT_FLAG, ALPHA, THRESHOLD, MAX_COUNT);

        model_graph.initialize(edge_connections_pure, cost_edge_variable_pure, cost_matrix_pure, edge_ec_to_degree_constraint_m_pure, edge_degree_constraint_cost_pure);

        model_graph.iterate_relaxed_graph();

        if (!model_graph.subtourConstraintsSatisfiedFlag && AUGMENTATION_FLAG)
        {
            model_graph.perform_augmentation(MAX_AUGM_COUNT);

        }

        fstream file_output("../ut_results/intrinsic_out_edge_ec_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                    file_output.write((char*)(&model_graph.outputs.IntrinsicOutputEdgeEc[i]), sizeof(model_graph.outputs.IntrinsicOutputEdgeEc[i]));}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}


    }

}


void ut_model_graph_ksat(string&);

void ut_model_graph_ksat(string& ut_name){

    const char *n_variables_bash=getenv("NUM_VARIABLES");
    if(n_variables_bash == NULL)
    {cout << "n_variables_bash not available\n";}

    const char *k_variable_bash=getenv("K_VARIABLE");
    if(k_variable_bash == NULL)
    {cout << "k_variable_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}

    const char *n_iter_bash=getenv("N_ITER");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";}

    const int NUM_VARIABLES = atoi(n_variables_bash);  
    const int K_VARIABLE = atoi(k_variable_bash);
    const int N_ITER = atoi(n_iter_bash);
    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");
    // cout<<"FILT_FLAG: "<<FILT_FLAG<<"\n";

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input_num_constraints = cnpy::npy_load("../ut_inputs/num_constraints.npy");
    int* num_constraints_pure = input_num_constraints.data<int>();
    const int NUM_CONSTRAINTS = *num_constraints_pure;

    cnpy::NpyArray input_num_used_variables = cnpy::npy_load("../ut_inputs/size_used_variables_pure.npy");
    int* size_used_variables_pure = input_num_used_variables.data<int>();
    const int NUM_USED_VARIABLES = *size_used_variables_pure;

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/used_variables_pure.npy");
    int* used_variables_pure = input1.data<int>(); 

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/variables_connections_pure.npy");
    int* variables_connections_pure = input2.data<int>(); 

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/sizes_of_variables_connections_pure.npy");
    int* sizes_of_variables_connections_pure = input3.data<int>(); 

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/constraints_connections_pure.npy");
    int* constraints_connections_pure = input4.data<int>(); 

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/constraints_connections_type_pure.npy");
    int* constraints_connections_type_pure = input5.data<int>(); 

    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/incoming_metrics_cost_pure.npy");
    impalib_type* incoming_metrics_cost_pure = input6.data<impalib_type>(); 

    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/variable_ec_to_ksat_constraint_m_pure.npy");
    impalib_type* variable_ec_to_ksat_constraint_m_pure = input7.data<impalib_type>(); 

    if (ut_name == "Iterate"){

        cout<<"-------"<<"\n";
        cout<<"C++"<<"\n";

        GraphicalModelKsat model_graph(N_ITER, NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILT_FLAG, ALPHA, NUM_USED_VARIABLES);

        model_graph.initialize(used_variables_pure, variables_connections_pure, sizes_of_variables_connections_pure, constraints_connections_pure, constraints_connections_type_pure, \
                                incoming_metrics_cost_pure, variable_ec_to_ksat_constraint_m_pure);

        model_graph.iterate();

        fstream file_output("../ut_results/extrinsic_out_variable_ec_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<NUM_VARIABLES; i++){
                    file_output.write((char*)(&model_graph.outputs.ExtrinsicOutputVariableEc[i]), sizeof(model_graph.outputs.ExtrinsicOutputVariableEc[i]));}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}

    }

}


void ut_model_graph_mobarp(string&);

void ut_model_graph_mobarp(string& ut_name){

    cout<<"-------"<<"\n";
    cout<<"C++"<<"\n";

    const char *n_fixed_tx_bash=getenv("NUM_FIXED_TX");
    if(n_fixed_tx_bash == NULL)
    {cout << "n_fixed_tx_bash not available\n";}

    const char *n_mobile_tx_bash=getenv("NUM_MOBILE_TX");
    if(n_mobile_tx_bash == NULL)
    {cout << "n_mobile_tx_bash not available\n";}

    const char *n_bands_bash=getenv("NUM_BANDS");
    if(n_bands_bash == NULL)
    {cout << "n_bands_bash not available\n";}

    const char *n_time_steps_bash=getenv("NUM_TIME_STEPS");
    if(n_time_steps_bash == NULL)
    {cout << "n_time_steps_bash not available\n";}

    const char *n_rx_locs_bash=getenv("NUM_RX_LOCS");
    if(n_rx_locs_bash == NULL)
    {cout << "n_rx_locs_bash not available\n";}

    const char *n_mobile_tx_locs_bash=getenv("NUM_MOBILE_TX_LOCS");
    if(n_mobile_tx_locs_bash == NULL)
    {cout << "n_mobile_tx_locs_bash not available\n";}

    const char *exclude_capac_flag_bash=getenv("EXCLUDE_CAP_FLAG");
    if(exclude_capac_flag_bash == NULL)
    {cout << "exclude_capac_flag_bash not available\n";}

    const char *n_iter_bash=getenv("NUM_ITERATIONS");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";}

    const int NUM_FIXED_TX = atoi(n_fixed_tx_bash);  
    const int NUM_MOBILE_TX = atoi(n_mobile_tx_bash); 
    const int NUM_BANDS = atoi(n_bands_bash); 
    const int NUM_TIME_STEPS = atoi(n_time_steps_bash); 
    const int NUM_RX_LOCS = atoi(n_rx_locs_bash); 
    const int NUM_MOBILE_TX_LOCS = atoi(n_mobile_tx_locs_bash);

    // const bool EXCLUDE_CAP_FLAG(exclude_capac_flag_bash);
    const bool EXCLUDE_CAP_FLAG = (exclude_capac_flag_bash != NULL && std::string(exclude_capac_flag_bash) == "1");
    const int NUM_ITERATIONS = atoi(n_iter_bash);

    // cout<<"exclude_capac_flag_bash: "<<exclude_capac_flag_bash<<", EXCLUDE_CAP_FLAG: "<<EXCLUDE_CAP_FLAG<<"\n";

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
        {cout << "filt_flag_bash not available\n";}
    
    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");
    // cout<<"filt_flag_bash: "<<filt_flag_bash<<" "<<", FILT_FLAG: "<<FILT_FLAG<<"\n";

    cout << "num_iterations: " << NUM_ITERATIONS << "\n";
    cout << "num_fixed_tx: " << NUM_FIXED_TX << "\n";
    cout << "num_mobile_tx: " << NUM_MOBILE_TX << "\n";
    cout << "num_bands: " << NUM_BANDS << "\n";
    cout << "num_time_steps: " << NUM_TIME_STEPS << "\n";
    cout << "num_rx_locs: " << NUM_RX_LOCS << "\n";
    cout << "num_mobile_tx_locs: " << NUM_MOBILE_TX_LOCS << "\n";
    cout << "filtering_flag: " << (FILT_FLAG ? "True" : "False") << "\n";
    cout << "alpha: " << ALPHA << "\n";
    cout << "exclude_cap_flag: " << (EXCLUDE_CAP_FLAG ? "True" : "False")<< "\n";

    if (ut_name == "Iterate"){

        GraphicalModelMOBARP model_graph(NUM_ITERATIONS, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, ALPHA, FILT_FLAG, EXCLUDE_CAP_FLAG);

        cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/fixed_x_costs.npy");
        impalib_type* fixed_x_costs_pure = input1.data<impalib_type>();

        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/mobile_x_costs.npy");
        impalib_type* mobile_x_costs_pure = input2.data<impalib_type>();

        cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/z_costs.npy");
        impalib_type* z_costs_pure = input3.data<impalib_type>();

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/r_costs.npy");
        impalib_type* r_costs_pure = input4.data<impalib_type>(); 

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/connectivity_mobile_tx.npy");
        int* connectivity_mobile_tx_pure = input5.data<int>();

        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/connectivity_fixed_tx.npy");
        int* connectivity_fixed_tx_pure = input6.data<int>();

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/fixed_capacity_constraints.npy");
        int* fixed_capacity_constraints_pure = input7.data<int>();

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/mobile_capacity_constraints.npy");
        int* mobile_capacity_constraints_pure = input8.data<int>();

        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy");
        int* conx_fixed_tx_per_num_rx_locs_pure = input9.data<int>();

        cnpy::NpyArray input10 = cnpy::npy_load("../ut_inputs/conx_mob_tx_rx.npy");
        int* conx_mob_tx_rx_pure = input10.data<int>();

        cnpy::NpyArray input11 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input11.data<int>();
        
        cnpy::NpyArray input12 = cnpy::npy_load("../ut_inputs/fixed_x_eq_const_to_fixed_capac_const_m.npy");
        impalib_type* fixed_x_eq_const_to_fixed_capac_const_m_pure = input12.data<impalib_type>();

        //yes
        cnpy::NpyArray input13 = cnpy::npy_load("../ut_inputs/mobile_x_eq_const_to_mobile_capac_const_m.npy");
        impalib_type* mobile_x_eq_const_to_mobile_capac_const_m_pure = input13.data<impalib_type>();

        //yes
        cnpy::NpyArray input14 = cnpy::npy_load("../ut_inputs/r_eq_const_to_auxiliary_const_m.npy");
        impalib_type* r_eq_const_to_auxiliary_const_m_pure = input14.data<impalib_type>();


        model_graph.initialize(fixed_x_eq_const_to_fixed_capac_const_m_pure, mobile_x_eq_const_to_mobile_capac_const_m_pure,
                                                    r_eq_const_to_auxiliary_const_m_pure, fixed_x_costs_pure, mobile_x_costs_pure,
                                                    z_costs_pure, r_costs_pure, connectivity_fixed_tx_pure,
                                                    connectivity_mobile_tx_pure, fixed_capacity_constraints_pure, mobile_capacity_constraints_pure,
                                                    conx_mob_tx_per_num_mob_tx_locs_pure, conx_fixed_tx_per_num_rx_locs_pure, conx_mob_tx_rx_pure);

        model_graph.iterate();

        // model_graph.process_ouputs(pExtrinsic_fixed_x, pExtrinsic_mobile_x, pExtrinsic_z, pExtrinsic_r)


        fstream file_output_1("../ut_results/extrinsic_fixed_x_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<model_graph.outputs.ExtrinsicFixedX.size(); i++){
                file_output_1.write((char*)(&model_graph.outputs.ExtrinsicFixedX[i]), sizeof(model_graph.outputs.ExtrinsicFixedX[i]));}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/extrinsic_mobile_x_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<model_graph.outputs.ExtrinsicMobileX.size(); i++){
                file_output_2.write((char*)(&model_graph.outputs.ExtrinsicMobileX[i]), sizeof(model_graph.outputs.ExtrinsicMobileX[i]));}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_3("../ut_results/extrinsic_r_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_3.is_open()) {
            for (int i=0; i<model_graph.outputs.ExtrinsicR.size(); i++){
                file_output_3.write((char*)(&model_graph.outputs.ExtrinsicR[i]), sizeof(model_graph.outputs.ExtrinsicR[i]));}
                file_output_3.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_4("../ut_results/extrinsic_z_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_4.is_open()) {
            for (int i=0; i<model_graph.outputs.ExtrinsicZ.size(); i++){
            for (int j=0; j<model_graph.outputs.ExtrinsicZ[0].size(); j++){
                file_output_4.write((char*)(&model_graph.outputs.ExtrinsicZ[i][j]), sizeof(model_graph.outputs.ExtrinsicZ[i][j]));}}
                file_output_4.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}


    }

}