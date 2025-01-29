// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_input_output_kc_mwm(string&);

void ut_input_output_kc_mwm(string& ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const int N_DEPARTMENT = atoi(n_departments_bash);  
    const int N_PROJECTS = atoi(n_projects_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    impalib_type* reward_project_pure = input2.data<impalib_type>();
    vector<vector<impalib_type>> reward_project(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int project_index=0; project_index<N_PROJECTS; project_index++){
    copy ( reward_project_pure + N_TEAMS*project_index, reward_project_pure+N_TEAMS*(project_index+1), reward_project[project_index].begin() );
    }

    OutputsKcMwm outputs(N_DEPARTMENT, N_TEAMS, N_PROJECTS);

    if (ut_name == "ExtrinsicOutputTeamUpdate"){
        
        cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/extrinsic_output_department_pure.npy");
        impalib_type* extrinsic_output_department_pure = input3.data<impalib_type>();
        vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENT, vector<impalib_type>(N_TEAMS,zero_value));

        for (int department_index=0; department_index<N_DEPARTMENT; department_index++){
        copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure+N_TEAMS*(department_index+1), extrinsic_output_department[department_index].begin() );
        }

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/oric_to_team_m_pure.npy");
        impalib_type* oric_to_team_m_pure = input4.data<impalib_type>();
        vector<impalib_type> oric_to_team_m(N_TEAMS,zero_value);
        copy(oric_to_team_m_pure, oric_to_team_m_pure + N_TEAMS, oric_to_team_m.begin());

        
        outputs.extrinsic_output_team_update(extrinsic_output_department, oric_to_team_m);

        fstream file_output("../ut_results/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_TEAMS; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputTeam[i]), sizeof(outputs.ExtrinsicOutputTeam[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

    if (ut_name == "IntrinsicOutMwmUpdate"){
        
        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/oric_to_eq_constraint_m_pure.npy");
        impalib_type* oric_to_eq_constraint_m_pure = input4.data<impalib_type>();
        vector<vector<impalib_type>> oric_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( oric_to_eq_constraint_m_pure + N_TEAMS*project_index, oric_to_eq_constraint_m_pure+N_TEAMS*(project_index+1), oric_to_eq_constraint_m[project_index].begin() );
        }

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/project_to_eq_constraint_m_pure.npy");
        impalib_type* project_to_eq_constraint_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( project_to_eq_constraint_m_pure + N_TEAMS*project_index, project_to_eq_constraint_m_pure+N_TEAMS*(project_index+1), project_to_eq_constraint_m[project_index].begin() );
        }

        outputs.intrinsic_out_mwm_update(oric_to_eq_constraint_m, project_to_eq_constraint_m, reward_project);

        fstream file_output("../ut_results/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_PROJECTS*N_TEAMS; i++){
                        file_output.write((char*)(&outputs.IntrinsicOutMwm[i]), sizeof(outputs.IntrinsicOutMwm[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}

    }

}

void ut_input_output_tsp(string&);

void ut_input_output_tsp(string& ut_name){

    const char *n_nodes_bash=getenv("N_NODES");
    if(n_nodes_bash == NULL)
    {cout << "n_nodes_bash not available\n";}

    const char *n_subtours_bash=getenv("N_SUBTOURS");
    if(n_subtours_bash == NULL)
    {cout << "n_subtours_bash not available\n";}

    const int N_NODES = atoi(n_nodes_bash);  
    const int N_EDGE_VARIABLES = N_NODES*N_NODES-N_NODES;
    const int N_SUBTOURS = atoi(n_subtours_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/degree_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* degree_constraint_to_eq_constraint_m_pure = input1.data<impalib_type>();

    vector<vector<impalib_type>> degree_constraint_to_eq_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (degree_constraint_to_eq_constraint_m_pure + N_NODES*edge_variable_index, degree_constraint_to_eq_constraint_m_pure+N_NODES*(edge_variable_index+1), degree_constraint_to_eq_constraint_m[edge_variable_index].begin() );
    }

    OutputsTsp outputs(N_NODES, N_EDGE_VARIABLES);

    if (ut_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate"){
        
        outputs.extrinsic_output_edge_ec_relaxed_graph_update(degree_constraint_to_eq_constraint_m);
        
        fstream file_output("../ut_results/extrinsic_output_edge_ec_relaxed_graph_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_EDGE_VARIABLES; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputEdgeEc[i]), sizeof(outputs.ExtrinsicOutputEdgeEc[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

    if (ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"){
        
        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input2.data<impalib_type>();
        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy ( subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        outputs.extrinsic_output_edge_ec_augmented_graph_update(degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m);

        fstream file_output("../ut_results/extrinsic_output_edge_ec_augmented_graph_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputEdgeEc[i]), sizeof(outputs.ExtrinsicOutputEdgeEc[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}
    }

}

void ut_input_output_ksat(string&);

void ut_input_output_ksat(string& ut_name){

    const char *n_variables_bash=getenv("NUM_VARIABLES");
    if(n_variables_bash == NULL)
    {cout << "n_variables_bash not available\n";}

    const char *n_constraints_bash=getenv("NUM_CONSTRAINTS");
    if(n_constraints_bash == NULL)
    {cout << "n_constraints_bash not available\n";}

    const char *k_variable_bash=getenv("K_VARIABLE");
    if(k_variable_bash == NULL)
    {cout << "k_variable_bash not available\n";}

    const int NUM_VARIABLES = atoi(n_variables_bash);  
    const int NUM_CONSTRAINTS = atoi(n_constraints_bash);
    const int K_VARIABLE = atoi(k_variable_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ksat_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* ksat_constraint_to_eq_constraint_m_pure = input1.data<impalib_type>();

    vector<vector<impalib_type>> ksat_constraint_to_eq_constraint_m(NUM_CONSTRAINTS, vector<impalib_type>(NUM_VARIABLES,zero_value));

    for (int constraint_index=0; constraint_index<NUM_CONSTRAINTS; constraint_index++){
    copy (ksat_constraint_to_eq_constraint_m_pure + NUM_VARIABLES*constraint_index, ksat_constraint_to_eq_constraint_m_pure+NUM_VARIABLES*(constraint_index+1), ksat_constraint_to_eq_constraint_m[constraint_index].begin() );
    }

    OutputsKsat outputs(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE);

    if (ut_name == "ExtrinsicOutputVariableEcUpdate"){
        
        outputs.update_extrinsic(ksat_constraint_to_eq_constraint_m);
        
        fstream file_output("../ut_results/extrinsic_output_variable_ec_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<NUM_VARIABLES; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputVariableEc[i]), sizeof(outputs.ExtrinsicOutputVariableEc[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

}


void ut_input_output_mobarp(string&);

void ut_input_output_mobarp(string& ut_name){
    
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

    const int NUM_FIXED_TX = atoi(n_fixed_tx_bash);  
    const int NUM_MOBILE_TX = atoi(n_mobile_tx_bash); 
    const int NUM_BANDS = atoi(n_bands_bash); 
    const int NUM_TIME_STEPS = atoi(n_time_steps_bash); 
    const int NUM_RX_LOCS = atoi(n_rx_locs_bash); 
    const int NUM_MOBILE_TX_LOCS = atoi(n_mobile_tx_locs_bash); 
    // const bool EXCLUDE_CAP_FLAG(exclude_capac_flag_bash);
    const bool EXCLUDE_CAP_FLAG = (exclude_capac_flag_bash != NULL && std::string(exclude_capac_flag_bash) == "1");

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/fixed_x_costs.npy");
    impalib_type* fixed_x_costs_pure = input1.data<impalib_type>();
    vector<vector<vector<impalib_type>>> FixedTxCosts;

    for (int i=0; i< NUM_FIXED_TX; i++){
        FixedTxCosts.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        for (int j=0; j< NUM_BANDS; j++){
            copy(fixed_x_costs_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, fixed_x_costs_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, FixedTxCosts[i][j].begin());
        }
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/mobile_x_costs.npy");
    impalib_type* mobile_x_costs_pure = input2.data<impalib_type>();
    vector<vector<vector<impalib_type>>> MobileTxCosts;

    for (int i=0; i< NUM_MOBILE_TX; i++){
        MobileTxCosts.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        for (int j=0; j< NUM_BANDS; j++){
            copy(mobile_x_costs_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, mobile_x_costs_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, MobileTxCosts[i][j].begin());
        }
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/r_costs.npy");
    impalib_type* r_costs_pure = input3.data<impalib_type>();
    vector<vector<impalib_type>> RCosts;

    for (int i=0; i< NUM_MOBILE_TX; i++){
        RCosts.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
        copy(r_costs_pure + NUM_MOBILE_TX_LOCS * i, r_costs_pure + NUM_MOBILE_TX_LOCS * (i + 1), RCosts[i].begin());
    }

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/connectivity_mobile_tx.npy");
    int* connectivity_mobile_tx_pure = input4.data<int>();
    vector<vector<vector<vector<int>>>> ConnectivityMobileTx;

    for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
        ConnectivityMobileTx.push_back(vector<vector<vector<int>>>(NUM_BANDS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_RX_LOCS, 0))));
        for (int j=0; j<NUM_BANDS; j++){
            for (int k=0; k< NUM_TIME_STEPS; k++){
                copy(connectivity_mobile_tx_pure + NUM_RX_LOCS*k + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*n, connectivity_mobile_tx_pure + NUM_RX_LOCS*(k+1) + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*n, ConnectivityMobileTx[n][j][k].begin());
        }
        }
    }

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/connectivity_fixed_tx.npy");
    int* connectivity_fixed_tx_pure = input5.data<int>();
    vector<vector<vector<vector<int>>>> ConnectivityFixedTx;

    for (int i=0; i< NUM_FIXED_TX; i++){
        ConnectivityFixedTx.push_back(vector<vector<vector<int>>>(NUM_BANDS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_RX_LOCS, 0))));
        for (int j=0; j<NUM_BANDS; j++){
            for (int k=0; k< NUM_TIME_STEPS; k++){
                copy(connectivity_fixed_tx_pure + NUM_RX_LOCS*k + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*i, connectivity_fixed_tx_pure + NUM_RX_LOCS*(k+1) + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*i, ConnectivityFixedTx[i][j][k].begin());
        }
        }
    }


    OutputsMOBARP outputs(NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, EXCLUDE_CAP_FLAG);

    if (ut_name == "ExtrinsicUpdate"){
        
        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/fixed_capac_const_to_fixed_x_eq_const_m_pure.npy");
        impalib_type* fixed_capac_const_to_fixed_x_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<vector<impalib_type>>> fixed_capac_const_to_fixed_x_eq_const_m;

        for (int i=0; i< NUM_FIXED_TX; i++){
            fixed_capac_const_to_fixed_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(fixed_capac_const_to_fixed_x_eq_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, fixed_capac_const_to_fixed_x_eq_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, fixed_capac_const_to_fixed_x_eq_const_m[i][j].begin());
            }
        }

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/mobile_capac_const_to_mobile_x_eq_const_m_pure.npy");
        impalib_type* mobile_capac_const_to_mobile_x_eq_const_m_pure = input7.data<impalib_type>();
        vector<vector<vector<impalib_type>>> mobile_capac_const_to_mobile_x_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            mobile_capac_const_to_mobile_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(mobile_capac_const_to_mobile_x_eq_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, mobile_capac_const_to_mobile_x_eq_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, mobile_capac_const_to_mobile_x_eq_const_m[i][j].begin());
            }
        }

        
        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_mobile_x_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_mobile_x_eq_const_m_pure = input8.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_mobile_x_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_mobile_x_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_mobile_x_eq_const_m[i].begin());
        }

        
        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_fixed_x_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_fixed_x_eq_const_m_pure = input9.data<impalib_type>();
        vector<vector<vector<impalib_type>>> set_cover_ineq_const_to_fixed_x_eq_const_m;

        for (int i=0; i< NUM_TIME_STEPS; i++){
            set_cover_ineq_const_to_fixed_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0)));
            for (int j=0; j< NUM_RX_LOCS; j++){
                copy(set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * j + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * (j + 1) + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m[i][j].begin());
            }
        }

        cnpy::NpyArray input10 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_r_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_r_eq_const_m_pure = input10.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_r_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_r_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_r_eq_const_m[i].begin());
        }

        cnpy::NpyArray input11 = cnpy::npy_load("../ut_inputs/mobile_loc_eq_const_to_r_eq_const_m_pure.npy");
        impalib_type* mobile_loc_eq_const_to_r_eq_const_m_pure = input11.data<impalib_type>();
        vector<vector<impalib_type>> mobile_loc_eq_const_to_r_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            mobile_loc_eq_const_to_r_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(mobile_loc_eq_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, mobile_loc_eq_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), mobile_loc_eq_const_to_r_eq_const_m[i].begin());
        }

        cnpy::NpyArray input12 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_z_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_z_eq_const_m_pure = input12.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_z_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_z_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_z_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_z_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_z_eq_const_m[i].begin());
        }

        cnpy::NpyArray input13 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_z_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_z_eq_const_m_pure = input13.data<impalib_type>();
        vector<vector<vector<vector<impalib_type>>>> set_cover_ineq_const_to_z_eq_const_m;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            set_cover_ineq_const_to_z_eq_const_m.push_back(vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m[k][l][n].begin());
            }
            }
        }

        cnpy::NpyArray input14 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input14.data<int>();
        vector<vector<int>> conx_mob_tx_per_num_mob_tx_locs;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_mob_tx_per_num_mob_tx_locs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
            copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), conx_mob_tx_per_num_mob_tx_locs[i].begin());
        }

        cnpy::NpyArray input15 = cnpy::npy_load("../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy");
        int* conx_fixed_tx_per_num_rx_locs_pure = input15.data<int>();
        vector<vector<int>> conx_fixed_tx_per_num_rx_locs;

        for (int i=0; i< NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_fixed_tx_per_num_rx_locs.push_back(vector<int>(NUM_RX_LOCS, 0));
            copy(conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * i, conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * (i + 1), conx_fixed_tx_per_num_rx_locs[i].begin());
        }

        cnpy::NpyArray input16 = cnpy::npy_load("../ut_inputs/conx_mob_tx_rx.npy");
        int* conx_mob_tx_rx_pure = input16.data<int>();
        vector<vector<vector<vector<int>>>> conx_mob_tx_rx;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            conx_mob_tx_rx.push_back(vector<vector<vector<int>>>(NUM_RX_LOCS, vector<vector<int>>(NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx[k][l][n].begin());
            }
            }
        }



        vector<vector<vector<int>>> reshaped_1(NUM_MOBILE_TX, vector<vector<int>>(NUM_BANDS*NUM_TIME_STEPS, vector<int>(NUM_MOBILE_TX_LOCS, 0)));

        for (size_t i = 0; i < NUM_MOBILE_TX; i++) {
            for (size_t j_k = 0; j_k < NUM_BANDS * NUM_TIME_STEPS; j_k++) {
                for (size_t n = 0; n < NUM_MOBILE_TX_LOCS; n++) {
                    reshaped_1[i][j_k][n] = conx_mob_tx_per_num_mob_tx_locs[i * NUM_BANDS * NUM_TIME_STEPS + j_k][n];
                }
            }
        }

        vector<vector<vector<int>>> reshaped_2(NUM_MOBILE_TX, vector<vector<int>>(NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_TIME_STEPS, 0)));


        for (size_t i = 0; i < NUM_MOBILE_TX; i++) {
            for (size_t n = 0; n < NUM_MOBILE_TX_LOCS; n++) {
                for (size_t j_k = 0; j_k < NUM_BANDS * NUM_TIME_STEPS; j_k++) {
                    reshaped_2[i][n][j_k] = reshaped_1[i][j_k][n];
                }
            }
        }

        vector<vector<int>> ConxMobTxR(vector<vector<int>>(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_TIME_STEPS, 0)));

        for (size_t i = 0; i < NUM_MOBILE_TX; i++) {
            for (size_t n = 0; n < NUM_MOBILE_TX_LOCS; n++) {
                for (size_t j_k = 0; j_k < NUM_BANDS * NUM_TIME_STEPS; j_k++) {
                    ConxMobTxR[i * NUM_MOBILE_TX_LOCS + n][j_k] = reshaped_2[i][n][j_k];
                }
            }
        }


        outputs.extrinsic_update(fixed_capac_const_to_fixed_x_eq_const_m, mobile_capac_const_to_mobile_x_eq_const_m, auxiliary_const_to_mobile_x_eq_const_m, set_cover_ineq_const_to_fixed_x_eq_const_m, auxiliary_const_to_r_eq_const_m,
                                mobile_loc_eq_const_to_r_eq_const_m, auxiliary_const_to_z_eq_const_m, set_cover_ineq_const_to_z_eq_const_m, conx_mob_tx_per_num_mob_tx_locs, conx_fixed_tx_per_num_rx_locs,
                                conx_mob_tx_rx, ConxMobTxR);

        fstream file_output_1("../ut_results/extrinsic_fixed_x_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<outputs.ExtrinsicFixedX.size(); i++){
                file_output_1.write((char*)(&outputs.ExtrinsicFixedX[i]), sizeof(outputs.ExtrinsicFixedX[i]));}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/extrinsic_mobile_x_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<outputs.ExtrinsicMobileX.size(); i++){
                file_output_2.write((char*)(&outputs.ExtrinsicMobileX[i]), sizeof(outputs.ExtrinsicMobileX[i]));}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_3("../ut_results/extrinsic_r_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_3.is_open()) {
            for (int i=0; i<outputs.ExtrinsicR.size(); i++){
                file_output_3.write((char*)(&outputs.ExtrinsicR[i]), sizeof(outputs.ExtrinsicR[i]));}
                file_output_3.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_4("../ut_results/extrinsic_z_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_4.is_open()) {
            for (int i=0; i<outputs.ExtrinsicZ.size(); i++){
            for (int j=0; j<outputs.ExtrinsicZ[0].size(); j++){
                file_output_4.write((char*)(&outputs.ExtrinsicZ[i][j]), sizeof(outputs.ExtrinsicZ[i][j]));}}
                file_output_4.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

}

