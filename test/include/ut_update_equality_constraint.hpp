// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_eq_constraint_kc_mwm(string&);

void ut_eq_constraint_kc_mwm(string& ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}
    
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    const int N_DEPARTMENTS = atoi(n_departments_bash);  
    const int N_PROJECTS = atoi(n_projects_bash);

    EqualityConstraint modelEqConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS);

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    impalib_type* reward_project_pure = input2.data<impalib_type>();
    vector<vector<impalib_type>> reward_project(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

    for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( reward_project_pure + N_TEAMS*project_index, reward_project_pure + N_TEAMS*(project_index+1), reward_project[project_index].begin() );
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    impalib_type* reward_team_pure = input3.data<impalib_type>();
    vector<impalib_type> reward_team(N_TEAMS, zero_value);
    copy(reward_team_pure, reward_team_pure+N_TEAMS, reward_team.begin());

    if (ut_name == "TeamEc2OricUpdate"){

        vector<impalib_type> team_to_oric_m(N_TEAMS, zero_value);

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/extrinsic_output_department_pure.npy");
        impalib_type* extrinsic_output_department_pure = input4.data<impalib_type>();
        vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS, zero_value));

        for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
            copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure + N_TEAMS*(department_index+1), extrinsic_output_department[department_index].begin() );
        }

            modelEqConstraint.team_eq_constraint_to_oric_update(extrinsic_output_department, team_to_oric_m, reward_team);

            fstream file_output("../ut_results/team_to_oric_m_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_TEAMS; i++){
                        file_output.write((char*)(&team_to_oric_m[i]), sizeof(team_to_oric_m[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

    else if (ut_name == "ProjectEqConst2OricUpdate"){
        
        vector<vector<impalib_type>> eq_constraint_to_oric_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));
        
        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/project_to_eq_constraint_m_pure.npy");
        impalib_type* project_to_eq_constraint_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
            copy ( project_to_eq_constraint_m_pure + N_TEAMS*project_index, project_to_eq_constraint_m_pure + N_TEAMS*(project_index+1), project_to_eq_constraint_m[project_index].begin() );
        }
    
    modelEqConstraint.project_eq_constraint_to_oric_update(project_to_eq_constraint_m, eq_constraint_to_oric_m, reward_project);

    fstream file_output("../ut_results/eq_constraint_to_oric_m_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output.is_open()) {
            for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&eq_constraint_to_oric_m[i][j]), sizeof(eq_constraint_to_oric_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}

    }

}

void ut_equality_constraint_tsp(string&);

void ut_equality_constraint_tsp(string& ut_name){

    const char *n_nodes_bash=getenv("N_NODES");
    if(n_nodes_bash == NULL)
    {cout << "n_nodes_bash not available\n";}

    const char *n_subtours_bash=getenv("N_SUBTOURS");
    if(n_subtours_bash == NULL)
    {cout << "n_subtours_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const int N_NODES = atoi(n_nodes_bash);  
    const int N_SUBTOURS = atoi(n_subtours_bash);
    const int N_EDGE_VARIABLES = N_NODES*N_NODES-N_NODES;
    // const bool FILT_FLAG(filt_flag_bash);
    const bool FILT_FLAG = (filt_flag_bash != NULL && std::string(filt_flag_bash) == "1");

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;
    
    EqualityConstraint modelEqualityConstraint(N_NODES, N_EDGE_VARIABLES, FILT_FLAG, ALPHA);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/edge_connections_pure.npy");
    int* edge_connections_pure = input1.data<int>();

    int num_connections = 2;

    vector<vector<int>> edge_connections(N_EDGE_VARIABLES, vector<int>(num_connections,0));
    
    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (edge_connections_pure + num_connections*edge_variable_index, edge_connections_pure+num_connections*(edge_variable_index+1), edge_connections[edge_variable_index].begin() );
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/edge_degree_constraint_cost_pure.npy");
    impalib_type* edge_degree_constraint_cost_pure = input2.data<impalib_type>();

    vector<vector<impalib_type>> edge_degree_constraint_cost(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index < N_EDGE_VARIABLES; edge_variable_index++){
    copy (edge_degree_constraint_cost_pure + N_NODES*edge_variable_index, edge_degree_constraint_cost_pure+N_NODES*(edge_variable_index+1), edge_degree_constraint_cost[edge_variable_index].begin() );
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/degree_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* degree_constraint_to_eq_constraint_m_pure = input3.data<impalib_type>();

    vector<vector<impalib_type>> degree_constraint_to_eq_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (degree_constraint_to_eq_constraint_m_pure + N_NODES*edge_variable_index, degree_constraint_to_eq_constraint_m_pure+N_NODES*(edge_variable_index+1), degree_constraint_to_eq_constraint_m[edge_variable_index].begin() );
    }

    if (ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate"){
        
        vector<vector<impalib_type>> edge_ec_to_degree_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

        modelEqualityConstraint.edge_ec_to_degree_constraint_relaxed_graph_update(edge_connections, edge_degree_constraint_cost, degree_constraint_to_eq_constraint_m, edge_ec_to_degree_constraint_m);
        
        fstream file_output("../ut_results/edge_ec_to_degree_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                for (int j=0; j<N_NODES; j++){
                    file_output.write((char*)(&edge_ec_to_degree_constraint_m[i][j]), sizeof(edge_ec_to_degree_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
    }

    else if (ut_name == "EdgeEc2SubtourConstraintsUpdate"){

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/delta_S_indices_list_sizes_pure.npy");
        const int* delta_S_indices_list_sizes_pure = input4.data<int>();

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/delta_S_indices_list_pure.npy");
        const int* delta_S_indices_list_pure = input5.data<int>();

        vector<vector<int>> delta_S_indices_list;

        int idx=0;
        for (int subtour_index=0; subtour_index < N_SUBTOURS; subtour_index++){
            int size = delta_S_indices_list_sizes_pure[subtour_index];
            vector<int> temp_vector(delta_S_indices_list_pure + idx, delta_S_indices_list_pure + idx + size);
            delta_S_indices_list.push_back(temp_vector);
            idx += size;
        }
        
        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/cost_edge_variable_pure.npy");
        
        impalib_type* cost_edge_variable_pure = input6.data<impalib_type>();
        vector<impalib_type> cost_edge_variable(cost_edge_variable_pure, cost_edge_variable_pure + N_EDGE_VARIABLES);
        
        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input7.data<impalib_type>();

        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy (subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        vector<vector<impalib_type>> edge_ec_to_subtour_constraints_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES, zero_value));

        edge_ec_to_subtour_constraints_m = modelEqualityConstraint.edge_ec_to_subtour_constraints_update(delta_S_indices_list, cost_edge_variable, degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m, edge_connections);
        
        fstream file_output("../ut_results/edge_ec_to_subtour_constraints_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_SUBTOURS; i++){
                for (int j=0; j<N_EDGE_VARIABLES; j++){
                    file_output.write((char*)(&edge_ec_to_subtour_constraints_m[i][j]), sizeof(edge_ec_to_subtour_constraints_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
        
    }

    else if (ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate"){

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input8.data<impalib_type>();

        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy (subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        vector<vector<impalib_type>> edge_ec_to_degree_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

        modelEqualityConstraint.edge_ec_to_degree_constraint_augmented_graph_update(degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m, edge_connections, edge_degree_constraint_cost, edge_ec_to_degree_constraint_m);
        
        fstream file_output("../ut_results/edge_ec_to_degree_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                for (int j=0; j<N_NODES; j++){
                    file_output.write((char*)(&edge_ec_to_degree_constraint_m[i][j]), sizeof(edge_ec_to_degree_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
    
    }
}

void ut_equality_constraint_ksat(string&);

void ut_equality_constraint_ksat(string& ut_name){


    const char *n_variables_bash=getenv("NUM_VARIABLES");
    if(n_variables_bash == NULL)
    {cout << "n_variables_bash not available\n";}

    const char *k_variable_bash=getenv("K_VARIABLE");
    if(k_variable_bash == NULL)
    {cout << "k_variable_bash not available\n";}

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}

    const int NUM_VARIABLES = atoi(n_variables_bash);  
    const int K_VARIABLE = atoi(k_variable_bash);
    const bool FILT_FLAG(filt_flag_bash);

    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;

    cnpy::NpyArray input_num_constraints = cnpy::npy_load("../ut_inputs/num_constraints.npy");
    int* num_constraints_pure = input_num_constraints.data<int>();
    const int NUM_CONSTRAINTS = *num_constraints_pure;

    EqualityConstraint modelEqualityConstraint(NUM_VARIABLES, NUM_CONSTRAINTS, K_VARIABLE, FILT_FLAG, ALPHA);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/size_used_variables_pure.npy");
    int* size_used_variables_pure = input1.data<int>();
    int size_used_variables = *size_used_variables_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/used_variables_pure.npy");
    int* used_variables_pure = input2.data<int>(); 

    vector<int> used_variables;

    copy(used_variables_pure, used_variables_pure + size_used_variables, back_inserter(used_variables));

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/sizes_of_variables_connections_pure.npy");
    int* sizes_of_variables_connections_pure = input3.data<int>(); 
    vector<int> sizes_of_variables_connections;
    copy(sizes_of_variables_connections_pure, sizes_of_variables_connections_pure + NUM_VARIABLES, back_inserter(sizes_of_variables_connections));

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/variables_connections_pure.npy");
    int* variables_connections_pure = input4.data<int>(); 
    
    vector<vector<int>> variables_connections;

    int connections_size_old = 0;

    for (int variable_index = 0; variable_index < NUM_VARIABLES; variable_index++){
        
        if (find(used_variables.begin(), used_variables.end(), variable_index) != used_variables.end()) {
            
            int connections_size = sizes_of_variables_connections[variable_index];
            
            variables_connections.push_back(vector<int>(connections_size, 0));
            
            copy(variables_connections_pure + connections_size_old,
             variables_connections_pure + connections_size_old + connections_size,
             variables_connections[variable_index].begin());

            connections_size_old += connections_size;

        } else {
            variables_connections.push_back(vector<int>());
        }
    }

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/incoming_metrics_cost_pure.npy");
    impalib_type* incoming_metrics_cost_pure = input5.data<impalib_type>(); 
    vector<impalib_type> incoming_metrics_cost;
    copy(incoming_metrics_cost_pure, incoming_metrics_cost_pure + NUM_VARIABLES, back_inserter(incoming_metrics_cost));

    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ksat_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* ksat_constraint_to_eq_constraint_m_pure = input6.data<impalib_type>(); 
    vector<vector<impalib_type>> ksat_constraint_to_eq_constraint_m(NUM_CONSTRAINTS, vector<impalib_type>(NUM_VARIABLES,zero_value));

    for (int constraint_index=0; constraint_index<NUM_CONSTRAINTS; constraint_index++){
    copy (ksat_constraint_to_eq_constraint_m_pure + NUM_VARIABLES*constraint_index, ksat_constraint_to_eq_constraint_m_pure+NUM_VARIABLES*(constraint_index+1), ksat_constraint_to_eq_constraint_m[constraint_index].begin() );
    }

    cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/constraints_connections_pure.npy");
    int* constraints_connections_pure = input7.data<int>(); 
    vector<vector<int>> constraints_connections(NUM_CONSTRAINTS, vector<int>(K_VARIABLE,0));

    for (int constraint_index=0; constraint_index<NUM_CONSTRAINTS; constraint_index++){
    copy (constraints_connections_pure + K_VARIABLE*constraint_index, constraints_connections_pure+K_VARIABLE*(constraint_index+1), constraints_connections[constraint_index].begin() );
    }

    if (ut_name == "VariableEc2KsatConstraintUpdate"){
        
        vector<vector<impalib_type>> variable_ec_to_ksat_constraint_m(NUM_CONSTRAINTS, vector<impalib_type>(NUM_VARIABLES, zero_value));

        modelEqualityConstraint.variable_ec_to_ksat_constraint_update(ksat_constraint_to_eq_constraint_m, variable_ec_to_ksat_constraint_m, used_variables, incoming_metrics_cost, variables_connections);
        
        fstream file_output("../ut_results/variable_ec_to_ksat_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<NUM_CONSTRAINTS; i++){
                for (int j=0; j<NUM_VARIABLES; j++){
                    file_output.write((char*)(&variable_ec_to_ksat_constraint_m[i][j]), sizeof(variable_ec_to_ksat_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
    }

}

void ut_equality_constraint_mobarp(string&);

void ut_equality_constraint_mobarp(string& ut_name){
    
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

    const impalib_type ALPHA = 0.5;
    const bool FILTERING_FLAG = false;

    EqualityConstraint model_equality_constraint(NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, 
                                                    NUM_MOBILE_TX_LOCS, ALPHA, FILTERING_FLAG, EXCLUDE_CAP_FLAG);

    //yes
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/fixed_x_costs.npy");
    impalib_type* fixed_x_costs_pure = input1.data<impalib_type>();
    vector<vector<vector<impalib_type>>> FixedTxCosts;

    for (int i=0; i< NUM_FIXED_TX; i++){
        FixedTxCosts.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        for (int j=0; j< NUM_BANDS; j++){
            copy(fixed_x_costs_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, fixed_x_costs_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, FixedTxCosts[i][j].begin());
        }
    }

    //yes
    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/mobile_x_costs.npy");
    impalib_type* mobile_x_costs_pure = input2.data<impalib_type>();
    vector<vector<vector<impalib_type>>> MobileTxCosts;

    for (int i=0; i< NUM_MOBILE_TX; i++){
        MobileTxCosts.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        for (int j=0; j< NUM_BANDS; j++){
            copy(mobile_x_costs_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, mobile_x_costs_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, MobileTxCosts[i][j].begin());
        }
    }

    //yes
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


    

    if (ut_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate"){
        
        //yes
        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/fixed_capac_const_to_fixed_x_eq_const_m_pure.npy");
        impalib_type* fixed_capac_const_to_fixed_x_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<vector<impalib_type>>> fixed_capac_const_to_fixed_x_eq_const_m;

        for (int i=0; i< NUM_FIXED_TX; i++){
            fixed_capac_const_to_fixed_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(fixed_capac_const_to_fixed_x_eq_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, fixed_capac_const_to_fixed_x_eq_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, fixed_capac_const_to_fixed_x_eq_const_m[i][j].begin());
            }
        }

        //yes
        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/mobile_capac_const_to_mobile_x_eq_const_m_pure.npy");
        impalib_type* mobile_capac_const_to_mobile_x_eq_const_m_pure = input7.data<impalib_type>();
        vector<vector<vector<impalib_type>>> mobile_capac_const_to_mobile_x_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            mobile_capac_const_to_mobile_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(mobile_capac_const_to_mobile_x_eq_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, mobile_capac_const_to_mobile_x_eq_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, mobile_capac_const_to_mobile_x_eq_const_m[i][j].begin());
            }
        }

        //yes
        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_mobile_x_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_mobile_x_eq_const_m_pure = input8.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_mobile_x_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_mobile_x_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_mobile_x_eq_const_m[i].begin());
        }

        //yes
        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_fixed_x_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_fixed_x_eq_const_m_pure = input9.data<impalib_type>();
        vector<vector<vector<impalib_type>>> set_cover_ineq_const_to_fixed_x_eq_const_m;

        for (int i=0; i< NUM_TIME_STEPS; i++){
            set_cover_ineq_const_to_fixed_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0)));
            for (int j=0; j< NUM_RX_LOCS; j++){
                copy(set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * j + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * (j + 1) + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m[i][j].begin());
            }
        }

        //yes
        cnpy::NpyArray input10 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input10.data<int>();
        vector<vector<int>> conx_mob_tx_per_num_mob_tx_locs;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_mob_tx_per_num_mob_tx_locs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
            copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), conx_mob_tx_per_num_mob_tx_locs[i].begin());
        }

        //yes
        cnpy::NpyArray input11 = cnpy::npy_load("../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy");
        int* conx_fixed_tx_per_num_rx_locs_pure = input11.data<int>();
        vector<vector<int>> conx_fixed_tx_per_num_rx_locs;

        for (int i=0; i< NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_fixed_tx_per_num_rx_locs.push_back(vector<int>(NUM_RX_LOCS, 0));
            copy(conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * i, conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * (i + 1), conx_fixed_tx_per_num_rx_locs[i].begin());
        }

        vector<vector<impalib_type>> FixedXEqConst2SetCoverConstM(NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_RX_LOCS, 0));
        vector<vector<impalib_type>> MobileXEqConst2AuxiliaryConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        model_equality_constraint.x_eq_const_to_auxiliary_and_set_cover_const_update(fixed_capac_const_to_fixed_x_eq_const_m, mobile_capac_const_to_mobile_x_eq_const_m, auxiliary_const_to_mobile_x_eq_const_m, 
                                                set_cover_ineq_const_to_fixed_x_eq_const_m, FixedTxCosts,
                                                    MobileTxCosts, conx_mob_tx_per_num_mob_tx_locs, MobileXEqConst2AuxiliaryConstM, conx_fixed_tx_per_num_rx_locs, FixedXEqConst2SetCoverConstM);


        fstream file_output_1("../ut_results/FixedXEqConst2SetCoverConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<FixedXEqConst2SetCoverConstM.size(); i++){
            for (int j=0; j<FixedXEqConst2SetCoverConstM[0].size(); j++){
                file_output_1.write((char*)(&FixedXEqConst2SetCoverConstM[i][j]), sizeof(FixedXEqConst2SetCoverConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/MobileXEqConst2AuxiliaryConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<MobileXEqConst2AuxiliaryConstM.size(); i++){
            for (int j=0; j<MobileXEqConst2AuxiliaryConstM[0].size(); j++){
                file_output_2.write((char*)(&MobileXEqConst2AuxiliaryConstM[i][j]), sizeof(MobileXEqConst2AuxiliaryConstM[i][j]));}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

    else if (ut_name == "ReqConstActivation"){
    
        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_r_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_r_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_r_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_r_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_r_eq_const_m[i].begin());
        }

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/mobile_loc_eq_const_to_r_eq_const_m_pure.npy");
        impalib_type* mobile_loc_eq_const_to_r_eq_const_m_pure = input7.data<impalib_type>();
        vector<vector<impalib_type>> mobile_loc_eq_const_to_r_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            mobile_loc_eq_const_to_r_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(mobile_loc_eq_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, mobile_loc_eq_const_to_r_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), mobile_loc_eq_const_to_r_eq_const_m[i].begin());
        }

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input8.data<int>();
        vector<vector<int>> conx_mob_tx_per_num_mob_tx_locs;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_mob_tx_per_num_mob_tx_locs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
            copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), conx_mob_tx_per_num_mob_tx_locs[i].begin());
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

        vector<vector<impalib_type>> REqConst2AuxiliaryConstM(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_TIME_STEPS, 0));
        vector<vector<impalib_type>> REqConst2MobileLocEqConstM(NUM_MOBILE_TX, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        model_equality_constraint.r_eq_const_activation(auxiliary_const_to_r_eq_const_m, mobile_loc_eq_const_to_r_eq_const_m, REqConst2AuxiliaryConstM,
                                                    REqConst2MobileLocEqConstM, ConxMobTxR, RCosts);

        fstream file_output_1("../ut_results/REqConst2AuxiliaryConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<REqConst2AuxiliaryConstM.size(); i++){
            for (int j=0; j<REqConst2AuxiliaryConstM[0].size(); j++){
                file_output_1.write((char*)(&REqConst2AuxiliaryConstM[i][j]), sizeof(REqConst2AuxiliaryConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}


        fstream file_output_2("../ut_results/REqConst2MobileLocEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<REqConst2MobileLocEqConstM.size(); i++){
            for (int j=0; j<REqConst2MobileLocEqConstM[0].size(); j++){
                file_output_2.write((char*)(&REqConst2MobileLocEqConstM[i][j]), sizeof(REqConst2MobileLocEqConstM[i][j]));}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

    else if (ut_name == "ZEqConst2AuxiliaryConstUpdate"){

        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_z_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_z_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<vector<vector<impalib_type>>>> set_cover_ineq_const_to_z_eq_const_m;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            set_cover_ineq_const_to_z_eq_const_m.push_back(vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m[k][l][n].begin());
            }
            }
        }

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/conx_mob_tx_rx.npy");
        int* conx_mob_tx_rx_pure = input7.data<int>();
        vector<vector<vector<vector<int>>>> conx_mob_tx_rx;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            conx_mob_tx_rx.push_back(vector<vector<vector<int>>>(NUM_RX_LOCS, vector<vector<int>>(NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx[k][l][n].begin());
            }
            }
        }

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input8.data<int>();
        vector<vector<int>> conx_mob_tx_per_num_mob_tx_locs;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_mob_tx_per_num_mob_tx_locs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
            copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), conx_mob_tx_per_num_mob_tx_locs[i].begin());
        }


        vector<vector<impalib_type>> ZEqConst2AuxiliaryConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/z_costs.npy");
        impalib_type* z_costs_pure = input9.data<impalib_type>();
        vector<vector<impalib_type>> ZCosts;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            ZCosts.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(z_costs_pure + NUM_MOBILE_TX_LOCS * i, z_costs_pure + NUM_MOBILE_TX_LOCS * (i + 1), ZCosts[i].begin());
        }

        model_equality_constraint.z_eq_const_to_auxiliary_const_update(set_cover_ineq_const_to_z_eq_const_m, conx_mob_tx_rx, conx_mob_tx_per_num_mob_tx_locs,
                        ZEqConst2AuxiliaryConstM, ZCosts);

        
        fstream file_output_1("../ut_results/ZEqConst2AuxiliaryConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<ZEqConst2AuxiliaryConstM.size(); i++){
            for (int j=0; j<ZEqConst2AuxiliaryConstM[0].size(); j++){
                file_output_1.write((char*)(&ZEqConst2AuxiliaryConstM[i][j]), sizeof(ZEqConst2AuxiliaryConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }


    else if (ut_name == "XEqConstActivation"){


        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_mobile_x_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_mobile_x_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_mobile_x_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_mobile_x_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_mobile_x_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_mobile_x_eq_const_m[i].begin());
        }

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_fixed_x_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_fixed_x_eq_const_m_pure = input7.data<impalib_type>();
        vector<vector<vector<impalib_type>>> set_cover_ineq_const_to_fixed_x_eq_const_m;

        for (int i=0; i< NUM_TIME_STEPS; i++){
            set_cover_ineq_const_to_fixed_x_eq_const_m.push_back(vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0)));
            for (int j=0; j< NUM_RX_LOCS; j++){
                copy(set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * j + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m_pure + NUM_FIXED_TX*NUM_BANDS * (j + 1) + NUM_RX_LOCS*NUM_FIXED_TX*NUM_BANDS*i, set_cover_ineq_const_to_fixed_x_eq_const_m[i][j].begin());
            }
        }

        vector<vector<vector<impalib_type>>> MobileXEqConst2MobileCapacConstM(NUM_MOBILE_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        vector<vector<vector<impalib_type>>> FixedXEqConst2FixedCapacConstM(NUM_FIXED_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
        int* conx_mob_tx_per_num_mob_tx_locs_pure = input8.data<int>();
        vector<vector<int>> conx_mob_tx_per_num_mob_tx_locs;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_mob_tx_per_num_mob_tx_locs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
            copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), conx_mob_tx_per_num_mob_tx_locs[i].begin());
        }

        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy");
        int* conx_fixed_tx_per_num_rx_locs_pure = input9.data<int>();
        vector<vector<int>> conx_fixed_tx_per_num_rx_locs;

        for (int i=0; i< NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_fixed_tx_per_num_rx_locs.push_back(vector<int>(NUM_RX_LOCS, 0));
            copy(conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * i, conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * (i + 1), conx_fixed_tx_per_num_rx_locs[i].begin());
        }
        model_equality_constraint.x_eq_const_activation(auxiliary_const_to_mobile_x_eq_const_m, set_cover_ineq_const_to_fixed_x_eq_const_m,
                            MobileXEqConst2MobileCapacConstM, FixedXEqConst2FixedCapacConstM,
                            conx_mob_tx_per_num_mob_tx_locs, conx_fixed_tx_per_num_rx_locs, MobileTxCosts, FixedTxCosts);
    
        fstream file_output_1("../ut_results/MobileXEqConst2MobileCapacConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<NUM_MOBILE_TX; i++){
            for (int j=0; j<NUM_BANDS; j++){
                for (int k=0; k<NUM_TIME_STEPS; k++){
                file_output_1.write((char*)(&MobileXEqConst2MobileCapacConstM[i][j][k]), sizeof(MobileXEqConst2MobileCapacConstM[i][j][k]));}}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/FixedXEqConst2FixedCapacConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<NUM_FIXED_TX; i++){
            for (int j=0; j<NUM_BANDS; j++){
                for (int k=0; k<NUM_TIME_STEPS; k++){
                file_output_2.write((char*)(&FixedXEqConst2FixedCapacConstM[i][j][k]), sizeof(FixedXEqConst2FixedCapacConstM[i][j][k]));}}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";} 

    }

    else if (ut_name == "ZEqConst2SetCoverIneqConstUpdate"){

        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/auxiliary_const_to_z_eq_const_m_pure.npy");
        impalib_type* auxiliary_const_to_z_eq_const_m_pure = input6.data<impalib_type>();
        vector<vector<impalib_type>> auxiliary_const_to_z_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            auxiliary_const_to_z_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(auxiliary_const_to_z_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, auxiliary_const_to_z_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), auxiliary_const_to_z_eq_const_m[i].begin());
        }

        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/set_cover_ineq_const_to_z_eq_const_m_pure.npy");
        impalib_type* set_cover_ineq_const_to_z_eq_const_m_pure = input7.data<impalib_type>();
        vector<vector<vector<vector<impalib_type>>>> set_cover_ineq_const_to_z_eq_const_m;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            set_cover_ineq_const_to_z_eq_const_m.push_back(vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, set_cover_ineq_const_to_z_eq_const_m[k][l][n].begin());
            }
            }
        }

        vector<vector<impalib_type>> ZEqConst2SetCoverIneqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS*NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_RX_LOCS, 0));

        vector<vector<vector<vector<int>>>> TransposedConxMobTxRx(NUM_BANDS*NUM_MOBILE_TX, vector<vector<vector<int>>>(NUM_RX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_MOBILE_TX_LOCS, 0))));

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/conx_mob_tx_rx.npy");
        int* conx_mob_tx_rx_pure = input8.data<int>();
        vector<vector<vector<vector<int>>>> conx_mob_tx_rx;

        for (int k=0; k< NUM_TIME_STEPS; k++){
            conx_mob_tx_rx.push_back(vector<vector<vector<int>>>(NUM_RX_LOCS, vector<vector<int>>(NUM_MOBILE_TX_LOCS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))));
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    copy(conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*n + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx_pure + NUM_BANDS*NUM_MOBILE_TX*(n+1) + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*l + NUM_BANDS*NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS*NUM_RX_LOCS*k, conx_mob_tx_rx[k][l][n].begin());
            }
            }
        }

        for (int k=0; k< NUM_TIME_STEPS; k++){
        for (int l=0; l<NUM_RX_LOCS; l++){
            for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                for (int j_i=0; j_i <NUM_BANDS*NUM_MOBILE_TX; j_i++){
                    TransposedConxMobTxRx[j_i][l][k][n] = conx_mob_tx_rx[k][l][n][j_i];
                }
            }
        }
        }

        cnpy::NpyArray input9 = cnpy::npy_load("../ut_inputs/z_costs.npy");
        impalib_type* z_costs_pure = input9.data<impalib_type>();
        vector<vector<impalib_type>> ZCosts;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            ZCosts.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(z_costs_pure + NUM_MOBILE_TX_LOCS * i, z_costs_pure + NUM_MOBILE_TX_LOCS * (i + 1), ZCosts[i].begin());
        }

        model_equality_constraint.z_eq_const_to_set_cover_ineq_const_update(auxiliary_const_to_z_eq_const_m, 
                                        set_cover_ineq_const_to_z_eq_const_m, ZEqConst2SetCoverIneqConstM,
                                        TransposedConxMobTxRx, ZCosts);

        fstream file_output_1("../ut_results/ZEqConst2SetCoverIneqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<ZEqConst2SetCoverIneqConstM.size(); i++){
            for (int j=0; j<ZEqConst2SetCoverIneqConstM[0].size(); j++){
                file_output_1.write((char*)(&ZEqConst2SetCoverIneqConstM[i][j]), sizeof(ZEqConst2SetCoverIneqConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

}