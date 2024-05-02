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
    const bool FILT_FLAG(filt_flag_bash);

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