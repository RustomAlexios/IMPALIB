// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_eq_constraint_kc_mwm(string);

void ut_eq_constraint_kc_mwm(string ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}
    
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_EqConstraint/ut_"+ ut_name+"/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    const int N_DEPARTMENTS = atoi(n_departments_bash);  
    const int N_PROJECTS = atoi(n_projects_bash);

    EqualityConstraintKcMwm modelEqConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS);

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_EqConstraint/ut_"+ ut_name+"/reward_project_pure.npy");
    impalib_type* reward_project_pure = input2.data<impalib_type>();
    vector<vector<impalib_type>> reward_project(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

    for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( reward_project_pure + N_TEAMS*project_index, reward_project_pure + N_TEAMS*(project_index+1), reward_project[project_index].begin() );
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_EqConstraint/ut_"+ ut_name+"/reward_team_pure.npy");
    impalib_type* reward_team_pure = input3.data<impalib_type>();
    vector<impalib_type> reward_team(N_TEAMS, zero_value);
    copy(reward_team_pure, reward_team_pure+N_TEAMS, reward_team.begin());

    if (ut_name == "TeamEc2OricUpdate"){

        vector<impalib_type> team_to_oric_m(N_TEAMS, zero_value);

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_EqConstraint/ut_"+ ut_name+"/extrinsic_output_department_pure.npy");
        impalib_type* extrinsic_output_department_pure = input4.data<impalib_type>();
        vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS, zero_value));

        for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
            copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure + N_TEAMS*(department_index+1), extrinsic_output_department[department_index].begin() );
        }

            modelEqConstraint.team_eq_constraint_to_oric_update(extrinsic_output_department, team_to_oric_m, reward_team);

            fstream file_output("../ut_results/ut_EqConstraint/ut_"+ ut_name+"/team_to_oric_m_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_TEAMS; i++){
                        file_output.write((char*)(&team_to_oric_m[i]), sizeof(team_to_oric_m[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;}
    
    }

    else if (ut_name == "ProjectEqConst2OricUpdate"){
        
        vector<vector<impalib_type>> eq_constraint_to_oric_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));
        
        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_EqConstraint/ut_"+ ut_name+"/project_to_eq_constraint_m_pure.npy");
        impalib_type* project_to_eq_constraint_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
            copy ( project_to_eq_constraint_m_pure + N_TEAMS*project_index, project_to_eq_constraint_m_pure + N_TEAMS*(project_index+1), project_to_eq_constraint_m[project_index].begin() );
        }
    
    modelEqConstraint.project_eq_constraint_to_oric_update(project_to_eq_constraint_m, eq_constraint_to_oric_m, reward_project);

    fstream file_output("../ut_results/ut_EqConstraint/ut_"+ ut_name+"/eq_constraint_to_oric_m_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output.is_open()) {
            for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&eq_constraint_to_oric_m[i][j]), sizeof(eq_constraint_to_oric_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}

    }

}

void ut_equality_constraint_tsp(string);

void ut_equality_constraint_tsp(string ut_name){

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


    cnpy::NpyArray input_alpha = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/alpha.npy");
    impalib_type* alpha_pure = input_alpha.data<impalib_type>();
    const impalib_type ALPHA = *alpha_pure;
    
    EqualityConstraintTsp modelEqualityConstraint(N_NODES, N_EDGE_VARIABLES, FILT_FLAG, ALPHA);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/edge_connections_pure.npy");
    int* edge_connections_pure = input1.data<int>();

    int num_connections = 2;

    vector<vector<int>> edge_connections(N_EDGE_VARIABLES, vector<int>(num_connections,0));
    
    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (edge_connections_pure + num_connections*edge_variable_index, edge_connections_pure+num_connections*(edge_variable_index+1), edge_connections[edge_variable_index].begin() );
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/edge_degree_constraint_cost_pure.npy");
    impalib_type* edge_degree_constraint_cost_pure = input2.data<impalib_type>();

    vector<vector<impalib_type>> edge_degree_constraint_cost(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index < N_EDGE_VARIABLES; edge_variable_index++){
    copy (edge_degree_constraint_cost_pure + N_NODES*edge_variable_index, edge_degree_constraint_cost_pure+N_NODES*(edge_variable_index+1), edge_degree_constraint_cost[edge_variable_index].begin() );
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/degree_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* degree_constraint_to_eq_constraint_m_pure = input3.data<impalib_type>();

    vector<vector<impalib_type>> degree_constraint_to_eq_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (degree_constraint_to_eq_constraint_m_pure + N_NODES*edge_variable_index, degree_constraint_to_eq_constraint_m_pure+N_NODES*(edge_variable_index+1), degree_constraint_to_eq_constraint_m[edge_variable_index].begin() );
    }

    if (ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate"){
        
        vector<vector<impalib_type>> edge_ec_to_degree_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

        modelEqualityConstraint.edge_ec_to_degree_constraint_relaxed_graph_update(edge_connections, edge_degree_constraint_cost, degree_constraint_to_eq_constraint_m, edge_ec_to_degree_constraint_m);
        
        fstream file_output("../ut_results/ut_EqualityConstraint/ut_"+ ut_name+"/edge_ec_to_degree_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                for (int j=0; j<N_NODES; j++){
                    file_output.write((char*)(&edge_ec_to_degree_constraint_m[i][j]), sizeof(edge_ec_to_degree_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
    }

    else if (ut_name == "EdgeEc2SubtourConstraintsUpdate"){

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/delta_S_indices_list_sizes_pure.npy");
        const int* delta_S_indices_list_sizes_pure = input4.data<int>();

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/delta_S_indices_list_pure.npy");
        const int* delta_S_indices_list_pure = input5.data<int>();

        vector<vector<int>> delta_S_indices_list;

        int idx=0;
        for (int subtour_index=0; subtour_index < N_SUBTOURS; subtour_index++){
            int size = delta_S_indices_list_sizes_pure[subtour_index];
            vector<int> temp_vector(delta_S_indices_list_pure + idx, delta_S_indices_list_pure + idx + size);
            delta_S_indices_list.push_back(temp_vector);
            idx += size;
        }
        
        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/cost_edge_variable_pure.npy");
        
        impalib_type* cost_edge_variable_pure = input6.data<impalib_type>();
        vector<impalib_type> cost_edge_variable(cost_edge_variable_pure, cost_edge_variable_pure + N_EDGE_VARIABLES);
        
        cnpy::NpyArray input7 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input7.data<impalib_type>();

        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy (subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        vector<vector<impalib_type>> edge_ec_to_subtour_constraints_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES, zero_value));

        edge_ec_to_subtour_constraints_m = modelEqualityConstraint.edge_ec_to_subtour_constraints_update(delta_S_indices_list, cost_edge_variable, degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m, edge_connections);
        

        fstream file_output("../ut_results/ut_EqualityConstraint/ut_"+ ut_name+"/edge_ec_to_subtour_constraints_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_SUBTOURS; i++){
                for (int j=0; j<N_EDGE_VARIABLES; j++){
                    file_output.write((char*)(&edge_ec_to_subtour_constraints_m[i][j]), sizeof(edge_ec_to_subtour_constraints_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
        
    }

    else if (ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate"){

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/ut_EqualityConstraint/ut_"+ ut_name+"/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input8.data<impalib_type>();

        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy (subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        vector<vector<impalib_type>> edge_ec_to_degree_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES, zero_value));

        modelEqualityConstraint.edge_ec_to_degree_constraint_augmented_graph_update(degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m, edge_connections, edge_degree_constraint_cost, edge_ec_to_degree_constraint_m);
        
        fstream file_output("../ut_results/ut_EqualityConstraint/ut_"+ ut_name+"/edge_ec_to_degree_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                for (int j=0; j<N_NODES; j++){
                    file_output.write((char*)(&edge_ec_to_degree_constraint_m[i][j]), sizeof(edge_ec_to_degree_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
    
    }
}
