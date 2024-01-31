// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_input_output_kc_mwm(string);

void ut_input_output_kc_mwm(string ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const int N_DEPARTMENT = atoi(n_departments_bash);  
    const int N_PROJECTS = atoi(n_projects_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/reward_project_pure.npy");
    impalib_type* reward_project_pure = input2.data<impalib_type>();
    vector<vector<impalib_type>> reward_project(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int project_index=0; project_index<N_PROJECTS; project_index++){
    copy ( reward_project_pure + N_TEAMS*project_index, reward_project_pure+N_TEAMS*(project_index+1), reward_project[project_index].begin() );
    }

    OutputsKcMwm outputs(N_DEPARTMENT, N_TEAMS, N_PROJECTS);

    if (ut_name == "ExtrinsicOutputTeamUpdate"){
        
        cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/extrinsic_output_department_pure.npy");
        impalib_type* extrinsic_output_department_pure = input3.data<impalib_type>();
        vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENT, vector<impalib_type>(N_TEAMS,zero_value));

        for (int department_index=0; department_index<N_DEPARTMENT; department_index++){
        copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure+N_TEAMS*(department_index+1), extrinsic_output_department[department_index].begin() );
        }

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/oric_to_team_m_pure.npy");
        impalib_type* oric_to_team_m_pure = input4.data<impalib_type>();
        vector<impalib_type> oric_to_team_m(N_TEAMS,zero_value);
        copy(oric_to_team_m_pure, oric_to_team_m_pure + N_TEAMS, oric_to_team_m.begin());

        
        outputs.extrinsic_output_team_update(extrinsic_output_department, oric_to_team_m);

        fstream file_output("../ut_results/ut_InputOutput/ut_"+ ut_name+"/extrinsic_output_team_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_TEAMS; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputTeam[i]), sizeof(outputs.ExtrinsicOutputTeam[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;}
    
    }

    if (ut_name == "IntrinsicOutMwmUpdate"){
        
        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/oric_to_eq_constraint_m_pure.npy");
        impalib_type* oric_to_eq_constraint_m_pure = input4.data<impalib_type>();
        vector<vector<impalib_type>> oric_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( oric_to_eq_constraint_m_pure + N_TEAMS*project_index, oric_to_eq_constraint_m_pure+N_TEAMS*(project_index+1), oric_to_eq_constraint_m[project_index].begin() );
        }

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/project_to_eq_constraint_m_pure.npy");
        impalib_type* project_to_eq_constraint_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS,zero_value));

        for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( project_to_eq_constraint_m_pure + N_TEAMS*project_index, project_to_eq_constraint_m_pure+N_TEAMS*(project_index+1), project_to_eq_constraint_m[project_index].begin() );
        }

        outputs.intrinsic_out_mwm_update(oric_to_eq_constraint_m, project_to_eq_constraint_m, reward_project);

        fstream file_output("../ut_results/ut_InputOutput/ut_"+ ut_name+"/intrinsic_out_mwm_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_PROJECTS*N_TEAMS; i++){
                        file_output.write((char*)(&outputs.IntrinsicOutMwm[i]), sizeof(outputs.IntrinsicOutMwm[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;}

    }

}

void ut_input_output_tsp(string);

void ut_input_output_tsp(string ut_name){

    const char *n_nodes_bash=getenv("N_NODES");
    if(n_nodes_bash == NULL)
    {cout << "n_nodes_bash not available\n";}

    const char *n_subtours_bash=getenv("N_SUBTOURS");
    if(n_subtours_bash == NULL)
    {cout << "n_subtours_bash not available\n";}

    const int N_NODES = atoi(n_nodes_bash);  
    const int N_EDGE_VARIABLES = N_NODES*N_NODES-N_NODES;
    const int N_SUBTOURS = atoi(n_subtours_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/degree_constraint_to_eq_constraint_m_pure.npy");
    impalib_type* degree_constraint_to_eq_constraint_m_pure = input1.data<impalib_type>();

    vector<vector<impalib_type>> degree_constraint_to_eq_constraint_m(N_EDGE_VARIABLES, vector<impalib_type>(N_NODES,zero_value));

    for (int edge_variable_index=0; edge_variable_index<N_EDGE_VARIABLES; edge_variable_index++){
    copy (degree_constraint_to_eq_constraint_m_pure + N_NODES*edge_variable_index, degree_constraint_to_eq_constraint_m_pure+N_NODES*(edge_variable_index+1), degree_constraint_to_eq_constraint_m[edge_variable_index].begin() );
    }

    OutputsTsp outputs(N_NODES, N_EDGE_VARIABLES);

    if (ut_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate"){
        
        outputs.extrinsic_output_edge_ec_relaxed_graph_update(degree_constraint_to_eq_constraint_m);
        
        fstream file_output("../ut_results/ut_InputOutput/ut_"+ ut_name+"/extrinsic_output_edge_ec_relaxed_graph_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_EDGE_VARIABLES; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputEdgeEc[i]), sizeof(outputs.ExtrinsicOutputEdgeEc[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;}
    
    }

    if (ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"){
        
        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_InputOutput/ut_"+ ut_name+"/subtour_constraints_to_edge_ec_m_pure.npy");
        impalib_type* subtour_constraints_to_edge_ec_m_pure = input2.data<impalib_type>();
        vector<vector<impalib_type>> subtour_constraints_to_edge_ec_m(N_SUBTOURS, vector<impalib_type>(N_EDGE_VARIABLES,zero_value));

        for (int subtour_index=0; subtour_index<N_SUBTOURS; subtour_index++){
        copy ( subtour_constraints_to_edge_ec_m_pure + N_EDGE_VARIABLES*subtour_index, subtour_constraints_to_edge_ec_m_pure+N_EDGE_VARIABLES*(subtour_index+1), subtour_constraints_to_edge_ec_m[subtour_index].begin() );
        }

        outputs.extrinsic_output_edge_ec_augmented_graph_update(degree_constraint_to_eq_constraint_m, subtour_constraints_to_edge_ec_m);

        fstream file_output("../ut_results/ut_InputOutput/ut_"+ ut_name+"/extrinsic_output_edge_ec_augmented_graph_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_EDGE_VARIABLES; i++){
                        file_output.write((char*)(&outputs.ExtrinsicOutputEdgeEc[i]), sizeof(outputs.ExtrinsicOutputEdgeEc[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;}
    }

}


