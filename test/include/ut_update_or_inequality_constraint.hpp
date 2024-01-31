// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"

void ut_oric(string);


void ut_oric(string ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}

    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const int N_DEPARTMENTS = atoi(n_departments_bash);  
    const int N_PROJECTS = atoi(n_projects_bash);

    //cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_ORIC/N_TEAMS_pure.npy");
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;

    OrInequalityConstraint modelOric(N_DEPARTMENTS, N_TEAMS, N_PROJECTS);

    //cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_ORIC/eq_constraint_to_oric_m_pure.npy");
    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/eq_constraint_to_oric_m_pure.npy");
    //cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/ut_ORIC/reward_project_pure.npy");
    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/reward_project_pure.npy");
    impalib_type* eq_constraint_to_oric_m_pure = input2.data<impalib_type>();
    impalib_type* reward_project_pure = input3.data<impalib_type>();
    vector<vector<impalib_type>> eq_constraint_to_oric_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));
    vector<vector<impalib_type>> reward_project(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

    for (int project_index=0; project_index<N_PROJECTS; project_index++){
        copy ( eq_constraint_to_oric_m_pure + N_TEAMS*project_index, eq_constraint_to_oric_m_pure+N_TEAMS*(project_index+1), eq_constraint_to_oric_m[project_index].begin() );
        copy ( reward_project_pure + N_TEAMS*project_index, reward_project_pure + N_TEAMS*(project_index+1), reward_project[project_index].begin() );
    }

    if (ut_name == "Oric2TeamUpdate"){
        vector<impalib_type> oric_to_team_m(N_TEAMS, zero_value);

        modelOric.oric_to_team_update(eq_constraint_to_oric_m, oric_to_team_m);

        //fstream file_output1("../ut_results/ut_ORIC/oric_to_team_m_wrapper", ios::out | ios::binary | ios:: trunc);
        fstream file_output1("../ut_results/oric_to_team_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output1.is_open()) {
                for (int i=0; i<N_TEAMS; i++){
                    file_output1.write((char*)(&oric_to_team_m[i]), sizeof(oric_to_team_m[i]));}
                    file_output1.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
    
    }

    else if (ut_name == "Oric2ProjectEcUpdate"){

        vector<impalib_type> team_to_oric_m(N_TEAMS, zero_value);
        //cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/ut_ORIC/team_to_oric_m_pure.npy");
        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/team_to_oric_m_pure.npy");
        impalib_type* team_to_oric_m_pure = input4.data<impalib_type>();
        copy ( team_to_oric_m_pure, team_to_oric_m_pure + N_TEAMS, team_to_oric_m.begin() );

        vector<vector<impalib_type>> oric_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));
        vector<vector<impalib_type>> eq_constraint_to_project_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));   
        
        modelOric.oric_to_project_eq_constraint_update(eq_constraint_to_oric_m, team_to_oric_m, oric_to_eq_constraint_m, eq_constraint_to_project_m, reward_project);

            //fstream file_output2("../ut_results/ut_ORIC/oric_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            fstream file_output2("../ut_results/oric_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output2.is_open()) {
                for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output2.write((char*)(&oric_to_eq_constraint_m[i][j]), sizeof(oric_to_eq_constraint_m[i][j]));}}
                    file_output2.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
            
            //fstream file_output3("../ut_results/ut_ORIC/eq_constraint_to_project_m_wrapper", ios::out | ios::binary | ios:: trunc);
            fstream file_output3("../ut_results/eq_constraint_to_project_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output3.is_open()) {
                for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output3.write((char*)(&eq_constraint_to_project_m[i][j]), sizeof(eq_constraint_to_project_m[i][j]));}}
                    file_output3.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
        
    }

}