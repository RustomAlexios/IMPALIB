// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_eq_constraint(string);

void ut_eq_constraint(string ut_name){

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