    // Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_project_ineq_constraint(string);

void ut_project_ineq_constraint(string ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}
    
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/ut_projectIneqConstraint/ut_"+ ut_name+"/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;
    
    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const int N_DEPARTMENTS = atoi(n_departments_bash); 
    const int N_PROJECTS = atoi(n_projects_bash);

    InequalityConstraint projectIneqConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS);

    if (ut_name == "ProjectInequalityConstraintUpdate"){

            cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/ut_projectIneqConstraint/ut_"+ ut_name+"/eq_constraint_to_project_m_pure.npy");
            impalib_type* eq_constraint_to_project_m_pure = input2.data<impalib_type>();
            vector<vector<impalib_type>> eq_constraint_to_project_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

            for (int project_index=0; project_index<N_PROJECTS; project_index++){
                copy ( eq_constraint_to_project_m_pure + N_TEAMS*project_index, eq_constraint_to_project_m_pure + N_TEAMS*(project_index+1), eq_constraint_to_project_m[project_index].begin() );
            }
            
            vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

            projectIneqConstraint.project_inequality_constraint_update(eq_constraint_to_project_m, project_to_eq_constraint_m);

            fstream file_output("../ut_results/ut_projectIneqConstraint/ut_"+ ut_name+"/project_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&project_to_eq_constraint_m[i][j]), sizeof(project_to_eq_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << endl;}
    
    }

}