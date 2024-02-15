// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib_unit_tests.hpp"
#include "ut_utils.hpp"

void ut_forward_backward(string&);
void ut_extrinsic_output_department(string&);
void ut_team_to_knapsack_update(string&);
void ut_process_extrinsic_output_department(string&);

void ut_forward_backward(string& ut_name){

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}

    const bool FILT_FLAG(filt_flag_bash);

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
    const int* non_zero_weight_indices_sizes_pure = input2.data<int>();

    int max_size_non_zero_weight = *max_element(non_zero_weight_indices_sizes_pure , non_zero_weight_indices_sizes_pure + N_DEPARTMENTS);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    int* teams_weights_per_department_pure = input3.data<int>();
    vector<vector<int>> teams_weights_per_department(N_DEPARTMENTS, vector<int>(N_TEAMS,0));
    
    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( teams_weights_per_department_pure + N_TEAMS*department_index, teams_weights_per_department_pure+N_TEAMS*(department_index+1), teams_weights_per_department[department_index].begin() );
    }

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    int* non_zero_weight_indices_arr_pure = input4.data<int>();
    vector<vector<int>> non_zero_weight_indices;

        for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
            non_zero_weight_indices.push_back(vector<int>(non_zero_weight_indices_sizes_pure[department_index],0));
        
        copy ( non_zero_weight_indices_arr_pure + max_size_non_zero_weight*department_index, 
        non_zero_weight_indices_arr_pure+non_zero_weight_indices_sizes_pure[department_index]+
            max_size_non_zero_weight*department_index, non_zero_weight_indices[department_index].begin() );
    }

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* team_to_knapsack_m_pure = input5.data<impalib_type>();
    vector<vector<impalib_type>> team_to_knapsack_m(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( team_to_knapsack_m_pure + N_TEAMS*department_index, team_to_knapsack_m_pure+N_TEAMS*(department_index+1), team_to_knapsack_m[department_index].begin() );
    }
    
    Knapsack modelKnapsacks(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA);

    for (int department_index = 0; department_index<N_DEPARTMENTS; department_index++){
        int max_state_department = Nu[department_index];
        if (ut_name == "KnapsackForward"){
            vector<vector<impalib_type>> stage_forward_messages(N_TEAMS+1, vector<impalib_type>(max_state_department+1,zero_value));
            modelKnapsacks.forward(department_index, stage_forward_messages, max_state_department,
            non_zero_weight_indices, non_zero_weight_indices_sizes_pure,
            teams_weights_per_department, team_to_knapsack_m);
            fstream file_output("../ut_results/forward_wrapper"+ std::to_string(department_index), ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_TEAMS+1; i++){
                for (int j=0; j<max_state_department+1; j++){
                    file_output.write((char*)(&stage_forward_messages[i][j]), sizeof(stage_forward_messages[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
        
        }

        else if (ut_name == "KnapsackBackward"){
            vector<vector<impalib_type>> stage_backward_messages(N_TEAMS+1, vector<impalib_type>(max_state_department+1,zero_value));
            modelKnapsacks.backward(department_index, stage_backward_messages, max_state_department,
            non_zero_weight_indices, non_zero_weight_indices_sizes_pure,
            teams_weights_per_department, team_to_knapsack_m);
            
            fstream file_output("../ut_results/backward_wrapper"+ std::to_string(department_index), ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_TEAMS+1; i++){
                for (int j=0; j<max_state_department+1; j++){
                    file_output.write((char*)(&stage_backward_messages[i][j]), sizeof(stage_backward_messages[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
        }
    }
}

void ut_extrinsic_output_department(string& ut_name){

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const bool FILT_FLAG(filt_flag_bash);

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

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/teams_weights_per_department_pure.npy");
    int* teams_weights_per_department_pure = input2.data<int>();
    vector<vector<int>> teams_weights_per_department(N_DEPARTMENTS, vector<int>(N_TEAMS,0));
    
    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( teams_weights_per_department_pure + N_TEAMS*department_index, teams_weights_per_department_pure+N_TEAMS*(department_index+1), teams_weights_per_department[department_index].begin() );
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/team_to_knapsack_m_pure.npy");
    impalib_type* team_to_knapsack_m_pure = input3.data<impalib_type>();
    vector<vector<impalib_type>> team_to_knapsack_m(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( team_to_knapsack_m_pure + N_TEAMS*department_index, team_to_knapsack_m_pure+N_TEAMS*(department_index+1), team_to_knapsack_m[department_index].begin() );
    }
    
    Knapsack modelKnapsacks(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA);

    vector<vector<impalib_type>> extrinsicOutputDepartment(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int department_index = 0; department_index<N_DEPARTMENTS; department_index++){
        int max_state_department = Nu[department_index];
        vector<vector<impalib_type>> stage_forward_messages(N_TEAMS+1, vector<impalib_type>(max_state_department+1,zero_value));
        vector<vector<impalib_type>> stage_backward_messages(N_TEAMS+1, vector<impalib_type>(max_state_department+1,zero_value));
        cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/stage_forward_messages"+std::to_string(department_index)+".npy");
        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/stage_backward_messages"+std::to_string(department_index)+".npy");
        impalib_type* stage_forward_messages_pure = input1.data<impalib_type>(); 
        impalib_type* stage_backward_messages_pure = input2.data<impalib_type>();
        
        for (int i=0; i<N_TEAMS+1; i++){
            copy ( stage_forward_messages_pure + (max_state_department+1)*i, stage_forward_messages_pure + (max_state_department+1)*(i+1), stage_forward_messages[i].begin() );
            copy ( stage_backward_messages_pure + (max_state_department+1)*i, stage_backward_messages_pure + (max_state_department+1)*(i+1), stage_backward_messages[i].begin() );
        }

        modelKnapsacks.extrinsic_output_department_lhs(teams_weights_per_department,
                            stage_forward_messages, team_to_knapsack_m,
                            department_index, stage_backward_messages,
                            max_state_department, extrinsicOutputDepartment);
    }
        
        fstream file_output("../ut_results/extrinsic_output_department_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output.is_open()) {
            for (int i=0; i<N_DEPARTMENTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&extrinsicOutputDepartment[i][j]), sizeof(extrinsicOutputDepartment[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}

}


void ut_team_to_knapsack_update(string& ut_name){

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const bool FILT_FLAG(filt_flag_bash);

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
    const int* non_zero_weight_indices_sizes_pure = input2.data<int>();

    int max_size_non_zero_weight = *max_element(non_zero_weight_indices_sizes_pure , non_zero_weight_indices_sizes_pure + N_DEPARTMENTS);

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/non_zero_weight_indices_arr_pure.npy");
    int* non_zero_weight_indices_arr_pure = input3.data<int>();
    vector<vector<int>> non_zero_weight_indices;

        for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
            non_zero_weight_indices.push_back(vector<int>(non_zero_weight_indices_sizes_pure[department_index],0));
        
        copy ( non_zero_weight_indices_arr_pure + max_size_non_zero_weight*department_index, 
        non_zero_weight_indices_arr_pure+non_zero_weight_indices_sizes_pure[department_index]+
            max_size_non_zero_weight*department_index, non_zero_weight_indices[department_index].begin() );
        }
    
    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/reward_team_pure.npy");
    impalib_type* reward_team_pure = input4.data<impalib_type>();
    vector<impalib_type> reward_team(N_TEAMS,zero_value);
    copy ( reward_team_pure, reward_team_pure+N_TEAMS, reward_team.begin() );

    cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/extrinsic_output_department_pure.npy");
    impalib_type* extrinsic_output_department_pure = input5.data<impalib_type>();
    vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure+N_TEAMS*(department_index+1), extrinsic_output_department[department_index].begin() );
    }

    cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/oric_to_team_m_pure.npy");
    impalib_type* oric_to_team_m_pure = input6.data<impalib_type>();
    vector<impalib_type> oric_to_team_m(N_TEAMS,zero_value);
    copy ( oric_to_team_m_pure, oric_to_team_m_pure+N_TEAMS, oric_to_team_m.begin() );

    vector<vector<impalib_type>> team_to_knapsack_m(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    Knapsack modelKnapsacks(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA);
    modelKnapsacks.team_to_knapsack_update(non_zero_weight_indices,
                                team_to_knapsack_m, 
                                reward_team, extrinsic_output_department,
                                oric_to_team_m);
    
    fstream file_output("../ut_results/team_to_knapsack_m_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output.is_open()) {
            for (int i=0; i<N_DEPARTMENTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&team_to_knapsack_m[i][j]), sizeof(team_to_knapsack_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
}

void ut_process_extrinsic_output_department(string& ut_name){

    const char *filt_flag_bash=getenv("FILT_FLAG");
    if(filt_flag_bash == NULL)
    {cout << "filt_flag_bash not available\n";}
    
    const char *n_iter_bash=getenv("N_ITER");
    if(n_iter_bash == NULL)
    {cout << "n_iter_bash not available\n";} 

    const bool FILT_FLAG(filt_flag_bash);
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

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/extrinsic_output_department_dummy_pure.npy");
    impalib_type* extrinsic_output_department_pure = input2.data<impalib_type>();
    vector<vector<impalib_type>> extrinsic_output_department_dummy(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
        copy ( extrinsic_output_department_pure + N_TEAMS*department_index, extrinsic_output_department_pure+N_TEAMS*(department_index+1), extrinsic_output_department_dummy[department_index].begin() );
    }

    Knapsack modelKnapsacks(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA);
    
    vector<vector<impalib_type>> extrinsic_output_department(N_DEPARTMENTS, vector<impalib_type>(N_TEAMS,zero_value));

    for (int iter=0; iter<N_ITER; iter++){
        for (int department_index=0; department_index<N_DEPARTMENTS; department_index++){
            modelKnapsacks.process_extrinsic_output_department(department_index,
                                    iter, 
                                    extrinsic_output_department_dummy, extrinsic_output_department);
        }
    }

    fstream file_output("../ut_results/extrinsic_output_department_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output.is_open()) {
            for (int i=0; i<N_DEPARTMENTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&extrinsic_output_department[i][j]), sizeof(extrinsic_output_department[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}

}