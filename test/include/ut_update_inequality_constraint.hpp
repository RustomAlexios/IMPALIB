    // Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_project_ineq_constraint(string&);

void ut_project_ineq_constraint(string& ut_name){

    const char *n_departments_bash=getenv("N_DEPARTMENTS");
    if(n_departments_bash == NULL)
    {cout << "n_departments_bash not available\n";}
    
    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/N_TEAMS_pure.npy");
    int* n_teams_pure = input1.data<int>();
    const int N_TEAMS = *n_teams_pure;
    
    const char *n_projects_bash=getenv("N_PROJECTS");
    if(n_projects_bash == NULL)
    {cout << "n_projects_bash not available\n";}

    const int N_DEPARTMENTS = atoi(n_departments_bash); 
    const int N_PROJECTS = atoi(n_projects_bash);

    InequalityConstraintKcMwm projectIneqConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS);

    if (ut_name == "ProjectInequalityConstraintUpdate"){

            cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/eq_constraint_to_project_m_pure.npy");
            impalib_type* eq_constraint_to_project_m_pure = input2.data<impalib_type>();
            vector<vector<impalib_type>> eq_constraint_to_project_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

            for (int project_index=0; project_index<N_PROJECTS; project_index++){
                copy ( eq_constraint_to_project_m_pure + N_TEAMS*project_index, eq_constraint_to_project_m_pure + N_TEAMS*(project_index+1), eq_constraint_to_project_m[project_index].begin() );
            }
            
            vector<vector<impalib_type>> project_to_eq_constraint_m(N_PROJECTS, vector<impalib_type>(N_TEAMS, zero_value));

            projectIneqConstraint.project_inequality_constraint_update(eq_constraint_to_project_m, project_to_eq_constraint_m);

            fstream file_output("../ut_results/project_to_eq_constraint_m_wrapper", ios::out | ios::binary | ios:: trunc);
            if (file_output.is_open()) {
                for (int i=0; i<N_PROJECTS; i++){
                for (int j=0; j<N_TEAMS; j++){
                    file_output.write((char*)(&project_to_eq_constraint_m[i][j]), sizeof(project_to_eq_constraint_m[i][j]));}}
                    file_output.close();}
            else {cout << "Error! File cannot be opened!" << "\n";}
    
    }

}


void ut_ineq_const_mobarp_update(string&);

void ut_ineq_const_mobarp_update(string& ut_name){

// InequalityConstraintMOBARP(const int NUM_FIXED_TX, const int NUM_MOBILE_TX, const int NUM_BANDS, 
//                                                                 const int NUM_TIME_STEPS, const int NUM_RX_LOCS, const int NUM_MOBILE_TX_LOCS, const impalib_type ALPHA, const bool FILTERING_FLAG

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

    const int NUM_FIXED_TX = atoi(n_fixed_tx_bash);  
    const int NUM_MOBILE_TX = atoi(n_mobile_tx_bash); 
    const int NUM_BANDS = atoi(n_bands_bash); 
    const int NUM_TIME_STEPS = atoi(n_time_steps_bash); 
    const int NUM_RX_LOCS = atoi(n_rx_locs_bash); 
    const int NUM_MOBILE_TX_LOCS = atoi(n_mobile_tx_locs_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/connectivity_mobile_tx.npy");
    int* connectivity_mobile_tx_pure = input1.data<int>();
    vector<vector<vector<vector<int>>>> ConnectivityMobileTx;

    for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
        ConnectivityMobileTx.push_back(vector<vector<vector<int>>>(NUM_BANDS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_RX_LOCS, 0))));
        for (int j=0; j<NUM_BANDS; j++){
            for (int k=0; k< NUM_TIME_STEPS; k++){
                copy(connectivity_mobile_tx_pure + NUM_RX_LOCS*k + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*n, connectivity_mobile_tx_pure + NUM_RX_LOCS*(k+1) + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*n, ConnectivityMobileTx[n][j][k].begin());
        }
        }
    }

    cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/connectivity_fixed_tx.npy");
    int* connectivity_fixed_tx_pure = input2.data<int>();
    vector<vector<vector<vector<int>>>> ConnectivityFixedTx;

    for (int i=0; i< NUM_FIXED_TX; i++){
        ConnectivityFixedTx.push_back(vector<vector<vector<int>>>(NUM_BANDS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_RX_LOCS, 0))));
        for (int j=0; j<NUM_BANDS; j++){
            for (int k=0; k< NUM_TIME_STEPS; k++){
                copy(connectivity_fixed_tx_pure + NUM_RX_LOCS*k + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*i, connectivity_fixed_tx_pure + NUM_RX_LOCS*(k+1) + NUM_RX_LOCS*NUM_TIME_STEPS * j + NUM_RX_LOCS*NUM_BANDS*NUM_TIME_STEPS*i, ConnectivityFixedTx[i][j][k].begin());
        }
        }
    }

    cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/fixed_capacity_constraints.npy");
    int* fixed_capacity_constraints_pure = input3.data<int>();
    vector<int> FixedCapacityConstraints;

    copy(fixed_capacity_constraints_pure, fixed_capacity_constraints_pure + NUM_FIXED_TX, back_inserter(FixedCapacityConstraints));

    cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/mobile_capacity_constraints.npy");
    int* mobile_capacity_constraints_pure = input4.data<int>();
    vector<int> MobileCapacityConstraints;

    copy(mobile_capacity_constraints_pure, mobile_capacity_constraints_pure + NUM_MOBILE_TX, back_inserter(MobileCapacityConstraints));

    InequalityConstraintMOBARP modelIneqConstraint(NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, 0.0, false);

    if (ut_name == "IneqCapacConstUpdate"){

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/fixed_x_eq_const_to_fixed_capac_const_m.npy");
        impalib_type* fixed_x_eq_const_to_fixed_capac_const_m_pure = input5.data<impalib_type>();
        vector<vector<vector<impalib_type>>> fixed_x_eq_const_to_fixed_capac_const_m;

        for (int i=0; i< NUM_FIXED_TX; i++){
            fixed_x_eq_const_to_fixed_capac_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(fixed_x_eq_const_to_fixed_capac_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, fixed_x_eq_const_to_fixed_capac_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, fixed_x_eq_const_to_fixed_capac_const_m[i][j].begin());
            }
        }

        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/mobile_x_eq_const_to_mobile_capac_const_m.npy");
        impalib_type* mobile_x_eq_const_to_mobile_capac_const_m_pure = input6.data<impalib_type>();
        vector<vector<vector<impalib_type>>> mobile_x_eq_const_to_mobile_capac_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            mobile_x_eq_const_to_mobile_capac_const_m.push_back(vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
            for (int j=0; j< NUM_BANDS; j++){
                copy(mobile_x_eq_const_to_mobile_capac_const_m_pure + NUM_TIME_STEPS * j + NUM_BANDS*NUM_TIME_STEPS*i, mobile_x_eq_const_to_mobile_capac_const_m_pure + NUM_TIME_STEPS * (j + 1) + NUM_BANDS*NUM_TIME_STEPS*i, mobile_x_eq_const_to_mobile_capac_const_m[i][j].begin());
            }
        }

        vector<vector<vector<impalib_type>>> FixedCapacConst2FixedXEqConstDummyM(NUM_FIXED_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));
        vector<vector<vector<impalib_type>>> MobileCapacConst2MobileXEqConstDummyM(NUM_MOBILE_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0)));

        modelIneqConstraint.ineq_capac_const_update(fixed_x_eq_const_to_fixed_capac_const_m, mobile_x_eq_const_to_mobile_capac_const_m, FixedCapacityConstraints, 
                                            MobileCapacityConstraints, FixedCapacConst2FixedXEqConstDummyM, MobileCapacConst2MobileXEqConstDummyM);

        fstream file_output_1("../ut_results/FixedCapacConst2FixedXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<NUM_FIXED_TX; i++){
            for (int j=0; j<NUM_BANDS; j++){
                for (int k=0; k<NUM_TIME_STEPS; k++){
                file_output_1.write((char*)(&FixedCapacConst2FixedXEqConstDummyM[i][j][k]), sizeof(FixedCapacConst2FixedXEqConstDummyM[i][j][k]));}}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/MobileCapacConst2MobileXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<NUM_MOBILE_TX; i++){
            for (int j=0; j<NUM_BANDS; j++){
                for (int k=0; k<NUM_TIME_STEPS; k++){
                file_output_2.write((char*)(&MobileCapacConst2MobileXEqConstDummyM[i][j][k]), sizeof(MobileCapacConst2MobileXEqConstDummyM[i][j][k]));}}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

   else if (ut_name == "MobileLocEqConst2REqConstUpdate"){

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/r_eq_const_to_mobile_loc_eq_const_m.npy");
        impalib_type* r_eq_const_to_mobile_loc_eq_const_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> r_eq_const_to_mobile_loc_eq_const_m;

        for (int i=0; i< NUM_MOBILE_TX; i++){
            r_eq_const_to_mobile_loc_eq_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(r_eq_const_to_mobile_loc_eq_const_m_pure + NUM_MOBILE_TX_LOCS * i, r_eq_const_to_mobile_loc_eq_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), r_eq_const_to_mobile_loc_eq_const_m[i].begin());
        }

        vector<vector<impalib_type>> MobileLocEqConst2REqConstDummyM(NUM_MOBILE_TX, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        modelIneqConstraint.mobile_loc_eq_const_to_r_eq_const_update(r_eq_const_to_mobile_loc_eq_const_m, MobileLocEqConst2REqConstDummyM);

        fstream file_output_1("../ut_results/MobileLocEqConst2REqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<MobileLocEqConst2REqConstDummyM.size(); i++){
            for (int j=0; j<MobileLocEqConst2REqConstDummyM[0].size(); j++){
                file_output_1.write((char*)(&MobileLocEqConst2REqConstDummyM[i][j]), sizeof(MobileLocEqConst2REqConstDummyM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

    else if (ut_name == "SetCoverIneqConstUpdate"){

        cnpy::NpyArray input5 = cnpy::npy_load("../ut_inputs/z_eq_const_to_set_cover_ineq_const_m.npy");
        impalib_type* z_eq_const_to_set_cover_ineq_const_m_pure = input5.data<impalib_type>();
        vector<vector<impalib_type>> z_eq_const_to_set_cover_ineq_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS*NUM_MOBILE_TX_LOCS; i++){
            z_eq_const_to_set_cover_ineq_const_m.push_back(vector<impalib_type>(NUM_RX_LOCS, 0));
            copy(z_eq_const_to_set_cover_ineq_const_m_pure + NUM_RX_LOCS * i, z_eq_const_to_set_cover_ineq_const_m_pure + NUM_RX_LOCS * (i + 1), z_eq_const_to_set_cover_ineq_const_m[i].begin());
        }

        cnpy::NpyArray input6 = cnpy::npy_load("../ut_inputs/fixed_x_eq_const_to_set_cover_const_m.npy");
        impalib_type* fixed_x_eq_const_to_set_cover_const_m_pure = input6.data<impalib_type>();
        vector<vector<impalib_type>> fixed_x_eq_const_to_set_cover_const_m;

        for (int i=0; i< NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            fixed_x_eq_const_to_set_cover_const_m.push_back(vector<impalib_type>(NUM_RX_LOCS, 0));
            copy(fixed_x_eq_const_to_set_cover_const_m_pure + NUM_RX_LOCS * i, fixed_x_eq_const_to_set_cover_const_m_pure + NUM_RX_LOCS * (i + 1), fixed_x_eq_const_to_set_cover_const_m[i].begin());
        }

        vector<vector<vector<vector<int>>>> TempConxMobTxRx(NUM_RX_LOCS, vector<vector<vector<int>>>(NUM_MOBILE_TX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))));

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

        for (int k=0; k< NUM_TIME_STEPS; k++){
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    for (int j_i=0; j_i <NUM_BANDS*NUM_MOBILE_TX; j_i++){
                        TempConxMobTxRx[l][n][k][j_i] = conx_mob_tx_rx[k][l][n][j_i];
                    }
                }
            }
        }

        cnpy::NpyArray input8 = cnpy::npy_load("../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy");
        int* conx_fixed_tx_per_num_rx_locs_pure = input8.data<int>();
        vector<vector<int>> conx_fixed_tx_per_num_rx_locs;

        for (int i=0; i< NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            conx_fixed_tx_per_num_rx_locs.push_back(vector<int>(NUM_RX_LOCS, 0));
            copy(conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * i, conx_fixed_tx_per_num_rx_locs_pure + NUM_RX_LOCS * (i + 1), conx_fixed_tx_per_num_rx_locs[i].begin());
        }

        vector<vector<int>> TempReshapedConnectivityFixedTx(NUM_FIXED_TX*NUM_BANDS, vector<int>(NUM_TIME_STEPS*NUM_RX_LOCS, 0));

        vector<int> flatData;
        for (const auto& row : conx_fixed_tx_per_num_rx_locs) {
            for (int val : row) {
                flatData.push_back(val);
            }
        }

        int index = 0;
        for (int i = 0; i < NUM_FIXED_TX*NUM_BANDS; ++i) {
            for (int j = 0; j < NUM_TIME_STEPS*NUM_RX_LOCS; ++j) {
                TempReshapedConnectivityFixedTx[i][j] = flatData[index++];
            }
        }

        vector<vector<int>> TempReshapedConxMobTxRx(NUM_MOBILE_TX_LOCS*NUM_BANDS*NUM_MOBILE_TX, vector<int>(NUM_RX_LOCS*NUM_TIME_STEPS, 0));

        vector<vector<vector<vector<int>>>> swaped_temp_conx_mob_tx_rx(NUM_BANDS*NUM_MOBILE_TX, vector<vector<vector<int>>>(NUM_MOBILE_TX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_RX_LOCS, 0))));
        
        for (int k=0; k< NUM_TIME_STEPS; k++){
            for (int l=0; l<NUM_RX_LOCS; l++){
                for (int n=0; n< NUM_MOBILE_TX_LOCS; n++){
                    for (int j_i=0; j_i <NUM_BANDS*NUM_MOBILE_TX; j_i++){
                        swaped_temp_conx_mob_tx_rx[j_i][n][k][l] = TempConxMobTxRx[l][n][k][j_i];
                    }
                }
            }
        }

        for (int j_i = 0; j_i < NUM_BANDS*NUM_MOBILE_TX; ++j_i) {
            for (int n = 0; n < NUM_MOBILE_TX_LOCS; ++n) {
                int row_index = j_i * NUM_MOBILE_TX_LOCS + n;
                int col_index = 0;
                for (int k = 0; k < NUM_TIME_STEPS; ++k) {
                    for (int l = 0; l < NUM_RX_LOCS; ++l) {
                        TempReshapedConxMobTxRx[row_index][col_index] = swaped_temp_conx_mob_tx_rx[j_i][n][k][l];
                        ++col_index;
                    }
                }
            }
        }


        vector<vector<vector<impalib_type>>> SetCoverIneqConst2FixedXEqConstDummyM(NUM_TIME_STEPS, vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0)));
        vector<vector<vector<vector<impalib_type>>>> SetCoverIneqConst2ZEqConstDummyM(NUM_TIME_STEPS, vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))));

        modelIneqConstraint.set_cover_ineq_const_update(z_eq_const_to_set_cover_ineq_const_m, fixed_x_eq_const_to_set_cover_const_m, 
                                            TempConxMobTxRx, TempReshapedConnectivityFixedTx, TempReshapedConxMobTxRx, SetCoverIneqConst2FixedXEqConstDummyM, SetCoverIneqConst2ZEqConstDummyM);

        fstream file_output_1("../ut_results/SetCoverIneqConst2ZEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<SetCoverIneqConst2ZEqConstDummyM.size(); i++){
            for (int j=0; j<SetCoverIneqConst2ZEqConstDummyM[0].size(); j++){
                for (int k=0; k<SetCoverIneqConst2ZEqConstDummyM[0][0].size(); k++){
                    for (int l=0; l<SetCoverIneqConst2ZEqConstDummyM[0][0][0].size(); l++){
                        file_output_1.write((char*)(&SetCoverIneqConst2ZEqConstDummyM[i][j][k][l]), sizeof(SetCoverIneqConst2ZEqConstDummyM[i][j][k][l]));}}}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}


        fstream file_output_2("../ut_results/SetCoverIneqConst2FixedXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<SetCoverIneqConst2FixedXEqConstDummyM.size(); i++){
            for (int j=0; j<SetCoverIneqConst2FixedXEqConstDummyM[0].size(); j++){
                for (int k=0; k<SetCoverIneqConst2FixedXEqConstDummyM[0][0].size(); k++){
                        file_output_2.write((char*)(&SetCoverIneqConst2FixedXEqConstDummyM[i][j][k]), sizeof(SetCoverIneqConst2FixedXEqConstDummyM[i][j][k]));}}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

}
