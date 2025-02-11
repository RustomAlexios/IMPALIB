// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

void ut_auxiliary_const_update(string&);

void ut_auxiliary_const_update(string& ut_name){
    
    // AuxiliaryConstraint(const int NUM_MOBILE_TX, const int NUM_BANDS, const int NUM_TIME_STEPS, const int NUM_MOBILE_TX_LOCS)

    const char *n_mobile_tx_bash=getenv("NUM_MOBILE_TX");
    if(n_mobile_tx_bash == NULL)
    {cout << "n_mobile_tx_bash not available\n";}

    const char *n_bands_bash=getenv("NUM_BANDS");
    if(n_bands_bash == NULL)
    {cout << "n_bands_bash not available\n";}

    const char *n_time_steps_bash=getenv("NUM_TIME_STEPS");
    if(n_time_steps_bash == NULL)
    {cout << "n_time_steps_bash not available\n";}

    const char *n_mobile_tx_locs_bash=getenv("NUM_MOBILE_TX_LOCS");
    if(n_mobile_tx_locs_bash == NULL)
    {cout << "n_mobile_tx_locs_bash not available\n";}


    const int NUM_MOBILE_TX = atoi(n_mobile_tx_bash); 
    const int NUM_BANDS = atoi(n_bands_bash); 
    const int NUM_TIME_STEPS = atoi(n_time_steps_bash); 
    const int NUM_MOBILE_TX_LOCS = atoi(n_mobile_tx_locs_bash);

    cnpy::NpyArray input1 = cnpy::npy_load("../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy");
    int* conx_mob_tx_per_num_mob_tx_locs_pure = input1.data<int>();
    vector<vector<int>> ConxMobTxPerNumMobTxLocs;

    for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
        ConxMobTxPerNumMobTxLocs.push_back(vector<int>(NUM_MOBILE_TX_LOCS, 0));
        copy(conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * i, conx_mob_tx_per_num_mob_tx_locs_pure + NUM_MOBILE_TX_LOCS * (i + 1), ConxMobTxPerNumMobTxLocs[i].begin());
    }


    AuxiliaryConstraint modelAuxiliaryConstraint(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS);

    // if (ut_name == "AuxiliaryConst2ZEqConstUpdate"){

    //     cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy");
    //     impalib_type* mobile_x_eq_const_to_auxiliary_const_m_pure = input2.data<impalib_type>();
    //     vector<vector<impalib_type>> mobile_x_eq_const_to_auxiliary_const_m;

    //     for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
    //         mobile_x_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
    //         copy(mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), mobile_x_eq_const_to_auxiliary_const_m[i].begin());
    //     }

    //     cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/r_eq_const_to_auxiliary_const_m.npy");
    //     impalib_type* r_eq_const_to_auxiliary_const_m_pure = input3.data<impalib_type>();
    //     vector<vector<impalib_type>> r_eq_const_to_auxiliary_const_m;

    //     for (int i=0; i< NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS; i++){
    //         r_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_BANDS*NUM_TIME_STEPS, 0));
    //         copy(r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * i, r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * (i + 1), r_eq_const_to_auxiliary_const_m[i].begin());
    //     }

    //     vector<vector<impalib_type>> AuxiliaryConst2ZEqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

    //     modelAuxiliaryConstraint.auxiliary_const_to_z_eq_const_update(mobile_x_eq_const_to_auxiliary_const_m, r_eq_const_to_auxiliary_const_m, AuxiliaryConst2ZEqConstM, 
    //                                         ConxMobTxPerNumMobTxLocs);

    //     fstream file_output_1("../ut_results/AuxiliaryConst2ZEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    //     if (file_output_1.is_open()) {
    //         for (int i=0; i<AuxiliaryConst2ZEqConstM.size(); i++){
    //         for (int j=0; j<AuxiliaryConst2ZEqConstM[0].size(); j++){
    //             file_output_1.write((char*)(&AuxiliaryConst2ZEqConstM[i][j]), sizeof(AuxiliaryConst2ZEqConstM[i][j]));}}
    //             file_output_1.close();
    //             }
    //     else {cout << "Error! File cannot be opened!" << "\n";}

    // }

    // else if (ut_name == "AuxiliaryConst2REqConstUpdate"){
        
    //     cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/z_eq_const_to_auxiliary_const_m.npy");
    //     impalib_type* z_eq_const_to_auxiliary_const_m_pure = input2.data<impalib_type>();
    //     vector<vector<impalib_type>> z_eq_const_to_auxiliary_const_m;

    //     for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
    //         z_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
    //         copy(z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), z_eq_const_to_auxiliary_const_m[i].begin());
    //     }

    //     cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy");
    //     impalib_type* mobile_x_eq_const_to_auxiliary_const_m_pure = input3.data<impalib_type>();
    //     vector<vector<impalib_type>> mobile_x_eq_const_to_auxiliary_const_m;

    //     for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
    //         mobile_x_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
    //         copy(mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), mobile_x_eq_const_to_auxiliary_const_m[i].begin());
    //     }

    //     vector<vector<impalib_type>> AuxiliaryConst2REqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

    //     modelAuxiliaryConstraint.auxiliary_const_to_r_eq_const_update(z_eq_const_to_auxiliary_const_m, mobile_x_eq_const_to_auxiliary_const_m, ConxMobTxPerNumMobTxLocs, AuxiliaryConst2REqConstM);

    //     fstream file_output_1("../ut_results/AuxiliaryConst2REqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    //     if (file_output_1.is_open()) {
    //         for (int i=0; i<AuxiliaryConst2REqConstM.size(); i++){
    //         for (int j=0; j<AuxiliaryConst2REqConstM[0].size(); j++){
    //             file_output_1.write((char*)(&AuxiliaryConst2REqConstM[i][j]), sizeof(AuxiliaryConst2REqConstM[i][j]));}}
    //             file_output_1.close();
    //             }
    //     else {cout << "Error! File cannot be opened!" << "\n";}

    // }

    if (ut_name == "AuxiliaryConst2ZAndREqConstUpdate"){

        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy");
        impalib_type* mobile_x_eq_const_to_auxiliary_const_m_pure = input2.data<impalib_type>();
        vector<vector<impalib_type>> mobile_x_eq_const_to_auxiliary_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            mobile_x_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, mobile_x_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), mobile_x_eq_const_to_auxiliary_const_m[i].begin());
        }

        cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/r_eq_const_to_auxiliary_const_m.npy");
        impalib_type* r_eq_const_to_auxiliary_const_m_pure = input3.data<impalib_type>();
        vector<vector<impalib_type>> r_eq_const_to_auxiliary_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS; i++){
            r_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_BANDS*NUM_TIME_STEPS, 0));
            copy(r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * i, r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * (i + 1), r_eq_const_to_auxiliary_const_m[i].begin());
        }

        vector<vector<impalib_type>> AuxiliaryConst2ZEqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        // modelAuxiliaryConstraint.auxiliary_const_to_z_eq_const_update(mobile_x_eq_const_to_auxiliary_const_m, r_eq_const_to_auxiliary_const_m, AuxiliaryConst2ZEqConstM, 
        //                                     ConxMobTxPerNumMobTxLocs);

        cnpy::NpyArray input4 = cnpy::npy_load("../ut_inputs/z_eq_const_to_auxiliary_const_m.npy");
        impalib_type* z_eq_const_to_auxiliary_const_m_pure = input4.data<impalib_type>();
        vector<vector<impalib_type>> z_eq_const_to_auxiliary_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            z_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), z_eq_const_to_auxiliary_const_m[i].begin());
        }

        vector<vector<impalib_type>> AuxiliaryConst2REqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        // modelAuxiliaryConstraint.auxiliary_const_to_r_eq_const_update(z_eq_const_to_auxiliary_const_m, mobile_x_eq_const_to_auxiliary_const_m, ConxMobTxPerNumMobTxLocs, AuxiliaryConst2REqConstM);

        modelAuxiliaryConstraint.auxiliary_const_to_z_and_r_eq_const_update(z_eq_const_to_auxiliary_const_m, mobile_x_eq_const_to_auxiliary_const_m,
                                ConxMobTxPerNumMobTxLocs, r_eq_const_to_auxiliary_const_m, AuxiliaryConst2ZEqConstM, AuxiliaryConst2REqConstM);

        fstream file_output_1("../ut_results/AuxiliaryConst2ZEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<AuxiliaryConst2ZEqConstM.size(); i++){
            for (int j=0; j<AuxiliaryConst2ZEqConstM[0].size(); j++){
                file_output_1.write((char*)(&AuxiliaryConst2ZEqConstM[i][j]), sizeof(AuxiliaryConst2ZEqConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

        fstream file_output_2("../ut_results/AuxiliaryConst2REqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_2.is_open()) {
            for (int i=0; i<AuxiliaryConst2REqConstM.size(); i++){
            for (int j=0; j<AuxiliaryConst2REqConstM[0].size(); j++){
                file_output_2.write((char*)(&AuxiliaryConst2REqConstM[i][j]), sizeof(AuxiliaryConst2REqConstM[i][j]));}}
                file_output_2.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }
    
    else if (ut_name == "AuxiliaryConst2MobileXEqConstUpdate"){
        
        cnpy::NpyArray input2 = cnpy::npy_load("../ut_inputs/z_eq_const_to_auxiliary_const_m.npy");
        impalib_type* z_eq_const_to_auxiliary_const_m_pure = input2.data<impalib_type>();
        vector<vector<impalib_type>> z_eq_const_to_auxiliary_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS; i++){
            z_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));
            copy(z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * i, z_eq_const_to_auxiliary_const_m_pure + NUM_MOBILE_TX_LOCS * (i + 1), z_eq_const_to_auxiliary_const_m[i].begin());
        }

        cnpy::NpyArray input3 = cnpy::npy_load("../ut_inputs/r_eq_const_to_auxiliary_const_m.npy");
        impalib_type* r_eq_const_to_auxiliary_const_m_pure = input3.data<impalib_type>();
        vector<vector<impalib_type>> r_eq_const_to_auxiliary_const_m;

        for (int i=0; i< NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS; i++){
            r_eq_const_to_auxiliary_const_m.push_back(vector<impalib_type>(NUM_BANDS*NUM_TIME_STEPS, 0));
            copy(r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * i, r_eq_const_to_auxiliary_const_m_pure + NUM_BANDS*NUM_TIME_STEPS * (i + 1), r_eq_const_to_auxiliary_const_m[i].begin());
        }

        vector<vector<impalib_type>> AuxiliaryConst2MobileXEqConstM(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0));

        modelAuxiliaryConstraint.auxiliary_const_to_mobile_x_eq_const_update(z_eq_const_to_auxiliary_const_m, r_eq_const_to_auxiliary_const_m, ConxMobTxPerNumMobTxLocs, AuxiliaryConst2MobileXEqConstM);

        fstream file_output_1("../ut_results/AuxiliaryConst2MobileXEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
        if (file_output_1.is_open()) {
            for (int i=0; i<AuxiliaryConst2MobileXEqConstM.size(); i++){
            for (int j=0; j<AuxiliaryConst2MobileXEqConstM[0].size(); j++){
                file_output_1.write((char*)(&AuxiliaryConst2MobileXEqConstM[i][j]), sizeof(AuxiliaryConst2MobileXEqConstM[i][j]));}}
                file_output_1.close();
                }
        else {cout << "Error! File cannot be opened!" << "\n";}

    }

}
