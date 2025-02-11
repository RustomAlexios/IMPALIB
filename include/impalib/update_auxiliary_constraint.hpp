// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for AUXILIARY constraint
 */
class AuxiliaryConstraint {
    int numMobileTx_;
    int numBands_;
    int numTimeSteps_;
    int numMobileTxLocs_;

   public:
    // void auxiliary_const_to_z_eq_const_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>>&, vector<vector<impalib_type>>&, vector<vector<int>> &) const;

    vector<vector<impalib_type>> reshape_function(const vector<vector<impalib_type>> &) const;
    
    // void auxiliary_const_to_r_eq_const_update(vector<vector<impalib_type>>&, vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<impalib_type>> &) const;

    void auxiliary_const_to_mobile_x_eq_const_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<impalib_type>> & ) const;

    void auxiliary_const_to_z_and_r_eq_const_update(vector<vector<impalib_type>>&, vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<impalib_type>>&,
                                vector<vector<impalib_type>>&, vector<vector<impalib_type>> &) const;

    AuxiliaryConstraint(int NUM_MOBILE_TX, int NUM_BANDS, int NUM_TIME_STEPS, int NUM_MOBILE_TX_LOCS); ///< constructor
};


inline AuxiliaryConstraint::AuxiliaryConstraint(const int NUM_MOBILE_TX, const int NUM_BANDS, const int NUM_TIME_STEPS, const int NUM_MOBILE_TX_LOCS)
    : numMobileTx_(NUM_MOBILE_TX),
      numBands_(NUM_BANDS),
      numTimeSteps_(NUM_TIME_STEPS),
      numMobileTxLocs_(NUM_MOBILE_TX_LOCS){};

// inline void AuxiliaryConstraint::auxiliary_const_to_z_eq_const_update(vector<vector<impalib_type>> & rMobileXEqConst2AuxiliaryConstM, vector<vector<impalib_type>>& rREqConst2AuxiliaryConstM, 
//                                         vector<vector<impalib_type>>& rAuxiliaryConst2ZEqConstM, vector<vector<int>> & rConxMobTxPerNumMobTxLocs) const {
        
//         auto reshaped_r_eq_const_to_auxiliary_const = reshape_function(rREqConst2AuxiliaryConstM);

//         // fstream file_output("./ut_results/reshapedREqConst2AuxiliaryConstM_wrapper", ios::out | ios::binary | ios:: trunc);
//         // if (file_output.is_open()) {
//         //     for (int i=0; i<reshaped_r_eq_const_to_auxiliary_const.size(); i++){
//         //     for (int j=0; j<reshaped_r_eq_const_to_auxiliary_const[0].size(); j++){
//         //         file_output.write((char*)(&reshaped_r_eq_const_to_auxiliary_const[i][j]), sizeof(reshaped_r_eq_const_to_auxiliary_const[i][j]));}}
//         //         file_output.close();
//         //         }
//         // else {cout << "Error! File cannot be opened!" << "\n";}

//         for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++){
//             for (int n=0; n<numMobileTxLocs_; n++){
//                 if (rConxMobTxPerNumMobTxLocs[i][n] == 1){
//                     rAuxiliaryConst2ZEqConstM[i][n] = max(rMobileXEqConst2AuxiliaryConstM[i][n] + reshaped_r_eq_const_to_auxiliary_const[i][n], max(rMobileXEqConst2AuxiliaryConstM[i][n], reshaped_r_eq_const_to_auxiliary_const[i][n]));
//                 }
//             }
//         }

//     // cout << "Matrix rAuxiliaryConst2ZEqConstM_:" << "\n";
//     // for (const auto& row : rAuxiliaryConst2ZEqConstM_) {
//     //     for (const auto& val : row) {
//     //         cout << setw(10) << val;
//     //     }
//     //     cout << "\n";
//     // }

//     // exit(0);

//     // fstream file_output_1("./ut_results/rAuxiliaryConst2ZEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
//     // if (file_output_1.is_open()) {
//     //     for (int i=0; i<rAuxiliaryConst2ZEqConstM.size(); i++){
//     //     for (int j=0; j<rAuxiliaryConst2ZEqConstM[0].size(); j++){
//     //         file_output_1.write((char*)(&rAuxiliaryConst2ZEqConstM[i][j]), sizeof(rAuxiliaryConst2ZEqConstM[i][j]));}}
//     //         file_output_1.close();
//     //         }
//     // else {cout << "Error! File cannot be opened!" << "\n";}

// }

// inline void AuxiliaryConstraint::auxiliary_const_to_r_eq_const_update(vector<vector<impalib_type>>& rZEqConst2AuxiliaryConstM, vector<vector<impalib_type>> &rMobileXEqConst2AuxiliaryConstM, 
//                                                         vector<vector<int>> & rConxMobTxPerNumMobTxLocs, vector<vector<impalib_type>> &rAuxiliaryConst2REqConstM) const {

//         for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++){
//             for (int n=0; n<numMobileTxLocs_; n++){
//                 if (rConxMobTxPerNumMobTxLocs[i][n] == 1){
//                     rAuxiliaryConst2REqConstM[i][n] = min(max(0.0, -rMobileXEqConst2AuxiliaryConstM[i][n]), max(rZEqConst2AuxiliaryConstM[i][n], rZEqConst2AuxiliaryConstM[i][n] + rMobileXEqConst2AuxiliaryConstM[i][n]));
//                 }
//             }
//         }

//     // cout << "Matrix rAuxiliaryConst2REqConstM_:" << "\n";
//     // for (const auto& row : rAuxiliaryConst2REqConstM_) {
//     //     for (const auto& val : row) {
//     //         cout << setw(10) << val;
//     //     }
//     //     cout << "\n";
//     // }

//     // exit(0);

//     // fstream file_output_1("./ut_results/rAuxiliaryConst2REqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
//     // if (file_output_1.is_open()) {
//     //     for (int i=0; i<rAuxiliaryConst2REqConstM.size(); i++){
//     //     for (int j=0; j<rAuxiliaryConst2REqConstM[0].size(); j++){
//     //         file_output_1.write((char*)(&rAuxiliaryConst2REqConstM[i][j]), sizeof(rAuxiliaryConst2REqConstM[i][j]));}}
//     //         file_output_1.close();
//     //         }
//     // else {cout << "Error! File cannot be opened!" << "\n";}

// }


inline void AuxiliaryConstraint::auxiliary_const_to_z_and_r_eq_const_update(vector<vector<impalib_type>>& rZEqConst2AuxiliaryConstM, vector<vector<impalib_type>> & rMobileXEqConst2AuxiliaryConstM,
                                vector<vector<int>> &  rConxMobTxPerNumMobTxLocs, vector<vector<impalib_type>>& rREqConst2AuxiliaryConstM,
                                vector<vector<impalib_type>>& rAuxiliaryConst2ZEqConstM, vector<vector<impalib_type>> & rAuxiliaryConst2REqConstM)const{

        auto reshaped_r_eq_const_to_auxiliary_const = reshape_function(rREqConst2AuxiliaryConstM);

        for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++){
            for (int n=0; n<numMobileTxLocs_; n++){
                if (rConxMobTxPerNumMobTxLocs[i][n] == 1){
                    rAuxiliaryConst2ZEqConstM[i][n] = max(rMobileXEqConst2AuxiliaryConstM[i][n] + reshaped_r_eq_const_to_auxiliary_const[i][n], max(rMobileXEqConst2AuxiliaryConstM[i][n], reshaped_r_eq_const_to_auxiliary_const[i][n]));
                    rAuxiliaryConst2REqConstM[i][n] = min(max(0.0, -rMobileXEqConst2AuxiliaryConstM[i][n]), max(rZEqConst2AuxiliaryConstM[i][n], rZEqConst2AuxiliaryConstM[i][n] + rMobileXEqConst2AuxiliaryConstM[i][n]));
                }
            }
        }
}


inline vector<vector<impalib_type>> AuxiliaryConstraint::reshape_function(const vector<vector<impalib_type>> & rREqConst2AuxiliaryConstM) const
{
    //rREqConst2AuxiliaryConstM has shape (numMobileTx_*numMobileTxLocs_, numBands_*numTimeSteps_)

    // auto start_1 = chrono::high_resolution_clock::now();
    // vector<vector<vector<impalib_type>>> reshaped_1(numMobileTx_, vector<vector<impalib_type>>(numMobileTxLocs_, vector<impalib_type>(numBands_*numTimeSteps_, 0)));

    // vector<vector<vector<impalib_type>>> reshaped_2(numMobileTx_, vector<vector<impalib_type>>(numBands_*numTimeSteps_, vector<impalib_type>(numMobileTxLocs_, 0)));

    // vector<vector<impalib_type>> reshaped_3(numBands_*numTimeSteps_*numMobileTx_, vector<impalib_type>(numMobileTxLocs_, 0));

    // for (size_t i = 0; i < numMobileTx_; i++) {
    //     for (size_t n = 0; n < numMobileTxLocs_; n++) {
    //         for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
    //             reshaped_1[i][n][j_k] = rREqConst2AuxiliaryConstM[i * numMobileTxLocs_ + n][j_k];
    //         }
    //     }
    // }

    // for (size_t i = 0; i < numMobileTx_; i++) {
    //     for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
    //         for (size_t n = 0; n < numMobileTxLocs_; n++) {
    //             reshaped_2[i][j_k][n] = reshaped_1[i][n][j_k];
    //         }
    //     }
    // }

    // for (size_t i = 0; i < numMobileTx_; i++) {
    //     for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
    //         for (size_t n = 0; n < numMobileTxLocs_; n++) {
    //             reshaped_3[i * (numBands_ * numTimeSteps_) + j_k][n] = reshaped_2[i][j_k][n];
    //         }
    //     }
    // }

    // auto end_1 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_1 = end_1 - start_1;
    // cout << "Execution time unoptimized: " << elapsed_1.count() << " seconds\n";

    // auto start_2 = chrono::high_resolution_clock::now();
    vector<vector<impalib_type>> reshaped_3(numBands_*numTimeSteps_*numMobileTx_, vector<impalib_type>(numMobileTxLocs_, 0));
    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
            for (size_t n = 0; n < numMobileTxLocs_; n++) {
                reshaped_3[i * (numBands_ * numTimeSteps_) + j_k][n] = rREqConst2AuxiliaryConstM[i * numMobileTxLocs_ + n][j_k];
            }
        }
    }

    // auto end_2 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_2 = end_2 - start_2;
    // cout << "Execution time optimized: " << elapsed_2.count() << " seconds\n";

    // for (int i = 0; i < numBands_*numTimeSteps_*numMobileTx_; ++i) {
    //     for (int n = 0; n < numMobileTxLocs_; ++n) {
    //         if (std::abs(reshaped_3_new[i][n] - reshaped_3[i][n]) > 1e-7) {
    //             cout<<"reshaped_3_new[i][n]: "<<reshaped_3_new[i][n]<<"\n";
    //             cout<<"reshaped_3[i][n]: "<<reshaped_3_new[i][n]<<"\n";
    //             exit(EXIT_FAILURE);
    //         }
    //     }
    // }


    // for (size_t i = 0; i < reshaped_3.size(); ++i) {
    //     for (size_t j = 0; j < reshaped_3[i].size(); ++j) {
    //         cout << reshaped_3[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    // exit(0);
    return reshaped_3;

}


inline void AuxiliaryConstraint::auxiliary_const_to_mobile_x_eq_const_update(vector<vector<impalib_type>> & rZEqConst2AuxiliaryConstM, vector<vector<impalib_type>> & rREqConst2AuxiliaryConstM, 
                                    vector<vector<int>> & rConxMobTxPerNumMobTxLocs, vector<vector<impalib_type>> & rAuxiliaryConst2MobileXEqConstM) const {

    auto reshaped_r_eq_const_to_auxiliary_const = reshape_function(rREqConst2AuxiliaryConstM);

    for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++){
        for (int n=0; n<numMobileTxLocs_; n++){
            if (rConxMobTxPerNumMobTxLocs[i][n] == 1){
                rAuxiliaryConst2MobileXEqConstM[i][n] = min(max(0.0, -reshaped_r_eq_const_to_auxiliary_const[i][n]), max(rZEqConst2AuxiliaryConstM[i][n], rZEqConst2AuxiliaryConstM[i][n] + reshaped_r_eq_const_to_auxiliary_const[i][n]));
            }
        }
    }

    // cout<<"-------------------------------------------------------------------------------------------------"<<"\n";
    // for (int i=0; i<numBands_*numMobileTx_*numTimeSteps_; i++){
    //     for (int j=0; j<numMobileTxLocs_; j++){
    //         cout<<setw(13)<<rAuxiliaryConst2MobileXEqConstM_[i][j];
    //     }
    //     cout<<"\n";
    // }
    // // exit(0);

    // fstream file_output_1("./ut_results/rAuxiliaryConst2MobileXEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<rAuxiliaryConst2MobileXEqConstM.size(); i++){
    //     for (int j=0; j<rAuxiliaryConst2MobileXEqConstM[0].size(); j++){
    //         file_output_1.write((char*)(&rAuxiliaryConst2MobileXEqConstM[i][j]), sizeof(rAuxiliaryConst2MobileXEqConstM[i][j]));}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

}