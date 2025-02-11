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

    vector<vector<impalib_type>> reshape_function(const vector<vector<impalib_type>> &) const;

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

    vector<vector<impalib_type>> reshaped_3(numBands_*numTimeSteps_*numMobileTx_, vector<impalib_type>(numMobileTxLocs_, 0));
    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
            for (size_t n = 0; n < numMobileTxLocs_; n++) {
                reshaped_3[i * (numBands_ * numTimeSteps_) + j_k][n] = rREqConst2AuxiliaryConstM[i * numMobileTxLocs_ + n][j_k];
            }
        }
    }

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

}