# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib, zero_value, deepcopy, os

class AuxiliaryConstraintMOBARP:
    def __init__(self, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS, ALPHA, FILTERING_FLAG, connectivity_mobile_tx):
        self.num_mobile_tx = NUM_MOBILE_TX
        self.num_bands = NUM_BANDS
        self.num_time_steps = NUM_TIME_STEPS
        self.num_mobile_tx_locs = NUM_MOBILE_TX_LOCS
        self.connectivity_mobile_tx = connectivity_mobile_tx
        self.alpha = ALPHA
        self.filtering_flag = FILTERING_FLAG
        self.connectivity_mobile_tx = connectivity_mobile_tx
        self.conx_mob_tx_per_num_mob_tx_locs = np.vstack([np.clip(np.sum(self.connectivity_mobile_tx, axis=3).transpose(1,2,0).reshape(-1, self.num_mobile_tx_locs), 0, 1)] * self.num_mobile_tx)

    def auxiliary_const_to_z_eq_const_update(self, mobile_x_eq_const_to_auxiliary_const_m, r_eq_const_to_auxiliary_const_m):
        #calculate messages from auxiliary constraint to z equality constraint
        auxiliary_const_to_z_eq_const_m = np.zeros(mobile_x_eq_const_to_auxiliary_const_m.shape, dtype=np_impa_lib)
        reshaped_r_eq_const_to_auxiliary_const_m = r_eq_const_to_auxiliary_const_m.reshape(self.num_mobile_tx, self.num_mobile_tx_locs, -1).transpose(0, 2, 1).reshape(-1,self.num_mobile_tx_locs)
        
        mask = self.conx_mob_tx_per_num_mob_tx_locs ==1
        sums_inputs = np.zeros_like(self.conx_mob_tx_per_num_mob_tx_locs, dtype = np_impa_lib)
        sums_inputs[mask] = mobile_x_eq_const_to_auxiliary_const_m[mask] + reshaped_r_eq_const_to_auxiliary_const_m[mask]
        self.auxiliary_const_to_z_eq_const_m = np.where(self.conx_mob_tx_per_num_mob_tx_locs, np.maximum(sums_inputs, np.maximum(reshaped_r_eq_const_to_auxiliary_const_m, mobile_x_eq_const_to_auxiliary_const_m)), auxiliary_const_to_z_eq_const_m)

    def auxiliary_const_to_r_eq_const_update(self, z_eq_const_to_auxiliary_const_m, mobile_x_eq_const_to_auxiliary_const_m):
        #calculate messages from auxiliary constraint to r equality constraint
        mask = self.conx_mob_tx_per_num_mob_tx_locs ==1
        auxiliary_const_to_r_eq_const_m = np.zeros(z_eq_const_to_auxiliary_const_m.shape, dtype = np_impa_lib)
        sums_inputs = np.zeros_like(self.conx_mob_tx_per_num_mob_tx_locs, dtype = np_impa_lib)
        sums_inputs[mask] = z_eq_const_to_auxiliary_const_m[mask] + mobile_x_eq_const_to_auxiliary_const_m[mask]
        self.auxiliary_const_to_r_eq_const_m = np.where(self.conx_mob_tx_per_num_mob_tx_locs, np.minimum(np.maximum(zero_value, -mobile_x_eq_const_to_auxiliary_const_m), np.maximum(z_eq_const_to_auxiliary_const_m, sums_inputs)), auxiliary_const_to_r_eq_const_m)

    def auxiliary_const_to_mobile_x_eq_const_update(self, z_eq_const_to_auxiliary_const_m, r_eq_const_to_auxiliary_const_m):
        #calculate messages from auxiliary constraint to mobile x equality constraint
        reshaped_r_eq_const_to_auxiliary_const_m = r_eq_const_to_auxiliary_const_m.reshape(self.num_mobile_tx, self.num_mobile_tx_locs, -1).transpose(0, 2, 1).reshape(-1,self.num_mobile_tx_locs)
        mask = self.conx_mob_tx_per_num_mob_tx_locs ==1
        auxiliary_const_to_mobile_x_eq_const_m = np.zeros(z_eq_const_to_auxiliary_const_m.shape, dtype = np_impa_lib)
        sums_inputs = np.zeros_like(self.conx_mob_tx_per_num_mob_tx_locs, dtype = np_impa_lib)
        sums_inputs[mask] = z_eq_const_to_auxiliary_const_m[mask] + reshaped_r_eq_const_to_auxiliary_const_m[mask]
        # print("----------------------------")
        # print(self.conx_mob_tx_per_num_mob_tx_locs)
        self.auxiliary_const_to_mobile_x_eq_const_m = np.where(self.conx_mob_tx_per_num_mob_tx_locs, np.minimum(np.maximum(zero_value, -reshaped_r_eq_const_to_auxiliary_const_m), np.maximum(z_eq_const_to_auxiliary_const_m, sums_inputs)), auxiliary_const_to_mobile_x_eq_const_m)