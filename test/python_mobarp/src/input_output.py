# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib
import time

class OutputsMOBARP:
    def __init__(self, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_x_costs, mobile_x_costs, r_costs, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx, EXCLUDE_CAP_FLAG):
        self.num_fixed_tx = NUM_FIXED_TX
        self.num_mobile_tx = NUM_MOBILE_TX
        self.num_bands = NUM_BANDS
        self.num_time_steps = NUM_TIME_STEPS
        self.fixed_x_costs = fixed_x_costs
        self.mobile_x_costs = mobile_x_costs
        self.r_costs = r_costs
        self.num_rx_locs = NUM_RX_LOCS
        self.num_mobile_tx_locs = NUM_MOBILE_TX_LOCS
        self.connectivity_mobile_tx = connectivity_mobile_tx
        self.connectivity_fixed_tx = connectivity_fixed_tx
        self.conx_mob_tx_per_num_mob_tx_locs = np.vstack([np.clip(np.sum(self.connectivity_mobile_tx, axis=3).transpose(1,2,0).reshape(-1, self.num_mobile_tx_locs), 0, 1)] * self.num_mobile_tx)
        self.conx_fixed_tx_per_num_rx_locs =  self.connectivity_fixed_tx.reshape(-1, self.num_rx_locs)
        self.conx_mob_tx_rx = np.concatenate([np.transpose(self.connectivity_mobile_tx, (2,3,0,1))]*self.num_mobile_tx, axis=3)
        self.exclude_cap_flag = EXCLUDE_CAP_FLAG

        
    def extrinsic_update(self, fixed_capac_const_to_fixed_x_eq_const_m, mobile_capac_const_to_mobile_x_eq_const_m, 
                            auxiliary_const_to_mobile_x_eq_const_m, set_cover_ineq_const_to_fixed_x_eq_const_m, 
                            auxiliary_const_to_r_eq_const_m, mobile_loc_eq_const_to_r_eq_const_m, 
                            auxiliary_const_to_z_eq_const_m, set_cover_ineq_const_to_z_eq_const_m):

        #calculate extrinsic messages to get hard decision
        conx_rows_mobile = np.where(np.clip(np.sum(self.conx_mob_tx_per_num_mob_tx_locs, axis=1), 0, 1) ==1)[0]
        extrinsic_mobile_x = np.zeros_like(self.mobile_x_costs.flatten(), dtype = np_impa_lib)
        if (not self.exclude_cap_flag):
            extrinsic_mobile_x[conx_rows_mobile] = np.sum(auxiliary_const_to_mobile_x_eq_const_m, axis=1, where = self.conx_mob_tx_per_num_mob_tx_locs==1)[conx_rows_mobile] + mobile_capac_const_to_mobile_x_eq_const_m.flatten()[conx_rows_mobile]
        else:
            extrinsic_mobile_x[conx_rows_mobile] = np.sum(auxiliary_const_to_mobile_x_eq_const_m, axis=1, where = self.conx_mob_tx_per_num_mob_tx_locs==1)[conx_rows_mobile]

        conx_rows_fixed = np.where(np.clip(np.sum(self.conx_fixed_tx_per_num_rx_locs, axis=1), 0, 1) ==1)[0]
        extrinsic_fixed_x = np.zeros_like(self.fixed_x_costs.flatten(), dtype = np_impa_lib)
        reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m = set_cover_ineq_const_to_fixed_x_eq_const_m.transpose(2, 0, 1).reshape(-1, self.num_rx_locs)
        
        if (not self.exclude_cap_flag):
            extrinsic_fixed_x[conx_rows_fixed] = fixed_capac_const_to_fixed_x_eq_const_m.flatten()[conx_rows_fixed] + np.sum(reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m, axis=1, where = self.conx_fixed_tx_per_num_rx_locs==1)[conx_rows_fixed]
        else:
            extrinsic_fixed_x[conx_rows_fixed] = np.sum(reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m, axis=1, where = self.conx_fixed_tx_per_num_rx_locs==1)[conx_rows_fixed]

        reshaped_auxiliary_const_to_r_eq_const_m = auxiliary_const_to_r_eq_const_m.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)

        conx_mob_tx_per_num_mob_tx_locs = self.conx_mob_tx_per_num_mob_tx_locs.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)
        extrinsic_r = np.sum(reshaped_auxiliary_const_to_r_eq_const_m, axis=1, where = conx_mob_tx_per_num_mob_tx_locs==1) + mobile_loc_eq_const_to_r_eq_const_m.flatten()
        
        reshaped_set_cover_ineq_const_to_z_eq_const_m = np.reshape(np.transpose(np.sum(set_cover_ineq_const_to_z_eq_const_m, axis=1, where=self.conx_mob_tx_rx==1).reshape(-1, set_cover_ineq_const_to_z_eq_const_m.shape[-1])), (self.mobile_x_costs.size,self.num_mobile_tx_locs))
        
        extrinsic_z = reshaped_set_cover_ineq_const_to_z_eq_const_m + auxiliary_const_to_z_eq_const_m
        return extrinsic_fixed_x, extrinsic_mobile_x, extrinsic_r, extrinsic_z