# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib, os
            
class EqualityConstraintMOBARP:
    def __init__(self, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_x_costs, mobile_x_costs, NUM_MOBILE_TX_LOCS, z_costs, NUM_RX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx, r_costs, EXCLUDE_CAP_FLAG):
        self.num_fixed_tx = NUM_FIXED_TX
        self.num_mobile_tx = NUM_MOBILE_TX
        self.num_bands = NUM_BANDS
        self.num_time_steps = NUM_TIME_STEPS
        self.fixed_x_costs = fixed_x_costs
        self.mobile_x_costs = mobile_x_costs
        self.num_mobile_tx_locs = NUM_MOBILE_TX_LOCS
        self.z_costs = z_costs
        self.num_rx_locs = NUM_RX_LOCS
        self.connectivity_mobile_tx = connectivity_mobile_tx
        self.connectivity_fixed_tx = connectivity_fixed_tx
        self.conx_mob_tx_per_num_mob_tx_locs = np.vstack([np.clip(np.sum(self.connectivity_mobile_tx, axis=3).transpose(1,2,0).reshape(-1, self.num_mobile_tx_locs), 0, 1)] * self.num_mobile_tx)
        self.conx_fixed_tx_per_num_rx_locs =  self.connectivity_fixed_tx.reshape(-1, self.num_rx_locs)
        self.conx_mob_tx_rx = np.concatenate([np.transpose(self.connectivity_mobile_tx, (2,3,0,1))]*self.num_mobile_tx, axis=3)
        self.r_costs = r_costs
        self.exclude_cap_flag = EXCLUDE_CAP_FLAG

    def x_eq_const_to_auxiliary_and_set_cover_const_update(self, fixed_capac_const_to_fixed_x_eq_const_m, mobile_capac_const_to_mobile_x_eq_const_m, auxiliary_const_to_mobile_x_eq_const_m, set_cover_ineq_const_to_fixed_x_eq_const_m):
        #calculate messages from x to auxiliary constraints and set cover constraints
        conx_mob_tx_per_num_mob_tx_locs = self.conx_mob_tx_per_num_mob_tx_locs

        conx_columns = {col: np.where(conx_mob_tx_per_num_mob_tx_locs[:, col] == 1)[0] for col in range(self.num_mobile_tx_locs)}

        mobile_x_eq_const_to_auxiliary_const_m = np.zeros((self.mobile_x_costs.size, self.num_mobile_tx_locs), dtype=np_impa_lib)
        fixed_x_eq_const_to_set_cover_const_m = np.zeros((self.fixed_x_costs.size, self.num_rx_locs), dtype=np_impa_lib)
        reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m = set_cover_ineq_const_to_fixed_x_eq_const_m.transpose(2, 0, 1).reshape(-1, self.num_rx_locs)
        
        sums_mobile_x_eq_const = np.zeros_like(self.mobile_x_costs.flatten(), dtype = np_impa_lib)
        conx_rows = np.where(np.sum(conx_mob_tx_per_num_mob_tx_locs, axis=1) !=0)[0]
        
        if (not self.exclude_cap_flag):
            sums_mobile_x_eq_const[conx_rows] = np.sum(auxiliary_const_to_mobile_x_eq_const_m, axis=1, where= conx_mob_tx_per_num_mob_tx_locs==1)[conx_rows] + mobile_capac_const_to_mobile_x_eq_const_m.flatten()[conx_rows] + self.mobile_x_costs.flatten()[conx_rows]
        else:
            sums_mobile_x_eq_const[conx_rows] = np.sum(auxiliary_const_to_mobile_x_eq_const_m, axis=1, where= conx_mob_tx_per_num_mob_tx_locs==1)[conx_rows] + self.mobile_x_costs.flatten()[conx_rows]

        for index in range(self.num_mobile_tx_locs):
            mobile_x_eq_const_to_auxiliary_const_m[conx_columns[index], index] = sums_mobile_x_eq_const[conx_columns[index]] - auxiliary_const_to_mobile_x_eq_const_m[conx_columns[index], index]

        self.mobile_x_eq_const_to_auxiliary_const_m = mobile_x_eq_const_to_auxiliary_const_m

        conx_fixed_tx_per_num_rx_locs = self.conx_fixed_tx_per_num_rx_locs
        conx_fixed_tx_cols = {col: np.where(conx_fixed_tx_per_num_rx_locs[:, col] == 1)[0] for col in range(self.num_rx_locs)}
        conx_fixed_tx_rows = np.where(np.sum(conx_fixed_tx_per_num_rx_locs, axis=1) !=0)[0]
        sums_fixed_x_eq_const = np.zeros_like(self.fixed_x_costs.flatten(), dtype = np_impa_lib)
        
        if (not self.exclude_cap_flag):
            sums_fixed_x_eq_const[conx_fixed_tx_rows] = fixed_capac_const_to_fixed_x_eq_const_m.flatten()[conx_fixed_tx_rows] + self.fixed_x_costs.flatten()[conx_fixed_tx_rows] + np.sum(reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m, axis=1, where = conx_fixed_tx_per_num_rx_locs==1)[conx_fixed_tx_rows]
        else:
            sums_fixed_x_eq_const[conx_fixed_tx_rows] = self.fixed_x_costs.flatten()[conx_fixed_tx_rows] + np.sum(reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m, axis=1, where = conx_fixed_tx_per_num_rx_locs==1)[conx_fixed_tx_rows]
        
        for index in range(self.num_rx_locs):
            fixed_x_eq_const_to_set_cover_const_m[conx_fixed_tx_cols[index], index] = sums_fixed_x_eq_const[conx_fixed_tx_cols[index]] - reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m[conx_fixed_tx_cols[index], index]
        
        self.fixed_x_eq_const_to_set_cover_const_m = fixed_x_eq_const_to_set_cover_const_m
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rFixedXEqConst2SetCoverConstM_"
        # rFixedXEqConst2SetCoverConstM_ = self.fixed_x_eq_const_to_set_cover_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rFixedXEqConst2SetCoverConstM_.flatten())
    
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rMobileXEqConst2AuxiliaryConstM_"
        # rMobileXEqConst2AuxiliaryConstM_ = self.mobile_x_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        # np.save(f_input_2, rMobileXEqConst2AuxiliaryConstM_.flatten())


    def r_eq_const_activation(self, auxiliary_const_to_r_eq_const_m, mobile_loc_eq_const_to_r_eq_const_m):
        #activate r equality constraint: calculate messages from r equality constraint to auxiliary constraint and mobile-loc equality constraint
        reshaped_auxiliary_const_to_r_eq_const_m = auxiliary_const_to_r_eq_const_m.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)
        conx_mob_tx_r = self.conx_mob_tx_per_num_mob_tx_locs.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)

        conx_rows = np.where(np.sum(conx_mob_tx_r, axis=1) !=0)[0]
        sums_r_eq_const = np.zeros_like(self.r_costs.flatten(), dtype = np_impa_lib)
        sums_r_eq_const[conx_rows] = np.sum(reshaped_auxiliary_const_to_r_eq_const_m, axis=1, where = conx_mob_tx_r==1)[conx_rows] + \
                            self.r_costs.flatten()[conx_rows] + mobile_loc_eq_const_to_r_eq_const_m.flatten()[conx_rows]

        conx_cols = {col: np.where(conx_mob_tx_r[:, col] == 1)[0] for col in range(self.num_bands*self.num_time_steps)}

        r_eq_const_to_auxiliary_const_m = np.zeros(reshaped_auxiliary_const_to_r_eq_const_m.shape, dtype=np_impa_lib)

        for index in range(self.num_bands*self.num_time_steps):
            r_eq_const_to_auxiliary_const_m[conx_cols[index], index] = sums_r_eq_const[conx_cols[index]] - reshaped_auxiliary_const_to_r_eq_const_m[conx_cols[index], index]

        self.r_eq_const_to_auxiliary_const_m = r_eq_const_to_auxiliary_const_m

        r_eq_const_to_mobile_loc_eq_const_m = np.sum(reshaped_auxiliary_const_to_r_eq_const_m, axis=1, where=conx_mob_tx_r==1) + self.r_costs.flatten()
        
        self.r_eq_const_to_mobile_loc_eq_const_m = r_eq_const_to_mobile_loc_eq_const_m.reshape(self.num_mobile_tx, -1)
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rREqConst2AuxiliaryConstM_"
        # rREqConst2AuxiliaryConstM_ = self.r_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rREqConst2AuxiliaryConstM_.flatten())
        
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rREqConst2MobileLocEqConstM_"
        # rREqConst2MobileLocEqConstM_ = self.r_eq_const_to_mobile_loc_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_2, rREqConst2MobileLocEqConstM_.flatten())
        
    def z_eq_const_to_auxiliary_const_update(self, set_cover_ineq_const_to_z_eq_const_m):   
        #calculate messages from z equality constraint to auxiliary constraint
        reshaped_set_cover_ineq_const_to_z_eq_const_m = np.reshape(np.transpose(np.sum(set_cover_ineq_const_to_z_eq_const_m, axis=1, where=self.conx_mob_tx_rx==1).reshape(-1, set_cover_ineq_const_to_z_eq_const_m.shape[-1])), (self.mobile_x_costs.size,self.num_mobile_tx_locs))
        z_eq_const_to_auxiliary_const_m = np.zeros_like(reshaped_set_cover_ineq_const_to_z_eq_const_m, dtype = np_impa_lib)
        mask = self.conx_mob_tx_per_num_mob_tx_locs == 1
        z_eq_const_to_auxiliary_const_m[mask] = self.z_costs[mask] + reshaped_set_cover_ineq_const_to_z_eq_const_m[mask]
        self.z_eq_const_to_auxiliary_const_m = z_eq_const_to_auxiliary_const_m
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rZEqConst2AuxiliaryConstM_"
        # rZEqConst2AuxiliaryConstM_ = self.z_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rZEqConst2AuxiliaryConstM_.flatten())
        
    
    def x_eq_const_activation(self, auxiliary_const_to_mobile_x_eq_const_m, set_cover_ineq_const_to_fixed_x_eq_const_m):
        #calculate messages from fixed and mobile x to capacity constraints 
        conx_rows_mobile = np.where(np.clip(np.sum(self.conx_mob_tx_per_num_mob_tx_locs, axis=1), 0, 1) ==1)[0]
        mobile_x_eq_const_to_mobile_capac_const_m = np.zeros_like(self.mobile_x_costs.flatten(), dtype=np_impa_lib)

        mobile_x_eq_const_to_mobile_capac_const_m[conx_rows_mobile] = np.sum(auxiliary_const_to_mobile_x_eq_const_m, axis=1, where= self.conx_mob_tx_per_num_mob_tx_locs==1)[conx_rows_mobile] + \
                                            self.mobile_x_costs.flatten()[conx_rows_mobile]
        reshaped_mobile_x_eq_const_to_mobile_capac_const_m = mobile_x_eq_const_to_mobile_capac_const_m.reshape(self.num_mobile_tx, self.num_bands, self.num_time_steps)

        reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m = set_cover_ineq_const_to_fixed_x_eq_const_m.transpose(2, 0, 1).reshape(-1, self.num_rx_locs)
        conx_rows_cols = np.where(np.clip(np.sum(self.conx_fixed_tx_per_num_rx_locs, axis=1), 0, 1) ==1)[0]
        fixed_x_eq_const_to_fixed_capac_const_m = np.zeros_like(self.fixed_x_costs.flatten(), dtype = np_impa_lib)
        fixed_x_eq_const_to_fixed_capac_const_m[conx_rows_cols] = np.sum(reshaped_ineq_set_cover_const_to_fixed_x_eq_const_m, axis=1, where = self.conx_fixed_tx_per_num_rx_locs==1)[conx_rows_cols] + self.fixed_x_costs.flatten()[conx_rows_cols]
        reshaped_fixed_x_eq_const_to_fixed_capac_const_m = fixed_x_eq_const_to_fixed_capac_const_m.reshape(self.num_fixed_tx, self.num_bands, self.num_time_steps)
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rFixedXEqConst2FixedCapacConstM_"
        # rFixedXEqConst2FixedCapacConstM_ = reshaped_fixed_x_eq_const_to_fixed_capac_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rFixedXEqConst2FixedCapacConstM_.flatten())
        
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rMobileXEqConst2MobileCapacConstM_"
        # rMobileXEqConst2MobileCapacConstM_ = reshaped_mobile_x_eq_const_to_mobile_capac_const_m.astype(np_impa_lib)
        # np.save(f_input_2, rMobileXEqConst2MobileCapacConstM_.flatten())
        
        return reshaped_mobile_x_eq_const_to_mobile_capac_const_m, reshaped_fixed_x_eq_const_to_fixed_capac_const_m

    def z_eq_const_to_set_cover_ineq_const_update(self, auxiliary_const_to_z_eq_const_m, set_cover_ineq_const_to_z_eq_const_m):
        #calculate messages from z equality constraint to set cover inequality constraint
        # reshaped_auxiliary_const_to_z_eq_const_m = auxiliary_const_to_z_eq_const_m.reshape(self.num_mobile_tx, self.num_mobile_tx_locs, -1).transpose(0, 2, 1).reshape(-1,self.num_mobile_tx_locs)
        temp_set_cover_ineq_const_to_z_eq_const_m = set_cover_ineq_const_to_z_eq_const_m.transpose(3, 1, 0, 2)
        reshaped_set_cover_ineq_const_to_z_eq_const_m = np.array([temp_set_cover_ineq_const_to_z_eq_const_m[:, i].flatten() for i in range(self.num_rx_locs)]).transpose()
        
        temp_connectivity_mobile_tx_rx = self.conx_mob_tx_rx.transpose(3, 1, 0, 2)
        reshaped_temp_connectivity_mobile_tx_rx = np.array([temp_connectivity_mobile_tx_rx[:, i].flatten() for i in range(self.num_rx_locs)]).transpose()

        z_eq_const_to_set_cover_ineq_const_m = np.zeros(reshaped_set_cover_ineq_const_to_z_eq_const_m.shape, dtype=np_impa_lib)
        
        conx_cols = {col: np.where(reshaped_temp_connectivity_mobile_tx_rx[:, col] == 1)[0] for col in range(self.num_rx_locs)}

        conx_elems = np.where(self.conx_mob_tx_per_num_mob_tx_locs.flatten()==1)[0]
        
        sums_z_eq_const = np.zeros(auxiliary_const_to_z_eq_const_m.flatten().shape, dtype=np_impa_lib)
        sums_z_eq_const[conx_elems] = np.sum(reshaped_set_cover_ineq_const_to_z_eq_const_m, axis=1, where=reshaped_temp_connectivity_mobile_tx_rx==1)[conx_elems] + (auxiliary_const_to_z_eq_const_m).flatten()[conx_elems] + \
                                        + (self.z_costs).flatten()[conx_elems]
        
        for index in range(self.num_rx_locs):
            z_eq_const_to_set_cover_ineq_const_m[conx_cols[index], index] = sums_z_eq_const[conx_cols[index]] - reshaped_set_cover_ineq_const_to_z_eq_const_m[conx_cols[index], index]
        
        self.z_eq_const_to_set_cover_ineq_const_m = z_eq_const_to_set_cover_ineq_const_m
        
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rZEqConst2SetCoverIneqConstM_"
        # rZEqConst2SetCoverIneqConstM_ = self.z_eq_const_to_set_cover_ineq_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rZEqConst2SetCoverIneqConstM_.flatten())
        

        
        