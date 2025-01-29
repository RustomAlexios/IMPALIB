# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib, zero_value, deepcopy, os
        
class InequalityConstraintMOBARP:
    def __init__(self, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_capac_constraints, mobile_capac_constraints, ALPHA, FILTERING_FLAG, NUM_MOBILE_TX_LOCS, NUM_RX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx):
        self.num_fixed_tx = NUM_FIXED_TX
        self.num_mobile_tx = NUM_MOBILE_TX
        self.num_bands = NUM_BANDS
        self.num_time_steps = NUM_TIME_STEPS
        self.fixed_capac_constraints = fixed_capac_constraints
        self.mobile_capac_constraints = mobile_capac_constraints
        self.alpha = ALPHA
        self.filtering_flag = FILTERING_FLAG
        self.num_mobile_tx_locs = NUM_MOBILE_TX_LOCS
        self.num_rx_locs = NUM_RX_LOCS
        self.connectivity_mobile_tx = connectivity_mobile_tx
        self.connectivity_fixed_tx = connectivity_fixed_tx
        self.conx_fixed_tx_per_num_rx_locs =  self.connectivity_fixed_tx.reshape(-1, self.num_rx_locs)
        self.conx_mob_tx_rx = np.concatenate([np.transpose(self.connectivity_mobile_tx, (2,3,0,1))]*self.num_mobile_tx, axis=3)

    def ineq_capac_const_update(self, fixed_x_eq_const_to_fixed_capac_const_m, mobile_x_eq_const_to_mobile_capac_const_m):
        #calculate messages from capacity constraints to fixed & mobile x
        fixed_capac_const_to_fixed_x_eq_const_m = np.zeros((self.num_fixed_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
        mobile_capac_const_to_mobile_x_eq_const_m = np.zeros((self.num_mobile_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
        remaining_indices = np.arange(self.num_bands)
        for j in range(self.num_bands):
            mask = remaining_indices != j
            remaining_fixed_messages = fixed_x_eq_const_to_fixed_capac_const_m[:, mask, :]
            optimum_fixed_m = (np.sort(remaining_fixed_messages, axis=1))[np.arange(self.num_fixed_tx),self.fixed_capac_constraints-1,:]
            fixed_capac_const_to_fixed_x_eq_const_m[:, j, :] = np.maximum(-optimum_fixed_m, zero_value)
            
            remaining_mobile_messages = mobile_x_eq_const_to_mobile_capac_const_m[:, mask, :]
            optimum_mobile_m = (np.sort(remaining_mobile_messages, axis=1))[np.arange(self.num_mobile_tx),self.mobile_capac_constraints-1,:]
            mobile_capac_const_to_mobile_x_eq_const_m[:, j, :] = np.maximum(-optimum_mobile_m, zero_value)

        self.fixed_capac_const_to_fixed_x_eq_const_m_dummy = deepcopy(fixed_capac_const_to_fixed_x_eq_const_m)
        self.mobile_capac_const_to_mobile_x_eq_const_m_dummy = deepcopy(mobile_capac_const_to_mobile_x_eq_const_m)
        
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rFixedCapacConst2FixedXEqConstDummyM_"
        # rFixedCapacConst2FixedXEqConstDummyM_ = self.fixed_capac_const_to_fixed_x_eq_const_m_dummy.astype(np_impa_lib)
        # np.save(f_input_1, rFixedCapacConst2FixedXEqConstDummyM_.flatten())
    
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rMobileCapacConst2MobileXEqConstDummyM_"
        # rMobileCapacConst2MobileXEqConstDummyM_ = self.mobile_capac_const_to_mobile_x_eq_const_m_dummy.astype(np_impa_lib)
        # np.save(f_input_2, rMobileCapacConst2MobileXEqConstDummyM_.flatten())
        
        # exit()
        
    def process_filtering_capac_const(self, iter):
        #perform fitering on messages from capacity constraints to fixed & mobile x
        alpha = self.alpha
        filtering_flag = self.filtering_flag

        if iter == 0 and filtering_flag and alpha != zero_value:
            fixed_capac_const_to_fixed_x_eq_const_m = (1 - alpha) * self.fixed_capac_const_to_fixed_x_eq_const_m_dummy
            self.fixed_capac_constr_to_fixed_x_m_old = deepcopy(fixed_capac_const_to_fixed_x_eq_const_m)
            
            mobile_capac_const_to_mobile_x_eq_const_m = (1 - alpha) * self.mobile_capac_const_to_mobile_x_eq_const_m_dummy
            self.mobile_capac_constr_to_mobile_x_m_old = deepcopy(mobile_capac_const_to_mobile_x_eq_const_m)
            
        elif iter > 0 and filtering_flag and alpha != 0:
            fixed_capac_const_to_fixed_x_eq_const_m = alpha * self.fixed_capac_constr_to_fixed_x_m_old + (1 - alpha) * self.fixed_capac_const_to_fixed_x_eq_const_m_dummy
            self.fixed_capac_constr_to_fixed_x_m_old = deepcopy(fixed_capac_const_to_fixed_x_eq_const_m)
            
            mobile_capac_const_to_mobile_x_eq_const_m = alpha * self.mobile_capac_constr_to_mobile_x_m_old + (1 - alpha) * self.mobile_capac_const_to_mobile_x_eq_const_m_dummy
            self.mobile_capac_constr_to_mobile_x_m_old = deepcopy(mobile_capac_const_to_mobile_x_eq_const_m)
            
        elif not filtering_flag or alpha==0:
            fixed_capac_const_to_fixed_x_eq_const_m = deepcopy(self.fixed_capac_const_to_fixed_x_eq_const_m_dummy)
            mobile_capac_const_to_mobile_x_eq_const_m = deepcopy(self.mobile_capac_const_to_mobile_x_eq_const_m_dummy)

        self.fixed_capac_const_to_fixed_x_eq_const_m = fixed_capac_const_to_fixed_x_eq_const_m     
        self.mobile_capac_const_to_mobile_x_eq_const_m = mobile_capac_const_to_mobile_x_eq_const_m
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rFixedCapacConst2FixedXEqConstM_"
        # rFixedCapacConst2FixedXEqConstM_ = self.fixed_capac_const_to_fixed_x_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rFixedCapacConst2FixedXEqConstM_.flatten())
    
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rMobileCapacConst2MobileXEqConstM_"
        # rMobileCapacConst2MobileXEqConstM_ = self.mobile_capac_const_to_mobile_x_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_2, rMobileCapacConst2MobileXEqConstM_.flatten())
        
    def mobile_loc_eq_const_to_r_eq_const_update(self, r_eq_const_to_mobile_loc_eq_const_m):
        #calculate messages from mobile-loc equality constraint to r equality constraint
        mobile_loc_eq_const_to_r_eq_const_m = np.zeros_like(r_eq_const_to_mobile_loc_eq_const_m, dtype = np_impa_lib)
        
        remaining_indices = np.arange(self.num_mobile_tx_locs)
        for index in range(self.num_mobile_tx_locs):
            mask = remaining_indices != index
            remaining_msgs = r_eq_const_to_mobile_loc_eq_const_m[:, mask]
            mobile_loc_eq_const_to_r_eq_const_m[:, index] = -np.min(remaining_msgs, axis=1)
        self.mobile_loc_eq_const_to_r_eq_const_m_dummy = deepcopy(mobile_loc_eq_const_to_r_eq_const_m)
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rMobileLocEqConst2REqConstDummyM_"
        # rMobileLocEqConst2REqConstDummyM_ = self.mobile_loc_eq_const_to_r_eq_const_m_dummy.astype(np_impa_lib)
        # np.save(f_input_1, rMobileLocEqConst2REqConstDummyM_.flatten())

    def process_filtering_mobile_loc_eq(self, iter):
        #perform messages on messages from mobile-loc equality constraint to r equality constraint
        alpha = self.alpha
        filtering_flag = self.filtering_flag

        if iter == 0 and filtering_flag and alpha != zero_value:
            mobile_loc_eq_const_to_r_eq_const_m = (1 - alpha) * self.mobile_loc_eq_const_to_r_eq_const_m_dummy
            self.mobile_loc_eq_const_to_r_eq_const_m_old = deepcopy(mobile_loc_eq_const_to_r_eq_const_m)
            
        elif iter > 0 and filtering_flag and alpha != 0:
            mobile_loc_eq_const_to_r_eq_const_m = alpha * self.mobile_loc_eq_const_to_r_eq_const_m_old + (1 - alpha) * self.mobile_loc_eq_const_to_r_eq_const_m_dummy
            self.mobile_loc_eq_const_to_r_eq_const_m_old = deepcopy(mobile_loc_eq_const_to_r_eq_const_m)
            
        elif not filtering_flag or alpha==0:
            mobile_loc_eq_const_to_r_eq_const_m = deepcopy(self.mobile_loc_eq_const_to_r_eq_const_m_dummy)

        self.mobile_loc_eq_const_to_r_eq_const_m = mobile_loc_eq_const_to_r_eq_const_m  
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rMobileLocEqConst2REqConstM_"
        # rMobileLocEqConst2REqConstM_ = self.mobile_loc_eq_const_to_r_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rMobileLocEqConst2REqConstM_.flatten())
        
    def set_cover_ineq_const_update(self, z_eq_const_to_set_cover_ineq_const_m, fixed_x_eq_const_to_set_cover_const_m):
        #calculate messages from set cover inequality constraint to fixed and mobile x
        reshaped_fixed_x_eq_const_to_set_cover_const_m = fixed_x_eq_const_to_set_cover_const_m.reshape(self.num_fixed_tx*self.num_bands, -1)
        reshaped_z_eq_const_to_set_cover_ineq_const_m = z_eq_const_to_set_cover_ineq_const_m.reshape(self.num_bands*self.num_mobile_tx,self.num_time_steps,self.num_mobile_tx_locs,self.num_rx_locs).transpose()
        
        temp_reshaped_connectivity_fixed_tx = self.conx_fixed_tx_per_num_rx_locs.reshape(self.num_fixed_tx*self.num_bands, -1)
        
        set_cover_ineq_const_to_fixed_x_eq_const_m = np.zeros(reshaped_fixed_x_eq_const_to_set_cover_const_m.shape, dtype=np_impa_lib)
        
        temp_conx_mob_tx_rx = self.conx_mob_tx_rx.transpose(1,2,0,3)

        remaining_indices = np.arange(self.num_fixed_tx*self.num_bands)
        for index in range(self.num_fixed_tx*self.num_bands):
            mask = remaining_indices != index
            min_remaining_fixed_msgs = np.min(reshaped_fixed_x_eq_const_to_set_cover_const_m[mask, :], axis=0, where = temp_reshaped_connectivity_fixed_tx[mask, :]==1, initial=np.inf)
            min_mobile_msgs = np.min(reshaped_z_eq_const_to_set_cover_ineq_const_m, axis=(3,1), where = temp_conx_mob_tx_rx==1, initial=np.inf).T.flatten()
            set_cover_ineq_const_to_fixed_x_eq_const_m[index, :] = np.where(temp_reshaped_connectivity_fixed_tx[index, :], -np.maximum(zero_value, np.minimum(min_remaining_fixed_msgs, min_mobile_msgs)), set_cover_ineq_const_to_fixed_x_eq_const_m[index, :])
        
        self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy = set_cover_ineq_const_to_fixed_x_eq_const_m.T.reshape(self.num_time_steps, self.num_rx_locs, self.num_fixed_tx*self.num_bands)
        
        # for k in range(self.num_time_steps):
        #     print(f"Time Step: {k}")
        #     for l in range(self.num_rx_locs):
        #         print(f"  Rx Location: {l}")
        #         for i_j in range(self.num_fixed_tx * self.num_bands):
        #             print(f"{self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy[k][l][i_j]:10.4f}", end=" ")
        #         print()
        
        temp_reshaped_z_eq_const_to_set_cover_ineq_const_m = np.swapaxes(reshaped_z_eq_const_to_set_cover_ineq_const_m, 0, 3).reshape(-1,self.num_time_steps*self.num_rx_locs)
        
        temp_reshaped_conx_mob_tx_rx = np.swapaxes(temp_conx_mob_tx_rx, 0, 3).reshape(-1,self.num_time_steps*self.num_rx_locs)

        set_cover_ineq_const_to_z_eq_const_m = np.zeros(temp_reshaped_z_eq_const_to_set_cover_ineq_const_m.shape, dtype=np_impa_lib)
        
        min_fixed_msgs = np.min(reshaped_fixed_x_eq_const_to_set_cover_const_m, axis=0, where=temp_reshaped_connectivity_fixed_tx==1, initial=np.inf)

        remaining_indices = np.arange(self.num_bands*self.num_mobile_tx*self.num_mobile_tx_locs)
        for index in range(self.num_bands*self.num_mobile_tx*self.num_mobile_tx_locs):
            mask = remaining_indices != index
            min_remaining_mobile_msgs = np.min(temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[mask, :], axis=0, where = temp_reshaped_conx_mob_tx_rx[mask, :]==1, initial=np.inf)
            set_cover_ineq_const_to_z_eq_const_m[index, :] = np.where(temp_reshaped_conx_mob_tx_rx[index, :], -np.maximum(zero_value, np.minimum(min_remaining_mobile_msgs, min_fixed_msgs)), set_cover_ineq_const_to_z_eq_const_m[index, :])
        
        reshaped_set_cover_ineq_const_to_z_eq_const_m = np.zeros((self.num_time_steps*self.num_rx_locs, self.num_mobile_tx_locs, self.num_bands*self.num_mobile_tx))
        temp_reshaped_set_cover_ineq_const_to_z_eq_const_m = set_cover_ineq_const_to_z_eq_const_m.T
        for index in range(self.num_time_steps*self.num_rx_locs):
            for index_n in range(self.num_mobile_tx_locs):
                arr = temp_reshaped_set_cover_ineq_const_to_z_eq_const_m[index, index_n::self.num_mobile_tx_locs]
                reshaped_set_cover_ineq_const_to_z_eq_const_m[index, index_n] = arr
        
        temp_set_cover_ineq_const_to_z_eq_const_m = reshaped_set_cover_ineq_const_to_z_eq_const_m.reshape(self.num_time_steps, self.num_rx_locs, self.num_mobile_tx_locs, self.num_bands*self.num_mobile_tx)
        # temp_set_cover_ineq_const_to_z_eq_const_m = temp_set_cover_ineq_const_to_z_eq_const_m
        self.set_cover_ineq_const_to_z_eq_const_m_dummy = temp_set_cover_ineq_const_to_z_eq_const_m

        # for k in range(self.num_time_steps):
        #     print(f"Time Step: {k}")
        #     for l in range(self.num_rx_locs):
        #         print(f"  Rx Location: {l}")
        #         for n in range(self.num_mobile_tx_locs):
        #             print(f"    Mobile Tx Location: {n} -> ", end="")
        #             for j_i in range(self.num_bands * self.num_mobile_tx):
        #                 print(f"{self.set_cover_ineq_const_to_z_eq_const_m_dummy[k][l][n][j_i]:12.3f}", end=" ")
        #             print()
        
        # exit()      
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rSetCoverIneqConst2ZEqConstDummyM_"
        # rSetCoverIneqConst2ZEqConstDummyM_ = self.set_cover_ineq_const_to_z_eq_const_m_dummy.astype(np_impa_lib)
        # np.save(f_input_1, rSetCoverIneqConst2ZEqConstDummyM_.flatten())
        
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rSetCoverIneqConst2FixedXEqConstDummyM_"
        # rSetCoverIneqConst2FixedXEqConstDummyM_ = self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy.astype(np_impa_lib)
        # np.save(f_input_2, rSetCoverIneqConst2FixedXEqConstDummyM_.flatten())
            
    def process_filtering_set_cover_const(self, iter):
        #perform filtering on messages from set cover inequality constraint to fixed and mobile x
        alpha = self.alpha
        filtering_flag = self.filtering_flag
        
        if iter == 0 and filtering_flag and alpha != zero_value:
            set_cover_ineq_const_to_fixed_x_eq_const_m = (1 - alpha) * self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy
            self.set_cover_ineq_const_to_fixed_x_eq_const_m_old = deepcopy(set_cover_ineq_const_to_fixed_x_eq_const_m)
            
            set_cover_ineq_const_to_z_eq_const_m = (1 - alpha) * self.set_cover_ineq_const_to_z_eq_const_m_dummy
            self.set_cover_ineq_const_to_z_eq_const_m_old = deepcopy(set_cover_ineq_const_to_z_eq_const_m)
            
        elif iter > 0 and filtering_flag and alpha != 0:
            set_cover_ineq_const_to_fixed_x_eq_const_m = alpha * self.set_cover_ineq_const_to_fixed_x_eq_const_m_old + (1 - alpha) * self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy
            self.set_cover_ineq_const_to_fixed_x_eq_const_m_old = deepcopy(set_cover_ineq_const_to_fixed_x_eq_const_m)
            
            set_cover_ineq_const_to_z_eq_const_m = alpha * self.set_cover_ineq_const_to_z_eq_const_m_old + (1 - alpha) * self.set_cover_ineq_const_to_z_eq_const_m_dummy
            self.set_cover_ineq_const_to_z_eq_const_m_old = deepcopy(set_cover_ineq_const_to_z_eq_const_m)
            
        elif not filtering_flag or alpha==0:
            set_cover_ineq_const_to_fixed_x_eq_const_m = deepcopy(self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy)
            set_cover_ineq_const_to_z_eq_const_m = deepcopy(self.set_cover_ineq_const_to_z_eq_const_m_dummy)

        self.set_cover_ineq_const_to_fixed_x_eq_const_m = set_cover_ineq_const_to_fixed_x_eq_const_m     
        self.set_cover_ineq_const_to_z_eq_const_m = set_cover_ineq_const_to_z_eq_const_m
        
        # f_input_1 = os.getcwd() + "/../../../src/impa/ut_results/rSetCoverIneqConst2ZEqConstM_"
        # rSetCoverIneqConst2ZEqConstM_ = self.set_cover_ineq_const_to_z_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_1, rSetCoverIneqConst2ZEqConstM_.flatten())
        
        # f_input_2 = os.getcwd() + "/../../../src/impa/ut_results/rSetCoverIneqConst2FixedXEqConstM_"
        # rSetCoverIneqConst2FixedXEqConstM_ = self.set_cover_ineq_const_to_fixed_x_eq_const_m.astype(np_impa_lib)
        # np.save(f_input_2, rSetCoverIneqConst2FixedXEqConstM_.flatten())
        
        # for k in range(self.num_time_steps):
        #     print(f"Time Step: {k}")
        #     for l in range(self.num_rx_locs):
        #         print(f"  Rx Location: {l}")
        #         for i_j in range(self.num_fixed_tx * self.num_bands):
        #             print(f"{self.set_cover_ineq_const_to_fixed_x_eq_const_m[k][l][i_j]:10.4f}", end=" ")
        #         print()
        # exit()
        
        # for k in range(self.num_time_steps):
        #     print(f"Time Step: {k}")
        #     for l in range(self.num_rx_locs):
        #         print(f"  Rx Location: {l}")
        #         for n in range(self.num_mobile_tx_locs):
        #             print(f"    Mobile Tx Location: {n} -> ", end="")
        #             for j_i in range(self.num_bands * self.num_mobile_tx):
        #                 print(f"{self.set_cover_ineq_const_to_z_eq_const_m[k][l][n][j_i]:12.3f}", end=" ")
        #             print()
        
        # exit()
            