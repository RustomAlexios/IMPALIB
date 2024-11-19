# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib, deepcopy, math, itertools, combinations, defaultdict, Counter, product, time, pkl, os
from update_equality_constraint import EqualityConstraintMOBARP
from update_inequality_constraint import InequalityConstraintMOBARP
from input_output import OutputsMOBARP
from update_auxiliary_constraint import AuxiliaryConstraintMOBARP


class GraphicalModelMOBARP:
    def __init__(self, NUM_ITERATIONS, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, THRESHOLD, FILTERING_FLAG, ALPHA, RANDOM_TEST_FLAG, POST_PROCESS_FLAG, \
                        OVER_WRITE_CAP_FLAG, OVER_WRITE_CAP_VAL, EXCLUDE_CAP_FLAG, GET_SOL_APPROACH, CRITERIA_IM, PERCENTAGE_NEGATIVE_IM, OVERWRITE_IM):
        
        self.num_iterations = NUM_ITERATIONS
        self.num_fixed_tx = NUM_FIXED_TX
        self.num_mobile_tx = NUM_MOBILE_TX
        self.num_bands = NUM_BANDS
        self.num_time_steps = NUM_TIME_STEPS
        self.num_rx_locs = NUM_RX_LOCS
        self.num_mobile_tx_locs = NUM_MOBILE_TX_LOCS
        self.random_test_flag = RANDOM_TEST_FLAG
        self.post_process_flag = POST_PROCESS_FLAG
        self.alpha = ALPHA
        self.threshold = THRESHOLD
        self.filtering_flag = FILTERING_FLAG
        self.overwrite_cap_flag = OVER_WRITE_CAP_FLAG
        self.overwrite_cap_val = OVER_WRITE_CAP_VAL
        self.exclude_cap_flag = EXCLUDE_CAP_FLAG
        self.get_sol_approach = GET_SOL_APPROACH
        self.criteria_im = CRITERIA_IM
        self.percentage_neg_im = PERCENTAGE_NEGATIVE_IM
        self.overwrite_im = OVERWRITE_IM

    def initialize(self):

        input_load = self.input_load

        self.snr_threshold = -np.inf
        #retrieve parameters if not random
        if not self.random_test_flag:
            #read parameters
            print("Retrieving parameters")

            # if (input_load[10] != np.inf):
            #     self.snr_threshold = input_load[10]
            #     connectivity_fixed_tx = np.zeros(self.snr_fixed.shape, dtype=np.int64)
            #     connectivity_fixed_tx = np.where((self.snr_fixed >= self.snr_threshold), 1, connectivity_fixed_tx)
            #     connectivity_fixed_tx = np.transpose(connectivity_fixed_tx, (0, 3, 2, 1))
            #     connectivity_mobile_tx = np.zeros(self.snr_mobile.shape, dtype=np.int64)
            #     connectivity_mobile_tx = np.where((self.snr_mobile >= self.snr_threshold), 1, connectivity_mobile_tx)
            #     connectivity_mobile_tx = np.transpose(connectivity_mobile_tx, (0, 3, 2, 1))
            #     exit()
            # else:
            connectivity_fixed_tx = input_load[8]
            connectivity_mobile_tx = input_load[9]
            self.snr_threshold = input_load[10]
            
            self.connectivity_fixed_tx = connectivity_fixed_tx
            self.connectivity_mobile_tx = connectivity_mobile_tx

            assert connectivity_mobile_tx.shape[-1] == connectivity_fixed_tx.shape[-1], "Different Number of RX LOCS."
            assert connectivity_mobile_tx.shape[1] == connectivity_fixed_tx.shape[1], "Different Number of FREQS."
            assert connectivity_mobile_tx.shape[2] == connectivity_fixed_tx.shape[2], "Different Number of Time Steps."

            self.num_mobile_tx = input_load[1]
            self.num_mobile_tx_locs = input_load[5]
            
            assert self.num_mobile_tx_locs == connectivity_mobile_tx.shape[0], "Issue with NUM_MOBILE_TX_LOCS."

            self.num_fixed_tx, self.num_bands, self.num_time_steps, self.num_rx_locs = input_load[0], input_load[2], input_load[3], input_load[4]

            assert (self.num_fixed_tx, self.num_bands, self.num_time_steps, self.num_rx_locs) == connectivity_fixed_tx.shape

        #generate IM
        fixed_x_costs, mobile_x_costs, connectivity_fixed_tx, connectivity_mobile_tx, r_costs, z_costs = self.costs_connectivity_generation()
        
        if (not self.random_test_flag):
            self.exclude_cap_flag = input_load[11]
            connectivity_fixed_tx = self.connectivity_fixed_tx
            connectivity_mobile_tx = self.connectivity_mobile_tx
        else:
            self.connectivity_fixed_tx = connectivity_fixed_tx
            self.connectivity_mobile_tx = connectivity_mobile_tx
            
        #update some parameters if random or not
        #not random, include capacity
        if (not self.random_test_flag and not self.exclude_cap_flag):
            print("Reading Input, Including Capacity")
            # connectivity_fixed_tx = self.connectivity_fixed_tx
            # connectivity_mobile_tx = self.connectivity_mobile_tx
            fixed_capac_constraints = input_load[6]
            mobile_capac_constraints = input_load[7]
            assert len(fixed_capac_constraints) == self.num_fixed_tx, "Issue with Fixed Capacity Constraints."
            assert len(mobile_capac_constraints) == self.num_mobile_tx, "Issue with Mobile Capacity Constraints."
        #random, include capacity
        elif (self.random_test_flag and not self.exclude_cap_flag):
            print("Random Input, Including Capacity")
            #high is self.num_bands
            max_cap_val = self.num_bands
            if (self.overwrite_cap_flag):
                max_cap_val = self.overwrite_cap_val +1
            fixed_capac_constraints = np.random.randint(low = 1, high=max_cap_val, size=self.num_fixed_tx, dtype=int)
            mobile_capac_constraints = np.random.randint(low = 1, high=max_cap_val, size=self.num_mobile_tx, dtype=int)
        #do not include capacity
        else:
            print("Excluding Capacity")
            fixed_capac_constraints = self.num_bands*np.ones(self.num_fixed_tx, dtype=int)
            mobile_capac_constraints = self.num_bands*np.ones(self.num_mobile_tx, dtype=int)

        self.fixed_capac_constraints = fixed_capac_constraints
        self.mobile_capac_constraints = mobile_capac_constraints
        
        self.fixed_x_costs = fixed_x_costs
        self.mobile_x_costs = mobile_x_costs
        # self.connectivity_fixed_tx = connectivity_fixed_tx
        # self.connectivity_mobile_tx = connectivity_mobile_tx
        
        self.r_costs = r_costs
        self.z_costs = z_costs
        
        self.intrinsic_fixed_x = np.zeros((self.num_fixed_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
        self.intrinsic_mobile_x = np.zeros((self.num_mobile_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
        self.intrinsic_mobile_r = np.zeros((self.num_mobile_tx, self.num_mobile_tx_locs), dtype=np_impa_lib)

        #if capacity constraints are not excluded
        if (not self.exclude_cap_flag):
            self.fixed_capac_const_to_fixed_x_eq_const_m = np.zeros((self.num_fixed_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
            self.mobile_capac_const_to_mobile_x_eq_const_m = np.zeros((self.num_mobile_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
            self.fixed_x_eq_const_to_fixed_capac_const_m = deepcopy(fixed_x_costs)
            self.mobile_x_eq_const_to_mobile_capac_const_m = deepcopy(mobile_x_costs)

        print("num_fixed_tx: ", self.num_fixed_tx)
        print("num_mobile_tx: ", self.num_mobile_tx)
        print("num_bands: ", self.num_bands)
        print("num_time_steps: ", self.num_time_steps)
        print("num_rx_locs: ", self.num_rx_locs)
        print("num_mobile_tx_locs: ", self.num_mobile_tx_locs)
        print("filtering_flag: ", self.filtering_flag)
        if (self.filtering_flag):
            print("alpha: ", self.formatted_alpha)
        print("exclude_cap_flag: ", self.exclude_cap_flag)
        print("get_sol_approach: ", self.get_sol_approach)
        print("self.criteria_im: ", self.criteria_im)
        print("self.percentage_neg_im: ", self.percentage_neg_im)
        print("self.overwrite_im: ", self.overwrite_im)
        print("self.snr_threshold: ", self.snr_threshold)
        
        #construct equality constraint object
        self.model_eq_constraint = EqualityConstraintMOBARP(
            self.num_fixed_tx,
            self.num_mobile_tx,
            self.num_bands,
            self.num_time_steps,
            self.fixed_x_costs,
            self.mobile_x_costs,
            self.num_mobile_tx_locs, 
            self.z_costs,
            self.num_rx_locs,
            self.connectivity_mobile_tx,
            self.connectivity_fixed_tx,
            self.r_costs,
            self.exclude_cap_flag
        )
        
        #construct inequality constraint object
        self.model_ineq_constraint = InequalityConstraintMOBARP(
            self.num_fixed_tx,
            self.num_mobile_tx,
            self.num_bands,
            self.num_time_steps,
            self.fixed_capac_constraints,
            self.mobile_capac_constraints,
            self.alpha,
            self.filtering_flag,
            self.num_mobile_tx_locs,
            self.num_rx_locs,
            self.connectivity_mobile_tx, 
            self.connectivity_fixed_tx
        )
        
        #construct outputs object
        self.outputs = OutputsMOBARP(self.num_fixed_tx, self.num_mobile_tx, self.num_bands, self.num_time_steps, self.fixed_x_costs, self.mobile_x_costs, self.r_costs, 
                                    self.num_rx_locs, self.num_mobile_tx_locs, self.connectivity_mobile_tx, self.connectivity_fixed_tx, self.exclude_cap_flag)
        
        #construct auxiliary constraint object
        self.model_auxiliary_constraint = AuxiliaryConstraintMOBARP(
            self.num_mobile_tx, 
            self.num_bands, 
            self.num_time_steps, 
            self.num_mobile_tx_locs,
            self.alpha, 
            self.filtering_flag,
            self.connectivity_mobile_tx,
        )

        #initialize some messages to be used during iteration 0
        self.model_auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_m = np.zeros((self.mobile_x_costs.size, self.num_mobile_tx_locs), dtype= np_impa_lib)
        self.model_eq_constraint.r_eq_const_to_auxiliary_const_m = np.concatenate([self.r_costs.flatten()[..., np.newaxis]]*self.num_bands*self.num_time_steps, axis=1)
        self.model_eq_constraint.z_eq_const_to_auxiliary_const_m = np.zeros((self.mobile_x_costs.size, self.num_mobile_tx_locs), dtype= np_impa_lib)
        self.model_ineq_constraint.mobile_loc_eq_const_to_r_eq_const_m= np.zeros((self.num_mobile_tx, self.num_mobile_tx_locs), dtype=np_impa_lib)
        self.reshaped_connectivity_mobile_tx_aux = np.vstack([np.clip(np.sum(self.connectivity_mobile_tx, axis=3).transpose(1,2,0).reshape(-1, self.num_mobile_tx_locs), 0, 1)] * self.num_mobile_tx)
        reshaped_connectivity_mobile_tx_aux = self.reshaped_connectivity_mobile_tx_aux.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)
        self.model_eq_constraint.r_eq_const_to_auxiliary_const_m *=reshaped_connectivity_mobile_tx_aux
        self.model_ineq_constraint.set_cover_ineq_const_to_z_eq_const_m = np.zeros((self.num_time_steps, self.num_rx_locs, self.num_mobile_tx_locs, self.num_bands*self.num_mobile_tx), dtype=np_impa_lib)
        self.model_ineq_constraint.set_cover_ineq_const_to_fixed_x_eq_const_m = np.zeros((self.num_time_steps, self.num_rx_locs, self.num_fixed_tx*self.num_bands), dtype=np_impa_lib)
        self.conx_mob_tx_per_num_mob_tx_locs = np.vstack([np.clip(np.sum(self.connectivity_mobile_tx, axis=3).transpose(1,2,0).reshape(-1, self.num_mobile_tx_locs), 0, 1)] * self.num_mobile_tx)

        #if capacity constraints, initialize different messages for iteration 0
        if (self.exclude_cap_flag):
            self.model_ineq_constraint.fixed_capac_const_to_fixed_x_eq_const_m = np.zeros((self.num_fixed_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
            self.model_ineq_constraint.mobile_capac_const_to_mobile_x_eq_const_m = np.zeros((self.num_mobile_tx, self.num_bands, self.num_time_steps), dtype=np_impa_lib)
            mobile_x_eq_const_to_auxiliary_const_m = np.tile(self.mobile_x_costs.flatten()[:, np.newaxis], (1, self.num_mobile_tx_locs))
            mobile_x_eq_const_to_auxiliary_const_m*=self.conx_mob_tx_per_num_mob_tx_locs
            self.conx_fixed_tx_per_num_rx_locs =  self.connectivity_fixed_tx.reshape(-1, self.num_rx_locs)
            fixed_x_eq_const_to_set_cover_const_m = np.tile(self.fixed_x_costs.flatten()[:, np.newaxis], (1, self.num_rx_locs))
            fixed_x_eq_const_to_set_cover_const_m*=self.conx_fixed_tx_per_num_rx_locs

    def costs_connectivity_generation(self):
        #generate IM and connectivities
        normal_mean = 0
        normal_variance = 3
        self.normal_mean = normal_mean
        self.normal_variance = normal_variance
        
        num_fixed_tx = self.num_fixed_tx
        num_bands = self.num_bands
        num_mobile_tx_locs = self.num_mobile_tx_locs
        num_time_steps = self.num_time_steps
        num_mobile_tx = self.num_mobile_tx
        num_rx_locs = self.num_rx_locs
        
        if (not self.random_test_flag and not self.overwrite_im):
            fixed_x_costs = self.input_load[12]
            mobile_x_costs = self.input_load[13]
            r_costs = self.input_load[14]
            # print(f"np.min(fixed_x_costs): {np.min(fixed_x_costs)}, np.max(fixed_x_costs): {np.max(fixed_x_costs)}")
            # print(f"np.min(mobile_x_costs): {np.min(mobile_x_costs)}, np.max(mobile_x_costs): {np.max(mobile_x_costs)}")
        elif (self.criteria_im == 1): #normal
            fixed_x_costs = np.random.normal(normal_mean, normal_variance, size=(num_fixed_tx, num_bands, num_time_steps))
            mobile_x_costs = np.random.normal(normal_mean, normal_variance, size=(num_mobile_tx, num_bands, num_time_steps))
            r_costs = np.random.normal(normal_mean, normal_variance, size=(num_mobile_tx, num_mobile_tx_locs))
        elif (self.criteria_im == 2): #pos X, normal R
            fixed_x_costs = np.random.uniform(10, 100, size=(num_fixed_tx, num_bands, num_time_steps))
            mobile_x_costs = np.random.uniform(10, 100, size=(num_mobile_tx, num_bands, num_time_steps))
            r_costs = np.random.normal(normal_mean, normal_variance, size=(num_mobile_tx, num_mobile_tx_locs))
        
        percentage_neg = self.percentage_neg_im
        size_fixed = fixed_x_costs.size
        num_fixed_negatives = int(size_fixed * ((percentage_neg) / 100))
        num_fixed_positives = size_fixed - num_fixed_negatives
        fixed_array_sign = np.array([1] * num_fixed_positives + [-1] * num_fixed_negatives)
        np.random.shuffle(fixed_array_sign)
        fixed_x_costs = fixed_x_costs*fixed_array_sign.reshape((num_fixed_tx, num_bands, num_time_steps))
        
        size_mobile = mobile_x_costs.size
        num_mobile_negatives = int(size_mobile * ((percentage_neg) / 100))
        num_mobile_positives = size_mobile - num_mobile_negatives
        mobile_array_sign = np.array([1] * num_mobile_positives + [-1] * num_mobile_negatives)
        np.random.shuffle(mobile_array_sign)
        mobile_x_costs = mobile_x_costs*mobile_array_sign.reshape((num_mobile_tx, num_bands, num_time_steps))
        
        # size_r = r_costs.size
        # num_r_negatives = int(size_r * ((percentage) / 100))
        # num_r_positives = size_r - num_r_negatives
        # r_array_sign = np.array([1] * num_r_positives + [-1] * num_r_negatives)
        # np.random.shuffle(r_array_sign)
        # r_costs = r_costs*r_array_sign.reshape((num_mobile_tx, num_mobile_tx_locs))

        z_costs = np.zeros((mobile_x_costs.size, num_mobile_tx_locs), dtype = np_impa_lib)

        if (not self.random_test_flag):
            #set to extreme if reading from a file
            connectivity_fixed_tx = -10000*np.ones((num_fixed_tx, num_bands, num_time_steps, num_rx_locs), dtype=int)
            connectivity_mobile_tx = -10000*np.ones((num_mobile_tx_locs, num_bands, num_time_steps, num_rx_locs), dtype=int)
        else:
            connectivity_fixed_tx = np.random.randint(2, size=(num_fixed_tx, num_bands, num_time_steps, num_rx_locs), dtype=int)
            connectivity_mobile_tx = np.random.randint(2, size=(num_mobile_tx_locs, num_bands, num_time_steps, num_rx_locs), dtype=int)

        return fixed_x_costs, mobile_x_costs, connectivity_fixed_tx, connectivity_mobile_tx, r_costs, z_costs
        
        
    def run_impa(self):
        
        prev_extrinsic_fixed_x = None
        counter_wait = 0
        hard_decision_data = []

        for iter in range(0, self.num_iterations,):
            
            #only run when capacity constraints are included
            if (not self.exclude_cap_flag):
                #update messages from capacity constraint to x equality constraint
                self.model_ineq_constraint.ineq_capac_const_update(self.fixed_x_eq_const_to_fixed_capac_const_m, self.mobile_x_eq_const_to_mobile_capac_const_m)
                #perform filtering on messages from capacity constraint to x equality constraint
                self.model_ineq_constraint.process_filtering_capac_const(iter)
            #update messages from x equality constraint to auxiliary and set cover constraints
            self.model_eq_constraint.x_eq_const_to_auxiliary_and_set_cover_const_update(self.model_ineq_constraint.fixed_capac_const_to_fixed_x_eq_const_m, self.model_ineq_constraint.mobile_capac_const_to_mobile_x_eq_const_m, self.model_auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_m, self.model_ineq_constraint.set_cover_ineq_const_to_fixed_x_eq_const_m)
            #update messages from auxiliary constraint to z equality constraint
            self.model_auxiliary_constraint.auxiliary_const_to_z_eq_const_update(iter, self.model_eq_constraint.mobile_x_eq_const_to_auxiliary_const_m, self.model_eq_constraint.r_eq_const_to_auxiliary_const_m)
            #update messages from auxiliary constraint to r equality constraint
            self.model_auxiliary_constraint.auxiliary_const_to_r_eq_const_update(self.model_eq_constraint.z_eq_const_to_auxiliary_const_m, self.model_eq_constraint.mobile_x_eq_const_to_auxiliary_const_m)
            #update messages from z equlaity constraint to set cover inequality constraint
            self.model_eq_constraint.z_eq_const_to_set_cover_ineq_const_update(self.model_auxiliary_constraint.auxiliary_const_to_z_eq_const_m, self.model_ineq_constraint.set_cover_ineq_const_to_z_eq_const_m)
            #activate r equality constraint
            self.model_eq_constraint.r_eq_const_activation(self.model_auxiliary_constraint.auxiliary_const_to_r_eq_const_m, self.model_ineq_constraint.mobile_loc_eq_const_to_r_eq_const_m)
            #update messages from mobile loc equality constraint to r equality constraint
            self.model_ineq_constraint.mobile_loc_eq_const_to_r_eq_const_update(self.model_eq_constraint.r_eq_const_to_mobile_loc_eq_const_m)
            #perform filtering on messages from mobile loc equality constraint to r equality constraint
            self.model_ineq_constraint.process_filtering_mobile_loc_eq(iter)
            #update messages from set cover to z and fixed x equality constraints
            self.model_ineq_constraint.set_cover_ineq_const_update(self.model_eq_constraint.z_eq_const_to_set_cover_ineq_const_m, self.model_eq_constraint.fixed_x_eq_const_to_set_cover_const_m)
            #perform filtering on messages from set cover to z and fixed x equality constraints
            self.model_ineq_constraint.process_filtering_set_cover_const(iter)
            #update messages from z equality constraint to auxiliary constraint
            self.model_eq_constraint.z_eq_const_to_auxiliary_const_update(self.model_ineq_constraint.set_cover_ineq_const_to_z_eq_const_m)
            #update messages from auxiliary constraint to mobile x equality constraint
            self.model_auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_update(self.model_eq_constraint.z_eq_const_to_auxiliary_const_m, self.model_eq_constraint.r_eq_const_to_auxiliary_const_m)
            
            #only do this update when capacity constraints are included
            if (not self.exclude_cap_flag):
                #calculate messages from x equality constraint to capacity constraint
                self.mobile_x_eq_const_to_mobile_capac_const_m, self.fixed_x_eq_const_to_fixed_capac_const_m = self.model_eq_constraint.x_eq_const_activation(self.model_auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_m, self.model_ineq_constraint.set_cover_ineq_const_to_fixed_x_eq_const_m)
        
            #calculate extrinsic messages
            extrinsic_fixed_x, extrinsic_mobile_x, extrinsic_r, extrinsic_z = self.outputs.extrinsic_update(self.model_ineq_constraint.fixed_capac_const_to_fixed_x_eq_const_m, self.model_ineq_constraint.mobile_capac_const_to_mobile_x_eq_const_m, 
                                                                                                self.model_auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_m, self.model_ineq_constraint.set_cover_ineq_const_to_fixed_x_eq_const_m, 
                                                                                                self.model_auxiliary_constraint.auxiliary_const_to_r_eq_const_m, self.model_ineq_constraint.mobile_loc_eq_const_to_r_eq_const_m, 
                                                                                     self.model_auxiliary_constraint.auxiliary_const_to_z_eq_const_m, self.model_ineq_constraint.set_cover_ineq_const_to_z_eq_const_m)
            
            if prev_extrinsic_fixed_x is not None:
                changes_in_iteration = []
                
                if not np.allclose(prev_extrinsic_fixed_x, extrinsic_fixed_x, atol=1e-3):
                    change_size = np.linalg.norm(extrinsic_fixed_x - prev_extrinsic_fixed_x)
                    changes_in_iteration.append(f"extrinsic_fixed_x (change size: {change_size:.4f})")
                    flag_changes_fixed_tx = True
                else:
                    flag_changes_fixed_tx = False
                    
                if not np.allclose(prev_extrinsic_mobile_x, extrinsic_mobile_x, atol=1e-3):
                    change_size = np.linalg.norm(extrinsic_mobile_x - prev_extrinsic_mobile_x)
                    changes_in_iteration.append(f"extrinsic_mobile_x (change size: {change_size:.4f})")
                    flag_changes_mobile_tx = True
                else:
                    flag_changes_mobile_tx = False
                    
                if not np.allclose(prev_extrinsic_r, extrinsic_r, atol=1e-3):
                    change_size = np.linalg.norm(extrinsic_r - prev_extrinsic_r)
                    changes_in_iteration.append(f"extrinsic_r (change size: {change_size:.4f})")
                    flag_changes_r = True
                else:
                    flag_changes_r = False
                
                if not np.allclose(prev_extrinsic_z, extrinsic_z, atol=1e-3):
                    change_size = np.linalg.norm(extrinsic_z - prev_extrinsic_z)
                    changes_in_iteration.append(f"extrinsic_z (change size: {change_size:.4f})")
                    flag_changes_z = True
                else:
                    flag_changes_z = False
                
                if (flag_changes_fixed_tx or flag_changes_mobile_tx or flag_changes_r or flag_changes_z):
                    counter_wait=0
                else:
                    counter_wait+=1
            
            prev_extrinsic_fixed_x, prev_extrinsic_mobile_x, prev_extrinsic_r, prev_extrinsic_z = extrinsic_fixed_x, extrinsic_mobile_x, extrinsic_r, extrinsic_z
        
            self.end_time = time.time()
            self.impa_runtime = self.end_time - self.start_time
            # print(f"IMPA run_time: {self.impa_runtime}")
            
            self.extrinsic_fixed_x = extrinsic_fixed_x
            self.extrinsic_mobile_x = extrinsic_mobile_x
            self.extrinsic_r = extrinsic_r
            self.extrinsic_z = extrinsic_z
            
            #calculate intrinsic messages
            intrinsic_fixed_x = extrinsic_fixed_x + self.fixed_x_costs.flatten()
            intrinsic_mobile_x = extrinsic_mobile_x + self.mobile_x_costs.flatten()   
            intrinsic_r = extrinsic_r + self.r_costs.flatten()
            intrinsic_z = extrinsic_z + self.z_costs  
            
            self.intrinsic_fixed_x = intrinsic_fixed_x
            self.intrinsic_mobile_x = intrinsic_mobile_x
            self.intrinsic_r = intrinsic_r
            self.intrinsic_z = intrinsic_z
            
            #get hard decision on variables based on intrinsic messages
            self.hard_decision_analysis()
            
            data_iter = {
                    "iter": iter,
                    "hard_decision_fixed_x": self.hard_decision_fixed_x,
                    "hard_decision_mobile_x": self.hard_decision_mobile_x,
                    "hard_decision_r": self.hard_decision_r,
                    "hard_decision_z": self.hard_decision_z
                }
            
            hard_decision_data.append(data_iter)
            
            self.stopping_at_iter = iter #last iteration that is done
            if (counter_wait==5):
                self.stopping_at_iter = iter-counter_wait+1
                break
                
        self.hard_decision_data = hard_decision_data
        
        data_iter_fixed_tx_hd = [data["hard_decision_fixed_x"].tolist() for data in hard_decision_data]
        stable_fixed_tx_hd_iter = self.find_stable_hard_decision(data_iter_fixed_tx_hd)
        
        data_iter_mobile_tx_hd = [data["hard_decision_mobile_x"].tolist() for data in hard_decision_data]
        stable_mobile_tx_hd_iter = self.find_stable_hard_decision(data_iter_mobile_tx_hd)

        data_iter_r_hd = [data["hard_decision_r"].tolist() for data in hard_decision_data]
        stable_r_hd_iter = self.find_stable_hard_decision(data_iter_r_hd)
        
        data_iter_z_hd = [data["hard_decision_z"].tolist() for data in hard_decision_data]
        stable_z_hd_iter = self.find_stable_hard_decision(data_iter_z_hd)
        
        print(f"Stopping at iteration: {self.stopping_at_iter}")
        self.stable_hd_iter = {"hard_decision_fixed_x": stable_fixed_tx_hd_iter, "hard_decision_mobile_x": stable_mobile_tx_hd_iter, "hard_decision_r": stable_r_hd_iter, "hard_decision_z": stable_z_hd_iter}
        print(f"Stable HD Iter: {self.stable_hd_iter}")

    def find_stable_hard_decision(self,data_history):
        if (data_history[0] == [-100]*len(data_history[0])):
            return -100
        stable_value = data_history[0]
        stable_hard_decision_iter = 0
        for i in range(1, len(data_history)):
            if data_history[i] != stable_value:
                stable_value = data_history[i]
                stable_hard_decision_iter = i
            elif all(val == stable_value for val in data_history[i:]):
                return stable_hard_decision_iter
        return self.num_iterations
        # return None
    
    def hard_decision_analysis(self,):
        
        intrinsic_fixed_x = self.intrinsic_fixed_x
        intrinsic_mobile_x = self.intrinsic_mobile_x
        intrinsic_r = self.intrinsic_r
        intrinsic_z = self.intrinsic_z

        hard_decision_fixed_x = np.array(deepcopy(intrinsic_fixed_x))
        hard_decision_fixed_x[intrinsic_fixed_x > self.threshold] = 0
        hard_decision_fixed_x[intrinsic_fixed_x <= self.threshold] = 1
        
        conx_fixed_tx_per_num_rx_locs =  self.connectivity_fixed_tx.reshape(-1, self.num_rx_locs)
        excluded_fixed_tx = np.where(np.clip(np.sum(conx_fixed_tx_per_num_rx_locs, axis=1), 0,1) !=1)[0]
        #do not consider variables that have zero connectivity
        hard_decision_fixed_x[excluded_fixed_tx] = -100

        # if (len(excluded_fixed_tx) !=0):
        #     print("Some Fixed TX have zero-connectivity.")

        hard_decision_mobile_x = np.array(deepcopy(intrinsic_mobile_x))
        hard_decision_mobile_x[intrinsic_mobile_x > self.threshold] = 0
        hard_decision_mobile_x[intrinsic_mobile_x <= self.threshold] = 1

        excluded_mobile_tx = np.where(np.clip(np.sum(self.conx_mob_tx_per_num_mob_tx_locs, axis=1), 0,1) !=1)[0]
        #do not consider variables that have zero connectivity
        hard_decision_mobile_x[excluded_mobile_tx] = -100
        # if (len(excluded_mobile_tx) !=0):
        #     print("Some Mobile TX have zero-connectivity.")

        conx_mob_tx_r = self.conx_mob_tx_per_num_mob_tx_locs.reshape(self.num_mobile_tx, -1, self.num_mobile_tx_locs).transpose(0, 2, 1).reshape(-1, self.num_bands*self.num_time_steps)
        excluded_r = np.where(np.sum(conx_mob_tx_r, axis=1) ==0)[0]

        # if (len(excluded_r) !=0):
        #     print("Some r-eq. const. have zero-connectivity.")
            #not supposed to have this issue, if that is not the case, it is the issue with connectivity

        hard_decision_r = np.array(deepcopy(intrinsic_r))
        hard_decision_r[intrinsic_r > self.threshold] = 0
        hard_decision_r[intrinsic_r <= self.threshold] = 1

        #do not consider variables that have zero connectivity
        hard_decision_r[excluded_r] = -100
        
        hard_decision_z = np.array(deepcopy(intrinsic_z))
        hard_decision_z[intrinsic_z > self.threshold] = 0
        hard_decision_z[intrinsic_z <= self.threshold] = 1

        self.hard_decision_fixed_x = hard_decision_fixed_x
        self.hard_decision_mobile_x = hard_decision_mobile_x
        
        self.hard_decision_r = hard_decision_r
        self.hard_decision_z = hard_decision_z
        
        # # print(self.mobile_x_costs[0,0])

    def run_analysis(self,):
        
        # check fixed capacity constraints
        self.reshaped_hard_decision_fixed_x = self.hard_decision_fixed_x.reshape(self.num_fixed_tx, self.num_bands, self.num_time_steps)
        mask_fixed = self.reshaped_hard_decision_fixed_x !=-100
        used_fixed_tx_capacities = np.sum(self.reshaped_hard_decision_fixed_x , axis=1, where=mask_fixed ==1).astype(np.int64)
        reshaped_fixed_tx_capacities = np.concatenate([self.fixed_capac_constraints[..., np.newaxis]]*self.num_time_steps, axis=1)
        flag_fixed_tx_capacities = np.all(used_fixed_tx_capacities<=reshaped_fixed_tx_capacities)
        indices_violated_fixed_tx_capacities = np.where(used_fixed_tx_capacities>reshaped_fixed_tx_capacities)
        
        print('-----------')
        print("Fixed TX Capacity Flag: ", flag_fixed_tx_capacities)

        self.flag_fixed_tx_capacities = flag_fixed_tx_capacities
        self.used_fixed_tx_capacities = used_fixed_tx_capacities
        self.indices_violated_fixed_tx_capacities = indices_violated_fixed_tx_capacities

        print("Used Capacity\t\t\t\t       Max Capacity")
        for row1, row2 in zip(used_fixed_tx_capacities, self.fixed_capac_constraints[..., np.newaxis]):
            print("{:<50}{}".format(str(row1),str(row2)))

        #check mobile capacity constraints
        self.reshaped_hard_decision_mobile_x = self.hard_decision_mobile_x.reshape(self.num_mobile_tx, self.num_bands, self.num_time_steps)
        mask_mobile = self.reshaped_hard_decision_mobile_x !=-100
        used_mobile_tx_capacities = np.sum(self.reshaped_hard_decision_mobile_x, axis=1, where = mask_mobile==1).astype(np.int64)
        reshaped_mobile_tx_capacities = np.concatenate([self.mobile_capac_constraints[..., np.newaxis]]*self.num_time_steps, axis=1)
        flag_mobile_tx_capacities = np.all(used_mobile_tx_capacities<=reshaped_mobile_tx_capacities)
        indices_violated_mobile_tx_capacities = np.where(used_mobile_tx_capacities>reshaped_mobile_tx_capacities)

        print('-----------')
        print("Mobile TX Capacity Flag: ", flag_mobile_tx_capacities)

        self.flag_mobile_tx_capacities = flag_mobile_tx_capacities
        self.used_mobile_tx_capacities = used_mobile_tx_capacities
        self.indices_violated_mobile_tx_capacities = indices_violated_mobile_tx_capacities

        print("Used Capacity\t\t\t\t       Max Capacity")
        for row1, row2 in zip(used_mobile_tx_capacities, self.mobile_capac_constraints[..., np.newaxis]):
            print("{:<50}{}".format(str(row1), str(row2)))

        #check mobile-loc assignment constraints
        self.reshaped_hard_decision_r = self.hard_decision_r.reshape(self.num_mobile_tx, self.num_mobile_tx_locs)
        mask_r = self.reshaped_hard_decision_r !=-100
        sum_used_mobile_loc = np.sum(self.reshaped_hard_decision_r, axis=1, where=mask_r==1).astype(np.int64)
        flag_sum_used_mobile_loc = np.all(sum_used_mobile_loc==1)
        indices_violated_tx_loc_assignment = np.where(sum_used_mobile_loc!=1)

        print('-----------')
        print("Mobile TX-Location Assignment")
        print('-----------')
        print("Mobile TX-Location Assignment Flag: ", flag_sum_used_mobile_loc)
        
        self.flag_sum_used_mobile_loc = flag_sum_used_mobile_loc
        self.sum_used_mobile_loc = sum_used_mobile_loc
        self.indices_violated_tx_loc_assignment = indices_violated_tx_loc_assignment

        self.reshaped_hard_decision_z = self.hard_decision_z.reshape(self.num_mobile_tx, self.num_bands, self.num_time_steps, self.num_mobile_tx_locs).transpose(0,3,1,2)
        
        indices_activated_r = np.where(self.reshaped_hard_decision_r == 1)

        configurations_activated_r = list(zip(*indices_activated_r))

        results_r = defaultdict(list)
        for (i_prime, n) in configurations_activated_r:
            results_r[i_prime].append(n)

        self.configurations_activated_r = configurations_activated_r
        
        for i_prime in range(self.num_mobile_tx):
            print(f"Mobile TX {i_prime} is at location {results_r[i_prime]}")

        missing_mobile_tx_to_loc_assignment = [key for key in range(self.num_mobile_tx) if not results_r[key]]

        # if missing_mobile_tx_to_loc_assignment:
        #     print(f"Missing Assignments for Mobile TX: {missing_mobile_tx_to_loc_assignment}")
        
        print('-----------')
        print("Set Cover Constraints")
        print('-----------')
        
        #only select activated z that are satisfied to get selected mobile TX
        indices_activated_z = np.where(self.reshaped_hard_decision_z==1)
        configurations_activated_z = list(zip(*indices_activated_z))
        print(f"{len(configurations_activated_z)} activated auxiliary constraints")
        set_z_i_prime_n = [(config[0], config[1]) for config in configurations_activated_z]
        matching_z_configurations = [configurations_activated_z[index] for index, config in enumerate(set_z_i_prime_n) if config in configurations_activated_r]
        set_z_i_prime_j_k = [(config[0], config[2], config[3]) for config in matching_z_configurations]
        indices_activated_mobile_x = np.where(self.reshaped_hard_decision_mobile_x==1)
        configurations_activated_mobile_x = list(zip(*indices_activated_mobile_x))
        matching_z_configurations = [matching_z_configurations[index] for index, config in enumerate(set_z_i_prime_j_k) if config in configurations_activated_mobile_x]
        print(f"{len(matching_z_configurations)} valid activated auxiliary constraints")
        results_mobile_tx = defaultdict(list)
        for matching_z_config in matching_z_configurations:
            assigned_rxs = np.where(self.connectivity_mobile_tx[matching_z_config[1], matching_z_config[2], matching_z_config[3]])
            candidate_config = (matching_z_config[0], matching_z_config[2], matching_z_config[3])
            for index in assigned_rxs[0]:
                if (candidate_config not in results_mobile_tx[index]):
                    results_mobile_tx[index].append(candidate_config)
        results_mobile_tx = {key: results_mobile_tx[key] for key in sorted(results_mobile_tx)}

        #select fixed TX
        indices_activated_fixed_tx = np.where(self.reshaped_hard_decision_fixed_x==1)
        configurations_activated_fixed_tx = list(zip(*indices_activated_fixed_tx))
        results_fixed_tx = defaultdict(list)
        for fixed_config in configurations_activated_fixed_tx:
            assigned_rxs = np.where(self.connectivity_fixed_tx[fixed_config[0], fixed_config[1], fixed_config[2]])
            for index in assigned_rxs[0]:
                if (fixed_config not in results_fixed_tx[index]):
                    results_fixed_tx[index].append(fixed_config)
        results_fixed_tx = {key: results_fixed_tx[key] for key in sorted(results_fixed_tx)}

        # print('----')
        #get approach for finding solution, efficient one is 2 (default)
        approach = self.get_sol_approach
        counts_k_l = defaultdict(list)

        if (approach==1):
            activated_x = defaultdict(lambda: defaultdict(list))
            self.get_activations_x_approach_1(results_mobile_tx, activated_x, counts_k_l, is_mobile=True)
            self.get_activations_x_approach_1(results_fixed_tx, activated_x, counts_k_l, is_mobile=False)
        elif (approach==2):
            activated_x = defaultdict(list)
            self.get_activations_x_approach_2(results_mobile_tx, activated_x, counts_k_l, is_mobile=True)
            self.get_activations_x_approach_2(results_fixed_tx, activated_x, counts_k_l, is_mobile=False)
        
        #investigate Set Cover Constraints
        self.set_cover_investigation(counts_k_l)
        self.activated_x = activated_x
        self.counts_k_l = counts_k_l
        self.results_mobile_tx = results_mobile_tx
        self.results_fixed_tx = results_fixed_tx

        
        #get solution that minimizes the objective function
        if (approach==1):
            x_assignment, used_x_list = self.get_solution_approach_1(results_mobile_tx, results_fixed_tx)
        elif (approach==2):
            x_assignment, used_x_list = self.get_solution_approach_2(results_mobile_tx, results_fixed_tx)
            
        #exit()
        self.x_assignment = x_assignment
        self.used_x_list = used_x_list

        self.objective_cost = len(used_x_list)
        #for key, value in x_assignment.items():
        #    print(f"(k,l) = {key}: {value}")
        
        count_m = sum(1 for item in used_x_list if item[0] == 'm')
        count_f = sum(1 for item in used_x_list if item[0] == 'f')
        self.count_m = count_m
        self.count_f = count_f

        print('-----------')
        print(f"Objective Value: {len(used_x_list)}/{self.fixed_x_costs.size + self.mobile_x_costs.size}")
        print(f"Number of Mobile TX: {count_m}")
        print(f"Number of Fixed TX: {count_f}")
        print('-----------')
        if (self.flag_fixed_tx_capacities and self.flag_mobile_tx_capacities and self.flag_sum_used_mobile_loc and self.flag_k_l):
            print("All constraints are satisfied.")
        else:
            not_satisfied = []
            if not self.flag_fixed_tx_capacities:
                not_satisfied.append("Fixed Capacities")
            if not self.flag_mobile_tx_capacities:
                not_satisfied.append("Mobile Capacities")
            if not self.flag_sum_used_mobile_loc:
                not_satisfied.append("Mobile TX-LOC Assignment")
            if not self.flag_k_l:
                not_satisfied.append("Set Cover Constraints")
            print("Constraints not satisfied:", ", ".join(not_satisfied))

        # print('-----------')
        activated_x_before_pp = len(activated_x)
        activated_x_after_pp = len(used_x_list)
        print(f"Number of activated TX before post-processing: {activated_x_before_pp}/{self.fixed_x_costs.size + self.mobile_x_costs.size}")
        print(f"Number of selected TX after post-processing: {activated_x_after_pp}/{self.fixed_x_costs.size + self.mobile_x_costs.size}")
        print("self.used_x_list: \n", self.used_x_list)
        
        if (self.save_flag):
            self.results_composed = [
                    self.filtering_flag, #0
                    self.alpha, #1
                    self.impa_runtime, #2
                    self.objective_cost, #3
                    dict(self.activated_x), #4
                    dict(self.sorted_activated_x), #5
                    self.count_m, #6
                    self.count_f, #7
                    self.used_x_list, #8
                    self.configurations_activated_r, #9
                    self.counts_k_l, #10
                    self.violated_k_l, #11
                    self.flag_fixed_tx_capacities, #12
                    self.flag_mobile_tx_capacities, #13
                    self.flag_sum_used_mobile_loc, #14
                    self.flag_k_l, #15
                    self.used_fixed_tx_capacities, #16
                    self.used_mobile_tx_capacities, #17
                    self.sum_used_mobile_loc, #18
                    self.indices_violated_fixed_tx_capacities, #19
                    self.indices_violated_mobile_tx_capacities, #20
                    self.indices_violated_tx_loc_assignment, #21
                    self.fixed_x_costs, #22
                    self.mobile_x_costs, #23
                    self.r_costs, #24
                    self.stopping_at_iter, #25
                    self.hard_decision_data, #26
                    self.stable_hd_iter, #27
            ]

    def get_activations_x_approach_1(self, results, x, count_dict, is_mobile):
        #get activated x according to approach 1 (dictionary with keys being configurations, and values being lists for mobile and fixed coverage)
        for key, value_list in results.items():
            for value in value_list:
                index_tx, j, k = value
                if (k, key) in count_dict:
                    count_dict[(k, key)] += 1
                else:
                    count_dict[(k, key)] = 1
                if is_mobile:
                    if (key not in x[value]['m']):
                        x[value]['m'].append(key)
                else:
                    if (key not in x[value]['f']):
                        x[value]['f'].append(key)
    
    def get_activations_x_approach_2(self, results, x, count_dict, is_mobile):
        #get activated x according to approach 2 (dictionary of with keys being configuration, and values are list of covered RX). Unlike approach 1, fixed and mobile coverage
        #are decoupled
        for key, value_list in results.items():
            for value in value_list:
                index_tx, j, k = value
                if (k, key) in count_dict:
                    count_dict[(k, key)] += 1
                else:
                    count_dict[(k, key)] = 1
                
                if is_mobile:
                    config = ('m', value)
                    if (key not in x[config]):
                        x[config].append(key)
                else:
                    config = ('f', value)
                    if (key not in x[config]):
                        x[config].append(key)

    def set_cover_investigation(self, counts):
        #statistics on set cover constraints
        all_k_l_combinations = {(k, l) for k in range(self.num_time_steps) for l in range(self.num_rx_locs)}
        violated_k_l = all_k_l_combinations - counts.keys()

        self.violated_k_l = violated_k_l

        if (len(violated_k_l)):
            self.flag_k_l = False
        else:
            self.flag_k_l = True

    def get_solution_approach_1(self, results_mobile_tx, results_fixed_tx):
        #get solution by looping over the dictionary obtained from approach 1, and assigning configurations to RX
        activated_x = self.activated_x
        sorted_activated_x_items = sorted(activated_x.items(), key=self.sum_fixed_mobile_per_x, reverse=True)
        sorted_activated_x = defaultdict(lambda: defaultdict(list), sorted_activated_x_items)
        self.sorted_activated_x = sorted_activated_x

        x_assignment = defaultdict(list)
        used_x_list = []
        k_l_combinations = [(k, l) for k in range(self.num_time_steps) for l in range(self.num_rx_locs)]
        for config in sorted_activated_x:
            first, second = ('m', 'f') if len(sorted_activated_x[config]['m']) > len(sorted_activated_x[config]['f']) else ('f', 'm')
            for element in sorted_activated_x[config][first]:
                if ((config[-1], element) not in k_l_combinations):
                    continue
                used_x = (first, config)
                x_assignment[(config[-1], element)] = used_x
                k_l_combinations.remove((config[-1], element))
                if (used_x not in used_x_list):
                    used_x_list.append(used_x)
            
            for element in sorted_activated_x[config][second]:
                if ((config[-1], element) not in k_l_combinations):
                    continue
                used_x = (second, config)
                x_assignment[(config[-1], element)] = used_x
                k_l_combinations.remove((config[-1], element))
                if (used_x not in used_x_list):
                    used_x_list.append(used_x)
        if (k_l_combinations):
            print("Set Cover Constraints Flag: ", self.flag_k_l)
            print(f"{len(k_l_combinations)}/{self.num_time_steps*self.num_rx_locs} violated set cover constraints.")
        else:
            print("Set Cover Constraints Flag: ", self.flag_k_l)
            print(f"All {self.num_time_steps*self.num_rx_locs} Set Cover Constraints are satisfied.")

        return x_assignment, used_x_list

    def get_solution_approach_2(self, results_mobile_tx, results_fixed_tx):
        #get solution by looping over the dictionary obtained from approach 2, and assigning configurations to RX (more efficient)
        activated_x = self.activated_x
        sorted_activated_x_items = sorted(activated_x.items(), key=self.length_coverage, reverse=True)
        sorted_activated_x = defaultdict(lambda: defaultdict(list), sorted_activated_x_items)
        self.sorted_activated_x = sorted_activated_x

        x_assignment = defaultdict(list)
        used_x_list = []
        k_l_combinations = [(k, l) for k in range(self.num_time_steps) for l in range(self.num_rx_locs)]
        for config in sorted_activated_x:
            for element in sorted_activated_x[config]:
                if ((config[1][-1], element) not in k_l_combinations):
                    continue
                x_assignment[(config[1][-1], element)] = config
                k_l_combinations.remove((config[1][-1], element))
                if (config not in used_x_list):
                    used_x_list.append(config)
        if (k_l_combinations):
            print("Set Cover Constraints Flag: ", self.flag_k_l)
            print(f"{len(k_l_combinations)}/{self.num_time_steps*self.num_rx_locs} violated set cover constraints.")
        else:
            print("Set Cover Constraints Flag: ", self.flag_k_l)
            print(f"All {self.num_time_steps*self.num_rx_locs} Set Cover Constraints are satisfied.")

        return x_assignment, used_x_list

    def sum_fixed_mobile_per_x(self, item):
        #used to sort the dictionary obtained in approach 1
        m_len = len(item[1].get('m', []))
        f_len = len(item[1].get('f', []))
        return m_len + f_len

    def length_coverage(self, item):
        #used to sort the dictionary obtained in approach 2
        return len(item[1])

    def save_outputs(self):
        #save outputs if save_flag is True
        output_file = os.path.join(self.folder_outputs, "outputs_set"+str(self.test_file)+".pkl")
        with open(output_file, "wb") as f:
        #with open(f"{self.folder_outputs}/outputs_set{self.test_file}.pkl", "wb") as f:
            pkl.dump(self.results_composed, f)