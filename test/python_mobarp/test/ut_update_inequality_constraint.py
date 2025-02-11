# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
import update_inequality_constraint
from ut_utils import *


def ut_ineq_capac_const_update(ut_name, num_fixed_tx, num_mobile_tx, num_bands, num_time_steps, fixed_capac_constraints, mobile_capac_constraints, alpha, filtering_flag, num_mobile_tx_locs, num_rx_locs,
                                connectivity_mobile_tx, connectivity_fixed_tx ):
    
    NUM_FIXED_TX = num_fixed_tx
    NUM_MOBILE_TX = num_mobile_tx
    NUM_BANDS = num_bands
    NUM_TIME_STEPS = num_time_steps
    NUM_RX_LOCS = num_rx_locs
    NUM_MOBILE_TX_LOCS = num_mobile_tx_locs
    ALPHA = alpha
    FILTERING_FLAG = filtering_flag
    
    f_input1 = os.getcwd() + "/../ut_inputs/fixed_capacity_constraints.npy"
    fixed_capac_constraints = fixed_capac_constraints.astype(np.int32)
    np.save(f_input1, fixed_capac_constraints.flatten())
    
    f_input2 = os.getcwd() + "/../ut_inputs/mobile_capacity_constraints.npy"
    mobile_capac_constraints = mobile_capac_constraints.astype(np.int32)
    np.save(f_input2, mobile_capac_constraints.flatten())
    
    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)
    
    f_input3 = os.getcwd() + "/../ut_inputs/connectivity_mobile_tx.npy"
    connectivity_mobile_tx = connectivity_mobile_tx.astype(np.int32)
    np.save(f_input3, connectivity_mobile_tx.flatten())
    
    f_input4 = os.getcwd() + "/../ut_inputs/connectivity_fixed_tx.npy"
    connectivity_fixed_tx = connectivity_fixed_tx.astype(np.int32)
    np.save(f_input4, connectivity_fixed_tx.flatten())
    
    inequality_constraint = update_inequality_constraint.InequalityConstraintMOBARP(NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, 
                                                    fixed_capac_constraints, mobile_capac_constraints, 
                                                    ALPHA, FILTERING_FLAG, NUM_MOBILE_TX_LOCS, NUM_RX_LOCS, 
                                                    connectivity_mobile_tx, connectivity_fixed_tx)

    f_input5 = os.getcwd() + "/../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy"
    conx_fixed_tx_per_num_rx_locs = inequality_constraint.conx_fixed_tx_per_num_rx_locs.astype(np.int32)
    np.save(f_input5, conx_fixed_tx_per_num_rx_locs.flatten())
    
    f_input6 = os.getcwd() + "/../ut_inputs/conx_mob_tx_rx.npy"
    conx_mob_tx_rx = inequality_constraint.conx_mob_tx_rx.astype(np.int32)
    np.save(f_input6, conx_mob_tx_rx.flatten())
    
    normal_mean = 0
    normal_variance = 100
        
    if ut_name == "IneqCapacConstUpdate":
    
        f_input1 = os.getcwd() + "/../ut_inputs/fixed_x_eq_const_to_fixed_capac_const_m.npy"
        # fixed_x_eq_const_to_fixed_capac_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))
        fixed_x_eq_const_to_fixed_capac_const_m_pure = np.random.uniform(low=-100, high=-1, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))
        fixed_x_eq_const_to_fixed_capac_const_m_pure = fixed_x_eq_const_to_fixed_capac_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, fixed_x_eq_const_to_fixed_capac_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_mobile_capac_const_m.npy"
        # mobile_x_eq_const_to_mobile_capac_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))
        mobile_x_eq_const_to_mobile_capac_const_m_pure = np.random.uniform(low=-100, high=-1, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))
        mobile_x_eq_const_to_mobile_capac_const_m_pure = mobile_x_eq_const_to_mobile_capac_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, mobile_x_eq_const_to_mobile_capac_const_m_pure.flatten())
    
        inequality_constraint.ineq_capac_const_update(fixed_x_eq_const_to_fixed_capac_const_m_pure, mobile_x_eq_const_to_mobile_capac_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/fixed_capac_const_to_fixed_x_eq_const_m_dummy_pure"
        output_file_1 = open(f_output_1, "wb")
        fixed_capac_const_to_fixed_x_eq_const_m_dummy_pure = inequality_constraint.fixed_capac_const_to_fixed_x_eq_const_m_dummy.astype(np_impa_lib)
        np.save(output_file_1, fixed_capac_const_to_fixed_x_eq_const_m_dummy_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/mobile_capac_const_to_mobile_x_eq_const_m_dummy_pure"
        mobile_capac_const_to_mobile_x_eq_const_m_dummy_pure = inequality_constraint.mobile_capac_const_to_mobile_x_eq_const_m_dummy.astype(np_impa_lib)
        output_file_2 = open(f_output_2, "wb")
        np.save(output_file_2, mobile_capac_const_to_mobile_x_eq_const_m_dummy_pure.flatten(), allow_pickle=True)
        output_file_2.close()
        
    elif ut_name == "MobileLocEqConst2REqConstUpdate":
        
        f_input1 = os.getcwd() + "/../ut_inputs/r_eq_const_to_mobile_loc_eq_const_m.npy"
        r_eq_const_to_mobile_loc_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS))
        r_eq_const_to_mobile_loc_eq_const_m_pure = r_eq_const_to_mobile_loc_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, r_eq_const_to_mobile_loc_eq_const_m_pure.flatten())
        
        inequality_constraint.mobile_loc_eq_const_to_r_eq_const_update(r_eq_const_to_mobile_loc_eq_const_m_pure)

        f_output_1 = os.getcwd() + "/../ut_results/mobile_loc_eq_const_to_r_eq_const_m_dummy_pure"
        output_file_1 = open(f_output_1, "wb")
        mobile_loc_eq_const_to_r_eq_const_m_dummy_pure = inequality_constraint.mobile_loc_eq_const_to_r_eq_const_m_dummy.astype(np_impa_lib)
        np.save(output_file_1, mobile_loc_eq_const_to_r_eq_const_m_dummy_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
    elif ut_name == "SetCoverIneqConstUpdate":
        # set_cover_ineq_const_update(self, z_eq_const_to_set_cover_ineq_const_m, fixed_x_eq_const_to_set_cover_const_m)
        
        f_input1 = os.getcwd() + "/../ut_inputs/z_eq_const_to_set_cover_ineq_const_m.npy"
        z_eq_const_to_set_cover_ineq_const_m_pure = np.random.uniform(low=1, high=100, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS*NUM_MOBILE_TX_LOCS, NUM_RX_LOCS))
        z_eq_const_to_set_cover_ineq_const_m_pure = z_eq_const_to_set_cover_ineq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, z_eq_const_to_set_cover_ineq_const_m_pure.flatten())
        
        
        f_input2 = os.getcwd() + "/../ut_inputs/fixed_x_eq_const_to_set_cover_const_m.npy"
        fixed_x_eq_const_to_set_cover_const_m_pure = np.random.uniform(low=1, high=100, size=(NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_RX_LOCS))
        fixed_x_eq_const_to_set_cover_const_m_pure = fixed_x_eq_const_to_set_cover_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, fixed_x_eq_const_to_set_cover_const_m_pure.flatten())
        
        inequality_constraint.set_cover_ineq_const_update(z_eq_const_to_set_cover_ineq_const_m_pure, fixed_x_eq_const_to_set_cover_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/set_cover_ineq_const_to_fixed_x_eq_const_m_dummy_pure"
        output_file_1 = open(f_output_1, "wb")
        set_cover_ineq_const_to_fixed_x_eq_const_m_dummy_pure = inequality_constraint.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy.astype(np_impa_lib)
        np.save(output_file_1, set_cover_ineq_const_to_fixed_x_eq_const_m_dummy_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/set_cover_ineq_const_to_z_eq_const_m_dummy_pure"
        output_file_2 = open(f_output_2, "wb")
        set_cover_ineq_const_to_z_eq_const_m_dummy_pure = inequality_constraint.set_cover_ineq_const_to_z_eq_const_m_dummy.astype(np_impa_lib)
        np.save(output_file_2, set_cover_ineq_const_to_z_eq_const_m_dummy_pure.flatten(), allow_pickle=True)
        output_file_2.close()