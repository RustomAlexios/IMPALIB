# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
import update_auxiliary_constraint
from ut_utils import *


def ut_auxiliary_constraint(ut_name, num_mobile_tx, num_bands, num_time_steps, num_mobile_tx_locs, alpha, filtering_flag, connectivity_mobile_tx):
    
    NUM_MOBILE_TX = num_mobile_tx
    NUM_BANDS = num_bands
    NUM_TIME_STEPS = num_time_steps
    NUM_MOBILE_TX_LOCS = num_mobile_tx_locs
    ALPHA = alpha
    FILTERING_FLAG = filtering_flag
    
    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)
    
    f_input1 = os.getcwd() + "/../ut_inputs/connectivity_mobile_tx.npy"
    connectivity_mobile_tx = connectivity_mobile_tx.astype(np.int32)
    np.save(f_input1, connectivity_mobile_tx.flatten())
    
    auxiliary_constraint = update_auxiliary_constraint.AuxiliaryConstraintMOBARP(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS, ALPHA, FILTERING_FLAG, connectivity_mobile_tx)

    f_input2 = os.getcwd() + "/../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy"
    conx_mob_tx_per_num_mob_tx_locs = auxiliary_constraint.conx_mob_tx_per_num_mob_tx_locs.astype(np.int32)
    np.save(f_input2, conx_mob_tx_per_num_mob_tx_locs.flatten())
    
    normal_mean = 0
    normal_variance = 100
        
    # if ut_name == "AuxiliaryConst2ZEqConstUpdate":
    
    #     f_input1 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy"
    #     mobile_x_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
    #     mobile_x_eq_const_to_auxiliary_const_m_pure = mobile_x_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
    #     np.save(f_input1, mobile_x_eq_const_to_auxiliary_const_m_pure.flatten())
        
    #     f_input2 = os.getcwd() + "/../ut_inputs/r_eq_const_to_auxiliary_const_m.npy"
    #     r_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, NUM_BANDS*NUM_TIME_STEPS))
    #     r_eq_const_to_auxiliary_const_m_pure = r_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
    #     np.save(f_input2, r_eq_const_to_auxiliary_const_m_pure.flatten())
    
    #     auxiliary_constraint.auxiliary_const_to_z_eq_const_update(mobile_x_eq_const_to_auxiliary_const_m_pure, r_eq_const_to_auxiliary_const_m_pure)
        
    #     f_output_1 = os.getcwd() + "/../ut_results/auxiliary_const_to_z_eq_const_m_pure"
    #     output_file_1 = open(f_output_1, "wb")
    #     auxiliary_const_to_z_eq_const_m_pure = auxiliary_constraint.auxiliary_const_to_z_eq_const_m.astype(np_impa_lib)
    #     np.save(output_file_1, auxiliary_const_to_z_eq_const_m_pure.flatten(), allow_pickle=True)
    #     output_file_1.close()
        
    # elif ut_name == "AuxiliaryConst2REqConstUpdate":

    #     f_input1 = os.getcwd() + "/../ut_inputs/z_eq_const_to_auxiliary_const_m.npy"
    #     z_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
    #     z_eq_const_to_auxiliary_const_m_pure = z_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
    #     np.save(f_input1, z_eq_const_to_auxiliary_const_m_pure.flatten())
        
    #     f_input2 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy"
    #     mobile_x_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
    #     mobile_x_eq_const_to_auxiliary_const_m_pure = mobile_x_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
    #     np.save(f_input2, mobile_x_eq_const_to_auxiliary_const_m_pure.flatten())
        
    #     auxiliary_constraint.auxiliary_const_to_r_eq_const_update(z_eq_const_to_auxiliary_const_m_pure, mobile_x_eq_const_to_auxiliary_const_m_pure)
        
    #     f_output_1 = os.getcwd() + "/../ut_results/auxiliary_const_to_r_eq_const_m_pure"
    #     output_file_1 = open(f_output_1, "wb")
    #     auxiliary_const_to_r_eq_const_m_pure = auxiliary_constraint.auxiliary_const_to_r_eq_const_m.astype(np_impa_lib)
    #     np.save(output_file_1, auxiliary_const_to_r_eq_const_m_pure.flatten(), allow_pickle=True)
    #     output_file_1.close()
    
    if ut_name == "AuxiliaryConst2ZAndREqConstUpdate":
        
        f_input1 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_auxiliary_const_m.npy"
        mobile_x_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
        mobile_x_eq_const_to_auxiliary_const_m_pure = mobile_x_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, mobile_x_eq_const_to_auxiliary_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/r_eq_const_to_auxiliary_const_m.npy"
        r_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, NUM_BANDS*NUM_TIME_STEPS))
        r_eq_const_to_auxiliary_const_m_pure = r_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, r_eq_const_to_auxiliary_const_m_pure.flatten())
    
        auxiliary_constraint.auxiliary_const_to_z_eq_const_update(mobile_x_eq_const_to_auxiliary_const_m_pure, r_eq_const_to_auxiliary_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/auxiliary_const_to_z_eq_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        auxiliary_const_to_z_eq_const_m_pure = auxiliary_constraint.auxiliary_const_to_z_eq_const_m.astype(np_impa_lib)
        np.save(output_file_1, auxiliary_const_to_z_eq_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_input3 = os.getcwd() + "/../ut_inputs/z_eq_const_to_auxiliary_const_m.npy"
        z_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
        z_eq_const_to_auxiliary_const_m_pure = z_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
        np.save(f_input3, z_eq_const_to_auxiliary_const_m_pure.flatten())
        
        auxiliary_constraint.auxiliary_const_to_r_eq_const_update(z_eq_const_to_auxiliary_const_m_pure, mobile_x_eq_const_to_auxiliary_const_m_pure)
        
        f_output_2 = os.getcwd() + "/../ut_results/auxiliary_const_to_r_eq_const_m_pure"
        output_file_2 = open(f_output_2, "wb")
        auxiliary_const_to_r_eq_const_m_pure = auxiliary_constraint.auxiliary_const_to_r_eq_const_m.astype(np_impa_lib)
        np.save(output_file_2, auxiliary_const_to_r_eq_const_m_pure.flatten(), allow_pickle=True)
        output_file_2.close()
        
    elif ut_name == "AuxiliaryConst2MobileXEqConstUpdate":
        
        f_input1 = os.getcwd() + "/../ut_inputs/z_eq_const_to_auxiliary_const_m.npy"
        z_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS))
        z_eq_const_to_auxiliary_const_m_pure = z_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, z_eq_const_to_auxiliary_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/r_eq_const_to_auxiliary_const_m.npy"
        r_eq_const_to_auxiliary_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS, NUM_BANDS*NUM_TIME_STEPS))
        r_eq_const_to_auxiliary_const_m_pure = r_eq_const_to_auxiliary_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, r_eq_const_to_auxiliary_const_m_pure.flatten())
        
        auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_update(z_eq_const_to_auxiliary_const_m_pure, r_eq_const_to_auxiliary_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/auxiliary_const_to_mobile_x_eq_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        auxiliary_const_to_mobile_x_eq_const_m_pure = auxiliary_constraint.auxiliary_const_to_mobile_x_eq_const_m.astype(np_impa_lib)
        np.save(output_file_1, auxiliary_const_to_mobile_x_eq_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        