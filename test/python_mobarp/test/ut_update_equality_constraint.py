# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
import update_equality_constraint as equality_constraint
from ut_utils import *


def ut_equality_constraint(ut_name,num_fixed_tx, num_mobile_tx, num_bands, num_time_steps, fixed_x_costs, mobile_x_costs, num_mobile_tx_locs, 
                                z_costs, num_rx_locs, connectivity_mobile_tx, connectivity_fixed_tx, r_costs, exclude_cap_flag):
    
    NUM_FIXED_TX = num_fixed_tx
    NUM_MOBILE_TX = num_mobile_tx
    NUM_BANDS = num_bands
    NUM_TIME_STEPS = num_time_steps
    NUM_RX_LOCS = num_rx_locs
    NUM_MOBILE_TX_LOCS = num_mobile_tx_locs
    EXCLUDE_CAP_FLAG = exclude_cap_flag
    
    model_equality_constraint = equality_constraint.EqualityConstraintMOBARP(NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_x_costs, mobile_x_costs, NUM_MOBILE_TX_LOCS, z_costs, 
                                                                            NUM_RX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx, r_costs, EXCLUDE_CAP_FLAG)
    
    f_input1 = os.getcwd() + "/../ut_inputs/fixed_x_costs.npy"
    fixed_x_costs = fixed_x_costs.astype(np_impa_lib)
    np.save(f_input1, fixed_x_costs.flatten())
    
    f_input2 = os.getcwd() + "/../ut_inputs/mobile_x_costs.npy"
    mobile_x_costs = mobile_x_costs.astype(np_impa_lib)
    np.save(f_input2, mobile_x_costs.flatten())
    
    f_input3 = os.getcwd() + "/../ut_inputs/z_costs.npy"
    z_costs = z_costs.astype(np_impa_lib)
    np.save(f_input3, z_costs.flatten())
    
    f_input4 = os.getcwd() + "/../ut_inputs/r_costs.npy"
    r_costs = r_costs.astype(np_impa_lib)
    np.save(f_input4, r_costs.flatten())
    
    f_input5 = os.getcwd() + "/../ut_inputs/connectivity_mobile_tx.npy"
    connectivity_mobile_tx = connectivity_mobile_tx.astype(np.int32)
    np.save(f_input5, connectivity_mobile_tx.flatten())
    
    f_input6 = os.getcwd() + "/../ut_inputs/connectivity_fixed_tx.npy"
    connectivity_fixed_tx = connectivity_fixed_tx.astype(np.int32)
    np.save(f_input6, connectivity_fixed_tx.flatten())
    
    f_input7 = os.getcwd() + "/../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy"
    conx_mob_tx_per_num_mob_tx_locs = model_equality_constraint.conx_mob_tx_per_num_mob_tx_locs.astype(np.int32)
    np.save(f_input7, conx_mob_tx_per_num_mob_tx_locs.flatten())
    
    f_input8 = os.getcwd() + "/../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy"
    conx_fixed_tx_per_num_rx_locs = model_equality_constraint.conx_fixed_tx_per_num_rx_locs.astype(np.int32)
    np.save(f_input8, conx_fixed_tx_per_num_rx_locs.flatten())
    
    f_input9 = os.getcwd() + "/../ut_inputs/conx_mob_tx_rx.npy"
    conx_mob_tx_rx = model_equality_constraint.conx_mob_tx_rx.astype(np.int32)
    np.save(f_input9, conx_mob_tx_rx.flatten())
    
    normal_mean = 0
    normal_variance = 100
        
    if (ut_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate"):
        
        f_input1 = os.getcwd() + "/../ut_inputs/fixed_capac_const_to_fixed_x_eq_const_m_pure.npy"
        fixed_capac_const_to_fixed_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))
        fixed_capac_const_to_fixed_x_eq_const_m_pure = fixed_capac_const_to_fixed_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, fixed_capac_const_to_fixed_x_eq_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/mobile_capac_const_to_mobile_x_eq_const_m_pure.npy"
        mobile_capac_const_to_mobile_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))
        mobile_capac_const_to_mobile_x_eq_const_m_pure = mobile_capac_const_to_mobile_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, mobile_capac_const_to_mobile_x_eq_const_m_pure.flatten())
        
        f_input3 = os.getcwd() + "/../ut_inputs/auxiliary_const_to_mobile_x_eq_const_m_pure.npy"
        auxiliary_const_to_mobile_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(mobile_x_costs.size, NUM_MOBILE_TX_LOCS))
        auxiliary_const_to_mobile_x_eq_const_m_pure = auxiliary_const_to_mobile_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input3, auxiliary_const_to_mobile_x_eq_const_m_pure.flatten())
        
        f_input4 = os.getcwd() + "/../ut_inputs/set_cover_ineq_const_to_fixed_x_eq_const_m_pure.npy"
        set_cover_ineq_const_to_fixed_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_TIME_STEPS, NUM_RX_LOCS, NUM_FIXED_TX*NUM_BANDS))
        set_cover_ineq_const_to_fixed_x_eq_const_m_pure = set_cover_ineq_const_to_fixed_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input4, set_cover_ineq_const_to_fixed_x_eq_const_m_pure.flatten())
        
        model_equality_constraint.x_eq_const_to_auxiliary_and_set_cover_const_update(fixed_capac_const_to_fixed_x_eq_const_m_pure, mobile_capac_const_to_mobile_x_eq_const_m_pure, auxiliary_const_to_mobile_x_eq_const_m_pure, set_cover_ineq_const_to_fixed_x_eq_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/mobile_x_eq_const_to_auxiliary_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        mobile_x_eq_const_to_auxiliary_const_m_pure = model_equality_constraint.mobile_x_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        np.save(output_file_1, mobile_x_eq_const_to_auxiliary_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/fixed_x_eq_const_to_set_cover_const_m_pure"
        output_file_2 = open(f_output_2, "wb")
        fixed_x_eq_const_to_set_cover_const_m_pure = model_equality_constraint.fixed_x_eq_const_to_set_cover_const_m.astype(np_impa_lib)
        np.save(output_file_2, fixed_x_eq_const_to_set_cover_const_m_pure.flatten(), allow_pickle=True)
        output_file_2.close()
        
    elif (ut_name == "ReqConstActivation"):
        
        f_input1 = os.getcwd() + "/../ut_inputs/auxiliary_const_to_r_eq_const_m_pure.npy"
        auxiliary_const_to_r_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(mobile_x_costs.size, NUM_MOBILE_TX_LOCS))
        auxiliary_const_to_r_eq_const_m_pure = auxiliary_const_to_r_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, auxiliary_const_to_r_eq_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/mobile_loc_eq_const_to_r_eq_const_m_pure.npy"
        mobile_loc_eq_const_to_r_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS))
        mobile_loc_eq_const_to_r_eq_const_m_pure = mobile_loc_eq_const_to_r_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, mobile_loc_eq_const_to_r_eq_const_m_pure.flatten())
        
        model_equality_constraint.r_eq_const_activation(auxiliary_const_to_r_eq_const_m_pure, mobile_loc_eq_const_to_r_eq_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/r_eq_const_to_auxiliary_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        r_eq_const_to_auxiliary_const_m_pure = model_equality_constraint.r_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        np.save(output_file_1, r_eq_const_to_auxiliary_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/r_eq_const_to_mobile_loc_eq_const_m_pure"
        output_file_2 = open(f_output_2, "wb")
        r_eq_const_to_mobile_loc_eq_const_m_pure = model_equality_constraint.r_eq_const_to_mobile_loc_eq_const_m.astype(np_impa_lib)
        np.save(output_file_2, r_eq_const_to_mobile_loc_eq_const_m_pure.flatten(), allow_pickle=True)
        output_file_2.close()
        
    elif ut_name == "ZEqConst2AuxiliaryConstUpdate":
        
        f_input1 = os.getcwd() + "/../ut_inputs/set_cover_ineq_const_to_z_eq_const_m_pure.npy"
        set_cover_ineq_const_to_z_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, NUM_BANDS*NUM_MOBILE_TX))
        set_cover_ineq_const_to_z_eq_const_m_pure = set_cover_ineq_const_to_z_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, set_cover_ineq_const_to_z_eq_const_m_pure.flatten())
        
        model_equality_constraint.z_eq_const_to_auxiliary_const_update(set_cover_ineq_const_to_z_eq_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/z_eq_const_to_auxiliary_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        z_eq_const_to_auxiliary_const_m_pure = model_equality_constraint.z_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        np.save(output_file_1, z_eq_const_to_auxiliary_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()

    elif ut_name == "XEqConstActivation":
        
        f_input1 = os.getcwd() + "/../ut_inputs/auxiliary_const_to_mobile_x_eq_const_m_pure.npy"
        auxiliary_const_to_mobile_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(mobile_x_costs.size, NUM_MOBILE_TX_LOCS))
        auxiliary_const_to_mobile_x_eq_const_m_pure = auxiliary_const_to_mobile_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, auxiliary_const_to_mobile_x_eq_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/set_cover_ineq_const_to_fixed_x_eq_const_m_pure.npy"
        set_cover_ineq_const_to_fixed_x_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_TIME_STEPS, NUM_RX_LOCS, NUM_FIXED_TX*NUM_BANDS))
        set_cover_ineq_const_to_fixed_x_eq_const_m_pure = set_cover_ineq_const_to_fixed_x_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, set_cover_ineq_const_to_fixed_x_eq_const_m_pure.flatten())
        
        mobile_x_eq_const_to_mobile_capac_const_m, fixed_x_eq_const_to_fixed_capac_const_m = model_equality_constraint.x_eq_const_activation(auxiliary_const_to_mobile_x_eq_const_m_pure, set_cover_ineq_const_to_fixed_x_eq_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/mobile_x_eq_const_to_mobile_capac_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        mobile_x_eq_const_to_mobile_capac_const_m_pure = mobile_x_eq_const_to_mobile_capac_const_m.astype(np_impa_lib)
        np.save(output_file_1, mobile_x_eq_const_to_mobile_capac_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/fixed_x_eq_const_to_fixed_capac_const_m_pure"
        output_file_2 = open(f_output_2, "wb")
        fixed_x_eq_const_to_fixed_capac_const_m_pure = fixed_x_eq_const_to_fixed_capac_const_m.astype(np_impa_lib)
        np.save(output_file_2, fixed_x_eq_const_to_fixed_capac_const_m_pure.flatten(), allow_pickle=True)
        output_file_2.close()
        
    elif ut_name == "ZEqConst2SetCoverIneqConstUpdate":
        
        f_input1 = os.getcwd() + "/../ut_inputs/auxiliary_const_to_z_eq_const_m_pure.npy"
        auxiliary_const_to_z_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(mobile_x_costs.size, NUM_MOBILE_TX_LOCS))
        auxiliary_const_to_z_eq_const_m_pure = auxiliary_const_to_z_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input1, auxiliary_const_to_z_eq_const_m_pure.flatten())
        
        f_input2 = os.getcwd() + "/../ut_inputs/set_cover_ineq_const_to_z_eq_const_m_pure.npy"
        set_cover_ineq_const_to_z_eq_const_m_pure = np.random.normal(normal_mean, normal_variance, size=(NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, NUM_BANDS*NUM_MOBILE_TX))
        set_cover_ineq_const_to_z_eq_const_m_pure = set_cover_ineq_const_to_z_eq_const_m_pure.astype(np_impa_lib)
        np.save(f_input2, set_cover_ineq_const_to_z_eq_const_m_pure.flatten())
        
        model_equality_constraint.z_eq_const_to_set_cover_ineq_const_update(auxiliary_const_to_z_eq_const_m_pure, set_cover_ineq_const_to_z_eq_const_m_pure)
        
        f_output_1 = os.getcwd() + "/../ut_results/z_eq_const_to_set_cover_ineq_const_m_pure"
        output_file_1 = open(f_output_1, "wb")
        z_eq_const_to_set_cover_ineq_const_m_pure = model_equality_constraint.z_eq_const_to_set_cover_ineq_const_m.astype(np_impa_lib)
        np.save(output_file_1, z_eq_const_to_set_cover_ineq_const_m_pure.flatten(), allow_pickle=True)
        output_file_1.close()
        