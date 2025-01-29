# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
import graphical_model as model_graph
from ut_utils import *


def ut_model_graph(ut_name, NUM_ITERATIONS,
                NUM_FIXED_TX,
                NUM_MOBILE_TX,
                NUM_BANDS,
                NUM_TIME_STEPS,
                NUM_RX_LOCS,
                NUM_MOBILE_TX_LOCS,
                THRESHOLD,
                FILTERING_FLAG,
                ALPHA,
                RANDOM_TEST_FLAG,
                POST_PROCESS_FLAG,
                OVER_WRITE_CAP_FLAG, 
                OVER_WRITE_CAP_VAL,
                EXCLUDE_CAP_FLAG,
                GET_SOL_APPROACH,
                CRITERIA_IM,
                PERCENTAGE_NEGATIVE_IM,
                OVERWRITE_IM):
    
    print("-------")
    print("Python Pure")

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, ALPHA)

    f_input_threshold = os.getcwd() + "/../ut_inputs/threshold.npy"
    np.save(f_input_threshold, THRESHOLD)

    if ut_name == "Iterate":
        
        ModelIMPA = model_graph.GraphicalModelMOBARP(
            NUM_ITERATIONS,
            NUM_FIXED_TX,
            NUM_MOBILE_TX,
            NUM_BANDS,
            NUM_TIME_STEPS,
            NUM_RX_LOCS,
            NUM_MOBILE_TX_LOCS,
            THRESHOLD,
            FILTERING_FLAG,
            ALPHA,
            RANDOM_TEST_FLAG,
            POST_PROCESS_FLAG,
            OVER_WRITE_CAP_FLAG, 
            OVER_WRITE_CAP_VAL,
            EXCLUDE_CAP_FLAG,
            GET_SOL_APPROACH,
            CRITERIA_IM,
            PERCENTAGE_NEGATIVE_IM,
            OVERWRITE_IM
        )

        ModelIMPA.input_load = []
        formatted_alpha = "{:.1f}".format(ALPHA)
        ModelIMPA.formatted_alpha = formatted_alpha
        
        ModelIMPA.initialize()
        
        f_input1 = os.getcwd() + "/../ut_inputs/fixed_x_costs.npy"
        fixed_x_costs = ModelIMPA.fixed_x_costs.astype(np_impa_lib)
        np.save(f_input1, fixed_x_costs.flatten())
    
        f_input2 = os.getcwd() + "/../ut_inputs/mobile_x_costs.npy"
        mobile_x_costs = ModelIMPA.mobile_x_costs.astype(np_impa_lib)
        np.save(f_input2, mobile_x_costs.flatten())
        
        f_input3 = os.getcwd() + "/../ut_inputs/z_costs.npy"
        z_costs = ModelIMPA.z_costs.astype(np_impa_lib)
        np.save(f_input3, z_costs.flatten())
    
        f_input4 = os.getcwd() + "/../ut_inputs/r_costs.npy"
        r_costs = ModelIMPA.r_costs.astype(np_impa_lib)
        np.save(f_input4, r_costs.flatten())
        
        f_input5 = os.getcwd() + "/../ut_inputs/connectivity_mobile_tx.npy"
        connectivity_mobile_tx = ModelIMPA.connectivity_mobile_tx.astype(np.int32)
        np.save(f_input5, connectivity_mobile_tx.flatten())
        
        f_input6 = os.getcwd() + "/../ut_inputs/connectivity_fixed_tx.npy"
        connectivity_fixed_tx = ModelIMPA.connectivity_fixed_tx.astype(np.int32)
        np.save(f_input6, connectivity_fixed_tx.flatten())
        
        f_input7 = os.getcwd() + "/../ut_inputs/conx_mob_tx_per_num_mob_tx_locs.npy"
        conx_mob_tx_per_num_mob_tx_locs = ModelIMPA.outputs.conx_mob_tx_per_num_mob_tx_locs.astype(np.int32)
        np.save(f_input7, conx_mob_tx_per_num_mob_tx_locs.flatten())
        
        f_input8 = os.getcwd() + "/../ut_inputs/conx_fixed_tx_per_num_rx_locs.npy"
        conx_fixed_tx_per_num_rx_locs = ModelIMPA.outputs.conx_fixed_tx_per_num_rx_locs.astype(np.int32)
        np.save(f_input8, conx_fixed_tx_per_num_rx_locs.flatten())
        
        f_input9 = os.getcwd() + "/../ut_inputs/conx_mob_tx_rx.npy"
        conx_mob_tx_rx = ModelIMPA.outputs.conx_mob_tx_rx.astype(np.int32)
        np.save(f_input9, conx_mob_tx_rx.flatten())
    
        f_input10 = os.getcwd() + "/../ut_inputs/fixed_capacity_constraints.npy"
        fixed_capac_constraints = ModelIMPA.fixed_capac_constraints.astype(np.int32)
        np.save(f_input10, fixed_capac_constraints.flatten())
        
        f_input11 = os.getcwd() + "/../ut_inputs/mobile_capacity_constraints.npy"
        mobile_capac_constraints = ModelIMPA.mobile_capac_constraints.astype(np.int32)
        np.save(f_input11, mobile_capac_constraints.flatten())
        
        if (not EXCLUDE_CAP_FLAG):
            f_input12 = os.getcwd() + "/../ut_inputs/fixed_x_eq_const_to_fixed_capac_const_m.npy"
            fixed_x_eq_const_to_fixed_capac_const_m_pure = ModelIMPA.fixed_x_eq_const_to_fixed_capac_const_m.astype(np_impa_lib)
            np.save(f_input12, fixed_x_eq_const_to_fixed_capac_const_m_pure.flatten())
            
            f_input13 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_mobile_capac_const_m.npy"
            mobile_x_eq_const_to_mobile_capac_const_m_pure = ModelIMPA.mobile_x_eq_const_to_mobile_capac_const_m.astype(np_impa_lib)
            np.save(f_input13, mobile_x_eq_const_to_mobile_capac_const_m_pure.flatten())
        else:
            f_input12 = os.getcwd() + "/../ut_inputs/fixed_x_eq_const_to_fixed_capac_const_m.npy"
            fixed_x_eq_const_to_fixed_capac_const_m_pure = np.zeros((NUM_FIXED_TX*NUM_BANDS*NUM_TIME_STEPS), dtype = np_impa_lib)
            np.save(f_input12, fixed_x_eq_const_to_fixed_capac_const_m_pure.flatten())
            f_input13 = os.getcwd() + "/../ut_inputs/mobile_x_eq_const_to_mobile_capac_const_m.npy"
            mobile_x_eq_const_to_mobile_capac_const_m_pure = np.zeros((NUM_MOBILE_TX*NUM_BANDS*NUM_BANDS), dtype = np_impa_lib)
            np.save(f_input13, mobile_x_eq_const_to_mobile_capac_const_m_pure.flatten())
            
        f_input14 = os.getcwd() + "/../ut_inputs/r_eq_const_to_auxiliary_const_m.npy"
        r_eq_const_to_auxiliary_const_m_pure = ModelIMPA.model_eq_constraint.r_eq_const_to_auxiliary_const_m.astype(np_impa_lib)
        np.save(f_input14, r_eq_const_to_auxiliary_const_m_pure.flatten())
        
        ModelIMPA.start_time = time.time()
        
        ModelIMPA.run_impa()
        
        f_output_1 = os.getcwd() + "/../ut_results/extrinsic_fixed_x_pure"
        output_file_1 = open(f_output_1, "wb")
        extrinsic_fixed_x_py = ModelIMPA.extrinsic_fixed_x.astype(np_impa_lib)
        np.save(output_file_1, extrinsic_fixed_x_py.flatten(), allow_pickle=True)
        output_file_1.close()
        
        f_output_2 = os.getcwd() + "/../ut_results/extrinsic_mobile_x_pure"
        extrinsic_mobile_x_py = ModelIMPA.extrinsic_mobile_x.astype(np_impa_lib)
        output_file_2 = open(f_output_2, "wb")
        np.save(output_file_2, extrinsic_mobile_x_py.flatten(), allow_pickle=True)
        output_file_2.close()
        
        f_output_3 = os.getcwd() + "/../ut_results/extrinsic_r_pure"
        extrinsic_r_py = ModelIMPA.extrinsic_r.astype(np_impa_lib)
        output_file_3 = open(f_output_3, "wb")
        np.save(output_file_3, extrinsic_r_py.flatten(), allow_pickle=True)
        output_file_3.close()
        
        f_output_4 = os.getcwd() + "/../ut_results/extrinsic_z_pure"
        extrinsic_z_py = ModelIMPA.extrinsic_z.astype(np_impa_lib)
        output_file_4 = open(f_output_4, "wb")
        np.save(output_file_4, extrinsic_z_py.flatten(), allow_pickle=True)
        output_file_4.close()
