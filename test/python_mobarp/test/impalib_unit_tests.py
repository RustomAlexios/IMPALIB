# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ut_update_auxiliary_constraint
import ut_update_equality_constraint
import ut_update_inequality_constraint
import ut_graphical_model
import ut_input_output
from environmentModule import *

def costs_connectivity_generation(num_fixed_tx, num_bands, num_mobile_tx_locs, num_time_steps, num_mobile_tx, num_rx_locs, percentage_neg):
    #generate IM and connectivities
    normal_mean = 0
    normal_variance = 3
    
    fixed_x_costs = np.random.uniform(10, 100, size=(num_fixed_tx, num_bands, num_time_steps))
    mobile_x_costs = np.random.uniform(10, 100, size=(num_mobile_tx, num_bands, num_time_steps))
    r_costs = np.random.normal(normal_mean, normal_variance, size=(num_mobile_tx, num_mobile_tx_locs))
    
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

    z_costs = np.zeros((mobile_x_costs.size, num_mobile_tx_locs), dtype = np_impa_lib)

    connectivity_fixed_tx = np.random.randint(2, size=(num_fixed_tx, num_bands, num_time_steps, num_rx_locs), dtype=int)
    connectivity_mobile_tx = np.random.randint(2, size=(num_mobile_tx_locs, num_bands, num_time_steps, num_rx_locs), dtype=int)

    return fixed_x_costs, mobile_x_costs, connectivity_fixed_tx, connectivity_mobile_tx, r_costs, z_costs
      
parser = argparse.ArgumentParser()
parser.add_argument("--sub_test_num", type=int, default=1, help="Sub-Test of Unit Test")
parser.add_argument("--sub_tests_total", type=int, default=1, help="Total Number of Sub-Tests of Unit Test")
parser.add_argument("--ut_name", type=str, help="Unit Test Name")
parser.add_argument("--nITER", type=int, default=1, help="Number of Iterations of IMPA")
parser.add_argument("--nFTX", type=int, default=2, help="Number of fixed TX")
parser.add_argument("--nMTX", type=int, default=2, help="Number of mobile TX")
parser.add_argument("--nBands", type=int, default=3, help="Number of bands")
parser.add_argument("--nTimeSteps", type=int, default=4, help="Number of discrete time steps")
parser.add_argument("--nRX", type=int, default=3, help="Number of RX")
parser.add_argument("--nMTXLoc", type=int, default=4, help="Number of mobile TX locations")
parser.add_argument("--filteringFlag", type=bool, help="Activate Filtering or not")
parser.add_argument("--overWriteCapFlag", type=bool, default=False, help="Over Write Max Capacity or not")
parser.add_argument("--overWriteCapVal", type=int, default=4, help="Over Write Max Capacity Value")
parser.add_argument("--alpha", type=np_impa_lib, default=0.0, help="Filtering Rate [0,1]")
parser.add_argument("--testFile", type=int, default=9000, help="Test File Index")
parser.add_argument("--saveFlag", type=bool, default=False, help="Save Outputs or not")
parser.add_argument("--excludeCapFlag", type=bool, default=False, help="Exclude capacities or not")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument("--getSolApproach", type=int, default=2, help="Approach for getting a solution")
parser.add_argument("--criteriaIM", type=int, default=2, help="1: normal, 2: pos X, normal R")
parser.add_argument("--percNegIM", type=np_impa_lib, default=0, help="Percentage of Negative samples in IM")
parser.add_argument("--overWriteIM", type=bool, default=False, help="overwrite IM")

if __name__ == "__main__":
    args = parser.parse_args()
    
    sub_test_num = args.sub_test_num
    total_sub_tests = args.sub_tests_total
    ut_name = args.ut_name
    
    NUM_ITERATIONS = args.nITER # number of iterations of IMPA
    NUM_FIXED_TX = args.nFTX # number of fixed TX
    NUM_MOBILE_TX = args.nMTX # number of mobile TX
    NUM_BANDS = args.nBands # number of bands
    NUM_TIME_STEPS = args.nTimeSteps # number of time steps
    NUM_RX_LOCS = args.nRX # number of RX
    NUM_MOBILE_TX_LOCS = args.nMTXLoc # number of mobile TX locs
    THRESHOLD = args.threshold # threshold for hard decision analysis
    FILTERING_FLAG = args.filteringFlag # perform filtering on messages from degree constraints and subtour elimination constraints
    OVER_WRITE_CAP_FLAG = args.overWriteCapFlag # overwrite capacities flag
    OVER_WRITE_CAP_VAL = args.overWriteCapVal # overwrite capacities value
    ALPHA = args.alpha # filtering parameter
    EXCLUDE_CAP_FLAG = args.excludeCapFlag # exclude capacity constraints flag
    GET_SOL_APPROACH = args.getSolApproach # approach (for now, 1 or 2) to get solution with minimum objective (2 is more efficient)
    CRITERIA_IM = args.criteriaIM
    PERCENTAGE_NEGATIVE_IM = args.percNegIM
    OVERWRITE_IM = args.overWriteIM
    
    # print("FILTERING_FLAG: ", FILTERING_FLAG)

    if (ut_name != "Iterate"):
        fixed_x_costs, mobile_x_costs, connectivity_fixed_tx, connectivity_mobile_tx, r_costs, z_costs = costs_connectivity_generation(NUM_FIXED_TX, NUM_BANDS, NUM_MOBILE_TX_LOCS, NUM_TIME_STEPS, NUM_MOBILE_TX, NUM_RX_LOCS, PERCENTAGE_NEGATIVE_IM)
            
        if (not EXCLUDE_CAP_FLAG):
            max_cap_val = NUM_BANDS
            fixed_capac_constraints = np.random.randint(low = 1, high=max_cap_val, size=NUM_FIXED_TX, dtype=int)
            mobile_capac_constraints = np.random.randint(low = 1, high=max_cap_val, size=NUM_MOBILE_TX, dtype=int)
        else:
            fixed_capac_constraints = NUM_BANDS*np.ones(NUM_FIXED_TX, dtype=int)
            mobile_capac_constraints = NUM_BANDS*np.ones(NUM_MOBILE_TX, dtype=int)
            
    if ut_name == "ExtrinsicUpdate":
        ut_input_output.ut_io(
            ut_name=ut_name,
            num_fixed_tx = NUM_FIXED_TX,
            num_mobile_tx = NUM_MOBILE_TX,
            num_bands = NUM_BANDS,
            num_time_steps = NUM_TIME_STEPS,
            num_rx_locs = NUM_RX_LOCS,
            num_mobile_tx_locs = NUM_MOBILE_TX_LOCS,
            exclude_cap_flag = EXCLUDE_CAP_FLAG,
            fixed_x_costs = fixed_x_costs,
            mobile_x_costs = mobile_x_costs,
            connectivity_fixed_tx = connectivity_fixed_tx,
            connectivity_mobile_tx = connectivity_mobile_tx,
            r_costs = r_costs,
            z_costs = z_costs,
        )
    elif ut_name == "IneqCapacConstUpdate" or ut_name == "MobileLocEqConst2REqConstUpdate" or ut_name == "SetCoverIneqConstUpdate":
        ut_update_inequality_constraint.ut_ineq_capac_const_update(
            ut_name, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_capac_constraints, mobile_capac_constraints, 
            ALPHA, FILTERING_FLAG, NUM_MOBILE_TX_LOCS, NUM_RX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx)
    # elif ut_name == "AuxiliaryConst2ZEqConstUpdate" or ut_name == "AuxiliaryConst2REqConstUpdate" or ut_name == "AuxiliaryConst2MobileXEqConstUpdate" or ut_name == "AuxiliaryConst2ZAndREqConstUpdate":
    elif ut_name == "AuxiliaryConst2ZAndREqConstUpdate" or ut_name == "AuxiliaryConst2MobileXEqConstUpdate":
        ut_update_auxiliary_constraint.ut_auxiliary_constraint(
            ut_name, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_MOBILE_TX_LOCS, ALPHA, FILTERING_FLAG, connectivity_mobile_tx)
    elif ut_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate" or ut_name == "ReqConstActivation" or ut_name == "ZEqConst2AuxiliaryConstUpdate" or ut_name == "XEqConstActivation" or ut_name == "ZEqConst2SetCoverIneqConstUpdate":
        ut_update_equality_constraint.ut_equality_constraint(
                    ut_name,NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, fixed_x_costs, mobile_x_costs, NUM_MOBILE_TX_LOCS, 
                        z_costs, NUM_RX_LOCS, connectivity_mobile_tx, connectivity_fixed_tx, r_costs, EXCLUDE_CAP_FLAG)
    elif ut_name == "Iterate":
        RANDOM_TEST_FLAG = True
        POST_PROCESS_FLAG = False
        ut_graphical_model.ut_model_graph(ut_name, NUM_ITERATIONS, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS,
                                NUM_MOBILE_TX_LOCS, THRESHOLD, FILTERING_FLAG, ALPHA, RANDOM_TEST_FLAG, POST_PROCESS_FLAG, OVER_WRITE_CAP_FLAG, OVER_WRITE_CAP_VAL, EXCLUDE_CAP_FLAG, GET_SOL_APPROACH,
                                    CRITERIA_IM, PERCENTAGE_NEGATIVE_IM, OVERWRITE_IM)