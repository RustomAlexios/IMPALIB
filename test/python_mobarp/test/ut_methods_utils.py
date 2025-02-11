# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + "/../src")
from ut_utils import *


parser = argparse.ArgumentParser()
parser.add_argument("--sub_test_num", type=int, default=1, help="Sub-Test of Unit Test")
parser.add_argument("--sub_tests_total", type=int, default=1, help="Total Number of Sub-Tests of Unit Test")
parser.add_argument("--ut_name", type=str, help="Unit Test Name")

if __name__ == "__main__":
    args = parser.parse_args()
    sub_test_num = args.sub_test_num
    total_sub_tests = args.sub_tests_total
    ut_name = args.ut_name
    rtol = 5e-02
    atol = 1e-04  # user will specify these tolerances

    if ut_name == "ExtrinsicUpdate":
        f_path_pure_list = ["../ut_results/extrinsic_fixed_x_pure", "../ut_results/extrinsic_mobile_x_pure", "../ut_results/extrinsic_r_pure", "../ut_results/extrinsic_z_pure"]
        f_path_wrapper_list = ["../ut_results/extrinsic_fixed_x_wrapper", "../ut_results/extrinsic_mobile_x_wrapper", "../ut_results/extrinsic_r_wrapper", "../ut_results/extrinsic_z_wrapper"]
    elif ut_name == "IneqCapacConstUpdate":
        f_path_pure_list = ["../ut_results/fixed_capac_const_to_fixed_x_eq_const_m_dummy_pure", "../ut_results/mobile_capac_const_to_mobile_x_eq_const_m_dummy_pure"]
        f_path_wrapper_list = ["../ut_results/FixedCapacConst2FixedXEqConstDummyM_wrapper", "../ut_results/MobileCapacConst2MobileXEqConstDummyM_wrapper"]
    elif ut_name == "MobileLocEqConst2REqConstUpdate":
        f_path_pure = "../ut_results/mobile_loc_eq_const_to_r_eq_const_m_dummy_pure"
        f_path_wrapper = "../ut_results/MobileLocEqConst2REqConstDummyM_wrapper"
    elif ut_name == "SetCoverIneqConstUpdate":
        f_path_pure_list = ["../ut_results/set_cover_ineq_const_to_fixed_x_eq_const_m_dummy_pure", "../ut_results/set_cover_ineq_const_to_z_eq_const_m_dummy_pure"]
        f_path_wrapper_list = ["../ut_results/SetCoverIneqConst2FixedXEqConstDummyM_wrapper", "../ut_results/SetCoverIneqConst2ZEqConstDummyM_wrapper"]
    elif ut_name == "AuxiliaryConst2ZEqConstUpdate":
        f_path_pure = "../ut_results/auxiliary_const_to_z_eq_const_m_pure"
        f_path_wrapper = "../ut_results/AuxiliaryConst2ZEqConstM_wrapper" 
    elif ut_name == "AuxiliaryConst2REqConstUpdate":
        f_path_pure = "../ut_results/auxiliary_const_to_r_eq_const_m_pure"
        f_path_wrapper = "../ut_results/AuxiliaryConst2REqConstM_wrapper" 
    elif ut_name == "AuxiliaryConst2ZAndREqConstUpdate":
        f_path_pure_list = ["../ut_results/auxiliary_const_to_z_eq_const_m_pure", "../ut_results/auxiliary_const_to_r_eq_const_m_pure"]
        f_path_wrapper_list = ["../ut_results/AuxiliaryConst2ZEqConstM_wrapper", "../ut_results/AuxiliaryConst2REqConstM_wrapper"] 
    elif ut_name == "AuxiliaryConst2MobileXEqConstUpdate":
        f_path_pure = "../ut_results/auxiliary_const_to_mobile_x_eq_const_m_pure"
        f_path_wrapper = "../ut_results/AuxiliaryConst2MobileXEqConstM_wrapper" 
    elif ut_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate":
        f_path_pure_list = ["../ut_results/mobile_x_eq_const_to_auxiliary_const_m_pure", "../ut_results/fixed_x_eq_const_to_set_cover_const_m_pure"]
        f_path_wrapper_list = ["../ut_results/MobileXEqConst2AuxiliaryConstM_wrapper", "../ut_results/FixedXEqConst2SetCoverConstM_wrapper"]
    elif ut_name == "ReqConstActivation":
        f_path_pure_list = ["../ut_results/r_eq_const_to_auxiliary_const_m_pure", "../ut_results/r_eq_const_to_mobile_loc_eq_const_m_pure"]
        f_path_wrapper_list = ["../ut_results/REqConst2AuxiliaryConstM_wrapper", "../ut_results/REqConst2MobileLocEqConstM_wrapper"]
    elif ut_name == "ZEqConst2AuxiliaryConstUpdate":
        f_path_pure = "../ut_results/z_eq_const_to_auxiliary_const_m_pure"
        f_path_wrapper = "../ut_results/ZEqConst2AuxiliaryConstM_wrapper" 
    elif ut_name == "XEqConstActivation":
        f_path_pure_list = ["../ut_results/mobile_x_eq_const_to_mobile_capac_const_m_pure", "../ut_results/fixed_x_eq_const_to_fixed_capac_const_m_pure"]
        f_path_wrapper_list = ["../ut_results/MobileXEqConst2MobileCapacConstM_wrapper", "../ut_results/FixedXEqConst2FixedCapacConstM_wrapper"]
    elif ut_name == "ZEqConst2SetCoverIneqConstUpdate":
        f_path_pure = "../ut_results/z_eq_const_to_set_cover_ineq_const_m_pure"
        f_path_wrapper = "../ut_results/ZEqConst2SetCoverIneqConstM_wrapper" 
    elif ut_name == "Iterate":
        f_path_pure_list = ["../ut_results/extrinsic_fixed_x_pure", "../ut_results/extrinsic_mobile_x_pure", "../ut_results/extrinsic_r_pure", "../ut_results/extrinsic_z_pure"]
        f_path_wrapper_list = ["../ut_results/extrinsic_fixed_x_wrapper", "../ut_results/extrinsic_mobile_x_wrapper", "../ut_results/extrinsic_r_wrapper", "../ut_results/extrinsic_z_wrapper"]
        
    if (
        ut_name == "ExtrinsicUpdate"
        or ut_name == "IneqCapacConstUpdate"
        or ut_name == "SetCoverIneqConstUpdate"
        or ut_name == "AuxiliaryConst2ZAndREqConstUpdate"
        or ut_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate"
        or ut_name == "ReqConstActivation"
        or ut_name == "XEqConstActivation"
        or ut_name == "Iterate"
    ):
        
        for i in range(len(f_path_pure_list)):
            f_path_pure = f_path_pure_list[i]
            f_path_wrapper = f_path_wrapper_list[i]
            
            file_array_pure = open(f_path_pure, "rb")
            y_pure = np.load(file_array_pure)
            y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib)
            # print("y_pure: ", y_pure)
            # print("y_wrapper: ", y_wrapper)
            # y_wrapper = np.fromfile(f_path_wrapper, dtype=np.int32); #for integer y_wrapper
            check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=rtol, atol=atol)
            
    # elif (ut_name == "MobileLocEqConst2REqConstUpdate" or ut_name == "AuxiliaryConst2ZEqConstUpdate" or ut_name == "AuxiliaryConst2REqConstUpdate" or ut_name == "AuxiliaryConst2MobileXEqConstUpdate"
    #       or ut_name == "ZEqConst2AuxiliaryConstUpdate" or ut_name == "ZEqConst2SetCoverIneqConstUpdate"):
    elif (ut_name == "MobileLocEqConst2REqConstUpdate" or ut_name == "AuxiliaryConst2MobileXEqConstUpdate"
          or ut_name == "ZEqConst2AuxiliaryConstUpdate" or ut_name == "ZEqConst2SetCoverIneqConstUpdate"):
            file_array_pure = open(f_path_pure, "rb")
            y_pure = np.load(file_array_pure)
            y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib)
            # print("y_pure: ", y_pure)
            # print("y_wrapper: ", y_wrapper)
            # y_wrapper = np.fromfile(f_path_wrapper, dtype=np.int32); #for integer y_wrapper
            check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=rtol, atol=atol)
