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

    if ut_name == "ExtrinsicOutputVariableEcUpdate":
        f_path_pure = "../ut_results/extrinsic_output_variable_ec_pure"
        f_path_wrapper = "../ut_results/extrinsic_output_variable_ec_wrapper"
    elif ut_name == "VariableEc2KsatConstraintUpdate":
        f_path_pure = "../ut_results/variable_ec_to_ksat_constraint_m_pure"
        f_path_wrapper = "../ut_results/variable_ec_to_ksat_constraint_m_wrapper"
    elif ut_name == "KsatConstraint2VariableEcUpdate":
        f_path_pure = "../ut_results/ksat_constraint_to_eq_constraint_m_pure"
        f_path_wrapper = "../ut_results/ksat_constraint_to_eq_constraint_m_wrapper"
    elif ut_name == "Iterate":
        f_path_pure = "../ut_results/extrinsic_out_variable_ec_pure"
        f_path_wrapper = "../ut_results/extrinsic_out_variable_ec_wrapper"
        
    if (
        ut_name == "ExtrinsicOutputVariableEcUpdate"
        or ut_name == "VariableEc2KsatConstraintUpdate"
        or ut_name == "KsatConstraint2VariableEcUpdate"
        or ut_name == "Iterate"
    ):
        file_array_pure = open(f_path_pure, "rb")
        y_pure = np.load(file_array_pure)
        y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib)
        #print("y_pure: ", y_pure)
        #print("y_wrapper: ", y_wrapper)
        # y_wrapper = np.fromfile(f_path_wrapper, dtype=np.int32); #for integer y_wrapper
        check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=rtol, atol=atol)
