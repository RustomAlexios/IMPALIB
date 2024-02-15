# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")
from ut_utils import check_agreement, check_forward_backward_results
from environmentModule import argparse, np, np_impa_lib

parser = argparse.ArgumentParser()
parser.add_argument("--sub_test_num", type=int, default=1, help="Sub-Test of Unit Test")
parser.add_argument(
    "--sub_tests_total",
    type=int,
    default=1,
    help="Total Number of Sub-Tests of Unit Test",
)
parser.add_argument("--ut_name", type=str, help="Unit Test Name")

if __name__ == "__main__":
    args = parser.parse_args()
    sub_test_num = args.sub_test_num
    total_sub_tests = args.sub_tests_total
    ut_name = args.ut_name
    rtol = 5e-02
    atol = 1e-08  # user will specify these tolerances
    if ut_name == "KnapsackForward" or ut_name == "KnapsackBackward":
        check_forward_backward_results(sub_test_num, total_sub_tests, ut_name, rtol=rtol, atol=atol)

    elif ut_name == "ExtrinsicOutputDepartment":
        f_path_pure = "../ut_results/ut_Knapsack/ut_" + ut_name + "/extrinsic_output_department_pure"
        f_path_wrapper = "../ut_results/ut_Knapsack/ut_" + ut_name + "/extrinsic_output_department_wrapper"
    elif ut_name == "Team2KnapsackUpdate":
        f_path_pure = "../ut_results/ut_Knapsack/ut_" + ut_name + "/team_to_knapsack_m_pure"
        f_path_wrapper = "../ut_results/ut_Knapsack/ut_" + ut_name + "/team_to_knapsack_m_wrapper"
    elif ut_name == "ProcessExtrinsicOutputDepartment":
        f_path_pure = "../ut_results/ut_Knapsack/ut_" + ut_name + "/extrinsic_output_department_pure"
        f_path_wrapper = "../ut_results/ut_Knapsack/ut_" + ut_name + "/extrinsic_output_department_wrapper"
    elif ut_name == "Oric2TeamUpdate":
        f_path_pure = "../ut_results/ut_ORIC/ut_" + ut_name + "/oric_to_team_m_pure"
        f_path_wrapper = "../ut_results/ut_ORIC/ut_" + ut_name + "/oric_to_team_m_wrapper"
    elif ut_name == "Oric2ProjectEcUpdate":
        f_path_pure1 = "../ut_results/ut_ORIC/ut_" + ut_name + "/eq_constraint_to_project_m_pure"
        f_path_wrapper1 = "../ut_results/ut_ORIC/ut_" + ut_name + "/eq_constraint_to_project_m_wrapper"
        file_array_pure1 = open(f_path_pure1, "rb")
        y_pure1 = np.load(file_array_pure1)
        y_wrapper1 = np.fromfile(f_path_wrapper1, dtype=np_impa_lib)
        check_agreement(
            sub_test_num,
            total_sub_tests,
            ut_name,
            y_pure1,
            y_wrapper1,
            rtol=rtol,
            atol=atol,
        )
        f_path_pure2 = "../ut_results/ut_ORIC/ut_" + ut_name + "/oric_to_eq_constraint_m_pure"
        f_path_wrapper2 = "../ut_results/ut_ORIC/ut_" + ut_name + "/oric_to_eq_constraint_m_wrapper"
        file_array_pure2 = open(f_path_pure2, "rb")
        y_pure2 = np.load(file_array_pure2)
        y_wrapper2 = np.fromfile(f_path_wrapper2, dtype=np_impa_lib)
        check_agreement(
            sub_test_num,
            total_sub_tests,
            ut_name,
            y_pure2,
            y_wrapper2,
            rtol=rtol,
            atol=atol,
        )
    elif ut_name == "ExtrinsicOutputTeamUpdate":
        f_path_pure = "../ut_results/ut_InputOutput/ut_" + ut_name + "/extrinsic_output_team_pure"
        f_path_wrapper = "../ut_results/ut_InputOutput/ut_" + ut_name + "/extrinsic_output_team_wrapper"
    elif ut_name == "IntrinsicOutMwmUpdate":
        f_path_pure = "../ut_results/ut_InputOutput/ut_" + ut_name + "/intrinsic_out_mwm_pure"
        f_path_wrapper = "../ut_results/ut_InputOutput/ut_" + ut_name + "/intrinsic_out_mwm_wrapper"
    elif ut_name == "TeamEc2OricUpdate":
        f_path_pure = "../ut_results/ut_EqConstraint/ut_" + ut_name + "/team_to_oric_m_pure"
        f_path_wrapper = "../ut_results/ut_EqConstraint/ut_" + ut_name + "/team_to_oric_m_wrapper"
    elif ut_name == "ProjectEqConst2OricUpdate":
        f_path_pure = "../ut_results/ut_EqConstraint/ut_" + ut_name + "/eq_constraint_to_oric_m_pure"
        f_path_wrapper = "../ut_results/ut_EqConstraint/ut_" + ut_name + "/eq_constraint_to_oric_m_wrapper"
    elif ut_name == "ProjectInequalityConstraintUpdate":
        f_path_pure = "../ut_results/ut_projectIneqConstraint/ut_" + ut_name + "/project_to_eq_constraint_m_pure"
        f_path_wrapper = "../ut_results/ut_projectIneqConstraint/ut_" + ut_name + "/project_to_eq_constraint_m_wrapper"
    elif ut_name == "Iterate" or ut_name == "IterateSampleGraph":
        f_path_pure1 = "../ut_results/ut_GraphicalModel/ut_" + ut_name + "/extrinsic_output_team_pure"
        f_path_wrapper1 = "../ut_results/ut_GraphicalModel/ut_" + ut_name + "/extrinsic_output_team_wrapper"
        file_array_pure1 = open(f_path_pure1, "rb")
        y_pure1 = np.load(file_array_pure1)
        y_wrapper1 = np.fromfile(f_path_wrapper1, dtype=np_impa_lib)
        check_agreement(
            sub_test_num,
            total_sub_tests,
            ut_name,
            y_pure1,
            y_wrapper1,
            rtol=rtol,
            atol=atol,
        )
        f_path_pure2 = "../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_mwm_pure"
        f_path_wrapper2 = "../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_mwm_wrapper"
        file_array_pure2 = open(f_path_pure2, "rb")
        y_pure2 = np.load(file_array_pure2)
        y_wrapper2 = np.fromfile(f_path_wrapper2, dtype=np_impa_lib)
        check_agreement(
            sub_test_num,
            total_sub_tests,
            ut_name,
            y_pure2,
            y_wrapper2,
            rtol=rtol,
            atol=atol,
        )

    if (
        ut_name == "Oric2TeamUpdate"
        or ut_name == "ExtrinsicOutputTeamUpdate"
        or ut_name == "IntrinsicOutMwmUpdate"
        or ut_name == "TeamEc2OricUpdate"
        or ut_name == "ProjectEqConst2OricUpdate"
        or ut_name == "ProjectInequalityConstraintUpdate"
        or ut_name == "ExtrinsicOutputDepartment"
        or ut_name == "Team2KnapsackUpdate"
        or ut_name == "ProcessExtrinsicOutputDepartment"
    ):
        file_array_pure = open(f_path_pure, "rb")
        y_pure = np.load(file_array_pure)
        y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib)
        # y_wrapper = np.fromfile(f_path_wrapper, dtype=np.int32); #for integer y_wrapper
        check_agreement(
            sub_test_num,
            total_sub_tests,
            ut_name,
            y_pure,
            y_wrapper,
            rtol=rtol,
            atol=atol,
        )
