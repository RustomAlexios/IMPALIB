// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

//#define IMPALIB_FLOAT_TYPE double

#include "impalib_unit_tests.hpp"

int main() {
    char *test_name_bash = getenv("test_name");
    if (test_name_bash == NULL) {
        cout << "test_name_bash not available\n";
    }

    string test_name(test_name_bash);

    if (test_name == "Oric2TeamUpdate" || test_name == "Oric2ProjectEcUpdate") {
        ut_oric(test_name);
    }
    if (test_name == "ExtrinsicOutputTeamUpdate" || test_name == "IntrinsicOutMwmUpdate") {
        ut_input_output_kc_mwm(test_name);
    }

    if (test_name == "TeamEc2OricUpdate" || test_name == "ProjectEqConst2OricUpdate") {
        ut_eq_constraint_kc_mwm(test_name);
    }

    if (test_name == "ProjectInequalityConstraintUpdate") {
        ut_project_ineq_constraint(test_name);
    }

    if (test_name == "KnapsackForward" || test_name == "KnapsackBackward") {
        ut_forward_backward(test_name);
    }

    if (test_name == "ExtrinsicOutputDepartment") {
        ut_extrinsic_output_department(test_name);
    }

    if (test_name == "Team2KnapsackUpdate") {
        ut_team_to_knapsack_update(test_name);
    }

    if (test_name == "ProcessExtrinsicOutputDepartment") {
        ut_process_extrinsic_output_department(test_name);
    }

    if (test_name == "Iterate") {
        ut_iterate_kc_mwm(test_name);
    }

    if (test_name == "IterateSampleGraph") {
        ut_iterate_sample_graph_kc_mwm(test_name);
    }
}