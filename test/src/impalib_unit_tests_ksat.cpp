// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "../include/impalib_unit_tests.hpp"

int main(){

    char *test_name_bash = getenv("test_name");
    if(test_name_bash == NULL){cout << "test_name_bash not available\n";}

    string test_name(test_name_bash);
    
    if (test_name == "ExtrinsicOutputVariableEcUpdate"){
        ut_input_output_ksat(test_name);
    }

    if (test_name == "VariableEc2KsatConstraintUpdate"){
        ut_equality_constraint_ksat(test_name);
    }

    if (test_name == "KsatConstraint2VariableEcUpdate"){
        ut_ksat_constraint(test_name);
    }

    if (test_name == "Iterate"){
        ut_model_graph_ksat(test_name);
    }


}