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
    
    if (test_name == "ExtrinsicUpdate"){
         ut_input_output_mobarp(test_name);
    }

    if (test_name == "IneqCapacConstUpdate" || test_name =="MobileLocEqConst2REqConstUpdate" || test_name =="SetCoverIneqConstUpdate"){
        ut_ineq_const_mobarp_update(test_name);
    }

    if (test_name == "AuxiliaryConst2ZEqConstUpdate" || test_name == "AuxiliaryConst2REqConstUpdate" || test_name == "AuxiliaryConst2MobileXEqConstUpdate"){
        ut_auxiliary_const_update(test_name);
    }

    if (test_name == "XEqConst2AuxiliaryAndSetCoverConstUpdate" || test_name == "ReqConstActivation" || test_name == "ZEqConst2AuxiliaryConstUpdate" || test_name == "XEqConstActivation" || test_name == "ZEqConst2SetCoverIneqConstUpdate"){
        ut_equality_constraint_mobarp(test_name);
    }

    if (test_name == "Iterate"){
         ut_model_graph_mobarp(test_name);
    }
}