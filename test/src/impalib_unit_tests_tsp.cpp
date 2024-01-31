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
    
    if (test_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate" || test_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"){
        ut_input_output_tsp(test_name);
    }
    
    if (test_name == "DegreeConstraint2EdgeEcUpdate"){
        ut_degree_constraint(test_name);
    }

    if (test_name == "SubtourConstraints2EdgeEcUpdate"){
        ut_subtour_constraint(test_name);
    }

    if (test_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate" or test_name == "EdgeEc2SubtourConstraintsUpdate" or test_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate"){
        ut_equality_constraint_tsp(test_name);
    }

    if (test_name == "IterateRelaxedGraph" or test_name == "IterateAugmentedGraph"){
        ut_model_graph_tsp(test_name);
    }


}