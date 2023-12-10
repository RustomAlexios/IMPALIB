# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + '/../src')
from ut_utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--sub_test_num', type=int,default=1,help="Sub-Test of Unit Test")
parser.add_argument('--sub_tests_total', type=int,default=1,help="Total Number of Sub-Tests of Unit Test")
parser.add_argument('--ut_name', type=str,help="Unit Test Name")
        
if __name__ == '__main__':
    args = parser.parse_args() 
    sub_test_num = args.sub_test_num 
    total_sub_tests = args.sub_tests_total 
    ut_name = args.ut_name
    rtol=5e-02; atol = 1e-04 #user will specify these tolerances
        
    if (ut_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate"):
        f_path_pure  = '../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_relaxed_graph_pure'
        f_path_wrapper  = '../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_relaxed_graph_wrapper'
    elif (ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate"):
        f_path_pure  = '../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_augmented_graph_pure'
        f_path_wrapper  = '../ut_results/ut_InputOutput/ut_'+ ut_name+'/extrinsic_output_edge_ec_augmented_graph_wrapper'
    elif (ut_name == "DegreeConstraint2EdgeEcUpdate"):
        f_path_pure  = '../ut_results/ut_DegreeConstraint/ut_'+ ut_name+'/degree_constraint_to_eq_constraint_m_pure'
        f_path_wrapper  = '../ut_results/ut_DegreeConstraint/ut_'+ ut_name+'/degree_constraint_to_eq_constraint_m_wrapper'
    elif (ut_name == "SubtourConstraints2EdgeEcUpdate"):
        f_path_pure  = '../ut_results/ut_SubtourConstraint/ut_'+ ut_name+'/subtour_constraints_to_edge_ec_m_pure'
        f_path_wrapper  = '../ut_results/ut_SubtourConstraint/ut_'+ ut_name+'/subtour_constraints_to_edge_ec_m_wrapper'
    elif (ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate"):
        f_path_pure  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_degree_constraint_m_pure'
        f_path_wrapper  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_degree_constraint_m_wrapper'
    elif (ut_name == "EdgeEc2SubtourConstraintsUpdate"):
        f_path_pure  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_subtour_constraints_m_pure'
        f_path_wrapper  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_subtour_constraints_m_wrapper'
    elif (ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate"):
        f_path_pure  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_degree_constraint_m_pure'
        f_path_wrapper  = '../ut_results/ut_EqualityConstraint/ut_'+ ut_name+'/edge_ec_to_degree_constraint_m_wrapper'
    elif (ut_name == "IterateRelaxedGraph" or ut_name == "IterateAugmentedGraph"):
        f_path_pure  = '../ut_results/ut_GraphicalModel/ut_'+ ut_name+'/intrinsic_out_edge_ec_pure'
        f_path_wrapper  = '../ut_results/ut_GraphicalModel/ut_'+ ut_name+'/intrinsic_out_edge_ec_wrapper'
    
    if (ut_name == "ExtrinsicOutputEdgeEcRelaxedGraphUpdate" or ut_name == "ExtrinsicOutputEdgeEcAugmentedGraphUpdate" or ut_name == "DegreeConstraint2EdgeEcUpdate" \
                    or ut_name == "SubtourConstraints2EdgeEcUpdate" or ut_name == "EdgeEc2DegreeConstraintRelaxedGraphUpdate" or ut_name == "EdgeEc2SubtourConstraintsUpdate" \
                        or ut_name == "EdgeEc2DegreeConstraintAugmentedGraphUpdate" or ut_name == "IterateRelaxedGraph" or ut_name == "IterateAugmentedGraph"):
        file_array_pure = open(f_path_pure, 'rb');
        y_pure = np.load(file_array_pure);
        y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib);
        #print("y_pure: ", y_pure)
        #print("y_wrapper: ", y_wrapper)
        #y_wrapper = np.fromfile(f_path_wrapper, dtype=np.int32); #for integer y_wrapper
        check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=rtol, atol = atol)  

