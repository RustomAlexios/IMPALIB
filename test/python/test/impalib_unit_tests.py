# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ut_update_or_inequality_constraint 
import ut_update_equality_constraint
import ut_project_inequality_constraint
import ut_input_output
import ut_knapsack
import ut_graphical_model
from environmentModule import *


parser = argparse.ArgumentParser()
parser.add_argument('--sub_test_num', type=int,default=1,help="Sub-Test of Unit Test")
parser.add_argument('--sub_tests_total', type=int,default=1,help="Total Number of Sub-Tests of Unit Test")
parser.add_argument('--ut_name', type=str,help="Unit Test Name")
parser.add_argument('--nDepartments', type=int, help="Number of Departments")
parser.add_argument('--nTeams', type=int,help="Number of Teams")
parser.add_argument('--nProjects', type=int,help="Number of Projects")
parser.add_argument('--Nu', help='delimited Number of Units list', type=str)
parser.add_argument('--pTypes', help='Items Type', type=str)
parser.add_argument('--alpha', help='ALPHA', type=np_impa_lib)
parser.add_argument('--fFlag', help='FILTERING FLAG', type=bool)
parser.add_argument('--nIter', help='Number of Iterations', type=int)
parser.add_argument('--ppFlag', help='POST PROCESSING FLAG', type=bool)
parser.add_argument('--threshold', help='THRESHOLD', type=np_impa_lib)

if __name__ == '__main__':
    args = parser.parse_args() 
    sub_test_num = args.sub_test_num 
    total_sub_tests = args.sub_tests_total 
    ut_name = args.ut_name
    N_DEPARTMENTS = args.nDepartments
    N_TEAMS = args.nTeams
    N_PROJECTS = args.nProjects
    N_ITER = args.nIter
    if (N_TEAMS != N_PROJECTS):
        unbalanced_flag = True
    else:
        unbalanced_flag = False

    
    if (ut_name == "Oric2TeamUpdate" or ut_name == "Oric2ProjectEcUpdate"):
        ut_update_or_inequality_constraint.ut_oric(ut_name = ut_name, n_departments = N_DEPARTMENTS, n_teams = N_TEAMS, n_projects=N_PROJECTS, unbalanced_flag=unbalanced_flag);
    elif (ut_name ==  "ExtrinsicOutputTeamUpdate" or ut_name == "IntrinsicOutMwmUpdate"):
        ut_input_output.ut_io(ut_name = ut_name, n_departments = N_DEPARTMENTS, n_teams = N_TEAMS, n_projects=N_PROJECTS)
    elif (ut_name == "TeamEc2OricUpdate" or ut_name == "ProjectEqConst2OricUpdate"):
        ut_update_equality_constraint.ut_eq_consraint(ut_name = ut_name, n_departments = N_DEPARTMENTS, n_teams = N_TEAMS, n_projects=N_PROJECTS)
    elif (ut_name == "ProjectInequalityConstraintUpdate"):
        ut_project_inequality_constraint.ut_project_ineq_constraint(ut_name = ut_name, n_departments = N_DEPARTMENTS, n_teams = N_TEAMS, n_projects=N_PROJECTS, unbalanced_flag=unbalanced_flag);
    elif (ut_name == "KnapsackForward" or ut_name == "KnapsackBackward"):
        N_u = np.array([int(item) for item in args.Nu.split(',')]); team_types = [int(item) for item in args.pTypes.split(',')]; filtering_flag = args.fFlag; alpha = args.alpha
        ut_knapsack.ut_forward_backward(ut_name = ut_name, N_u = N_u, team_types = team_types, filtering_flag = filtering_flag,alpha=alpha)
    elif (ut_name == "ExtrinsicOutputDepartment"):
        N_u = np.array([int(item) for item in args.Nu.split(',')]); team_types = [int(item) for item in args.pTypes.split(',')]; filtering_flag = args.fFlag; alpha = args.alpha
        ut_knapsack.ut_extrinsic_output_department(ut_name = ut_name, N_u = N_u, team_types = team_types, filtering_flag = filtering_flag,alpha=alpha)
    elif (ut_name == "Team2KnapsackUpdate"):
        N_u = np.array([int(item) for item in args.Nu.split(',')]); team_types = [int(item) for item in args.pTypes.split(',')]; filtering_flag = args.fFlag; alpha = args.alpha
        ut_knapsack.ut_team_to_knapsack_update(ut_name = ut_name, N_u = N_u, team_types = team_types, filtering_flag = filtering_flag,alpha=alpha)
    elif (ut_name == "ProcessExtrinsicOutputDepartment"):
        N_u = np.array([int(item) for item in args.Nu.split(',')]); team_types = [int(item) for item in args.pTypes.split(',')]; filtering_flag = args.fFlag; alpha = args.alpha
        ut_knapsack.ut_process_extrinsic_output_department(ut_name = ut_name, N_u = N_u, team_types = team_types, filtering_flag = filtering_flag,alpha=alpha, N_ITER=N_ITER)
    elif (ut_name == "Iterate"):
        N_u = np.array([int(item) for item in args.Nu.split(',')]); team_types = [int(item) for item in args.pTypes.split(',')]; filtering_flag = args.fFlag; alpha = args.alpha
        ppFlag = args.ppFlag; threshold = args.threshold
        ut_graphical_model.ut_iterate(ut_name = ut_name, N_u = N_u, team_types = team_types, filtering_flag = filtering_flag,alpha=alpha, N_ITER=N_ITER, ppFlag = ppFlag, n_projects = N_PROJECTS, threshold = threshold)        
    elif (ut_name == "IterateSampleGraph"):
        filtering_flag = args.fFlag; alpha = args.alpha
        ppFlag = args.ppFlag; threshold = args.threshold
        ut_graphical_model.ut_iterate_sample_graph(ut_name = ut_name, filtering_flag = filtering_flag, alpha=alpha, N_ITER=N_ITER, ppFlag = ppFlag, threshold = threshold)      