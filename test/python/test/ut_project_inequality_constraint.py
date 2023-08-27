# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
sys.path.append(sys.path[0] + '/../src')

import project_inequality_constraint as ineq_constraint_project
from environmentModule import *

def ut_project_ineq_constraint(ut_name, n_departments, n_teams , n_projects, unbalanced_flag):
    N_DEPARTMENTS = n_departments
    N_TEAMS = n_teams
    N_PROJECTS = n_projects

    f_input1 = os.getcwd() + '/../ut_inputs/ut_projectIneqConstraint/ut_'+ ut_name+'/N_TEAMS_pure.npy'
    np.save(f_input1, N_TEAMS)
    
    project_ineq_constraint = ineq_constraint_project.InequalityConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, unbalanced_flag)

    f_input_path  = os.getcwd() + '/../ut_inputs/ut_projectIneqConstraint/ut_'+ ut_name+'/eq_constraint_to_project_m_pure.npy'
    eq_constraint_to_project_m_pure = np.random.uniform(-300, 300, size=(N_PROJECTS, N_TEAMS))
    eq_constraint_to_project_m_pure = eq_constraint_to_project_m_pure.astype(np_impa_lib)
    np.save(f_input_path, eq_constraint_to_project_m_pure.flatten())
        
    if (ut_name == "ProjectInequalityConstraintUpdate"):

        f_output_path  = os.getcwd() + '/../ut_results/ut_projectIneqConstraint/ut_'+ ut_name+'/project_to_eq_constraint_m_pure'
        project_to_eq_constraint_m_pure = project_ineq_constraint.project_inequality_constraint_update(eq_constraint_to_project_m_pure)
        output_file_python_pure = open(f_output_path, 'wb')
        np.save(output_file_python_pure, project_to_eq_constraint_m_pure, allow_pickle=True)
        output_file_python_pure.close()
    