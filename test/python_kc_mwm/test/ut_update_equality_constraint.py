# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import update_equality_constraint as eq_constraint
import sys
from environmentModule import np, os, np_impa_lib

sys.path.append(sys.path[0] + "/../src")

def ut_eq_consraint(ut_name, n_departments, n_teams, n_projects):
    N_DEPARTMENTS = n_departments
    N_TEAMS = n_teams
    N_PROJECTS = n_projects

    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    reward_project_pure = np.random.uniform(-300, 300, size=(N_PROJECTS, N_TEAMS))
    reward_project_pure = reward_project_pure.astype(np_impa_lib)
    f_input2 = os.getcwd() + "/../ut_inputs/reward_project_pure.npy"
    np.save(f_input2, reward_project_pure.flatten())

    reward_team_pure = np.random.uniform(-300, 300, N_TEAMS)
    reward_team_pure = reward_team_pure.astype(np_impa_lib)
    f_reward_team_path = os.getcwd() + "/../ut_inputs/reward_team_pure.npy"
    np.save(f_reward_team_path, reward_team_pure.flatten())

    model_eq_constraint = eq_constraint.EqualityConstraintKcMwm(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, reward_team_pure, reward_project_pure)

    if ut_name == "TeamEc2OricUpdate":
        extrinsic_output_department_pure = np.random.uniform(-300, 300, size=(N_DEPARTMENTS, N_TEAMS))
        extrinsic_output_department_pure = extrinsic_output_department_pure.astype(np_impa_lib)
        f_input_path = os.getcwd() + "/../ut_inputs/extrinsic_output_department_pure"
        np.save(f_input_path, extrinsic_output_department_pure.flatten())

        f_output_path = os.getcwd() + "/../ut_results/team_to_oric_m_pure"
        team_to_oric_m_pure = model_eq_constraint.team_ec_to_oric_update(extrinsic_output_department_pure)
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, team_to_oric_m_pure, allow_pickle=True)
        output_file_python_pure.close()

    elif ut_name == "ProjectEqConst2OricUpdate":
        f_input_path = os.getcwd() + "/../ut_inputs/project_to_eq_constraint_m_pure.npy"
        project_to_eq_constraint_m_pure = np.random.uniform(-300, 300, size=(N_PROJECTS, N_TEAMS))
        project_to_eq_constraint_m_pure = project_to_eq_constraint_m_pure.astype(np_impa_lib)
        np.save(f_input_path, project_to_eq_constraint_m_pure.flatten())

        f_output_path = os.getcwd() + "/../ut_results/eq_constraint_to_oric_m_pure"
        eq_constraint_to_oric_m_pure = model_eq_constraint.project_eq_const_to_oric_update(project_to_eq_constraint_m_pure)
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, eq_constraint_to_oric_m_pure, allow_pickle=True)
        output_file_python_pure.close()
