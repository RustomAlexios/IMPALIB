# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)


import update_or_inequality_constraint as oric
import sys
from environmentModule import os, np, np_impa_lib

sys.path.append(sys.path[0] + "/../src")

def ut_oric(ut_name, n_departments, n_teams, n_projects, unbalanced_flag):
    N_DEPARTMENTS = n_departments
    N_TEAMS = n_teams
    N_PROJECTS = n_projects

    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    reward_project_pure = np.random.uniform(-300, 300, size=(N_PROJECTS, N_TEAMS))
    reward_project_pure = reward_project_pure.astype(np_impa_lib)
    f_input2 = os.getcwd() + "/../ut_inputs/reward_project_pure.npy"
    np.save(f_input2, reward_project_pure.flatten())

    model_oric = oric.OrInequalityConstraint(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, unbalanced_flag, reward_project_pure)

    eq_constraint_to_oric_m_pure = np.random.uniform(-300, 300, size=(N_PROJECTS, N_TEAMS))
    eq_constraint_to_oric_m_pure = eq_constraint_to_oric_m_pure.astype(np_impa_lib)
    f_input3 = os.getcwd() + "/../ut_inputs/eq_constraint_to_oric_m_pure.npy"
    np.save(f_input3, eq_constraint_to_oric_m_pure.flatten())

    if ut_name == "Oric2TeamUpdate":
        f_output_path = os.getcwd() + "/../ut_results/oric_to_team_m_pure"
        oric_to_team_m_pure = model_oric.oric_to_team_update(eq_constraint_to_oric_m_pure)
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, oric_to_team_m_pure, allow_pickle=True)
        output_file_python_pure.close()

    elif ut_name == "Oric2ProjectEcUpdate":
        team_to_oric_m_pure = np.random.uniform(-300, 300, N_TEAMS)
        team_to_oric_m_pure = team_to_oric_m_pure.astype(np_impa_lib)
        f_input = os.getcwd() + "/../ut_inputs/team_to_oric_m_pure.npy"
        np.save(f_input, team_to_oric_m_pure.flatten())

        (
            oric_to_eq_constraint_m_pure,
            eq_constraint_to_project_m_pure,
        ) = model_oric.oric_to_project_ec_update(eq_constraint_to_oric_m_pure, team_to_oric_m_pure)

        f_output_path1 = os.getcwd() + "/../ut_results/oric_to_eq_constraint_m_pure"
        output_file_python1 = open(f_output_path1, "wb")
        np.save(output_file_python1, oric_to_eq_constraint_m_pure, allow_pickle=True)
        output_file_python1.close()

        f_output_path2 = os.getcwd() + "/../ut_results/eq_constraint_to_project_m_pure"
        output_file_python2 = open(f_output_path2, "wb")
        np.save(output_file_python2, eq_constraint_to_project_m_pure, allow_pickle=True)
        output_file_python2.close()
