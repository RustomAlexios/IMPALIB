# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys
import input_output
from environmentModule import np, os, np_impa_lib
# from ut_utils import *

sys.path.append(sys.path[0] + "/../src")


def ut_io(ut_name, n_departments, n_teams, n_projects):
    N_DEPARTMENTS = n_departments
    N_TEAMS = n_teams
    N_PROJECTS = n_projects

    f_reward_project_path = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/reward_project_pure.npy"
    reward_project_pure = np.random.uniform(10, 500, size=(N_PROJECTS, N_TEAMS))
    reward_project_pure = reward_project_pure.astype(np_impa_lib)
    np.save(f_reward_project_path, reward_project_pure.flatten())

    outputs = input_output.OutputsKcMwm(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, reward_project_pure)

    f_input1 = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    if ut_name == "ExtrinsicOutputTeamUpdate":
        extrinsic_output_department_pure = np.random.uniform(10, 500, size=(N_DEPARTMENTS, N_TEAMS))
        extrinsic_output_department_pure = extrinsic_output_department_pure.astype(np_impa_lib)
        f_input2 = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/extrinsic_output_department_pure.npy"
        np.save(f_input2, extrinsic_output_department_pure.flatten())

        oric_to_team_m_pure = np.random.uniform(10, 500, N_TEAMS)
        oric_to_team_m_pure = oric_to_team_m_pure.astype(np_impa_lib)
        f_input3 = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/oric_to_team_m_pure.npy"
        np.save(f_input3, oric_to_team_m_pure.flatten())

        extrinsic_output_team_pure = outputs.extrinsic_output_team_update(
            extrinsic_output_department_pure, oric_to_team_m_pure
        )
        f_output_path = os.getcwd() + "/../ut_results/ut_InputOutput/ut_" + ut_name + "/extrinsic_output_team_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(
            output_file_python_pure,
            extrinsic_output_team_pure.flatten(),
            allow_pickle=True,
        )
        output_file_python_pure.close()

    elif ut_name == "IntrinsicOutMwmUpdate":
        oric_to_eq_constraint_m_pure = np.random.uniform(10, 500, size=(N_PROJECTS, N_TEAMS))
        oric_to_eq_constraint_m_pure = oric_to_eq_constraint_m_pure.astype(np_impa_lib)
        f_input2 = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/oric_to_eq_constraint_m_pure.npy"
        np.save(f_input2, oric_to_eq_constraint_m_pure.flatten())

        project_to_eq_constraint_m_pure = np.random.uniform(10, 500, size=(N_PROJECTS, N_TEAMS))
        project_to_eq_constraint_m_pure = project_to_eq_constraint_m_pure.astype(np_impa_lib)
        f_input3 = os.getcwd() + "/../ut_inputs/ut_InputOutput/ut_" + ut_name + "/project_to_eq_constraint_m_pure.npy"
        np.save(f_input3, project_to_eq_constraint_m_pure.flatten())

        intrinsic_out_mwm_pure = outputs.intrinsic_out_mwm_update(
            oric_to_eq_constraint_m_pure, project_to_eq_constraint_m_pure
        )
        f_output_path = os.getcwd() + "/../ut_results/ut_InputOutput/ut_" + ut_name + "/intrinsic_out_mwm_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, intrinsic_out_mwm_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
