# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")
import knapsack as knapsack_methods
from ut_utils import team_reward_generation, prune_teams
from environmentModule import np, os, np_impa_lib


def ut_forward_backward(ut_name, N_u, team_types, filtering_flag, alpha):
    N_DEPARTMENTS = N_u.size
    available_combinations = prune_teams(N_u)
    teams_weights_per_department, teams_types_per_department = team_reward_generation(available_combinations, N_u, team_types)
    teams_types = teams_types_per_department[0]

    teams_weights_per_department = np.array(teams_weights_per_department, dtype=np.int32)  # added to make it array with type int32

    N_TEAMS = len(teams_types)

    non_zero_weight_indices = []
    for department_index in range(0, N_DEPARTMENTS):
        indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
        non_zero_weight_indices = non_zero_weight_indices + [indices]

    non_zero_weight_indices_sizes = [len(indices) for indices in non_zero_weight_indices]  # input to C++
    max_size_nonzero_weights = max(non_zero_weight_indices_sizes)
    non_zero_weight_indices_arr = np.zeros((N_DEPARTMENTS, max_size_nonzero_weights), dtype=np.int32)  # input to C++ # added np.int32

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, alpha)

    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    f_input2 = os.getcwd() + "/../ut_inputs/non_zero_weight_indices_sizes_pure.npy"
    np.save(f_input2, np.array(non_zero_weight_indices_sizes, dtype=np.int32).flatten())

    f_input3 = os.getcwd() + "/../ut_inputs/teams_weights_per_department_pure"
    np.save(f_input3, teams_weights_per_department.flatten())

    for i in range(len(non_zero_weight_indices)):
        for j in range(len(non_zero_weight_indices[i])):
            non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

    f_input4 = os.getcwd() + "/../ut_inputs/non_zero_weight_indices_arr_pure.npy"
    np.save(f_input4, non_zero_weight_indices_arr.flatten())

    team_to_knapsack_m = np.zeros((N_DEPARTMENTS, N_TEAMS), dtype=np_impa_lib)
    reward_team = -np.random.uniform(10, 1000, N_TEAMS)
    reward_team = reward_team.astype(np_impa_lib)

    for i in range(0, N_DEPARTMENTS):
        for j in range(0, N_TEAMS):
            if teams_weights_per_department[i][j] != 0:
                team_to_knapsack_m[i][j] = reward_team[j]

    f_input5 = os.getcwd() + "/../ut_inputs/team_to_knapsack_m_pure.npy"
    np.save(f_input5, team_to_knapsack_m.flatten())

    model_knapsacks = knapsack_methods.Knapsack(N_DEPARTMENTS, N_TEAMS, filtering_flag, alpha, reward_team)

    for department_index in range(0, N_DEPARTMENTS):
        if ut_name == "KnapsackForward":
            model_knapsacks.forward(
                N_u[department_index],
                team_to_knapsack_m[department_index],
                teams_weights_per_department[department_index],
            )
            f_pure_path = os.getcwd() + "/../ut_results/forward_pure" + str(department_index)
            output_file_python_pure = open(f_pure_path, "wb")
            np.save(
                output_file_python_pure,
                model_knapsacks.stage_forward_messages,
                allow_pickle=True,
            )
            output_file_python_pure.close()

        if ut_name == "KnapsackBackward":
            model_knapsacks.backward(
                N_u[department_index],
                team_to_knapsack_m[department_index],
                teams_weights_per_department[department_index],
            )
            f_pure_path = os.getcwd() + "/../ut_results/backward_pure" + str(department_index)
            output_file_python_pure = open(f_pure_path, "wb")
            np.save(
                output_file_python_pure,
                model_knapsacks.stage_backward_messages,
                allow_pickle=True,
            )
            output_file_python_pure.close()


def ut_extrinsic_output_department(ut_name, N_u, team_types, filtering_flag, alpha):
    N_DEPARTMENTS = N_u.size
    available_combinations = prune_teams(N_u)
    teams_weights_per_department, teams_types_per_department = team_reward_generation(available_combinations, N_u, team_types)
    teams_types = teams_types_per_department[0]

    teams_weights_per_department = np.array(teams_weights_per_department, dtype=np.int32)  # added to make it array with type int32

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, alpha)

    N_TEAMS = len(teams_types)
    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    f_input2 = os.getcwd() + "/../ut_inputs/teams_weights_per_department_pure"
    np.save(f_input2, teams_weights_per_department.flatten())

    team_to_knapsack_m = np.zeros((N_DEPARTMENTS, N_TEAMS), dtype=np_impa_lib)
    reward_team = -np.random.uniform(10, 1000, N_TEAMS)
    reward_team = reward_team.astype(np_impa_lib)

    for i in range(0, N_DEPARTMENTS):
        for j in range(0, N_TEAMS):
            if teams_weights_per_department[i][j] != 0:
                team_to_knapsack_m[i][j] = reward_team[j]

    f_input3 = os.getcwd() + "/../ut_inputs/team_to_knapsack_m_pure.npy"
    np.save(f_input3, team_to_knapsack_m.flatten())

    model_knapsacks = knapsack_methods.Knapsack(N_DEPARTMENTS, N_TEAMS, filtering_flag, alpha, reward_team)

    for department_index in range(0, N_DEPARTMENTS):
        stage_forward_messages = np.random.uniform(low=-300, high=300, size=(N_TEAMS + 1, N_u[department_index] + 1))
        stage_forward_messages = stage_forward_messages.astype(np_impa_lib)
        stage_backward_messages = np.random.uniform(low=-300, high=300, size=(N_TEAMS + 1, N_u[department_index] + 1))
        stage_backward_messages = stage_backward_messages.astype(np_impa_lib)

        f_input_fv = os.getcwd() + "/../ut_inputs/stage_forward_messages" + str(department_index) + ".npy"
        np.save(f_input_fv, stage_forward_messages.flatten())
        f_input_bv = os.getcwd() + "/../ut_inputs/stage_backward_messages" + str(department_index) + ".npy"
        np.save(f_input_bv, stage_backward_messages.flatten())

        model_knapsacks.stage_forward_messages = stage_forward_messages
        model_knapsacks.stage_backward_messages = stage_backward_messages
        model_knapsacks.extrinsic_output_department_lhs(
            department_index,
            N_u[department_index],
            team_to_knapsack_m[department_index],
            teams_weights_per_department[department_index],
        )

    f_pure_path = os.getcwd() + "/../ut_results/extrinsic_output_department_pure"
    output_file_python_pure = open(f_pure_path, "wb")
    np.save(
        output_file_python_pure,
        model_knapsacks.extrinsic_output_department_dummy.flatten(),
        allow_pickle=True,
    )
    output_file_python_pure.close()


def ut_team_to_knapsack_update(ut_name, N_u, team_types, filtering_flag, alpha):
    N_DEPARTMENTS = N_u.size
    available_combinations = prune_teams(N_u)
    teams_weights_per_department, teams_types_per_department = team_reward_generation(available_combinations, N_u, team_types)
    teams_types = teams_types_per_department[0]

    teams_weights_per_department = np.array(teams_weights_per_department, dtype=np.int32)  # added to make it array with type int32

    N_TEAMS = len(teams_types)

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, alpha)

    non_zero_weight_indices = []
    for department_index in range(0, N_DEPARTMENTS):
        indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
        non_zero_weight_indices = non_zero_weight_indices + [indices]

    non_zero_weight_indices_sizes = [len(indices) for indices in non_zero_weight_indices]  # input to C++
    max_size_nonzero_weights = max(non_zero_weight_indices_sizes)
    non_zero_weight_indices_arr = np.zeros((N_DEPARTMENTS, max_size_nonzero_weights), dtype=np.int32)  # input to C++ # added np.int32

    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    f_input2 = os.getcwd() + "/../ut_inputs/non_zero_weight_indices_sizes_pure.npy"
    np.save(f_input2, np.array(non_zero_weight_indices_sizes, dtype=np.int32).flatten())

    for i in range(len(non_zero_weight_indices)):
        for j in range(len(non_zero_weight_indices[i])):
            non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

    f_input3 = os.getcwd() + "/../ut_inputs/non_zero_weight_indices_arr_pure.npy"
    np.save(f_input3, non_zero_weight_indices_arr.flatten())

    team_to_knapsack_m = np.zeros((N_DEPARTMENTS, N_TEAMS), dtype=np_impa_lib)
    reward_team = -np.random.uniform(10, 1000, N_TEAMS)
    reward_team = reward_team.astype(np_impa_lib)

    team_to_knapsack_m = np.random.uniform(10, 1000, size=(N_DEPARTMENTS, N_TEAMS))
    team_to_knapsack_m = team_to_knapsack_m.astype(np_impa_lib)

    f_input4 = os.getcwd() + "/../ut_inputs/team_to_knapsack_m_pure.npy"
    np.save(f_input4, team_to_knapsack_m.flatten())

    f_input5 = os.getcwd() + "/../ut_inputs/reward_team_pure.npy"
    np.save(f_input5, reward_team.flatten())

    extrinsic_output_department = np.random.uniform(-300, 300, size=(N_DEPARTMENTS, N_TEAMS))
    extrinsic_output_department = extrinsic_output_department.astype(np_impa_lib)
    f_input6 = os.getcwd() + "/../ut_inputs/extrinsic_output_department_pure.npy"
    np.save(f_input6, extrinsic_output_department.flatten())

    oric_to_team_m = np.random.uniform(-300, 300, N_TEAMS)
    oric_to_team_m = oric_to_team_m.astype(np_impa_lib)
    f_input7 = os.getcwd() + "/../ut_inputs/oric_to_team_m_pure.npy"
    np.save(f_input7, oric_to_team_m.flatten())

    model_knapsacks = knapsack_methods.Knapsack(N_DEPARTMENTS, N_TEAMS, filtering_flag, alpha, reward_team)
    model_knapsacks.extrinsic_output_department = extrinsic_output_department

    team_to_knapsack_m = model_knapsacks.team_to_department_update(teams_weights_per_department, oric_to_team_m)

    f_pure_path = os.getcwd() + "/../ut_results/team_to_knapsack_m_pure"
    output_file_python_pure = open(f_pure_path, "wb")
    np.save(output_file_python_pure, team_to_knapsack_m, allow_pickle=True)
    output_file_python_pure.close()


def ut_process_extrinsic_output_department(ut_name, N_u, team_types, filtering_flag, alpha, N_ITER):
    N_DEPARTMENTS = N_u.size
    available_combinations = prune_teams(N_u)
    teams_weights_per_department, teams_types_per_department = team_reward_generation(available_combinations, N_u, team_types)
    teams_types = teams_types_per_department[0]

    N_TEAMS = len(teams_types)

    f_input_alpha = os.getcwd() + "/../ut_inputs/alpha.npy"
    np.save(f_input_alpha, alpha)

    f_input1 = os.getcwd() + "/../ut_inputs/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    reward_team = -np.random.uniform(10, 1000, N_TEAMS)
    reward_team = reward_team.astype(np_impa_lib)

    extrinsic_output_department_dummy = np.random.uniform(-300, 300, size=(N_DEPARTMENTS, N_TEAMS))
    extrinsic_output_department_dummy = extrinsic_output_department_dummy.astype(np_impa_lib)
    f_input2 = os.getcwd() + "/../ut_inputs/extrinsic_output_department_dummy_pure.npy"
    np.save(f_input2, extrinsic_output_department_dummy.flatten())

    model_knapsacks = knapsack_methods.Knapsack(N_DEPARTMENTS, N_TEAMS, filtering_flag, alpha, reward_team)
    model_knapsacks.extrinsic_output_department_dummy = extrinsic_output_department_dummy
    for iter in range(0, N_ITER):
        model_knapsacks.process_extrinsic_output_department(iter)

    f_pure_path = os.getcwd() + "/../ut_results/extrinsic_output_department_pure"
    output_file_python_pure = open(f_pure_path, "wb")
    np.save(
        output_file_python_pure,
        model_knapsacks.extrinsic_output_department,
        allow_pickle=True,
    )
    output_file_python_pure.close()
