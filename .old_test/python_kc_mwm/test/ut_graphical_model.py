# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")

import graphical_model as model_graph
from ut_utils import prune_teams, team_reward_generation
from environmentModule import np, os, np_impa_lib, pkl


def ut_iterate(
    ut_name,
    N_u,
    team_types,
    filtering_flag,
    alpha,
    N_ITER,
    ppFlag,
    n_projects,
    threshold,
):
    np.random.seed(17)
    N_DEPARTMENTS = len(N_u)
    available_combinations = prune_teams(N_u)
    teams_weights_per_department, teams_types_per_department = team_reward_generation(
        available_combinations, N_u, team_types
    )
    teams_types = teams_types_per_department[0]

    N_TEAMS = len(teams_types)

    non_zero_weight_indices = []
    for department_index in range(0, N_DEPARTMENTS):
        indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
        non_zero_weight_indices = non_zero_weight_indices + [indices]

    non_zero_weight_indices_sizes = [len(indices) for indices in non_zero_weight_indices]
    max_size_nonzero_weights = max(non_zero_weight_indices_sizes)
    non_zero_weight_indices_arr = np.zeros((N_DEPARTMENTS, max_size_nonzero_weights), dtype=np.int32)

    f_input_alpha = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/alpha.npy"
    np.save(f_input_alpha, alpha)

    f_input1 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    f_input2 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/non_zero_weight_indices_sizes_pure.npy"
    np.save(f_input2, np.array(non_zero_weight_indices_sizes, dtype=np.int32).flatten())

    for department_index in range(0, N_DEPARTMENTS):
        teams_weights_per_department[department_index] = [
            np.int32(x) for x in teams_weights_per_department[department_index]
        ]  # input to C++

    f_input3 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/teams_weights_per_department_pure"
    np.save(f_input3, np.array(teams_weights_per_department, dtype=np.int32).flatten())

    # print('np.array(teams_weights_per_department, dtype=np.int32).flatten(): \n', np.array(teams_weights_per_department, dtype=np.int32).flatten())

    for i in range(len(non_zero_weight_indices)):
        for j in range(len(non_zero_weight_indices[i])):
            non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

    # print('non_zero_weight_indices_arr: ', non_zero_weight_indices_arr.flatten())

    f_input4 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/non_zero_weight_indices_arr_pure.npy"
    np.save(f_input4, non_zero_weight_indices_arr.flatten())
    reward_team = np.random.uniform(-300, 0, N_TEAMS)
    reward_team = reward_team.astype(np_impa_lib)
    reward_project = np.random.uniform(0, 300, size=(n_projects, N_TEAMS))
    reward_project = reward_project.astype(np_impa_lib)

    # print('reward_team: \n', reward_team)
    # print('reward_project: \n', reward_project)

    team_to_knapsack_m = np.zeros((N_DEPARTMENTS, N_TEAMS), dtype=np_impa_lib)
    for i in range(0, N_DEPARTMENTS):
        for j in range(0, N_TEAMS):
            if teams_weights_per_department[i][j] != 0:
                team_to_knapsack_m[i][j] = reward_team[j]

    # print('team_to_knapsack_m: \n', team_to_knapsack_m)

    f_input5 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/team_to_knapsack_m_pure.npy"
    np.save(f_input5, team_to_knapsack_m.flatten())
    input_load = []
    N_u = N_u.astype("int32")
    input_load.append(N_u)
    input_load.append(team_types)

    input_load.append(reward_team)
    input_load.append(reward_project.T)
    ModelIMPA = model_graph.GraphicalModelKcMwm(N_ITER, filtering_flag, ppFlag, alpha, THRESHOLD=threshold)
    ModelIMPA.initialize(input_load, test_flag=True)

    f_input6 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/reward_team_pure.npy"
    np.save(f_input6, reward_team.flatten())
    f_input7 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/reward_project_pure.npy"
    np.save(f_input7, reward_project.flatten())

    f_input8 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/max_state_py.npy"
    np.save(f_input8, N_u)

    ModelIMPA.iterate()

    f_output_path1 = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/extrinsic_output_team_pure"
    output_file_python1 = open(f_output_path1, "wb")
    np.save(
        output_file_python1,
        ModelIMPA.extrinsic_output_team.flatten(),
        allow_pickle=True,
    )
    output_file_python1.close()

    f_output_path2 = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_mwm_pure"
    output_file_python2 = open(f_output_path2, "wb")
    np.save(output_file_python2, ModelIMPA.intrinsic_out_mwm.flatten(), allow_pickle=True)
    output_file_python2.close()


def ut_iterate_sample_graph(ut_name, filtering_flag, alpha, N_ITER, ppFlag, threshold):
    setfile = np.random.randint(0, 2)
    # folder_inputs = '../../data/inputs_1000'
    folder_inputs = "../../data/inputs_random_params_1000"
    print("Graphical Model of Test Set: ", setfile)
    with open(str(folder_inputs) + "/inputs_set" + str(setfile) + ".pkl", "rb") as f:
        input_load = pkl.load(f)

    ModelIMPA = model_graph.GraphicalModelKcMwm(N_ITER, filtering_flag, ppFlag, alpha, THRESHOLD=threshold)
    ModelIMPA.initialize(input_load, test_flag=False)
    N_u = ModelIMPA.max_state
    N_DEPARTMENTS = len(N_u)
    teams_types = ModelIMPA.teams_types
    N_TEAMS = len(teams_types)

    teams_weights_per_department = ModelIMPA.teams_weights_per_department

    non_zero_weight_indices = []
    for department_index in range(0, N_DEPARTMENTS):
        indices = [i for i, e in enumerate(teams_weights_per_department[department_index]) if e != 0]
        non_zero_weight_indices = non_zero_weight_indices + [indices]

    non_zero_weight_indices_sizes = [len(indices) for indices in non_zero_weight_indices]
    max_size_nonzero_weights = max(non_zero_weight_indices_sizes)
    non_zero_weight_indices_arr = np.zeros((N_DEPARTMENTS, max_size_nonzero_weights), dtype=np.int32)

    f_input_projects = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/N_PROJECTS_pure.npy"
    np.save(f_input_projects, ModelIMPA.num_projects)

    f_input_departments = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/N_DEPARTMENTS_pure.npy"
    np.save(f_input_departments, N_DEPARTMENTS)

    f_input_alpha = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/alpha.npy"
    np.save(f_input_alpha, alpha)

    f_input1 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/N_TEAMS_pure.npy"
    np.save(f_input1, N_TEAMS)

    f_input2 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/non_zero_weight_indices_sizes_pure.npy"
    np.save(f_input2, np.array(non_zero_weight_indices_sizes, dtype=np.int32).flatten())

    for department_index in range(0, N_DEPARTMENTS):
        teams_weights_per_department[department_index] = [
            np.int32(x) for x in teams_weights_per_department[department_index]
        ]  # input to C++

    f_input3 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/teams_weights_per_department_pure"
    np.save(f_input3, np.array(teams_weights_per_department, dtype=np.int32).flatten())

    for i in range(len(non_zero_weight_indices)):
        for j in range(len(non_zero_weight_indices[i])):
            non_zero_weight_indices_arr[i, j] = non_zero_weight_indices[i][j]

    f_input4 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/non_zero_weight_indices_arr_pure.npy"
    np.save(f_input4, non_zero_weight_indices_arr.flatten())

    reward_team = ModelIMPA.reward_team
    reward_project = ModelIMPA.reward_project

    team_to_knapsack_m = ModelIMPA.team_to_knapsack_m

    f_input5 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/team_to_knapsack_m_pure.npy"
    np.save(f_input5, team_to_knapsack_m.flatten())

    f_input6 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/reward_team_pure.npy"
    np.save(f_input6, reward_team.flatten())

    f_input7 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/reward_project_pure.npy"
    np.save(f_input7, reward_project.flatten())

    N_u = N_u.astype("int32")
    f_input8 = os.getcwd() + "/../ut_inputs/ut_GraphicalModel/ut_" + ut_name + "/max_state_py.npy"
    np.save(f_input8, N_u)

    ModelIMPA.iterate()

    f_output_path1 = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/extrinsic_output_team_pure"
    output_file_python1 = open(f_output_path1, "wb")
    np.save(
        output_file_python1,
        ModelIMPA.extrinsic_output_team.flatten(),
        allow_pickle=True,
    )
    output_file_python1.close()

    f_output_path2 = os.getcwd() + "/../ut_results/ut_GraphicalModel/ut_" + ut_name + "/intrinsic_out_mwm_pure"
    output_file_python2 = open(f_output_path2, "wb")
    np.save(output_file_python2, ModelIMPA.intrinsic_out_mwm.flatten(), allow_pickle=True)
    output_file_python2.close()
