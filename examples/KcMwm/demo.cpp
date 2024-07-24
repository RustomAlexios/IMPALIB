// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"

int main() {
    // Problem parameters
    const int N_PROJECTS = 2;                                          ///< number of projects
    const bool FILT_FLAG = true;                                       ///< whether filtering is activated or not
    const int N_ITER = 400;                                            ///< number of iterations of IMPA
    const int N_DEPARTMENTS = 2;                                       ///< number of departments
    array<int, N_DEPARTMENTS> max_state = {3, 3};                      ///< size of each department
    const impalib_type ALPHA = 0.9;                                    ///< filtering parameter
    const int N_TEAMS = 5;                                             ///< number of teams
    array<int, N_DEPARTMENTS> non_zero_weight_indices_sizes = {4, 4};  ///< sizes of non-zero connections per department

    // Extracting data pointers
    const int *pNON_ZERO_WEIGHT_INDICES_SIZES = non_zero_weight_indices_sizes.data();

    // Calculate max size of non-zero weight
    const auto max_size_non_zero_weight_iter = max_element(non_zero_weight_indices_sizes.begin(), non_zero_weight_indices_sizes.end());
    int max_size_non_zero_weight = *max_size_non_zero_weight_iter;

    // Initialize the graphical model
    GraphicalModelKcMwm model_graph(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, max_size_non_zero_weight, N_ITER, FILT_FLAG, ALPHA);

    // Define rewards for teams and projects
    array<impalib_type, N_TEAMS> reward_team = {-22, -31, -68, -39, -84};
    array<impalib_type, N_PROJECTS *N_TEAMS> reward_project = {44, 1, 41, 10, 3, 7, 17, 56, 98, 63};

    // Extracting data pointers
    const impalib_type *pREWARD_PROJECT = reward_project.data();
    const impalib_type *pREWARD_TEAM = reward_team.data();

    // Defining and extracting variables
    array<impalib_type, N_DEPARTMENTS *N_TEAMS> transition_model = {-22, -31, -68, -39, 0, 0, -31, -68, -39, -84};
    impalib_type *pTransition_model = transition_model.data();

    array<int, N_DEPARTMENTS *N_TEAMS> teams_weights_per_department = {2, 1, 1, 1, 0, 0, 1, 1, 1, 2};

    const int *pTEAMS_WEIGHTS_PER_DEPARTMENT = teams_weights_per_department.data();

    array<int, N_DEPARTMENTS *N_TEAMS> non_zero_weight_indices = {0, 1, 2, 3, 1, 2, 3, 4};
    const int *p_NON_ZERO_WEIGHT_INDICES = non_zero_weight_indices.data();

    const int *pMAX_STATE = max_state.data();

    // Initialize and iterate through the model
    model_graph.initialize(pREWARD_TEAM, pTransition_model, pTEAMS_WEIGHTS_PER_DEPARTMENT, pNON_ZERO_WEIGHT_INDICES_SIZES, p_NON_ZERO_WEIGHT_INDICES, pREWARD_PROJECT, pMAX_STATE);
    model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES);

    // Print results
    for (int i = 0; i < N_TEAMS; i++) {
        cout << model_graph.outputs.ExtrinsicOutputTeam[i] + model_graph.modelInputs_.RewardTeam[i] << '\n';
    }
}