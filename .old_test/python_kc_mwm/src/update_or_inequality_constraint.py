# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np, np_impa_lib


class OrInequalityConstraint:
    def __init__(self, N_DEPARTMENTS, N_TEAMS, N_PROJECTS, unbalanced_flag, reward_project):
        self.num_departments = N_DEPARTMENTS
        self.num_teams = N_TEAMS
        self.num_projects = N_PROJECTS

        all_projects = np.array(range(0, self.num_projects))
        self.all_projects = all_projects
        self.unbalanced_flag = unbalanced_flag
        self.reward_project = reward_project

    def oric_to_project_ec_update(self, eq_constraint_to_oric_m, team_to_oric_m):
        num_projects = self.num_projects
        num_teams = self.num_teams

        oric_to_eq_constraint_m = np.zeros((num_projects, num_teams), dtype=np_impa_lib)
        eq_constraint_to_project_m = np.zeros((num_projects, num_teams), dtype=np_impa_lib)

        for project_index in range(0, num_projects):
            oric_to_eq_constraint_m[project_index, :] = self.oric_activation(
                project_index, eq_constraint_to_oric_m, team_to_oric_m
            )

        eq_constraint_to_project_m = oric_to_eq_constraint_m + self.reward_project
        self.oric_to_eq_constraint_m = oric_to_eq_constraint_m
        self.eq_constraint_to_project_m = eq_constraint_to_project_m
        return oric_to_eq_constraint_m, eq_constraint_to_project_m

    def oric_activation(self, project_index, eq_constraint_to_oric_m, team_to_oric_m):
        all_projects = self.all_projects

        remaining_projects = np.setdiff1d(all_projects, project_index)

        unbalanced_flag = self.unbalanced_flag
        or_input = team_to_oric_m
        ineq_constraint_input = eq_constraint_to_oric_m[remaining_projects, :]

        if unbalanced_flag:
            min_ic_input = ineq_constraint_input.min(axis=0)
            output_messages = -np.minimum(min_ic_input, -or_input)
        else:
            min_ic_input = ineq_constraint_input.min(axis=0)
            output_messages = -np.minimum(min_ic_input, -or_input)
        return output_messages

    def oric_to_team_update(self, eq_constraint_to_oric_m):
        oric_to_team_m = eq_constraint_to_oric_m.min(axis=0)
        self.oric_to_team_m = oric_to_team_m
        return oric_to_team_m
