# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np

class OutputsKcMwm:
    def __init__(self, N_DEPARTMENTS, N_TEAMS, N_PROJECTS, reward_project):
        self.num_departments = N_DEPARTMENTS
        self.num_teams = N_TEAMS
        self.num_projects = N_PROJECTS
        self.reward_project = reward_project

    def intrinsic_out_mwm_update(self, oric_to_eq_constraint_m, project_to_eq_constraint_m):
        intrinsic_out_mwm = oric_to_eq_constraint_m + project_to_eq_constraint_m + self.reward_project
        return intrinsic_out_mwm

    def extrinsic_output_team_update(self, extrinsic_output_department, oric_to_team_m):
        extrinsic_output_team = np.sum(extrinsic_output_department, axis=0) + oric_to_team_m
        return extrinsic_output_team
