# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *

class EqualityConstraintKcMwm:
    def __init__(self, N_DEPARTMENTS, N_TEAMS, N_PROJECTS, reward_team, reward_project):
        self.num_bins = N_DEPARTMENTS
        self.num_items = N_TEAMS
        self.num_targets = N_PROJECTS
        self.reward_team = reward_team
        self.reward_project = reward_project
        
    def team_ec_to_oric_update(self, extrinsic_output_department):
        
        team_to_oric_m = np.sum(extrinsic_output_department, axis=0) + self.reward_team
        
        return team_to_oric_m
    
    def project_eq_const_to_oric_update(self, project_to_eq_constraint_m):
        eq_constraint_to_oric_m = project_to_eq_constraint_m + self.reward_project
        return eq_constraint_to_oric_m