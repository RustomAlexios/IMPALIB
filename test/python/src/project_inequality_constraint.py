# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *

class InequalityConstraint:
    def __init__(self, NUM_DEPARTMENTS, NUM_TEAMS, NUM_PROJECTS, unbalanced_flag):
        self.num_departments = NUM_DEPARTMENTS
        self.num_teams = NUM_TEAMS
        self.num_projects = NUM_PROJECTS
        self.unbalanced_flag = unbalanced_flag
        
    def project_inequality_constraint_update(self, A):
        return self.exactly_one_update(A, 0)

    def exactly_one_update(self, A, cols):
        
        B = np.zeros(A.shape, dtype = np_impa_lib)

        if cols:
            B[0] = self.constraint_rule(np.min(A[1:,:], axis  = 0))
            for j in range(1, A.shape[0] - 1):
                l = np.min(A[:j,:], axis  = 0)
                h = np.min(A[j+1:,:], axis  = 0)
                B[j] = self.constraint_rule(np.minimum(l, h))
            B[A.shape[0]-1] = self.constraint_rule(np.min(A[:A.shape[0] - 1,:], axis  = 0))
        else:
            B[:,0] = self.constraint_rule(np.min(A[:,1:], axis  = 1))
            for j in range(1, A.shape[1] - 1):
                l = np.min(A[:,:j], axis  = 1)
                h = np.min(A[:,j+1:], axis  = 1)
                B[:,j] = self.constraint_rule(np.minimum(l, h))
            B[:, A.shape[1]-1] = self.constraint_rule(np.min(A[:,:A.shape[1] - 1], axis  = 1))
        return -B #target_to_eq_constraint_m
        
    
    def constraint_rule(self, min_array):
        unbalanced_flag = self.unbalanced_flag
        if (not unbalanced_flag):
            min_message = min_array
        else:
            min_message = np.minimum(min_array, np.zeros(len(min_array), dtype = np_impa_lib))
        return min_message