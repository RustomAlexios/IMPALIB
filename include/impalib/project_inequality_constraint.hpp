// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for the inequality constraint for the Knapsack-MWM problem
 */
class InequalityConstraint {
   private:
    int numProjects_;     ///< number of projects
    int numTeams_;        ///< number of teams
    int numDepartments_;  ///< number of departments
    int maxStateIc_ = 1;  ///< maximum value of project inequality constraint (<=1)

   public:
    vector<vector<impalib_type>> messages_to_equality(const vector<vector<impalib_type>> &eq2proj);  ///< calculate messages from project inequality constraint to project equality constraint
    InequalityConstraint(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);                            ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numDepartments_: N_DEPARTMENTS
 *
 */

InequalityConstraint::InequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){};

/**
 * Calculate messages from project inequality constraint to project equality constraint for the Knapsack-MWM problem
 *
 * @param[in] eq2proj: messages from project equality constraint to project inequality constraint
 * @param[out] rProject2EqConstraintM: messages from project inequality constraint to project equality constraint
 *
 */

vector<vector<impalib_type>> InequalityConstraint::messages_to_equality(const vector<vector<impalib_type>> &eq2proj) {
    vector<vector<impalib_type>> proj2eq(numProjects_, vector<impalib_type>(numTeams_, 0));
    vector<vector<impalib_type>> forward(numTeams_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> backward(numTeams_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int i = 0; i < proj2eq.size(); i++) {
        // Initialize forward messages
        vector<impalib_type> forward0(maxStateIc_ + 1, zero_value), backward0(maxStateIc_ + 1, zero_value);
        fill(forward0.begin() + 1, forward0.end(), value_inf);

        forward[0] = forward0;

        for (int s = 0; s < numTeams_; s++) {
            forward[s + 1][0] = forward[s][0];
            forward[s + 1][1] = min(forward[s][1], forward[s][0] + eq2proj[i][s]);
        }

        backward[numTeams_] = backward0;

        for (int s = numTeams_ - 1; s >= 0; s--) {
            backward[s][0] = min(backward[s + 1][0], backward[s + 1][1] + eq2proj[i][s]);
            backward[s][1] = backward[s + 1][1];
        }

        // Update project to equality constraint messages
        for (int j = 0; j < numTeams_; j++) {
            impalib_type minimumValue = zero_value;
            minimumValue = min(forward[j][1], backward[j + 1][0]);
            proj2eq[i][j] = -min(minimumValue, zero_value);
        }
    }
    return proj2eq;
}