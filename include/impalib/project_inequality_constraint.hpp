// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "common.hpp"

/**
 * Represents a class for the inequality constraint for the Knapsack-MWM problem
 */
class InequalityConstraint {
   public:
    InequalityConstraint(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);  ///< constructor
    vector<vector<impalib_type>> project_inequality_constraint_update(
        const vector<vector<impalib_type>> &) const;                       ///< calculate messages from project inequality constraint to project equality constraint

   private:
    int nProj_;           ///< number of projects
    int nTeams_;          ///< number of teams
    int nDept_;           ///< number of departments
    int maxStateIc_ = 1;  ///< maximum value of project inequality constraint (<=1)
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 *
 */

inline InequalityConstraint::InequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : nProj_(N_PROJECTS), nTeams_(N_TEAMS), nDept_(N_DEPARTMENTS){};

/**
 * Calculate messages from project inequality constraint to project equality constraint for the Knapsack-MWM problem
 *
 * @param[in] eq2ProjectM: messages from project equality constraint to project inequality constraint
 * @returns: messages from project inequality constraint to project equality constraint
 *
 */

inline vector<vector<impalib_type>> InequalityConstraint::project_inequality_constraint_update(const vector<vector<impalib_type>> &eq2ProjectM) const {
    auto project2EqM = eq2ProjectM;
    vector<vector<impalib_type>> forward(nTeams_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> backward(nTeams_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int project = 0; project < project2EqM.size(); project++) {
        // Initialize forward messages
        vector<impalib_type> forward0(maxStateIc_ + 1, zero_value);
        fill(forward0.begin() + 1, forward0.end(), value_inf);

        forward[0] = forward0;

        for (int stage = 0; stage < nTeams_; stage++) {
            forward[stage + 1][0] = forward[stage][0];
            forward[stage + 1][1] = min(forward[stage][1], forward[stage][0] + eq2ProjectM[project][stage]);
        }

        vector<impalib_type> backward0(maxStateIc_ + 1, zero_value);
        backward[nTeams_] = backward0;

        for (int stage = nTeams_ - 1; stage >= 0; stage--) {
            backward[stage][0] = min(backward[stage + 1][0], backward[stage + 1][1] + eq2ProjectM[project][stage]);
            backward[stage][1] = backward[stage + 1][1];
        }

        // Update project to equality constraint messages
        for (int team = 0; team < nTeams_; team++) {
            impalib_type minimumValue = min(forward[team][1], backward[team + 1][0]);
            project2EqM[project][team] = -min(minimumValue, zero_value);
        }
    }

    return project2EqM;
}