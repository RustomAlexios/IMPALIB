// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for the ORIC for Knapsack-MWM problem
 */
class OrInequalityConstraint {
   private:
    int nTeams_;          ///< number of teams
    int nDept_;           ///< number of departments
    int nProj_;           ///< number of projects
    int maxStateIc_ = 1;  ///< maximum state of inequality constraint (>=1)

   public:
    void oric_to_project_eq_constraint_update(const vector<vector<impalib_type>> &, const vector<impalib_type> &, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              const vector<vector<impalib_type>> &) const;  ///< update messages from ORIC to project equality constraint

    static void oric_to_team_update(const vector<vector<impalib_type>> &, vector<impalib_type> &);  ///< calculate messages from team ORIC to team equality constraint

    OrInequalityConstraint(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS);  ///< constructor
};

/**
 * Construct ORIC object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 *
 */

inline OrInequalityConstraint::OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : nProj_(N_PROJECTS), nTeams_(N_TEAMS), nDept_(N_DEPARTMENTS){};

/**
 * Calculate messages from ORIC to project equality constraints for the Knapsack-MWM problem
 *
 * @param[in] eq2OricM: messages from project equality constraint to ORIC
 * @param[in] team2OricM: messages from teams to ORIC
 * @param[out] oric2EqM: messages from ORIC to project equality constraints
 * @param[in] eq2ProjectM: messages from equality constraints to project inequality constraints
 * @param[in] rewards: rewards for project equality constraints
 *
 */

inline void OrInequalityConstraint::oric_to_project_eq_constraint_update(const vector<vector<impalib_type>> &eq2OricM, const vector<impalib_type> &team2OricM, vector<vector<impalib_type>> &oric2EqM,
                                                                         vector<vector<impalib_type>> &eq2ProjectM, const vector<vector<impalib_type>> &rewards) const {
    vector<vector<impalib_type>> forward(nProj_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> backward(nProj_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int team = 0; team < nTeams_; team++) {
        vector<impalib_type> forward0(maxStateIc_ + 1, zero_value);
        vector<impalib_type> backward0(maxStateIc_ + 1, zero_value);

        fill(forward0.begin() + 1, forward0.end(), value_inf);
        forward[0] = forward0;

        // Calculate forward messages
        for (int stage = 0; stage < nProj_; stage++) {
            forward[stage + 1][0] = forward[stage][0];
            forward[stage + 1][1] = min(forward[stage][1], forward[stage][0] + eq2OricM[stage][team] + team2OricM[team]);
        }

        // Set initial backward messages
        backward[nProj_] = backward0;

        // Calculate backward messages
        for (int stage = nProj_ - 1; stage >= 0; stage--) {
            backward[stage][0] = min(backward[stage + 1][0], backward[stage + 1][1] + eq2OricM[stage][team] + team2OricM[team]);
            backward[stage][1] = backward[stage + 1][1];
        }

        for (int project = 0; project < nProj_; project++) {
            impalib_type minimumValue = zero_value;
            minimumValue = min(forward[project][1], backward[project + 1][0]);
            minimumValue = min(minimumValue, zero_value);
            oric2EqM[project][team] = team2OricM[team] - minimumValue;
            eq2ProjectM[project][team] = oric2EqM[project][team] + rewards[project][team];
        }
    }
}

/**
 * Calculate messages from ORIC to team equality constraints for the Knapsack-MWM problem
 *
 * @param[in] eq2OricM: messages from project equality constraint to ORIC
 * @param[out] oric2TeamM: messages from ORIC to teams
 *
 */

inline void OrInequalityConstraint::oric_to_team_update(const vector<vector<impalib_type>> &eq2OricM, vector<impalib_type> &oric2TeamM) {
    for (int team = 0; team < oric2TeamM.size(); team++) {
        impalib_type minValue = 1000000;

        for (int project = 0; project < eq2OricM.size(); project++) {
            if (eq2OricM[project][team] < minValue) {
                minValue = eq2OricM[project][team];
            }
        }
        oric2TeamM[team] = minValue;
    }
}