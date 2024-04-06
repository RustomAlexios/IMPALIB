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
    int numTeams_;        ///< number of teams
    int numDepartments_;  ///< number of departments
    int numProjects_;     ///< number of projects
    int maxStateIc_ = 1;  ///< maximum state of inequality constraint (>=1)

   public:
    void messages_to_project_eq(const vector<vector<impalib_type>> &eq2oric, const vector<impalib_type> &team2oric, vector<vector<impalib_type>> &oric2eq, vector<vector<impalib_type>> &eq2proj,
                                const vector<vector<impalib_type>> &reward_project);  ///< update messages from ORIC to project equality constraint

    vector<impalib_type> messages_to_team_eq(const vector<vector<impalib_type>> &eq2oric);  ///< calculate messages from team ORIC to team equality constraint

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

inline OrInequalityConstraint::OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS) : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){};

/**
 * Calculate messages from ORIC to project equality constraints for the Knapsack-MWM problem
 *
 * @param[in] eq2oric: messages from project equality constraint to ORIC
 * @param[in] team2oric: messages from teams to ORIC
 * @param[out] oric2eq: messages from ORIC to project equality constraints
 * @param[in] eq2proj: messages from equality constraints to project inequality constraints
 * @param[in] reward_project: rewards for project equality constraints
 *
 */
inline void OrInequalityConstraint::messages_to_project_eq(const vector<vector<impalib_type>> &eq2oric, const vector<impalib_type> &team2oric, vector<vector<impalib_type>> &oric2eq,
                                                    vector<vector<impalib_type>> &eq2proj, const vector<vector<impalib_type>> &reward_project) {
    vector<vector<impalib_type>> forward(numProjects_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> backward(numProjects_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int i = 0; i < numTeams_; i++) {
        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value), initial_backward_messages(maxStateIc_ + 1, zero_value);

        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);
        forward[0] = initial_forward_messages;

        // Calculate forward messages
        for (int s = 0; s < numProjects_; s++) {
            forward[s + 1][0] = forward[s][0];
            forward[s + 1][1] = min(forward[s][1], forward[s][0] + eq2oric[s][i] + team2oric[i]);
        }

        // Set initial backward messages
        backward[numProjects_] = initial_backward_messages;

        // Calculate backward messages
        for (int s = numProjects_ - 1; s >= 0; s--) {
            backward[s][0] = min(backward[s + 1][0], backward[s + 1][1] + eq2oric[s][i] + team2oric[i]);
            backward[s][1] = backward[s + 1][1];
        }

        for (int j = 0; j < numProjects_; j++) {
            impalib_type minimumValue = zero_value;
            minimumValue = min(forward[j][1], backward[j + 1][0]);
            minimumValue = min(minimumValue, zero_value);
            oric2eq[j][i] = team2oric[i] - minimumValue;
            eq2proj[j][i] = oric2eq[j][i] + reward_project[j][i];
        }
    }
}

/**
 * Calculate messages from ORIC to team equality constraints for the Knapsack-MWM problem
 *
 * @param[in] eq2oric: messages from project equality constraint to ORIC
 * @returns : messages from ORIC to teams
 *
 */
inline vector<impalib_type> OrInequalityConstraint::messages_to_team_eq(const vector<vector<impalib_type>> &eq2oric) {
    vector<impalib_type> oric2team(numTeams_);
    for (int i = 0; i < oric2team.size(); i++) {
        impalib_type minValue = 1000000;

        for (int j = 0; j < eq2oric.size(); j++) {
            if (eq2oric[j][i] < minValue) {
                minValue = eq2oric[j][i];
            }
        }
        oric2team[i] = minValue;
    }
    return oric2team;
}