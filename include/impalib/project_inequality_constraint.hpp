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
class InequalityConstraint
{
private:
    int numProjects_; ///< number of projects
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int maxStateIc_ = 1; ///< maximum value of project inequality constraint (<=1)

public:
    void project_inequality_constraint_update(const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &) const; ///< calculate messages from project inequality constraint to project equality constraint
    InequalityConstraint(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS); ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * 
 */

inline InequalityConstraint::InequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                    };

/**
 * Calculate messages from project inequality constraint to project equality constraint for the Knapsack-MWM problem
 *
 * @param[in] rEqConstraint2ProjectM: messages from project equality constraint to project inequality constraint
 * @param[out] rProject2EqConstraintM: messages from project inequality constraint to project equality constraint
 * 
 */

inline void InequalityConstraint::project_inequality_constraint_update(const vector<vector<impalib_type>> &rEqConstraint2ProjectM,
                                                                vector<vector<impalib_type>> &rProject2EqConstraintM) const
{

    vector<vector<impalib_type>> stage_forward_messages_project_EC(numTeams_ + 1,
                                                                   vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages_project_EC(numTeams_ + 1,
                                                                    vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int project_index = 0; project_index < rProject2EqConstraintM.size(); project_index++)
    {
        // Initialize forward messages
        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value),
            initial_backward_messages(maxStateIc_ + 1, zero_value);
        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages_project_EC[0] = initial_forward_messages;
        
        for (int stage = 0; stage < numTeams_; stage++)
        {
            stage_forward_messages_project_EC[stage + 1][0] = stage_forward_messages_project_EC[stage][0];
            stage_forward_messages_project_EC[stage + 1][1] =
                min(stage_forward_messages_project_EC[stage][1],
                    stage_forward_messages_project_EC[stage][0] + rEqConstraint2ProjectM[project_index][stage]);
        }

        stage_backward_messages_project_EC[numTeams_] = initial_backward_messages;

        for (int stage = numTeams_ - 1; stage >= 0; stage--)
        {
            stage_backward_messages_project_EC[stage][0] =
                min(stage_backward_messages_project_EC[stage + 1][0],
                    stage_backward_messages_project_EC[stage + 1][1] + rEqConstraint2ProjectM[project_index][stage]);
            stage_backward_messages_project_EC[stage][1] = stage_backward_messages_project_EC[stage + 1][1];
        }

        // Update project to equality constraint messages
        for (int team_index = 0; team_index < numTeams_; team_index++)
        {
            impalib_type minimumValue                         = zero_value;
            minimumValue                                      = min(stage_forward_messages_project_EC[team_index][1],
                                                                    stage_backward_messages_project_EC[team_index + 1][0]);
            rProject2EqConstraintM[project_index][team_index] = -min(minimumValue, zero_value);
        }
    }
}