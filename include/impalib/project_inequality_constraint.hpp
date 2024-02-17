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
    void project_inequality_constraint_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &); ///< calculate messages from project inequality constraint to project equality constraint
    InequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS); ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numDepartments_: N_DEPARTMENTS
 * 
 */

InequalityConstraint::InequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                    };

/**
 * Calculate messages from project inequality constraint to project equality constraint for the Knapsack-MWM problem
 *
 * @param[in] rEqConstraint2ProjectM: messages from project equality constraint to project inequality constraint
 * @param[out] rProject2EqConstraintM: messages from project inequality constraint to project equality constraint
 * 
 */

void InequalityConstraint::project_inequality_constraint_update(vector<vector<impalib_type>> &rEqConstraint2ProjectM,
                                                                vector<vector<impalib_type>> &rProject2EqConstraintM)
{

    vector<vector<impalib_type>> stage_forward_messages_project_EC(numTeams_ + 1,
                                                                   vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages_project_EC(numTeams_ + 1,
                                                                    vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int i = 0; i < rProject2EqConstraintM.size(); i++)
    {
        // Initialize forward messages
        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value),
            initial_backward_messages(maxStateIc_ + 1, zero_value);
        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages_project_EC[0] = initial_forward_messages;
        
        for (int s = 0; s < numTeams_; s++)
        {
            stage_forward_messages_project_EC[s + 1][0] = stage_forward_messages_project_EC[s][0];
            stage_forward_messages_project_EC[s + 1][1] =
                min(stage_forward_messages_project_EC[s][1],
                    stage_forward_messages_project_EC[s][0] + rEqConstraint2ProjectM[i][s]);
        }

        stage_backward_messages_project_EC[numTeams_] = initial_backward_messages;

        for (int s = numTeams_ - 1; s >= 0; s--)
        {
            stage_backward_messages_project_EC[s][0] =
                min(stage_backward_messages_project_EC[s + 1][0],
                    stage_backward_messages_project_EC[s + 1][1] + rEqConstraint2ProjectM[i][s]);
            stage_backward_messages_project_EC[s][1] = stage_backward_messages_project_EC[s + 1][1];
        }

        // Update project to equality constraint messages
        for (int j = 0; j < numTeams_; j++)
        {
            impalib_type minimumValue                         = zero_value;
            minimumValue                                      = min(stage_forward_messages_project_EC[j][1],
                                                                    stage_backward_messages_project_EC[j + 1][0]);
            rProject2EqConstraintM[i][j] = -min(minimumValue, zero_value);
        }
    }
}