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
class OrInequalityConstraint
{
private:
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int numProjects_; ///< number of projects
    int maxStateIc_ = 1; ///< maximum state of inequality constraint (>=1)

public:
    void oric_to_project_eq_constraint_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                              vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &); ///< update messages from ORIC to project equality constraint

    void oric_to_team_update(vector<vector<impalib_type>> &, vector<impalib_type> &); ///< calculate messages from team ORIC to team equality constraint

    OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS); ///< constructor
};

/**
 * Construct ORIC object for the Knapsack-MWM problem
 *
 * @param[in] N_DEPARTMENTS: number of departments
 * @param[in] N_TEAMS: number of teams
 * @param[in] N_PROJECTS: number of projects
 * @param[out] numProjects_: N_PROJECTS
 * @param[out] numTeams_: N_TEAMS
 * @param[out] numDepartments_: N_DEPARTMENTS
 * 
 */

OrInequalityConstraint::OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                    };

/**
 * Calculate messages from ORIC to project equality constraints for the Knapsack-MWM problem
 *
 * @param[in] rEqConstraint2OricM: messages from project equality constraint to ORIC
 * @param[in] mTeam2ORIC: messages from teams to ORIC
 * @param[out] rOric2EqConstraintM: messages from ORIC to project equality constraints
 * @param[in] rEqConstraint2ProjectM: messages from equality constraints to project inequality constraints
 * @param[in] rRewardProject: rewards for project equality constraints
 * 
 */

void OrInequalityConstraint::oric_to_project_eq_constraint_update(vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                                  vector<impalib_type>         &mTeam2ORIC,
                                                                  vector<vector<impalib_type>> &rOric2EqConstraintM,
                                                                  vector<vector<impalib_type>> &rEqConstraint2ProjectM,
                                                                  vector<vector<impalib_type>> &rRewardProject)
{
    // Initialize forward and backward messages
    vector<vector<impalib_type>> stage_forward_messages_ORIC_project(numProjects_ + 1,
                                                                     vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages_ORIC_project(
        numProjects_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    // Iterate over each team
    for (int team_index = 0; team_index < numTeams_; team_index++)
    {

        // Initialize forward and backward messages
        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value),
            initial_backward_messages(maxStateIc_ + 1, zero_value);
        
        // Set initial forward messages
        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages_ORIC_project[0] = initial_forward_messages;
        
        // Calculate forward messages
        for (int stage = 0; stage < numProjects_; stage++)
        {
            stage_forward_messages_ORIC_project[stage + 1][0] = stage_forward_messages_ORIC_project[stage][0];
            stage_forward_messages_ORIC_project[stage + 1][1] =
                min(stage_forward_messages_ORIC_project[stage][1], stage_forward_messages_ORIC_project[stage][0]
                                                                       + rEqConstraint2OricM[stage][team_index]
                                                                       + mTeam2ORIC[team_index]);
        }

        // Set initial backward messages
        stage_backward_messages_ORIC_project[numProjects_] = initial_backward_messages;

        // Calculate backward messages
        for (int stage = numProjects_ - 1; stage >= 0; stage--)
        {
            stage_backward_messages_ORIC_project[stage][0] =
                min(stage_backward_messages_ORIC_project[stage + 1][0],
                    stage_backward_messages_ORIC_project[stage + 1][1] + rEqConstraint2OricM[stage][team_index]
                        + mTeam2ORIC[team_index]);
            stage_backward_messages_ORIC_project[stage][1] = stage_backward_messages_ORIC_project[stage + 1][1];
        }

        // Update messages from ORIC to project equality constraints and from project equality constraints to project inequality constraints
        for (int project_index = 0; project_index < numProjects_; project_index++)
        {
            impalib_type minimumValue                      = zero_value;
            // Calculate the minimum value between forward and backward messages
            minimumValue                                   = min(stage_forward_messages_ORIC_project[project_index][1],
                                                                 stage_backward_messages_ORIC_project[project_index + 1][0]);
            // Calculate the minimum value between previous minimum and zero
            minimumValue                                   = min(minimumValue, zero_value);
            // Update messages from ORIC to project equality constraints
            rOric2EqConstraintM[project_index][team_index] = mTeam2ORIC[team_index] - minimumValue;
            // Update messages from project equality constraints to project inequality constraints
            rEqConstraint2ProjectM[project_index][team_index] =
                rOric2EqConstraintM[project_index][team_index] + rRewardProject[project_index][team_index];
        }
    }
}

/**
 * Calculate messages from ORIC to team equality constraints for the Knapsack-MWM problem
 *
 * @param[in] rEqConstraint2OricM: messages from project equality constraint to ORIC
 * @param[out] rOric2TeamM: messages from ORIC to teams
 * 
 */

void OrInequalityConstraint::oric_to_team_update(vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                 vector<impalib_type>         &rOric2TeamM)
{
    // Iterate over each team
    for (int team_index = 0; team_index < rOric2TeamM.size(); team_index++)
    {
        // Initialize minimum value
        impalib_type minValue = 1000000;

        // Find the minimum value among messages from project equality constraints to ORIC for this team
        for (int project_index = 0; project_index < rEqConstraint2OricM.size(); project_index++)
        {
            if (rEqConstraint2OricM[project_index][team_index] < minValue)
            {
                minValue = rEqConstraint2OricM[project_index][team_index];
            }
        }
        // Update message from ORIC to team with the minimum value
        rOric2TeamM[team_index] = minValue;
    }
}