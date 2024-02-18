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

    vector<impalib_type> oric_to_team_update(vector<vector<impalib_type>> &); ///< calculate messages from team ORIC to team equality constraint

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
    vector<vector<impalib_type>> stage_forward_messages_ORIC_project(numProjects_ + 1,
                                                                     vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages_ORIC_project(
        numProjects_ + 1, vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int i = 0; i < numTeams_; i++)
    {

        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value),
            initial_backward_messages(maxStateIc_ + 1, zero_value);

        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);
        stage_forward_messages_ORIC_project[0] = initial_forward_messages;
        
        // Calculate forward messages
        for (int s = 0; s < numProjects_; s++)
        {
            stage_forward_messages_ORIC_project[s + 1][0] = stage_forward_messages_ORIC_project[s][0];
            stage_forward_messages_ORIC_project[s + 1][1] =
                min(stage_forward_messages_ORIC_project[s][1], stage_forward_messages_ORIC_project[s][0]
                                                                       + rEqConstraint2OricM[s][i]
                                                                       + mTeam2ORIC[i]);
        }

        // Set initial backward messages
        stage_backward_messages_ORIC_project[numProjects_] = initial_backward_messages;

        // Calculate backward messages
        for (int s = numProjects_ - 1; s >= 0; s--)
        {
            stage_backward_messages_ORIC_project[s][0] =
                min(stage_backward_messages_ORIC_project[s + 1][0],
                    stage_backward_messages_ORIC_project[s + 1][1] + rEqConstraint2OricM[s][i]
                        + mTeam2ORIC[i]);
            stage_backward_messages_ORIC_project[s][1] = stage_backward_messages_ORIC_project[s + 1][1];
        }

        for (int j = 0; j < numProjects_; j++)
        {
            impalib_type minimumValue                      = zero_value;
            minimumValue                                   = min(stage_forward_messages_ORIC_project[j][1],
                                                                 stage_backward_messages_ORIC_project[j + 1][0]);
            minimumValue                                   = min(minimumValue, zero_value);
            rOric2EqConstraintM[j][i] = mTeam2ORIC[i] - minimumValue;
            rEqConstraint2ProjectM[j][i] =
                rOric2EqConstraintM[j][i] + rRewardProject[j][i];
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

vector<impalib_type> OrInequalityConstraint::oric_to_team_update(vector<vector<impalib_type>> &rEqConstraint2OricM)
{
    vector<impalib_type> rOric2TeamM(numTeams_);
    for (int i = 0; i < rOric2TeamM.size(); i++)
    {
        impalib_type minValue = 1000000;

        for (int j = 0; j < rEqConstraint2OricM.size(); j++)
        {
            if (rEqConstraint2OricM[j][i] < minValue)
            {
                minValue = rEqConstraint2OricM[j][i];
            }
        }
        rOric2TeamM[i] = minValue;
    }
    return rOric2TeamM;
}