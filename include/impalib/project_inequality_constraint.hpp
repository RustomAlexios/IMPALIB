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
 * @param[in] eq2ProjectM: messages from project equality constraint to project inequality constraint
 * @param[out] project2EqM: messages from project inequality constraint to project equality constraint
 * 
 */

inline void InequalityConstraint::project_inequality_constraint_update(const vector<vector<impalib_type>> &eq2ProjectM,
                                                                vector<vector<impalib_type>> &project2EqM) const
{

    vector<vector<impalib_type>> forward(numTeams_ + 1,
                                                                   vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> backward(numTeams_ + 1,
                                                                    vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int project_index = 0; project_index < project2EqM.size(); project_index++)
    {
        // Initialize forward messages
        vector<impalib_type> forward0(maxStateIc_ + 1, zero_value),
            backward0(maxStateIc_ + 1, zero_value);
        fill(forward0.begin() + 1, forward0.end(), value_inf);

        forward[0] = forward0;

        for (int stage = 0; stage < numTeams_; stage++)
        {
            forward[stage + 1][0] = forward[stage][0];
            forward[stage + 1][1] =
                min(forward[stage][1],
                    forward[stage][0] + eq2ProjectM[project_index][stage]);
        }

        backward[numTeams_] = backward0;

        for (int stage = numTeams_ - 1; stage >= 0; stage--)
        {
            backward[stage][0] =
                min(backward[stage + 1][0],
                    backward[stage + 1][1] + eq2ProjectM[project_index][stage]);
            backward[stage][1] = backward[stage + 1][1];
        }

        // Update project to equality constraint messages
        for (int team = 0; team < numTeams_; team++)
        {
            impalib_type minimumValue                         = zero_value;
            minimumValue                                      = min(forward[team][1],
                                                                    backward[team + 1][0]);
            project2EqM[project_index][team] = -min(minimumValue, zero_value);
        }
    }
}