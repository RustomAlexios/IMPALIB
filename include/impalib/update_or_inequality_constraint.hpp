// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class OrInequalityConstraint{
    private:
        int numTeams_;
        int numDepartments_;
        int numProjects_;
        int maxStateIc_; 

    public:
        void oric_to_project_eq_constraint_update(vector<vector<impalib_type>>&, vector<impalib_type>&,
                        vector<vector<impalib_type>>&, vector<vector<impalib_type>>&, vector<vector<impalib_type>>&);
        
        void oric_to_team_update(vector<vector<impalib_type>>&, vector<impalib_type>&);
        

    OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS);
};

OrInequalityConstraint::OrInequalityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS){
    numProjects_ = N_PROJECTS;
    numTeams_ = N_TEAMS;
    numDepartments_ = N_DEPARTMENTS;
    maxStateIc_ = 1;

};

void OrInequalityConstraint::oric_to_project_eq_constraint_update(vector<vector<impalib_type>> &rEqConstraint2OricM, vector<impalib_type> &mTeam2ORIC,
                        vector<vector<impalib_type>> &rOric2EqConstraintM, vector<vector<impalib_type>> &rEqConstraint2ProjectM, vector<vector<impalib_type>>& rRewardProject){

        vector<vector<impalib_type>> stage_forward_messages_ORIC_project(numProjects_+1, vector<impalib_type>(maxStateIc_+1,zero_value));
        vector<vector<impalib_type>> stage_backward_messages_ORIC_project(numProjects_+1, vector<impalib_type>(maxStateIc_+1,zero_value));

        for (int team_index=0; team_index<numTeams_; team_index++){

            vector<impalib_type> initial_forward_messages(maxStateIc_+1, zero_value), initial_backward_messages(maxStateIc_+1, zero_value);
            fill(initial_forward_messages.begin()+1, initial_forward_messages.end(), value_inf);

            stage_forward_messages_ORIC_project[0] = initial_forward_messages;

            for (int stage = 0; stage <numProjects_; stage++){
                stage_forward_messages_ORIC_project[stage+1][0] = stage_forward_messages_ORIC_project[stage][0];
                stage_forward_messages_ORIC_project[stage+1][1] = min(stage_forward_messages_ORIC_project[stage][1], stage_forward_messages_ORIC_project[stage][0]+rEqConstraint2OricM[stage][team_index]+mTeam2ORIC[team_index]);
            }

            stage_backward_messages_ORIC_project[numProjects_] = initial_backward_messages;

            for (int stage = numProjects_-1; stage >=0; stage--){
                stage_backward_messages_ORIC_project[stage][0] = min(stage_backward_messages_ORIC_project[stage+1][0], stage_backward_messages_ORIC_project[stage+1][1]+rEqConstraint2OricM[stage][team_index]+mTeam2ORIC[team_index]);
                stage_backward_messages_ORIC_project[stage][1] = stage_backward_messages_ORIC_project[stage+1][1];
            }

            for (int project_index = 0; project_index <numProjects_; project_index++){
                impalib_type minimumValue;
                minimumValue = min(stage_forward_messages_ORIC_project[project_index][1],stage_backward_messages_ORIC_project[project_index+1][0]);
                minimumValue = min(minimumValue,zero_value);
                rOric2EqConstraintM[project_index][team_index] = mTeam2ORIC[team_index]-minimumValue;
                rEqConstraint2ProjectM[project_index][team_index] = rOric2EqConstraintM[project_index][team_index] + rRewardProject[project_index][team_index];
            }
        }


        }

void OrInequalityConstraint::oric_to_team_update(vector<vector<impalib_type>> &rEqConstraint2OricM, vector<impalib_type> &rOric2TeamM){

            for (int team_index = 0; team_index < rOric2TeamM.size(); team_index++){
                impalib_type minValue = 1000000;
                for (int project_index=0; project_index < rEqConstraint2OricM.size(); project_index++){
                    if (rEqConstraint2OricM[project_index][team_index]<minValue){minValue = rEqConstraint2OricM[project_index][team_index];}
                }
                rOric2TeamM[team_index] = minValue;
                }
        }