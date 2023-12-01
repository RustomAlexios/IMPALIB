// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class EqualityConstraintKcMwm{
    private:
        int numTeams_;
        int numProjects_;
        int numDepartments_;

    public: 
        void team_eq_constraint_to_oric_update(vector<vector<impalib_type>>&, vector<impalib_type>&, vector<impalib_type>&);

        void project_eq_constraint_to_oric_update(vector<vector<impalib_type>>&, vector<vector<impalib_type>>&, vector<vector<impalib_type>>&);

    EqualityConstraintKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS);
};

EqualityConstraintKcMwm::EqualityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS){
    numProjects_ = N_PROJECTS;
    numTeams_ = N_TEAMS;
    numDepartments_ = N_DEPARTMENTS;

};

    void EqualityConstraintKcMwm::team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rExtrinsicOutputDepartment, 
                                vector<impalib_type> &rTeam2OricM, vector<impalib_type>& rewardTeam){

        vector <impalib_type> intermediate_team_to_oric_m(numTeams_,0);
        
        for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++){
            transform(rExtrinsicOutputDepartment[department_index].begin(), rExtrinsicOutputDepartment[department_index].end(), 
                                    intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
        }
        transform(intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.end(), 
                                    rewardTeam.begin(), rTeam2OricM.begin(), std::plus<impalib_type>());
        }

    void EqualityConstraintKcMwm::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rProject2EqConstraintM,
                                        vector<vector<impalib_type>> &rEqConstraint2OricM, vector<vector<impalib_type>>& rewardProject){
        for (int project_index=0; project_index < rProject2EqConstraintM.size(); project_index++){
            transform(rProject2EqConstraintM[project_index].begin(), rProject2EqConstraintM[project_index].end(), 
                                    rewardProject[project_index].begin(), rEqConstraint2OricM[project_index].begin(), std::plus<impalib_type>());
        }

        }