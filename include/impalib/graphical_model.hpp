// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib.hpp"
class GraphicalModel{
        private:
                int numProjects_;
                int numTeams_;
                int numDepartments_;
                int maxSizeNonZeroWeights_;
                int numIterations_;
                bool filteringFlag_;
                impalib_type alpha_;
                vector<vector<impalib_type>> extrinsicOutputDepartmentDummy_;
                vector<vector<impalib_type>> extrinsicOutputDepartment_;
                vector<impalib_type> oric2PackageM_;
                vector<vector<impalib_type>> eqConstraint2OricM_;
                vector<vector<impalib_type>> oric2EqConstraintM_;
                vector<vector<impalib_type>> eqConstraint2ProjectM_;
                vector<vector<impalib_type>> project2EqConstraintM_;
                vector<impalib_type> team2OricM_;
                Knapsack modelKnapsacks_;
                InequalityConstraint projectIneqConstraint_;
                EqualityConstraint modelEqConstraint_; 
                OrInequalityConstraint modelOric_;

        public: 
                OutputsImpa outputs;
                InputsImpa modelInputs_;
                void initialize(const impalib_type*, impalib_type*, const int*,
                        const int*, const int*, const impalib_type*, 
                        const int*);
                void iterate(const int*);
                GraphicalModel(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG, const impalib_type ALPHA);
};

GraphicalModel::GraphicalModel(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG, const impalib_type ALPHA)
                                : modelInputs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS)
                                , outputs(N_DEPARTMENTS, N_TEAMS, N_PROJECTS), modelKnapsacks_(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA), projectIneqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS)
                                , modelEqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS), modelOric_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS) {

        numProjects_ = N_PROJECTS;
        numTeams_ = N_TEAMS;
        numDepartments_ = N_DEPARTMENTS;
        maxSizeNonZeroWeights_ = MAX_SIZE_NON_ZERO_WEIGHTS;
        numIterations_ = N_ITERATIONS;
        filteringFlag_ = FILT_FLAG;
        alpha_ = ALPHA;

        oric2PackageM_.resize(numTeams_);
        team2OricM_.resize(numTeams_);
        fill(oric2PackageM_.begin(), oric2PackageM_.begin()+numTeams_, zero_value);
        fill(team2OricM_.begin(), team2OricM_.begin()+numTeams_, zero_value);

        eqConstraint2OricM_.reserve(numProjects_);
        oric2EqConstraintM_.reserve(numProjects_);
        eqConstraint2ProjectM_.reserve(numProjects_);
        project2EqConstraintM_.reserve(numProjects_);

        for (int project_index = 0; project_index < numProjects_; project_index++){
                eqConstraint2OricM_.push_back(vector<impalib_type>(numTeams_,zero_value));
                oric2EqConstraintM_.push_back(vector<impalib_type>(numTeams_,zero_value));
                eqConstraint2ProjectM_.push_back(vector<impalib_type>(numTeams_,zero_value));
                project2EqConstraintM_.push_back(vector<impalib_type>(numTeams_,zero_value));
                }

        
        extrinsicOutputDepartment_.reserve(numDepartments_);
        extrinsicOutputDepartmentDummy_.reserve(numDepartments_);
        
        for (int department_index = 0; department_index < numDepartments_; department_index++){
                extrinsicOutputDepartment_.push_back(vector<impalib_type>(numTeams_,zero_value));
                extrinsicOutputDepartmentDummy_.push_back(vector<impalib_type>(numTeams_,zero_value));
                }
        
};

void GraphicalModel::initialize(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                    const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, 
                    const int *pMAX_STATE_PY){

        modelInputs_.process_inputs(pREWARD_TEAM_PY, pTransition_model_py, pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY,
                        pMAX_STATE_PY);
}


void GraphicalModel::iterate(const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY){
        for (int iter = 0; iter <numIterations_; iter++){

                for (int department_index = 0; department_index<numDepartments_; department_index++){
                
                int max_state_department = modelInputs_.MaxState[department_index];
                
                vector<vector<impalib_type>> stage_forward_messages(numTeams_+1, vector<impalib_type>(max_state_department+1,zero_value));
                vector<vector<impalib_type>> stage_backward_messages(numTeams_+1, vector<impalib_type>(max_state_department+1,zero_value));
                
                modelKnapsacks_.forward(department_index, stage_forward_messages, max_state_department,
                        modelInputs_.NonZeroWeightIndices, pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                        modelInputs_.TeamsWeightsPerDepartment, modelInputs_.Team2KnapsackM);
                
                modelKnapsacks_.backward(department_index, stage_backward_messages, max_state_department,
                                modelInputs_.NonZeroWeightIndices, pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
                        modelInputs_.TeamsWeightsPerDepartment, modelInputs_.Team2KnapsackM);
                
                modelKnapsacks_.extrinsic_output_department_lhs(modelInputs_.TeamsWeightsPerDepartment,
                                stage_forward_messages, modelInputs_.Team2KnapsackM,
                                department_index, stage_backward_messages,
                                max_state_department, extrinsicOutputDepartmentDummy_);

                modelKnapsacks_.process_extrinsic_output_department(department_index, iter, extrinsicOutputDepartmentDummy_, extrinsicOutputDepartment_);
                }

                modelEqConstraint_.team_eq_constraint_to_oric_update(extrinsicOutputDepartment_, team2OricM_, modelInputs_.RewardTeam);

                modelOric_.oric_to_project_eq_constraint_update(eqConstraint2OricM_, team2OricM_,
                                oric2EqConstraintM_, eqConstraint2ProjectM_, modelInputs_.RewardProject);

                projectIneqConstraint_.project_inequality_constraint_update(eqConstraint2ProjectM_, project2EqConstraintM_);

                modelEqConstraint_.project_eq_constraint_to_oric_update(project2EqConstraintM_, eqConstraint2OricM_, modelInputs_.RewardProject);

                modelOric_.oric_to_team_update(eqConstraint2OricM_, oric2PackageM_);

                outputs.intrinsic_out_mwm_update(oric2EqConstraintM_, project2EqConstraintM_, modelInputs_.RewardProject);

                modelKnapsacks_.team_to_knapsack_update(modelInputs_.NonZeroWeightIndices,
                                        modelInputs_.Team2KnapsackM,
                                        modelInputs_.RewardTeam, extrinsicOutputDepartment_,
                                        oric2PackageM_);
        }
                outputs.extrinsic_output_team_update(extrinsicOutputDepartment_, oric2PackageM_);

}