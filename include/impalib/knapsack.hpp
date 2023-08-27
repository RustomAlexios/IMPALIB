// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class Knapsack{
    private:
        int numDepartments_;
        int numTeams_;
        bool filteringFlag_;
        impalib_type alpha_;
        vector<vector<impalib_type>> extrinsicOutputDepartmentOld_;

    public:
        void forward(int , vector<vector<impalib_type>>&, int , 
                    vector<vector<int>>&, const int*, 
                    vector<vector<int>>&, vector<vector<impalib_type>>&  );
        
        void backward(int , vector<vector<impalib_type>>& , int , 
                    vector<vector<int>>&, const int*, 
                    vector<vector<int>>&, vector<vector<impalib_type>>&  );
        
        void extrinsic_output_department_lhs(vector<vector<int>>&,
                        vector<vector<impalib_type>>&, vector<vector<impalib_type>>&,
                        int, vector<vector<impalib_type>>&,
                        int, vector<vector<impalib_type>>&);
        
        void team_to_knapsack_update(vector<vector<int>>&,
                                vector<vector<impalib_type>>&,
                                vector<impalib_type>&, vector<vector<impalib_type>>&,
                                vector<impalib_type>&);
        void process_extrinsic_output_department(int, int, vector<vector<impalib_type>>&, vector<vector<impalib_type>>&);

    Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA);
};

Knapsack::Knapsack(const int N_DEPARTMENTS, const int N_TEAMS, const bool FILT_FLAG, const impalib_type ALPHA){
        numDepartments_ = N_DEPARTMENTS;
        numTeams_ = N_TEAMS;
        filteringFlag_ = FILT_FLAG;
        alpha_ = ALPHA;

        extrinsicOutputDepartmentOld_.reserve(numDepartments_);

        for (int department_index = 0; department_index < numDepartments_; department_index++){
            extrinsicOutputDepartmentOld_.push_back(vector<impalib_type>(numTeams_,zero_value));
        }

};

void Knapsack::forward(int department_index, vector<vector<impalib_type>> &rStageForwardMessages, int max_state_department,
            vector<vector<int>> &rNonZeroWeightIndices, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
            vector<vector<int>> &rTeamsWeightsPerDepartment, vector<vector<impalib_type>> &rTeam2KnapsackM){
            
            vector<int>::iterator upper;
            
            vector<impalib_type> initial_forward_messages(max_state_department+1, zero_value);
            fill(initial_forward_messages.begin()+1, initial_forward_messages.end(), value_inf);

            rStageForwardMessages[0] = initial_forward_messages;

            if (rNonZeroWeightIndices[department_index][0]!=0){
                    for (int j=0; j<rNonZeroWeightIndices[department_index][0]; j++){
                        rStageForwardMessages[j+1] = initial_forward_messages;
                    }
            }
            
            for (int l=0; l<pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]; l++){
                int t = rNonZeroWeightIndices[department_index][l];
                for (int a = 0; a <= max_state_department; a++){
                    if (a-rTeamsWeightsPerDepartment[department_index][t]>=0 && (! (t==0 && a!=rTeamsWeightsPerDepartment[department_index][t]))){
                        rStageForwardMessages[t+1][a] = min(initial_forward_messages[a], initial_forward_messages[a-rTeamsWeightsPerDepartment[department_index][t]]+rTeam2KnapsackM[department_index][t]);
                    } else {
                        rStageForwardMessages[t+1][a] = initial_forward_messages[a];
                    }
                }
                
                initial_forward_messages = rStageForwardMessages[t+1];
                
                if (rNonZeroWeightIndices[department_index].size() != numTeams_){
                    if (!binary_search(rNonZeroWeightIndices[department_index].begin(), rNonZeroWeightIndices[department_index].end(), t+1)){
                        if (t+1 >=rNonZeroWeightIndices[department_index].back() && t+1<numTeams_){
                                for (int j=t+1; j<numTeams_; j++){
                                    rStageForwardMessages[j+1] = initial_forward_messages;
                                }
                        }else if (t+1 < rNonZeroWeightIndices[department_index].back()){
                            upper = upper_bound(rNonZeroWeightIndices[department_index].begin(), rNonZeroWeightIndices[department_index].end(), t);
                            int next_index = upper - rNonZeroWeightIndices[department_index].begin();
                                for (int j=t+1; j<rNonZeroWeightIndices[department_index][next_index]; j++){
                                    rStageForwardMessages[j+1] = initial_forward_messages;
                                }
                        }
                    }
                }
            }

}


void Knapsack::backward(int department_index, vector<vector<impalib_type>> &rStageBackwardMessages, int max_state_department,
            vector<vector<int>> &rNonZeroWeightIndices, const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY,
            vector<vector<int>> &rTeamsWeightsPerDepartment, vector<vector<impalib_type>> &rTeam2KnapsackM){

            vector<int>::iterator upper;
            vector<impalib_type> initial_backward_messages(max_state_department+1, zero_value);
            
            rStageBackwardMessages[numTeams_] = initial_backward_messages;
            
            if (rNonZeroWeightIndices[department_index][pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]-1]!= numTeams_-1){
                    for (int j=numTeams_-1; j>rNonZeroWeightIndices[department_index][0]; j--){
                        rStageBackwardMessages[j] = initial_backward_messages;
                    }
            }

            for (int l=pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]-1; l>=0; l--){
                int t = rNonZeroWeightIndices[department_index][l];
                for (int a = 0; a<=max_state_department; a++){
                    if ((t>0 && a+rTeamsWeightsPerDepartment[department_index][t]<=max_state_department) || (a==0 && t==0)){
                        rStageBackwardMessages[t][a] = min(initial_backward_messages[a], initial_backward_messages[a+rTeamsWeightsPerDepartment[department_index][t]]+rTeam2KnapsackM[department_index][t]);
                    } else {
                        rStageBackwardMessages[t][a] = initial_backward_messages[a];
                    }
                }
                
                initial_backward_messages = rStageBackwardMessages[t];

                if (pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]-1 != numTeams_){
                    if (!binary_search(rNonZeroWeightIndices[department_index].begin(), rNonZeroWeightIndices[department_index].end(), t-1)){
                        if ((t-1 <= rNonZeroWeightIndices[department_index][0]) && t-1>=0){
                                for (int j=t-1; j>=0; j--){
                                    rStageBackwardMessages[j] = initial_backward_messages;
                                }
                        }else if (t-1 > rNonZeroWeightIndices[department_index][0]){
                            upper = upper_bound(rNonZeroWeightIndices[department_index].begin(), rNonZeroWeightIndices[department_index].end(), t-1);
                            int next_index = upper - rNonZeroWeightIndices[department_index].begin();
                                for (int j=t-1; j>rNonZeroWeightIndices[department_index][next_index-1]; j--){
                                    rStageBackwardMessages[j] = initial_backward_messages;
                                }
                        }
                    }
                }
            }

            }

void Knapsack::extrinsic_output_department_lhs(vector<vector<int>> &rTeamsWeightsPerDepartment,
                            vector<vector<impalib_type>> &rStageForwardMessages, vector<vector<impalib_type>> &rTeam2KnapsackM,
                            int department_index, vector<vector<impalib_type>> &rStageBackwardMessages,
                            int max_state_department, vector<vector<impalib_type>> &rExtrinsicOutputDepartment){
        
        vector<impalib_type> metric_path_solid, metric_path_dash;

        for (int team_index = 0; team_index < rTeamsWeightsPerDepartment[department_index].size(); team_index++){
            metric_path_dash.clear(); metric_path_solid.clear();
            if (rTeamsWeightsPerDepartment[department_index][team_index] != 0){
                if (team_index==0){
                    metric_path_solid.push_back(rStageForwardMessages[team_index][team_index] + rStageBackwardMessages[team_index+1][rTeamsWeightsPerDepartment[department_index][team_index]]
                                                            +rTeam2KnapsackM[department_index][team_index]);
                    metric_path_dash.push_back(rStageForwardMessages[team_index][team_index]+rStageBackwardMessages[team_index+1][0]);
                }else{
                    for (int i = 0; i <= max_state_department; i++){
                        if (i+rTeamsWeightsPerDepartment[department_index][team_index] <= max_state_department){
                            metric_path_solid.push_back(rStageForwardMessages[team_index][i] +rStageBackwardMessages[team_index+1][i+rTeamsWeightsPerDepartment[department_index][team_index]]
                                                        +rTeam2KnapsackM[department_index][team_index]);    
                        }
                        metric_path_dash.push_back(rStageForwardMessages[team_index][i]+rStageBackwardMessages[team_index+1][i]);
                    }
                }
                rExtrinsicOutputDepartment[department_index][team_index] = *min_element(metric_path_solid.begin(), metric_path_solid.end())
                                                        - *min_element(metric_path_dash.begin(), metric_path_dash.end())
                                                        - rTeam2KnapsackM[department_index][team_index];
            }else{
                rExtrinsicOutputDepartment[department_index][team_index] = zero_value;
            }
            }

        }

void Knapsack::team_to_knapsack_update(vector<vector<int>> &rNonZeroWeightIndices,
                                vector<vector<impalib_type>> &rTeam2KnapsackM, 
                                vector<impalib_type> &rRewardTeam, vector<vector<impalib_type>> &rExtrinsicOutputDepartment,
                                vector<impalib_type> &mORIC2Team){
        
        for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++){
            vector<int> remaining_departments(numDepartments_);
            iota(remaining_departments.begin(), remaining_departments.end(),0);
            vector<int> unique_edge_department(rNonZeroWeightIndices[department_index]);
            remaining_departments.erase(remaining_departments.begin()+department_index);

            for (auto t : remaining_departments){
                
                vector<int> intersection(numTeams_);
                vector<int>::iterator it;
                it=set_intersection (rNonZeroWeightIndices[department_index].begin(), rNonZeroWeightIndices[department_index].end(),
                            rNonZeroWeightIndices[t].begin(), rNonZeroWeightIndices[t].end(), intersection.begin());
                intersection.resize(it-intersection.begin()); 
                for (auto l: intersection){
                    
                    rTeam2KnapsackM[department_index][l] =  rRewardTeam[l] + rExtrinsicOutputDepartment[t][l] + mORIC2Team[l];
                    vector<int>::iterator position = std::find(unique_edge_department.begin(), unique_edge_department.end(), l);
                    unique_edge_department.erase(position);
                }
                }
            
            for (auto u: unique_edge_department){
                rTeam2KnapsackM[department_index][u] = rRewardTeam[u] + mORIC2Team[u];
            }
            } 
}


void Knapsack:: process_extrinsic_output_department(int department_index, int iter, vector<vector<impalib_type>> &rExtrinsicOutputDepartmentDummy, vector<vector<impalib_type>> &rExtrinsicOutputDepartment){
            
        if ((filteringFlag_) and (alpha_ != zero_value)){
            
            vector<impalib_type> intermediate_dummy(rExtrinsicOutputDepartmentDummy[department_index]), intermediate_old(extrinsicOutputDepartmentOld_[department_index]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1-alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](impalib_type &c){return c*w_2;});
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](impalib_type &c){ return c*w_1; });
                if (iter==0 ){
                    copy(intermediate_dummy.begin(), intermediate_dummy.end(), rExtrinsicOutputDepartment[department_index].begin());
                } else {
                    transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), std::back_inserter(intermediate_extrinsic), std::plus<impalib_type>());
                    copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rExtrinsicOutputDepartment[department_index].begin());
                }
            
            copy(rExtrinsicOutputDepartment[department_index].begin(), rExtrinsicOutputDepartment[department_index].end(), extrinsicOutputDepartmentOld_[department_index].begin());
        }

        else{
            copy(rExtrinsicOutputDepartmentDummy[department_index].begin(), rExtrinsicOutputDepartmentDummy[department_index].end(), rExtrinsicOutputDepartment[department_index].begin());
        }
}