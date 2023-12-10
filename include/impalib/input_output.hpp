// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class InputsKcMwm{
    private:
        int numTeams_;
        int numDepartments_;
        int numProjects_;
        int maxSizeNonzeroWeights_;
    public: 
        vector<impalib_type> RewardTeam;
        vector<vector<impalib_type>> Team2KnapsackM;
        vector<vector<int>> TeamsWeightsPerDepartment;
        vector<vector<impalib_type>> RewardProject;
        vector<int> MaxState;
        vector<vector<int>> NonZeroWeightIndices;
        
        void process_inputs(const impalib_type*,
                    impalib_type*, const int*,
                    const int*, const int*, const impalib_type*, 
                    const int*);
        
        InputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS);
};

class OutputsKcMwm{
    private:
        int numTeams_;
        int numDepartments_;
        int numProjects_;

    public:
        vector<impalib_type> ExtrinsicOutputTeam;
        vector<impalib_type> IntrinsicOutMwm;
        void intrinsic_out_mwm_update(vector<vector<impalib_type>>&, vector<vector<impalib_type>>&, vector<vector<impalib_type>>&);
        void extrinsic_output_team_update(vector<vector<impalib_type>>&, vector<impalib_type>&);
        OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS);
};

InputsKcMwm::InputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS,const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS){
        
        numDepartments_ = N_DEPARTMENTS;
        numTeams_ = N_TEAMS;
        numProjects_ = N_PROJECTS;
        maxSizeNonzeroWeights_ = MAX_SIZE_NON_ZERO_WEIGHTS;

        TeamsWeightsPerDepartment.reserve(numDepartments_);
        NonZeroWeightIndices.reserve(numDepartments_);
        RewardProject.reserve(numProjects_);
        Team2KnapsackM.reserve(numDepartments_);

        for (int department_index = 0; department_index < numDepartments_; department_index++){
            Team2KnapsackM.push_back(vector<impalib_type>(numTeams_,zero_value));
            TeamsWeightsPerDepartment.push_back(vector<int>(numTeams_,0));
        }

        for (int project_index=0; project_index < numProjects_; project_index++){
            RewardProject.push_back(vector<impalib_type>(numTeams_,zero_value));
        }
};

void InputsKcMwm:: process_inputs(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                    const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, 
                    const int *pMAX_STATE_PY){

    copy(pREWARD_TEAM_PY, pREWARD_TEAM_PY + numTeams_,back_inserter(RewardTeam));
    copy(pMAX_STATE_PY, pMAX_STATE_PY + numDepartments_,back_inserter(MaxState));

    for (int department_index = 0; department_index < numDepartments_; department_index++){
        NonZeroWeightIndices.push_back(vector<int>(pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index],0));
        copy ( pTransition_model_py + numTeams_*department_index, pTransition_model_py+numTeams_*(department_index+1), Team2KnapsackM[department_index].begin() );
        copy ( pTEAMS_WEIGHTS_PER_DEPARTMENT_PY + numTeams_*department_index, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY+numTeams_*(department_index+1), TeamsWeightsPerDepartment[department_index].begin() );
        copy ( p_NON_ZERO_WEIGHT_INDICES_PY + maxSizeNonzeroWeights_*department_index, 
        p_NON_ZERO_WEIGHT_INDICES_PY+pNON_ZERO_WEIGHT_INDICES_SIZES_PY[department_index]+ maxSizeNonzeroWeights_*department_index, NonZeroWeightIndices[department_index].begin() );
    }

    
    for (int project_index = 0; project_index < numProjects_; project_index++){
        copy ( pREWARD_PROJECT_PY + numTeams_*project_index, pREWARD_PROJECT_PY+numTeams_*(project_index+1), RewardProject[project_index].begin() );
    }
    
}

OutputsKcMwm::OutputsKcMwm(const int N_DEPARTMENTS, const int N_TEAMS,const int N_PROJECTS){
        
        numDepartments_ = N_DEPARTMENTS;
        numTeams_ = N_TEAMS;
        numProjects_ = N_PROJECTS;
        
        ExtrinsicOutputTeam.reserve(numTeams_);
        ExtrinsicOutputTeam.resize(numTeams_);
        fill(ExtrinsicOutputTeam.begin(), ExtrinsicOutputTeam.begin()+numTeams_, zero_value);

        IntrinsicOutMwm.reserve(numProjects_*numTeams_);
        IntrinsicOutMwm.resize(numProjects_*numTeams_);
        fill(IntrinsicOutMwm.begin(), IntrinsicOutMwm.begin()+numProjects_*numTeams_, zero_value);

};

void OutputsKcMwm::intrinsic_out_mwm_update(vector<vector<impalib_type>> &rOric2EqConstraintM, vector<vector<impalib_type>> &rProject2EqConstraintM, vector<vector<impalib_type>>& rRewardProject){

        for (int project_index = 0; project_index < rRewardProject.size(); project_index++){
            for (int team_index = 0; team_index < rRewardProject[project_index].size(); team_index++){
                IntrinsicOutMwm[project_index+team_index+project_index*(numTeams_-1)] = rOric2EqConstraintM[project_index][team_index] + 
                                                    rProject2EqConstraintM[project_index][team_index] + rRewardProject[project_index][team_index];
            }
        }
    }


void OutputsKcMwm::extrinsic_output_team_update( vector<vector<impalib_type>> &rExtrinsicOutputDepartment,
                                    vector<impalib_type> &rOric2TeamM){

        copy(rOric2TeamM.begin(), rOric2TeamM.end(), ExtrinsicOutputTeam.begin());

        for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++){
            transform(rExtrinsicOutputDepartment[department_index].begin(), rExtrinsicOutputDepartment[department_index].end(), 
                                    ExtrinsicOutputTeam.begin(), ExtrinsicOutputTeam.begin(), std::plus<impalib_type>());
        }
}

class InputsTsp{
    private:
        int numNodes_;
        int numEdgeVariables_;
        int numNodesPerEdge_;
    public: 
        vector<vector<int>> EdgeConnections;
        vector<impalib_type> CostEdgeVariable;
        vector<vector<impalib_type>> CostMatrix;
        vector<vector<impalib_type>> EdgeDegreeConstraintCost;
        vector<vector<impalib_type>> EdgeEc2DegreeConstraintM;
        
        void process_inputs(const int*, const impalib_type*, const impalib_type*, impalib_type*, const impalib_type*);
        
        InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES);
};


InputsTsp::InputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES){
        
        numNodes_ = NUM_NODES;
        numEdgeVariables_ = NUM_EDGE_VARIABLES;
        numNodesPerEdge_ = 2;

};

class OutputsTsp{
    private:
        int numNodes_;
        int numEdgeVariables_;
    public:
        vector<impalib_type> ExtrinsicOutputEdgeEc;
        vector<impalib_type> IntrinsicOutputEdgeEc;
        void extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>>&);
        void extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>>&, vector<vector<impalib_type>>&);
        void intrinsic_output_edge_ec_update(vector<impalib_type>&);
        OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES);
};

OutputsTsp::OutputsTsp(const int NUM_NODES, const int NUM_EDGE_VARIABLES){
        
        numNodes_ = NUM_NODES;
        numEdgeVariables_ = NUM_EDGE_VARIABLES;

        ExtrinsicOutputEdgeEc.reserve(numEdgeVariables_);
        ExtrinsicOutputEdgeEc.resize(numEdgeVariables_);
        IntrinsicOutputEdgeEc.reserve(numEdgeVariables_);
        IntrinsicOutputEdgeEc.resize(numEdgeVariables_);
        fill(ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.begin()+numEdgeVariables_, zero_value);
        fill(IntrinsicOutputEdgeEc.begin(), IntrinsicOutputEdgeEc.begin()+numEdgeVariables_, zero_value);

};


void InputsTsp:: process_inputs(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY, impalib_type *pEdge_ec_to_degree_constraint_m_py, \
                                    const impalib_type* pEDGE_DEGREE_CONSTRAINT_COST_PY){

    copy(pCOST_EDGE_VARIABLE_PY, pCOST_EDGE_VARIABLE_PY + numEdgeVariables_,back_inserter(CostEdgeVariable));
    
    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++){
        EdgeConnections.push_back(vector<int>(numNodesPerEdge_,0));
        copy ( pEDGE_CONNECTIONS_PY + numNodesPerEdge_*edge_variable_index, pEDGE_CONNECTIONS_PY+numNodesPerEdge_*(edge_variable_index+1), EdgeConnections[edge_variable_index].begin() );
        EdgeEc2DegreeConstraintM.push_back(vector<impalib_type>(numNodes_,zero_value));
        copy ( pEdge_ec_to_degree_constraint_m_py + numNodes_*edge_variable_index, pEdge_ec_to_degree_constraint_m_py+numNodes_*(edge_variable_index+1), EdgeEc2DegreeConstraintM[edge_variable_index].begin() );
        EdgeDegreeConstraintCost.push_back(vector<impalib_type>(numNodes_,zero_value));
        copy ( pEDGE_DEGREE_CONSTRAINT_COST_PY + numNodes_*edge_variable_index, pEDGE_DEGREE_CONSTRAINT_COST_PY+numNodes_*(edge_variable_index+1), EdgeDegreeConstraintCost[edge_variable_index].begin() );
    }

    for (int node_index = 0; node_index < numNodes_; node_index++){
        CostMatrix.push_back(vector<impalib_type>(numNodes_,zero_value));
        copy ( pCOST_MATRIX_PY + numNodes_*node_index, pCOST_MATRIX_PY+numNodes_*(node_index+1), CostMatrix[node_index].begin() );
    }
}

void OutputsTsp::extrinsic_output_edge_ec_relaxed_graph_update(vector<vector<impalib_type>>& rDegreeConstraint2EqConstraintM){
    
        for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size(); edge_variable_index++){
            ExtrinsicOutputEdgeEc[edge_variable_index] = accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(), rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
        }

}

void OutputsTsp::extrinsic_output_edge_ec_augmented_graph_update(vector<vector<impalib_type>>& rDegreeConstraint2EqConstraintM, vector<vector<impalib_type>>& rSubtourConstraints2EdgeEcM){
    
        for (int edge_variable_index = 0; edge_variable_index < rDegreeConstraint2EqConstraintM.size(); edge_variable_index++){
            ExtrinsicOutputEdgeEc[edge_variable_index] = accumulate(rDegreeConstraint2EqConstraintM[edge_variable_index].begin(), rDegreeConstraint2EqConstraintM[edge_variable_index].end(), zero_value);
        }

        for (int subtour_constraint_index = 0; subtour_constraint_index<rSubtourConstraints2EdgeEcM.size(); subtour_constraint_index++){
            transform(rSubtourConstraints2EdgeEcM[subtour_constraint_index].begin(), rSubtourConstraints2EdgeEcM[subtour_constraint_index].end(), 
                                    ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.begin(), plus<impalib_type>());  
        }
}

void OutputsTsp::intrinsic_output_edge_ec_update(vector<impalib_type>& rCostEdgeVariable){
    
    transform(ExtrinsicOutputEdgeEc.begin(), ExtrinsicOutputEdgeEc.end(), rCostEdgeVariable.begin(), IntrinsicOutputEdgeEc.begin(), plus<impalib_type>());
}