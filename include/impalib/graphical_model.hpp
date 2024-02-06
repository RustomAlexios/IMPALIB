// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"
class GraphicalModelKcMwm{
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
                EqualityConstraintKcMwm modelEqConstraint_; 
                OrInequalityConstraint modelOric_;

        public: 
                OutputsKcMwm outputs;
                InputsKcMwm modelInputs_;
                void initialize(const impalib_type*, impalib_type*, const int*,
                        const int*, const int*, const impalib_type*, 
                        const int*);
                void iterate(const int*);
                GraphicalModelKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG, const impalib_type ALPHA);
};

GraphicalModelKcMwm::GraphicalModelKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS, const int MAX_SIZE_NON_ZERO_WEIGHTS, const int N_ITERATIONS, const bool FILT_FLAG, const impalib_type ALPHA)
                                : modelInputs_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS)
                                , outputs(N_DEPARTMENTS, N_TEAMS, N_PROJECTS), modelKnapsacks_(N_DEPARTMENTS, N_TEAMS, FILT_FLAG, ALPHA), projectIneqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS)
                                , modelEqConstraint_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS), modelOric_(N_DEPARTMENTS, N_TEAMS, N_PROJECTS), \
                                numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS), maxSizeNonZeroWeights_(MAX_SIZE_NON_ZERO_WEIGHTS), numIterations_(N_ITERATIONS), \
                                filteringFlag_(FILT_FLAG), alpha_(ALPHA){

        //numProjects_ = N_PROJECTS;
        //numTeams_ = N_TEAMS;
        //numDepartments_ = N_DEPARTMENTS;
        //maxSizeNonZeroWeights_ = MAX_SIZE_NON_ZERO_WEIGHTS;
        //numIterations_ = N_ITERATIONS;
        //filteringFlag_ = FILT_FLAG;
        //alpha_ = ALPHA;

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

void GraphicalModelKcMwm::initialize(const impalib_type *pREWARD_TEAM_PY, impalib_type *pTransition_model_py, const int *pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                    const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, const impalib_type *pREWARD_PROJECT_PY, 
                    const int *pMAX_STATE_PY){

        modelInputs_.process_inputs(pREWARD_TEAM_PY, pTransition_model_py, pITEMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY,
                        pMAX_STATE_PY);
}


void GraphicalModelKcMwm::iterate(const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY){
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


class GraphicalModelTsp{
        private:
                int numNodes_;
                int numIterations_;
                int numEdgeVariables_;
                bool resetFlag_;
                bool filteringFlag_;
                bool augmentationFlag_;
                vector<int> hard_decision;
                impalib_type alpha_;
                impalib_type threshold_;
                EqualityConstraintTsp modelEqConstraint_; 
                DegreeConstraint modelDegreeConstraint_;
                SubtourEliminationConstraint modelSubtourEliminationConstraint_;
                vector<vector<impalib_type>> DegreeConstraint2EqConstraintDummyM_;
                vector<vector<impalib_type>> DegreeConstraint2EqConstraintM_;
                vector<vector<int>> selected_edges_old_;
                vector<vector<int>> selected_edges_old_old_;
                int maxCount_;
                bool tourImpaFlag_ = false;
                vector<vector<int>> delta_S_indices_list;
                vector<vector<impalib_type>> subtourConstraints2EdgeEcDummyM_;
                vector<vector<impalib_type>> edgeEc2SubtourConstraintsM_;
                void iterate_augmented_graph();
                void subtour_elimination_constraints_analysis(unordered_map<int, vector<int>>&, vector<vector<int>>&);
                void hard_decision_analysis(vector<vector<int>> &);
                bool isSubsequence(const vector<int>&, const vector<int>&, int);
                vector<vector<int>> get_closed_loops(unordered_map<int, vector<int>>&, vector<vector<int>>&);
                vector<int> find_closed_loop(unordered_map<int, vector<int>>&, int, int, unordered_set<int>, vector<int>, vector<int>&);
                InputsTsp modelInputs_;
                vector<vector<int>> selectedEdges_;
                int numAugmentations_;
                int noConsClosedLoopsCount_;
                int solOscCount_;
                bool noConsClosedLoopsCountExcFlag_;
                int noImprovSolCount_;
                bool noImprovSolCountExcFlag_;
                bool solOscCountExcFlag_;
                vector<int> tourImpa_;
                vector<vector<int>> subtourPaths_;
                vector<int> closedPathsSize_;
                impalib_type costImpa_;
                vector<vector<impalib_type>> subtourConstraints2EdgeEcM_;
        public:
                bool subtourConstraintsSatisfiedFlag = false;
                OutputsTsp outputs;
                void initialize(const int*, const impalib_type*, const impalib_type*, impalib_type*, const impalib_type*);
                void iterate_relaxed_graph();
                void perform_augmentation(const int);
                void process_ouputs(impalib_type*, int*, int*, int*, impalib_type*, bool *, bool *, bool *, int*, int*, int*, int*);
                GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG, const impalib_type ALPHA, \
                        const impalib_type THRESHOLD, const int MAX_COUNT);
};


GraphicalModelTsp::GraphicalModelTsp(const int NUM_ITERATIONS, const int NUM_NODES, const int NUM_EDGE_VARIABLES, const bool AUGMENTATION_FLAG, const bool RESET_FLAG, const bool FILTERING_FLAG, const impalib_type ALPHA, \
                                const impalib_type THRESHOLD, const int MAX_COUNT)
                                : modelDegreeConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
                                modelEqConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
                                modelSubtourEliminationConstraint_(NUM_NODES, NUM_EDGE_VARIABLES, FILTERING_FLAG, ALPHA),
                                modelInputs_(NUM_NODES, NUM_EDGE_VARIABLES),
                                outputs(NUM_NODES, NUM_EDGE_VARIABLES), numIterations_(NUM_ITERATIONS), filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), \
                                numNodes_(NUM_NODES), augmentationFlag_(AUGMENTATION_FLAG), resetFlag_(RESET_FLAG), numEdgeVariables_(NUM_EDGE_VARIABLES), threshold_(THRESHOLD), \
                                maxCount_(MAX_COUNT)
        {
        //numIterations_ = NUM_ITERATIONS;
        //filteringFlag_ = FILTERING_FLAG;
        //alpha_ = ALPHA;
        //numNodes_ = NUM_NODES;
        //augmentationFlag_ = AUGMENTATION_FLAG;
        //resetFlag_ = RESET_FLAG;
        //numEdgeVariables_ = NUM_EDGE_VARIABLES;
        //threshold_ = THRESHOLD;
        //maxCount_ = MAX_COUNT;

        DegreeConstraint2EqConstraintDummyM_.reserve(numEdgeVariables_);
        DegreeConstraint2EqConstraintM_.reserve(numEdgeVariables_);

        //tourImpaFlag = false;
        //subtourConstraintsSatisfiedFlag = false;

        for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++){
                DegreeConstraint2EqConstraintDummyM_.push_back(vector<impalib_type>(numNodes_,zero_value));
                DegreeConstraint2EqConstraintM_.push_back(vector<impalib_type>(numNodes_,zero_value));
        }
        numAugmentations_ = 0;
        costImpa_ = zero_value;
        noConsClosedLoopsCount_ = 0;
        noImprovSolCount_ = 0;
        solOscCount_ = 0;

        noConsClosedLoopsCountExcFlag_ = false;
        noImprovSolCountExcFlag_ = false;
        solOscCountExcFlag_ = false;
        
};


void GraphicalModelTsp::initialize(const int *pEDGE_CONNECTIONS_PY, const impalib_type *pCOST_EDGE_VARIABLE_PY, const impalib_type *pCOST_MATRIX_PY, impalib_type *pEdge_ec_to_degree_constraint_m_py, \
                                const impalib_type* pEDGE_DEGREE_CONSTRAINT_COST_PY){
        modelInputs_.process_inputs(pEDGE_CONNECTIONS_PY, pCOST_EDGE_VARIABLE_PY, pCOST_MATRIX_PY, pEdge_ec_to_degree_constraint_m_py, pEDGE_DEGREE_CONSTRAINT_COST_PY);
}


void GraphicalModelTsp::iterate_relaxed_graph(){

        for (int iter = 0; iter <numIterations_; iter++){
                
                modelDegreeConstraint_.degree_constraint_to_edge_ec_update(modelInputs_.EdgeEc2DegreeConstraintM, modelInputs_.EdgeConnections, DegreeConstraint2EqConstraintDummyM_);

                modelDegreeConstraint_.process_filtering(iter, DegreeConstraint2EqConstraintDummyM_, DegreeConstraint2EqConstraintM_);
                
                modelEqConstraint_.edge_ec_to_degree_constraint_relaxed_graph_update(modelInputs_.EdgeConnections, modelInputs_.EdgeDegreeConstraintCost ,DegreeConstraint2EqConstraintM_, modelInputs_.EdgeEc2DegreeConstraintM);
        }
        
        outputs.extrinsic_output_edge_ec_relaxed_graph_update(DegreeConstraint2EqConstraintM_);
        outputs.intrinsic_output_edge_ec_update(modelInputs_.CostEdgeVariable);

        vector<vector<int>> selected_edges;
        hard_decision_analysis(selected_edges);

        unordered_map<int, vector<int>> graph;
        
        subtour_elimination_constraints_analysis(graph, selected_edges); 
        selectedEdges_ = selected_edges;
}


void GraphicalModelTsp::perform_augmentation(const int MAX_AUGM_COUNT){

        if (delta_S_indices_list.size() > 0)
        {
            vector<vector<impalib_type>> temp(delta_S_indices_list.size(),
                                              vector<impalib_type>(numEdgeVariables_, zero_value));
            subtourConstraints2EdgeEcM_.insert(subtourConstraints2EdgeEcM_.end(), temp.begin(),
                                                          temp.end());
            subtourConstraints2EdgeEcDummyM_ = subtourConstraints2EdgeEcM_;
        }

    while (!subtourConstraintsSatisfiedFlag && augmentationFlag_ && !tourImpaFlag_)
    {

        if (noConsClosedLoopsCountExcFlag_ || noImprovSolCountExcFlag_
            || solOscCountExcFlag_)
        {
            cout << "noConsClosedLoopsCountExcFlag_: " << noConsClosedLoopsCountExcFlag_
                 << '\n';
            cout << "noImprovSolCountExcFlag_: " << noImprovSolCountExcFlag_ << '\n';
            cout << "solOscCountExcFlag_: " << solOscCountExcFlag_ << '\n';
            break;
        }
        if (costImpa_ == zero_value)
        {
            cout << "Cost is zero" << '\n';
            cout << "Possibly Nan" << '\n';
            exit(0);
        }

        iterate_augmented_graph();

        if (numAugmentations_ == MAX_AUGM_COUNT)
        {
            cout << "MAX_AUGM_COUNT reached" << '\n';
            break;
        }

        if (subtourConstraints2EdgeEcM_.size() != delta_S_indices_list.size())
        {
            size_t numLists2Add =
                delta_S_indices_list.size() - subtourConstraints2EdgeEcM_.size();
            vector<vector<impalib_type>> temp(numLists2Add, vector<impalib_type>(numEdgeVariables_, zero_value));
            subtourConstraints2EdgeEcM_.insert(subtourConstraints2EdgeEcM_.end(), temp.begin(),
                                                          temp.end());
            subtourConstraints2EdgeEcDummyM_ = subtourConstraints2EdgeEcM_;
        }
    }
}

void GraphicalModelTsp::process_ouputs(impalib_type* pExtrinsic_output_edge_ec, int* pNum_augmentations, int* pNum_added_constraints, int* pTour_impa, impalib_type* pCost_impa, \
                                        bool *pNo_improv_sol_count_exc_flag, bool *pNo_cons_loops_count_exc_flag, bool *pSol_osc_count_exc_flag, int* pSelected_edges, int* pSelected_edges_size, \
                                        int* pSubtour_paths, int* pSubtour_paths_size){
        
        copy(outputs.ExtrinsicOutputEdgeEc.begin(),
                outputs.ExtrinsicOutputEdgeEc.begin() + numEdgeVariables_, pExtrinsic_output_edge_ec);

        *pNum_augmentations     = numAugmentations_;

        *pNum_added_constraints = static_cast<int>(subtourConstraints2EdgeEcM_.size());

        copy(tourImpa_.begin(), tourImpa_.begin() + static_cast<int>(tourImpa_.size()),
                pTour_impa);

        *pCost_impa                    = costImpa_;
        *pNo_improv_sol_count_exc_flag = noImprovSolCountExcFlag_;
        *pNo_cons_loops_count_exc_flag = noConsClosedLoopsCountExcFlag_;
        *pSol_osc_count_exc_flag       = solOscCountExcFlag_;

        vector<int> flattened_selected_edges =
                accumulate(selectedEdges_.begin(), selectedEdges_.end(), vector<int>{},
                        [](vector<int> &acc, const vector<int> &inner)
                        {
                        acc.insert(acc.end(), inner.begin(), inner.end());
                        return acc;
                        });

        copy(flattened_selected_edges.begin(),
                flattened_selected_edges.begin() + static_cast<int>(flattened_selected_edges.size()), pSelected_edges);
        *pSelected_edges_size = static_cast<int>(flattened_selected_edges.size());

        vector<int> flattened_closed_paths =
                accumulate(subtourPaths_.begin(), subtourPaths_.end(), vector<int>{},
                        [](vector<int> &acc, const vector<int> &inner)
                        {
                        acc.insert(acc.end(), inner.begin(), inner.end());
                        return acc;
                        });
        copy(flattened_closed_paths.begin(),
                flattened_closed_paths.begin() + static_cast<int>(flattened_closed_paths.size()), pSubtour_paths);
        copy(closedPathsSize_.begin(),
                closedPathsSize_.begin() + static_cast<int>(closedPathsSize_.size()),
                pSubtour_paths_size);

}

void GraphicalModelTsp::iterate_augmented_graph(){
        
        if (resetFlag_){
                for (size_t i=0; i<modelInputs_.EdgeConnections.size(); i++){
                        auto connection = modelInputs_.EdgeConnections[i];
                        impalib_type cost = modelInputs_.CostMatrix[connection[0]][connection[1]];
                        modelInputs_.EdgeEc2DegreeConstraintM[i][connection[0]] = cost;
                        modelInputs_.EdgeEc2DegreeConstraintM[i][connection[1]] = cost;
                }
        }

        cout << "-----------------------" << '\n';
        cout << "iterate_augmented_graph" << '\n';
        cout << "-----------------------" << '\n';

        numAugmentations_+=1;
        cout<<"Augmentation count: "<<numAugmentations_<<'\n';
        cout<<"delta_S_indices_list.size(): "<<delta_S_indices_list.size()<<'\n';

        modelSubtourEliminationConstraint_.subtourConstraints2EdgeEcOld_.reserve(delta_S_indices_list.size());

        for (int index_subtour_constraint = 0; index_subtour_constraint < delta_S_indices_list.size(); index_subtour_constraint++) {
                modelSubtourEliminationConstraint_.subtourConstraints2EdgeEcOld_.push_back(vector<impalib_type>(numEdgeVariables_,zero_value));
        }
        
        for (int iter = 0; iter <numIterations_; iter++){
                modelDegreeConstraint_.degree_constraint_to_edge_ec_update(modelInputs_.EdgeEc2DegreeConstraintM, modelInputs_.EdgeConnections, DegreeConstraint2EqConstraintDummyM_);
                
                modelDegreeConstraint_.process_filtering(iter, DegreeConstraint2EqConstraintDummyM_, DegreeConstraint2EqConstraintM_);
                
                edgeEc2SubtourConstraintsM_ = modelEqConstraint_.edge_ec_to_subtour_constraints_update(delta_S_indices_list, modelInputs_.CostEdgeVariable, DegreeConstraint2EqConstraintM_, subtourConstraints2EdgeEcM_, modelInputs_.EdgeConnections);

                modelSubtourEliminationConstraint_.subtour_constraints_to_edge_ec_update(edgeEc2SubtourConstraintsM_, delta_S_indices_list, subtourConstraints2EdgeEcDummyM_);

                modelSubtourEliminationConstraint_.process_filtering(iter, subtourConstraints2EdgeEcDummyM_, subtourConstraints2EdgeEcM_, delta_S_indices_list);
                
                modelEqConstraint_.edge_ec_to_degree_constraint_augmented_graph_update(DegreeConstraint2EqConstraintM_, subtourConstraints2EdgeEcM_, modelInputs_.EdgeConnections, modelInputs_.EdgeDegreeConstraintCost, modelInputs_.EdgeEc2DegreeConstraintM);          
        }


        bool flag_nan=false;
        for (int i = 0; i < modelInputs_.EdgeEc2DegreeConstraintM.size(); i++) {
                for (int j = 0; j < modelInputs_.EdgeEc2DegreeConstraintM[i].size(); j++) {
                        if (isnan(modelInputs_.EdgeEc2DegreeConstraintM[i][j])) {
                                flag_nan=true;
                        }
                }
        }

        if (flag_nan){
                cout << "NaN detected" << '\n';
                exit(0);
        }

        outputs.extrinsic_output_edge_ec_augmented_graph_update(DegreeConstraint2EqConstraintM_, subtourConstraints2EdgeEcM_);
        outputs.intrinsic_output_edge_ec_update(modelInputs_.CostEdgeVariable);
        
        selectedEdges_.clear();
        vector<vector<int>> selected_edges;
        hard_decision_analysis(selected_edges);
        if (!selected_edges_old_.empty()){
                bool areEqual = (selected_edges == selected_edges_old_);
                if (areEqual){
                        cout << "Exited: No Improvement of IMPA Solution" << '\n'; 
                        noImprovSolCount_+=1;
                        if (noImprovSolCount_ > maxCount_){
                        cout << "Exited: noImprovSolCount_ Exceeded maximum allowable number " << maxCount_<<'\n';
                        noImprovSolCountExcFlag_ = true;
                }
                }
                else{
                        noImprovSolCount_=0;
                }

        }


        if (!selected_edges_old_old_.empty()){
                bool oscillation_flag = ((selected_edges == selected_edges_old_old_) && (selected_edges_old_old_ != selected_edges_old_));
                if (oscillation_flag){
                        cout << "Exited: Solution is oscillating" << '\n'; 
                        solOscCount_+=1;
                        if (solOscCount_ > maxCount_){
                        cout << "Exited: solOscCount_ Exceeded maximum allowable number " << maxCount_<<'\n';
                        solOscCountExcFlag_ = true;
                }
                }
                else{
                        solOscCount_=0;
                }

        }

        if (!selected_edges_old_.empty()){
                selected_edges_old_old_ = selected_edges_old_;
        }

        selected_edges_old_ = selected_edges;
        selectedEdges_ = selected_edges;

        unordered_map<int, vector<int>> graph;
        
        subtour_elimination_constraints_analysis(graph, selected_edges);
}

void GraphicalModelTsp::hard_decision_analysis(vector<vector<int>> &rSelectedEdges){

        hard_decision.resize(numEdgeVariables_);
        fill(hard_decision.begin(), hard_decision.begin()+numEdgeVariables_, numeric_limits<int>::max());

        transform(outputs.IntrinsicOutputEdgeEc.begin(), outputs.IntrinsicOutputEdgeEc.end(), hard_decision.begin(),
                        [&](impalib_type value) { return value > threshold_ ? zero_value : 1; });
        
        for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++){
                if (hard_decision[edge_variable_index]  == 1) {
                rSelectedEdges.push_back(modelInputs_.EdgeConnections[edge_variable_index]);
                }
        }

        set<int> uniqueNodes;

        cout << "selected_edges: [";
        for(const vector<int>& edge : rSelectedEdges) {

                for(int node : edge) {
                        if (uniqueNodes.find(node) == uniqueNodes.end()) {
                                uniqueNodes.insert(node);
                        }
                }
                cout << "[";
                for(size_t i = 0; i < edge.size(); i++) {
                        cout << edge[i];
                        if(i != edge.size() - 1) {
                        cout << ", ";
                        }
                }

        cout << "]";

        if(&edge != &rSelectedEdges.back()) {
                cout << ", ";
        }
        }
        cout << "]" << '\n';
        //cout << "Number of activated nodes: " << uniqueNodes.size()  << '\n';

        impalib_type cost_impa = zero_value;
        for(size_t i = 0; i < hard_decision.size(); ++i) {
                if(hard_decision[i] == 1) {
                        cost_impa += modelInputs_.CostEdgeVariable[i];
                }
        }

        cout<<"C++ cost_impa: "<<cost_impa<<'\n';
        costImpa_ = cost_impa;
}


void GraphicalModelTsp:: subtour_elimination_constraints_analysis(unordered_map<int, vector<int>>& rGraph, vector<vector<int>>& rSelectedEdges){
        
        closedPathsSize_.clear(); //just store it at the end if applicable
        vector<vector<int>> loops_list = get_closed_loops(rGraph, rSelectedEdges);
        subtourPaths_.clear();
        subtourPaths_ = loops_list;
        
        
        if (loops_list.empty()){
                cout << "Exited: Cannot find tours using get_closed_loops()" << '\n';
                noConsClosedLoopsCount_+=1;
                if (noConsClosedLoopsCount_ > maxCount_){
                        cout << "Exited: noConsClosedLoopsCount_ Exceeded maximum allowable number: " << maxCount_<<'\n';
                        noConsClosedLoopsCountExcFlag_ = true;
                }
        }
        else if (loops_list.size()==1 && loops_list[0].size()==numNodes_){
                tourImpaFlag_ = true;
                subtourConstraintsSatisfiedFlag = true;
                tourImpa_ = loops_list[0];
                tourImpa_.push_back(loops_list[0][0]);
                cout << "Tour found" << '\n';
                
        cout << "tour_impa: [";
        for (size_t i = 0; i < tourImpa_.size(); ++i) {
                cout << tourImpa_[i];
                if (i != tourImpa_.size() - 1) {
                cout << ", ";
                }
        }
        cout << "]" << '\n';
        }

        else{
                noConsClosedLoopsCount_=0;
                for (const auto& loop : loops_list) {
                        closedPathsSize_.push_back(static_cast<int>(loop.size()));
                        cout << "Subtour of size " << loop.size() << " detected @: ";
                        cout << "[";
                        for (size_t i = 0; i < loop.size(); ++i) {
                                cout << loop[i];
                                if (i != loop.size() - 1) {
                                        cout << ", ";
                                }
                        }
                        cout << "]" << '\n';

                        vector<int> delta_S_indices;
                        for (size_t i = 0; i < modelInputs_.EdgeConnections.size(); ++i) {
                                const auto& connection = modelInputs_.EdgeConnections[i];
                                if (find(loop.begin(), loop.end(), connection[0]) != loop.end() && find(loop.begin(), loop.end(), connection[1]) == loop.end()) {
                                        delta_S_indices.push_back(static_cast<int>(i));    
                                }
                        }

                        if (delta_S_indices_list.empty() && !delta_S_indices.empty()){ //added !delta_S_indices.empty() to make sure if a full tour was detected so it is not be added with subtours
                                delta_S_indices_list.push_back(delta_S_indices);
                        }
                        else{
                                bool found = false;
                                for (const auto& existing_indices: delta_S_indices_list){
                                        if (existing_indices == delta_S_indices){
                                                found =true;
                                                break;
                                        }
                                }
                                if (!found && !delta_S_indices.empty()){//added !delta_S_indices.empty() to make sure if a full tour was detected so it is not be added with subtours
                                        delta_S_indices_list.push_back(delta_S_indices);
                                }
                        }
                }
        }

}

vector<vector<int>> GraphicalModelTsp:: get_closed_loops(unordered_map<int, vector<int>> &rGraph, vector<vector<int>>& rSelectedEgdes){
        
        for (const auto& connection : rSelectedEgdes) {
                if (rGraph.find(connection[0]) != rGraph.end()) {
                        rGraph[connection[0]].push_back(connection[1]);
                } else {
                        rGraph[connection[0]] = {connection[1]};
                }
        }

        vector<vector<int>> loops_list;
        vector<int> visited_nodes;

        for (const auto& connection : rGraph) {
                int start_node = connection.first;
                const vector<int>& end_nodes = connection.second;
                if (find(visited_nodes.begin(), visited_nodes.end(), start_node) != visited_nodes.end()) {
                        continue;
                } else {
                for (int j = 0; j < end_nodes.size(); j++) {
                        unordered_set<int> visited_set;
                        vector<int> path;
                        vector<int> closed_loop = find_closed_loop(rGraph, start_node, end_nodes[j], visited_set, path, visited_nodes);
                        if (!closed_loop.empty()){
                                loops_list.push_back(closed_loop);}
                }
                }
        }


        vector<int> list_indices_to_remove;
        vector<vector<int>> new_loops_list;

        for (int i=0; i<static_cast<int>(loops_list.size()); i++){
                const auto& list_1 = loops_list[i];
                if (list_1.empty()){
                        list_indices_to_remove.push_back(i);
                        continue;       
                }
                else if (find(list_indices_to_remove.begin(), list_indices_to_remove.end(), i) != list_indices_to_remove.end()){
                        continue;
                }
                vector<int> double_list_1 = list_1;
                double_list_1.insert(double_list_1.end(), list_1.begin(), list_1.end());
                for (int j=i+1; j<static_cast<int>(loops_list.size()); j++){
                        const auto& list_2 = loops_list[j];
                        
                        if (list_2.empty()){
                                list_indices_to_remove.push_back(j);
                        }

                        if (double_list_1.size() <list_2.size()){ //since if seq.size() - subseq.size()<0, this would give an arbitrary large number and core dumped
                                continue;
                        }

                        else if (isSubsequence(double_list_1, list_2, j)){
                                list_indices_to_remove.push_back(j);
                        }
                }
        }
        
        
        for (size_t i = 0; i < loops_list.size(); ++i) {
                if (find(list_indices_to_remove.begin(), list_indices_to_remove.end(), i) == list_indices_to_remove.end()) {
                        new_loops_list.push_back(loops_list[i]);
                }
        }

        return new_loops_list;
}

bool GraphicalModelTsp::isSubsequence(const vector<int>& seq, const vector<int>& subseq, int j) {

        for (int i = 0; i < static_cast<int>(seq.size() - subseq.size()); ++i) {
                if (equal(seq.begin() + i, seq.begin() + i + static_cast<int>(subseq.size()), subseq.begin())) {
                        return true;
                }
        }
        return false;
}        

vector<int> GraphicalModelTsp:: find_closed_loop(unordered_map<int, vector<int>> &rGraph, int start_node, int current, unordered_set<int> visited, vector<int> path, vector<int>& visited_nodes){
        
        visited.insert(current);
        path.push_back(current);
        
        if (current == start_node) {
                return path;
        }

        if (rGraph.find(current) == rGraph.end()) {
                vector<int> new_path = {start_node};
                new_path.insert(new_path.end(), path.begin(), path.end());
                return vector<int>();
        }

        for (int node : rGraph.at(current)) {
                if (visited.find(node) == visited.end()) {
                        vector<int> new_path = find_closed_loop(rGraph, start_node, node, visited, path, visited_nodes);
                        if (!new_path.empty()) {
                        return new_path;
                }
                }
        }
        
        return vector<int>();

}