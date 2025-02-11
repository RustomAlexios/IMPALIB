// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class EqualityConstraint
{
private:
    int numTeams_; ///< number of teams
    int numProjects_; ///< number of projects
    int numDepartments_; ///< number of departments
    int          numNodes_; ///< number of nodes
    int          numEdgeVariables_; ///< number of edge connections
    bool         filteringFlag_; ///< filtering flag
    impalib_type alpha_; ///< filtering parameter
    int numVariables_;
    int numConstraints_;
    int kvariable_;
    int numFixedTx_;
    int numMobileTx_;
    int numBands_;
    int numTimeSteps_;
    int numRxLocs_;
    int numMobileTxLocs_;
    bool excludeCapFlag_;
public:

    EqualityConstraint(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS), 
        filteringFlag_(false), alpha_(zero_value), numNodes_(0), numEdgeVariables_(0),
        numVariables_(0), numConstraints_(0), kvariable_(0), 
        numFixedTx_(0), numMobileTx_(0), numBands_(0), numTimeSteps_(0),
        numRxLocs_(0), numMobileTxLocs_(0), excludeCapFlag_(0){};

    EqualityConstraint(const int NUM_NODES, const int NUM_EDGE_VARIABLES,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : numProjects_(0), numTeams_(0), numDepartments_(0),
      filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(NUM_NODES),
      numEdgeVariables_(NUM_EDGE_VARIABLES), numVariables_(0), numConstraints_(0), kvariable_(0),
      numFixedTx_(0), numMobileTx_(0), numBands_(0), numTimeSteps_(0),
      numRxLocs_(0), numMobileTxLocs_(0), excludeCapFlag_(0){};
    
    EqualityConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE,
                                             const bool FILTERING_FLAG, const impalib_type ALPHA)
    : numProjects_(0), numTeams_(0), numDepartments_(0),
      filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numNodes_(0),
      numEdgeVariables_(0), numVariables_(NUM_VARIABLES), numConstraints_(NUM_CONSTRAINTS), kvariable_(K_VARIABLE),
      numFixedTx_(0), numMobileTx_(0), numBands_(0), numTimeSteps_(0),
      numRxLocs_(0), numMobileTxLocs_(0), excludeCapFlag_(0){};

    EqualityConstraint(const int NUM_FIXED_TX, const int NUM_MOBILE_TX, const int NUM_BANDS, const int NUM_TIME_STEPS, const int NUM_RX_LOCS, 
                        const int NUM_MOBILE_TX_LOCS, const impalib_type ALPHA, const bool FILTERING_FLAG, const bool EXCLUDE_CAP_FLAG)
    :   numFixedTx_(NUM_FIXED_TX), numMobileTx_(NUM_MOBILE_TX), numBands_(NUM_BANDS), numTimeSteps_(NUM_TIME_STEPS),
        numRxLocs_(NUM_RX_LOCS), numMobileTxLocs_(NUM_MOBILE_TX_LOCS),
        alpha_(ALPHA), filteringFlag_(FILTERING_FLAG), excludeCapFlag_(EXCLUDE_CAP_FLAG),
        numProjects_(0), numTeams_(0), numDepartments_(0),
        numNodes_(0), numEdgeVariables_(0), numVariables_(0), numConstraints_(0), kvariable_(0){};

    void team_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<impalib_type> &,
                                           vector<impalib_type> &) const; ///< calculate messages from team equality constraint to ORIC

    static void project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &,
                                              vector<vector<impalib_type>> &); ///< calculate messages from project equality constraint to ORIC

    void edge_ec_to_degree_constraint_relaxed_graph_update(const vector<vector<int>> &, vector<vector<impalib_type>> &,
                                                           const vector<vector<impalib_type>> &,
                                                           vector<vector<impalib_type>> &) const; ///< calculate messages from edge to degree constraints for relaxed TSP
    
    static void flip_matrix(const vector<vector<impalib_type>> &, const vector<vector<int>> &, vector<vector<impalib_type>> &); ///< flip matrix
    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_update(const vector<vector<int>> &, const vector<impalib_type> &,
                                                                       const vector<vector<impalib_type>> &,
                                                                       const vector<vector<impalib_type>> &,
                                                                       const vector<vector<int>> &) const; ///< calculate message from edge to subtour constraint
    void                         edge_ec_to_degree_constraint_augmented_graph_update(const vector<vector<impalib_type>> &,
                                                                                     const vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                                                                     const vector<vector<impalib_type>> &,
                                                                                     vector<vector<impalib_type>> &) const; ///< calculate messages from edge to degree constraints for augmented TSP

    void variable_ec_to_ksat_constraint_update(const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<int> &, const vector<impalib_type> &, const vector<vector<int>> &) const;

    void x_eq_const_to_auxiliary_and_set_cover_const_update(vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &, vector<vector<impalib_type>> &, 
                                                            vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &,
                                                            vector<vector<int>> &, vector<vector<impalib_type>> &, vector<vector<int>> &,
                                                            vector<vector<impalib_type>> & ) const;
    
    vector<vector<impalib_type>> transpose_reshape(vector<vector<vector<impalib_type>>>) const;

    void z_eq_const_to_set_cover_ineq_const_update(vector<vector<impalib_type>>&, vector<vector<vector<vector<impalib_type>>>>&, vector<vector<impalib_type>>&,
                                                vector<vector<vector<vector<int>>>>&, vector<vector<impalib_type>>&) const;


    void r_eq_const_activation(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &,
                            vector<vector<impalib_type>> &) const;

    void z_eq_const_to_auxiliary_const_update(vector<vector<vector<vector<impalib_type>>>> &, vector<vector<vector<vector<int>>>> &, vector<vector<int>> &, vector<vector<impalib_type>> &,
                                            vector<vector<impalib_type>> &) const;

    void x_eq_const_activation(vector<vector<impalib_type>>&, vector<vector<vector<impalib_type>>>& ,
                        vector<vector<vector<impalib_type>>>&, vector<vector<vector<impalib_type>>>&,
                        vector<vector<int>>&, vector<vector<int>>&, vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &) const;                                        
};

/**
 * Calculate messages from team equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] rExtrinsicOutputDepartment: messages from departments to teams
 * @param[out] rTeam2OricM: messages from team equality constraint to ORIC
 * @param[in] rewardTeam: rewards of teams
 */

inline void EqualityConstraint::team_eq_constraint_to_oric_update(
    vector<vector<impalib_type>> &rExtrinsicOutputDepartment, vector<impalib_type> &rTeam2OricM,
    vector<impalib_type> &rewardTeam) const
{
    vector<impalib_type> intermediate_team_to_oric_m(numTeams_, 0);

    for (int department_index = 0; department_index < rExtrinsicOutputDepartment.size(); department_index++)
    {
        transform(rExtrinsicOutputDepartment[department_index].begin(),
                  rExtrinsicOutputDepartment[department_index].end(), intermediate_team_to_oric_m.begin(),
                  intermediate_team_to_oric_m.begin(), std::plus<impalib_type>());
    }
    transform(intermediate_team_to_oric_m.begin(), intermediate_team_to_oric_m.end(), rewardTeam.begin(),
              rTeam2OricM.begin(), std::plus<impalib_type>());
}

/**
 * Calculate messages from project equality constraint to ORIC for the Knapsack-MWM problem
 *
 * @param[in] rProject2EqConstraintM: messages from projects inequality constraints to project equality constraint
 * @param[out] rEqConstraint2OricM: messages from project equality constraints to ORIC
 * @param[in] rewardProject: rewards of teams-projects combinations
 */

inline void EqualityConstraint::project_eq_constraint_to_oric_update(vector<vector<impalib_type>> &rProject2EqConstraintM,
                                                                   vector<vector<impalib_type>> &rEqConstraint2OricM,
                                                                   vector<vector<impalib_type>> &rewardProject)
{
    for (int project_index = 0; project_index < rProject2EqConstraintM.size(); project_index++)
    {
        transform(rProject2EqConstraintM[project_index].begin(), rProject2EqConstraintM[project_index].end(),
                  rewardProject[project_index].begin(), rEqConstraint2OricM[project_index].begin(),
                  std::plus<impalib_type>());
    }
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the relaxed TSP
 *
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[in] rEdgeDegreeConstraintCost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to equality constraints
 * @param[out] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * 
 */

inline void EqualityConstraint::edge_ec_to_degree_constraint_relaxed_graph_update(
    const vector<vector<int>> &rEdgeConnections, vector<vector<impalib_type>> &rEdgeDegreeConstraintCost,
    const vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM) const
{

    vector<vector<impalib_type>> flipped_degree_constraint_to_eq_constraint_m = rDegreeConstraint2EqConstraintM;
    flip_matrix(rDegreeConstraint2EqConstraintM, rEdgeConnections, flipped_degree_constraint_to_eq_constraint_m);

    for (int edge_variable_index = 0; edge_variable_index < numEdgeVariables_; edge_variable_index++)
    {
        transform(flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].begin(),
                  flipped_degree_constraint_to_eq_constraint_m[edge_variable_index].end(),
                  rEdgeDegreeConstraintCost[edge_variable_index].begin(),
                  rEdgeEc2DegreeConstraintM[edge_variable_index].begin(), plus<impalib_type>());
    }
}


/**
 * Calculate messages from edge equality constraints to subtour elimination constraints for the TSP
 *
 * @param[in] rDeltaSIndicesList: list of edge indices forming the subtour elimination constraints
 * @param[in] rCostEdgeVaribale: costs of each edge variable
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraint to equality constraint
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour constraints to edge equality constraints
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @return edge_ec_to_subtour_constraints_list: messages from edge equality constraint to subtour elimination constraints
 * 
 */

inline vector<vector<impalib_type>> EqualityConstraint::edge_ec_to_subtour_constraints_update(
    const vector<vector<int>> &rDeltaSIndicesList, const vector<impalib_type> &rCostEdgeVaribale,
    const vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    const vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, const vector<vector<int>> &rEdgeConnections) const
{

    vector<vector<impalib_type>> edge_ec_to_subtour_constraints_list;

    if (rDeltaSIndicesList.size() == 1)
    {
        vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
        }

        for (size_t i = 0; i < rDeltaSIndicesList[0].size(); i++)
        {
            edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[0][i]] =
                combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[0][i]]
                + rCostEdgeVaribale[rDeltaSIndicesList[0][i]];
        }
        edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
    }
    
    else
    {
        vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
        for (const auto &row : rSubtourConstraints2EdgeEcM)
        {
            transform(combined_subtour_constraints_to_edge_ec_m.begin(),
                      combined_subtour_constraints_to_edge_ec_m.end(), row.begin(),
                      combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
        }

        vector<impalib_type> combined_degree_constraint_to_eq_constraint_m(numEdgeVariables_, zero_value);
        for (size_t index_edge_variable = 0; index_edge_variable < numEdgeVariables_; ++index_edge_variable)
        {
            combined_degree_constraint_to_eq_constraint_m[index_edge_variable] =
                rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][0]]
                + rDegreeConstraint2EqConstraintM[index_edge_variable][rEdgeConnections[index_edge_variable][1]];
        }

        for (size_t index_subtour_constraint = 0; index_subtour_constraint < rDeltaSIndicesList.size();
             index_subtour_constraint++)
        {
            vector<impalib_type> edge_ec_to_subtour_constraints_m(numEdgeVariables_, zero_value);

            for (size_t i = 0; i < rDeltaSIndicesList[index_subtour_constraint].size(); i++)
            {
                edge_ec_to_subtour_constraints_m[rDeltaSIndicesList[index_subtour_constraint][i]] =
                    combined_subtour_constraints_to_edge_ec_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + combined_degree_constraint_to_eq_constraint_m[rDeltaSIndicesList[index_subtour_constraint][i]]
                    + rCostEdgeVaribale[rDeltaSIndicesList[index_subtour_constraint][i]]
                    - rSubtourConstraints2EdgeEcM[index_subtour_constraint]
                                                 [rDeltaSIndicesList[index_subtour_constraint][i]];
            }
            edge_ec_to_subtour_constraints_list.push_back(edge_ec_to_subtour_constraints_m);
        }
    }
    return edge_ec_to_subtour_constraints_list;
}

/**
 * Calculate messages from edge equality constraints to degree constraints for the augmented TSP
 *
 * @param[in] rDegreeConstraint2EqConstraintM: messages from degree constraints to equality constraints
 * @param[in] rSubtourConstraints2EdgeEcM: messages from subtour elimination constraints to edge equality constraints
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[in] rEdgeDegreeConstraintCost: cost matrix of edges that has size function of number of edges and number of nodes
 * @param[out] rEdgeEc2DegreeConstraintM: messages from edge equality constraints to degree constraints
 * 
 */

inline void EqualityConstraint::edge_ec_to_degree_constraint_augmented_graph_update(
    const vector<vector<impalib_type>> &rDegreeConstraint2EqConstraintM,
    const vector<vector<impalib_type>> &rSubtourConstraints2EdgeEcM, const vector<vector<int>> &rEdgeConnections,
    const vector<vector<impalib_type>> &rEdgeDegreeConstraintCost, vector<vector<impalib_type>> &rEdgeEc2DegreeConstraintM) const
{

    vector<impalib_type> combined_subtour_constraints_to_edge_ec_m(numEdgeVariables_, zero_value);
    for (const auto &row : rSubtourConstraints2EdgeEcM)
    {
        transform(combined_subtour_constraints_to_edge_ec_m.begin(), combined_subtour_constraints_to_edge_ec_m.end(),
                  row.begin(), combined_subtour_constraints_to_edge_ec_m.begin(), std::plus<impalib_type>());
    }

    for (size_t i = 0; i < rEdgeConnections.size(); i++)
    {

        // Update the message for the first node of the edge
        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][0]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][1]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][0]];

        // Update the message for the second node of the edge
        rEdgeEc2DegreeConstraintM[i][rEdgeConnections[i][1]] =
            combined_subtour_constraints_to_edge_ec_m[i] + rDegreeConstraint2EqConstraintM[i][rEdgeConnections[i][0]]
            + rEdgeDegreeConstraintCost[i][rEdgeConnections[i][1]];
    }
}

/**
 * Flip a matrix to facilitate message updates for the TSP. Matrices will be in the same format during IMPA
 *
 * @param[in] rMatrix: matrix that requires flipping
 * @param[in] rEdgeConnections: list of connections for each edge equality constraint
 * @param[out] rFlippedMatrix: flipped matrix
 * 
 */

inline void EqualityConstraint::flip_matrix(const vector<vector<impalib_type>> &rMatrix, const vector<vector<int>> &rEdgeConnections,
                                        vector<vector<impalib_type>> &rFlippedMatrix)
{
    // Iterate over each edge connection
    for (size_t l = 0; l < rEdgeConnections.size(); ++l)
    {
        // Row index
        int row                = rEdgeConnections[l][0];
        // Column index
        int col                = rEdgeConnections[l][1];
        rFlippedMatrix[l][row] = rMatrix[l][col];
        rFlippedMatrix[l][col] = rMatrix[l][row];
    }
}

/**
 * Calculate messages from variable equality constraints to k-sat constraints for the K-SAT problem
 *
 * @param[in] rKsatConstraint2EqConstraintM_: messages from k-sat constraints to variable equality constraints
 * @param[out] rVariableEc2KsatConstraintM: messages variable equality constraints to from k-sat constraints
 * @param[in] rUsedVariables: used variables in creating the constraints
 * @param[in] rVariablesConnections: constraints connections for each variable
 *
 */

inline void EqualityConstraint::variable_ec_to_ksat_constraint_update(const vector<vector<impalib_type>> &rKsatConstraint2EqConstraintM_, vector<vector<impalib_type>> &rVariableEc2KsatConstraintM, vector<int> &rUsedVariables, const vector<impalib_type> &rIncomingMetricsCost, const vector<vector<int>> &rVariablesConnections) const
{

    for(auto& row : rVariableEc2KsatConstraintM) {
        row.assign(row.size(), zero_value);
    }

    vector<impalib_type> used_incoming_metrics_cost(numVariables_, zero_value);
    
    for_each(rUsedVariables.begin(), rUsedVariables.end(), [&](int n) {
        used_incoming_metrics_cost[n] = rIncomingMetricsCost[n];
    });

    vector<impalib_type> sum_messages(numVariables_, zero_value);

    for (int i = 0; i < numVariables_; ++i) {
        for (int j = 0; j < rVariablesConnections[i].size(); ++j) {
            sum_messages[i] += rKsatConstraint2EqConstraintM_[rVariablesConnections[i][j]][i];
        }
        sum_messages[i] +=used_incoming_metrics_cost[i];
    }

    for (int index_variable = 0; index_variable < rUsedVariables.size(); ++index_variable) {
        int variable = rUsedVariables[index_variable];
        for (int i = 0; i < rVariablesConnections[variable].size(); ++i) {
            int constraint = rVariablesConnections[variable][i];
            // This if statement check was added to account for the fact that a constraint can have the same variable more than once,
            // like in the benchmarks datasets. However, in practical cases, a variable cannot appear more than once in a constraint
            // and thus this if statement check can be dropped
            if (abs(rVariableEc2KsatConstraintM[constraint][variable])<abs(sum_messages[variable] - rKsatConstraint2EqConstraintM_[constraint][variable])){
                rVariableEc2KsatConstraintM[constraint][variable] = sum_messages[variable] - rKsatConstraint2EqConstraintM_[constraint][variable];
            }
        }
    }
}


inline void EqualityConstraint::x_eq_const_to_auxiliary_and_set_cover_const_update(vector<vector<vector<impalib_type>>> &rFixedCapacConst2FixedXEqConstM, vector<vector<vector<impalib_type>>> &rMobileCapacConst2MobileXEqConstM,
                                                         vector<vector<impalib_type>> & rAuxiliaryConst2MobileXEqConstM, vector<vector<vector<impalib_type>>> & rSetCoverIneqConst2FixedXEqConstM, 
                                                         vector<vector<vector<impalib_type>>> & rFixedTxCosts, vector<vector<vector<impalib_type>>> &rMobileTxCosts,
                                                         vector<vector<int>> & rConxMobTxPerNumMobTxLocs, vector<vector<impalib_type>> & rMobileXEqConst2AuxiliaryConstM, vector<vector<int>> &rConxFixedTxPerNumRXLocs,
                                                         vector<vector<impalib_type>> & rFixedXEqConst2SetCoverConstM) const

{
    vector<impalib_type> sums_mobile_x_eq_const(numMobileTx_*numBands_*numTimeSteps_, 0);
    vector<int> sums_conx_per_row(rConxMobTxPerNumMobTxLocs.size(), 0);
    vector<int> conx_rows;

    vector<impalib_type> flattened_mobile_tx_costs;

    std::vector<impalib_type> flattened_mobile_capac_const_2_mobile_x_eq_const;

    for (int i=0; i<numMobileTx_; i++){
        for (int j=0; j<numBands_; j++){
                flattened_mobile_tx_costs.insert(flattened_mobile_tx_costs.end(), rMobileTxCosts[i][j].begin(), rMobileTxCosts[i][j].end());
                flattened_mobile_capac_const_2_mobile_x_eq_const.insert(flattened_mobile_capac_const_2_mobile_x_eq_const.end(), rMobileCapacConst2MobileXEqConstM[i][j].begin(), rMobileCapacConst2MobileXEqConstM[i][j].end());
        }
    }

    for (int i=0; i<rConxMobTxPerNumMobTxLocs.size(); i++){
        sums_conx_per_row[i] = accumulate(rConxMobTxPerNumMobTxLocs[i].begin(), rConxMobTxPerNumMobTxLocs[i].end(), 0);
        if (sums_conx_per_row[i] !=0) {conx_rows.push_back(i);}
    }

    for (int i=0; i<conx_rows.size(); i++){
        
        impalib_type sum_elements = 0;
        
        for (size_t l = 0; l < rAuxiliaryConst2MobileXEqConstM[conx_rows[i]].size(); l++) {
            if (rConxMobTxPerNumMobTxLocs[conx_rows[i]][l] == 1) {
                sum_elements += rAuxiliaryConst2MobileXEqConstM[conx_rows[i]][l];
            }
        }

        if (! excludeCapFlag_){
            sums_mobile_x_eq_const[conx_rows[i]] = sum_elements + flattened_mobile_capac_const_2_mobile_x_eq_const[conx_rows[i]] + flattened_mobile_tx_costs[conx_rows[i]];
        }
        else
        {
            sums_mobile_x_eq_const[conx_rows[i]] = sum_elements + flattened_mobile_tx_costs[conx_rows[i]];
        }
        
    }

    for (int i=0; i<conx_rows.size(); i++){
        for (int n = 0; n< numMobileTxLocs_; n++){
            if (rConxMobTxPerNumMobTxLocs[conx_rows[i]][n] == 1){
                rMobileXEqConst2AuxiliaryConstM[conx_rows[i]][n] = sums_mobile_x_eq_const[conx_rows[i]] - rAuxiliaryConst2MobileXEqConstM[conx_rows[i]][n];
            }
    }

    }

    vector<impalib_type> sums_fixed_x_eq_const(numFixedTx_*numBands_*numTimeSteps_, 0);

    vector<int> sums_fixed_conx_per_row(rConxFixedTxPerNumRXLocs.size(), 0);
    vector<int> fixed_conx_rows;

    vector<impalib_type> flattened_fixed_tx_costs;

    std::vector<impalib_type> flattened_fixed_capac_const_2_fixed_x_eq_const;

    for (int i=0; i<numFixedTx_; i++){
        for (int j=0; j<numBands_; j++){
            flattened_fixed_tx_costs.insert(flattened_fixed_tx_costs.end(), rFixedTxCosts[i][j].begin(), rFixedTxCosts[i][j].end());
            flattened_fixed_capac_const_2_fixed_x_eq_const.insert(flattened_fixed_capac_const_2_fixed_x_eq_const.end(), rFixedCapacConst2FixedXEqConstM[i][j].begin(), rFixedCapacConst2FixedXEqConstM[i][j].end());
        }
    }

    for (int i=0; i<rConxFixedTxPerNumRXLocs.size(); i++){
        sums_fixed_conx_per_row[i] = accumulate(rConxFixedTxPerNumRXLocs[i].begin(), rConxFixedTxPerNumRXLocs[i].end(), 0);
        if (sums_fixed_conx_per_row[i] !=0) {fixed_conx_rows.push_back(i);}
    }

    auto reshaped_set_cover_ineq_const_to_fixed_x_eq_const = transpose_reshape(rSetCoverIneqConst2FixedXEqConstM);

    for (int i=0; i<fixed_conx_rows.size(); i++){
        
        impalib_type fixed_sum_elements = 0;
        for (size_t l = 0; l < reshaped_set_cover_ineq_const_to_fixed_x_eq_const[fixed_conx_rows[i]].size(); l++) {
            if (rConxFixedTxPerNumRXLocs[fixed_conx_rows[i]][l] == 1) {
                fixed_sum_elements += reshaped_set_cover_ineq_const_to_fixed_x_eq_const[fixed_conx_rows[i]][l];
            }
        }

        if (! excludeCapFlag_){
            sums_fixed_x_eq_const[fixed_conx_rows[i]] = fixed_sum_elements + flattened_fixed_capac_const_2_fixed_x_eq_const[fixed_conx_rows[i]] + flattened_fixed_tx_costs[fixed_conx_rows[i]];
        }
        else
        {
            sums_fixed_x_eq_const[fixed_conx_rows[i]] = fixed_sum_elements + flattened_fixed_tx_costs[fixed_conx_rows[i]];
        }
        
    }

    for (int i=0; i<fixed_conx_rows.size(); i++){
        for (int l = 0; l< numRxLocs_; l++){
            if (rConxFixedTxPerNumRXLocs[fixed_conx_rows[i]][l] == 1){
                rFixedXEqConst2SetCoverConstM[fixed_conx_rows[i]][l] = sums_fixed_x_eq_const[fixed_conx_rows[i]] - reshaped_set_cover_ineq_const_to_fixed_x_eq_const[fixed_conx_rows[i]][l];
            }
    }

    }

}


inline vector<vector<impalib_type>> EqualityConstraint::transpose_reshape(vector<vector<vector<impalib_type>>> rSetCoverIneqConst2FixedXEqConstM) const
{

    vector<vector<vector<impalib_type>>> transposed_matrix(numFixedTx_*numBands_, vector<vector<impalib_type>>(numTimeSteps_, vector<impalib_type>(numRxLocs_, 0)));

    vector<vector<impalib_type>> reshaped_matrix;

    for (int k = 0; k < numFixedTx_*numBands_; ++k) {
        for (int i = 0; i < numTimeSteps_; ++i) {
            vector<impalib_type> row;
            row.reserve(numRxLocs_);
            for (int j = 0; j < numRxLocs_; ++j) {
                row.push_back(rSetCoverIneqConst2FixedXEqConstM[i][j][k]);
            }
            reshaped_matrix.push_back(row);
        }
    }
    
    return reshaped_matrix;
}

inline void EqualityConstraint::z_eq_const_to_set_cover_ineq_const_update(vector<vector<impalib_type>>& rAuxiliaryConst2ZEqConstM, 
                                vector<vector<vector<vector<impalib_type>>>>& rSetCoverIneqConst2ZEqConstM, vector<vector<impalib_type>>& rZEqConst2SetCoverIneqConstM,
                                vector<vector<vector<vector<int>>>>& rTransposedConxMobTxRx, vector<vector<impalib_type>>& rZcosts) const{

        vector<vector<impalib_type>> reshaped_set_cover_ineq_const_to_z_eq_const_m;
        vector<vector<int>> reshaped_connectivity_mobile_tx_rx;

        for (size_t l=0; l<numRxLocs_; l++){
            vector<impalib_type> flattened;
            vector<int> flattened_conx;
            for (size_t j_i =0; j_i< numMobileTx_*numBands_; j_i++){
                for (size_t k=0; k<numTimeSteps_; k++){
                    for (size_t n=0; n<numMobileTxLocs_; n++){
                        flattened.push_back(rSetCoverIneqConst2ZEqConstM[k][l][n][j_i]);
                        flattened_conx.push_back(rTransposedConxMobTxRx[j_i][l][k][n]);
                    }
                }
            }
            reshaped_set_cover_ineq_const_to_z_eq_const_m.push_back(flattened); //transpose was not applied
            reshaped_connectivity_mobile_tx_rx.push_back(flattened_conx);
        }

    vector<vector<int>> transposed_reshaped_connectivity_mobile_tx_rx(numMobileTx_*numBands_*numTimeSteps_*numMobileTxLocs_, vector<int>(numRxLocs_, 0));

    for (size_t i = 0; i < numRxLocs_; ++i) {
        for (size_t j = 0; j < numMobileTx_*numBands_*numTimeSteps_*numMobileTxLocs_; ++j) {
            transposed_reshaped_connectivity_mobile_tx_rx[j][i] = reshaped_connectivity_mobile_tx_rx[i][j];
        }
    }

    vector<impalib_type> sums_z_eq_const(numMobileTx_*numBands_*numTimeSteps_*numMobileTxLocs_, 0);

    vector<int> sums_conx_per_row(numMobileTx_*numBands_*numTimeSteps_*numMobileTxLocs_, 0);
    vector<int> conx_rows;

    for (int i=0; i<sums_conx_per_row.size(); i++){
        sums_conx_per_row[i] = accumulate(transposed_reshaped_connectivity_mobile_tx_rx[i].begin(), transposed_reshaped_connectivity_mobile_tx_rx[i].end(), 0);
        if (sums_conx_per_row[i] !=0) {conx_rows.push_back(i);}
    }

    vector<impalib_type> flattened_auxiliary_const_to_z_eq_const_m;

    for (const auto& row : rAuxiliaryConst2ZEqConstM) {
        flattened_auxiliary_const_to_z_eq_const_m.insert(flattened_auxiliary_const_to_z_eq_const_m.end(), row.begin(), row.end());
    }

    vector<impalib_type> flattened_z_costs;

    for (const auto& row : rZcosts) {
        flattened_z_costs.insert(flattened_z_costs.end(), row.begin(), row.end());
    }

    for (int i=0; i<conx_rows.size(); i++){
        
        impalib_type sum_elements = 0;
        
        for (size_t l = 0; l < numRxLocs_; l++) {
            if (transposed_reshaped_connectivity_mobile_tx_rx[conx_rows[i]][l] == 1) {
                sum_elements += reshaped_set_cover_ineq_const_to_z_eq_const_m[l][conx_rows[i]];
            }
        }
        
        sums_z_eq_const[conx_rows[i]] = sum_elements + flattened_auxiliary_const_to_z_eq_const_m[conx_rows[i]] + flattened_z_costs[conx_rows[i]];
    }


    for (int i=0; i<conx_rows.size(); i++){
        for (int l = 0; l< numRxLocs_; l++){
            if (transposed_reshaped_connectivity_mobile_tx_rx[conx_rows[i]][l] == 1){
                rZEqConst2SetCoverIneqConstM[conx_rows[i]][l] = sums_z_eq_const[conx_rows[i]] - reshaped_set_cover_ineq_const_to_z_eq_const_m[l][conx_rows[i]];
            }
    }
    }

}

inline void EqualityConstraint::r_eq_const_activation(vector<vector<impalib_type>> & rAuxiliaryConst2REqConstM, vector<vector<impalib_type>> & rMobileLocEqConst2REqConstM, vector<vector<impalib_type>> & rREqConst2AuxiliaryConstM,
                                           vector<vector<impalib_type>> & rREqConst2MobileLocEqConstM, vector<vector<int>> &rConxMobTxR, vector<vector<impalib_type>> & rRCosts) const {

    vector<vector<impalib_type>> reshaped_auxiliary_const_to_r_eq_const_m(numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numBands_*numTimeSteps_, 0));

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t n = 0; n < numMobileTxLocs_; n++) {
            for (size_t j_k = 0; j_k < numBands_ * numTimeSteps_; j_k++) {
                reshaped_auxiliary_const_to_r_eq_const_m[i * numMobileTxLocs_ + n][j_k] = rAuxiliaryConst2REqConstM[i * numBands_ * numTimeSteps_ + j_k][n];
            }
        }
    }

    vector<impalib_type> sums_r_eq_const(numMobileTx_*numMobileTxLocs_, 0);

    vector<int> sums_conx_per_row(numMobileTx_*numMobileTxLocs_, 0);
    vector<int> conx_rows;

    for (int i=0; i<sums_conx_per_row.size(); i++){
        sums_conx_per_row[i] = accumulate(rConxMobTxR[i].begin(), rConxMobTxR[i].end(), 0);
        if (sums_conx_per_row[i] !=0) {conx_rows.push_back(i);}
    }

    vector<impalib_type> flattened_r_costs;

    for (const auto& row : rRCosts) {
        flattened_r_costs.insert(flattened_r_costs.end(), row.begin(), row.end());
    }

    vector<impalib_type> flattened_mobile_loc_eq_const_to_r_eq_const_m;

    for (const auto& row : rMobileLocEqConst2REqConstM) {
        flattened_mobile_loc_eq_const_to_r_eq_const_m.insert(flattened_mobile_loc_eq_const_to_r_eq_const_m.end(), row.begin(), row.end());
    }

    for (int i=0; i<conx_rows.size(); i++){
        
        impalib_type sum_elements = 0;
        
        for (size_t l = 0; l < numBands_*numTimeSteps_; l++) {
            if (rConxMobTxR[conx_rows[i]][l] == 1) {
                sum_elements += reshaped_auxiliary_const_to_r_eq_const_m[conx_rows[i]][l];
            }
        }
        
        sums_r_eq_const[conx_rows[i]] = sum_elements + flattened_mobile_loc_eq_const_to_r_eq_const_m[conx_rows[i]] + flattened_r_costs[conx_rows[i]];
    }


    for (int i=0; i<conx_rows.size(); i++){
        for (int l = 0; l< numBands_*numTimeSteps_; l++){
            if (rConxMobTxR[conx_rows[i]][l] == 1){
                rREqConst2AuxiliaryConstM[conx_rows[i]][l] = sums_r_eq_const[conx_rows[i]] - reshaped_auxiliary_const_to_r_eq_const_m[conx_rows[i]][l];
            }
    }
    }

    vector<impalib_type> flattened_r_eq_const_to_mobile_loc_eq_const_m(numMobileTx_*numMobileTxLocs_, 0);
    
    for (size_t i=0; i<numMobileTx_*numMobileTxLocs_; i++){
        impalib_type sum_elements = 0;
        for (size_t l = 0; l < numBands_*numTimeSteps_; l++) {
            if (rConxMobTxR[conx_rows[i]][l] == 1) {
                sum_elements += reshaped_auxiliary_const_to_r_eq_const_m[conx_rows[i]][l];
            }
        }
        flattened_r_eq_const_to_mobile_loc_eq_const_m[conx_rows[i]] = sum_elements + flattened_r_costs[conx_rows[i]];
    }

    for (size_t i = 0; i < numMobileTx_; i++) {
        for (size_t j = 0; j < numMobileTxLocs_; j++) {
            rREqConst2MobileLocEqConstM[i][j] = flattened_r_eq_const_to_mobile_loc_eq_const_m[i * numMobileTxLocs_ + j];
        }
    }

    }


inline void EqualityConstraint::z_eq_const_to_auxiliary_const_update(vector<vector<vector<vector<impalib_type>>>> & rSetCoverIneqConst2ZEqConstM, vector<vector<vector<vector<int>>>> & rConxMobTxRx, vector<vector<int>> & rConxMobTxPerNumMobTxLocs,
                                vector<vector<impalib_type>> & rZEqConst2AuxiliaryConstM, vector<vector<impalib_type>> &rZCosts) const {

    for (int j_i=0; j_i< numBands_*numMobileTx_; j_i++){
        for (int k=0; k< numTimeSteps_; k++){
            for (int n=0; n< numMobileTxLocs_; n++){
                impalib_type temp_sum = 0;
                for (int l=0; l< numRxLocs_; l++){
                    if (rConxMobTxRx[k][l][n][j_i] == 1){
                        temp_sum += rSetCoverIneqConst2ZEqConstM[k][l][n][j_i];
                    }
                } 
            if (rConxMobTxPerNumMobTxLocs[j_i*numTimeSteps_ + k][n] ==1){
                rZEqConst2AuxiliaryConstM[j_i*numTimeSteps_ + k][n] = temp_sum + rZCosts[j_i*numTimeSteps_ + k][n];
            }
            }
        }
    }
}


inline void EqualityConstraint::x_eq_const_activation(vector<vector<impalib_type>>& rAuxiliaryConst2MobileXEqConstM, vector<vector<vector<impalib_type>>>& rSetCoverIneqConst2FixedXEqConstM,
                    vector<vector<vector<impalib_type>>>& rMobileXEqConst2MobileCapacConstM, vector<vector<vector<impalib_type>>>& rFixedXEqConst2FixedCapacConstM,
                    vector<vector<int>>& rConxMobTxPerNumMobTxLocs, vector<vector<int>>& rConxFixedTxPerNumRXLocs, vector<vector<vector<impalib_type>>> &rMobileTxCosts,
                    vector<vector<vector<impalib_type>>> & rFixedTxCosts) const {

    vector<int> conx_rows_mobile;

    for (int i=0; i<numMobileTx_*numBands_*numTimeSteps_; i++) {
        int sum = accumulate(rConxMobTxPerNumMobTxLocs[i].begin(), rConxMobTxPerNumMobTxLocs[i].end(), 0);
        if (sum !=0){
            conx_rows_mobile.push_back(i);
        }
    }

    vector<impalib_type> flat_data_1(numMobileTx_*numBands_*numTimeSteps_, 0);

    vector<impalib_type> flattened_mobile_tx_costs;

    for (auto& matrix : rMobileTxCosts) {
        for (auto& row : matrix) {
            flattened_mobile_tx_costs.insert(flattened_mobile_tx_costs.end(), row.begin(), row.end());
        }
    }
    

    for (int i=0; i<conx_rows_mobile.size(); i++) {
        impalib_type row_sum = 0;
        for (int j = 0; j < numMobileTxLocs_; j++) {
            if (rConxMobTxPerNumMobTxLocs[conx_rows_mobile[i]][j] == 1) {
                row_sum += rAuxiliaryConst2MobileXEqConstM[conx_rows_mobile[i]][j];
            }
        }
        flat_data_1[conx_rows_mobile[i]] = row_sum + flattened_mobile_tx_costs[conx_rows_mobile[i]];
    }

    int flat_index_1 = 0;
    for (int i=0; i<numMobileTx_; i++){
        for (int j=0; j<numBands_; j++){
            for (int k=0; k<numTimeSteps_; k++){
                rMobileXEqConst2MobileCapacConstM[i][j][k] = flat_data_1[flat_index_1++];
            }
        }
    }

    auto reshaped_set_cover_ineq_const_to_fixed_x_eq_const = transpose_reshape(rSetCoverIneqConst2FixedXEqConstM);

    vector<int> conx_rows_fixed;

    for (int i=0; i<numFixedTx_*numBands_*numTimeSteps_; i++) {
        int sum = std::accumulate(rConxFixedTxPerNumRXLocs[i].begin(), rConxFixedTxPerNumRXLocs[i].end(), 0);
        if (sum !=0){
            conx_rows_fixed.push_back(i);
        }
    }

    vector<impalib_type> flat_data_2(numFixedTx_*numBands_*numTimeSteps_, 0);
    vector<impalib_type> flattened_fixed_tx_costs;

    for (auto& matrix : rFixedTxCosts) {
        for (auto& row : matrix) {
            flattened_fixed_tx_costs.insert(flattened_fixed_tx_costs.end(), row.begin(), row.end());
        }
    }

    for (int i=0; i<conx_rows_fixed.size(); i++) {
        impalib_type row_sum = 0;
        for (int j = 0; j < numRxLocs_; j++) {
            if (rConxFixedTxPerNumRXLocs[conx_rows_fixed[i]][j] == 1) {
                row_sum += reshaped_set_cover_ineq_const_to_fixed_x_eq_const[conx_rows_fixed[i]][j];
            }
        }
        flat_data_2[conx_rows_fixed[i]] = row_sum + flattened_fixed_tx_costs[conx_rows_fixed[i]];
    }

    int flat_index_2 = 0;
    for (int i=0; i<numFixedTx_; i++){
        for (int j=0; j<numBands_; j++){
            for (int k=0; k<numTimeSteps_; k++){
                rFixedXEqConst2FixedCapacConstM[i][j][k] = flat_data_2[flat_index_2++];
            }
        }
    }


}