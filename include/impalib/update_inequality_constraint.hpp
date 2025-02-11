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
class InequalityConstraintKcMwm
{
private:
    int numProjects_; ///< number of projects
    int numTeams_; ///< number of teams
    int numDepartments_; ///< number of departments
    int maxStateIc_ = 1; ///< maximum value of project inequality constraint (<=1)

public:
    void project_inequality_constraint_update(const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &) const; ///< calculate messages from project inequality constraint to project equality constraint
    InequalityConstraintKcMwm(int N_DEPARTMENTS, int N_TEAMS, int N_PROJECTS); ///< constructor
};

/**
 * Construct InequalityConstraint object for the Knapsack-MWM problem
 *
 * 
 */

inline InequalityConstraintKcMwm::InequalityConstraintKcMwm(const int N_DEPARTMENTS, const int N_TEAMS, const int N_PROJECTS)
    : numProjects_(N_PROJECTS), numTeams_(N_TEAMS), numDepartments_(N_DEPARTMENTS){
                                                    };

/**
 * Calculate messages from project inequality constraint to project equality constraint for the Knapsack-MWM problem
 *
 * @param[in] rEqConstraint2ProjectM: messages from project equality constraint to project inequality constraint
 * @param[out] rProject2EqConstraintM: messages from project inequality constraint to project equality constraint
 * 
 */

inline void InequalityConstraintKcMwm::project_inequality_constraint_update(const vector<vector<impalib_type>> &rEqConstraint2ProjectM,
                                                                vector<vector<impalib_type>> &rProject2EqConstraintM) const
{

    vector<vector<impalib_type>> stage_forward_messages_project_EC(numTeams_ + 1,
                                                                   vector<impalib_type>(maxStateIc_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages_project_EC(numTeams_ + 1,
                                                                    vector<impalib_type>(maxStateIc_ + 1, zero_value));

    for (int project_index = 0; project_index < rProject2EqConstraintM.size(); project_index++)
    {
        // Initialize forward messages
        vector<impalib_type> initial_forward_messages(maxStateIc_ + 1, zero_value),
            initial_backward_messages(maxStateIc_ + 1, zero_value);
        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages_project_EC[0] = initial_forward_messages;
        
        for (int stage = 0; stage < numTeams_; stage++)
        {
            stage_forward_messages_project_EC[stage + 1][0] = stage_forward_messages_project_EC[stage][0];
            stage_forward_messages_project_EC[stage + 1][1] =
                min(stage_forward_messages_project_EC[stage][1],
                    stage_forward_messages_project_EC[stage][0] + rEqConstraint2ProjectM[project_index][stage]);
        }

        stage_backward_messages_project_EC[numTeams_] = initial_backward_messages;

        for (int stage = numTeams_ - 1; stage >= 0; stage--)
        {
            stage_backward_messages_project_EC[stage][0] =
                min(stage_backward_messages_project_EC[stage + 1][0],
                    stage_backward_messages_project_EC[stage + 1][1] + rEqConstraint2ProjectM[project_index][stage]);
            stage_backward_messages_project_EC[stage][1] = stage_backward_messages_project_EC[stage + 1][1];
        }

        // Update project to equality constraint messages
        for (int team_index = 0; team_index < numTeams_; team_index++)
        {
            impalib_type minimumValue                         = zero_value;
            minimumValue                                      = min(stage_forward_messages_project_EC[team_index][1],
                                                                    stage_backward_messages_project_EC[team_index + 1][0]);
            rProject2EqConstraintM[project_index][team_index] = -min(minimumValue, zero_value);
        }
    }
}


/**
 * Represents a class for the inequality constraint for the MOBARP problem
 */
class InequalityConstraintMOBARP
{
private:
    int numFixedTx_;
    int numMobileTx_;
    int numBands_;
    int numTimeSteps_;
    int numRxLocs_;
    int numMobileTxLocs_;
    impalib_type alpha_;
    bool filteringFlag_;
    vector<vector<vector<impalib_type>>> FixedCapacConst2FixedXEqConstOld_;
    vector<vector<vector<impalib_type>>> MobileCapacConst2MobileXEqConstOld_;
    vector<vector<impalib_type>> MobileLocEqConst2REqConstOld_;
    vector<vector<vector<impalib_type>>> SetCoverIneqConst2FixedXEqConstOld_;
    vector<vector<vector<vector<impalib_type>>>> SetCoverIneqConst2ZEqConstOld_;
    impalib_type                 initial_forward_message_; ///< initial forward message of forward-backward algorithm
    impalib_type                 initial_backward_message_; ///< initial backward message of forward-backward algorithm
    int maxState_=1;

public:
    void ineq_capac_const_update(const vector<vector<vector<impalib_type>>> &, const vector<vector<vector<impalib_type>>> &, const vector<int> &, const vector<int> &, 
                                vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &) const;
    void process_filtering(int, vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &,
                            vector<vector<vector<impalib_type>>> &, vector<vector<vector<impalib_type>>> &);  ///< perform filtering
    
    void mobile_loc_eq_const_to_r_eq_const_update(vector<vector<impalib_type>> &,vector<vector<impalib_type>> &) const;

    void process_filtering_mobile_loc_eq(int , vector<vector<impalib_type>> &, vector<vector<impalib_type>> &);

    void set_cover_ineq_const_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<vector<vector<int>>>> &, vector<vector<int>> &, vector<vector<int>> &,
                                    vector<vector<vector<impalib_type>>> &, vector<vector<vector<vector<impalib_type>>>> &);

    void process_filtering_set_cover_const(int , vector<vector<vector<impalib_type>>> &, vector<vector<vector<vector<impalib_type>>>> &, vector<vector<vector<impalib_type>>> &, vector<vector<vector<vector<impalib_type>>>> &);

    InequalityConstraintMOBARP(int NUM_FIXED_TX, int NUM_MOBILE_TX, int NUM_BANDS, 
                                int NUM_TIME_STEPS, int NUM_RX_LOCS, int NUM_MOBILE_TX_LOCS, impalib_type ALPHA, bool FILTERING_FLAG); ///< constructor
};

inline InequalityConstraintMOBARP::InequalityConstraintMOBARP(const int NUM_FIXED_TX, const int NUM_MOBILE_TX, const int NUM_BANDS, 
                                                                const int NUM_TIME_STEPS, const int NUM_RX_LOCS, const int NUM_MOBILE_TX_LOCS, const impalib_type ALPHA, const bool FILTERING_FLAG)
    : numFixedTx_(NUM_FIXED_TX), numMobileTx_(NUM_MOBILE_TX), numBands_(NUM_BANDS), numTimeSteps_(NUM_TIME_STEPS), numRxLocs_(NUM_RX_LOCS), numMobileTxLocs_(NUM_MOBILE_TX_LOCS), 
        alpha_(ALPHA), filteringFlag_(FILTERING_FLAG), 
        FixedCapacConst2FixedXEqConstOld_(NUM_FIXED_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0))),
        MobileCapacConst2MobileXEqConstOld_(NUM_MOBILE_TX, vector<vector<impalib_type>>(NUM_BANDS, vector<impalib_type>(NUM_TIME_STEPS, 0))),
        MobileLocEqConst2REqConstOld_(NUM_MOBILE_TX, vector<impalib_type>(NUM_MOBILE_TX_LOCS, 0)),
        SetCoverIneqConst2FixedXEqConstOld_(NUM_TIME_STEPS, vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0))),
        SetCoverIneqConst2ZEqConstOld_(vector<vector<vector<vector<impalib_type>>>>(NUM_TIME_STEPS, vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))))),
        initial_forward_message_(value_inf), initial_backward_message_(value_inf){};


inline void InequalityConstraintMOBARP::ineq_capac_const_update(const vector<vector<vector<impalib_type>>> &rFixedXEqConst2FixedCapacConstM, const vector<vector<vector<impalib_type>>> & rMobileXEqConst2MobileCapacConstM, 
                                const vector<int> &rFixedCapacConstraints, const vector<int> & rMobileCapacConstraints, 
                                vector<vector<vector<impalib_type>>> & rFixedCapacConst2FixedXEqConstDummyM, vector<vector<vector<impalib_type>>> & rMobileCapacConst2MobileXEqConstDummyM) const
{

    // for (int i=0; i<numFixedTx_; i++){
    //     for (int k=0; k<numTimeSteps_; k++){

    //         for (int j=0; j<numBands_; j++){
    //            vector<impalib_type> fixed_result;
    //            fixed_result.reserve(numBands_-1);

    //            for (int indx_j = 0; indx_j< numBands_; indx_j++){
    //             if (indx_j != j){
    //                 fixed_result.push_back(rFixedXEqConst2FixedCapacConstM[i][indx_j][k]);
    //             }
    //            }
    //             sort(fixed_result.begin(), fixed_result.end());
    //             rFixedCapacConst2FixedXEqConstDummyM[i][j][k] = max(-fixed_result[rFixedCapacConstraints[i]-1], zero_value);
    //         }
    // }

    // }

    // cout<<"Fixed: \n";
    for (int i=0; i<numFixedTx_; i++){
        // cout<<"maxState_ Calculation:\n";
        int maxState_ = rFixedCapacConstraints[i];
        // cout<<"maxState_: "<<maxState_<<"\n";

        vector<vector<impalib_type>> stage_forward_messages(numBands_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
        vector<vector<impalib_type>> stage_backward_messages(numBands_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

        for (int k=0; k<numTimeSteps_; k++){
            // cout<<"i"<<i<<", k"<<k<<"\n";
            vector<impalib_type> initial_forward_messages( maxState_+ 1, zero_value);
            vector<impalib_type>  initial_backward_messages(maxState_ + 1, zero_value);
            fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);
            stage_forward_messages[0] = initial_forward_messages;
            stage_backward_messages[numBands_] = initial_backward_messages;

            for (int j=0; j< numBands_; j++){
                // cout<<"rFixedXEqConst2FixedCapacConstM[i][j][k]: "<<rFixedXEqConst2FixedCapacConstM[i][j][k]<<"\n";
                // for (int a=0; a<=min(maxState_, j+1); a++){
                for (int a=0; a<=maxState_; a++){
                    // if (a - 1>=0){
                        if (a-1>=0 && a<=min(maxState_, j+1)){
                        stage_forward_messages[j + 1][a] =
                            min(initial_forward_messages[a],
                                initial_forward_messages[a - 1]
                                    + rFixedXEqConst2FixedCapacConstM[i][j][k]);
                        }
                        else{
                            stage_forward_messages[j + 1][a] = initial_forward_messages[a];
                        }
                }

                initial_forward_messages = stage_forward_messages[j + 1];
            }

            // for (int index=0; index<numBands_; index++){
            //     for (int a=0; a<=maxState_; a++){
            //         cout<<stage_forward_messages[index][a]<<" ";
            //     }
            //     cout<<"\n";
            // }
            // cout<<"\n";

            for (int j = numBands_ - 1; j >= 0; j--)
            {
                for (int a = 0; a <= min(j, maxState_); a++)
                {
                    if (a == maxState_){
                        stage_backward_messages[j][a] = initial_backward_messages[a];
                    }
                    else{
                        stage_backward_messages[j][a] =
                            min(initial_backward_messages[a],
                                initial_backward_messages[a + 1]
                                    + rFixedXEqConst2FixedCapacConstM[i][j][k]);
                    }

                    }
                initial_backward_messages = stage_backward_messages[j];

            }

            // for (int index=0; index<numBands_; index++){
            //     for (int a=0; a<=maxState_; a++){
            //         cout<<stage_backward_messages[index][a]<<" ";
            //     }
            //     cout<<"\n";
            // }
            // cout<<"\n";

            vector<impalib_type> metric_path_solid, metric_path_dash;

            for (int j = 0; j < numBands_; j++)
            {
                // cout<<"j: "<<j<<"\n";
                metric_path_dash.clear();
                metric_path_solid.clear();

                if (j == 0)
                {
                    metric_path_solid.push_back(
                        stage_forward_messages[j][0]
                        + stage_backward_messages[j + 1][1]
                        + rFixedXEqConst2FixedCapacConstM[i][j][k]);
                    
                    metric_path_dash.push_back(stage_forward_messages[0][0]
                                            + stage_backward_messages[j + 1][0]);
                }
                else
                {
                    for (int a = 0; a <= min(j, maxState_-1); a++)
                    {
                        metric_path_solid.push_back(stage_forward_messages[j][a]+ stage_backward_messages[j + 1][a + 1] + rFixedXEqConst2FixedCapacConstM[i][j][k]);
                    }

                    for (int a = 0; a <= min(j,maxState_); a++)
                    {
                        metric_path_dash.push_back(stage_forward_messages[j][a] + stage_backward_messages[j + 1][a]);
                    }

                }

                // cout<<"j: "<<j<<"\n";
                // cout<<"metric_path_solid\n";
                // for (auto& e: metric_path_solid){
                //     cout<<e<<" ";
                // }
                // cout<<"\nmetric_path_dash\n";
                // for (auto& e: metric_path_dash){
                //     cout<<e<<" ";
                // }
                // cout<<"\n";
                rFixedCapacConst2FixedXEqConstDummyM[i][j][k]=
                    *min_element(metric_path_solid.begin(), metric_path_solid.end())
                    - *min_element(metric_path_dash.begin(), metric_path_dash.end())
                    - rFixedXEqConst2FixedCapacConstM[i][j][k];
                // cout<<"rFixedCapacConst2FixedXEqConstDummyM[i][j][k]: "<<rFixedCapacConst2FixedXEqConstDummyM[i][j][k]<<"\n";
            }
        }
    }


    // cout<<"Mobile: \n";
    for (int i=0; i<numMobileTx_; i++){
        // cout<<"maxState_ Calculation:\n";
        int maxState_ = rMobileCapacConstraints[i];
        // cout<<"maxState_: "<<maxState_<<"\n";

        vector<vector<impalib_type>> stage_forward_messages(numBands_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
        vector<vector<impalib_type>> stage_backward_messages(numBands_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

        for (int k=0; k<numTimeSteps_; k++){
            // cout<<"i"<<i<<", k"<<k<<"\n";
            
            vector<impalib_type> initial_forward_messages( maxState_+ 1, zero_value);
            vector<impalib_type>  initial_backward_messages(maxState_ + 1, zero_value);

            fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);
            stage_forward_messages[0] = initial_forward_messages;
            stage_backward_messages[numBands_] = initial_backward_messages;

            for (int j=0; j< numBands_; j++){
                for (int a=0; a<=maxState_; a++){
                    // if (a - 1>=0){
                        if (a-1>=0 && a<=min(maxState_, j+1)){
                        stage_forward_messages[j + 1][a] =
                            min(initial_forward_messages[a],
                                initial_forward_messages[a - 1]
                                    + rMobileXEqConst2MobileCapacConstM[i][j][k]);
                        }
                        else{
                            stage_forward_messages[j + 1][a] = initial_forward_messages[a];
                        }
                }

                initial_forward_messages = stage_forward_messages[j + 1];
            }

            for (int j = numBands_ - 1; j >= 0; j--)
            {
                for (int a = 0; a <= min(j, maxState_); a++)
                {
                    if (a==maxState_){
                        stage_backward_messages[j][a]  = initial_backward_messages[a];
                    }
                    else{
                        stage_backward_messages[j][a] =
                            min(initial_backward_messages[a],
                                initial_backward_messages[a + 1]
                                    + rMobileXEqConst2MobileCapacConstM[i][j][k]);
                    }

                    }
                initial_backward_messages = stage_backward_messages[j];

            }

            vector<impalib_type> metric_path_solid, metric_path_dash;

            for (int j = 0; j < numBands_; j++)
            {
                // cout<<"j: "<<j<<"\n";
                metric_path_dash.clear();
                metric_path_solid.clear();

                if (j == 0)
                {
                    metric_path_solid.push_back(
                        stage_forward_messages[j][0]
                        + stage_backward_messages[j + 1][1]
                        + rMobileXEqConst2MobileCapacConstM[i][j][k]);
                    
                    metric_path_dash.push_back(stage_forward_messages[0][0]
                                            + stage_backward_messages[j + 1][0]);
                }
                else
                {
                    for (int a = 0; a <= min(j, maxState_-1); a++)
                    {
                        metric_path_solid.push_back(stage_forward_messages[j][a]+ stage_backward_messages[j + 1][a + 1] + rMobileXEqConst2MobileCapacConstM[i][j][k]);
                    }

                    for (int a = 0; a <= min(j,maxState_); a++)
                    {
                        metric_path_dash.push_back(stage_forward_messages[j][a] + stage_backward_messages[j + 1][a]);
                    }

                }
                rMobileCapacConst2MobileXEqConstDummyM[i][j][k]=
                    *min_element(metric_path_solid.begin(), metric_path_solid.end())
                    - *min_element(metric_path_dash.begin(), metric_path_dash.end())
                    - rMobileXEqConst2MobileCapacConstM[i][j][k];
            }
        }
    }

    // for (int i=0; i<numMobileTx_; i++){
    //     for (int k=0; k<numTimeSteps_; k++){

    //         for (int j=0; j<numBands_; j++){
    //            vector<impalib_type> mobile_result;
    //            mobile_result.reserve(numBands_-1);

    //            for (int indx_j = 0; indx_j< numBands_; indx_j++){
    //             if (indx_j != j){
    //                 mobile_result.push_back(rMobileXEqConst2MobileCapacConstM[i][indx_j][k]);
    //             }
    //            }
    //             sort(mobile_result.begin(), mobile_result.end());
    //             rMobileCapacConst2MobileXEqConstDummyM[i][j][k] = max(-mobile_result[rMobileCapacConstraints[i]-1], zero_value);
    //         }
    // }
    // }

    // cout<<"Printing"<<"\n";
    // for (int i = 0; i < numFixedTx_; ++i) {
    //     for (int j = 0; j < numBands_; ++j) {
    //         for (int k = 0; k < numTimeSteps_; ++k) {
    //             cout << rFixedCapacConst2FixedXEqConstDummyM[i][j][k]<<" ";
    //         }
    //     }
    // }

    // cout<<"\n";
    
    // fstream file_output_1("./ut_results/rFixedCapacConst2FixedXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<numFixedTx_; i++){
    //     for (int j=0; j<numBands_; j++){
    //         for (int k=0; k<numTimeSteps_; k++){
    //         file_output_1.write((char*)(&rFixedCapacConst2FixedXEqConstDummyM[i][j][k]), sizeof(rFixedCapacConst2FixedXEqConstDummyM[i][j][k]));}}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

    // fstream file_output_2("./ut_results/rMobileCapacConst2MobileXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_2.is_open()) {
    //     for (int i=0; i<numMobileTx_; i++){
    //     for (int j=0; j<numBands_; j++){
    //         for (int k=0; k<numTimeSteps_; k++){
    //         file_output_2.write((char*)(&rMobileCapacConst2MobileXEqConstDummyM[i][j][k]), sizeof(rMobileCapacConst2MobileXEqConstDummyM[i][j][k]));}}}
    //         file_output_2.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

// exit(0);
}

inline void InequalityConstraintMOBARP::process_filtering(const int iter, vector<vector<vector<impalib_type>>> &rFixedCapacConst2FixedXEqConstDummyM, vector<vector<vector<impalib_type>>> &rMobileCapacConst2MobileXEqConstDummyM,
                                                            vector<vector<vector<impalib_type>>> &rFixedCapacConst2FixedXEqConstM, vector<vector<vector<impalib_type>>> &rMobileCapacConst2MobileXEqConstM) {
    
    for (int j=0; j<numBands_; j++){

        for (int i = 0; i < numFixedTx_; i++) {

            if ((filteringFlag_) and (alpha_ != zero_value)) {
                vector<impalib_type> intermediate_dummy(rFixedCapacConst2FixedXEqConstDummyM[i][j]), intermediate_old(FixedCapacConst2FixedXEqConstOld_[i][j]), intermediate_extrinsic;

                impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
                transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

                if (iter == 0) {
                    copy(intermediate_dummy.begin(), intermediate_dummy.end(), rFixedCapacConst2FixedXEqConstM[i][j].begin());
                } else {
                    transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                    copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rFixedCapacConst2FixedXEqConstM[i][j].begin());
                }
                copy(rFixedCapacConst2FixedXEqConstM[i][j].begin(), rFixedCapacConst2FixedXEqConstM[i][j].end(), FixedCapacConst2FixedXEqConstOld_[i][j].begin());
            }

            else {
                copy(rFixedCapacConst2FixedXEqConstDummyM[i][j].begin(), rFixedCapacConst2FixedXEqConstDummyM[i][j].end(), rFixedCapacConst2FixedXEqConstM[i][j].begin());
            }
        }

        for (int i = 0; i < numMobileTx_; i++) {

            if ((filteringFlag_) and (alpha_ != zero_value)) {
                vector<impalib_type> intermediate_dummy(rMobileCapacConst2MobileXEqConstDummyM[i][j]), intermediate_old(MobileCapacConst2MobileXEqConstOld_[i][j]), intermediate_extrinsic;

                impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
                transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

                if (iter == 0) {
                    copy(intermediate_dummy.begin(), intermediate_dummy.end(), rMobileCapacConst2MobileXEqConstM[i][j].begin());
                } else {
                    transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                    copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rMobileCapacConst2MobileXEqConstM[i][j].begin());
                }
                copy(rMobileCapacConst2MobileXEqConstM[i][j].begin(), rMobileCapacConst2MobileXEqConstM[i][j].end(), MobileCapacConst2MobileXEqConstOld_[i][j].begin());
            }

            else {
                copy(rMobileCapacConst2MobileXEqConstDummyM[i][j].begin(), rMobileCapacConst2MobileXEqConstDummyM[i][j].end(), rMobileCapacConst2MobileXEqConstM[i][j].begin());
            }
        }

    }

    // for (int i = 0; i < numFixedTx_; i++) {
    //     cout << "Layer " << i + 1 << ":\n";
    //     for (int j = 0; j < numBands_; j++) {
    //         for (int k = 0; k < numTimeSteps_; k++) {
    //             cout << rFixedCapacConst2FixedXEqConstM_[i][j][k] << " ";
    //         }
    //         cout << "\n";
    //     }
    //     cout << "\n";
    // }

    // for (int i = 0; i < numMobileTx_; i++) {
    //     cout << "Layer " << i + 1 << ":\n";
    //     for (int j = 0; j < numBands_; j++) {
    //         for (int k = 0; k < numTimeSteps_; k++) {
    //             cout << rMobileCapacConst2MobileXEqConstM_[i][j][k] << " ";
    //         }
    //         cout << "\n";
    //     }
    //     cout << "\n";
    // }
    // exit(0);

    // fstream file_output_1("./ut_results/rFixedCapacConst2FixedXEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<numFixedTx_; i++){
    //     for (int j=0; j<numBands_; j++){
    //         for (int k=0; k<numTimeSteps_; k++){
    //         file_output_1.write((char*)(&rFixedCapacConst2FixedXEqConstM[i][j][k]), sizeof(rFixedCapacConst2FixedXEqConstM[i][j][k]));}}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

    // fstream file_output_2("./ut_results/rMobileCapacConst2MobileXEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_2.is_open()) {
    //     for (int i=0; i<numMobileTx_; i++){
    //     for (int j=0; j<numBands_; j++){
    //         for (int k=0; k<numTimeSteps_; k++){
    //         file_output_2.write((char*)(&rMobileCapacConst2MobileXEqConstM[i][j][k]), sizeof(rMobileCapacConst2MobileXEqConstM[i][j][k]));}}}
    //         file_output_2.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}

}

inline void InequalityConstraintMOBARP::mobile_loc_eq_const_to_r_eq_const_update(vector<vector<impalib_type>> & rREqConst2MobileLocEqConstM,vector<vector<impalib_type>> & rMobileLocEqConst2REqConstDummyM) const {

    //need to do it using FBA
    // auto start_1 = chrono::high_resolution_clock::now();
    // for (size_t n = 0; n < numMobileTxLocs_; n++) {
        
    //     for (size_t i = 0; i < numMobileTx_; i++) {
    //         impalib_type min_value = value_inf;
    //         for (size_t col = 0; col < numMobileTxLocs_; ++col) {
    //             if (col != n) {
    //                 min_value = min(min_value, rREqConst2MobileLocEqConstM[i][col]);
    //             }
    //         }
    //         rMobileLocEqConst2REqConstDummyM[i][n] = -min_value;
    //     }
    // }
    // auto end_1 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_1 = end_1 - start_1;
    // cout << "Execution time unoptimized: " << elapsed_1.count() << " seconds\n";

    // auto start_2 = chrono::high_resolution_clock::now();
    // vector<vector<impalib_type>> MobileLocEqConst2REqConstDummyMNew(numMobileTx_, vector<impalib_type>(numMobileTxLocs_, zero_value));

    vector<vector<impalib_type>> stage_forward_messages(numMobileTxLocs_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages(numMobileTxLocs_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

    for (size_t i=0; i< numMobileTx_; i++){
            vector<impalib_type> initial_forward_messages(maxState_ + 1, zero_value), initial_backward_messages(maxState_ + 1, zero_value);
            fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

            stage_forward_messages[0] = initial_forward_messages;
            
            for (int stage = 0; stage < numMobileTxLocs_; stage++)
            {
                stage_forward_messages[stage + 1][0] = stage_forward_messages[stage][0];
                stage_forward_messages[stage + 1][1] = min(stage_forward_messages[stage][1], stage_forward_messages[stage][0] + rREqConst2MobileLocEqConstM[i][stage]);
            }

            stage_backward_messages[numMobileTxLocs_] = initial_backward_messages;

            for (int stage = numMobileTxLocs_ - 1; stage >= 0; stage--)
            {
                if (stage == numMobileTxLocs_ - 1){
                    stage_backward_messages[stage][0] = rREqConst2MobileLocEqConstM[i][stage];
                }
                else{
                    stage_backward_messages[stage][0] = min(stage_backward_messages[stage + 1][0], stage_backward_messages[stage + 1][1] + rREqConst2MobileLocEqConstM[i][stage]);
                }
                stage_backward_messages[stage][1] = stage_backward_messages[stage + 1][1];
            }

        impalib_type min_dashed_edges = zero_value;
        impalib_type min_solid_edges = zero_value;

        for (int n = 0; n < numMobileTxLocs_; n++)
        {
            min_solid_edges = stage_forward_messages[n][0] + stage_backward_messages[n+1][1] + rREqConst2MobileLocEqConstM[i][n];
            
            if (n == numMobileTxLocs_-1){
                min_dashed_edges = stage_forward_messages[n][1] + stage_backward_messages[n+1][1];
            }
            else {
                min_dashed_edges = min(stage_forward_messages[n][0] + stage_backward_messages[n+1][0], stage_forward_messages[n][1] + stage_backward_messages[n+1][1]);
            }

            rMobileLocEqConst2REqConstDummyM[i][n] = (min_solid_edges  - min_dashed_edges - rREqConst2MobileLocEqConstM[i][n]);
        }

        }

    // auto end_2 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_2 = end_2 - start_2;
    // cout << "Execution time optimized: " << elapsed_2.count() << " seconds\n";

    // for (int i = 0; i < numMobileTx_; ++i) {
    //     for (int n = 0; n < numMobileTxLocs_; ++n) {
    //         if (std::abs(MobileLocEqConst2REqConstDummyMNew[i][n] - rMobileLocEqConst2REqConstDummyM[i][n]) > 1e-7) {
    //             cout << "Error: Values are not equal at index [" << i << "][" << n << "]\n";
    //             cout<<"MobileLocEqConst2REqConstDummyMNew[i][n]: "<<MobileLocEqConst2REqConstDummyMNew[i][n]<<"\n";
    //             cout<<"rMobileLocEqConst2REqConstDummyM[i][n]: "<<rMobileLocEqConst2REqConstDummyM[i][n]<<"\n";
    //             exit(EXIT_FAILURE);
    //         }
    //     }
    // }


    // fstream file_output_1("./ut_results/rMobileLocEqConst2REqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<rMobileLocEqConst2REqConstDummyM.size(); i++){
    //     for (int j=0; j<rMobileLocEqConst2REqConstDummyM[0].size(); j++){
    //         file_output_1.write((char*)(&rMobileLocEqConst2REqConstDummyM[i][j]), sizeof(rMobileLocEqConst2REqConstDummyM[i][j]));}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}
    
    // cout << "rMobileLocEqConst2REqConstDummyM_:" << "\n";
    // for (size_t i = 0; i < numMobileTx_; ++i) {
    //     for (size_t n = 0; n < numMobileTxLocs_; ++n) {
    //         cout << setw(8) << rMobileLocEqConst2REqConstDummyM_[i][n] << " "; 
    //     }
    //     cout << "\n"; // Newline after each row
    // }

}

inline void InequalityConstraintMOBARP::process_filtering_mobile_loc_eq(int iter, vector<vector<impalib_type>> & rMobileLocEqConst2REqConstDummyM, vector<vector<impalib_type>> & rMobileLocEqConst2REqConstM) {

 for (int i = 0; i < numMobileTx_; i++)
    {
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            // Calculate weighted values for current and old messages
            vector<impalib_type> intermediate_dummy(rMobileLocEqConst2REqConstDummyM[i]),intermediate_old(MobileLocEqConst2REqConstOld_[i]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(), rMobileLocEqConst2REqConstM[i].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rMobileLocEqConst2REqConstM[i].begin());
            }
            copy(rMobileLocEqConst2REqConstM[i].begin(), rMobileLocEqConst2REqConstM[i].end(), MobileLocEqConst2REqConstOld_[i].begin());
        }

        else
        {
            copy(rMobileLocEqConst2REqConstDummyM[i].begin(), rMobileLocEqConst2REqConstDummyM[i].end(), rMobileLocEqConst2REqConstM[i].begin());
        }
    }

    // cout << "Matrix MobileLocEqConst2REqConstM_:" << "\n";
    // for (size_t i = 0; i < numMobileTx_; ++i) {
    //     for (size_t n = 0; n < numMobileTxLocs_; ++n) {
    //         cout << setw(8) << rMobileLocEqConst2REqConstM_[i][n] << " "; 
    //     }
    //     cout << "\n";
    // }

    // fstream file_output_1("./ut_results/rMobileLocEqConst2REqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<rMobileLocEqConst2REqConstM.size(); i++){
    //     for (int j=0; j<rMobileLocEqConst2REqConstM[0].size(); j++){
    //         file_output_1.write((char*)(&rMobileLocEqConst2REqConstM[i][j]), sizeof(rMobileLocEqConst2REqConstM[i][j]));}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}


}


inline void InequalityConstraintMOBARP::set_cover_ineq_const_update(vector<vector<impalib_type>> & rZEqConst2SetCoverIneqConstM, vector<vector<impalib_type>> &rFixedXEqConst2SetCoverConstM, vector<vector<vector<vector<int>>>> & rTempConxMobTxRx,
                                                   vector<vector<int>> & rTempReshapedConnectivityFixedTx, vector<vector<int>> & rTempReshapedConxMobTxRx, vector<vector<vector<impalib_type>>> &rSetCoverIneqConst2FixedXEqConstDummyM,
                                                   vector<vector<vector<vector<impalib_type>>>> &rSetCoverIneqConst2ZEqConstDummyM){

        vector<vector<impalib_type>> temp_set_cover_ineq_const_to_fixed_x_eq_const_m(numFixedTx_*numBands_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));

        vector<vector<impalib_type>> reshaped_fixed_x_eq_const_to_set_cover_const_m(numFixedTx_*numBands_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));

        vector<impalib_type> flat_data_1;
        for (const auto& row : rFixedXEqConst2SetCoverConstM) {
            for (auto &val : row) {
                flat_data_1.push_back(val);
            }
        }

        int index_1 = 0;
        for (int i = 0; i < numFixedTx_*numBands_; ++i) {
            for (int j = 0; j < numTimeSteps_*numRxLocs_; ++j) {
                reshaped_fixed_x_eq_const_to_set_cover_const_m[i][j] = flat_data_1[index_1++];
            }
        }

        //reshaped_z_eq_const_to_set_cover_ineq_const_m = z_eq_const_to_set_cover_ineq_const_m.reshape(self.num_bands*self.num_mobile_tx,self.num_time_steps,self.num_mobile_tx_locs,self.num_rx_locs).transpose()
        vector<vector<vector<vector<impalib_type>>>> reshaped_z_eq_const_to_set_cover_ineq_const_m(numRxLocs_, vector<vector<vector<impalib_type>>>(numMobileTxLocs_, vector<vector<impalib_type>>(numTimeSteps_, vector<impalib_type>(numMobileTx_*numBands_, 0))));


        vector<impalib_type> flat_data_2;
        for (const auto& row : rZEqConst2SetCoverIneqConstM) {
            for (auto &val : row) {
                flat_data_2.push_back(val);
            }
        }

        // vector<vector<vector<vector<impalib_type>>>> temp_1(numBands_*numMobileTx_, vector<vector<vector<impalib_type>>>(numTimeSteps_, vector<vector<impalib_type>>(numMobileTxLocs_, vector<impalib_type>(numRxLocs_, 0))));
        
        int index_2 = 0;
        for (int j_i = 0; j_i < numBands_*numMobileTx_; j_i++){
            for (int k=0; k<numTimeSteps_; k++){
                for (int n=0;n<numMobileTxLocs_; n++){
                    for (int l=0; l<numRxLocs_; l++){
                        // temp_1[j_i][k][n][l] = flat_data_2[index_2++];
                        // reshaped_z_eq_const_to_set_cover_ineq_const_m[l][n][k][j_i] = temp_1[j_i][k][n][l];
                        reshaped_z_eq_const_to_set_cover_ineq_const_m[l][n][k][j_i] = flat_data_2[index_2++];
                    }
                }
            }
        }


        //TempConxMobTxRx(vector<vector<vector<vector<int>>>>(NUM_RX_LOCS, vector<vector<vector<int>>>(NUM_MOBILE_TX_LOCS, vector<vector<int>>(NUM_TIME_STEPS, vector<int>(NUM_BANDS*NUM_MOBILE_TX, 0))))),
        vector<impalib_type> min_mobile_msgs(rTempConxMobTxRx[0][0].size()*rTempConxMobTxRx.size(), numeric_limits<impalib_type>::infinity());

        int flat_index_1 = 0;
        //min_mobile_msgs = np.min(reshaped_z_eq_const_to_set_cover_ineq_const_m, axis=(3,1), where = temp_conx_mob_tx_rx==1, initial=np.inf).T.flatten()
        //flipping first for loops accounts for the transpose
        for (size_t k = 0; k < numTimeSteps_; k++) {
            for (size_t l = 0; l < numRxLocs_; ++l) {
                
                impalib_type min_val = std::numeric_limits<impalib_type>::infinity();
                for (size_t n = 0; n < numMobileTxLocs_; n++) {
                    for (size_t j_i = 0; j_i < numBands_*numMobileTx_; j_i++) {
                        if (rTempConxMobTxRx[l][n][k][j_i] == 1) {
                            min_val = min(min_val, reshaped_z_eq_const_to_set_cover_ineq_const_m[l][n][k][j_i]);
                        }

                    }
                }
                min_mobile_msgs[flat_index_1] = min(min_mobile_msgs[flat_index_1], min_val);
                flat_index_1++;
            }
        }

        // // auto start_1 = chrono::high_resolution_clock::now();
        // //TempReshapedConnectivityFixedTx(vector<vector<int>>(NUM_FIXED_TX*NUM_BANDS, vector<int>(NUM_TIME_STEPS*NUM_RX_LOCS, 0)))
        // for (int index = 0; index<numFixedTx_*numBands_; index++){
        //     vector<bool> mask(numFixedTx_*numBands_, true);
        //     mask[index] = false;
        //     vector<impalib_type> min_remaining_fixed_msgs(rTempReshapedConnectivityFixedTx[0].size(), numeric_limits<impalib_type>::infinity());

        //     for (int i=0; i<numFixedTx_*numBands_; i++){
        //         if (!mask[i]) continue;
        //         for (int j=0; j<rTempReshapedConnectivityFixedTx[i].size(); j++){
        //             if (rTempReshapedConnectivityFixedTx[i][j] == 1){
        //                 min_remaining_fixed_msgs[j] = min(min_remaining_fixed_msgs[j], reshaped_fixed_x_eq_const_to_set_cover_const_m[i][j]);
        //             }
        //         }
        //     }

        //     for (int j=0; j<numTimeSteps_*numRxLocs_; j++){
        //         if (rTempReshapedConnectivityFixedTx[index][j] ==1){
        //             impalib_type min_val = 0.0;
        //             min_val = min(min_remaining_fixed_msgs[j], min_mobile_msgs[j]);
        //             temp_set_cover_ineq_const_to_fixed_x_eq_const_m[index][j] = -max(0.0, min_val);
                    
        //         }
        //     }

        // }

    // auto end_1 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_1 = end_1 - start_1;
    // cout << "Execution time unoptimized: " << elapsed_1.count() << " seconds\n";

    // auto start_2 = chrono::high_resolution_clock::now();

    vector<impalib_type> stage_forward_messages_1(numFixedTx_*numBands_ + 1, zero_value);
    vector<impalib_type> stage_backward_messages_1(numFixedTx_*numBands_  + 1, zero_value);

    // vector<vector<impalib_type>> temp_set_cover_ineq_const_to_fixed_x_eq_const_m_new(numFixedTx_*numBands_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));

    for (int k_l=0; k_l < numTimeSteps_*numRxLocs_; k_l++){

        vector<int> connections;

        for (int i_j = 0; i_j < rTempReshapedConnectivityFixedTx.size(); i_j++)
        {
            if (rTempReshapedConnectivityFixedTx[i_j][k_l]==1){
                connections.push_back(i_j);
            }
        }

        if (connections.size()==0){
            continue;
        }

        // Calculate forward messages
        stage_forward_messages_1[connections[0]] = initial_forward_message_;

        for (int stage = 1; stage < connections.size(); stage++)
        {

            stage_forward_messages_1[connections[stage]] =
                min(stage_forward_messages_1[connections[stage - 1]],
                    reshaped_fixed_x_eq_const_to_set_cover_const_m[connections[stage - 1]][k_l]);
        }

        // Calculate backward messages
        stage_backward_messages_1[connections[connections.size() - 1] + 1] = initial_backward_message_;
        
        for (size_t stage = connections.size() - 1; stage >= 1; stage--)
        {
            stage_backward_messages_1[connections[stage - 1] + 1] =
                min(stage_backward_messages_1[connections[stage] + 1],
                    reshaped_fixed_x_eq_const_to_set_cover_const_m[connections[stage]][k_l]);
        }

        for (int conx_index = 0; conx_index < connections.size(); conx_index++)
        {
            impalib_type minimumValue = min(stage_forward_messages_1[connections[conx_index]],
                                            stage_backward_messages_1[connections[conx_index] + 1]);
            minimumValue = min(minimumValue, min_mobile_msgs[k_l]);                       
            temp_set_cover_ineq_const_to_fixed_x_eq_const_m[connections[conx_index]][k_l] = -max(minimumValue, 0.0);
        }
    }
   
    // auto end_2 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_2 = end_2 - start_2;
    // cout << "Execution time optimized: " << elapsed_2.count() << " seconds\n";

    // for (int i = 0; i < numFixedTx_*numBands_; ++i) {
    //     for (int n = 0; n < numTimeSteps_*numRxLocs_; ++n) {
    //         if (std::abs(temp_set_cover_ineq_const_to_fixed_x_eq_const_m_new[i][n] - temp_set_cover_ineq_const_to_fixed_x_eq_const_m[i][n]) > 1e-7) {
    //             cout << "Error: Values are not equal at index [" << i << "][" << n << "]\n";
    //             cout<<"temp_set_cover_ineq_const_to_fixed_x_eq_const_m_new[i][n]: "<<temp_set_cover_ineq_const_to_fixed_x_eq_const_m_new[i][n]<<"\n";
    //             cout<<"temp_set_cover_ineq_const_to_fixed_x_eq_const_m[i][n]: "<<temp_set_cover_ineq_const_to_fixed_x_eq_const_m[i][n]<<"\n";
    //             exit(EXIT_FAILURE);
    //         }
    //     }
    // }

    // exit(0);
    //self.set_cover_ineq_const_to_fixed_x_eq_const_m_dummy = set_cover_ineq_const_to_fixed_x_eq_const_m.T.reshape(self.num_time_steps, self.num_rx_locs, self.num_fixed_tx*self.num_bands)
    vector<vector<impalib_type>> temp_2(numTimeSteps_*numRxLocs_, vector<impalib_type>(numFixedTx_*numBands_, 0));

    for (int index = 0; index<numFixedTx_*numBands_; index++){
        for (int j=0; j<numTimeSteps_*numRxLocs_; j++){
            temp_2[j][index] = temp_set_cover_ineq_const_to_fixed_x_eq_const_m[index][j];
        }
    }


    vector<impalib_type> flat_data_3;
    for (const auto& row : temp_2) {
        for (auto &val : row) {
            flat_data_3.push_back(val);
        }
    }

    int index_3 = 0;
    for (int k=0; k<numTimeSteps_; k++){
        for (int l=0; l<numRxLocs_; l++){
            for (int i_j=0; i_j<numFixedTx_*numBands_; i_j++){
                rSetCoverIneqConst2FixedXEqConstDummyM[k][l][i_j] = flat_data_3[index_3++];
            }
        }
    }

    // for (int k = 0; k < numTimeSteps_; ++k) {
    //     cout << "Time Step: " << k << "\n";
    //     for (int l = 0; l < numRxLocs_; ++l) {
    //         cout << "  Rx Location: " << l << "\n";
    //         for (int i_j = 0; i_j < numFixedTx_ * numBands_; ++i_j) {
    //             cout << fixed << setprecision(4) << setw(10) 
    //                  << rSetCoverIneqConst2FixedXEqConstDummyM_[k][l][i_j] << " ";
    //         }
    //         cout << "\n";
    //     }
    // }

    //vector<vector<impalib_type>> reshaped_fixed_x_eq_const_to_set_cover_const_m(numFixedTx_*numBands_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));
    //min_fixed_msgs = np.min(reshaped_fixed_x_eq_const_to_set_cover_const_m, axis=0, where=temp_reshaped_connectivity_fixed_tx==1, initial=np.inf)
    //TempReshapedConnectivityFixedTx(vector<vector<int>>(NUM_FIXED_TX*NUM_BANDS, vector<int>(NUM_TIME_STEPS*NUM_RX_LOCS, 0))){};

    vector<impalib_type> min_fixed_msgs(numTimeSteps_*numRxLocs_, 0);
    for (int k_l=0; k_l<numTimeSteps_*numRxLocs_; k_l++){
        impalib_type min_value = numeric_limits<impalib_type>::infinity();
        for (int i_j=0; i_j<numFixedTx_*numBands_; i_j++){
            if (rTempReshapedConnectivityFixedTx[i_j][k_l] == 1){
                min_value = min(min_value, reshaped_fixed_x_eq_const_to_set_cover_const_m[i_j][k_l]);
            }
        min_fixed_msgs[k_l] = min_value;
        }
    }

    //temp_reshaped_z_eq_const_to_set_cover_ineq_const_m = np.swapaxes(reshaped_z_eq_const_to_set_cover_ineq_const_m, 0, 3).reshape(-1,self.num_time_steps*self.num_rx_locs)
    // vector<vector<vector<vector<impalib_type>>>> swaped_reshaped_z_eq_const_to_set_cover_ineq_const_m(numBands_*numMobileTx_, vector<vector<vector<impalib_type>>>(numMobileTxLocs_, vector<vector<impalib_type>>(numTimeSteps_, vector<impalib_type>(numRxLocs_, 0))));
    vector<vector<impalib_type>> temp_reshaped_z_eq_const_to_set_cover_ineq_const_m(numBands_*numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));
    
    // for (int k=0; k< numTimeSteps_; k++){
    //     for (int l=0; l<numRxLocs_; l++){
    //         for (int n=0; n< numMobileTxLocs_; n++){
    //             for (int j_i=0; j_i <numBands_*numMobileTx_; j_i++){
    //                 swaped_reshaped_z_eq_const_to_set_cover_ineq_const_m[j_i][n][k][l] = reshaped_z_eq_const_to_set_cover_ineq_const_m[l][n][k][j_i];
    //             }
    //         }
    //     }
    // }

    for (int j_i = 0; j_i < numBands_*numMobileTx_; ++j_i) {
        for (int n = 0; n < numMobileTxLocs_; ++n) {
            int row_index = j_i * numMobileTxLocs_ + n;
            int col_index = 0;
            for (int k = 0; k < numTimeSteps_; ++k) {
                for (int l = 0; l < numRxLocs_; ++l) {
                    // temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[row_index][col_index] = swaped_reshaped_z_eq_const_to_set_cover_ineq_const_m[j_i][n][k][l];
                    temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[row_index][col_index] = reshaped_z_eq_const_to_set_cover_ineq_const_m[l][n][k][j_i];
                    ++col_index;
                }
            }
        }
    }

    vector<vector<impalib_type>> set_cover_ineq_const_to_z_eq_const_m(numBands_*numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));

    // auto start_3 = chrono::high_resolution_clock::now();

    // for (int index = 0; index<numBands_*numMobileTx_*numMobileTxLocs_; index++){
        
    //     vector<bool> mask(numBands_*numMobileTx_*numMobileTxLocs_, true);
    //     mask[index] = false;
    //     //rTempReshapedConxMobTxRx(vector<vector<int>>(NUM_MOBILE_TX_LOCS*NUM_BANDS*NUM_MOBILE_TX, vector<int>(NUM_RX_LOCS*NUM_TIME_STEPS, 0))),
    //     vector<impalib_type> min_remaining_mobile_msgs(rTempReshapedConxMobTxRx[0].size(), numeric_limits<impalib_type>::infinity());

    //     for (int i=0; i<numBands_*numMobileTx_*numMobileTxLocs_; i++){
    //         if (!mask[i]) continue;
    //         for (int j=0; j<rTempReshapedConxMobTxRx[i].size(); j++){
    //             if (rTempReshapedConxMobTxRx[i][j] == 1){
    //                 min_remaining_mobile_msgs[j] = min(min_remaining_mobile_msgs[j], temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[i][j]);
    //             }
    //         }
    //     }

    //     // set_cover_ineq_const_to_z_eq_const_m[index, :] = np.where(temp_reshaped_conx_mob_tx_rx[index, :], -np.maximum(zero_value, np.minimum(min_remaining_mobile_msgs, min_fixed_msgs)), set_cover_ineq_const_to_z_eq_const_m[index, :])

    //     for (int j=0; j<numTimeSteps_*numRxLocs_; j++){
    //         if (rTempReshapedConxMobTxRx[index][j] ==1){
    //             impalib_type min_val = 0.0;
    //             min_val = min(min_fixed_msgs[j], min_remaining_mobile_msgs[j]);
    //             set_cover_ineq_const_to_z_eq_const_m[index][j] = -max(0.0, min_val);
                
    //         }
    //     }
    // }


    // auto end_3 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_3 = end_3 - start_3;
    // cout << "Execution time unoptimized: " << elapsed_3.count() << " seconds\n";


    // auto start_4 = chrono::high_resolution_clock::now();

    vector<impalib_type> stage_forward_messages_2(numMobileTx_*numMobileTxLocs_*numBands_ + 1, zero_value);
    vector<impalib_type> stage_backward_messages_2(numMobileTx_*numMobileTxLocs_*numBands_  + 1, zero_value);

    // vector<vector<impalib_type>> set_cover_ineq_const_to_z_eq_const_m_new(numBands_*numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));

    for (int k_l=0; k_l < numTimeSteps_*numRxLocs_; k_l++){

        vector<int> connections;

        for (int i_j = 0; i_j < rTempReshapedConxMobTxRx.size(); i_j++)
        {
            if (rTempReshapedConxMobTxRx[i_j][k_l]==1){
                connections.push_back(i_j);
                // cout<<i_j<<" ";
            }
        }

        if (connections.size()==0){
            continue;
        }
        // cout<<"\n";
        // Calculate forward messages
        // cout<<"connections.size(): "<<connections.size()<<"\n";
        // cout<<"connections[0]: "<<connections[0]<<"\n";
        // cout<<"stage_forward_messages_2.size(): "<<stage_forward_messages_2.size()<<"\n";

        stage_forward_messages_2[connections[0]] = initial_forward_message_;

        // cout<<"DONE"<<"\n";
        // exit(0);

        // cout<<"temp_reshaped_z_eq_const_to_set_cover_ineq_const_m.size(): "<<temp_reshaped_z_eq_const_to_set_cover_ineq_const_m.size()<<"\n";
        // cout<<"temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[0].size(): "<<temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[0].size()<<"\n";
        // cout<<"k_l: "<<k_l<<"\n";

        for (int stage = 1; stage < connections.size(); stage++)
        {
            // cout<<"connections[stage]: "<<connections[stage]<<"\n";
            stage_forward_messages_2[connections[stage]] =
                min(stage_forward_messages_2[connections[stage - 1]],
                    temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[connections[stage - 1]][k_l]);
        }

        // Calculate backward messages
        stage_backward_messages_2[connections[connections.size() - 1] + 1] = initial_backward_message_;

        for (size_t stage = connections.size() - 1; stage >= 1; stage--)
        {
            stage_backward_messages_2[connections[stage - 1] + 1] =
                min(stage_backward_messages_2[connections[stage] + 1],
                    temp_reshaped_z_eq_const_to_set_cover_ineq_const_m[connections[stage]][k_l]);
        }

        for (int conx_index = 0; conx_index < connections.size(); conx_index++)
        {
            impalib_type minimumValue = min(stage_forward_messages_2[connections[conx_index]],
                                            stage_backward_messages_2[connections[conx_index] + 1]);
            minimumValue = min(minimumValue, min_fixed_msgs[k_l]);                       
            set_cover_ineq_const_to_z_eq_const_m[connections[conx_index]][k_l] = -max(minimumValue, 0.0);
        }
    }

    // auto end_4 = chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed_4 = end_4 - start_4;
    // cout << "Execution time optimized: " << elapsed_4.count() << " seconds\n";

    // for (int i = 0; i < numBands_*numMobileTx_*numMobileTxLocs_; ++i) {
    //     for (int n = 0; n < numTimeSteps_*numRxLocs_; ++n) {
    //         if (std::abs(set_cover_ineq_const_to_z_eq_const_m_new[i][n] - set_cover_ineq_const_to_z_eq_const_m[i][n]) > 1e-7) {
    //             cout << "Error: Values are not equal at index [" << i << "][" << n << "]\n";
    //             cout<<"set_cover_ineq_const_to_z_eq_const_m_new[i][n]: "<<set_cover_ineq_const_to_z_eq_const_m_new[i][n]<<"\n";
    //             cout<<"set_cover_ineq_const_to_z_eq_const_m[i][n]: "<<set_cover_ineq_const_to_z_eq_const_m[i][n]<<"\n";
    //             exit(EXIT_FAILURE);
    //         }
    //     }
    // }

    // cout<<"DONE"<<"\n";
    // exit(0);

    //vector<vector<impalib_type>> set_cover_ineq_const_to_z_eq_const_m(numBands_*numMobileTx_*numMobileTxLocs_, vector<impalib_type>(numTimeSteps_*numRxLocs_, 0));
    // vector<vector<impalib_type>> temp_transpose(numTimeSteps_*numRxLocs_, vector<impalib_type>(numBands_*numMobileTx_*numMobileTxLocs_, 0));

    // for (size_t i = 0; i < set_cover_ineq_const_to_z_eq_const_m.size(); ++i) { //numBands_*numMobileTx_*numMobileTxLocs_
    //     for (size_t j = 0; j < set_cover_ineq_const_to_z_eq_const_m[0].size(); ++j) {//numTimeSteps_*numRxLocs_
    //         temp_transpose[j][i] = set_cover_ineq_const_to_z_eq_const_m[i][j];
    //     }
    // }

    // vector<vector<vector<impalib_type>>> reshaped_set_cover_ineq_const_to_z_eq_const_m(numTimeSteps_ * numRxLocs_, vector<vector<impalib_type>>(numMobileTxLocs_, vector<impalib_type>(numBands_ * numMobileTx_, 0)));

    vector<impalib_type> flat_data_4;
    
    for (int k_l = 0; k_l < numTimeSteps_ * numRxLocs_; ++k_l) {
        for (int n = 0; n < numMobileTxLocs_; ++n) {
            for (int j_i = 0; j_i < numBands_ * numMobileTx_; ++j_i) {
                // flat_data_4.push_back(temp_transpose[k_l][n + j_i * numMobileTxLocs_]);
                flat_data_4.push_back(set_cover_ineq_const_to_z_eq_const_m[n + j_i * numMobileTxLocs_][k_l]);
                // reshaped_set_cover_ineq_const_to_z_eq_const_m[k_l][n][j_i] = 
                //     temp_transpose[k_l][n + j_i * numMobileTxLocs_];
            }
        }
    }

    // vector<impalib_type> flat_data_5;
    
    // for (int k_l = 0; k_l < numTimeSteps_ * numRxLocs_; ++k_l) {
    //     for (int n = 0; n < numMobileTxLocs_; ++n) {
    //         for (int j_i = 0; j_i < numBands_ * numMobileTx_; ++j_i) {
    //             flat_data_5.push_back(temp_transpose[k_l][n + j_i * numMobileTxLocs_]);
    //         }
    //     }
    // }

    // for (int i = 0; i < flat_data_5.size(); ++i) {
    //         // cout<<flat_data_5[i]<< ", "<<flat_data_4[i]<<"\n";
    //         if (std::abs(flat_data_5[i] - flat_data_4[i]) > 1e-7) {
    //             cout << "Error: Values are not equal at index [" << i <<"\n";
    //             cout<<"flat_data_4[i]: "<<flat_data_4[i]<<"\n";
    //             cout<<"flat_data_5[i]: "<<flat_data_5[i]<<"\n";
    //             exit(EXIT_FAILURE);
    //         }
    // }

    //temp_set_cover_ineq_const_to_z_eq_const_m = reshaped_set_cover_ineq_const_to_z_eq_const_m.reshape(self.num_time_steps, self.num_rx_locs, self.num_mobile_tx_locs, self.num_bands*self.num_mobile_tx)

    int flat_index_2 = 0;
    for (int k=0; k<numTimeSteps_; k++){
        for (int l=0; l<numRxLocs_; l++){
            for (int n=0; n<numMobileTxLocs_; n++){
                for (int j_i=0; j_i< numBands_*numMobileTx_; j_i++){
                    rSetCoverIneqConst2ZEqConstDummyM[k][l][n][j_i] = flat_data_4[flat_index_2++];
                }
            }
        }
    }

    // cout<<"DONE \n";
    // for (int k = 0; k < numTimeSteps_; ++k) {
    //     cout << "Time Step: " << k << "\n";
    //     for (int l = 0; l < numRxLocs_; ++l) {
    //         cout << "  Rx Location: " << l << "\n";
    //         for (int n = 0; n < numMobileTxLocs_; ++n) {
    //             cout << "    Mobile Tx Location: " << n << " -> ";
    //             for (int j_i = 0; j_i < numBands_ * numMobileTx_; ++j_i) {
    //                 cout << setw(12)<<setprecision(5)<<rSetCoverIneqConst2ZEqConstDummyM_[k][l][n][j_i] << " ";
    //             }
    //             cout << "\n";
    //         }
    //     }
    // }

    // exit(0);


    // fstream file_output_1("./ut_results/rSetCoverIneqConst2ZEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<rSetCoverIneqConst2ZEqConstDummyM.size(); i++){
    //     for (int j=0; j<rSetCoverIneqConst2ZEqConstDummyM[0].size(); j++){
    //         for (int k=0; k<rSetCoverIneqConst2ZEqConstDummyM[0][0].size(); k++){
    //             for (int l=0; l<rSetCoverIneqConst2ZEqConstDummyM[0][0][0].size(); l++){
    //                 file_output_1.write((char*)(&rSetCoverIneqConst2ZEqConstDummyM[i][j][k][l]), sizeof(rSetCoverIneqConst2ZEqConstDummyM[i][j][k][l]));}}}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}


    // fstream file_output_2("./ut_results/rSetCoverIneqConst2FixedXEqConstDummyM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_2.is_open()) {
    //     for (int i=0; i<rSetCoverIneqConst2FixedXEqConstDummyM.size(); i++){
    //     for (int j=0; j<rSetCoverIneqConst2FixedXEqConstDummyM[0].size(); j++){
    //         for (int k=0; k<rSetCoverIneqConst2FixedXEqConstDummyM[0][0].size(); k++){
    //                 file_output_2.write((char*)(&rSetCoverIneqConst2FixedXEqConstDummyM[i][j][k]), sizeof(rSetCoverIneqConst2FixedXEqConstDummyM[i][j][k]));}}}
    //         file_output_2.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}
    
    // cout<<"DONE"<<"\n";
}


inline void InequalityConstraintMOBARP::process_filtering_set_cover_const(int iter, vector<vector<vector<impalib_type>>> &rSetCoverIneqConst2FixedXEqConstDummyM, 
                                vector<vector<vector<vector<impalib_type>>>> &rSetCoverIneqConst2ZEqConstDummyM, vector<vector<vector<impalib_type>>> &rSetCoverIneqConst2FixedXEqConstM, 
                                vector<vector<vector<vector<impalib_type>>>> &rSetCoverIneqConst2ZEqConstM){

    //SetCoverIneqConst2FixedXEqConstOld_(NUM_TIME_STEPS, vector<vector<impalib_type>>(NUM_RX_LOCS, vector<impalib_type>(NUM_FIXED_TX*NUM_BANDS, 0))),
    //SetCoverIneqConst2ZEqConstOld_(vector<vector<vector<vector<impalib_type>>>>(NUM_TIME_STEPS, vector<vector<vector<impalib_type>>>(NUM_RX_LOCS, vector<vector<impalib_type>>(NUM_MOBILE_TX_LOCS, vector<impalib_type>(NUM_BANDS*NUM_MOBILE_TX, 0))))){};

 for (int k = 0; k < numTimeSteps_; k++){

    for (int l=0; l<numRxLocs_; l++){
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            // Calculate weighted values for current and old messages
            vector<impalib_type> intermediate_dummy(rSetCoverIneqConst2FixedXEqConstDummyM[k][l]),intermediate_old(SetCoverIneqConst2FixedXEqConstOld_[k][l]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(), rSetCoverIneqConst2FixedXEqConstM[k][l].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rSetCoverIneqConst2FixedXEqConstM[k][l].begin());
            }
            copy(rSetCoverIneqConst2FixedXEqConstM[k][l].begin(), rSetCoverIneqConst2FixedXEqConstM[k][l].end(), SetCoverIneqConst2FixedXEqConstOld_[k][l].begin());
        }

        else
        {
            copy(rSetCoverIneqConst2FixedXEqConstDummyM[k][l].begin(), rSetCoverIneqConst2FixedXEqConstDummyM[k][l].end(), rSetCoverIneqConst2FixedXEqConstM[k][l].begin());
        }

        for (int n=0; n<numMobileTxLocs_; n++){

            if ((filteringFlag_) and (alpha_ != zero_value))
            {
                // Calculate weighted values for current and old messages
                vector<impalib_type> intermediate_dummy(rSetCoverIneqConst2ZEqConstDummyM[k][l][n]),intermediate_old(SetCoverIneqConst2ZEqConstOld_[k][l][n]), intermediate_extrinsic;

                impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
                transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

                if (iter == 0)
                {
                    copy(intermediate_dummy.begin(), intermediate_dummy.end(), rSetCoverIneqConst2ZEqConstM[k][l][n].begin());
                }
                else
                {
                    transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                    copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rSetCoverIneqConst2ZEqConstM[k][l][n].begin());
                }
                copy(rSetCoverIneqConst2ZEqConstM[k][l][n].begin(), rSetCoverIneqConst2ZEqConstM[k][l][n].end(), SetCoverIneqConst2ZEqConstOld_[k][l][n].begin());
            }

            else
            {
                copy(rSetCoverIneqConst2ZEqConstDummyM[k][l][n].begin(), rSetCoverIneqConst2ZEqConstDummyM[k][l][n].end(), rSetCoverIneqConst2ZEqConstM[k][l][n].begin());
            }
            
        }
    }


    }

    // fstream file_output_1("./ut_results/rSetCoverIneqConst2ZEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_1.is_open()) {
    //     for (int i=0; i<rSetCoverIneqConst2ZEqConstM.size(); i++){
    //     for (int j=0; j<rSetCoverIneqConst2ZEqConstM[0].size(); j++){
    //         for (int k=0; k<rSetCoverIneqConst2ZEqConstM[0][0].size(); k++){
    //             for (int l=0; l<rSetCoverIneqConst2ZEqConstM[0][0][0].size(); l++){
    //                 file_output_1.write((char*)(&rSetCoverIneqConst2ZEqConstM[i][j][k][l]), sizeof(rSetCoverIneqConst2ZEqConstM[i][j][k][l]));}}}}
    //         file_output_1.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}


    // fstream file_output_2("./ut_results/rSetCoverIneqConst2FixedXEqConstM_wrapper", ios::out | ios::binary | ios:: trunc);
    // if (file_output_2.is_open()) {
    //     for (int i=0; i<rSetCoverIneqConst2FixedXEqConstM.size(); i++){
    //     for (int j=0; j<rSetCoverIneqConst2FixedXEqConstM[0].size(); j++){
    //         for (int k=0; k<rSetCoverIneqConst2FixedXEqConstM[0][0].size(); k++){
    //                 file_output_2.write((char*)(&rSetCoverIneqConst2FixedXEqConstM[i][j][k]), sizeof(rSetCoverIneqConst2FixedXEqConstM[i][j][k]));}}}
    //         file_output_2.close();
    //         }
    // else {cout << "Error! File cannot be opened!" << "\n";}


    // for (int k = 0; k < numTimeSteps_; ++k) {
    //     cout << "Time Step: " << k << "\n";
    //     for (int l = 0; l < numRxLocs_; ++l) {
    //         cout << "  Rx Location: " << l << "\n";
    //         for (int i_j = 0; i_j < numFixedTx_ * numBands_; ++i_j) {
    //             cout << fixed << setprecision(4) << setw(10) 
    //                  << rSetCoverIneqConst2FixedXEqConstM_[k][l][i_j] << " ";
    //         }
    //         cout << "\n";
    //     }
    // }

    // for (int k = 0; k < numTimeSteps_; ++k) {
    //     cout << "Time Step: " << k << "\n";
    //     for (int l = 0; l < numRxLocs_; ++l) {
    //         cout << "  Rx Location: " << l << "\n";
    //         for (int n = 0; n < numMobileTxLocs_; ++n) {
    //             cout << "    Mobile Tx Location: " << n << " -> ";
    //             for (int j_i = 0; j_i < numBands_ * numMobileTx_; ++j_i) {
    //                 cout << setw(12)<<setprecision(5)<<rSetCoverIneqConst2ZEqConstM_[k][l][n][j_i] << " ";
    //             }
    //             cout << "\n";
    //         }
    //     }
    // }  

}