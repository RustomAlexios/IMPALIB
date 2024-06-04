// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

/**
 * Represents a class for K-SAT constraint
 */
class KsatConstraint {
   private:
    int numVariables_;                                            ///< total number of variables
    int numConstraints_;                                          ///< number of constraints
    int kVariable_;                                               ///< number of variables per constraint
    bool filteringFlag_;                                          ///< filtering flag
    impalib_type alpha_;                                          ///< filtering parameter
    vector<vector<impalib_type>> ksatConstraint2EqConstraintOld;  ///< messages from k-sat constraints to equality constraints before filtering
    //impalib_type                 initial_forward_message_ = value_inf; ///< initial forward message of forward-backward algorithm
    //impalib_type                 initial_backward_message_ = zero_value; ///< initial backward message of forward-backward algorithm
    int maxState_ = 1; 

   public:
    void ksat_constraint_to_variable_ec_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &,
                                               vector<vector<int>> &);                            ///< update messages from k-sat constraints to equality constraints
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &);  ///< perform filtering

    KsatConstraint(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE, bool FILTERING_FLAG,
                   impalib_type ALPHA);  ///< constructor
};

/**
 * Construct KsatConstraint object for the K-SAT problem
 *
 * @param[in] NUM_VARIABLES: total number of variables
 * @param[in] NUM_CONSTRAINTS: number of constraints
 * @param[in] K_VARIABLE: number of variables per constraint
 * @param[in] FILTERING_FLAG: filtering on or off
 * @param[in] ALPHA: filtering parameter value (between 0 and 1)
 * @param[out] filteringFlag_: FILTERING_FLAG
 * @param[out] alpha_: ALPHA
 * @param[out] numVariables_: NUM_VARIABLES
 * @param[out] numConstraints_: NUM_CONSTRAINTS
 * @param[out] kVariable_: K_VARIABLE
 */

KsatConstraint::KsatConstraint(int NUM_VARIABLES, int NUM_CONSTRAINTS, int K_VARIABLE, bool FILTERING_FLAG, impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG),
      alpha_(ALPHA),
      numVariables_(NUM_VARIABLES),
      numConstraints_(NUM_CONSTRAINTS),
      kVariable_(K_VARIABLE),
      ksatConstraint2EqConstraintOld(numConstraints_, vector<impalib_type>(numVariables_, zero_value)){};

/**
 * Calculate messages from k-sat constraints to variable equality constraints for the K-SAT problem
 *
 * @param[in] rVariableEc2KsatConstraintM: messages from variable equality constraints to k-sat constraints
 * @param[out] rKsatConstraint2EqConstraintDummyM_: messages from k-sat constraints to variable equality constraints before filtering
 * @param[in] rConstraintsConnections: connections to variables for each constraint
 * @param[in] rConstraintsConnectionsType: types of connections to variables for each constraint
 *
 */

void KsatConstraint::ksat_constraint_to_variable_ec_update(vector<vector<impalib_type>> &rVariableEc2KsatConstraintM, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_,
                                                           vector<vector<int>> &rConstraintsConnections, vector<vector<int>> &rConstraintsConnectionsType) {
    
    bool optimized_flag = true;

    if (optimized_flag){
    for(auto& row : rKsatConstraint2EqConstraintDummyM_) {
        row.assign(row.size(), zero_value);
    }

    vector<vector<impalib_type>> stage_forward_messages(kVariable_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages(kVariable_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

    for (int c = 0; c < numConstraints_; ++c) {

        auto &conx = rConstraintsConnections[c];
        auto &conx_type = rConstraintsConnectionsType[c];

        vector<impalib_type> initial_forward_messages(maxState_ + 1, zero_value),
            initial_backward_messages(maxState_ + 1, zero_value);

        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages[0] = initial_forward_messages;

        for (int stage = 0; stage < kVariable_; stage++)
        {
            stage_forward_messages[stage + 1][0] = stage_forward_messages[stage][0];
            stage_forward_messages[stage + 1][1] =
                min(stage_forward_messages[stage][1] + min(zero_value,rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage]), stage_forward_messages[stage][0]+rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage]);
        }

        stage_backward_messages[kVariable_] = initial_backward_messages;

        for (int stage = kVariable_ - 1; stage >= 0; stage--)
        {
            if (stage == kVariable_ - 1){
                stage_backward_messages[stage][0] = rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage];
            }
            else{
                stage_backward_messages[stage][0] = min(stage_backward_messages[stage + 1][0], stage_backward_messages[stage + 1][1]+rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage]);
            }
            stage_backward_messages[stage][1] = min(stage_backward_messages[stage+1][1], stage_backward_messages[stage+1][1]+ rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage]);
        }

        impalib_type min_dashed_edges = zero_value;
        impalib_type min_solid_edges = zero_value;

        for (int stage = 0; stage < kVariable_; stage++)
        {
            min_solid_edges = min(stage_forward_messages[stage][0] + stage_backward_messages[stage+1][1] + rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage], stage_forward_messages[stage][1] + rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage] + stage_backward_messages[stage+1][1]);
            
            if (stage == kVariable_-1){
                min_dashed_edges = stage_forward_messages[stage][1] + stage_backward_messages[stage+1][1];
            }
            else {
                min_dashed_edges = min(stage_forward_messages[stage][0] + stage_backward_messages[stage+1][0], stage_forward_messages[stage][1] + stage_backward_messages[stage+1][1]);
            }
            
            // The if statements checks (abs(rKsatConstraint2EqConstraintDummyM_[c][conx[stage]])< abs(((min_solid_edges  - min_dashed_edges - rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage])*conx_type[stage])))
            // were added to account for the fact that a variable could occur multiple times in a constraint,
            // like the case of benchmarks, but in practice, a constraint should not have the same variable appearing more than once
            // and thus these if statements checks can be removed
            if (abs(rKsatConstraint2EqConstraintDummyM_[c][conx[stage]])< abs(((min_solid_edges  - min_dashed_edges - rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage])*conx_type[stage]))){
            rKsatConstraint2EqConstraintDummyM_[c][conx[stage]] = (min_solid_edges  - min_dashed_edges - rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage])*conx_type[stage];
            }
        }
    }

    }

    else{

    // if you want to use unoptimized version of the code, uncomment this part and comment rKsatConstraint2EqConstraintDummyMImproved initialization
    //for(auto& row : rKsatConstraint2EqConstraintDummyM_) {
    //    row.assign(row.size(), zero_value);
    //}
    vector<vector<impalib_type>> rKsatConstraint2EqConstraintDummyMImproved(numConstraints_, vector<impalib_type>(numVariables_, zero_value));

    for (int c = 0; c < numConstraints_; ++c) {
        // Get connections and connection types for the current constraint
        auto &conx = rConstraintsConnections[c];
        auto &conx_type = rConstraintsConnectionsType[c];

        for (size_t v = 0; v < conx.size(); ++v) {
            // Get the current connected variable
            auto variable = conx[v];

            vector<int> conx_solid;
            vector<int> conx_dashed;

            // Separate solid and dashed connections
            for (size_t i = 0; i < conx.size(); ++i) {
                if (conx_type[i] == 1){// && conx[i] != variable) {
                    // if (i!=v) was added && conx[i] != variable) was commented
                    // to take into consideration that a variable can appear more than once
                    // in a constraint (like bencmarks), but this is not common in practice
                    if (i!=v) conx_solid.push_back(conx[i]);
                } else if (conx_type[i] == -1){// && conx[i] != variable) {
                    if (i!=v) conx_dashed.push_back(conx[i]);
                }
            }

            impalib_type opt_solid_msgs = zero_value, opt_dashed_msgs = zero_value;

            // Calculate optimum messages for solid connections
            if (conx_solid.empty()) {
                opt_solid_msgs = (conx_type[v] == 1) ? -value_inf : value_inf;
            } else {
                impalib_type minValue = value_inf;

                for (int j = 0; j < conx_solid.size(); j++) {
                    if (rVariableEc2KsatConstraintM[c][conx_solid[j]] < minValue) {
                        minValue = rVariableEc2KsatConstraintM[c][conx_solid[j]];
                    }
                }

                if (conx_type[v] == 1) {
                    opt_solid_msgs = -minValue;
                } else {
                    opt_solid_msgs = minValue;
                }
            }

            // Calculate optimum messages for dashed connections
            if (conx_dashed.empty()) {
                opt_dashed_msgs = (conx_type[v] == 1) ? -value_inf : value_inf;
            } else {
                impalib_type maxValue = -value_inf;

                for (int j = 0; j < conx_dashed.size(); j++) {
                    if (rVariableEc2KsatConstraintM[c][conx_dashed[j]] > maxValue) {
                        maxValue = rVariableEc2KsatConstraintM[c][conx_dashed[j]];
                    }
                }

                if (conx_type[v] == 1) {
                    opt_dashed_msgs = maxValue;
                } else {
                    opt_dashed_msgs = -maxValue;
                }
            }

            // Calculate the final message for the current variable
            // The if statements checks (abs(rKsatConstraint2EqConstraintDummyMImproved[c][variable])<abs(min(zero_value, max(opt_solid_msgs, opt_dashed_msgs))))
            // were added to account for the fact that a variable could occur multiple times in a constraint,
            // like the case of benchmarks, but in practice, a constraint should not have the same variable appearing more than once
            // and thus these if statements can be removed
            if (conx_type[v] == 1) {
                if (abs(rKsatConstraint2EqConstraintDummyMImproved[c][variable])<abs(min(zero_value, max(opt_solid_msgs, opt_dashed_msgs)))){
                    //comment rKsatConstraint2EqConstraintDummyMImproved to use unimproved version of the code, and uncomment rKsatConstraint2EqConstraintDummyM_ to use improved version of the code
                    rKsatConstraint2EqConstraintDummyMImproved[c][variable] = min(zero_value, max(opt_solid_msgs, opt_dashed_msgs));
                    //rKsatConstraint2EqConstraintDummyM_[c][variable] = min(zero_value, max(opt_solid_msgs, opt_dashed_msgs));
                }
            } else {
                if (abs(rKsatConstraint2EqConstraintDummyMImproved[c][variable])<abs(max(zero_value, min(opt_solid_msgs, opt_dashed_msgs)))){
                //comment rKsatConstraint2EqConstraintDummyMImproved to use unimproved version of the code, and uncomment rKsatConstraint2EqConstraintDummyM_ to use improved version of the code
                rKsatConstraint2EqConstraintDummyMImproved[c][variable] = max(zero_value, min(opt_solid_msgs, opt_dashed_msgs));
                //rKsatConstraint2EqConstraintDummyM_[c][variable] = max(zero_value, min(opt_solid_msgs, opt_dashed_msgs));
                }
            }
        }
    }
    }
}

/**
 * Perform filtering on messages from k-sat constraints to variable equality constraints for the K-SAT problem
 *
 * @param[in] iter: iteration index of IMPA
 * @param[in] rKsatConstraint2EqConstraintDummyM_: messages from k-sat constraints to variable equality constraints before filtering
 * @param[in] rKsatConstraint2EqConstraintM_: messages from k-sat constraints to variable equality constraints after filtering
 *
 */

void KsatConstraint::process_filtering(int iter, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintM_) {
    for (int c = 0; c < numConstraints_; c++) {
        if ((filteringFlag_) and (alpha_ != zero_value)) {
            vector<impalib_type> intermediate_dummy(rKsatConstraint2EqConstraintDummyM_[c]), intermediate_old(ksatConstraint2EqConstraintOld[c]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](impalib_type &c) { return c * w_1; });

            if (iter == 0) {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(), rKsatConstraint2EqConstraintM_[c].begin());
            } else {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(), back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(), rKsatConstraint2EqConstraintM_[c].begin());
            }
            copy(rKsatConstraint2EqConstraintM_[c].begin(), rKsatConstraint2EqConstraintM_[c].end(), ksatConstraint2EqConstraintOld[c].begin());
        }

        else {
            copy(rKsatConstraint2EqConstraintDummyM_[c].begin(), rKsatConstraint2EqConstraintDummyM_[c].end(), rKsatConstraint2EqConstraintM_[c].begin());
        }
    }
}
