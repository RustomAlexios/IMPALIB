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

   public:
    void ksat_constraint_to_variable_ec_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &,
                                               vector<vector<int>> &);                            ///< update messages from k-sat constraints to equality constraints
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &);  ///< perform filtering

    KsatConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG,
                   const impalib_type ALPHA);  ///< constructor
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

KsatConstraint::KsatConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG, const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numVariables_(NUM_VARIABLES), numConstraints_(NUM_CONSTRAINTS), kVariable_(K_VARIABLE), ksatConstraint2EqConstraintOld(numConstraints_, vector<impalib_type>(numVariables_, zero_value)) {
    // Reserve memory
    // ksatConstraint2EqConstraintOld.reserve(numConstraints_);

    //for (int i = 0; i < numConstraints_; i++) {
    //    // Initialize a vector of size numVariables_ with zero_value and add it to ksatConstraint2EqConstraintOld vector
    //    ksatConstraint2EqConstraintOld.push_back(vector<impalib_type>(numVariables_, zero_value));
    //}
};

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
                if (conx_type[i] == 1 && conx[i] != variable) {
                    conx_solid.push_back(conx[i]);
                } else if (conx_type[i] == -1 && conx[i] != variable) {
                    conx_dashed.push_back(conx[i]);
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
            if (conx_type[v] == 1) {
                rKsatConstraint2EqConstraintDummyM_[c][variable] = min(zero_value, max(opt_solid_msgs, opt_dashed_msgs));
            } else {
                rKsatConstraint2EqConstraintDummyM_[c][variable] = max(zero_value, min(opt_solid_msgs, opt_dashed_msgs));
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
