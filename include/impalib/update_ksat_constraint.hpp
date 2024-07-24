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
    // impalib_type                 initial_forward_message_ = value_inf; ///< initial forward message of forward-backward algorithm
    // impalib_type                 initial_backward_message_ = zero_value; ///< initial backward message of forward-backward algorithm
    int maxState_ = 1;

   public:
    void ksat_constraint_to_variable_ec_update(const vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, const vector<vector<int>> &,
                                               const vector<vector<int>> &) const;                ///< update messages from k-sat constraints to equality constraints
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
 */

inline KsatConstraint::KsatConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG, const impalib_type ALPHA)
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

inline void KsatConstraint::ksat_constraint_to_variable_ec_update(const vector<vector<impalib_type>> &rVariableEc2KsatConstraintM, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_,
                                                                  const vector<vector<int>> &rConstraintsConnections, const vector<vector<int>> &rConstraintsConnectionsType) const {
    for (auto &row : rKsatConstraint2EqConstraintDummyM_) {
        row.assign(row.size(), zero_value);
    }

    vector<vector<impalib_type>> stage_forward_messages(kVariable_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
    vector<vector<impalib_type>> stage_backward_messages(kVariable_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

    for (int c = 0; c < numConstraints_; ++c) {
        auto &conx = rConstraintsConnections[c];
        auto &conx_type = rConstraintsConnectionsType[c];

        vector<impalib_type> initial_forward_messages(maxState_ + 1, zero_value), initial_backward_messages(maxState_ + 1, zero_value);

        fill(initial_forward_messages.begin() + 1, initial_forward_messages.end(), value_inf);

        stage_forward_messages[0] = initial_forward_messages;

        for (int stage = 0; stage < kVariable_; stage++) {
            stage_forward_messages[stage + 1][0] = stage_forward_messages[stage][0];
            stage_forward_messages[stage + 1][1] = min(stage_forward_messages[stage][1] + min(zero_value, rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]),
                                                       stage_forward_messages[stage][0] + rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]);
        }

        stage_backward_messages[kVariable_] = initial_backward_messages;

        for (int stage = kVariable_ - 1; stage >= 0; stage--) {
            if (stage == kVariable_ - 1) {
                stage_backward_messages[stage][0] = rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage];
            } else {
                stage_backward_messages[stage][0] = min(stage_backward_messages[stage + 1][0], stage_backward_messages[stage + 1][1] + rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]);
            }
            stage_backward_messages[stage][1] = min(stage_backward_messages[stage + 1][1], stage_backward_messages[stage + 1][1] + rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]);
        }

        impalib_type min_dashed_edges = zero_value;
        impalib_type min_solid_edges = zero_value;

        for (int stage = 0; stage < kVariable_; stage++) {
            min_solid_edges = min(stage_forward_messages[stage][0] + stage_backward_messages[stage + 1][1] + rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage],
                                  stage_forward_messages[stage][1] + rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage] + stage_backward_messages[stage + 1][1]);

            if (stage == kVariable_ - 1) {
                min_dashed_edges = stage_forward_messages[stage][1] + stage_backward_messages[stage + 1][1];
            } else {
                min_dashed_edges = min(stage_forward_messages[stage][0] + stage_backward_messages[stage + 1][0], stage_forward_messages[stage][1] + stage_backward_messages[stage + 1][1]);
            }

            // The if statements checks (abs(rKsatConstraint2EqConstraintDummyM_[c][conx[stage]])< abs(((min_solid_edges  - min_dashed_edges -
            // rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage])*conx_type[stage]))) were added to account for the fact that a variable could occur multiple times in a constraint, like the
            // case of benchmarks, but in practice, a constraint should not have the same variable appearing more than once and thus these if statements checks can be removed
            if (abs(rKsatConstraint2EqConstraintDummyM_[c][conx[stage]]) <
                abs(((min_solid_edges - min_dashed_edges - rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]) * conx_type[stage]))) {
                rKsatConstraint2EqConstraintDummyM_[c][conx[stage]] = (min_solid_edges - min_dashed_edges - rVariableEc2KsatConstraintM[c][conx[stage]] * conx_type[stage]) * conx_type[stage];
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

inline void KsatConstraint::process_filtering(const int iter, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintM_) {
    for (int c = 0; c < numConstraints_; c++) {
        if ((filteringFlag_) and (alpha_ != zero_value)) {
            vector<impalib_type> intermediate_dummy(rKsatConstraint2EqConstraintDummyM_[c]), intermediate_old(ksatConstraint2EqConstraintOld[c]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(), [w_2](const impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(), [w_1](const impalib_type &c) { return c * w_1; });

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
