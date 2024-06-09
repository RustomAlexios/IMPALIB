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
    int nVars_;                                 ///< total number of variables
    int nConstraints_;                          ///< number of constraints
    int k_;                                     ///< number of variables per constraint
    bool doFilter_;                             ///< filtering flag
    impalib_type alpha_;                        ///< filtering parameter
    vector<vector<impalib_type>> ksat2EqOldM_;  ///< messages from k-sat constraints to equality constraints before filtering
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
    : doFilter_(FILTERING_FLAG), alpha_(ALPHA), nVars_(NUM_VARIABLES), nConstraints_(NUM_CONSTRAINTS), k_(K_VARIABLE), ksat2EqOldM_(nConstraints_, vector<impalib_type>(nVars_, zero_value)){};

/**
 * Calculate messages from k-sat constraints to variable equality constraints for the K-SAT problem
 *
 * @param[in] var2KsatM: messages from variable equality constraints to k-sat constraints
 * @param[out] ksat2EqPreM: messages from k-sat constraints to variable equality constraints before filtering
 * @param[in] connections: connections to variables for each constraint
 * @param[in] types: types of connections to variables for each constraint
 *
 */

inline void KsatConstraint::ksat_constraint_to_variable_ec_update(const vector<vector<impalib_type>> &var2KsatM, vector<vector<impalib_type>> &ksat2EqPreM, const vector<vector<int>> &connections,
                                                                  const vector<vector<int>> &types) const {
    for (auto &row : ksat2EqPreM) {
        row.assign(row.size(), zero_value);
    }

    vector<vector<impalib_type>> forward(k_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));
    vector<vector<impalib_type>> backward(k_ + 1, vector<impalib_type>(maxState_ + 1, zero_value));

    for (int c = 0; c < nConstraints_; ++c) {
        auto &conx = connections[c];
        auto &conx_type = types[c];

        vector<impalib_type> forward0(maxState_ + 1, zero_value);
        vector<impalib_type> backward0(maxState_ + 1, zero_value);

        fill(forward0.begin() + 1, forward0.end(), value_inf);

        forward[0] = forward0;

        for (int stage = 0; stage < k_; stage++) {
            forward[stage + 1][0] = forward[stage][0];
            forward[stage + 1][1] = min(forward[stage][1] + min(zero_value, var2KsatM[c][conx[stage]] * conx_type[stage]), forward[stage][0] + var2KsatM[c][conx[stage]] * conx_type[stage]);
        }

        backward[k_] = backward0;

        for (int stage = k_ - 1; stage >= 0; stage--) {
            if (stage == k_ - 1) {
                backward[stage][0] = var2KsatM[c][conx[stage]] * conx_type[stage];
            } else {
                backward[stage][0] = min(backward[stage + 1][0], backward[stage + 1][1] + var2KsatM[c][conx[stage]] * conx_type[stage]);
            }
            backward[stage][1] = min(backward[stage + 1][1], backward[stage + 1][1] + var2KsatM[c][conx[stage]] * conx_type[stage]);
        }

        impalib_type min_dashed_edges = zero_value;
        impalib_type min_solid_edges = zero_value;

        for (int stage = 0; stage < k_; stage++) {
            min_solid_edges = min(forward[stage][0] + backward[stage + 1][1] + var2KsatM[c][conx[stage]] * conx_type[stage],
                                  forward[stage][1] + var2KsatM[c][conx[stage]] * conx_type[stage] + backward[stage + 1][1]);

            if (stage == k_ - 1) {
                min_dashed_edges = forward[stage][1] + backward[stage + 1][1];
            } else {
                min_dashed_edges = min(forward[stage][0] + backward[stage + 1][0], forward[stage][1] + backward[stage + 1][1]);
            }

            // The if statements checks (abs(rKsatConstraint2EqConstraintDummyM_[c][conx[stage]])< abs(((min_solid_edges  - min_dashed_edges -
            // rVariableEc2KsatConstraintM[c][conx[stage]]*conx_type[stage])*conx_type[stage]))) were added to account for the fact that a variable could occur multiple times in a constraint, like the
            // case of benchmarks, but in practice, a constraint should not have the same variable appearing more than once and thus these if statements checks can be removed
            if (abs(ksat2EqPreM[c][conx[stage]]) < abs(((min_solid_edges - min_dashed_edges - var2KsatM[c][conx[stage]] * conx_type[stage]) * conx_type[stage]))) {
                ksat2EqPreM[c][conx[stage]] = (min_solid_edges - min_dashed_edges - var2KsatM[c][conx[stage]] * conx_type[stage]) * conx_type[stage];
            }
        }
    }
}

/**
 * Perform filtering on messages from k-sat constraints to variable equality constraints for the K-SAT problem
 *
 * @param[in] iter: iteration index of IMPA
 * @param[in] ksat2EqPreM: messages from k-sat constraints to variable equality constraints before filtering
 * @param[in] ksat2EqM: messages from k-sat constraints to variable equality constraints after filtering
 *
 */

inline void KsatConstraint::process_filtering(const int iter, vector<vector<impalib_type>> &ksat2EqPreM, vector<vector<impalib_type>> &ksat2EqM) {
    if (!doFilter_) {
        ksat2EqM = ksat2EqPreM;
        return;
    }

    // Calculate weighted values for current and old messages
    for (int c = 0; c < nConstraints_; c++) {
        if (iter == 0) {
            ksat2EqM[c] = ksat2EqPreM[c];
        } else {
            vector<impalib_type> weighted(ksat2EqPreM[c].size());
            for (int i=0; i<weighted.size(); ++i) {
                weighted[i] = alpha_*ksat2EqOldM_[c][i] + (1-alpha_)*ksat2EqPreM[c][i];
            }
            ksat2EqM[c] = weighted;
        }
        ksat2EqOldM_[c] = ksat2EqM[c];
    }
}
