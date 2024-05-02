// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "impalib/impalib.hpp"

class KsatConstraint
{
private:
    int                          numVariables_;
    int                          numConstraints_;
    int                          kVariable_;
    bool                         filteringFlag_;
    impalib_type                 alpha_;
    vector<vector<impalib_type>> ksatConstraint2EqConstraintOld;

public:
    void ksat_constraint_to_variable_ec_update(vector<vector<impalib_type>> &, vector<vector<impalib_type>> &, vector<vector<int>> &, vector<vector<int>> &);
    void process_filtering(int, vector<vector<impalib_type>> &, vector<vector<impalib_type>> &);

    KsatConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG,
                     const impalib_type ALPHA);
};


KsatConstraint::KsatConstraint(const int NUM_VARIABLES, const int NUM_CONSTRAINTS, const int K_VARIABLE, const bool FILTERING_FLAG,
                                   const impalib_type ALPHA)
    : filteringFlag_(FILTERING_FLAG), alpha_(ALPHA), numVariables_(NUM_VARIABLES), numConstraints_(NUM_CONSTRAINTS), kVariable_(K_VARIABLE)
{
    ksatConstraint2EqConstraintOld.reserve(numConstraints_);

    for (int constraint_index = 0; constraint_index < numConstraints_; constraint_index++)
    {
        ksatConstraint2EqConstraintOld.push_back(vector<impalib_type>(numVariables_, zero_value));
    }
};

void KsatConstraint::ksat_constraint_to_variable_ec_update(vector<vector<impalib_type>> &rVariableEc2KsatConstraintM, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_, vector<vector<int>> &rConstraintsConnections, vector<vector<int>> &rConstraintsConnectionsType)
{   
    for (int constraint_index = 0; constraint_index < numConstraints_; ++constraint_index) {
            auto& connections = rConstraintsConnections[constraint_index];
            auto& connections_type = rConstraintsConnectionsType[constraint_index];
            
            for (size_t variable_index = 0; variable_index < connections.size(); ++variable_index) {
                
                auto variable = connections[variable_index];
                
                vector<int> connections_solid;
                vector<int> connections_dashed;

                for (size_t i = 0; i < connections.size(); ++i) {
                    if (connections_type[i] == 1 && connections[i] != variable) {
                        connections_solid.push_back(connections[i]);
                    } else if (connections_type[i] == -1 && connections[i] != variable) {
                        connections_dashed.push_back(connections[i]);
                    }
                }

                impalib_type optimum_solid_messages = zero_value, optimum_dashed_messages = zero_value;

                if (connections_solid.empty()) {
                    optimum_solid_messages = (connections_type[variable_index] == 1) ? -value_inf : value_inf;
                } else {

                    impalib_type minimumValue = value_inf;

                    for (int j=0; j<connections_solid.size(); j++){
                        if (rVariableEc2KsatConstraintM[constraint_index][connections_solid[j]]<minimumValue){
                            minimumValue = rVariableEc2KsatConstraintM[constraint_index][connections_solid[j]];
                        }
                    }

                    if (connections_type[variable_index] == 1) {
                            optimum_solid_messages = -minimumValue;
                    }  
                        else {
                            optimum_solid_messages = minimumValue;
                        }
                    }
                
                if (connections_dashed.empty()) {
                    optimum_dashed_messages = (connections_type[variable_index] == 1) ? -value_inf : value_inf;
                } else {

                    impalib_type maximumValue = -value_inf;

                    for (int j=0; j<connections_dashed.size(); j++){
                        if (rVariableEc2KsatConstraintM[constraint_index][connections_dashed[j]]>maximumValue){
                            maximumValue = rVariableEc2KsatConstraintM[constraint_index][connections_dashed[j]];
                        }
                    }
                    
                    if (connections_type[variable_index] == 1) {
                            optimum_dashed_messages = maximumValue;
                    }  
                        else {
                            optimum_dashed_messages = -maximumValue;
                        }
                    }

                if (connections_type[variable_index] == 1) {
                    rKsatConstraint2EqConstraintDummyM_[constraint_index][variable] = min(zero_value, max(optimum_solid_messages, optimum_dashed_messages));
                } else {
                    rKsatConstraint2EqConstraintDummyM_[constraint_index][variable] = max(zero_value, min(optimum_solid_messages, optimum_dashed_messages));
                }
            }
    }
}

void KsatConstraint::process_filtering(int iter, vector<vector<impalib_type>> &rKsatConstraint2EqConstraintDummyM_,
                                         vector<vector<impalib_type>> &rKsatConstraint2EqConstraintM_)
{   

    for (int constraint_index = 0; constraint_index < numConstraints_; constraint_index++)
    {
        if ((filteringFlag_) and (alpha_ != zero_value))
        {
            vector<impalib_type> intermediate_dummy(rKsatConstraint2EqConstraintDummyM_[constraint_index]),
                intermediate_old(ksatConstraint2EqConstraintOld[constraint_index]), intermediate_extrinsic;

            impalib_type w_1 = alpha_, w_2 = 1 - alpha_;
            transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_dummy.begin(),
                      [w_2](impalib_type &c) { return c * w_2; });
            transform(intermediate_old.begin(), intermediate_old.end(), intermediate_old.begin(),
                      [w_1](impalib_type &c) { return c * w_1; });

            if (iter == 0)
            {
                copy(intermediate_dummy.begin(), intermediate_dummy.end(),
                     rKsatConstraint2EqConstraintM_[constraint_index].begin());
            }
            else
            {
                transform(intermediate_dummy.begin(), intermediate_dummy.end(), intermediate_old.begin(),
                          back_inserter(intermediate_extrinsic), plus<impalib_type>());
                copy(intermediate_extrinsic.begin(), intermediate_extrinsic.end(),
                     rKsatConstraint2EqConstraintM_[constraint_index].begin());
            }

            copy(rKsatConstraint2EqConstraintM_[constraint_index].begin(),
                 rKsatConstraint2EqConstraintM_[constraint_index].end(),
                 ksatConstraint2EqConstraintOld[constraint_index].begin());
        }

        else
        {
            copy(rKsatConstraint2EqConstraintDummyM_[constraint_index].begin(),
                 rKsatConstraint2EqConstraintDummyM_[constraint_index].end(),
                 rKsatConstraint2EqConstraintM_[constraint_index].begin());
        }
    }
}
