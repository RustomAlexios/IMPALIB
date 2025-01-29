// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "impalib/impalib.hpp"
#include <iomanip>

int main()
{

    //run main_mobarp.py by setting np.random.seed(17) to reproduce the same experiment here
    //python3 main_mobarp.py --randomTestFlag=True --nFTX=2 --nMTX=2 --nBands=2 --nTimeSteps=2 --nRX=2 --nMTXLo=2 --filteringFlag=True --alpha=0.1 --percNegIM=100 --nITER=200
    
    const int  NUM_ITERATIONS           = 200;  ///< number of iterations of IMPA
    const int NUM_FIXED_TX = 2; ///< number of fixed TX
    const int NUM_MOBILE_TX = 2; ///< number of mobile TX
    const int NUM_BANDS = 2; ///< number of bands
    const int NUM_TIME_STEPS = 2; ///< number of time steps
    const int NUM_RX_LOCS = 2; ///< number of rx
    const int NUM_MOBILE_TX_LOCS = 2; ///< number of mobile tx locs
    const bool FILTERING_FLAG        = true; ///< whether filtering is activated or not
    const impalib_type ALPHA             = 0.5;     ///< filtering parameter
    const bool EXCLUDE_CAP_FLAG = false; ///< exclude capacities or not
    const impalib_type THRESHOLD         = -0.0001; ///< threshold parameter

    GraphicalModelMOBARP model_graph(NUM_ITERATIONS, NUM_FIXED_TX, NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS, NUM_MOBILE_TX_LOCS, ALPHA, FILTERING_FLAG, EXCLUDE_CAP_FLAG);
    
    array<array<array<impalib_type, NUM_TIME_STEPS>, NUM_BANDS>, NUM_FIXED_TX> fixed_x_costs = {{
        {{{{-36.51985024, -57.752808}}, {{-27.23687083, -16.11103224}}}},
        {{{{-80.8286914, -69.07001696}},{{-67.37688064, -61.80426044}}}}
        }};

    array<array<array<impalib_type, NUM_TIME_STEPS>, NUM_BANDS>, NUM_MOBILE_TX> mobile_x_costs = 
    {{
        {{
            {{-13.51566246, -42.2032244}}, {{-95.11148681, -15.40402123}}
        }},
        {{
            {{-87.76378932, -88.95614735}},{{-14.60742991, -68.71767539}}
        }}
    }};

    array<array<impalib_type, NUM_MOBILE_TX_LOCS>, NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS> z_costs = {0};

    array<array<impalib_type, NUM_MOBILE_TX_LOCS>, NUM_MOBILE_TX> r_costs = {{{6.51377099, 3.45693075},{-5.45643702, -0.41414802}}};

    array<array<array<array<int, NUM_RX_LOCS>, NUM_TIME_STEPS>, NUM_BANDS>, NUM_MOBILE_TX_LOCS> connectivity_mobile_tx = 
        {{
            {{
                {{
                    {0, 1}, {0, 0}
                }},
                {{
                    {0, 0}, {1, 1}
                }},

            }},

            {{
                {{
                    {0, 0}, {0, 1}
                }},
                {{
                    {1, 1}, {1, 1}
                }},
            }}
        }};

    array<array<array<array<int, NUM_RX_LOCS>, NUM_TIME_STEPS>, NUM_BANDS>, NUM_MOBILE_TX_LOCS> connectivity_fixed_tx = 
        {{
            {{
                {{
                    {1, 1}, {0, 1}
                }},
                {{
                    {1, 1}, {1, 0}
                }},

            }},

            {{
                {{
                    {1, 0}, {1, 1}
                }},
                {{
                    {0, 0}, {0, 1}
                }},
            }}
        }};

    array<int, NUM_FIXED_TX>             fixed_capacity_constraints    = {1, 1};
    array<int, NUM_MOBILE_TX>             mobile_capacity_constraints    = {1, 1};

    //pre-processed version of connectivity_fixed_tx (performed in Python)
    array<array<int, NUM_RX_LOCS>, NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS> conx_fixed_tx_per_num_rx_locs = {{{1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, 0}, {1, 1}, {0, 0}, {0, 1}}};

    //pre-processed version of connectivity_mobile_tx (performed in Python)
    array<array<array<array<int, NUM_BANDS*NUM_MOBILE_TX>, NUM_MOBILE_TX_LOCS>, NUM_RX_LOCS>, NUM_TIME_STEPS> conx_mob_tx_rx = 
        {{
            {{
                {{
                    {0, 0, 0, 0}, {0, 1, 0, 1},
                }},
                {{
                    {1, 0, 1, 0}, {0, 1, 0, 1},
                }},
            }},

            {{
                {{
                    {0, 1, 0, 1}, {0, 1, 0, 1},
                }},
                {{
                    {0, 1, 0, 1}, {1, 1, 1, 1},
                }},
            }},
        }};
    
    //pre-processed version of connectivity_mobile_tx (performed in Python)
    array<array<int, NUM_MOBILE_TX_LOCS>, NUM_MOBILE_TX*NUM_BANDS*NUM_TIME_STEPS> conx_mob_tx_per_num_mob_tx_locs = {{{1, 0}, {0, 1}, {0, 1}, {1, 1}, {1, 0}, {0, 1}, {0, 1}, {1, 1}}};

    array<array<array<impalib_type, NUM_TIME_STEPS>, NUM_BANDS>, NUM_FIXED_TX> fixed_x_eq_const_to_fixed_capac_const_m = 
    {{
        {{
            {{-36.51985024, -57.752808}}, {{-27.23687083, -16.11103224}}
        }},
        {{
            {{-80.8286914, -69.07001696}},{{-67.37688064, -61.80426044}}
        }}
    }};

    array<array<array<impalib_type, NUM_TIME_STEPS>, NUM_BANDS>, NUM_MOBILE_TX> mobile_x_eq_const_to_mobile_capac_const_m = 
    {{
        {{
            {{-13.51566246, -42.2032244}}, {{-95.11148681, -15.40402123}}
        }},
        {{
            {{-87.76378932, -88.95614735}},{{-14.60742991, -68.71767539}}
        }}
    }};

    array<array<impalib_type, NUM_BANDS*NUM_TIME_STEPS>, NUM_MOBILE_TX*NUM_MOBILE_TX_LOCS> r_eq_const_to_auxiliary_const_m = 
                    {{{6.51377099, 0, 0, 6.51377099}, {0, 3.45693075, 3.45693075, 3.45693075}, {-5.45643702, 0, 0, -5.45643702}, {0, -0.41414802, -0.41414802, -0.41414802}}};

    
    impalib_type *fixed_x_eq_const_to_fixed_capac_const_m_pure = addressof(get<0>(fixed_x_eq_const_to_fixed_capac_const_m[0][0]));
    impalib_type *mobile_x_eq_const_to_mobile_capac_const_m_pure = addressof(get<0>(mobile_x_eq_const_to_mobile_capac_const_m[0][0]));
    impalib_type *r_eq_const_to_auxiliary_const_m_pure = addressof(get<0>(r_eq_const_to_auxiliary_const_m[0]));
    const impalib_type *fixed_x_costs_pure = addressof(get<0>(fixed_x_costs[0][0]));
    const impalib_type *mobile_x_costs_pure = addressof(get<0>(mobile_x_costs[0][0]));
    const impalib_type *z_costs_pure = addressof(get<0>(z_costs[0]));
    const impalib_type *r_costs_pure = addressof(get<0>(r_costs[0]));
    const int *connectivity_fixed_tx_pure = addressof(get<0>(connectivity_fixed_tx[0][0][0]));
    const int *connectivity_mobile_tx_pure = addressof(get<0>(connectivity_mobile_tx[0][0][0]));
    const int *fixed_capacity_constraints_pure = fixed_capacity_constraints.data();
    const int *mobile_capacity_constraints_pure = mobile_capacity_constraints.data();
    const int *conx_mob_tx_per_num_mob_tx_locs_pure = addressof(get<0>(conx_mob_tx_per_num_mob_tx_locs[0]));
    const int *conx_fixed_tx_per_num_rx_locs_pure = addressof(get<0>(conx_fixed_tx_per_num_rx_locs[0]));
    const int *conx_mob_tx_rx_pure = addressof(get<0>(conx_mob_tx_rx[0][0][0]));


    model_graph.initialize(fixed_x_eq_const_to_fixed_capac_const_m_pure, mobile_x_eq_const_to_mobile_capac_const_m_pure,
                                                r_eq_const_to_auxiliary_const_m_pure, fixed_x_costs_pure, mobile_x_costs_pure,
                                                z_costs_pure, r_costs_pure, connectivity_fixed_tx_pure,
                                                connectivity_mobile_tx_pure, fixed_capacity_constraints_pure, mobile_capacity_constraints_pure,
                                                conx_mob_tx_per_num_mob_tx_locs_pure, conx_fixed_tx_per_num_rx_locs_pure, conx_mob_tx_rx_pure);

    model_graph.iterate(); 

    for (const auto& val : model_graph.outputs.ExtrinsicFixedX) {
        cout << setw(12) << val << " ";
    }
    cout << "\n--------\n";

    cout << "ExtrinsicMobileX:\n";
    for (const auto& val : model_graph.outputs.ExtrinsicMobileX) {
        cout << setw(12) << val << " ";
    }
    cout << "\n--------\n";

    cout << "ExtrinsicR:\n";
    for (const auto& val : model_graph.outputs.ExtrinsicR) {
        cout << setw(12) << val << " ";
    }
    cout << "\n--------\n";

    cout << "ExtrinsicZ:\n";
    for (const auto& row : model_graph.outputs.ExtrinsicZ) {
        for (const auto& val : row) {
            cout << setw(12) << val << " ";
        }
        cout << "\n";
    }

}