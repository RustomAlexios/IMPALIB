// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#include "../include/impalib/impalib.hpp"


//g++ BCJR_wrapper_optimized.cpp -O3 -march=native -fPIC -shared -o ../build/shared_library/libCfunc.so

//clang++ -std=c++11 -stdlib=libc++ -arch x86_64 BCJR_wrapper_optimized.cpp -shared -o ../build/shared_library/libCfunc.so

//g++ BCJR_wrapper.cpp -lm -O3 -o BCJR


extern "C" void BcjrWrapper(const int NUM_ITERATIONS, const int NUM_DEPARTMENTS, const int NUM_TEAMS, const int NUM_PROJECTS, 
                impalib_type *pTransition_model_py, const int *pMAX_STATE_PY, const impalib_type *pREWARD_TEAM_PY, 
                const impalib_type *pREWARD_PROJECT_PY, const int *pTEAMS_WEIGHTS_PER_DEPARTMENT_PY, const int *p_NON_ZERO_WEIGHT_INDICES_PY, 
                const int *pNON_ZERO_WEIGHT_INDICES_SIZES_PY, const int MAX_SIZE_NON_ZERO_WEIGHTS, impalib_type *pExtrinsic_output_team, impalib_type *pIntrinsic_out_mwm, \
                const bool FILTERING_FLAG, const impalib_type ALPHA){ 

GraphicalModel model_graph(NUM_DEPARTMENTS, NUM_TEAMS, NUM_PROJECTS, MAX_SIZE_NON_ZERO_WEIGHTS, NUM_ITERATIONS, FILTERING_FLAG, ALPHA);

model_graph.initialize(pREWARD_TEAM_PY, pTransition_model_py, pTEAMS_WEIGHTS_PER_DEPARTMENT_PY,
                        pNON_ZERO_WEIGHT_INDICES_SIZES_PY, p_NON_ZERO_WEIGHT_INDICES_PY, pREWARD_PROJECT_PY, 
                        pMAX_STATE_PY);

model_graph.iterate(pNON_ZERO_WEIGHT_INDICES_SIZES_PY);

copy(model_graph.outputs.ExtrinsicOutputTeam.begin(), model_graph.outputs.ExtrinsicOutputTeam.begin() + NUM_TEAMS, pExtrinsic_output_team);
copy(model_graph.outputs.IntrinsicOutMwm.begin(), model_graph.outputs.IntrinsicOutMwm.begin() + NUM_TEAMS*NUM_PROJECTS, pIntrinsic_out_mwm);
}