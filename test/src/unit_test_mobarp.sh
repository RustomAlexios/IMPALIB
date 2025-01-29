# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

echo "Unit Tests"

declare -i total_sub_tests=1

NUM_FIXED_TX=4
NUM_MOBILE_TX=5
NUM_BANDS=6
NUM_TIME_STEPS=5
NUM_RX_LOCS=7
NUM_MOBILE_TX_LOCS=7
EXCLUDE_CAP_FLAG=0
FILT_FLAG=1 #1 or ''
ALPHA=0.9
threshold=-0.0001
PERCENTAGE_NEGATIVE_IM=10
NUM_ITERATIONS=200

test_counter=1
ut="TEST"

if [ "$EXCLUDE_CAP_FLAG" -eq 1 ]; then
    unit_tests=("ExtrinsicUpdate" "MobileLocEqConst2REqConstUpdate" "SetCoverIneqConstUpdate" "AuxiliaryConst2ZEqConstUpdate" "AuxiliaryConst2REqConstUpdate" "AuxiliaryConst2MobileXEqConstUpdate" "XEqConst2AuxiliaryAndSetCoverConstUpdate" "ReqConstActivation" "ZEqConst2AuxiliaryConstUpdate" "ZEqConst2SetCoverIneqConstUpdate" "Iterate")
else
    unit_tests=("ExtrinsicUpdate" "IneqCapacConstUpdate" "MobileLocEqConst2REqConstUpdate" "SetCoverIneqConstUpdate" "AuxiliaryConst2ZEqConstUpdate" "AuxiliaryConst2REqConstUpdate" "AuxiliaryConst2MobileXEqConstUpdate" "XEqConst2AuxiliaryAndSetCoverConstUpdate" "ReqConstActivation" "ZEqConst2AuxiliaryConstUpdate" "XEqConstActivation" "ZEqConst2SetCoverIneqConstUpdate" "Iterate")
fi

for test_name in ${unit_tests[@]}; do
    mkdir -p ../ut_inputs
    mkdir -p ../ut_results
    echo "$ut: $test_counter"
    export test_name; export NUM_ITERATIONS; export NUM_FIXED_TX; export NUM_MOBILE_TX; export FILT_FLAG; export NUM_BANDS; export NUM_TIME_STEPS; export NUM_RX_LOCS; export NUM_MOBILE_TX_LOCS; export EXCLUDE_CAP_FLAG; export PERCENTAGE_NEGATIVE_IM
    for sub_test_number in $( seq 1 $total_sub_tests )
    do  
        #exclude capacities: EXCLUDE_CAP_FLAG=1
        if [ "$EXCLUDE_CAP_FLAG" -eq 1 ]; then
            python3 ../python_mobarp/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nITER=$NUM_ITERATIONS --nFTX=$NUM_FIXED_TX --nMTX=$NUM_MOBILE_TX --filteringFlag=$FILT_FLAG --nBands=$NUM_BANDS --alpha=$ALPHA --nTimeSteps=$NUM_TIME_STEPS --threshold=$threshold --nMTXLoc=$NUM_MOBILE_TX_LOCS --nRX=$NUM_RX_LOCS --excludeCapFlag=$EXCLUDE_CAP_FLAG --percNegIM=$PERCENTAGE_NEGATIVE_IM
        #do not exclude capacities: EXCLUDE_CAP_FLAG=0
        else
            python3 ../python_mobarp/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nITER=$NUM_ITERATIONS --nFTX=$NUM_FIXED_TX --nMTX=$NUM_MOBILE_TX --filteringFlag=$FILT_FLAG --nBands=$NUM_BANDS --alpha=$ALPHA --nTimeSteps=$NUM_TIME_STEPS --threshold=$threshold --nMTXLoc=$NUM_MOBILE_TX_LOCS --nRX=$NUM_RX_LOCS --percNegIM=$PERCENTAGE_NEGATIVE_IM
        fi 
        ../../impalib_unit_tests_mobarp
        python3 ../python_mobarp/test/ut_methods_utils.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name
    done
    test_counter=$(($test_counter+1))
    rm ../ut_inputs/*
    rm ../ut_results/*
    done
