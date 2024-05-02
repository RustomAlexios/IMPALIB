# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

echo "Unit Tests"

declare -i total_sub_tests=1

NUM_VARIABLES=8
NUM_CONSTRAINTS=3
K_VARIABLE=3
FILT_FLAG=true
ALPHA=0.5
threshold=-0.0001
N_ITER=200
TYPE_METRICS=0
PP_ELEMENTS=2
IM_VARIANCE=1

test_counter=1
ut="TEST"

unit_tests=("ExtrinsicOutputVariableEcUpdate" "VariableEc2KsatConstraintUpdate" "KsatConstraint2VariableEcUpdate" "Iterate")

for test_name in ${unit_tests[@]}; do
    mkdir -p ../ut_inputs
    mkdir -p ../ut_results
    echo "$ut: $test_counter"
    export test_name; export NUM_VARIABLES; export K_VARIABLE; export N_ITER; export FILT_FLAG; export NUM_CONSTRAINTS
    for sub_test_number in $( seq 1 $total_sub_tests )
    do
        python3 ../python_ksat/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nVariables=$NUM_VARIABLES --nConstraints=$NUM_CONSTRAINTS --filteringFlag=$FILT_FLAG --kVariable=$K_VARIABLE --alpha=$ALPHA --typeMetrics=$TYPE_METRICS --threshold=$threshold --nITER=$N_ITER --PPElements=$PP_ELEMENTS --var=$IM_VARIANCE
        ../../impalib_unit_tests_ksat
        python3 ../python_ksat/test/ut_methods_utils.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name
    done
    test_counter=$(($test_counter+1))
    rm ../ut_inputs/*
    rm ../ut_results/*
    done
