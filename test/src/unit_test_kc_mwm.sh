# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

echo "Unit Tests"

declare -i total_sub_tests=1

N_DEPARTMENTS=5
Nu="20,30,10,7,8"
N_ITER=400
pTypes="1,2,3,4,5"
N_TEAMS=300
N_PROJECTS=100
FILT_FLAG=true
POST_PROCESS_FLAG=true
ALPHA=0.9
threshold=-0.0001

test_counter=1
ut="TEST"

unit_tests=("Oric2TeamUpdate" "Oric2ProjectEcUpdate" "ExtrinsicOutputTeamUpdate" "IntrinsicOutMwmUpdate" "TeamEc2OricUpdate" "ProjectEqConst2OricUpdate" "ProjectInequalityConstraintUpdate" "KnapsackForward" "KnapsackBackward" "ExtrinsicOutputDepartment" "Team2KnapsackUpdate" "ProcessExtrinsicOutputDepartment" "Iterate" "IterateSampleGraph")

for test_name in ${unit_tests[@]}; do
    mkdir -p ../ut_inputs
    mkdir -p ../ut_results
    echo "$ut: $test_counter"
    export test_name; export N_DEPARTMENTS; export N_TEAMS; export N_PROJECTS; export FILT_FLAG; export Nu; export N_ITER
    for sub_test_number in $( seq 1 $total_sub_tests )
    do
        python3 ../python_kc_mwm/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nDepartments=$N_DEPARTMENTS --nTeams=$N_TEAMS --nProjects=$N_PROJECTS --Nu=$Nu --pTypes=$pTypes --fFlag=$FILT_FLAG --alpha=$ALPHA --nIter=$N_ITER --ppFlag=$POST_PROCESS_FLAG --threshold=$threshold
        ../../impalib_unit_tests_kc_mwm
        python3 ../python_kc_mwm/test/ut_methods_utils.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name
    done
    test_counter=$(($test_counter+1))
    rm ../ut_inputs/*
    rm ../ut_results/*
    done
