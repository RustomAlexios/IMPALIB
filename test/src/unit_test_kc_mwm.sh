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

mkdir -p ../ut_inputs/ut_EqConstraint/ut_TeamEc2OricUpdate
mkdir -p ../ut_inputs/ut_EqConstraint/ut_ProjectEqConst2OricUpdate
mkdir -p ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputTeamUpdate
mkdir -p ../ut_inputs/ut_InputOutput/ut_IntrinsicOutMwmUpdate
mkdir -p ../ut_inputs/ut_Knapsack/ut_ExtrinsicOutputDepartment
mkdir -p ../ut_inputs/ut_Knapsack/ut_KnapsackBackward
mkdir -p ../ut_inputs/ut_Knapsack/ut_KnapsackForward
mkdir -p ../ut_inputs/ut_Knapsack/ut_Team2KnapsackUpdate
mkdir -p ../ut_inputs/ut_Knapsack/ut_ProcessExtrinsicOutputDepartment
mkdir -p ../ut_inputs/ut_ORIC/ut_Oric2TeamUpdate
mkdir -p ../ut_inputs/ut_ORIC/ut_Oric2ProjectEcUpdate
mkdir -p ../ut_inputs/ut_projectIneqConstraint/ut_ProjectInequalityConstraintUpdate
mkdir -p ../ut_inputs/ut_GraphicalModel/ut_Iterate
mkdir -p ../ut_inputs/ut_GraphicalModel/ut_IterateSampleGraph
mkdir -p ../ut_results/ut_EqConstraint/ut_TeamEc2OricUpdate
mkdir -p ../ut_results/ut_EqConstraint/ut_ProjectEqConst2OricUpdate
mkdir -p ../ut_results/ut_InputOutput/ut_ExtrinsicOutputTeamUpdate
mkdir -p ../ut_results/ut_InputOutput/ut_IntrinsicOutMwmUpdate
mkdir -p ../ut_results/ut_Knapsack/ut_ExtrinsicOutputDepartment
mkdir -p ../ut_results/ut_Knapsack/ut_KnapsackBackward
mkdir -p ../ut_results/ut_Knapsack/ut_KnapsackForward
mkdir -p ../ut_results/ut_Knapsack/ut_Team2KnapsackUpdate
mkdir -p ../ut_results/ut_Knapsack/ut_ProcessExtrinsicOutputDepartment
mkdir -p ../ut_results/ut_ORIC/ut_Oric2TeamUpdate
mkdir -p ../ut_results/ut_ORIC/ut_Oric2ProjectEcUpdate
mkdir -p ../ut_results/ut_projectIneqConstraint/ut_ProjectInequalityConstraintUpdate
mkdir -p ../ut_results/ut_GraphicalModel/ut_Iterate
mkdir -p ../ut_results/ut_GraphicalModel/ut_IterateSampleGraph

unit_tests=("Oric2TeamUpdate" "Oric2ProjectEcUpdate" "ExtrinsicOutputTeamUpdate" "IntrinsicOutMwmUpdate" "TeamEc2OricUpdate" "ProjectEqConst2OricUpdate" "ProjectInequalityConstraintUpdate" "KnapsackForward" "KnapsackBackward" "ExtrinsicOutputDepartment" "Team2KnapsackUpdate" "ProcessExtrinsicOutputDepartment" "Iterate" "IterateSampleGraph")

for test_name in ${unit_tests[@]}; do
    echo "$ut: $test_counter"
    export test_name; export N_DEPARTMENTS; export N_TEAMS; export N_PROJECTS; export FILT_FLAG; export Nu; export N_ITER
    for sub_test_number in $( seq 1 $total_sub_tests )
    do
        python3 ../python_kc_mwm/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nDepartments=$N_DEPARTMENTS --nTeams=$N_TEAMS --nProjects=$N_PROJECTS --Nu=$Nu --pTypes=$pTypes --fFlag=$FILT_FLAG --alpha=$ALPHA --nIter=$N_ITER --ppFlag=$POST_PROCESS_FLAG --threshold=$threshold
        ../../impalib_unit_tests_kc_mwm
        python3 ../python_kc_mwm/test/ut_methods_utils.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name
    done
    test_counter=$(($test_counter+1))
    done
#<<comment
rm ../ut_inputs/ut_EqConstraint/ut_TeamEc2OricUpdate/*
rm ../ut_inputs/ut_EqConstraint/ut_ProjectEqConst2OricUpdate/*
rm ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputTeamUpdate/*
rm ../ut_inputs/ut_InputOutput/ut_IntrinsicOutMwmUpdate/*
rm ../ut_inputs/ut_Knapsack/ut_ExtrinsicOutputDepartment/*
rm ../ut_inputs/ut_Knapsack/ut_KnapsackBackward/*
rm ../ut_inputs/ut_Knapsack/ut_KnapsackForward/*
rm ../ut_inputs/ut_Knapsack/ut_Team2KnapsackUpdate/*
rm ../ut_inputs/ut_Knapsack/ut_ProcessExtrinsicOutputDepartment/*
rm ../ut_inputs/ut_ORIC/ut_Oric2TeamUpdate/*
rm ../ut_inputs/ut_ORIC/ut_Oric2ProjectEcUpdate/*
rm ../ut_inputs/ut_projectIneqConstraint/ut_ProjectInequalityConstraintUpdate/*
rm ../ut_inputs/ut_GraphicalModel/ut_Iterate/*
rm ../ut_inputs/ut_GraphicalModel/ut_IterateSampleGraph/*
rm ../ut_results/ut_EqConstraint/ut_TeamEc2OricUpdate/*
rm ../ut_results/ut_EqConstraint/ut_ProjectEqConst2OricUpdate/*
rm ../ut_results/ut_InputOutput/ut_ExtrinsicOutputTeamUpdate/*
rm ../ut_results/ut_InputOutput/ut_IntrinsicOutMwmUpdate/*
rm ../ut_results/ut_Knapsack/ut_ExtrinsicOutputDepartment/*
rm ../ut_results/ut_Knapsack/ut_KnapsackBackward/*
rm ../ut_results/ut_Knapsack/ut_KnapsackForward/*
rm ../ut_results/ut_Knapsack/ut_Team2KnapsackUpdate/*
rm ../ut_results/ut_Knapsack/ut_ProcessExtrinsicOutputDepartment/*
rm ../ut_results/ut_ORIC/ut_Oric2TeamUpdate/*
rm ../ut_results/ut_ORIC/ut_Oric2ProjectEcUpdate/*
rm ../ut_results/ut_projectIneqConstraint/ut_ProjectInequalityConstraintUpdate/*
rm ../ut_results/ut_GraphicalModel/ut_Iterate/*
rm ../ut_results/ut_GraphicalModel/ut_IterateSampleGraph/*
#comment