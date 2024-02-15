# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

echo "Unit Tests"

declare -i total_sub_tests=1

N_NODES=9
N_SUBTOURS=50
SYM_FLAG=true
FILT_FLAG=true
AUGM_FLAG=true
ALPHA=0.5
threshold=-0.0001
N_ITER=200

test_counter=1
ut="TEST"

mkdir -p ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputEdgeEcRelaxedGraphUpdate
mkdir -p ../ut_results/ut_InputOutput/ut_ExtrinsicOutputEdgeEcRelaxedGraphUpdate

mkdir -p ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputEdgeEcAugmentedGraphUpdate
mkdir -p ../ut_results/ut_InputOutput/ut_ExtrinsicOutputEdgeEcAugmentedGraphUpdate

mkdir -p ../ut_inputs/ut_DegreeConstraint/ut_DegreeConstraint2EdgeEcUpdate
mkdir -p ../ut_results/ut_DegreeConstraint/ut_DegreeConstraint2EdgeEcUpdate

mkdir -p ../ut_inputs/ut_SubtourConstraint/ut_SubtourConstraints2EdgeEcUpdate
mkdir -p ../ut_results/ut_SubtourConstraint/ut_SubtourConstraints2EdgeEcUpdate

mkdir -p ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintRelaxedGraphUpdate
mkdir -p ../ut_results/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintRelaxedGraphUpdate

mkdir -p ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2SubtourConstraintsUpdate
mkdir -p ../ut_results/ut_EqualityConstraint/ut_EdgeEc2SubtourConstraintsUpdate

mkdir -p ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintAugmentedGraphUpdate
mkdir -p ../ut_results/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintAugmentedGraphUpdate

mkdir -p ../ut_inputs/ut_GraphicalModel/ut_IterateRelaxedGraph
mkdir -p ../ut_results/ut_GraphicalModel/ut_IterateRelaxedGraph

mkdir -p ../ut_inputs/ut_GraphicalModel/ut_IterateAugmentedGraph
mkdir -p ../ut_results/ut_GraphicalModel/ut_IterateAugmentedGraph

unit_tests=("ExtrinsicOutputEdgeEcRelaxedGraphUpdate" "ExtrinsicOutputEdgeEcAugmentedGraphUpdate" "DegreeConstraint2EdgeEcUpdate" "SubtourConstraints2EdgeEcUpdate" "EdgeEc2DegreeConstraintRelaxedGraphUpdate" "EdgeEc2SubtourConstraintsUpdate" "EdgeEc2DegreeConstraintAugmentedGraphUpdate" "IterateRelaxedGraph" "IterateAugmentedGraph")

for test_name in ${unit_tests[@]}; do
    echo "$ut: $test_counter"
    export test_name; export N_NODES; export N_SUBTOURS; export FILT_FLAG; export SYM_FLAG; export N_ITER
    for sub_test_number in $( seq 1 $total_sub_tests )
    do
        if [ "$SYM_FLAG" = true ]; then
            python3 ../python_tsp/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nNodes=$N_NODES --nSubtours=$N_SUBTOURS --fFlag=$FILT_FLAG --symFlag=$SYM_FLAG --alpha=$ALPHA --augmFlag=$AUGM_FLAG --threshold=$threshold --nIter=$N_ITER
        else
            python3 ../python_tsp/test/impalib_unit_tests.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name --nNodes=$N_NODES --nSubtours=$N_SUBTOURS --fFlag=$FILT_FLAG --alpha=$ALPHA --augmFlag=$AUGM_FLAG --threshold=$threshold --nIter=$N_ITER
        fi
        ../../impalib_unit_tests_tsp
        python3 ../python_tsp/test/ut_methods_utils.py --sub_test_num=$sub_test_number --sub_tests_total=$total_sub_tests --ut_name=$test_name
    done
    test_counter=$(($test_counter+1))
    done

#<<'COMMENT'
rm ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputEdgeEcRelaxedGraphUpdate/*
rm ../ut_results/ut_InputOutput/ut_ExtrinsicOutputEdgeEcRelaxedGraphUpdate/*

rm ../ut_inputs/ut_InputOutput/ut_ExtrinsicOutputEdgeEcAugmentedGraphUpdate/*
rm ../ut_results/ut_InputOutput/ut_ExtrinsicOutputEdgeEcAugmentedGraphUpdate/*

rm ../ut_inputs/ut_DegreeConstraint/ut_DegreeConstraint2EdgeEcUpdate/*
rm ../ut_results/ut_DegreeConstraint/ut_DegreeConstraint2EdgeEcUpdate/*

rm ../ut_inputs/ut_SubtourConstraint/ut_SubtourConstraints2EdgeEcUpdate/*
rm ../ut_results/ut_SubtourConstraint/ut_SubtourConstraints2EdgeEcUpdate/*

rm ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintRelaxedGraphUpdate/*
rm ../ut_results/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintRelaxedGraphUpdate/*

rm ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2SubtourConstraintsUpdate/*
rm ../ut_results/ut_EqualityConstraint/ut_EdgeEc2SubtourConstraintsUpdate/*

rm ../ut_inputs/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintAugmentedGraphUpdate/*
rm ../ut_results/ut_EqualityConstraint/ut_EdgeEc2DegreeConstraintAugmentedGraphUpdate/*

rm ../ut_inputs/ut_GraphicalModel/ut_IterateRelaxedGraph/*
rm ../ut_results/ut_GraphicalModel/ut_IterateRelaxedGraph/*

rm ../ut_inputs/ut_GraphicalModel/ut_IterateAugmentedGraph/*
rm ../ut_results/ut_GraphicalModel/ut_IterateAugmentedGraph/*
#COMMENT