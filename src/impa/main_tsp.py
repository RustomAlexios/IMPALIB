# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impa.environmentModule import time, argparse, np_impa_lib, os, pkl
from impa.Impa import GraphicalModelTsp

parser = argparse.ArgumentParser()
parser.add_argument("--nITER", type=int, default=400, help="Number of Iterations of IMPA")
parser.add_argument("--nNodes", type=int, default=10, help="Number of nodes in the graphical model")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument("--symFlag", type=bool, default=False, help="Symmetric or asymmetric TSP")
parser.add_argument("--augmFlag", type=bool, default=False, help="Augmentative MPA or not")
parser.add_argument("--exactSolFlag", type=bool, default=False, help="Run exact solver or not")
parser.add_argument("--lkhSolFlag", type=bool, default=True, help="Run LKH-heuristic solver or not")
parser.add_argument(
    "--simAnnSolFlag",
    type=bool,
    default=False,
    help="Run Simulated Annealing-heuristic solver or not",
)
parser.add_argument(
    "--resetFlag",
    type=bool,
    default=False,
    help="Number of Sampled combinations per size for augmentation",
)
parser.add_argument("--filteringFlag", type=bool, default=False, help="Activate Filtering or not")
parser.add_argument("--alpha", type=np_impa_lib, default=0.7, help="Filtering Rate [0,1]")
parser.add_argument("--testFile", type=int, default=9000, help="Test File Index")
parser.add_argument("--saveFlag", type=bool, default=False, help="Save Outputs or not")
parser.add_argument(
    "--randomTestFlag",
    type=bool,
    default=False,
    help="Generate random test or use test files",
)
parser.add_argument(
    "--inputPath",
    metavar="path",
    default="../../data/impa_input",
    type=str,
    help="path to output directory",
)
parser.add_argument(
    "--outputPath",
    metavar="path",
    default="../../data/impa_output",
    type=str,
    help="path to output directory",
)
parser.add_argument("--maxCount", type=int, default=10, help="Max count of consecutive failures")
parser.add_argument("--PPFlag", type=bool, default=True, help="Activate Post-Processing or not")
parser.add_argument("--KOPTFlag", type=bool, default=False, help="Perform k-opt")
parser.add_argument("--maxAugmCount", type=int, default=200, help="Maximum augmentation count")

if __name__ == "__main__":
    print("MAIN_WRAPPER_TSP")

    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER
    NUM_NODES = args.nNodes
    THRESHOLD = args.threshold
    SYMMETRIC_FLAG = args.symFlag
    AUGMENTATION_FLAG = args.augmFlag
    EXACT_SOLVER_FLAG = args.exactSolFlag
    LKH_SOLVER_FLAG = args.lkhSolFlag
    SIM_ANNEALING_FLAG = args.simAnnSolFlag
    RESET_FLAG = args.resetFlag
    FILTERING_FLAG = args.filteringFlag
    ALPHA = args.alpha
    test_file = args.testFile
    SAVE_FLAG = args.saveFlag
    RANDOM_TEST_FLAG = args.randomTestFlag
    input_path = args.inputPath
    output_path = args.outputPath
    MAX_COUNT = args.maxCount
    POST_PROCESS_FLAG = args.PPFlag
    K_OPT_FLAG = args.KOPTFlag
    MAX_AUGM_COUNT = args.maxAugmCount

    formatted_alpha = "{:.1f}".format(ALPHA)

    ModelIMPA = GraphicalModelTsp(
        NUM_ITERATIONS,
        NUM_NODES,
        SYMMETRIC_FLAG,
        AUGMENTATION_FLAG,
        EXACT_SOLVER_FLAG,
        LKH_SOLVER_FLAG,
        SIM_ANNEALING_FLAG,
        RESET_FLAG,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
        MAX_COUNT,
        POST_PROCESS_FLAG,
        K_OPT_FLAG,
        MAX_AUGM_COUNT,
    )
    folder_inputs = "../../data/" + input_path
    if FILTERING_FLAG:
        ModelIMPA.folder_outputs = "../../data/" + output_path + f"_alpha{formatted_alpha}"
    else:
        ModelIMPA.folder_outputs = "../../data/" + output_path

    if POST_PROCESS_FLAG:
        ModelIMPA.folder_outputs += "_pp"

    ModelIMPA.save_flag = SAVE_FLAG

    if not (os.path.exists(f"{ModelIMPA.folder_outputs}")) and ModelIMPA.save_flag:
        os.makedirs(f"{ModelIMPA.folder_outputs}")

    if not RANDOM_TEST_FLAG:
        print(f"Test File: {test_file}")
    else:
        print("Random testing, no input is used.")

    print(f"Filtering Flag: {FILTERING_FLAG}")

    if FILTERING_FLAG:
        print(f"Alpha: {formatted_alpha}")
    print(f"Reset Flag: {RESET_FLAG}")

    ModelIMPA.input_load = []
    if not RANDOM_TEST_FLAG:
        ModelIMPA.test_file = test_file
        with open(str(folder_inputs) + "/inputs_set" + str(test_file) + ".pkl", "rb") as f:
            ModelIMPA.input_load = pkl.load(f)

    ModelIMPA.initialize()

    ModelIMPA.start_time = time.time()
    ModelIMPA.run_impa()
    ModelIMPA.run_pre_analysis()

    if ModelIMPA.save_flag and not ModelIMPA.random_test_flag:
        ModelIMPA.save_outputs()
