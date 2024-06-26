# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *
from graphical_model import *

parser = argparse.ArgumentParser()
parser.add_argument(
    "--nITER",
    type=int,
    default=400,
    help="Number of Iterations of IMPA",
)
parser.add_argument(
    "--nNodes",
    type=int,
    default=10,
    help="Number of nodes in the graphical model",
)
parser.add_argument(
    "--threshold",
    type=np_impa_lib,
    default=-0.0001,
    help="Threshold on hard decision",
)
parser.add_argument(
    "--symFlag",
    type=bool,
    default=False,
    help="Symmetric or asymmetric TSP",
)
parser.add_argument(
    "--augmFlag",
    type=bool,
    default=False,
    help="Augmentative MPA or not",
)
parser.add_argument(
    "--exactSolFlag",
    type=bool,
    default=False,
    help="Run exact solver or not",
)
parser.add_argument(
    "--lkhSolFlag",
    type=bool,
    default=True,
    help="Run LKH-heuristic solver or not",
)
parser.add_argument(
    "--simAnnSolFlag",
    type=bool,
    default=False,
    help="Run Simulated Annealing-heuristic solver or not",
)
parser.add_argument(
    "--nSamples",
    type=int,
    default=1000,
    help="Number of Sampled combinations per size for augmentation",
)
parser.add_argument(
    "--resetFlag",
    type=bool,
    default=False,
    help="Number of Sampled combinations per size for augmentation",
)
parser.add_argument(
    "--filteringFlag",
    type=bool,
    default=False,
    help="Activate Filtering or not",
)
parser.add_argument(
    "--alpha",
    type=np_impa_lib,
    default=0.5,
    help="Filtering Rate [0,1]",
)
parser.add_argument(
    "--testFile",
    type=int,
    default=9000,
    help="Test File Index",
)
parser.add_argument(
    "--saveFlag",
    type=bool,
    default=False,
    help="Save Outputs or not",
)
parser.add_argument(
    "--randomTestFlag",
    type=bool,
    default=False,
    help="Generate random test or use test files",
)

if __name__ == "__main__":
    print("MAIN_PYTHON_TSP")

    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER
    NUM_NODES = args.nNodes
    THRESHOLD = args.threshold
    SYMMETRIC_FLAG = args.symFlag
    AUGMENTATION_FLAG = args.augmFlag
    EXACT_SOLVER_FLAG = args.exactSolFlag
    LKH_SOLVER_FLAG = args.lkhSolFlag
    SIM_ANNEALING_FLAG = args.simAnnSolFlag
    NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE = args.nSamples
    RESET_FLAG = args.resetFlag
    FILTERING_FLAG = args.filteringFlag
    ALPHA = args.alpha
    test_file = args.testFile
    SAVE_FLAG = args.saveFlag
    RANDOM_TEST_FLAG = args.randomTestFlag

    formatted_alpha = "{:.1f}".format(ALPHA)

    ModelIMPA = GraphicalModel(
        NUM_ITERATIONS,
        NUM_NODES,
        SYMMETRIC_FLAG,
        AUGMENTATION_FLAG,
        EXACT_SOLVER_FLAG,
        LKH_SOLVER_FLAG,
        SIM_ANNEALING_FLAG,
        NUM_AUGMENTED_SAMPLES_PER_SUBTOUR_SIZE,
        RESET_FLAG,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
    )
    folder_inputs = "../../../src/data/inputs_fixed_1000_nNodes15"
    ModelIMPA.folder_outputs = f"../../../src/data/outputs_fixed_1000_alpha{formatted_alpha}"
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
        with open(
            str(folder_inputs) + "/inputs_set" + str(test_file) + ".pkl",
            "rb",
        ) as f:
            ModelIMPA.input_load = pkl.load(f)

    ModelIMPA.initialize()

    ModelIMPA.start_time = time.time()
    ModelIMPA.iterate_relaxed_graph()
    ModelIMPA.pre_analysis()

    if ModelIMPA.save_flag and not ModelIMPA.random_test_flag:
        ModelIMPA.save_outputs()
