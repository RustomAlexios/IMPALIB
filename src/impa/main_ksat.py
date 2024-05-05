# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impa.environmentModule import time, argparse, np_impa_lib, os, pkl, sys
from impa.Impa import GraphicalModelKsat

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nITER", type=int, default=200, help="Number of Iterations of IMPA")
parser.add_argument("--nVariables", type=int, default=10, help="Number of variables in k-SAT")
parser.add_argument("--nConstraints", type=int, default=3, help="Number of constraints in k-SAT")
parser.add_argument("--kVariable", type=int, default=3, help="Number of variables per constraint")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument("--PPElements", type=int, default=2, help="Number of elements in a PP combination")
parser.add_argument("--filteringFlag", type=bool, default=True, help="Activate Filtering or not")
parser.add_argument("--alpha", type=np_impa_lib, default=0.5, help="Filtering Rate [0,1]")
parser.add_argument("--testFile", type=int, default=9000, help="Test File Index")
parser.add_argument("--saveFlag", type=bool, default=False, help="Save Outputs or not")
parser.add_argument("--typeMetrics", type=int, default=0, help="1: Correctly Biased, 0: Random, 2: Zero Mean")
parser.add_argument("--var", type=np_impa_lib, default=3, help="Variance of incoming metrics")
parser.add_argument("--overwrite", type=bool, default=False, help="Overwrite Incoming metrics from file")
parser.add_argument("--PPFlag", type=bool, default=False, help="Activate Post-Processing or not")
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

if __name__ == "__main__":
    print("MAIN_WRAPPER_kSAT")

    # Parse command-line arguments
    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER  # number of iterations of IMPA
    NUM_VARIABLES = args.nVariables  # number of variables of k-SAT
    NUM_CONSTRAINTS = args.nConstraints  # number of constraints of k-SAT
    K_VARIABLE = args.kVariable
    THRESHOLD = args.threshold  # threshold for hard decision analysis
    FILTERING_FLAG = args.filteringFlag  # perform filtering on messages from degree constraints and subtour elimination constraints
    ALPHA = args.alpha  # filtering parameter
    test_file = args.testFile  # test file number example if randomTestFlag is false
    SAVE_FLAG = args.saveFlag  # save results to analyze
    RANDOM_TEST_FLAG = args.randomTestFlag  # perform random graph analysis
    input_path = args.inputPath  # input path of test file if randomTestFlag is false
    output_path = args.outputPath  # output path for saving results
    POST_PROCESS_FLAG = args.PPFlag  # post processing flag
    TYPE_METRICS = args.typeMetrics
    PP_ELEMENTS = args.PPElements
    IM_VARIANCE = args.var
    OVERWRITE = args.overwrite

    if K_VARIABLE >= NUM_VARIABLES:
        print("The value of K_VARIABLE has exceeded or reached the limit of NUM_VARIABLES. Exiting code.")
        sys.exit()

    # Format alpha for output folder naming
    formatted_alpha = "{:.1f}".format(ALPHA)

    # Initialize GraphicalModelTsp object with provided parameters
    ModelIMPA = GraphicalModelKsat(
        NUM_ITERATIONS,
        NUM_VARIABLES,
        NUM_CONSTRAINTS,
        K_VARIABLE,
        FILTERING_FLAG,
        ALPHA,
        THRESHOLD,
        RANDOM_TEST_FLAG,
        POST_PROCESS_FLAG,
        TYPE_METRICS,
        PP_ELEMENTS,
        IM_VARIANCE,
        OVERWRITE,
    )
    # Set folder paths for inputs and outputs
    folder_inputs = "../../data/" + input_path

    print(f"folder_inputs: {folder_inputs}")

    if FILTERING_FLAG:
        ModelIMPA.folder_outputs = "../../data/" + output_path + f"_alpha{formatted_alpha}"
    else:
        ModelIMPA.folder_outputs = "../../data/" + output_path

    if POST_PROCESS_FLAG:
        ModelIMPA.folder_outputs += "_pp"

    ModelIMPA.save_flag = SAVE_FLAG

    # Create output folder if it doesn't exist and saving is enabled
    if not (os.path.exists(f"{ModelIMPA.folder_outputs}")) and ModelIMPA.save_flag:
        os.makedirs(f"{ModelIMPA.folder_outputs}")

    if not RANDOM_TEST_FLAG:
        print(f"Test File: {test_file}")
    else:
        print("Random testing, no input is used.")

    print(f"Filtering Flag: {FILTERING_FLAG}")

    if FILTERING_FLAG:
        print(f"Alpha: {formatted_alpha}")

    # Load input data if not doing random testing
    ModelIMPA.input_load = []
    if not RANDOM_TEST_FLAG:
        ModelIMPA.test_file = test_file
        with open(str(folder_inputs) + "/inputs_set" + str(test_file) + ".pkl", "rb") as f:
            ModelIMPA.input_load = pkl.load(f)

    # Initialize the model
    ModelIMPA.initialize()

    # Start IMPA algorithm and pre-analysis
    ModelIMPA.start_time = time.time()
    ModelIMPA.run_impa()
    # ModelIMPA.run_pre_analysis()

    # Save outputs if saving is enabled and not doing random testing
    if ModelIMPA.save_flag and not ModelIMPA.random_test_flag:
        ModelIMPA.save_outputs()
