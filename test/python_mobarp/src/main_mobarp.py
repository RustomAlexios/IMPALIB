# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *
from graphical_model import *

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nITER", type=int, default=200, help="Number of Iterations of IMPA")
parser.add_argument("--nFTX", type=int, default=5, help="Number of fixed TX")
parser.add_argument("--nMTX", type=int, default=5, help="Number of mobile TX")
parser.add_argument("--nBands", type=int, default=5, help="Number of bands")
parser.add_argument("--nTimeSteps", type=int, default=5, help="Number of discrete time steps")
parser.add_argument("--nRX", type=int, default=5, help="Number of RX")
parser.add_argument("--nMTXLoc", type=int, default=5, help="Number of mobile TX locations")
parser.add_argument("--filteringFlag", type=bool, default=False, help="Activate Filtering or not")
parser.add_argument("--overWriteCapFlag", type=bool, default=False, help="Over Write Max Capacity or not")
parser.add_argument("--overWriteCapVal", type=int, default=3, help="Over Write Max Capacity Value")
parser.add_argument("--alpha", type=np_impa_lib, default=0.0, help="Filtering Rate [0,1]")
parser.add_argument("--testFile", type=int, default=900000, help="Test File Index")
parser.add_argument("--saveFlag", type=bool, default=False, help="Save Outputs or not")
parser.add_argument("--excludeCapFlag", type=bool, default=False, help="Save Outputs or not")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument("--getSolApproach", type=int, default=2, help="Approach for getting a solution")
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
parser.add_argument("--PPFlag", type=bool, default=False, help="Activate Post-Processing or not")

if __name__ == "__main__":
    print("MAIN_PURE_MOBARP")

    # Parse command-line arguments
    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER # number of iterations of IMPA
    NUM_FIXED_TX = args.nFTX # number of fixed TX
    NUM_MOBILE_TX = args.nMTX # number of mobile TX
    NUM_BANDS = args.nBands # number of bands
    NUM_TIME_STEPS = args.nTimeSteps # number of time steps
    NUM_RX_LOCS = args.nRX # number of RX
    NUM_MOBILE_TX_LOCS = args.nMTXLoc # number of mobile TX locs
    THRESHOLD = args.threshold # threshold for hard decision analysis
    FILTERING_FLAG = args.filteringFlag # perform filtering on messages from degree constraints and subtour elimination constraints
    OVER_WRITE_CAP_FLAG = args.overWriteCapFlag # overwrite capacities flag
    OVER_WRITE_CAP_VAL = args.overWriteCapVal # overwrite capacities value
    ALPHA = args.alpha # filtering parameter
    test_file = args.testFile # test file number example if randomTestFlag is false
    SAVE_FLAG = args.saveFlag # save results to analyze
    RANDOM_TEST_FLAG = args.randomTestFlag # perform random graph analysis
    input_path = args.inputPath # input path of test file if randomTestFlag is false
    output_path = args.outputPath # output path for saving results
    POST_PROCESS_FLAG = args.PPFlag # post processing flag
    EXCLUDE_CAP_FLAG = args.excludeCapFlag # exclude capacity constraints flag
    GET_SOL_APPROACH = args.getSolApproach # approach (for now, 1 or 2) to get solution with minimum objective (2 is more efficient)

    if OVER_WRITE_CAP_FLAG and OVER_WRITE_CAP_VAL >= NUM_BANDS:
        raise ValueError(f"OVER_WRITE_CAP_VAL ({OVER_WRITE_CAP_VAL}) cannot be greater than or equal to NUM_BANDS ({NUM_BANDS})")

    # Format alpha for output folder naming
    formatted_alpha = "{:.1f}".format(ALPHA)

    # Initialize GraphicalModelTsp object with provided parameters
    ModelIMPA = GraphicalModelMOBARP(
        NUM_ITERATIONS,
        NUM_FIXED_TX,
        NUM_MOBILE_TX,
        NUM_BANDS,
        NUM_TIME_STEPS,
        NUM_RX_LOCS,
        NUM_MOBILE_TX_LOCS,
        THRESHOLD,
        FILTERING_FLAG,
        ALPHA,
        RANDOM_TEST_FLAG,
        POST_PROCESS_FLAG,
        OVER_WRITE_CAP_FLAG, 
        OVER_WRITE_CAP_VAL,
        EXCLUDE_CAP_FLAG,
        GET_SOL_APPROACH
    )
    
    # Set folder paths for inputs and outputs
    folder_inputs = "../../../data/" + input_path
    ModelIMPA.formatted_alpha = formatted_alpha
    if FILTERING_FLAG:
        ModelIMPA.folder_outputs = "../../../data/" + output_path + f"_alpha{formatted_alpha}"
    else:
        ModelIMPA.folder_outputs = "../../../data/"+ output_path

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

    # Load input data if not doing random testing
    ModelIMPA.input_load = []
    if not RANDOM_TEST_FLAG:
        ModelIMPA.test_file = test_file
        with open(str(folder_inputs) + "/inputs_set" + str(test_file) + ".pkl", "rb") as f:
            ModelIMPA.input_load = pkl.load(f)
        #snr_fixed = np.load(str(folder_inputs) + "/fixed_test_set" + str(test_file) + ".npy")
        #snr_mobile = np.load(str(folder_inputs) + "/mobile_test_set" + str(test_file) + ".npy")
        #ModelIMPA.snr_fixed = snr_fixed 
        #ModelIMPA.snr_mobile = snr_mobile 

    # Initialize the model
    ModelIMPA.initialize()

    # Start IMPA algorithm and pre-analysis
    ModelIMPA.start_time = time.time()
    
    ModelIMPA.run_impa()

    # Save outputs if saving is enabled and not doing random testing
    if ModelIMPA.save_flag and not ModelIMPA.random_test_flag:
        ModelIMPA.save_outputs()