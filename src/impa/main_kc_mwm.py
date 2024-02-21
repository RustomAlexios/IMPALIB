#!/usr/bin/env python3
# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

# Import necessary modules
from impa.environmentModule import argparse, np_impa_lib, pkl, time
from impa.Impa import ImpaKcMwm

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nITER", type=int, default=400, help="Number of Iterations of IMPA")
parser.add_argument("--filteringFlag", type=bool, default=False, help="Activate Filtering or not")
parser.add_argument("--alpha", type=np_impa_lib, default=0.0, help="Filtering Rate [0,1]")
parser.add_argument("--PPFlag", type=bool, default=True, help="Activate Post-Processing or not")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument(
    "--PPOption",
    type=int,
    default=1,
    help="Post-processing technique (1-departments, 2-teams)",
)
parser.add_argument("--path", metavar="path", default="impa", type=str, help="path to directory")

# Main block
if __name__ == "__main__":
    print("MAIN_WRAPPER_OPTIMIZED")

    # Define index for file iteration
    index_file = 0
    end_file = index_file + 1

    # Parse command-line arguments
    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER # number of iterations of IMPA
    FILTERING_FLAG = args.filteringFlag # filtering flag on messages from knapsack to team constraints
    POST_PROCESS_FLAG = args.PPFlag # post processing flag
    ALPHA = args.alpha # filtering parameter
    THRESHOLD = args.threshold # threshold for hard-decision analysis
    PP_OPTION = args.PPOption # post processing option (on teams or knapsacks)
    output_path = args.path # output path for saving results

    # if not(os.path.exists(f'../impa/{output_path}')):
    #    os.makedirs(f'../impa/{output_path}')

    # Initialize the ImpaKcMwm model with provided parameters
    ModelKcMwm = ImpaKcMwm(NUM_ITERATIONS, FILTERING_FLAG, POST_PROCESS_FLAG, ALPHA, THRESHOLD, PP_OPTION)

    targs = 500  # Example value
    # folder_inputs = '../data/targ_vary/inputs_targs_'+str(targs)
    # folder_inputs = '../data/inputs_random_params_1000'
    folder_inputs = "../../data/inputs_1000"

    for setfile in range(index_file, end_file):
        print("SetFile: ", setfile)

        # Load input data from file
        with open(str(folder_inputs) + "/inputs_set" + str(setfile) + ".pkl", "rb") as f:
            input_load = pkl.load(f)

        # Initialize model with loaded data
        ModelKcMwm.initialize(input_load)
        start_time = time.time()
        # Perform iterations
        ModelKcMwm.iterate()
        runtime = time.time() - start_time

        # Perform pre-analysis
        ModelKcMwm.pre_analysis()

        # Perform post-analysis
        ModelKcMwm.post_analysis()
        # print(ModelKcMwm.intrinsic_output)
        print(f"Time: {runtime}")
        # with open(f'../impa/{output_path}/outputs_set{setfile}.pkl', 'wb') as f:
        #    pkl.dump((ModelKcMwm.results_composed, runtime), f)
