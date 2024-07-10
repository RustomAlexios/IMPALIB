# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import pkl, argparse, np_impa_lib
from graphical_model import GraphicalModelKcMwm

parser = argparse.ArgumentParser()
parser.add_argument("--nITER", type=int, default=400, help="Number of Iterations of IMPA")
parser.add_argument("--filteringFlag", type=bool, default=False, help="Activate Filtering or not")
parser.add_argument("--alpha", type=np_impa_lib, default=0.0, help="Filtering Rate [0,1]")
parser.add_argument("--PPFlag", type=bool, default=True, help="Activate Post-Processing or not")
parser.add_argument("--threshold", type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")

if __name__ == "__main__":
    print("MAIN_PYTHON_PURE_OPTIMIZED")

    targs = 500
    # folder_inputs = '../targ_vary/inputs_targs_'+str(targs)
    # folder_inputs = '../../../data/inputs_random_params_1000'
    folder_inputs = "../../../data/inputs_1000"

    index_file = 0
    end_file = index_file + 1

    args = parser.parse_args()
    NUM_ITERATIONS = args.nITER
    FILTERING_FLAG = args.filteringFlag
    POST_PROCESS_FLAG = args.PPFlag
    ALPHA = args.alpha
    THRESHOLD = args.threshold

    ModelIMPA = GraphicalModelKcMwm(NUM_ITERATIONS, FILTERING_FLAG, POST_PROCESS_FLAG, ALPHA, THRESHOLD)

    for setfile in range(index_file, end_file):
        with open(str(folder_inputs) + "/inputs_set" + str(setfile) + ".pkl", "rb") as f:
            input_load = pkl.load(f)

        ModelIMPA.initialize(input_load)

        ModelIMPA.iterate()

        ModelIMPA.pre_analysis()

        ModelIMPA.post_analysis()
