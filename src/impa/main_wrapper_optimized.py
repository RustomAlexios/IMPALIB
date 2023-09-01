#!/usr/bin/env python3
# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from impa.environmentModule import *
from impa.IMPA import *

parser = argparse.ArgumentParser()
parser.add_argument('--nITER', type=int,default=400,help="Number of Iterations of IMPA")
parser.add_argument('--filterFlag', type=bool, default=False, help="Activate Filtering or not")
parser.add_argument('--alpha', type=np_impa_lib, default=0.0, help="Filtering Rate [0,1]")
parser.add_argument('--PPFlag', type=bool, default=True, help="Activate Post-Processing or not")
parser.add_argument('--threshold', type=np_impa_lib, default=-0.0001, help="Threshold on hard decision")
parser.add_argument('--PPOption', type=int, default=1, help="Post-processing technique (1-departments, 2-teams)")
parser.add_argument('--path', metavar='path', default = 'impa', type=str, help='path to directory')

if __name__ == '__main__':
    print('MAIN_WRAPPER_OPTIMIZED')
    index_file = 0
    end_file = index_file+1
    
    args = parser.parse_args()  
    NUM_ITERATIONS = args.nITER
    FILTERING_FLAG = args.filterFlag
    POST_PROCESS_FLAG = args.PPFlag
    ALPHA = args.alpha
    THRESHOLD = args.threshold
    PP_OPTION = args.PPOption
    output_path = args.path
    
    #if not(os.path.exists(f'../impa/{output_path}')):
    #    os.makedirs(f'../impa/{output_path}')
    
    ModelIMPA = IMPA(NUM_ITERATIONS, FILTERING_FLAG, POST_PROCESS_FLAG, ALPHA, THRESHOLD, PP_OPTION)
    
    targs = 500
    #folder_inputs = '../data/targ_vary/inputs_targs_'+str(targs)
    #folder_inputs = '../data/inputs_random_params_1000'
    folder_inputs = '../data/inputs_1000'
    
    for setfile in range(index_file, end_file):
        print('SetFile: ', setfile)
        with open(str(folder_inputs)+'/inputs_set'+str(setfile)+'.pkl', 'rb') as f:
            input_load = pkl.load(f)
        
        ModelIMPA.initialize(input_load)
        start_time = time.time()
        ModelIMPA.iterate()
        runtime = time.time() - start_time
        
        ModelIMPA.pre_analysis()
        
        ModelIMPA.post_analysis()
        #print(ModelIMPA.intrinsic_output)
        print(f'Time: {runtime}')
        #with open(f'../impa/{output_path}/outputs_set{setfile}.pkl', 'wb') as f:
        #    pkl.dump((ModelIMPA.results_composed, runtime), f)
        
        
