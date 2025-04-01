from solvers import human, parser, composed
import pickle as pkl 
import numpy as np 
import pdb
import sys
import os 
from comparison_utils import get_objs
import pandas as pd 
import time 
import re 


if __name__ == '__main__':
    #from mac
    inputdir = sys.argv[1]
    outdir = sys.argv[2]
    real = 'False'
    timedict = []

    if not(os.path.exists(f'../data/{outdir}')):
        os.makedirs(f'../data/{outdir}')

    # numbers = []
    # pattern = re.compile(r'(\d+)\.pkl$')

    # for filename in os.listdir(f"../data/{outdir}"):
    #     match = pattern.search(filename)
    #     if match:
    #         numbers.append(int(match.group(1)))
    
    for file in os.listdir(inputdir):
        i = int(re.search(r'\d+', file).group())
        infile = inputdir + f'/{file}'
        print(infile)
        with open(infile, 'rb') as f:
            inputs = pkl.load(f)

        if real == 'True':
            N_u, package_types, package_costs, target_costs, l_t, k_t, s_t = inputs
        
        else:
            N_u, package_types, package_costs, target_costs = inputs

        target_costs = np.array(target_costs)
        number_of_targets = target_costs.shape[1]

        print('number_of_targets: ', number_of_targets)
        print('number_of_items: ', len(package_costs))

        packages_dict, targets_dict, pt_dict, reward_package, reward_target, l_t, k_t, s_t = parser(N_u, package_types, number_of_targets, \
                                                                          package_costs, target_costs)

        #start_time = time.time()
        runtime, results_or = composed(pt_dict, N_u)
        #runtime = time.time() - start_time
        print(f'Time: {runtime}')
        timedict.append(runtime)

        # pd.DataFrame(timedict).to_csv(f'Outputs/{outdir}_time.csv')
        output_or = pd.DataFrame(results_or)
        l_o, r_o = get_objs(output_or, package_costs, target_costs)
        s_o = l_o + r_o 
        print(f'LHS Objective: {l_o}, RHS Objective: {r_o}, Total Objective: {s_o}')

        with open(f'../data/{outdir}/outputs_set{i}.pkl', 'wb') as f:
            #pkl.dump(output_or, f)
            pkl.dump((output_or, runtime), f)
        #exit()
        

