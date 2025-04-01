from logging import exception
from unittest import result
import numpy as np 
import pandas as pd
import pdb
import os
import pickle as pkl
from solvers import parser
import sys 
from pathlib import Path
import re 
from comparison_utils import generate_comp, compare_shared

in_dir = sys.argv[1]
csvname = sys.argv[2]
or_dir = sys.argv[3]
impa_dir = sys.argv[4]
heur_dir = 'False'


if not(os.path.exists(f'../data/Results')):
    os.makedirs(f'../data/Results')

results_dict = {}
k = 0 
for file in os.listdir(or_dir):
    index = int(re.search(r'\d+', file).group())
    print(index, end='\r')
    infile = in_dir + f'/inputs_set{index}.pkl'
    outfilename = f'outputs_set{index}.pkl'

    cap_fail, target_fail, activated_targets, or_lhs, or_rhs, or_sum, output_or, input, in_size, or_time = generate_comp(infile, or_dir + '/' + outfilename, 'or_tools')
    out_dict = dict(or_lhs=or_lhs, or_rhs=or_rhs, or_sum=or_sum)
    #print(or_sum)

    if not(impa_dir == 'False'): 
        try:
            cap_fail, target_fail, activated_targets, impa_lhs, impa_rhs, impa_sum, output_impa, input, in_size, impa_time = generate_comp(infile, impa_dir + '/' + outfilename, 'impa')
        except:
            print(f'Failure to generate comparison on {infile}')
            continue
        package_types = input[1]
        numshared_assig, percentshared_assig, \
        numshared_pack, percentshared_pack, \
        numshared_targets, percentshared_targets, \
        pt_or_total, pt_or_percent, pt_other_total, pt_other_percent = compare_shared(output_or, output_impa, package_types)
        error_lhs = np.abs((or_lhs-impa_lhs)/or_lhs)
        error_rhs = np.abs((or_rhs-impa_rhs)/or_rhs)
        error_total = np.abs((or_sum-impa_sum)/or_sum)
        if heur_dir == 'False':
            out_dict = dict(cap_fail = cap_fail,
                package_fail = target_fail, 
                numshared_assig = numshared_assig,
                percent_shared_assig = percentshared_assig,
                num_shared_pack = numshared_pack,
                percent_shared_pack = percentshared_pack,
                num_shared_targets = numshared_targets,
                percent_shared_targets = percentshared_targets,
                num_pt_impa = pt_other_total,
                num_pt_or = pt_or_total,
                percent_pt_impa = pt_other_percent,
                percent_pt_or = pt_or_percent,
                impa_lhs = impa_lhs,
                or_lhs = or_lhs,
                p_lhs = impa_lhs/or_lhs,
                p_error_lhs = error_lhs,
                #heur_lhs = h_lhs,
                impa_rhs = impa_rhs,
                or_rhs = or_rhs,
                p_rhs = impa_rhs/or_rhs,
                p_error_rhs = error_rhs,
                #heur_rhs = h_rhs,
                impa_total = impa_sum,
                or_total = or_sum,
                p_total = impa_sum/or_sum,
                p_error_total = error_total,
                #heur_total = h_sum,
                instance_size = in_size,
                impa_time = impa_time,
                or_time = or_time)


    if not(heur_dir == 'False'):
        h_cap_fail, h_target_fail, activated_targets, h_lhs, h_rhs, h_sum, output_heur, input, in_size = generate_comp(infile, heur_dir + '/' + outfilename, 'heur')
        print(h_sum)
        if impa_dir == 'False':
            out_dict = dict(or_lhs = or_lhs,
                heur_lhs = h_lhs,
                or_rhs = or_rhs,
                heur_rhs = h_rhs,
                or_total = or_sum,
                heur_total = h_sum,
                instance_size = in_size)

    if not(heur_dir == 'False') and not(impa_dir == 'False'):
        out_dict = dict(cap_fail = cap_fail,
                package_fail = target_fail, 
                numshared_assig = numshared_assig,
                percent_shared_assig = percentshared_assig,
                num_shared_pack = numshared_pack,
                percent_shared_pack = percentshared_pack,
                num_shared_targets = numshared_targets,
                percent_shared_targets = percentshared_targets,
                num_pt_impa = pt_other_total,
                num_pt_or = pt_or_total,
                percent_pt_impa = pt_other_percent,
                percent_pt_or = pt_other_percent,
                impa_lhs = impa_lhs,
                or_lhs = or_lhs,
                p_error_lhs = error_lhs,
                heur_lhs = h_lhs,
                impa_rhs = impa_rhs,
                or_rhs = or_rhs,
                p_error_rhs = error_rhs,
                heur_rhs = h_rhs,
                impa_total = impa_sum,
                or_total = or_sum,
                p_error_total = error_total,
                heur_total = h_sum,
                instance_size = in_size)
                
    results_dict[index] = out_dict
    if k == 0:
        df = pd.DataFrame(results_dict)
        df = df.transpose()

    results = pd.DataFrame(results_dict)
    results = results.transpose()
    results = results.sort_index()
    #df = pd.concat([df, results], ignore_index=True)
    k +=1 
results.to_csv(f'../data/Results/{csvname}.csv')





        











        
