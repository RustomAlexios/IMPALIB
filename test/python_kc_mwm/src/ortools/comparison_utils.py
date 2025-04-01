import numpy as np 
import pickle as pkl 
from solvers import parser
import pdb
import pandas as pd
import time 

def get_caps(results):
    units = set(results[0]).union(set(results[1]))
    temp = {}
    fighters_type = np.array((1, 1, 2, 2, 2))
    weasels_type = np.array((0, 1, 0, 1, 2))
    for unit in units:
        temp[unit] = 0 

    for assignment in results.values: 
        assignment_type = assignment[2] - 1
        temp[assignment[0]] += fighters_type[assignment_type]
        temp[assignment[1]] += weasels_type[assignment_type]
    
    return temp

def get_objs(results, package_costs, target_costs):
    lhs = []
    rhs = [] 
    for assignment in results.values.tolist():
        package = assignment[3]
        target = assignment[4]
        try:
            p_cost = package_costs[package-1]
            t_cost = target_costs[package-1, target-1]

        except:
            pdb.set_trace()
        lhs.append(p_cost)
        rhs.append(t_cost)
    return np.sum(lhs), np.sum(rhs)


def generate_comp(infile, outfile, name):

    with open(infile, 'rb') as f:
        input = pkl.load(f)

    with open(outfile, 'rb') as f:
        output, time = pkl.load(f)

    try:
        N_u, package_types, package_costs, target_costs, l_t, k_t, s_t = input
    
    except: 
        N_u, package_types, package_costs, target_costs = input

    target_costs = np.array(target_costs)
    number_of_targets = target_costs.shape[1]
    
    target_fail = False 
    if name == 'impa':
        output = [(asg[0], asg[1], asg[2], asg[3][0], asg[4]) for asg in output]
        temp = [] 
        target_fail_list = [] 
        for assignment in output:
            
            if len(assignment[4]) > 1:
                target_fail_list.append(assignment)
                target_fail = True
            #print('AHHAHAH')
            for target in assignment[4]:
                new_assignment = list(assignment[0:4])
                new_assignment.append(target)
                temp.append(tuple(new_assignment))

        output = temp

    output = pd.DataFrame(output)
    caps = get_caps(output)

    cap_fail = False
    for unit, usage in caps.items():
        if N_u[unit-1] < usage:
            cap_fail = True

    assignments = output.iloc[:,[0,1,2,4]]
    lhs, rhs = get_objs(output, package_costs, target_costs)
    sum = lhs + rhs
    activated_targets = len(set(output[4]))
    prob_len = target_costs.shape[0] * target_costs.shape[1] 
    return cap_fail, target_fail, activated_targets, lhs, rhs, sum, output, input, prob_len, time

def compare_shared(output_or, output_other, package_types): 
    # Compare assignments 
    assignments_other = output_other.iloc[:,[0,1,2,4]]
    assignments_or = output_or.iloc[:,[0,1,2,4]]
    shared_assig = (np.array(assignments_other)[:, None] == np.array(assignments_or)).all(-1).any(-1)
    numshared_assig = np.sum(shared_assig)
    percentshared_assig = 100*numshared_assig/len(assignments_or)

    # Comparing packages. 
    packages_other = output_other.iloc[:,[0,1,2]]
    packages_or = output_or.iloc[:,[0,1,2]]
    shared_pack = (np.array(packages_other)[:, None] == np.array(packages_or)).all(-1).any(-1)
    numshared_pack = np.sum(shared_pack)
    percentshared_pack = 100*numshared_pack/(len(packages_or))

    # Comparing targets. 
    targets_other = output_other.iloc[:,4].sort_values()
    targets_or = output_or.iloc[:,4].sort_values()

    shared_targets = np.isin(targets_other, targets_or)
    numshared_targets = np.sum(shared_targets)
    percentshared_targets = 100*numshared_targets/len(targets_or)

    # Comparing package types. 
    pt_other_total = [0,0,0,0,0]
    pt_or_total = [0,0,0,0,0]
    pt_other_percent = [0,0,0,0,0]
    pt_or_percent = [0,0,0,0,0]

    pt_other = output_other.iloc[:,2]
    pt_or = output_or.iloc[:,2]
    for pt in package_types:
        pt_other_total[pt-1] = np.sum(pt_other == pt)
        pt_or_total[pt-1] = np.sum(pt_or == pt)
        pt_other_percent[pt-1] = np.sum(pt_other == pt)/len(pt_other)
        pt_or_percent[pt-1] = np.sum(pt_or == pt)/len(pt_or)

    return numshared_assig, percentshared_assig, \
           numshared_pack, percentshared_pack, \
           numshared_targets, percentshared_targets, \
           pt_or_total, pt_or_percent, pt_other_total, pt_other_percent





if __name__ == '__main__':
    cap_fail, target_fail, activated_targets, lhs, rhs, sum, output_or, input, prob_len = generate_comp('Inputs\inputs_100\inputs_set0.pkl', 'Outputs/random_obj_fixed_params/OR-Tools/outputs_100_or/outputs_set0.pkl', 'or_tools')
    #cap_fail, target_fail, activated_targets, lhs, rhs, sum, output_impa, input, in_size = generate_comp('Inputs\inputs_100\inputs_set0.pkl', 'outputs_100_impa\outputs_set0.pkl', 'impa')
    #pt_or_total, pt_or_percent, pt_other_total, pt_other_percent = compare_shared(output_or, output_impa, package_types)
