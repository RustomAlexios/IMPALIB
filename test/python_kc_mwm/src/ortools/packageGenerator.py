from copy import deepcopy
from cmath import inf
from turtle import distance

import numpy as np
import math
import pdb

#PACKAGES & REWARD GENERATION
fighters_type = np.array((1, 1, 2, 2, 2))
weasels_type = np.array((0, 1, 0, 1, 2))

K_1 = 0
K_3 = 1000


def package_reward_generation(available_combinations, N_u, package_types, rho=1):

    distance_metric = rho*np.array(([0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                            [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                            [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                            [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
                            [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                            [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                            [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
                            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                            [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2],
                            [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2]))


    packages_dict = {} 
    items_weights_per_bin = []
    items_types_per_bin = []
    reward_package = []

    for l in range(0,len(N_u)):
        items_weights_per_bin.append([])
        items_types_per_bin.append([])
    
    combination_package_size = []
    sum_packages = []
    combination_units = []
    combination_package = []

    F_u = np.array((0, 0, 0, 4, 0, 0, 4, 0, 0, 1, 1), dtype = int)
    W_u = np.array((2, 2, 2, 0, 2, 2, 0, 1, 1, 0, 0), dtype = int)

    all_bins  = np.array(range(1,len(N_u)+1))
    index_map = []
    i = 0 
    for u in range(1,len(N_u)+1):
        for v in range(1, len(N_u)+1):
            combination = (u,v)
            if (combination not in available_combinations):
                continue
            #package_array = np.zeros(len(package_types), dtype = int)
            package_array = []
            for type in package_types:
                if (type == 1):
                    if (u == v):
                        package_size = N_u[u-1]
                        for k in range(0,package_size):
                            items_weights_per_bin[u-1].append(fighters_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, v)
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1]))
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1]))
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] - K_1)
                            packages_dict[(u,v,type,i+1)] = (1, 0, distance_metric[v-1][u-1] + F_u[u-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i += 1 
                    else:
                        package_size = 0
                elif (type == 2):
                    if (u == v):
                        package_size = N_u[u-1]/2
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1] + weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, v)
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1]))
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            packages_dict[(u,v,type,i+1)] = (1, 1, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i+=1
                    else:
                        package_size = min(N_u[u-1], N_u[v-1])
                        for k in range(0,package_size):
                            items_weights_per_bin[u-1].append(fighters_type[type-1])
                            items_weights_per_bin[v-1].append(weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            items_types_per_bin[v-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, [u, v])
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1])) 
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1) 
                            packages_dict[(u,v,type,i+1)] = (1, 1, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)  
                            index_map.append((u,v,type,k+1))
                            i+=1                           
                elif (type == 3):
                    if (u == v):
                        package_size = N_u[u-1]/2
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, v)
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1]))
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1]))
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] - K_1)
                            packages_dict[(u,v,type,i+1)] = (1, 1, distance_metric[v-1][u-1] + F_u[u-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i+=1 
                    else:
                        package_size = 0
                elif (type == 4):
                    if (u == v):
                        package_size = N_u[u-1]/3
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1]+ weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, v)
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1])) 
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1) 
                            packages_dict[(u,v,type,i+1)] = (2, 1, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)  
                            index_map.append((u,v,type,k+1))
                            i+=1       
                    else:
                        package_size = min(N_u[u-1]/2, N_u[v-1])
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1])
                            items_weights_per_bin[v-1].append(weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            items_types_per_bin[v-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, [u, v])
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1])) 
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            packages_dict[(u,v,type,i+1)] = (2, 1, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i+=1 
                elif (type == 5):
                    if (u == v):
                        package_size = N_u[u-1]/4 
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1]+ weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)   
                            remaining_bins = np.setdiff1d(all_bins, v)
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1]))
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1) 
                            packages_dict[(u,v,type,i+1)] = (2, 2, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i+=1            
                    else:
                        package_size = min(N_u[u-1]/2, N_u[v-1]/2)
                        for k in range(0,math.floor(package_size)):
                            items_weights_per_bin[u-1].append(fighters_type[type-1])
                            items_weights_per_bin[v-1].append(weasels_type[type-1])
                            items_types_per_bin[u-1].append(type)
                            items_types_per_bin[v-1].append(type)
                            remaining_bins = np.setdiff1d(all_bins, [u, v])
                            for l in remaining_bins:
                                items_weights_per_bin[l-1].append(0)
                                items_types_per_bin[l-1].append(type)  
                            #reward_package.append(-K_1*math.exp(-distance_metric[v-1][u-1] - F_u[u-1] - W_u[v-1]))
                            reward_package.append(distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            packages_dict[(u,v,type,i+1)] = (2, 2, distance_metric[v-1][u-1] + F_u[u-1] + W_u[v-1] - K_1)
                            index_map.append((u,v,type,k+1))
                            i+=1 
                if (type in range(1,5+1)):
                    package_size = math.floor(package_size)
                    #package_array[type-1] = package_size
                    package_array.append(package_size)
            combination_package_size.append(np.sum(package_array)) 
            combination_package.append(package_array)
            sum_packages.append(sum(package_array))
            combination_units.append((u,v))
            #print(' u: ', u, '  v: ', v, '  Packages Size:  ', package_array, 'with sum: ', sum(package_array), '& last index: ', np.sum(sum_packages)-1) 

    return packages_dict, reward_package, items_weights_per_bin, items_types_per_bin, combination_units, \
    combination_package_size, combination_package, sum_packages, index_map


def package_to_target_rewards(last_index, combination_units, items_types, number_of_targets, rendezvouz_points, alpha=2, beta=2):
    
    number_of_packages = len(items_types)
    priority_targets = np.array(list(range(1,number_of_targets+1)))
    P_target = alpha*priority_targets
    D_package_target  = np.zeros((number_of_packages , number_of_targets))
    K_package_target = np.zeros((number_of_packages , number_of_targets))
    S_package_target = np.zeros((number_of_packages , number_of_targets))
    reward_target = np.zeros((number_of_packages , number_of_targets))

    #np.random.seed(100) # Delete later
    l_target = np.random.uniform(0, 1, number_of_targets)
    k_target = np.random.randint(min(fighters_type), max(fighters_type)+1, number_of_targets)
    s_target = np.random.randint(min(weasels_type),max(weasels_type)+1 , number_of_targets)



    distance_metric_rendezvouz = np.array(([rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0],rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0],rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0],rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0],rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[1], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]],
                                            [rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[1], rendezvouz_points[1], \
                                                rendezvouz_points[1], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0], rendezvouz_points[0]]))

    for i in range(0, number_of_packages):
        index = i
        #print('-------')
        #print('Package: ', index)
        location_package  = next(x[0] for x in enumerate(last_index) if x[1] >= index)
        units = combination_units[location_package]
        fighters_unit = units[0]
        weasels_unit = units[1]
        #print('fighters_unit: ', fighters_unit, 'weasels_unit: ', weasels_unit)
        for j in range(0,number_of_targets):
            if (distance_metric_rendezvouz[weasels_unit-1][fighters_unit-1] == rendezvouz_points[0]): #East Rendezvouz
                D_package_target[i,j] = beta*(1-l_target[j])
            elif (distance_metric_rendezvouz[weasels_unit-1][fighters_unit-1] == rendezvouz_points[1]): #West Rendezvouz
                D_package_target[i,j] = beta*l_target[j]
                #print('D_package_target: \n', D_package_target)
            if (items_types[index] == 1):
                f_p = 1; w_p = 0 
            elif (items_types[index] == 2):
                f_p = 1; w_p = 1
            elif (items_types[index] == 3): 
                f_p = 2; w_p = 0
            elif (items_types[index] == 4): 
                f_p = 2; w_p = 1
            elif (items_types[index] == 5): 
                f_p = 2; w_p = 2
            #print('f_p: ', f_p, 'w_p: ', w_p)
            K_package_target[i,j] = K_3*max(k_target[j]-f_p,0)
            #print('K_package_target: \n', K_package_target)
            S_package_target[i,j] = K_3*max(s_target[j]-w_p,0) 
            #print('S_package_target: \n', S_package_target)
            #reward_target[i,j] = K_2*math.exp(-P_target[j]-D_package_target[i,j]\
            #-K_package_target[i,j]-S_package_target[i,j])
            reward_target[i,j] = (P_target[j] + D_package_target[i,j]\
            + K_package_target[i,j] + S_package_target[i,j])
    return reward_target, l_target, k_target, s_target

'''
#Failure Set 1 EXCEEDING CAPACITY
#N_u = np.array((20,10), dtype = int)
#package_types = [4]
reward_package = np.array((-95.76636096, -49.87315831, -5.81116413, -30.68214207, -55.29666546,
                            -47.97233243, -13.23280859, -58.61278067, -85.1950681, -54.43160649,
                            -2.82328854, -29.40863014, -52.60179514, -56.10501736, -13.51966673, 
                            -12.43682453, -43.63107882, -95.74963336, -85.20236773, -42.01690685,
                            -78.13231954, -85.87414019, -72.78118831, -13.56985145))
'''

'''
#Failure Set 2 EXCEEDING CAPACITY
#N_u = np.array((20,10), dtype = int)
package_types = [4]
reward_package = np.array((-19.2314878, -13.59366678, -32.46565792, -13.2121892, -38.95372275,
                            -23.56397947, -22.73242766, -5.08886663, -37.9914633, -6.88904801,
                            -18.26892251, -23.87887244, -17.32570221, -42.93494897, -7.03703582,
                            -13.47742751, -29.30734718, -13.30069987, -45.69877777, -35.60170799,
                            -47.11943431, -5.64879255, -32.61117326, -4.56240884))
'''

'''
#Failure Set 3 BELOW CAPACITY
#N_u = np.array((20,10), dtype = int)
#package_types = [4]
reward_package = np.array((-21.20028248, -39.98253874, -15.24504692, -49.18181998, -25.71815791,
                            -44.24706712, -2.87870096, -33.84159098, -37.32201823, -17.22807573,
                            -27.52110849, -40.54436894, -3.71761878, -1.25336952, -45.379379,
                            -43.01789557, -40.06256304, -8.47026759, -34.32956508, -27.62230759,
                            -24.52062999, -11.82766268, -28.22580518, -14.07774011))
'''

'''
#Failure Set 4 EXCEEDING & BELOW CAPACITY
#N_u = np.array((20,10,10,4), dtype = int)
#package_types = [4]
reward_package = np.array((-22.78904233, -18.73342656, -5.83030304, -18.66100268, -1.88063542,
                            -15.22815364, -30.48562483, -25.26897682, -42.47984982, -47.7415164,
                            -15.73670265, -44.48415638, -25.2776985, -31.46045618, -35.44103014,
                            -46.97737577, -4.98577863, -2.77448358, -18.43864212, -13.85231812,
                            -41.76140327, -29.90405845, -42.58719704, -16.53038141, -23.15037636,
                            -48.77388907, -34.41675911, -46.93765022, -36.41251917, -11.45982193,
                            -41.84399401, -15.43142097, -3.88323574, -26.28799021, -47.64708855,
                            -41.15927458, -12.74540504, -27.2123815, -16.94303437, -35.3858961,
                            -32.90364703, -27.05213609, -24.02238212, -27.59280921, -18.51252675,
                            -48.4670415, -36.93336504, -38.60249075, -17.18634782, -16.69555278,
                            -11.67736148, -49.44374619, -8.40594946, -5.66471829, -20.33723588,
                            -41.34947988, -22.54867504, -43.51113001, -41.61537011, -32.49487297,
                            -11.47340542, -1.07161417, -16.2116931, -31.4247045, -4.95583063))
'''

