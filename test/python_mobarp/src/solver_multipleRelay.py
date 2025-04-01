from ortools.math_opt.python import mathopt
import pandas as pd
import numpy as np
import os
import pickle as pkl
import time


def create_vars(model_relay, NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS):
    dec_vars_sites = []
    #x_i,k,j
    for i in range(NUM_FIXED_TX):
        dec_vars = []
        for k in range(NUM_TIME_STEPS):
            temp = []
            for j in range(NUM_BANDS):
                temp.append(model_relay.add_binary_variable(name=f"x{i},{k},{j}"))
            dec_vars.append(temp)
        dec_vars_sites.append(dec_vars)
    
    #r_i_prime,n
    dec_vars_relay_locs = []
    for i_prime in range(NUM_MOBILE_TX):
        temp = []
        for n in range(NUM_MOBILE_TX_LOCS):
            temp.append(model_relay.add_binary_variable(name=f"r{i_prime},{n}"))
        dec_vars_relay_locs.append(temp)
    
    #x_i_prime,k,j   
    dec_vars_relay_broadcast = []
    for i_prime in range(NUM_MOBILE_TX):
        dec_vars = []
        for k in range(NUM_TIME_STEPS):
            temp = []
            for j in range(NUM_BANDS):
                temp.append(model_relay.add_binary_variable(name=f"xr{i_prime},{k},{j}"))
            dec_vars.append(temp)
        dec_vars_relay_broadcast.append(dec_vars)
    
    #z_i_prime,n,k,j
    dec_vars_aux = []
    for i_prime in range(NUM_MOBILE_TX):
        dec_vars_actual = []
        for n in range(NUM_MOBILE_TX_LOCS):
            dec_vars = []
            for k in range(NUM_TIME_STEPS):
                temp = []
                for j in range(NUM_BANDS):
                    temp.append(model_relay.add_binary_variable(name=f"z{i_prime},{n},{k},{j}"))
                dec_vars.append(temp)
            dec_vars_actual.append(dec_vars)
        dec_vars_aux.append(dec_vars_actual)
    
    return model_relay, dec_vars_sites, dec_vars_relay_locs, dec_vars_relay_broadcast, dec_vars_aux

def create_set_cover_constraints(model_relay, NUM_RX_LOCS, NUM_TIME_STEPS, NUM_FIXED_TX, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, connectivity_fixed_tx, connectivity_mobile_tx, \
                                dec_vars_sites, dec_vars_aux):
    
    #connectivity_fixed_tx (i, l, k, j)
    #connectivity_mobile_tx (n, l, k, j)
    
    # Covering constraints
    for l in range(NUM_RX_LOCS):
        for k in range(NUM_TIME_STEPS):
            terms = []
            
            for i in range(NUM_FIXED_TX):
                for j in range(NUM_BANDS):
                    terms.append(connectivity_fixed_tx[i, l, k, j] * dec_vars_sites[i][k][j])
            
            for i_prime in range(NUM_MOBILE_TX):
                for n in range(NUM_MOBILE_TX_LOCS):
                    for j in range(NUM_BANDS):
                        terms.append(connectivity_mobile_tx[n, l, k, j] * dec_vars_aux[i_prime][n][k][j])
            model_relay.add_linear_constraint(sum(terms) >= 1)

    return model_relay
   
def create_capacity_constraints(model_relay, NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, dec_vars_sites, dec_vars_relay_broadcast, capacity_fixed, capacity_mobile):
    for i in range(NUM_FIXED_TX):
        for k in range(NUM_TIME_STEPS):
            terms = []
            for j in range(NUM_BANDS):
                terms.append(dec_vars_sites[i][k][j])
            model_relay.add_linear_constraint(sum(terms) <= capacity_fixed[i])
        
    for i_prime in range(NUM_MOBILE_TX): 
        for k in range(NUM_TIME_STEPS):
            terms = []
            for j in range(NUM_BANDS):
                terms.append(dec_vars_relay_broadcast[i_prime][k][j])
            model_relay.add_linear_constraint(sum(terms) <= capacity_mobile[i_prime])
            
    return  model_relay

def create_auxiliary_constraints(model_relay, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, NUM_TIME_STEPS, NUM_BANDS, dec_vars_aux, dec_vars_relay_locs, dec_vars_relay_broadcast):
    # Auxiliary constraints
    for i_prime in range(NUM_MOBILE_TX):
        for n in range(NUM_MOBILE_TX_LOCS):
            for k in range(NUM_TIME_STEPS):
                for j in range(NUM_BANDS):
                    # print(i_prime, n, k, j)
                    model_relay.add_linear_constraint(dec_vars_aux[i_prime][n][k][j] <= dec_vars_relay_locs[i_prime][n]) #z_i_prime,n,k,j <= r_i_prime,n
                    model_relay.add_linear_constraint(dec_vars_aux[i_prime][n][k][j] <= dec_vars_relay_broadcast[i_prime][k][j]) #z_i_prime,n,k,j <= x_i_prime,k,j
                    model_relay.add_linear_constraint(
                        dec_vars_aux[i_prime][n][k][j] >= dec_vars_relay_locs[i_prime][n] + dec_vars_relay_broadcast[i_prime][k][j] - 1) #z_i_prime,n,k,j >= x_i_prime,k,j + r_i_prime,n - 1
    return model_relay

def create_relay_site_selection_constraints(model_relay, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, dec_vars_relay_locs):
    # Relay site selection constraint
    for i_prime in range(NUM_MOBILE_TX):
        terms = []
        for n in range(NUM_MOBILE_TX_LOCS):
            terms.append(dec_vars_relay_locs[i_prime][n])
        model_relay.add_linear_constraint(sum(terms) == 1)
    return model_relay

def get_solution(NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, result, dec_vars_sites, dec_vars_relay_locs, dec_vars_relay_broadcast):
    used_x_list = []
    configurations_activated_r = []
    
    for i in range(NUM_FIXED_TX):
        #print(f"Site {i} assignment:")
        for k in range(NUM_TIME_STEPS):
            for j in range(NUM_BANDS):
                if np.isclose(result.variable_values()[dec_vars_sites[i][k][j]], 1):
                    #print(f"\ttime {k}, band {j}")
                    used_x_list.append(('f', (i, j, k)))
                    
    for i_prime in range(NUM_MOBILE_TX):
        for n in range(NUM_MOBILE_TX_LOCS):
            if np.isclose(result.variable_values()[dec_vars_relay_locs[i_prime][n]], 1):
                # print(f"Relay {i_prime} @ node {n}")
                configurations_activated_r.append((i_prime, n))
                
    for i_prime in range(NUM_MOBILE_TX):
        #print(f"Relay {i_prime} assignment:")
        for k in range(NUM_TIME_STEPS):
            for j in range(NUM_BANDS):
                if np.isclose(result.variable_values()[dec_vars_relay_broadcast[i_prime][k][j]], 1):
                    #print(f"\ttime {k}, band {j}")
                    used_x_list.append(('m', (i_prime, j, k)))
    return used_x_list, configurations_activated_r

if __name__=="__main__":
    # input_file_name = "inputs_mobarp_random_cpsat"
    # output_file_name = "outputs_mobarp_random_cpsat"
    
    # input_file_name = "inputs_mobarp_fixed9_cpsat_im_avg1"
    # output_file_name = "outputs_mobarp_fixed9_cpsat_im_avg1"
    
    snr_threshold = 13
    fixed_size = 20
    
    # input_file_name = f"inputs_mobarp_fixed{fixed_size}_cpsat_snr_threshold{snr_threshold}"
    # output_file_name = f"outputs_mobarp_fixed{fixed_size}_cpsat_snr_threshold{snr_threshold}"
    
    set_number = 4
    type_sim = 'time_analysis'
    input_file_name = f"inputs_mobarp_optimized_random_{type_sim}_cpsat_snr_threshold{snr_threshold}/set{set_number}"
    
    output_file_name = f"outputs_mobarp_optimized_random_{type_sim}_cpsat_snr_threshold{snr_threshold}/set{set_number}"
    
    x = np.load("../../../data/inputs_mobarp_master_real/rx_SNR_relay_master.npy") 
    y = np.load("../../../data/inputs_mobarp_master_real/rx_SNR_master.npy")

    master_connectivity_mobile_tx = np.load("../../../data/inputs_mobarp_master_real/rx_SNR_relay_master.npy") >= snr_threshold
    master_connectivity_fixed_tx = np.load("../../../data/inputs_mobarp_master_real/rx_SNR_master.npy") >= snr_threshold 
    
    # print(f"Pre-sampled master_connectivity_mobile_tx shape: {master_connectivity_mobile_tx.shape}")
    # print(f"Pre-sampled master_connectivity_fixed_tx shape: {master_connectivity_fixed_tx.shape}")

    index_sample = 0
    n_samples = 50
    save_flag = True
    
    criteria_im = 2
    # low_im = 10
    # high_im = 100
    low_im = 10
    high_im = 100
    normal_mean = 0
    normal_variance = 3
    
    selected_indices_list = []
    
    while (index_sample<n_samples):
        print("trying")
        # capacity_fixed = 5
        # capacity_mobile = 5
        
        # min_num_mobile_tx = 2
        # max_num_mobile_tx = 7
        
        # max_band_value = 10
        
        # max_num_time_steps = 12 #24
        
        # #(20, 60), (10, 40), (5, 15)
        # min_num_rx_locs = 5 #10 #20
        # max_num_rx_locs = 15 #40 #60
        
        # #(20, 60), (10, 40), (5, 15)
        # min_num_mobile_tx_locs = 5 #10 #20
        # max_num_mobile_tx_locs = 15 #40 #60
        
        NUM_MOBILE_TX = 6 #np.random.randint(low = min_num_mobile_tx, high = max_num_mobile_tx+1, dtype = int) #np.random.randint(low = 1, high = max_num_mobile_tx+1, dtype = int) # fixed_size 
        NUM_BANDS =  10#2*NUM_MOBILE_TX #np.random.randint(low = 5, high = max_band_value+1, dtype = int) #np.random.randint(low = 2, high = max_band_value+1, dtype = int) # np.minimum(fixed_size, max_band_value)
        NUM_TIME_STEPS = 3*NUM_MOBILE_TX #max_num_time_steps #np.random.randint(low = 1, high = max_num_time_steps+1, dtype = int)  # fixed_size
        NUM_RX_LOCS =  7*NUM_MOBILE_TX #np.random.randint(low = min_num_rx_locs, high = max_num_rx_locs+1, dtype = int) #np.random.randint(low = 1, high = max_num_rx_locs+1, dtype = int) # fixed_size
        NUM_MOBILE_TX_LOCS =  8*NUM_MOBILE_TX #np.random.randint(low = min_num_mobile_tx_locs, high = max_num_mobile_tx_locs+1, dtype = int) #np.random.randint(low = 2, high = max_num_mobile_tx_locs+1, dtype = int) # fixed_size
        
        bands_sel = np.sort(np.random.choice(np.arange(master_connectivity_fixed_tx.shape[3]), NUM_BANDS, replace=False))
        rx_sel = np.sort(np.random.choice(np.arange(master_connectivity_fixed_tx.shape[1]), NUM_RX_LOCS, replace=False))
        timesteps_sel = np.sort(np.random.choice(np.arange(master_connectivity_fixed_tx.shape[2]), NUM_TIME_STEPS, replace=False))
        mobile_tx_locs = np.sort(np.random.choice(np.arange(master_connectivity_mobile_tx.shape[0]), NUM_MOBILE_TX_LOCS, replace=False))
        
        connectivity_mobile_tx = master_connectivity_mobile_tx[mobile_tx_locs][:, rx_sel][:, :, timesteps_sel][:, :, :, bands_sel] #(n, l, k, j)
        connectivity_fixed_tx = master_connectivity_fixed_tx[:, rx_sel][:, :, timesteps_sel][:, :, :, bands_sel] #(i, l, k, j)
        # print(f"Post-sampled connectivity_mobile_tx shape: {connectivity_mobile_tx.shape}")
        # print(f"Post-sampled connectivity_fixed_tx shape: {connectivity_fixed_tx.shape}")
        
        NUM_FIXED_TX = connectivity_fixed_tx.shape[0]  # shore sites
        
        selected_fixed_tx_indices = np.arange(NUM_FIXED_TX)
        selected_bands_indices = bands_sel
        selected_time_steps = timesteps_sel
        selected_rx_locs = rx_sel
        selected_mobile_tx_locs = mobile_tx_locs 
        
        selected_indices = [selected_fixed_tx_indices.tolist(), selected_bands_indices.tolist(), selected_time_steps.tolist(), selected_rx_locs.tolist(), selected_mobile_tx_locs.tolist()]

        if (selected_indices in selected_indices_list):
            print("already used selected_indices")
            exit()
            continue
        
        selected_indices_list.append(selected_indices)

        NUM_RX_LOCS = connectivity_fixed_tx.shape[1]  # receivers
        NUM_TIME_STEPS = connectivity_fixed_tx.shape[2]  # timesteps
        NUM_BANDS = connectivity_fixed_tx.shape[3]  # bands
        NUM_MOBILE_TX_LOCS = connectivity_mobile_tx.shape[0]  # relay locs

        EXCLUDE_CAP_FLAG = 0# np.random.randint(low = 0, high = 2, dtype = int)
        
        if EXCLUDE_CAP_FLAG:
            capacity_fixed = (NUM_BANDS*np.ones(NUM_FIXED_TX, dtype=int)).tolist()
            capacity_mobile = (NUM_BANDS*np.ones(NUM_MOBILE_TX, dtype=int)).tolist()
        else: 
            capacity_fixed = np.random.randint(low = 1, high=NUM_BANDS, size=NUM_FIXED_TX).tolist()
            capacity_mobile = np.random.randint(low = 1, high=NUM_BANDS, size=NUM_MOBILE_TX).tolist()
        
        start_time = time.time()
        model_relay = mathopt.Model(name="Fixed transmitters with taskable mobile transmitters")
        
        model_relay, dec_vars_sites, dec_vars_relay_locs, dec_vars_relay_broadcast, dec_vars_aux = create_vars(model_relay, NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS)

        model_relay = create_set_cover_constraints(model_relay, NUM_RX_LOCS, NUM_TIME_STEPS, NUM_FIXED_TX, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, connectivity_fixed_tx, connectivity_mobile_tx, \
                                        dec_vars_sites, dec_vars_aux)
                         
        model_relay = create_capacity_constraints(model_relay, NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, dec_vars_sites, dec_vars_relay_broadcast, capacity_fixed, capacity_mobile)

        model_relay = create_auxiliary_constraints(model_relay, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, NUM_TIME_STEPS, NUM_BANDS, dec_vars_aux, dec_vars_relay_locs, dec_vars_relay_broadcast)

        model_relay = create_relay_site_selection_constraints(model_relay, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, dec_vars_relay_locs)
            
        # Minimum-order objective
        model_relay.minimize(sum(sum(inner) for outer in dec_vars_sites for inner in outer) +
                            sum(sum(inner) for outer in dec_vars_relay_broadcast for inner in outer))

        params = mathopt.SolveParameters(enable_output=False, solution_pool_size=1)#, threads=1)
        # print(dir(params))
        # print("Number of threads set to:", params.threads)
        # exit()
        result = mathopt.solve(model_relay, mathopt.SolverType.CP_SAT, params=params)
        if result.termination.reason == mathopt.TerminationReason.OPTIMAL:
            print(f"Test file: {index_sample} & NUM_FIXED_TX: {NUM_FIXED_TX} & NUM_MOBILE_TX: {NUM_MOBILE_TX} & NUM_BANDS: {NUM_BANDS} & NUM_TIME_STEPS: {NUM_TIME_STEPS} & NUM_RX_LOCS: {NUM_RX_LOCS} & NUM_MOBILE_TX_LOCS: {NUM_MOBILE_TX_LOCS} & EXCLUDE_CAP_FLAG: {EXCLUDE_CAP_FLAG} & snr_threshold: {snr_threshold}")
            # print("Optimal solution found!")
            end_time = time.time()
            objective_value = result.objective_value()
            cpsat_time = end_time - start_time
            
            print(f"cpsat_time: {cpsat_time}")
            # print(result)
            
            # print("Input parameters:")
            # print(f"\t{NUM_FIXED_TX=}")
            # print(f"\t{NUM_MOBILE_TX=}")
            # print(f"\t{NUM_BANDS=}")
            # print(f"\t{NUM_TIME_STEPS=}")
            # print(f"\t{NUM_RX_LOCS=}")
            # print(f"\t{NUM_MOBILE_TX_LOCS=}")
            print("Objective value:", objective_value)
            
            # print(sum(sum(inner) for outer in dec_vars_sites for inner in outer) +
            #                 sum(sum(inner) for outer in dec_vars_relay_broadcast for inner in outer))

            SNR_THRESHOLD = snr_threshold #because connectivities were already thresholded       

            used_x_list, configurations_activated_r = get_solution(NUM_FIXED_TX, NUM_TIME_STEPS, NUM_BANDS, NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS, result, dec_vars_sites, dec_vars_relay_locs, dec_vars_relay_broadcast)
            # print(f"used_x_list: \n {used_x_list}")
            # print(f"configurations_activated_r: {configurations_activated_r}")
            # exit()
            # print("configurations_activated_r: ", configurations_activated_r)
            # exit()
            reshaped_connectivity_fixed_tx = np.transpose(connectivity_fixed_tx, (0, 3, 2, 1)) #connectivity_fixed_tx (i, l, k, j)
            reshaped_connectivity_mobile_tx = np.transpose(connectivity_mobile_tx, (0, 3, 2, 1)) #connectivity_mobile_tx (n, l, k, j)
            
            if (criteria_im == 1): #normal
                fixed_x_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))
                mobile_x_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))
                r_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS))
            elif (criteria_im == 2): #pos X, normal R
                fixed_x_costs = np.random.uniform(low_im, high_im, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))
                mobile_x_costs = np.random.uniform(low_im, high_im, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))
                r_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS))
            
            input = [
                NUM_FIXED_TX, #0
                NUM_MOBILE_TX, #1
                NUM_BANDS, #2
                NUM_TIME_STEPS, #3
                NUM_RX_LOCS, #4
                NUM_MOBILE_TX_LOCS, #5
                np.array(capacity_fixed, dtype=int), #6
                np.array(capacity_mobile, dtype=int), #7
                reshaped_connectivity_fixed_tx.astype(int), #8
                reshaped_connectivity_mobile_tx.astype(int), #9
                SNR_THRESHOLD, #10
                EXCLUDE_CAP_FLAG, #11
                fixed_x_costs, #12
                mobile_x_costs, #13
                r_costs #14
            ]
            
            if (save_flag):
                if not (os.path.isdir(f"../../../data/{input_file_name}")):
                    os.makedirs(f"../../../data/{input_file_name}")

                with open(f"../../../data/{input_file_name}/inputs_set{str(index_sample)}.pkl","wb",) as f:
                    pkl.dump(input, f)
                
                if not (os.path.isdir(f"../../../data/{output_file_name}")):
                    os.makedirs(f"../../../data/{output_file_name}")
                    
                results_composed = [objective_value, result.termination.reason, used_x_list, configurations_activated_r, cpsat_time]
                
                with open(f"../../../data/{output_file_name}/outputs_set{index_sample}.pkl", "wb") as f:
                    pkl.dump(results_composed, f)
            print('------')
            index_sample +=1
            
            
