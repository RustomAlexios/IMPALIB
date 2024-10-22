# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)
          
import numpy as np
import random
import os
import pickle as pkl

if __name__ == "__main__":
    
    np_impa_lib = np.float32
    zero_value = np_impa_lib(0)
    samples = 100
    
    NUM_FIXED_TX_MAX = 10
    NUM_MOBILE_TX_MAX = 10
    NUM_BANDS_MAX = 10
    NUM_TIME_STEPS_MAX = 10
    NUM_RX_LOCS_MAX = 10
    NUM_MOBILE_TX_LOCS_MAX = 10

    input_file_name = "inputs_mobarp_comparison"

    normal_mean = 0
    normal_variance = 3

    for index_sample in range(0, samples):
        
        NUM_FIXED_TX = np.random.randint(low = 1, high = NUM_FIXED_TX_MAX+1, dtype = int)
        NUM_MOBILE_TX = np.random.randint(low = 1, high = NUM_MOBILE_TX_MAX+1, dtype = int)
        NUM_BANDS = np.random.randint(low = 2, high = NUM_BANDS_MAX+1, dtype = int)
        NUM_TIME_STEPS = np.random.randint(low = 1, high = NUM_TIME_STEPS_MAX+1, dtype = int)
        NUM_RX_LOCS = np.random.randint(low = 1, high = NUM_RX_LOCS_MAX+1, dtype = int)
        NUM_MOBILE_TX_LOCS = np.random.randint(low = 2, high = NUM_MOBILE_TX_LOCS_MAX+1, dtype = int)

        SNR_THRESHOLD = np.inf
        
        EXCLUDE_CAP_FLAG = 0 #np.random.randint(low = 0, high = 2, dtype = int)

        fixed_capac_constraints = np.random.randint(low = 1, high=NUM_BANDS, size=NUM_FIXED_TX, dtype=int)

        mobile_capac_constraints = np.random.randint(low = 1, high=NUM_BANDS, size=NUM_MOBILE_TX, dtype=int)

        connectivity_fixed_tx = np.random.randint(2, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS), dtype=int)

        connectivity_mobile_tx = np.random.randint(2, size=(NUM_MOBILE_TX_LOCS, NUM_BANDS, NUM_TIME_STEPS, NUM_RX_LOCS), dtype=int)

        fixed_x_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_FIXED_TX, NUM_BANDS, NUM_TIME_STEPS))

        mobile_x_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_BANDS, NUM_TIME_STEPS))

        r_costs = np.random.normal(normal_mean, normal_variance, size=(NUM_MOBILE_TX, NUM_MOBILE_TX_LOCS))

        input = [
            NUM_FIXED_TX,
            NUM_MOBILE_TX,
            NUM_BANDS,
            NUM_TIME_STEPS,
            NUM_RX_LOCS, 
            NUM_MOBILE_TX_LOCS,
            fixed_capac_constraints,
            mobile_capac_constraints, 
            connectivity_fixed_tx,
            connectivity_mobile_tx,
            SNR_THRESHOLD, 
            EXCLUDE_CAP_FLAG,
            fixed_x_costs,
            mobile_x_costs,
            r_costs
        ]
        
        if not (os.path.isdir(f"../../../data/{input_file_name}")):
            os.makedirs(f"../../../data/{input_file_name}")

        with open(f"../../../data/{input_file_name}/inputs_set{str(index_sample)}.pkl","wb",) as f:
            print(f"Test file: {index_sample} & NUM_FIXED_TX: {NUM_FIXED_TX} & NUM_MOBILE_TX: {NUM_MOBILE_TX} & NUM_BANDS: {NUM_BANDS} & NUM_TIME_STEPS: {NUM_TIME_STEPS} & NUM_RX_LOCS: {NUM_RX_LOCS} & NUM_MOBILE_TX_LOCS: {NUM_MOBILE_TX_LOCS} & EXCLUDE_CAP_FLAG: {EXCLUDE_CAP_FLAG}")
            pkl.dump(input, f)