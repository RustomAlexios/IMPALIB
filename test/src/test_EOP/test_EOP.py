import numpy as np

np.set_printoptions(threshold=1)
np.set_printoptions(linewidth=np.inf)

N_BINS = 3
N_ITEMS = 100

np_impa_lib = np.float64

extrinsic_output_bin = np.random.uniform(-100, 100, size = (N_BINS, N_ITEMS))
extrinsic_output_bin = extrinsic_output_bin.astype(np_impa_lib)

oric_to_package_m = np.random.uniform(-100, 100, size = N_ITEMS)
oric_to_package_m = oric_to_package_m.astype(np_impa_lib)

f_input1 = 'extrinsic_output_bin.npy'
np.save(f_input1, extrinsic_output_bin)

f_input2 = 'oric_to_package_m.npy'
np.save(f_input2, oric_to_package_m)
    
extrinsic_output_package = np.sum(extrinsic_output_bin,axis=0) + oric_to_package_m
#print(extrinsic_output_bin)
#print(oric_to_package_m)
f_output_path  = 'extrinsic_output_package_pure'
output_file_python_pure = open(f_output_path, 'wb')
np.save(output_file_python_pure, extrinsic_output_package, allow_pickle=True)
output_file_python_pure.close()

extrinsic_output_package_dummy = np.zeros(N_ITEMS, dtype=np_impa_lib)

for bin_index in range(0,N_BINS):
    for item_index in range(0,N_ITEMS):
        extrinsic_output_package_dummy[item_index] = extrinsic_output_package_dummy[item_index] + extrinsic_output_bin[bin_index][item_index]
        
for item_index in range(0,N_ITEMS):
    extrinsic_output_package_dummy[item_index] = extrinsic_output_package_dummy[item_index] + oric_to_package_m[item_index]
    
print(np.array_equal(extrinsic_output_package_dummy, extrinsic_output_package))