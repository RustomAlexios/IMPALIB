import numpy as np

np.set_printoptions(threshold=1)
np.set_printoptions(linewidth=np.inf)

np_impa_lib = np.float64
f_path_pure  = 'extrinsic_output_package_pure'
f_path_wrapper  = 'extrinsic_output_package_wrapper'

file_array_pure = open(f_path_pure, 'rb');
y_pure = np.load(file_array_pure);

#file_array_wrapper = open(f_path_wrapper, 'rb');
#y_wrapper = np.load(file_array_wrapper, allow_pickle=False);

y_wrapper = np.fromfile(f_path_wrapper, dtype=np_impa_lib);

ut_failed = True
y_pure = y_pure.flatten()
y_wrapper = y_wrapper.flatten()
assert y_pure.shape == y_wrapper.shape, f'Shape mismatch: Python {y_pure.shape} and C++ {y_wrapper.shape}'

if (np.allclose(y_pure, y_wrapper)):
    ut_failed = False

max_error = np.max(abs(y_pure.flatten()-y_wrapper.flatten()))

print(max_error)
print(y_pure.flatten())
print(y_wrapper.flatten())