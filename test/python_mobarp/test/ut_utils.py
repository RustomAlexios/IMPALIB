# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import *

def check_agreement(sub_test_num, total_sub_tests, ut_name, y_pure, y_wrapper, rtol=1e-05, atol=1e-08):
    ut_failed = True
    y_pure = deepcopy(y_pure.flatten())
    y_wrapper = deepcopy(y_wrapper.flatten())
    assert y_pure.shape == y_wrapper.shape, f"Shape mismatch: Python {y_pure.shape} and C++ {y_wrapper.shape}"

    mask = ~((y_pure == -np.inf) & (y_wrapper == -np.inf)) #mask was added since both y_pure and y_wrapper had -inf (due to structure of problem: special case) in "SetCoverIneqConstUpdate" when NUM_FIXED_TX=2, NUM_MOBILE_TX=2, NUM_BANDS=2, NUM_TIME_STEPS=2, NUM_RX_LOCS=1, NUM_MOBILE_TX_LOCS=2, EXCLUDE_CAP_FLAG=0, FILT_FLAG=1 #1 or '', ALPHA=0.9, threshold=-0.0001, PERCENTAGE_NEGATIVE_IM=10, NUM_ITERATIONS=200
    max_absolute_error = np.max(abs(y_pure[mask] - y_wrapper[mask]))
    # max_relative_error = np.max(abs(y_pure - y_wrapper) / abs(np.max(y_pure) + 1e-30))
    max_relative_error = np.max(abs(y_pure[mask] - y_wrapper[mask]) / abs(np.max(y_pure[mask]) + 1e-4))
    
    # max_absolute_error = np.max(abs(y_pure - y_wrapper))
    # # # max_relative_error = np.max(abs(y_pure - y_wrapper) / abs(np.max(y_pure) + 1e-30))
    # max_relative_error = np.max(abs(y_pure - y_wrapper) / abs(np.max(y_pure) + 1e-4))
    
    if max_absolute_error > atol or max_relative_error > rtol:
        ut_failed = True
    else:
        ut_failed = False

    with open("mobarp.txt", "a") as file:
        if ut_failed:
            output = f"FAILED SUB-TEST {sub_test_num} out of {total_sub_tests}:: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}"
            print(output)
            file.write(output + "\n")
            print("y_p: ", y_pure)
            print("y_w: ", y_wrapper)
        else:
            output = f"PASSED SUB-TEST {sub_test_num} out of {total_sub_tests}:: Test Name: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}"
            print(output)
            file.write(output + "\n")
            # print('y_p: ', y_pure)
            # print('y_w: ', y_wrapper)
    file.close()
