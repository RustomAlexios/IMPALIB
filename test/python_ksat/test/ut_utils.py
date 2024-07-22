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

    max_absolute_error = np.max(abs(y_pure - y_wrapper))
    max_relative_error = np.max(abs(y_pure - y_wrapper) / abs(np.max(y_pure) + 1e-30))
    # index = np.argmax(abs(y_pure-y_wrapper)/abs(np.max(y_pure)+1e-30))
    # print('y_pure: ', y_pure[index])
    # print('y_wrapper: ', y_wrapper[index])
    # print('np.max(y_pure)', abs(np.max(y_pure)+1e-30))

    if max_absolute_error > atol or max_relative_error > rtol:
        ut_failed = True
    else:
        ut_failed = False

    with open("ksat.txt", "a") as file:
        if ut_failed:
            output = f"FAILED SUB-TEST {sub_test_num} out of {total_sub_tests}:: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}"
            print(output)
            file.write(output + "\n")
            #print("y_p: ", y_pure)
            #print("y_w: ", y_wrapper)
        else:
            output = f"PASSED SUB-TEST {sub_test_num} out of {total_sub_tests}:: Test Name: {ut_name}, Max. Abs. Error: {max_absolute_error:.4e}, Max. Rel. Error: {max_relative_error:.4e}"
            print(output)
            file.write(output + "\n")
            #print('y_p: ', y_pure)
            #print('y_w: ', y_wrapper)
    file.close()
