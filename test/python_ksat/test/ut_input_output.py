# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import sys

sys.path.append(sys.path[0] + "/../src")

import input_output
from ut_utils import *

def ut_io(ut_name, n_variables, n_constraints):
    
    NUM_VARIABLES = n_variables
    NUM_CONSTRAINTS = n_constraints

    f_input1 = os.getcwd() + "/../ut_inputs/ksat_constraint_to_eq_constraint_m_pure.npy"
    ksat_constraint_to_eq_constraint_m_pure = np.random.uniform(10, 500, size=(NUM_CONSTRAINTS, NUM_VARIABLES))
    
    ksat_constraint_to_eq_constraint_m_pure = ksat_constraint_to_eq_constraint_m_pure.astype(np_impa_lib)
    np.save(f_input1, ksat_constraint_to_eq_constraint_m_pure.flatten())

    outputs = input_output.OutputsImpa(NUM_VARIABLES, NUM_CONSTRAINTS, 0)

    if ut_name == "ExtrinsicOutputVariableEcUpdate":
        extrinsic_output_variable_ec_pure = outputs.extrinsic_output_variable_ec_update(ksat_constraint_to_eq_constraint_m_pure)
        f_output_path = os.getcwd() + "/../ut_results/extrinsic_output_variable_ec_pure"
        output_file_python_pure = open(f_output_path, "wb")
        np.save(output_file_python_pure, extrinsic_output_variable_ec_pure.flatten(), allow_pickle=True)
        output_file_python_pure.close()
