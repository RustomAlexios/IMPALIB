# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import pickle as pkl
import numpy as np
import os

np_impa_lib = np.float32
# np_impa_lib = np.float64

zero_value = np_impa_lib(0)

# num_nodes = 15
samples = 1000
type = "random"

for i in range(0, samples):
    # if (i<samples/2):
    #    symmetric_flag = True
    # else:
    #    symmetric_flag = False
    num_nodes = np.random.randint(20, 151)

    symmetric_flag = False

    upper_triangle = np.random.uniform(
        1,
        1000,
        size=(
            num_nodes,
            num_nodes,
        ),
    )

    if symmetric_flag:
        matrix = (
            np.triu(upper_triangle)
            + np.triu(
                upper_triangle,
                k=1,
            ).T
        )
        np.fill_diagonal(
            matrix,
            zero_value,
        )
    else:
        matrix = upper_triangle.copy()
        np.fill_diagonal(
            matrix,
            zero_value,
        )

    edge_connections = np.transpose(np.nonzero(matrix))

    input = [
        num_nodes,
        symmetric_flag,
        matrix,
        edge_connections,
    ]

    # if not(os.path.isdir(f'../../../data/inputs_{type}_{samples}_nNodes{num_nodes}')):
    #        os.makedirs(f'../../../data/inputs_{type}_{samples}_nNodes{num_nodes}')

    # with open(f'../../../data/inputs_{type}_{samples}_nNodes{num_nodes}/inputs_set{str(i)}.pkl', 'wb') as f:
    #    print(f'Test file: {i} & Symmetric Flag: {symmetric_flag}')
    #    pkl.dump(input, f)

    if not (os.path.isdir(f"../../../data/inputs_{type}_{samples}")):
        os.makedirs(f"../../../data/inputs_{type}_{samples}")

    with open(
        f"../../../data/inputs_{type}_{samples}/inputs_set{str(i)}.pkl",
        "wb",
    ) as f:
        print(f"Test file: {i} & num_nodes: {num_nodes} & Symmetric Flag: {symmetric_flag}")
        pkl.dump(input, f)
