# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from environmentModule import np

fighters_count = np.array((1, 1, 2, 2, 2))
weasels_count = np.array((0, 1, 0, 1, 2))

F_u = np.array((0, 0, 0, 4, 0, 0, 4, 0, 0, 1, 1), dtype=int)
W_u = np.array((2, 2, 2, 0, 2, 2, 0, 1, 1, 0, 0), dtype=int)

rho = 1
K_1 = 10

distance_metric = rho * np.array(
    (
        [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
        [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
        [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
        [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
        [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
        [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
        [1, 1, 1, 1, 0, 0, 0, 2, 2, 3, 3],
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2],
        [2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2],
    )
)


K_2 = 5
K_3 = 10
alpha = 2
beta = 2


rendezvouz_points = ["CV", "RR"]

# Define 0 for East and 1 West
distance_metric_rendezvouz = np.array(
    (
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
        [
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[1],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
            rendezvouz_points[0],
        ],
    )
)
