# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from itertools import combinations, permutations, chain, product, islice
from python_tsp.heuristics import solve_tsp_simulated_annealing
from python_tsp.exact import solve_tsp_brute_force
from copy import deepcopy
from tqdm import tqdm
from cmath import inf
import pickle as pkl
import numpy as np
import statistics
import itertools
import argparse
import random
import shutil
import elkai
import time
import math
import copy
import pdb
import sys
import re
import os

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=1)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(suppress=True)

# np_impa_lib = np.float32
np_impa_lib = np.float64
zero_value = np_impa_lib(0)

#np.random.seed(17)
#random.seed(17)