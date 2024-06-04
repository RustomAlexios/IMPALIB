# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from copy import deepcopy
from cmath import inf

# import bitstring
import argparse
import time
import numpy as np
import math
import itertools
import pickle as pkl
import pdb
import statistics
import os
import shutil
import re
import elkai
import copy
import random
import sys
from itertools import combinations, permutations, chain, product
from python_tsp.exact import solve_tsp_dynamic_programming, solve_tsp_brute_force
from python_tsp.heuristics import solve_tsp_simulated_annealing
from itertools import islice
from tqdm import tqdm
from ctypes import Structure, c_double
from statistics import mean
from collections import defaultdict

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=1)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(suppress=True)
# np_impa_lib = np.float32
np_impa_lib = np.float64
zero_value = np_impa_lib(0)

#np.random.seed(17)
#random.seed(17)
