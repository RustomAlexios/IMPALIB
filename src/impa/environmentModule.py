# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from copy import deepcopy
from cmath import inf

# import bitstring
import argparse

import numpy as np
import math
import itertools
from itertools import combinations, permutations, chain, product

np.set_printoptions(linewidth=np.inf)

import time

# np.random.seed(17)
from tqdm import tqdm

from multiprocessing import Process, Manager
import multiprocessing as mp

from multiprocessing.sharedctypes import Value, Array
from ctypes import Structure, c_double

import pickle as pkl
import pdb
# import pandas as pd

np.set_printoptions(threshold=1)
np.set_printoptions(threshold=np.inf)
import statistics
from statistics import mean
import os
import shutil
import re


from collections import defaultdict

np.set_printoptions(suppress=True)

from python_tsp.exact import solve_tsp_dynamic_programming, solve_tsp_brute_force
import elkai
from python_tsp.heuristics import solve_tsp_simulated_annealing
import random
from itertools import islice
import copy

# np_impa_lib = np.float32
np_impa_lib = np.float64

zero_value = np_impa_lib(0)
