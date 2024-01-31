# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import numpy as np
import math
import itertools
import argparse
import time
import pickle as pkl
import pdb
import os
import shutil
import statistics
import re
from copy import deepcopy
from cmath import inf
from itertools import combinations, permutations
from tqdm import tqdm
from ctypes import Structure, c_double
from statistics import mean
# import bitstring

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=np.inf)
# np_impa_lib = np.float32
np_impa_lib = np.float64

zero_value = np_impa_lib(0)
# np.random.seed(17)
