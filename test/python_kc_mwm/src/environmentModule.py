# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from itertools import combinations, permutations
from copy import deepcopy
from tqdm import tqdm
from cmath import inf
import pickle as pkl
import numpy as np
import statistics
import itertools
import argparse
import shutil
import time
import math
import pdb
import sys
import os
import re

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=np.inf)

# np_impa_lib = np.float32
np_impa_lib = np.float64
zero_value = np_impa_lib(0)

# np.random.seed(17)