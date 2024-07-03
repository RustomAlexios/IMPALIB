# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import numpy as np
import multiprocessing as mp
import pickle as pkl
import statistics
import itertools
import argparse
import random
import shutil
import elkai
import time
import math
import pdb
import os
import re
# import bitstring

from multiprocessing.sharedctypes import (
    Value,
    Array,
)
from multiprocessing import (
    Process,
    Manager,
)
from collections import (
    defaultdict,
)
from statistics import (
    mean,
)
from itertools import (
    combinations,
    permutations,
)
from ctypes import (
    Structure,
    c_double,
)
from cmath import (
    inf,
)
from copy import (
    deepcopy,
)
from tqdm import (
    tqdm,
)

# np.set_printoptions(threshold=1)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(suppress=True)

# np_impa_lib = np.float32
np_impa_lib = np.float64
zero_value = np_impa_lib(0)

# np.random.seed(33)
# random.seed(33)