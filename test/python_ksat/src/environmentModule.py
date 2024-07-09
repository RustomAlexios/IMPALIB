# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import pickle as pkl
import numpy as np
import statistics
import itertools
import argparse
import random
import shutil
import math
import time
import pdb
import re
import os

from collections import (
    defaultdict,
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
from tqdm import (
    tqdm,
)
from copy import (
    deepcopy,
)

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(suppress=True)
# np.set_printoptions(threshold=1)

# np_impa_lib = np.float32
np_impa_lib = np.float64
zero_value = np_impa_lib(0)

#np.random.seed(30)
#random.seed(30)