# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

from copy import deepcopy
from cmath import inf
#import bitstring
import argparse

import numpy as np
import math
import itertools
from itertools import combinations, permutations
#np.set_printoptions(linewidth=np.inf)

import time
#np.random.seed(17)
from tqdm import tqdm

from multiprocessing import Process, Manager
import multiprocessing as mp

from multiprocessing.sharedctypes import Value, Array
from ctypes import Structure, c_double

import pickle as pkl 
import pdb
#import pandas as pd

np.set_printoptions(threshold=1)
np.set_printoptions(threshold=np.inf)
import statistics
from statistics import mean
import os
import shutil
import re

#np_impa_lib = np.float32
np_impa_lib = np.float64