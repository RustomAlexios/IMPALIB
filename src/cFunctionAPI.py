# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ctypes
from environmentModule import *

c_lib = ctypes.CDLL("../build/shared_library/libCfunc.so")
c_int_p = ctypes.POINTER(ctypes.c_int)
c_double_p = ctypes.POINTER(ctypes.c_double)
#c_bool_p = ctypes.POINTER(ctypes.c_bool)
c_float_p = ctypes.POINTER(ctypes.c_float)

c_impa_lib_type_p = c_double_p #c_float_p
#c_impa_lib_type_p = c_float_p
c_impa_lib_type = ctypes.c_double #ctypes_cfloat
#c_impa_lib_type = ctypes.c_float
#np_impa_lib = np.double
#np_impa_lib = np.float32

BcjrWrapper = c_lib.BcjrWrapper

BcjrWrapper.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, c_impa_lib_type_p, c_int_p, c_impa_lib_type_p, \
    c_impa_lib_type_p, c_int_p, c_int_p, c_int_p, ctypes.c_int, c_impa_lib_type_p, c_impa_lib_type_p, ctypes.c_bool, c_impa_lib_type]