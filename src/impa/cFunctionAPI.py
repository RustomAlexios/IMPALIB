# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import ctypes
from impa.environmentModule import os

# Load the C library
c_lib = ctypes.CDLL(os.path.dirname(__file__) + "/lib_wrapper.so")
# c_lib_kc_mwm = ctypes.CDLL(os.path.dirname(__file__) + "/lib_wrapper.so")

# Define pointer types for C data types
c_int_p = ctypes.POINTER(ctypes.c_int)
c_double_p = ctypes.POINTER(ctypes.c_double)
c_bool_p = ctypes.POINTER(ctypes.c_bool)
c_float_p = ctypes.POINTER(ctypes.c_float)

c_impa_lib_type_p = c_double_p  # c_float_p
# c_impa_lib_type_p = c_float_p
c_impa_lib_type = ctypes.c_double  # ctypes_cfloat
# c_impa_lib_type = ctypes.c_float
# np_impa_lib = np.double
# np_impa_lib = np.float32

# Define function prototypes for the TSP wrapper
WrapperTsp = c_lib.WrapperTsp

WrapperTsp.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_bool,
    ctypes.c_bool,
    c_impa_lib_type,
    ctypes.c_bool,
    c_int_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_impa_lib_type,
    c_int_p,
    c_int_p,
    c_impa_lib_type_p,
    c_bool_p,
    c_bool_p,
    c_int_p,
    c_int_p,
    c_bool_p,
    c_int_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
]

# Define function prototypes for the KC-MWM wrapper
WrapperKcMwm = c_lib.WrapperKcMwm

WrapperKcMwm.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    c_impa_lib_type_p,
    c_int_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_int_p,
    c_int_p,
    c_int_p,
    ctypes.c_int,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    ctypes.c_bool,
    c_impa_lib_type,
]

WrapperKsat = c_lib.WrapperKsat

WrapperKsat.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    c_impa_lib_type,
    ctypes.c_bool,
    c_int_p,
    c_int_p,
    c_int_p,
    c_int_p,
    c_int_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    c_impa_lib_type_p,
    ctypes.c_int,
]