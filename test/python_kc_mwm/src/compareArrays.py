# Copyright 2023, Alexios Rustom.
# https://github.com/RustomAlexios/IMPALIB
# Distributed under the MIT License.
# (See accompanying LICENSE file or at
#  https://opensource.org/licenses/MIT)

import numpy as np
# import struct
# import bitstring

# differences: 102, 214, 264, 331, 415, 502, 585, 603, 882, 929
# for setfile in range(0,1000):
# print('SetFile: ', setfile)
setfile = 882
number_of_iterations = 400

number_of_items = 1894
number_of_targets = 100
number_of_bins = 11

metric = "IO"

# file_name = 'MEOB'
# output_file_pythonWrapper = open('binaryPythonWrapper'+ str(setfile) +'_Niter'+str(number_of_iterations), 'rb')
# output_file_pythonPureOriginal = open('binaryPythonPure'+ file_name+ str(setfile) +'_Niter'+str(number_of_iterations)+'_Original', 'rb')
# data_pythonPureOriginal = np.load(output_file_pythonPureOriginal)

# output_file_pythonPureIndustry = open('binaryPythonPure'+ file_name+ str(setfile) +'_Niter'+str(number_of_iterations)+'_Industry', 'rb')
# data_pythonPureIndustry = np.load(output_file_pythonPureIndustry)

# with open('binaryPythonPure' + file_name + str(setfile) +'_Niter'+str(number_of_iterations), mode='rb') as file: # b is important -> binary
#    fileContentPure = file.read()

# f1 = bitstring.BitArray(float=data_pythonPure[52], length=64)
# f1.hex

# with open("binaryPythonWrapper" + file_name, mode='rb') as file: # b is important -> binary
#    fileContentWrapper1 = file.read()

# output_original = open('testOriginal', 'rb')
# data_original = np.load(output_original)

# with open("testOriginal", mode='rb') as file: # b is important -> binary
#    fileContentWrapper2 = file.read()

"""
original_pure_iomwm = './binaryPythonPureOriginal/binaryPythonPureIOMWM'+ str(setfile) +'_Niter'+str(number_of_iterations)
optimized_pure_iomwm  = './binaryPythonPureOptimized/binaryPythonPureOptIOMWM'+ str(setfile) +'_Niter'+str(number_of_iterations)

original_pure_io = './binaryPythonPureOriginal/binaryPythonPureIO'+ str(setfile) +'_Niter'+str(number_of_iterations)
optimized_pure_io  = './binaryPythonPureOptimized/binaryPythonPureOptIO'+ str(setfile) +'_Niter'+str(number_of_iterations)

optimized_wrapper_io  = './binaryPythonWrapperOptimized/binaryPythonWrapperOptIO'+ str(setfile) +'_Niter'+str(number_of_iterations)
optimized_wrapper_iomwm  = './binaryPythonWrapperOptimized/binaryPythonWrapperOptIOMWM'+ str(setfile) +'_Niter'+str(number_of_iterations)
"""

# original_pure_iomwm = 'binaryPythonPureIOMWM'+ str(setfile) +'_Niter'+str(number_of_iterations)

optimized_pure_metric = "binaryPythonPureOpt" + metric + str(setfile) + "_Niter" + str(number_of_iterations)

# original_pure_io = 'binaryPythonPureIO'+ str(setfile) +'_Niter'+str(number_of_iterations)
# optimized_pure_io  = 'binaryPythonPureOptIO'+ str(setfile) +'_Niter'+str(number_of_iterations)

# optimized_wrapper_io  = 'binaryPythonWrapperOptIO'+ str(setfile) +'_Niter'+str(number_of_iterations)
optimized_wrapper_metric = "binaryPythonWrapperOpt" + metric + str(setfile) + "_Niter" + str(number_of_iterations)

# data_original = np.fromfile(original, dtype=np.float64) #for opening in C++

# f2 = bitstring.BitArray(float=data_pythonWrapper[52], length=64)
# f2.hex

# output_file_pythonWrapper = open('binaryPythonWrapperIO', 'rb') #for python only
# data_pythonWrapper = np.load(output_file_pythonWrapper) #for python only

# output_file_original_pure_io = open(original_pure_io, 'rb') #for python only
# data_original_pure_io = np.load(output_file_original_pure_io) #for python only
# output_file_optimized_pure_io = open(optimized_pure_io, 'rb') #for python only
# data_optimized_pure_io = np.load(output_file_optimized_pure_io) #for python only

"""
#output_file_original_pure_iomwm = open(original_pure_iomwm, 'rb') #for python only
#data_original_pure_iomwm = np.load(output_file_original_pure_iomwm) #for python only
output_file_optimized_pure_iomwm = open(optimized_pure_iomwm, 'rb') #for python only
data_optimized_pure_iomwm = np.load(output_file_optimized_pure_iomwm) #for python only

output_file_optimized_wrapper_iomwm = open(optimized_wrapper_iomwm, 'rb')
data_optimized_wrapper_iomwm = np.load(output_file_optimized_wrapper_iomwm)
"""

output_file_optimized_pure_metric = open(optimized_pure_metric, "rb")  # for python only
data_optimized_pure_metric = np.load(output_file_optimized_pure_metric)  # for python only

# data_optimized_wrapper_metric   = np.fromfile(optimized_wrapper_metric, dtype=np.float64) #for opening in C++
# data_optimized_wrapper_metric = data_optimized_wrapper_metric.reshape((number_of_bins, number_of_items))
output_file_optimized_wrapper_metric = open(optimized_wrapper_metric, "rb")
data_optimized_wrapper_metric = np.load(output_file_optimized_wrapper_metric)

# output_file_optimized_wrapper_io = open(optimized_wrapper_io, 'rb')
# data_optimized_wrapper_io = np.load(output_file_optimized_wrapper_io)

# data_pythonWrapperACF  = np.fromfile("binaryPythonWrapperMIOACF", dtype=np.float64)

# data = open("sample.bin", "rb").read()

# print(struct.unpack('d', bytearray(data.read(8))))
# print(struct.unpack('d',data[8:16]))

# equal_arrays = np.array_equal(data_pythonWrapper, data_pythonPure)
# print(data_original.shape)
# equal_arrays = np.array_equal(data_original.reshape((number_of_targets, number_of_items)), data_optimized.reshape((number_of_targets, number_of_items)))
# equal_arrays_2 = np.array_equal(data_pythonWrapperACF, data_pythonWrapper)
# equal_arrays = np.array_equal(data_original.reshape((number_of_bins, number_of_items)), data_optimized.reshape((number_of_bins, number_of_items)))

"""
print("--------------------------")
print("Pure: Pure vs. Optimized")
print("--------------------------")
"""

# equal_arrays_pure_io = np.array_equal(data_original_pure_io, data_optimized_pure_io)
# print("equal_arrays_pure_io: ", equal_arrays_pure_io)

# equal_arrays_pure_iomwm = np.array_equal(data_original_pure_iomwm, data_optimized_pure_iomwm)
# print("equal_arrays_pure_iomwm: ", equal_arrays_pure_iomwm)
# print("--------------------------")

"""
print("Optimized: Pure vs. Wrapper")
print("--------------------------")
"""

# equal_arrays_io = np.array_equal(data_optimized_wrapper_io, data_optimized_pure_io)
# print("equal_arrays_io: ", equal_arrays_io)

equal_arrays_metric = np.array_equal(data_optimized_wrapper_metric, data_optimized_pure_metric)
print("equal_arrays_metric: ", equal_arrays_metric)
print("--------------------------")

# if (not equal_arrays_pure_io or not equal_arrays_pure_iomwm):
# exit('Pure: Pure vs. Optimized Violated')

# if (not equal_arrays_io):
#    #print("data_optimized_wrapper_io: ", data_optimized_wrapper_io)
#    #print("data_optimized_pure_io: ", data_optimized_pure_io)
#    #comparison = [i == j for i,j in zip(data_optimized_wrapper_io,data_optimized_pure_io)]
#    #res = [i for i, val in enumerate(comparison) if not val]
#    print('Minimum Difference IO: ',np.min(abs(data_optimized_wrapper_io-data_optimized_pure_io)))
#    print('Maximum Difference IO: ',np.max(abs(data_optimized_wrapper_io-data_optimized_pure_io)))


if not equal_arrays_metric:
    # print("data_optimized_wrapper_io: ", data_optimized_wrapper_iomwm)
    # print("data_optimized_pure_io: ", data_optimized_pure_iomwm)
    print(
        "Minimum Difference Metric: ",
        np.min(np.min(abs(data_optimized_wrapper_metric - data_optimized_pure_metric))),
    )
    print(
        "Maximum Difference Metric: ",
        np.max(np.max(abs(data_optimized_wrapper_metric - data_optimized_pure_metric))),
    )

# print("equal_arrays_2: ", equal_arrays_2)

# if (not equal_arrays_metric):
#    comparison = [i == j for i,j in zip(data_optimized_wrapper_metric,data_optimized_pure_metric)]
#    res = [i for i, val in enumerate(comparison) if not val]
#    #print(max(abs(data_optimized_pure_metric[res]-data_optimized_wrapper_metric[res])))
#    print('indices of difference: ', res)
# print('Python Op: ', max(data_optimized[res]))
# print('Python Or: ', max(data_original[res]))
# print('Length: ', len(res))
