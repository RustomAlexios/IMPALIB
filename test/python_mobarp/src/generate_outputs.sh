#!/bin/bash

percNegIM=60
nSamples=200
alpha=0.9
# inputPath="inputs_mobarp_comparison"
# outputPath="outputs_mobarp_approach_$approach"

# inputPath="inputs_mobarp_fixed13_cpsat_snr_threshold16"
inputPath="inputs_mobarp_random_cpsat_snr_threshold18"
outputPath="outputs_mobarp_random_snr_threshold18_impa_stop_criteria_percNeg${percNegIM}"

# inputPath="inputs_mobarp_fixed5_cpsat_snr_threshold18"
# outputPath="outputs_mobarp_fixed5_snr_threshold18_impa_stop_criteria_percNeg${percNegIM}"

for NUMBER in {0..199}
do
    python3 main_mobarp.py \
        --threshold=-0.0001 \
        --nITER=200 \
        --inputPath="$inputPath"\
        --testFile=$NUMBER \
        --outputPath="$outputPath"\
        --saveFlag=True \
        --percNegIM=$percNegIM \
        --filteringFlag=True \
        --alpha=$alpha
done