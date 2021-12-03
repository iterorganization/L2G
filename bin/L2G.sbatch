#!/bin/bash

# Example run script for SLURM for submitting on ITER hpc.
# Missing SLURM directives intended!

pwd; hostname; date

# Check if L2G is in the environment
if [ ! -x "$(command -v runL2G)" ]; then
    echo "runL2G not found!"
    echo "L2G is not in the environment!"
    echo "Load the environment by either module load or sourcing the environment script!"
    exit 1
fi
export PYTHONUNBUFFERED=1 # So that we continuously get python output
# $1 is the json file

# -p is for disabling the use of carriage returns when loading meshes.
#    1 - return carriage progress bar
#    0 - normal print progress bar
runL2G -p 0 -i $1

date