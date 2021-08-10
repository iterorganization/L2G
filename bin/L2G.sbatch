#!/bin/bash

# Example SLURM batch script for running on ITER hpc. It takes one command line
# argument, which is the L2G json case file.

#SBATCH --job-name="L2G job"    # Job name
#SBATCH --partition=gen10_ib
#SBATCH --nodes=1                       # Run on one Node
#SBATCH --tasks-per-node=36
#SBATCH --output=l2g.log      # Standard output and error log

# Optional options.
##SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=%USER_ITER_EMAIL%        # Where to send mail
##SBATCH --mem-per-cpu=2gb                        # Job memory request
##SBATCH --time=02:00:00                  # Time limit hrs:min:sec

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

# -e is for outputing the equilibrium slices as EQDSK files in the /current/dir/EQDSKS
#    1 - generates
#    0 - does nothing
# -p is for disabling the use of carriage returns when loading meshes.
#    1 - return carriage progress bar
#    0 - normal print progress bar
runL2G -e 1 -p 0 -i $1

date