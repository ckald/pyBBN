#!/bin/env bash
##comment out lines by adding at least two `#' at the beginning
#SBATCH --account=boyarsky
#SBATCH --partition=lowmem
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=0-05:00:00
#SBATCH --mem=20000
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ckaldeg@gmail.com,filimonovastasya@gmail.com

################################################################################
# CAUTION: Copy this file to the test folder. To start run using `sbatch`
################################################################################

source ../../../env/bin/activate
cd ../../
PYTHONPATH=. PARALLELIZE=10 python3 ~-/__main__.py

# For future reference
#SBATCH --job-name=
