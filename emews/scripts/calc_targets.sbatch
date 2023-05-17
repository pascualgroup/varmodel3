#!/bin/bash

#SBATCH --job-name=calc_targets
#SBATCH --account=pi-pascualmm
#SBATCH --output=calc_targets_out.txt
#SBATCH --error=calc_targets_out.txt
#SBATCH --time=36:00:00
#SBATCH --partition=caslake  # accessible partitions listed by the sinfo command
#SBATCH --ntasks-per-node=1  # number of tasks per node
#SBATCH --cpus-per-task=1    # number of CPU cores per task
#SBATCH --exclusive
#SBATCH --mem=0

module load python/anaconda-2020.11

RESULT_AT=30540
EXP_DIR=/project2/pascualmm/ncollier/repos/varmodel3/emews/experiments/lhs_8K_03202023.ah.2.0
python3 /project2/pascualmm/ncollier/repos/varmodel3/emews/python/calc_targets.py $EXP_DIR $RESULT_AT 
