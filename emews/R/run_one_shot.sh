#! /bin/bash

set -eu

if [ "$#" -ne 1 ]; then
  script_name=$(basename $0)
  echo "Usage: ${script_name} max_iter"
  exit 1
fi


# 4 target files
# 4 results files
GIS=( "on" "off" )
N_TARGETS=( "4" "6" )
TARGETS=( "78" "84" "90" "96" )

TARGET_ROOT="/project/jozik/ncollier/repos/varmodel3/emews/data/parameters/targetsFiles"
PRIORS_ROOT="/project/jozik/ncollier/repos/varmodel3/emews/data/parameters"
OUTPUT_ROOT=/project/jozik/ncollier/repos/varmodel3/emews/scratch
MAX_ITER=$1

for GI in ${GIS[@]}; do
    for N_TARGET in ${N_TARGETS[@]}; do
        for TARGET in ${TARGETS[@]}; do
            TARGET_F=$TARGET_ROOT/targets_id_${TARGET}_11242025.csv
            PRIORS_F=$PRIORS_ROOT/gi_${GI}_priors_03262025.csv
            OUTPUT_D=$OUTPUT_ROOT/gi_${GI}_032025/${N_TARGET}t_${TARGET}
            mkdir -p $OUTPUT_D
            CFG_FILE=$OUTPUT_D/cfg.yaml
            cat > $CFG_FILE << EOF
targets: $TARGET_F
priors: $PRIORS_F
algo_param_file: /project/jozik/ncollier/repos/varmodel3/emews/data/algo_params/imabc_params.R
output_directory: $OUTPUT_D
results_files: 
  - /project/jozik/ncollier/repos/varmodel3/emews/results/SimulatedTargets_gi_${GI}_${N_TARGET}t_LHS2025.csv
n_targets: ${N_TARGET}
end_iter: $MAX_ITER
EOF
            Rscript imabc_one_shot.R $CFG_FILE > $OUTPUT_D/out.txt
            Rscript parms_rds2csv.R $OUTPUT_D/parms_to_run_${MAX_ITER}.RDS
            mv $OUTPUT_D/parms_to_run.csv $OUTPUT_D/params_to_run_gi_${GI}_${N_TARGET}t_${TARGET}_iter_${MAX_ITER}.csv
        done
    done
done





