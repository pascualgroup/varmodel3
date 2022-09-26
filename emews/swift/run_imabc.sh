#! /usr/bin/env bash

set -eu

if [ "$#" -ne 2 ]; then
  script_name=$(basename $0)
  echo "Usage: ${script_name} exp_id cfg_file"
  exit 1
fi

export TURBINE_LOG=0 TURBINE_DEBUG=0 ADLB_DEBUG=0
# export TURBINE_STDOUT=out-%%r.txt
export TURBINE_STDOUT=
export ADLB_TRACE=0
export EMEWS_PROJECT_ROOT=$( cd $( dirname $0 )/.. ; /bin/pwd )
# source some utility functions used by EMEWS in this script                                                                                 
source "${EMEWS_PROJECT_ROOT}/etc/emews_utils.sh"

export EXPID=$1
export TURBINE_OUTPUT=$EMEWS_PROJECT_ROOT/experiments/$EXPID
check_directory_exists

CFG_FILE=$2
source $CFG_FILE

echo "--------------------------"
echo "WALLTIME:              $CFG_WALLTIME"
echo "PROCS:                 $CFG_PROCS"
echo "PPN:                   $CFG_PPN"
echo "PRIORS:                $CFG_PRIORS"
echo "TARGETS:               $CFG_TARGETS"
echo "ALGO_PARAMS:           $CFG_ALGO_PARAMS_FILE"
echo "--------------------------"

export PROCS=$CFG_PROCS
export QUEUE=$CFG_QUEUE
export WALLTIME=$CFG_WALLTIME
export PPN=$CFG_PPN
export TURBINE_JOBNAME="${EXPID}_job"
export PROJECT=$CFG_PROJECT

# if R cannot be found, then these will need to be
# uncommented and set correctly.
# export R_HOME=/path/to/R
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$R_HOME/lib
# if python packages can't be found, then uncommited and set this
# PYTHONPATH="/lcrc/project/EMEWS/bebop/repos/probabilistic-sensitivity-analysis:"
# PYTHONPATH+="/lcrc/project/EMEWS/bebop/repos/panmodel-0.20.0:"
export PYTHONPATH+="$EMEWS_PROJECT_ROOT/python"
export PYTHONHOME="/software/python-anaconda-2021.05-el7-x86_64"
# export PYTHONPATH
# echo "PYTHONPATH: $PYTHONPATH"

export SITE=midway

# Resident task workers and ranks
export TURBINE_RESIDENT_WORK_WORKERS=1
export RESIDENT_WORK_RANKS=$(( PROCS - 2 ))

# EQ/R location
EQR=/project2/pascualmm/sfw/gcc-10.2.0/EQ-R


# set machine to your schedule type (e.g. pbs, slurm, cobalt etc.),
# or empty for an immediate non-queued unscheduled run
MACHINE="slurm"

if [ -n "$MACHINE" ]; then
  MACHINE="-m $MACHINE"
fi

mkdir -p $TURBINE_OUTPUT/instances
cp $CFG_FILE $TURBINE_OUTPUT/cfg.cfg

SRC_DEFAULT_PARAMS=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_DEFAULT_PARAMS
DST_DEFAULT_PARAMS=$TURBINE_OUTPUT/default_parameters.json
cp $SRC_DEFAULT_PARAMS $DST_DEFAULT_PARAMS

SRC_BRMS_FILE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_BITING_RATE_MULTIPLIERS
DST_BRMS_FILE=$TURBINE_OUTPUT/biting_rate_multipliers.txt
cp $SRC_BRMS_FILE $DST_BRMS_FILE

PRIORS_FILE_SOURCE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_PRIORS
PRIORS_FILE=$TURBINE_OUTPUT/priors.csv
cp $PRIORS_FILE_SOURCE $PRIORS_FILE

MEAS_FILE_SOURCE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_MEASUREMENT_FILE
MEAS_FILE=$TURBINE_OUTPUT/measurement.txt
cp $MEAS_FILE_SOURCE $MEAS_FILE

SRC_TARGETS_FILE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_TARGETS
DST_TARGETS_FILE=$TURBINE_OUTPUT/targets.csv
cp $SRC_TARGETS_FILE $DST_TARGETS_FILE

ALGO_R_FILE=$EMEWS_PROJECT_ROOT/R/imabc.R
AR_FILE=$TURBINE_OUTPUT/algo.R
cp $ALGO_R_FILE $AR_FILE

ALGO_PARAMS_FILE=$EMEWS_PROJECT_ROOT/data/algo_params/$CFG_ALGO_PARAMS_FILE
AP_FILE=$TURBINE_OUTPUT/algo.params.R
cp $ALGO_PARAMS_FILE $AP_FILE

CFG_EXTRA_FILES_TO_INCLUDE=${CFG_EXTRA_FILES_TO_INCLUDE:-}
for f in ${CFG_EXTRA_FILES_TO_INCLUDE[@]}; do
  tf="$(basename -- $f)"
  cp $EMEWS_PROJECT_ROOT/$f $TURBINE_OUTPUT/$tf
done

VARMODEL_X=$EMEWS_PROJECT_ROOT/../run.jl

CMD_LINE_ARGS="$* -algo_file=$AR_FILE -algo_params_file=$AP_FILE -varmodel_x=$VARMODEL_X "
CMD_LINE_ARGS+="-replicates=$CFG_REPLICATES -default_params_file=$DST_DEFAULT_PARAMS "
CMD_LINE_ARGS+="-biting_rate_multiplier_file=$DST_BRMS_FILE -measurement_file=$MEAS_FILE "
CMD_LINE_ARGS+="-result_at=$CFG_RESULT_AT -prev_result_file=$CFG_PREV_RESULT_FILE"

# Add any script variables that you want to log as
# part of the experiment meta data to the USER_VARS array,
# for example, USER_VARS=("VAR_1" "VAR_2")
USER_VARS=( )
# log variables and script to to TURBINE_OUTPUT directory

# export TURBINE_LAUNCHER=srun
# export TURBINE_SBATCH_ARGS="-c 18"
export TURBINE_SBATCH_ARGS="--mem-per-cpu=${CFG_MEM_PER_CPU}" # \n#SBATCH --exclusive"

log_script

# echo's anything following this standard out
# set -x

swift-t -n $PROCS $MACHINE -p \
    -r $EQR -I $EQR \
    -e EMEWS_PROJECT_ROOT \
    -e SITE \
    -e TURBINE_OUTPUT \
    -e TURBINE_LOG \
    -e TURBINE_DEBUG \
    -e ADLB_DEBUG \
    -e PYTHONPATH \
    -e PYTHONHOME \
    $EMEWS_PROJECT_ROOT/swift/imabc_workflow.swift $CMD_LINE_ARGS

chmod g+rw $TURBINE_OUTPUT/*.tic
