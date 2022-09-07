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
echo "UPF:                   $CFG_UPF"
echo "--------------------------"

export QUEUE=$CFG_QUEUE
export WALLTIME=$CFG_WALLTIME
export TURBINE_JOBNAME="${EXPID}_job"
export PROJECT=$CFG_PROJECT

EQ_SQL=$( readlink --canonicalize $EMEWS_PROJECT_ROOT/ext/EQ-SQL )

# if R cannot be found, then these will need to be
# uncommented and set correctly.
# export R_HOME=/path/to/R
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$R_HOME/lib
# if python packages can't be found, then uncommited and set this
# PYTHONPATH="/lcrc/project/EMEWS/bebop/repos/probabilistic-sensitivity-analysis:"
# PYTHONPATH+="/lcrc/project/EMEWS/bebop/repos/panmodel-0.20.0:"
export PYTHONPATH+="$EMEWS_PROJECT_ROOT/python:$EQ_SQL"
export PYTHONHOME="/software/python-anaconda-2021.05-el7-x86_64"
# export PYTHONPATH
# echo "PYTHONPATH: $PYTHONPATH"

export WORKER_POOL_NODES=$CFG_WORKER_POOL_NODES
export WORKER_POOL_PPN=$CFG_WORKER_POOL_PPN

EMEWS_EXT=$EMEWS_PROJECT_ROOT/ext/emews

export DB_HOST=$CFG_DB_HOST
export DB_USER=$CFG_DB_USER
export DB_PORT=$CFG_DB_PORT

export SITE=midway

export ME_NODES=$CFG_SWEEP_NODES
export ME_PPN=$CFG_SWEEP_PPN
export ME_PROCS=$(( ME_NODES * ME_PPN ))

ME_EXPORTS="PYTHONPATH=$EMEWS_PROJECT_ROOT/python:$EMEWS_PROJECT_ROOT/ext/EQ-SQL\n"
export ME_EXPORTS

# set machine to your schedule type (e.g. pbs, slurm, cobalt etc.),
# or empty for an immediate non-queued unscheduled run
MACHINE="slurm-multijob"

export TOTAL_PROCS=$(( WORKER_POOL_NODES * WORKER_POOL_PPN + ME_NODES * ME_PPN ))

# Set these to worker pool values because the worker pool
# runs under turbine and it needs these.
export PPN=$WORKER_POOL_PPN
export PROCS=$(( WORKER_POOL_NODES * WORKER_POOL_PPN ))
MEM_PER_CPU=$(( 56000 / PPN ))

if [ -n "$MACHINE" ]; then
  MACHINE="-m $MACHINE"
fi

mkdir -p $TURBINE_OUTPUT/instances
cp $CFG_FILE $TURBINE_OUTPUT/cfg.cfg

SRC_ME=$EMEWS_PROJECT_ROOT/python/sweep_me.py
TARGET_ME=$TURBINE_OUTPUT/sweep_me.py
cp $SRC_ME $TARGET_ME

SRC_DEFAULT_PARAMS=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_DEFAULT_PARAMS
DST_DEFAULT_PARAMS=$TURBINE_OUTPUT/default_parameters.json
cp $SRC_DEFAULT_PARAMS $DST_DEFAULT_PARAMS

SRC_BRMS_FILE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_BITING_RATE_MULTIPLIERS
DST_BRMS_FILE=$TURBINE_OUTPUT/biting_rate_multipliers.txt
cp $SRC_BRMS_FILE $DST_BRMS_FILE

MEAS_FILE_SOURCE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_MEASUREMENT_FILE
MEAS_FILE=$TURBINE_OUTPUT/measurement.txt
cp $MEAS_FILE_SOURCE $MEAS_FILE

UPF_SOURCE=$EMEWS_PROJECT_ROOT/data/parameters/$CFG_UPF
UPF_TARGET=$TURBINE_OUTPUT/upf.txt
cp $UPF_SOURCE $UPF_TARGET

CFG_EXTRA_FILES_TO_INCLUDE=${CFG_EXTRA_FILES_TO_INCLUDE:-}
for f in ${CFG_EXTRA_FILES_TO_INCLUDE[@]}; do
  tf="$(basename -- $f)"
  cp $EMEWS_PROJECT_ROOT/$f $TURBINE_OUTPUT/$tf
done

VARMODEL_X=$EMEWS_PROJECT_ROOT/../run.jl

CMD_LINE_ARGS="$*  -f=$UPF_TARGET -varmodel_x=$VARMODEL_X "
CMD_LINE_ARGS+="-default_params_file=$DST_DEFAULT_PARAMS "
CMD_LINE_ARGS+="-biting_rate_multiplier_file=$DST_BRMS_FILE -measurement_file=$MEAS_FILE "
CMD_LINE_ARGS+="-result_at=$CFG_RESULT_AT"

# Add any script variables that you want to log as
# part of the experiment meta data to the USER_VARS array,
# for example, USER_VARS=("VAR_1" "VAR_2")
USER_VARS=( )
# log variables and script to to TURBINE_OUTPUT directory

export TURBINE_LAUNCHER=srun
export TURBINE_SBATCH_ARGS="--exclusive\n#SBATCH --mail-type=BEGIN"
export TURBINE_LAUNCH_OPTIONS="--mem-per-cpu=${MEM_PER_CPU} " # \n#SBATCH --exclusive"
export ME_LAUNCH_OPTIONS="--mem-per-cpu=56000 "
ME_COMMAND="python3 $TARGET_ME -u $UPF_TARGET -m $MEAS_FILE "
ME_COMMAND+="-t $CFG_RESULT_AT -e $EXPID"
export ME_COMMAND

PGSQL_LIB=/project2/pascualmm/sfw/gcc-10.2.0/postgreql-14.2/lib

log_script

# echo's anything following this standard out
# set -x

swift-t -n $PROCS $MACHINE -p \
    -r $EQ_SQL -I $EQ_SQL \
    -I $EMEWS_EXT \
    -e EMEWS_PROJECT_ROOT \
    -e SITE \
    -e TURBINE_OUTPUT \
    -e TURBINE_LOG \
    -e TURBINE_DEBUG \
    -e ADLB_DEBUG \
    -e DB_HOST \
    -e DB_USER \
    -e DB_PORT \
    -e PYTHONPATH \
    -e PYTHONHOME \
    -e LD_LIBRARY_PATH=$PGSQL_LIB:$LD_LIBRARY_PATH \
    $EMEWS_PROJECT_ROOT/swift/sweep_worker_pool.swift $CMD_LINE_ARGS

chmod g+rw $TURBINE_OUTPUT/*.tic

# Dry Run Arg:    -t Y: \
