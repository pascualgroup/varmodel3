#!/bin/bash

set -eu

# Check for an optional timeout threshold in seconds. If the duration of the
# model run as executed below, takes longer that this threshhold
# then the run will be aborted. Note that the "timeout" command
# must be supported by executing OS.

# The timeout argument is optional. By default the "run_model" swift
# app fuction sends 3 arguments, and no timeout value is set. If there
# is a 4th (the TIMEOUT_ARG_INDEX) argument, we use that as the timeout value.

# !!! IF YOU CHANGE THE NUMBER OF ARGUMENTS PASSED TO THIS SCRIPT, YOU MUST
# CHANGE THE TIMEOUT_ARG_INDEX !!!
# TIMEOUT=""
# TIMEOUT_ARG_INDEX=6
# if [[ $# ==  $TIMEOUT_ARG_INDEX ]]
# then
# 	TIMEOUT=${!TIMEOUT_ARG_INDEX}
# fi

# TIMEOUT_CMD=""
# if [ -n "$TIMEOUT" ]; then
#   TIMEOUT_CMD="timeout $TIMEOUT"
# fi

source $EMEWS_PROJECT_ROOT/scripts/${SITE}_env.sh

VARMODEL_X=$1
INSTANCE_DIR=$2

cd $INSTANCE_DIR

arg_array=( "$VARMODEL_X" 
            "parameters.json" )

# MODEL_CMD="python -u ${arg_array[@]}"


julia --check-bounds=no -O3 "${arg_array[@]}"

#$TIMEOUT_CMD $COMMAND
# $? is the exit status of the most recently executed command (i.e the
# line above)
RES=$?
if [ "$RES" -ne 0 ]; then
	if [ "$RES" == 124 ]; then
    echo "---> Timeout error in $COMMAND"
  else
	  echo "---> Error in $COMMAND"
  fi
fi
