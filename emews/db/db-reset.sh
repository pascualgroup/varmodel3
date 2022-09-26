#!/bin/bash
set -eu

# DB RESET SH
# Deletes all the table data!
#         and resets the sequence to its start value.
# Environment variables:
# Set DB_CONFIRM=1   to retain    confirmations (default)
# Set DB_CONFIRM=0   to shorten   confirmations
# Set DB_CONFIRM=OFF to eliminate confirmations

if (( ${#} != 0 ))
then
  echo "Too many arguments!"
  exit 1
fi

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -v $*

DB_DELAY=5
if [[ ${DB_CONFIRM:-1} == 1 ]]
then
  echo
  echo "Deleting all table rows ... Enter to confirm ... Ctrl-C to cancel ..."
  read _
  echo "Deleting all table rows ..."
  echo
  # ignore non-zero exit code for no input:
  read -t $DB_DELAY _ || true
else
  echo "Deleting all table rows ..."
  if [[ ${DB_CONFIRM} == "OFF" ]]
  then
    DB_DELAY=0
  else
    DB_DELAY=1
  fi
  # ignore non-zero exit code for no input:
  read -t $DB_DELAY _ || true
fi
echo

sql <<EOF
\set ON_ERROR_STOP on
\dt
select pg_sleep($DB_DELAY);
delete from eq_exp_id_tasks;
delete from eq_tasks;
delete from emews_queue_OUT;
delete from emews_queue_IN;
alter sequence emews_id_generator restart;
EOF

echo "db-reset.sh: OK"
