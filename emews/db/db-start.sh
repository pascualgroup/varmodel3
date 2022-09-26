#!/bin/bash
set -eu

# DB START SH

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -vd || exit 1

# NOTE: We disabled fsync():
if ( set -x
     nice -n 19 \
          pg_ctl -D $DB_DATA \
                 -l $DB_DATA/db.log \
                 -o "-F -p $DB_PORT" \
                 start
   )
then
  : OK
else
  echo "db-start.sh: Could not start the DB!"
  log_hosts "FAILED TO START DB on" $(hostname) $DB_PORT
  exit 1
fi

log_hosts "START DB on" $(hostname) $DB_PORT
log_hosts "START PID" $( head -1 $DB_DATA/postmaster.pid )

TS=$(date +%s)
FNAME="db_env_vars_${TS}.txt"
echo
echo "Writing DB_SETTINGS to $FNAME"
env | grep DB > $FNAME
