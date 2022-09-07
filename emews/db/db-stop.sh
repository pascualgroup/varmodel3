#!/bin/bash

# DB STOP SH

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -d || exit 1

if [[ ! -f $DB_DATA/postmaster.pid ]]
then
  echo "Does not exist: $DB_DATA/postmaster.pid"
  echo "Server is probably not running!"
  exit 1
fi

log_hosts "STOP  PID" $( head -1 $DB_DATA/postmaster.pid )
log_hosts "STOP  DB on" $(hostname)

set -x
pg_ctl -D $DB_DATA stop
