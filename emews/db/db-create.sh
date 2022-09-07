#!/bin/bash
set -eu

# DB CREATE SH
# Use createdb to create the DB
# The DB must be running for this to work (cf. db-start.sh)

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -v $*

(
  set -x
  createdb -O polaris_user --host=$DB_HOST --port=$DB_PORT $DB_NAME
)

echo "db-create.sh: OK"
