#!/bin/bash
set -eu

# DB MK TABLES SH
# Uses workflow.sql to make the DB tables

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -v $*

sql --file $THIS/workflow.sql

echo
echo "db-mk-tables.sh: OK"
