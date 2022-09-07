#!/bin/bash
set -eu

# DB INIT SH
# Use pg_ctl to create the whole DB directory "cluster"
# Only used once
# This script insists that the directory not exist in advance
# To delete the "cluster", 'rm -r' the whole directory

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh -dv $*

if [[ -d $DB_DATA ]]
then
  echo "db-init.sh: ERROR: Already exists: DB_DATA=$DB_DATA"
  exit 1
fi

(
  set -x
  pg_ctl -D $DB_DATA initdb
)

echo "db-init.sh: pg_ctl: OK"
echo

# Add the permissive network settings from:
# https://www.alcf.anl.gov/support-center/theta/postgresql-and-sqlite

# COMMENT="# Added by db-init.sh at $( date "+%Y-%m-%d %H:%M" )"

# PG_CONF=$DB_DATA/postgresql.conf
# echo "db-init.sh: editing PG_CONF=$PG_CONF ..."
# cp --verbose --force --backup=numbered $PG_CONF $PG_CONF

# if grep -q "Added by db-init.sh" $PG_CONF
# then
#   echo "db-init.sh: ERROR: PG_CONF was already edited!"
#   exit 1
# fi

# LISTEN_ALL="listen_addresses = '*'"
# # Need many backslashes to escape the comment:
# sed -i "/#listen_addresses/{
# a\\\\$COMMENT
# a$LISTEN_ALL
# }" $PG_CONF
# echo "db-init.sh: edited PG_CONF: OK"

# echo

# HBA_CONF=$DB_DATA/pg_hba.conf
# echo "db-init.sh: editing HBA_CONF=$HBA_CONF ..."

# if grep -q "Added by db-init.sh" $HBA_CONF
# then
#   echo "db-init.sh: ERROR: HBA_CONF was already edited!"
#   exit 1
# fi

# cp --verbose --force --backup=numbered $HBA_CONF $HBA_CONF
# {
#   echo
#   echo $COMMENT
#   echo "host all all 0.0.0.0/0 trust"
# } >> $HBA_CONF
# echo "db-init.sh: edited HBA_CONF: OK"

# Done!
echo
echo "db-init.sh: use db-start.sh to start the DB"
