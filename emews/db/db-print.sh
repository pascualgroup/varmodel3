#!/bin/bash
set -eu

# DB PRINT SH
# Just print the tables for human inspection

if (( ${#} == 0 ))
then
  :
elif (( ${#} == 1 ))
then
  DB_NAME=$1
else
  echo "Too many arguments!"
  exit 1
fi

THIS=$( readlink --canonicalize $( dirname $0 ) )
source $THIS/db-settings.sh || exit 1

sql < $THIS/db-print.sql
