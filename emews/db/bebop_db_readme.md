# Bebop DB Notes

Postgres Executable is installed in `/lcrc/project/EMEWS/bebop/sfw/gcc-7.1.0/postgres-14.2/bin`

## Server Configuration

`env-bebop.sh` sets all the bebop specified variables. This
should be sourced before starting the database.
DB_HOST will be set to the current hostname (e.g., beboplogin1, etc.).

`export DB_HOST=${DB_HOST:-$(hostname --long)} `

All the bebop login nodes are set as listen addresses in postgresql.conf:

```
listen_addresses = 'beboplogin1, beboplogin2, beboplogin3, beboplogin4, beboplogin5, beboplogin6'
```

Assuming that the database is started on a bebop login node,
and `env-bebop.sh` has been sourced then
the DB should be listening correctly when clients try to connect
to `DB_HOST`. 

In pg_hba.conf, this line allows any connection from within Bebop,
so that worker pools can communicate with the DB from the 
compute nodes.

```
# Allow connection from within bebop
host    all             all             .lcrc.anl.gov          trust
```

## Client Connection

Client must `export DB_HOST=<where db is running>` before attempting
to connect to the database. Otherwise, running db-settings.sh sets
DB_HOST to the client's host which may not be the DB host. So for example,
if the DB is running beboplogin1, and the client is on beboplogin2, then
this works.

```
[collier@beboplogin2 db]$ export DB_HOST=beboplogin1.lcrc.anl.gov
[collier@beboplogin2 db]$ ./db-print.sh
...
```

Without the correct DB_HOST export, you will see something like this:

```
[collier@beboplogin2 db]$ ./db-print.sh
psql: error: connection to server at "beboplogin2.lcrc.anl.gov" (10.70.128.6), port 11219 failed: Connection refused
```

where you can see the client trying to use its own hostname.

## SSH Tunnel

See db-tunnel.sh.  

An example tunnel string assuming the db is running on beboplogin1.

```
ssh -L :11219:beboplogin1.lcrc.anl.gov:11219 username@beboplogin1.lcrc.anl.gov -N -C
```

This sets localhost:11219 to point to beboplogin1 11219

And set the appropriate exports:

```
$ export DB_HOST=localhost
$ export DB_USER=collier
```

`db-print.sh` run on the client will then print the contents of the DB on beboplogin1.

Without `DB_USER=collier`, db-settings.sh will set DB_USER to the os user name. And you
may get an error if the os user name doesn't match that of the db.

```
$ ./db-print.sh 
psql: error: FATAL:  role "nick" does not exist
```

I think collier works here because the DB is running under my user name `collier`
on Bebop. I think if I'd added a 'nick' as a user to the DB it might work.
