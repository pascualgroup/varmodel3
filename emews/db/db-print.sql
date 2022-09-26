
/*
 DB PRINT SQL
 Just dump the tables for human inspection
 See db-print.sh for usage
*/

\dt

\echo == EMEWS EXPID TASKS ==
select * from eq_exp_id_tasks;
\echo == EMEWS TASKS ==
select * from eq_tasks;

\echo == EMEWS QUEUE IN ==
select * from emews_queue_IN;
\echo == EMEWS QUEUE OUT ==
select * from emews_queue_OUT;
