
/**
    WORKFLOW SQL
    Initialize the PSQL DB for workflows
    See db-create.py for usage
*/

/* For SQLite: */
/* PRAGMA foreign_keys = ON; */

/* Each group here is a metadata collection of other groups or points
*/
create table eq_exp_id_tasks (
       /* e.g. 'experiment':'X-032' or 'iteration':421 */
       exp_id text,
       eq_task_id integer
);


/* Each row here is a model run
*/
create table eq_tasks (
       /* the task id; eq_id==0 is the dummy null point */
       eq_task_id integer PRIMARY KEY,
       /* See db_covid.py for valid status codes */
       eq_status integer,
       /* the task type, eq_type==0 means "any type" */
       eq_task_type integer,
       /* JSON-formatted payload from the OUT queue or status message, e.g., "EQ_STOP" */
       json_out text,
       /* JSON-formatted payload for the IN queue */
       json_in  text,
       /* time this task was created (json_out) */
       time_created timestamp,
       /* time this task started (json_out) */
       time_start timestamp,
       /* time this task finished (json_in) */
       time_stop  timestamp
);

/* This generator is just for the queues */
create sequence emews_id_generator start 1 no cycle;

create table emews_queue_OUT(
       /* the task type */
       eq_task_type integer,
       /* eq_id */
       eq_task_id integer,
       eq_priority integer
);

create table emews_queue_IN(
       /* the task type */
       eq_task_type integer,
       /*  eq_id */
       eq_task_id integer
);
