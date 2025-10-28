import location;
pragma worktypedef resident_work;

type message {
    int eq_task_id;
    string msg_type;
    string payload;
}

string code_get = """
import os
import eqsql_swift

eq_work_type = %i
worker_pool_id = '%s'
try:
    retry_threshold = int(os.environ.get('EQ_DB_RETRY_THRESHOLD', 10))
except ValueError as e:
    print("ENV VAR: EQ_DB_RETRY_THRESHOLD must be an integer")
    raise e

try:
    query_timeout = float(os.environ.get('EQ_QUERY_TASK_TIMEOUT', 120.0))
except ValueError as e:
    print("ENV VAR: EQ_QUERY_TASK_TIMEOUT must be a float")
    raise e

result_str = eqsql_swift.query_task(eq_work_type, worker_pool=worker_pool_id, query_timeout=query_timeout, retry_threshold=retry_threshold)
""";

(message msg) eq_task_query(int eq_type, string worker_pool_id) {
    string msg_string = python_persist(code_get % (eq_type, worker_pool_id), "result_str");
    // string msg_string = python_persist(code_parse_msg % result, "result_str");
    string msg_parts[] = split(msg_string, "|");
    msg.msg_type = msg_parts[0];
    msg.payload = msg_parts[1];
    if (size(msg_parts) == 3) {
        msg.eq_task_id = string2int(msg_parts[2]);
    } else {
        msg.eq_task_id = -1;
    }
}

string code_put = """
import os
import eqsql_swift

try:
    retry_threshold = int(os.environ.get('EQ_DB_RETRY_THRESHOLD', 10))
except ValueError as e:
    print("ENV VAR: EQ_DB_RETRY_THRESHOLD must be an integer")
    raise e

eq_task_id = %i
eq_type = %i
payload = r'%s'

eqsql_swift.report_task(eq_task_id, eq_type, payload, retry_threshold=retry_threshold)
""";

(void v) eq_task_report(int eq_task_id, int eq_type, string result_payload) {
    // trace("code: " + code_put % (eq_type, eq_ids));
    python_persist(code_put % (eq_task_id, eq_type, result_payload)) =>
        v = propagate();
}


@dispatch=resident_work
(void v) _void_py(string code, string expr="\"\"") "turbine" "0.1.0"
    [ "turbine::python 1 1 <<code>> <<expr>> "];


// TODO: Update repo version with timeout
@dispatch=resident_work
(string output) _string_py(string code, string expr) "turbine" "0.1.0"
    [ "set <<output>> [ turbine::python 1 1 <<code>> <<expr>> ]" ];

string init_querier_string = """
import eqsql_swift
import os

try:
    retry_threshold = int(os.environ.get('EQ_DB_RETRY_THRESHOLD', 10))
except ValueError as e:
    print("ENV VAR: EQ_DB_RETRY_THRESHOLD must be an integer")
    raise e

try:
    query_timeout = float(os.environ.get('EQ_QUERY_TASK_TIMEOUT', 120.0))
except ValueError as e:
    print("ENV VAR: EQ_QUERY_TASK_TIMEOUT must be a float")
    raise e

eqsql_swift.init_task_querier('%s', %d, %d, %d, query_timeout, retry_threshold)
""";

(void v) eq_init_batch_querier(location loc, string worker_pool, int batch_size, int threshold, int work_type) {
    string code = init_querier_string % (worker_pool, batch_size, threshold, work_type);
    // printf("Code is: \n%s", code);
    @location=loc _void_py(code) => v = propagate();
}

(void v) eq_stop_batch_querier(location loc){
    stop_string = "eqsql_swift.stop_task_querier()";
    @location=loc _void_py(stop_string) => v = propagate();
}


string get_string = "result = eqsql_swift.get_tasks_n()";

(message msgs[]) eq_batch_task_query(location loc) {
    //printf("eq_batch_task_query called");
    //printf("Code is: \n%s", code);
    result = @location=loc _string_py(get_string, "result");
    string msg_strs[] = split(result, ";");

    foreach msg_str, i in msg_strs {
        string msg_parts[] = split(msg_str, "|");
        message msg;
        msg.msg_type = msg_parts[0];
        msg.payload = msg_parts[1];
        if (size(msg_parts) == 3) {
            msg.eq_task_id = string2int(msg_parts[2]);
        } else {
            msg.eq_task_id = -1;
        }
        msgs[i] = msg;
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
