
/*
   EMEWS EQ.swift for SQL
*/

type message {
    int eq_task_id;
    string msg_type;
    string payload;
}

(string result) EQ_validate() {
    result = python_persist("import eq ; eq.init()", "eq.validate()");
}

string code_get = """
import sys
import eq
import json

eq.init(retry_threshold=200)
try:
    print('swift out_get', flush=True)
    # result is a msg map
    msg_map = eq.query_task(%i, timeout=7200.0)
    items = [msg_map['type'], msg_map['payload']]
    if msg_map['type'] == 'work':
        items.append(str(msg_map['eq_task_id']))
    result_str = '|'.join(items)
    # result_str = json.dumps(result)
    # result_str = '{}|{}'.format(eq_task_id, payload)

    print('swift out_get done', flush=True)
except Exception as e:
    import sys, traceback
    info = sys.exc_info()
    s = traceback.format_tb(info[2])
    print(str(e) + ' ... \\n' + ''.join(s))
    sys.stdout.flush()
    result_str = eq.ABORT_JSON_MSG
finally:
    eq.close()
""";

(message msg) eq_task_querier(int eq_type) {
    string msg_string = python_persist(code_get % eq_type, "result_str");
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
import sys
import eq
eq.init(retry_threshold=200)
try:
    eq_task_id = %i
    eq_type = %i
    # TODO this returns a ResultStatus, add FAILURE handling
    eq.report_task(eq_task_id, eq_type, r'%s')
finally:
    eq.close()
""";

(void v) eq_task_reporter(int eq_task_id, int eq_type, string result_payload) {
    // trace("code: " + code_put % (eq_type, eq_ids));
    python_persist(code_put % (eq_task_id, eq_type, result_payload)) =>
        v = propagate();

}

// Local Variables:
// c-basic-offset: 4
// End:
