import python;
import string;
import sys;

string split_json_template = """
import json

vals = json.loads('%s')
str_vals = []
for val in vals:
    str_vals.append(json.dumps(val))
params_str = ';'.join(str_vals)
""";

(string vals[]) parse_json_list(string json_list) {
    // Python stack trace can be passed in as the parameter list
    // if ME itself raises an Error. This checks for that
    // and prints that error.
    if (find(json_list, "Traceback (", 0, -1) == 0) {
        printf(json_list);
    }
    string code = split_json_template % json_list;
    string params_str = python_persist(code, "params_str");
    vals = split(params_str, ";");
}
