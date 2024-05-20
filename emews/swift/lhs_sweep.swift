import io;
import sys;
import files;
import string;
import python;
import launch;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

file upf = input(argv("f"));
file model_sh = input(emews_root+"/scripts/run_model.sh");
string varmodel_x = argv("varmodel_x");
int replicates = string2int(argv("replicates"));
string biting_rate_multiplier_file = argv("biting_rate_multiplier_file");
string default_params_file = argv("default_params_file");
string measurement_file = argv("measurement_file");
string result_at = argv("result_at");

string parse_params_template = """
import json
import numpy as np
import hashlib

params = '%s'
dp_file = '%s'
mos_pop = '%s'

p_map = json.loads(params)
# print(param_list, flush=True)

with open(dp_file) as fin:
    dp_map = json.load(fin)

bite_rate_mults = np.loadtxt(mos_pop)

imabc_seed = p_map.pop('seed')
h = hashlib.md5(str.encode(imabc_seed)).digest()[:6]
rng_seed = int.from_bytes(h, byteorder='little')
p_map['rng_seed'] = rng_seed

# update / add default params with p_map params
final_map = dp_map.copy()
final_map.update(p_map)
bite_rate = bite_rate_mults * p_map['biting_rate']
final_map['biting_rate'] = bite_rate.tolist()


err_a = final_map.pop('ectopic_recombination_rate_A')
err_b = final_map.pop('ectopic_recombination_rate_BC')
final_map['ectopic_recombination_rate'] = [err_a, err_b]

sr_a = final_map.pop('switching_rate_A')
sr_b = final_map.pop('switching_rate_BC')
final_map['switching_rate'] =[sr_a, sr_b]

fbc = final_map.pop('functionality_BC')
final_map['var_groups_functionality'] = [1.0, fbc]

# model expects this to be an int
final_map['n_genes_initial'] = int(round(final_map['n_genes_initial']))
# "n_alleles_per_locus_initial": "n_genes_initial/10",
final_map['n_alleles_per_locus_initial'] = int(round(final_map['n_genes_initial'] / 20))

instance = final_map.pop('instance')
params = '{}!{}'.format(instance, json.dumps(final_map))
""";

string stage_params = """
import json
import os
import pandas as pd
import pyarrow
import fastparquet

replicate = %d
instance_root = '%s'
instance_dir = f'{instance_root}_{replicate}'
os.makedirs(instance_dir)

params = json.loads('%s')
params['rng_seed'] = params['rng_seed'] + replicate
with open(f'{instance_dir}/parameters.json', 'w') as f_out:
    json.dump(params, f_out)
""";

string compute_result_code = """
import json
import reduce_db
import os

input_file = '%s'
out_f = '%s'
err_f = '%s'
result = reduce_db.run(input_file)
os.unlink(input_file)
os.unlink(out_f)
os.unlink(err_f)
result_j = json.dumps(result)
""";

app (file out, file err) run(string instance_dir) {
    "bash" model_sh varmodel_x instance_dir @stdout=out @stderr=err;
}

app (void o) rm(string filename) {
    "rm" filename;
}

(void v) run_model(string param_line) {
    string code = parse_params_template % (param_line, default_params_file, biting_rate_multiplier_file);
    string params_str = python_persist(code, "params");
    string ps[] = split(params_str, "!");
    string instance = ps[0];
    string params = ps[1];
    // submission script should create this directory
    string instance_root = "%s/instances/instance_%s" % (turbine_output, instance);
    // string dbs[];
    string db;
    //foreach i in [0:replicates-1:1] {
    int i = 1;
    stage_code = stage_params % (i, instance_root, params);
    instance_dir = python_persist(stage_code, "instance_dir");
    
    string out_f = "%s/out.txt" % instance_dir;
    string err_f = "%s/err.txt" % instance_dir;
    file out <out_f>;
    file err <err_f>;
    (out,err) = run(instance_dir) =>
    string db = "%s/output.sqlite" % instance_dir =>
    stringresult_code = compute_result_code % (db, out_f, err_f);
    python_persist(result_code, "result_j") =>
    v = propagate();
}

// call this to create any required directories
app (void o) make_dir(string dirname) {
    "mkdir" "-p" dirname;
}

main() {
    string upf_lines[] = file_lines(upf);
    foreach line, i in upf_lines {
        run_model(line);
    }
}
