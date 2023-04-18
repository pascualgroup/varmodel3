
/**
   EMEWS loop.swift
*/

import assert;
import io;
import python;
import string;
import sys;
import unix;
import EQ;
import emews;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

int job_id = string2int(getenv("SLURM_JOB_ID"));
printf("JOB ID: %d", job_id);

file upf = input(argv("f"));
file model_sh = input(emews_root+"/scripts/run_model.sh");
string varmodel_x = argv("varmodel_x");
int replicates = string2int(argv("replicates"));
string biting_rate_multiplier_file = argv("biting_rate_multiplier_file");
string default_params_file = argv("default_params_file");

string parse_params_template = """
import json
import numpy as np
import hashlib

p_map = json.loads('%s')
# print(param_list, flush=True)

dp_file = '%s'
with open(dp_file) as fin:
    dp_map = json.load(fin)

bite_rate_mults = np.loadtxt('%s')


imabc_seed = p_map.pop('seed')
h = hashlib.md5(str.encode(imabc_seed)).digest()[:6]
rng_seed = int.from_bytes(h, byteorder='little')
p_map['rng_seed'] = rng_seed

# update / add default params with p_map params
final_map = dp_map.copy()
final_map.update(p_map)
bite_rate = bite_rate_mults * p_map['biting_rate']
final_map['biting_rate'] = bite_rate.tolist()
# model expects this to be an int
final_map['n_genes_initial'] = int(round(final_map['n_genes_initial']))
# "n_alleles_per_locus_initial": "n_genes_initial/10",
final_map['n_alleles_per_locus_initial'] = int(round(final_map['n_genes_initial'] / 10))

instance = final_map.pop('instance')
params = '{}!{}'.format(instance, json.dumps(final_map))
""";

string stage_params = """
import json
import os

replicate = %d
instance_root = '%s'
instance_dir = f'{instance_root}_{replicate}'
os.makedirs(instance_dir)

params = json.loads('%s')
params['rng_seed'] = params['rng_seed'] + replicate
with open(f'{instance_dir}/parameters.json', 'w') as f_out:
    json.dump(params, f_out)
""";

app (file out, file err) run(string instance_dir) {
    "bash" model_sh varmodel_x instance_dir @stdout=out @stderr=err;
}

app (void o) rm(string filename) {
    "rm" filename;
}

(string data_path) run_model(string param_line) {
    string code = parse_params_template % (param_line, default_params_file, biting_rate_multiplier_file);
    string params_str = python_persist(code, "params");
    string ps[] = split(params_str, "!");
    string instance = ps[0];
    string params = ps[1];
    // submission script should create this directory
    string instance_root = "%s/instances/instance_%s" % (turbine_output, instance);
    //foreach i in [0:replicates-1:1] {
    int i = 1;
    stage_code = stage_params % (i, instance_root, params);
    instance_dir = python_persist(stage_code, "instance_dir");
    
    string out_f = "%s/out.txt" % instance_dir;
    string err_f = "%s/err.txt" % instance_dir;
    file out <out_f>;
    file err <err_f>;
    (out,err) = run(instance_dir) =>
    data_path = "%s/output.sqlite" % instance_dir;
}

(int num_runs) get_num_runs() {
  message msg = eq_task_querier(job_id);
  eq_task_reporter(msg.eq_task_id, job_id, "OK") =>
  num_runs = string2int(msg.payload);
  printf("NUM RUNS: %d", num_runs);
}

(void v) loop() {
  int num_runs = get_num_runs();
  string results[];
  foreach i in [0:num_runs-1:1] {
    message msg = eq_task_querier(job_id);
    int eq_task_id = msg.eq_task_id;
    // printf("loop.swift payload %s: ", msg.payload);
    string output_file = run_model(msg.payload);
    results[i] = output_file;
    eq_task_reporter(eq_task_id, job_id, output_file);
    // eq_task_reporter(eq_task_id, 1, "OK");
  }
  // query the stop message
  printf("num runs: %d, result size: %d", num_runs, size(results)) =>
  eq_task_querier(job_id) =>
  v = propagate();
}

loop() => printf("loop.swift: normal exit.");
