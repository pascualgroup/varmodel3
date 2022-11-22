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

string compute_result_code = """
import json
import calc_targets
import os

input_file = '%s'
instance_dir = '%s'
time = %s
measurement_file = '%s'
preval, MOIvar, pts, nbstrain, nbgene = calc_targets.calc_targets(input_file, time, measurement_file);
with open(os.path.join(instance_dir, f'results_{time}.csv'), 'w') as f_out:
    f_out.write('preval,MOIvar,pts,nbstrain,nbgene\\\\n')
    f_out.write(f'{preval},{MOIvar},{pts},{nbstrain},{nbgene}\\\\n')
os.remove(input_file)
result_json = json.dumps({'preval': preval, 'MOIvar': MOIvar, 'pts': pts, 'nbstrain': nbstrain, 'nbgene': nbgene})
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
    db = "%s/output.sqlite" % instance_dir;
    result_code = compute_result_code % (db, instance_dir, result_at, measurement_file);
    python_persist(result_code, "result_json") =>
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
