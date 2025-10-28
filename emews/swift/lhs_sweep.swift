import io;
import sys;
import files;
import string;
import python;
import R;
import location;
import unix;
import emews;

// deletes the specified directory
app (void o) rm_dir(string dirname) {
  "rm" "-rf" dirname;
}

// deletes the specified directories
app (void o) rm_dirs(file dirnames[]) {
  "rm" "-rf" dirnames;
}

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");
string site = getenv("SITE");

file upf = input(argv("f"));
string varmodel_x = argv("varmodel_x");
int replicates = string2int(argv("replicates"));
string biting_rate_multiplier_file = argv("biting_rate_multiplier_file");
string default_params_file = argv("default_params_file");
string measurement_file = argv("measurement_file");

string parse_params_template = """
import json
import numpy as np
import hashlib

turbine_output = '%s'
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

ect_bc = final_map.pop('ectopic_recombination_rate_BC')
final_map['ectopic_recombination_rate'] = [6.0e-5, ect_bc]

sr_a = final_map.pop('switching_rate_A')
sr_b = final_map.pop('switching_rate_BC')
final_map['switching_rate'] =[sr_a, sr_b]

var_groups_ratio_regional_pool_A = final_map.pop('var_groups_ratio_regional_pool_A')
final_map['var_groups_ratio_regional_pool'] = [var_groups_ratio_regional_pool_A, 1 - var_groups_ratio_regional_pool_A]

# model expects this to be an int
final_map['n_genes_initial'] = int(round(final_map['n_genes_initial']))
# "n_alleles_per_locus_initial": "n_genes_initial/10",
final_map['n_alleles_per_locus_initial'] = int(round(final_map['n_genes_initial'] / 20))

final_map["measurement_error_file_loc"] = f"{turbine_output}/misc/measurement_error_file_loc/"
final_map["MOI_estimation_info_file_loc"] = f"{turbine_output}/misc/MOI_estimation_info_file_loc/"

instance = final_map.pop('instance')
params = '{}!{}'.format(instance, json.dumps(final_map))
""";

string stage_params = """
import json
import os
import pandas as pd
import subprocess

replicate = %d
instance_root = '%s'
instance_dir = f'{instance_root}_{replicate}'
print("CREATING:", instance_dir, flush=True)
os.makedirs(instance_dir)

params = json.loads('%s')
params['rng_seed'] = params['rng_seed'] + replicate
with open(f'{instance_dir}/parameters.json', 'w') as f_out:
    json.dump(params, f_out)

os.environ['SITE'] = "%s"
emews_project_root = "%s"
os.environ['EMEWS_PROJECT_ROOT'] = emews_project_root

model_sh = f"{emews_project_root}/scripts/run_model.sh"
#  "bash" model_sh varmodel_x instance_dir @stdout=out @stderr=err;
varmodel_x = "%s"
cmd = ["bash", model_sh, varmodel_x, instance_dir]

try:
    proc = subprocess.run(cmd, cwd=instance_dir, capture_output=True, text=True)
    with open(f'{instance_dir}/out.txt', 'w') as fout:
        fout.write(proc.stdout)
    with open(f'{instance_dir}/err.txt', 'w') as fout:
        fout.write(proc.stderr)
except subprocess.CalledProcessError as e:
    with open(f'{instance_dir}/out.txt', 'w') as fout:
        fout.write(e.stdout)
    with open(f'{instance_dir}/err.txt', 'w') as fout:
        fout.write(e.stderr)
except OSError as e:
    print(f"OSError: {e}", flush=True)

sql_output = f"{instance_dir}/output.sqlite"
""";

(string result) run_task(string param_line) {
  string code = parse_params_template % (turbine_output, param_line, default_params_file, 
                                         biting_rate_multiplier_file);
  string params_str = python_persist(code, "params");
  string ps[] = split(params_str, "!");
  string instance = ps[0];
  string params = ps[1];

  // submission script should create this directory
  string instance_root = "%s/instances/instance_%s" % (turbine_output, instance);
  int i = 1;
  string stage_code = stage_params % (i, instance_root, params, site, emews_root, varmodel_x);
  result = python_persist(stage_code, "sql_output");
  // string db = "%s/output.sqlite" % instance_dir =>
  // string result_code = compute_result_code % (db, out_f, err_f);
  // result = python_persist(result_code, "result_j");
}


main() {
    string results[];
    string upf_lines[] = file_lines(upf);
    foreach line, i in upf_lines {
        results[i] = run_task(line);
    }
    int sz = size(results);
    printf("size: %d", sz);
}
