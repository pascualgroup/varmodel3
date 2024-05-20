import files;
import string;
import sys;
import io;
import python;
import R;
import location;
import unix;
import emews;

import EQSQL;

// deletes the specified directory
app (void o) rm_dir(string dirname) {
  "rm" "-rf" dirname;
}

// deletes the specified directories
app (void o) rm_dirs(file dirnames[]) {
  "rm" "-rf" dirnames;
}

int PROCESS_RESULT_TASK_TYPE = 100;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");
int resident_work_rank = string2int(getenv("RESIDENT_WORK_RANK"));

int TASK_TYPE = string2int(argv("task_type"));
printf("TASK TYPE: %s", TASK_TYPE);
int BATCH_SIZE = string2int(argv("batch_size"));
int BATCH_THRESHOLD = string2int(argv("batch_threshold", "1"));
string WORKER_POOL_ID = argv("worker_pool_id", "default");

file model_sh = input(emews_root+"/scripts/run_model.sh");
string varmodel_x = argv("varmodel_x");
int replicates = string2int(argv("replicates"));
string biting_rate_multiplier_file = argv("biting_rate_multiplier_file");
string default_params_file = argv("default_params_file");
string measurement_file = argv("measurement_file");

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

string compute_targets = """
import _targets
import json
import pandas as pd

payload = '%s'
df1_f, df2_f, df3_f, eq_id = json.loads(payload)
print(f'{eq_id} START', flush=True)
df1 = pd.read_parquet(df1_f)
df2 = pd.read_parquet(df2_f)
df3 = pd.read_parquet(df3_f)

moi_est_f = '/project/jozik/ncollier/repos/varmodel3/emews/python/MOIestObjs.pkl'
at = 35940
prop = 0.47
result = _targets.do_calc_targets(df1, df2, df3, at, prop, moi_est_f, 'mixtureDist')
jresult = json.dumps(result)
print(f'{eq_id} END {jresult}', flush=True)
""";

app (file out, file err) run_model(string instance_dir) {
    "bash" model_sh varmodel_x instance_dir @stdout=out @stderr=err;
}

app (void o) rm(string filename) {
    "rm" filename;
}

(string result) run_task(int task_id, string param_line) {
  if (TASK_TYPE == PROCESS_RESULT_TASK_TYPE) {
    string target_code = compute_targets % param_line;
    result = python_persist(target_code, "jresult");
  } else {
    string code = parse_params_template % (param_line, default_params_file, biting_rate_multiplier_file);
    string params_str = python_persist(code, "params");
    string ps[] = split(params_str, "!");
    string instance = ps[0];
    string params = ps[1];

    // submission script should create this directory
    string instance_root = "%s/instances/instance_%s" % (turbine_output, instance);
    int i = 1;
    stage_code = stage_params % (i, instance_root, params);
    instance_dir = python_persist(stage_code, "instance_dir");
    
    string out_f = "%s/out.txt" % instance_dir;
    string err_f = "%s/err.txt" % instance_dir;
    file out <out_f>;
    file err <err_f>;
    (out,err) = run_model(instance_dir) =>
    string db = "%s/output.sqlite" % instance_dir =>
    string result_code = compute_result_code % (db, out_f, err_f);
    result = python_persist(result_code, "result_j");
  }
}


run(message msgs[]) {
  // printf("MSGS SIZE: %d", size(msgs));
  foreach msg, i in msgs {
    result_payload = run_task(msg.eq_task_id, msg.payload);
    eq_task_report(msg.eq_task_id, TASK_TYPE, result_payload);
  }
}


(void v) loop(location querier_loc) {
  for (boolean b = true;
       b;
       b=c)
  {
    message msgs[] = eq_batch_task_query(querier_loc);
    boolean c;
    if (msgs[0].msg_type == "status") {
      if (msgs[0].payload == "EQ_STOP") {
        printf("loop.swift: STOP") =>
          v = propagate() =>
          c = false;
      } else {
        // sleep to give time for Python etc.
        // to flush messages
        sleep(5);
        printf("loop.swift: got %s: exiting!", msgs[0].payload) =>
        v = propagate() =>
        c = false;
      }
    } else {
      run(msgs);
      c = true;
    }
  }
}

(void o) start() {
  location querier_loc = locationFromRank(resident_work_rank);
  eq_init_batch_querier(querier_loc, WORKER_POOL_ID, BATCH_SIZE, BATCH_THRESHOLD, TASK_TYPE) =>
  printf("STARTING") =>
  loop(querier_loc) => {
    eq_stop_batch_querier(querier_loc);
    o = propagate();
  }
}

start() => printf("worker pool: normal exit.");