import io;
import sys;
import files;
import location;
import string;
import EQR;
import R;
import assert;
import python;
import unix;
import stats;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

string resident_work_ranks = getenv("RESIDENT_WORK_RANKS");
string r_ranks[] = split(resident_work_ranks,",");

string algo_file = argv("algo_file");
string algo_params_file = argv("algo_params_file");

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

param_list = json.loads('%s')
# print(param_list, flush=True)

dp_file = '%s'
with open(dp_file) as fin:
    dp_map = json.load(fin)

bite_rate_mults = np.loadtxt('%s')

all_params = ''
for p_map in param_list:
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
    if len(all_params) > 0:
      all_params = '{};{}'.format(all_params, params)
    else:
      all_params = params
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
time = %s
measurement_file = '%s'
preval, MOIvar, pts, nbstrain, nbgene = calc_targets.calc_targets(input_file, time, measurement_file);
result_json = json.dumps({'MOI': MOIvar, 'prevalence': preval, 'numberuniquegenes': nbgene, 
               'numberuniquestrains': nbstrain, 'PTS': pts})
os.remove(input_file)
""";

app (file out, file err) run(string instance_dir) {
    "bash" model_sh varmodel_x instance_dir @stdout=out @stderr=err;
}

app (void o) rm(string filename) {
    "rm" filename;
}

(string result) obj(string params_str) {
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
    // }
    result_code = compute_result_code % (db, result_at, measurement_file);
    result = python_persist(result_code, "result_json");
}

(void v) loop(location ME) {

    for (boolean b = true, int i = 1;
       b;
       b=c, i = i + 1)
    {
        string payload_str =  EQR_get(ME);
        // string payload[] = split(payload_str, "|");
        // string payload_type = payload[0];
        boolean c;

        if (payload_str == "DONE") {
            string finals =  EQR_get(ME);
            printf("Results: %s", finals) =>
            v = make_void() =>
            c = false;

        } else if (payload_str == "EQR_ABORT") {
            printf("EQR aborted: see output for R error") =>
            string why = EQR_get(ME);
            printf("%s", why) =>
            v = propagate() =>
            c = false;
        } else {
            string results[];
            string param_code = parse_params_template % (payload_str, default_params_file, 
                                                         biting_rate_multiplier_file);
            string json_params[] = split(python_persist(param_code, "all_params"), ";");
            foreach params, j in json_params {
                results[j] = obj(params);
            }
            string res = join(results, ";");
            EQR_put(ME, res) => c = true;
        }
    }
}


(void o) start(int ME_rank) {
    location ME = locationFromRank(ME_rank);
    //string algorithm = strcat(emews_root, "/R/easyabc.R");
    EQR_init_script(ME, algo_file) =>
    EQR_get(ME) =>
    EQR_put(ME, algo_params_file) =>
    loop(ME) => {
        EQR_stop(ME) =>
        EQR_delete_R(ME);
        o = propagate();
    }
}

// deletes the specified directory
app (void o) rm_dir(string dirname) {
    "rm" "-rf" dirname;
}

main() {

    assert(strlen(emews_root) > 0, "Set EMEWS_PROJECT_ROOT!");

    int ME_ranks[];
    foreach r_rank, i in r_ranks {
      ME_ranks[i] = toint(r_rank);
    }

    foreach ME_rank, i in ME_ranks {
    start(ME_rank) =>
        printf("End rank: %d", ME_rank);
    }
}
