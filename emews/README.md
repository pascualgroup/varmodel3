# Varmodel EMEWS Workflows #

This directory contains the code for running an EMEWS workflow
on the varmodel.

## Sweep ##

Given an unrolled parameter file (UPF), where each line
in the file is a set of parameters to run, the sweep workflow
will pass each line to a varmodel instance to run. 

To run the sweep, the workflow enviroment must be loaded.

```
source /scratch/cgsb/pascual/envs/swift-t-julia.sh 
```

This will perform some module loads and add some relevant
paths to the PATH and LD_LIBRARY_PATH variables. This only
needs to be sourced once per terminal instance. So for example,
if you login, and then do the above "source", you can submit
any number of workflows without sourcing again. You only need
to source the environment once you've logged out and in again.

Once the environment has been sourced, the workflow is submitted
using the `swift/greene_run_sweep.sh` script. This takes two arguments,

1. An experiment id
2. A configuration file

For example,

```
cd swift
./greene_run_sweep.sh lhs_10K_1.0 ../data/cfgs/greene_gi_off_sweep.cfg
```

This will submit the job to Greene's slurm scheduler, and create an
`experiments/<experiment id> directory in which the model runs will run.
In the above example, an `experiments/lhs_10K_1.0` directory is created.
File relevant to the workflow run will be copied in there (e.g., the default
parameters, the configuration file, etc.)

Once the job starts, the experiment directory will contain `instance` subdirectories
in which each model run occurs. For example,

```
gi_on_032025_lhs_sweep_503_3.0
├── biting_rate_multipliers.txt
├── cfg.cfg
├── default_parameters.json
├── instances
│   ├── instance_1063_1
│   │   ├── err.txt
│   │   ├── output.sqlite
│   │   ├── out.txt
│   │   └── parameters.json
│   ├── instance_1068_1
│   │   ├── err.txt
│   │   ├── output.sqlite
│   │   ├── out.txt
│   │   └── parameters.json
...
```

### Configuring a Run ###

The workflow runs are configured using a configuration file that specifies
the UPF to use, the walltime, number of nodes, and other model specific parameters.

The two files that can be used are `data/cfgs/greene_gi_off_sweep.cfg` and
`data/cfgs/greene_gi_on_sweep.cfg`. Their contents are the same except for
the default parameters specification: gi_off uses the "off" parameters and
gi_on uses the "on" parameters, defined in the CFG_DEFAULT_PARAMS configuration
entry.

The configuration file has the following entries:

```
# How long to run
CFG_WALLTIME=12:00:00
# Assigns the max memory per node and insures that only
# our job can use it.
CFG_SBATCH_ARGS="--mem=180G\n#SBATCH --exclusive"
# Processes per node: essentially the number of varmodel runs we can do on each node.
CFG_PPN=24
# The number of nodes to use
NODES=6

# The task type -- corresponds to what bucket the runs defined in the
# UPF will be.
CFG_TASK_TYPE=500
# The name of the UPF -- assumed to be the emews/data/upfs directory
CFG_UPF="${CFG_TASK_TYPE}_05072025_1.0_upf.txt"

# The total number of processes assigned to this job.
# Essentially the total number of varmodel runs that can
# be run at one time. Note that 1 is used for the workflow overhead.
CFG_PROCS=$(( NODES * CFG_PPN ))

# Varmodel repo location -- edit this to match yours
VARMODEL_ROOT=/scratch/cgsb/pascual/ncollier/repos/varmodel3

# Locations of files used in the workflow. These should not need to be changed
CFG_DEFAULT_PARAMS=$VARMODEL_ROOT/emews/data/parameters/default_params_GI_on_03262025.json
CFG_BITING_RATE_MULTIPLIERS=$VARMODEL_ROOT/emews/data/parameters/mosquito_population.txt
CFG_MEASUREMENT_FILE=$VARMODEL_ROOT/emews/data/parameters/measurement.txt
CFG_RESULT_AT=35940
CFG_VARMODEL_X=$VARMODEL_ROOT/run.jl
CFG_MEAS_ERROR=$VARMODEL_ROOT/emews/data/parameters/measurement_error_file_loc
CFG_MOI_INFO=$VARMODEL_ROOT/emews/data/parameters/MOI_estimation_info_file_loc
```

Of the configuration entries, the important ones are `CFG_PPN` and `CFG_TASK_TYPE`.
`CFG_PPN` determines how many varmodel runs will run on a single node. A lower PPN results
in more memory being avaiable to each varmodel run, and a higher PPN results in
less. Consequently, the long running varmodel runs need a lower PPN as they will
consume more memory. `CFG_TASK_TYPE` is a tag
identifying the predicted runtime bucket for a set of runs. 

| Task Type | Runtime | GI  | Midway PPN |
| --------- | ------- | --- | ---------- |
| 400       | Low     | off | 24         |
| 401       | Low Mid | off | 24         |
| 402       | High Mid| off | 15         |
| 403       | High    | off | 10         |
| 500       | Low     | on  | 24         |
| 501       | Low Mid | on  | 24         |
| 502       | High Mid| on  | 15         |
| 503       | High    | on  | 10         |

Given the task type tag, you can set the CFG_PPN appropriately.

## Sofware ##

The required software can be found in `/scratch/cgsb/pascual/sfw`.

* `julia` - the julia installation
* `gcc-11.3.1/tcl-8.6.17` - tcl installed used by the emews / swift-t workflow
* `openmpi-4.1.6/swift-t-10272025` - the swift-t emews install

### Julia ###

A shared julia has been installed in `/scratch/cgsb/pascual/sfw/julia`.
The easiest way to use it is to:

```
source emews/scripts/greene_env.sh
```

which will set the paths (JULIA_DEPOT_PATH and PATH) correctly. Once
source, `julia` will start the julia command line interpreter. Note
that on Greene doing anything remotely computational, like running a
julia script, on a login node seems to get you kicked off.
