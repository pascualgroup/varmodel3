import subprocess
import yaml
import os
import time
import sys

from eqsql.task_queues import local_queue
import _targets


MOI_EST_F = '/project/jozik/ncollier/repos/varmodel3/emews/python/MOIestObjs.pkl'
AT = 35940


def parse_instance(line):
    end = line.rfind('_')
    tag = 'instance_'
    start = line.rfind(tag)
    return int(line[start + len(tag): end])


def check_worker_pool():
    job_id = os.environ['SLURM_JOB_ID']
    #  print('JOB_ID: ', job_id, flush=True)
    result = subprocess.run(['sacct', '-P', '-j', job_id], capture_output=True)
    lines = result.stdout.decode('utf-8').split('\n')
    worker_pool_status = ''
    for line in lines:
        items = line.split('|')
        if len(items) > 2 and items[1].startswith('tclsh'):
            # print("WORKER POOL: ", items, flush=True)
            worker_pool_status = items[5]

    return worker_pool_status == 'RUNNING'


def compute_targets(input_file, prop):
    instance_dir = os.path.dirname(input_file)
    instance = parse_instance(input_file)
    result = _targets.CalculTargetsMeasFunc(input_file, AT, prop, MOI_EST_F, 'mixtureDist')
    preval, MOIvar, pts, ptsA, ptsBC, nbstrain, nbgene, nbgeneA, nbgeneBC = result
    # preval, MOIvar, pts, ptsA, ptsBC, nbstrain, nbgene, nbgeneA, nbgeneBC
    result_file = os.path.join(instance_dir, f'result_{instance}.csv')
    with open(result_file, 'w') as f_out:
        f_out.write('preval,MOIvar,pts,ptsA,ptsBC,nbstrain,nbgene,nbgeneA,nbgeneBC,prop,instance\n')
        f_out.write(f'{preval},{MOIvar},{pts},{ptsA},{ptsBC},{nbstrain},{nbgene},{nbgeneA},{nbgeneBC},{prop},{instance}\n')
    return result


def run(exp_id, cfg_file):
    with open(cfg_file) as fin:
        params = yaml.safe_load(fin)

    task_queue = None
    try:
        task_queue = local_queue.init_task_queue(params['db_host'], params['db_user'],
                                                 port=params['db_port'], db_name=params['db_name'])

        upf = params['upf']
        with open(upf) as fin:
            payloads = [line.strip() for line in fin]

        job_id = int(os.environ['SLURM_JOB_ID'])
        # use the job_id as the work type to match with the
        # worker pool running as part of the same job
        # higher priority so worker pool grabs it first
        _, fts = task_queue.submit_tasks(exp_id, job_id, payloads)

        # check every 10 minutes up to an hour if worker pool is still up
        # then poll for results.
        worker_pool_up = True
        for _ in range(6):
            time.sleep(600)
            worker_pool_up = check_worker_pool()
            if not worker_pool_up:
                break

        for ft in task_queue.as_completed(fts, sleep=10):
            _, db_file = ft.result()
            compute_targets(db_file, 0.47)
            # os.unlink(db_file)

    finally:
        if task_queue is not None:
            task_queue.close()


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
