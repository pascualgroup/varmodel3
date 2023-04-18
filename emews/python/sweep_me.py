import argparse
import eq
import time
import csv
import Targets
from multiprocessing import Pool, Manager
import os
import subprocess
import datetime


def create_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--exp_id', required=True, help="experiment id")
    parser.add_argument('-u', '--upf', required=True, help="Upf file")
    parser.add_argument('-m', '--measurement', required=True, help="measurement file")
    parser.add_argument('-t', '--time', type=int, required=True, help='Time to take result at')

    return parser


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


def compute_targets(input_file, result_at, measurement_file, props):
    instance_dir = os.path.dirname(input_file)
    with Manager() as manager:
        lock = manager.Lock()
        inputs = [(input_file, result_at, measurement_file, prop, lock)
                  for prop in props]
        with Pool(3) as p:
            map_results = p.starmap(Targets.calc_targets_meas, inputs)
    # preval, MOIvar, pts, nbstrain, nbgene = Targets.calc_targets_meas(input_file, result_at, measurement_file, prop)
    #  = calc_targets.calc_targets(input_file, result_at, measurement_file)
    # print(map_results)
    results = []
    for i, result in enumerate(map_results):
        prop = props[i]
        preval, MOIvar, pts, nbstrain, nbgene = result
        result_file = os.path.join(instance_dir, f'results_{prop}_{result_at}.csv')
        print(f'Writing result file: {result_file}', flush=True)
        instance = parse_instance(instance_dir)
        with open(result_file, 'w') as f_out:
            f_out.write('preval,MOIvar,pts,nbstrain,nbgene,instance\n')
            f_out.write(f'{preval},{MOIvar},{pts},{nbstrain},{nbgene},{instance}\n')
            results.append([preval, MOIvar, pts, nbstrain, nbgene, prop, instance])
    return results


def run(exp_id, upf, measurement, result_at):
    print("INITIALIZING DB", flush=True)
    eq.init(retry_threshold=200)
    print("DB INITIALIZED", flush=True)

    # get number of lines
    with open(upf) as f_in:
        num_runs = sum(1 for _ in f_in)
    print(f'NUM RUNS: {num_runs}', flush=True)

    job_id = int(os.environ['SLURM_JOB_ID'])
    # use the job_id as the work type to match with the
    # worker pool running as part of the same job
    # higher priority so worker pool grabs it first
    status, task_id = eq.submit_task(exp_id, job_id, f'{num_runs}', priority=10)
    # Get the 'OK' back
    status, result_str = eq.query_result(task_id, timeout=120.0)
    print(f'Init num runs: {status}, {result_str}', flush=True)

    task_ids = []
    with open(upf) as f_in:
        for line in f_in.readlines():
            line = line.strip()
            status, task_id = eq.submit_task(exp_id, job_id, line)
            if status != eq.ResultStatus.SUCCESS:
                print(f'Error submitting task {status}', flush=True)
            task_ids.append(task_id)

    # check every 10 minutes up to an hour if worker pool is still up
    # then poll for results.
    worker_pool_up = True
    for _ in range(6):
        time.sleep(600)
        worker_pool_up = check_worker_pool()
        if not worker_pool_up:
            break

    results = []
    go = True
    finished_tasks = set()
    while go and worker_pool_up:
        for task_id in task_ids:
            if task_id not in finished_tasks:
                status, result_str = eq.query_result(task_id, timeout=4.0)
                if status == eq.ResultStatus.SUCCESS:
                    print(f'{task_id}: {result_str}', flush=True)
                    start_t = datetime.datetime.now()
                    # props: 0.45, 0.50, 0.57, 0.60, and 0.65
                    # [0.45, 0.50, 0.57, 0.60, 0.65]
                    result = compute_targets(result_str, result_at, measurement, [
                                             0.45, 0.57, 0.65])
                    print(f'Deleting {result_str}', flush=True)
                    os.remove(result_str)
                    end_t = datetime.datetime.now()
                    print(f'{task_id} runtime: {(end_t - start_t).total_seconds()}')
                    for r in result:
                        results.append(r)
                    finished_tasks.add(task_id)
                    go = len(finished_tasks) != num_runs
                elif result_str == eq.EQ_ABORT:
                    print(f'Aborting while querying task {task_id}', flush=True)
                    go = False
                    break

        worker_pool_up = check_worker_pool()
        if not worker_pool_up:
            print('Aborting: Worker Pool Down')
        time.sleep(60)

    eq.stop_worker_pool(job_id)
    eq.close()

    with open('full_results.csv', 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['preval', 'MOIvar', 'pts', 'nbstrain', 'nbgene, prop, instance'])
        writer.writerows(results)
    print('SWEEP ME DONE', flush=True)


if __name__ == '__main__':
    print("ME STARTING", flush=True)
    parser = create_arg_parser()
    args = parser.parse_args()
    run(args.exp_id, args.upf, args.measurement, args.time)
