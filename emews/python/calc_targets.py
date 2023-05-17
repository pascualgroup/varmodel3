from multiprocessing import Manager
import os
import glob
import sys
sys.path.append('/project2/pascualmm/ncollier/repos/varmodel3/post_analysis')

import _targets as targets


def parse_instance(line):
    end = line.rfind('_')
    tag = 'instance_'
    start = line.rfind(tag)
    return int(line[start + len(tag): end])


def run(exp_dir, result_at):
    instances_dir = os.path.join(exp_dir, 'instances')
    measurement_file = os.path.join(exp_dir, 'measurement.txt')
    sqlites = glob.glob(f'{instances_dir}/instance*/output.sqlite')
    n = len(sqlites)
    for i, sqldb in enumerate(sqlites):
        instance_dir = os.path.dirname(sqldb)
        with Manager() as manager:
            lock = manager.Lock()
            for p in [0.47, 0.57, 0.67]:
                print(f"{i + 1} of {n}: {sqldb}, {p}", flush=True)
                preval, MOIvar, pts, nbstrain, nbgene = targets.calc_targets_meas(
                    sqldb, result_at, measurement_file, p, lock)
                result_file = os.path.join(
                    instance_dir, f'results_{p}_{result_at}.csv')
                print(f'Writing result file: {result_file}', flush=True)
                instance = parse_instance(instance_dir)
                with open(result_file, 'w') as f_out:
                    f_out.write('preval,MOIvar,pts,nbstrain,nbgene,instance\n')
                    f_out.write(
                        f'{preval},{MOIvar},{pts},{nbstrain},{nbgene},{instance}\n')

        print(f'Deleting {sqldb}', flush=True)
        os.remove(sqldb)


if __name__ == '__main__':
    exp_dir = sys.argv[1]
    result_at = int(sys.argv[2])
    run(exp_dir, result_at)
