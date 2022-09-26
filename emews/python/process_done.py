import db_tools
import logging
import os
import sweep_me
import csv
import json


def find_instance(line):
    end = line.rfind('_')
    tag = 'instance_'
    start = line.rfind(tag)
    return int(line[start + len(tag): end])


def run2():
    upf = '/project2/pascualmm/ncollier/repos/varmodel3/emews/experiments/lhs_8K.aa.2.0/upf.txt'
    params = {}
    with open(upf) as f_in:
        for line in f_in.readlines():
            line = line.strip()
            vals = json.loads(line)
            params[vals['instance']] = vals

    os.environ['DB_HOST'] = 'midway2-login2.rcc.local'
    os.environ['DB_PORT'] = '58577'
    db = db_tools.setup_db(log_level=logging.WARN, envs=True)
    db.connect()

    db.execute('select json_in from eq_tasks, emews_queue_in where eq_tasks.eq_task_id = emews_queue_in.eq_task_id;')
    rs = db.cursor.fetchall()
    db.close()
    for r in rs:
        instance = find_instance(r[0])
        params.pop(instance)
    
    update_upf = '/project2/pascualmm/ncollier/repos/varmodel3/emews/data/parameters/lhs_8K.aa.v2'
    with open(update_upf, 'w') as f_out:
        for p in params.values():
            f_out.write(f'{json.dumps(p)}\n')


def run():
    os.environ['DB_HOST'] = 'midway2-login2.rcc.local'
    os.environ['DB_PORT'] = '60926'
    db = db_tools.setup_db(log_level=logging.WARN, envs=True)
    db.connect()

    db.execute('select json_in from eq_tasks, emews_queue_in where eq_tasks.eq_task_id = emews_queue_in.eq_task_id;')
    rs = db.cursor.fetchall()
    print(rs)
    # db.close()
    # mf = '/project2/pascualmm/ncollier/repos/varmodel3/emews/experiments/lhs_8K.aa.2.0/measurement.txt'
    # results = []
    # for r in rs:
    #     input_file = r[0]
    #     vals = sweep_me.compute_targets(input_file, 30540, mf)
    #     results.append(vals)

    # with open('partial_8K_aa_results.csv', 'w', newline='') as fout:
    #     writer = csv.writer(fout)
    #     writer.writerow(['preval','MOIvar','pts','nbstrain','nbgene'])
    #     writer.writerows(results)

if __name__ == '__main__':
    run()