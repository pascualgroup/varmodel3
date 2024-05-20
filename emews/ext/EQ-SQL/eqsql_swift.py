""" Utility functions for interfacing swift code with database queues"""
import logging
import traceback
import threading
import time
import multiprocessing as mp
import os

from eqsql.task_queues import local_queue
from eqsql.task_queues.core import ABORT_MSG

password = None
check_password = True

def _create_eqsql(retry_threshold: int = 0, log_level=logging.WARN): 
    host = os.getenv('DB_HOST')
    user = os.getenv('DB_USER')
    if os.getenv('DB_PORT') is None or os.getenv('DB_PORT') == '':
        port = None
    else:
        port = int(os.getenv('DB_PORT'))

    global check_password, password
    if check_password:
        check_password = False
        if os.getenv('DB_PASSWORD_F') is not None and os.getenv('DB_PASSWORD_F') != '':
            password_f= os.getenv('DB_PASSWORD_F')
            with open(password_f) as fin:
                password = fin.readline().strip()

    db_name = os.getenv('DB_NAME')
    return local_queue.init_task_queue(host, user, port, db_name, password, retry_threshold, log_level)


def query_task(eq_work_type: int, worker_pool: str, query_timeout: float = 120.0,
               retry_threshold: int = 0, log_level=logging.WARN):
    eq_sql: local_queue.LocalTaskQueue = None
    try:
        eq_sql = _create_eqsql(retry_threshold, log_level)
        eq_sql.logger.debug('swift out_get')
        # result is a msg map
        msg_map = eq_sql.query_task(eq_work_type, worker_pool=worker_pool, timeout=query_timeout)
        items = [msg_map['type'], msg_map['payload']]
        if msg_map['type'] == 'work':
            items.append(str(msg_map['eq_task_id']))
        # result_str should be returned via swift's python persist
        eq_sql.logger.debug('swift out_get done')
        return '|'.join(items)
    except Exception:
        if eq_sql is None:
            print(f'eq_swift.query_task error {traceback.format_exc()}', flush=True)
        else:
            eq_sql.logger.error(f'eq_swift.query_task error {traceback.format_exc()}')
        # result_str returned via swift's python persist
        return ABORT_MSG
    finally:
        if eq_sql is not None:
            eq_sql.close()


def report_task(eq_task_id: int, eq_work_type: int, result_payload: str,
                retry_threshold: int = 0, log_level=logging.WARN):
    eq_sql = None
    try:
        eq_sql = _create_eqsql(retry_threshold, log_level)
        # TODO this returns a ResultStatus, add FAILURE handling
        eq_sql.report_task(eq_task_id, eq_work_type, result_payload)
    except Exception:
        if eq_sql is None:
            print(f'eq_swift.report_task error {traceback.format_exc()}', flush=True)
        else:
            eq_sql.logger.error(f'eq_swift.report_task error {traceback.format_exc()}')
    finally:
        if eq_sql is not None:
            eq_sql.close()


_q = mp.Queue(1)
_go = True

def query_tasks_n(batch_size: int, threshold: int, work_type: int, worker_pool: str,
                  timeout: float, retry_threshold: int, q: mp.Queue):
    running_task_ids = []
    wait = 0.25
    while _go:
        eq_sql = None
        try:
            eq_sql = _create_eqsql(retry_threshold)
            running_task_ids, tasks = eq_sql.query_more_tasks(work_type, running_task_ids,
                                                              batch_size=batch_size, threshold=threshold,
                                                              worker_pool=worker_pool, timeout=timeout)
        except Exception:
            if eq_sql is None:
                print(f'eq_swift.query_task_n error {traceback.format_exc()}', flush=True)
            else:
                eq_sql.logger.error(f'eq_swift.query_task_n error {traceback.format_exc()}')
            running_task_ids = []
            tasks = [ABORT_MSG]
        finally:
            if eq_sql is not None:
                eq_sql.close()

        n_tasks = len(tasks)
        # print("TASKS: ", tasks, flush=True)
        if n_tasks > 0:
            wait = 0.25
            if tasks[-1]['type'] == 'status':
                # Intention is that the stop / abort task
                # is pushed by itself.
                if n_tasks > 1:
                    q.put(tasks[:-1])
                    q.put([tasks[-1]])
                else:
                    q.put([tasks[0]])
            else:
                q.put(tasks)
        else:
            time.sleep(wait)
            if wait < 20:
                wait += 0.25


def init_task_querier(worker_pool: str, batch_size: int, threshold: int, work_type: int, 
                      timeout: float = 120, retry_threshold: int = 0):
    # wait_info = WaitInfo
    t = threading.Thread(target=query_tasks_n, args=(batch_size, threshold, work_type,
                         worker_pool, timeout, retry_threshold, _q))
    t.start()


def get_tasks_n(msg_delimiter: str = '|', list_delimiter: str = ';'):
    global _q
    msg_maps = _q.get(True)
    msgs = []
    for msg_map in msg_maps:
        items = [msg_map['type'], msg_map['payload']]
        if msg_map['type'] == 'work':
            items.append(str(msg_map['eq_task_id']))
        msgs.append(msg_delimiter.join(items))

    return list_delimiter.join(msgs)


def stop_task_querier():
    global _go
    _go = False
