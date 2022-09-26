
# EQ-SQL eq.py
from random import random
import threading
import traceback
import logging
import time
import json
from datetime import datetime, timezone
from enum import IntEnum
from typing import Tuple, Dict

import db_tools
from db_tools import Q


class ResultStatus(IntEnum):
    SUCCESS = 0
    FAILURE = 1


EQ_ABORT = 'EQ_ABORT'
EQ_TIMEOUT = 'EQ_TIMEOUT'
EQ_STOP = 'EQ_STOP'

ABORT_JSON_MSG = json.dumps({'type': 'status', 'payload': EQ_ABORT})

p = None
aborted = False
wait_info = None
# The psycopg2 handle:
DB = None
logger = None


class WaitInfo:

    def __init__(self):
        self.wait = 4

    def getWait(self):
        if self.wait < 60:
            self.wait += 1
        return self.wait


class ThreadRunner(threading.Thread):

    def __init__(self, runnable):
        threading.Thread.__init__(self)
        self.runnable = runnable
        self.exc = "Exited normally"

    def run(self):
        try:
            self.runnable.run()
        except BaseException:
            # tuple of type, value and traceback
            self.exc = traceback.format_exc()


def init(retry_threshold=0, log_level=logging.WARN):
    """Initializes the eq module by connecting to the DB,
    and setting up logging.

    Args:
        log_level: the logging threshold level.
    """
    global logger
    logger = db_tools.setup_log(__name__, log_level)

    global DB
    if DB is not None:
        return

    retries = 0
    while True:
        try:
            DB = db_tools.setup_db(log_level=log_level, envs=True)
            DB.connect()
            break
        except db_tools.ConnectionException as e:
            retries += 1
            if retries > retry_threshold:
                raise(e)
            time.sleep(random() * 4)

    return DB


def close():
    """Closes the DB connection. eq.init() is required to re-initilaize the connection
    before calling any other functions.
    """
    global DB
    if DB is not None:
        DB.close()
        DB = None


def validate():
    """ Connect to DB or die! """
    global DB
    # This code has no effect except to validate the connection:
    try:
        DB.execute("select * from emews_id_generator;")
        DB.get()
    except Exception:
        logger.error("ERROR: eq.validate() failed!")
        return None
    return "EQ-SQL:OK"


def _sql_pop_out_q(eq_type) -> str:
    """
    Generate sql for a queue pop from emews_queue_out
    From:
    https://www.2ndquadrant.com/en/blog/what-is-select-skip-locked-for-in-postgresql-9-5
    Can only select 1 column from the subquery,
    but we return * from the deleted row.
    See workflow.sql for the returned queue row
    """
    code = """
    DELETE FROM emews_queue_OUT
    WHERE  eq_task_id = (
    SELECT eq_task_id
    FROM emews_queue_OUT
    WHERE eq_task_type = {}
    ORDER BY eq_priority DESC, eq_task_id ASC
    FOR UPDATE SKIP LOCKED
    LIMIT 1
    )
    RETURNING *;
    """.format(eq_type)
    return code


def _sql_pop_in_q(eq_task_id) -> str:
    """
    Generate sql for a queue pop from emewws_queue_in
    From:
    https://www.2ndquadrant.com/en/blog/what-is-select-skip-locked-for-in-postgresql-9-5
    Can only select 1 column from the subquery,
    but we return * from the deleted row.
    See workflow.sql for the returned queue row
    """
    code = """
    DELETE FROM emews_queue_IN
    WHERE  eq_task_id = (
    SELECT eq_task_id
    FROM emews_queue_IN
    WHERE eq_task_id = {}
    ORDER BY eq_task_id
    FOR UPDATE SKIP LOCKED
    LIMIT 1
    )
    RETURNING *;
    """.format(eq_task_id)
    return code


def pop_out_queue(eq_type: int, delay, timeout) -> Tuple:
    """Pops the highest priority task of the specified work type off
    of the db out queue.

    This call repeatedly polls for a task of the specified type. The polling
    interval is specified by
    the delay such that the first interval is defined by the initial delay value
    which is increased exponentionally after the first poll. The polling will
    timeout after the amount of time specified by the timout value is has elapsed.

    Args:
        eq_type: the type of the work to pop from the queue
        delay: the initial polling delay value
        timeout: the duration after which this call will timeout
            and return.

    Returns: A two element tuple where the first elements is one of
        ResultStatus.SUCCESS or ResultStatus.FAILURE. On success the
        second element will be the popped eq task id. On failure, the second
        element will be one of EQ_ABORT or EQ_TIMEOUT depending on the
        cause of the failure.
    """
    sql_pop = _sql_pop_out_q(eq_type)
    res = _queue_pop(sql_pop, delay, timeout)
    logger.debug(f'pop_out_queue: {str(res)}')
    return res


def pop_in_queue(eq_task_id: int, delay, timeout) -> Tuple:
    """Pops the specified task off of the db in queue.

    This call repeatedly polls for a task with specified id. The polling
    interval is specified by
    the delay such that the first interval is defined by the initial delay value
    which is increased exponentionally after the first poll. The polling will
    timeout after the amount of time specified by the timout value is has elapsed.

    Args:
        eq_task_id: id of the task to pop
        delay: the initial polling delay value
        timeout: the duration after which this call will timeout
            and return.
    Returns: A two element tuple where the first elements is one of
        ResultStatus.SUCCESS or ResultStatus.FAILURE. On success the
        second element will be eq task id. On failure, the second
        element will be one of EQ_ABORT or EQ_TIMEOUT depending on the
        cause of the failure.
    """
    sql_pop = _sql_pop_in_q(eq_task_id)
    res = _queue_pop(sql_pop, delay, timeout)
    logger.debug(f'pop_in_queue: {str(res)}')
    return res


def _queue_pop(sql_pop: str, delay: float, timeout: float) -> Tuple[ResultStatus, str]:
    """Performs the actual queue pop as defined the sql string.

    This call repeatedly attempts the pop operation by executing sql until
    the operation completes or the timeout duration has passed. The polling
    interval is specified by
    the delay such that the first interval is defined by the initial delay value
    which is increased after the first poll. The polling will
    timeout after the amount of time specified by the timout value is has elapsed.

    Args:
        sql_pop: the sql query that defines the pop operation
        delay: the initial polling delay value
        timeout: the duration after which this call will timeout
            and return.

    Returns: A two element tuple where the first elements is one of
        ResultStatus.SUCCESS or ResultStatus.FAILURE. On success the
        second element will be the popped eq_task_id. On failure, the second
        element will be one of EQ_ABORT or EQ_TIMEOUT depending on the
        cause of the failure.
    """
    global DB
    start = time.time()
    try:
        while True:
            DB.execute(sql_pop)
            rs = DB.get()
            if rs is not None:
                break  # got good data
            if time.time() - start > timeout:
                return (ResultStatus.FAILURE, EQ_TIMEOUT)
            time.sleep(delay)
            if delay < 30:
                delay += 0.25
    except Exception as e:
        logger.error(f'queue_pop error: {e}')
        logger.error(f'queue_pop error {traceback.format_exc()}')
        return (ResultStatus.FAILURE, EQ_ABORT)

    return (ResultStatus.SUCCESS, rs[1])


def push_out_queue(eq_task_id, eq_type, priority=0) -> ResultStatus:
    """Pushes the specified task onto the output queue with
    the specified priority.

    Args:
        eq_task_id: the id of the task
        eq_type: the type of the task
        priority: the priority of the task
    Returns:
        ResultStatus.SUCCESS if the task was successfully pushed
        onto the output queue, otherwise ResultStatus.FAILURE.
    """
    try:
        # queue_push("emews_queue_OUT", eq_type, eq_task_id, priority)
        DB.insert('emews_queue_OUT', ["eq_task_type", "eq_task_id", "eq_priority"],
                  [eq_type, eq_task_id, priority])
        return ResultStatus.SUCCESS

    except Exception as e:
        logger.error(f'push_out_queue error: {e}')
        logger.error(f'push_out_queue error {traceback.format_exc()}')
        return ResultStatus.FAILURE


def push_in_queue(eq_task_id, eq_type) -> ResultStatus:
    """Pushes the specified task onto the input queue.

    Args:
        eq_task_id: the id of the task
        eq_type: the type of the task
    Returns:
        ResultStatus.SUCCESS if the task was successfully pushed
        onto the input queue, otherwise ResultStatus.FAILURE.
    """
    try:
        DB.insert('emews_queue_IN', ["eq_task_type", "eq_task_id"],
                  [eq_type, eq_task_id])
        return ResultStatus.FAILURE
    except Exception as e:
        logger.error(f'push_in_queue error: {e}')
        logger.error(f'push_in_queue error {traceback.format_exc()}')
        return ResultStatus.FAILURE


def insert_task(exp_id: str, eq_type: int, payload: str) -> Tuple:
    """Inserts the specified payload to the database, creating
    a task entry for it and returning its assigned task id

    Args:
        exp_id: the id of the experiment that this task is part of
        eq_type: the work type of this task
        payload: the task payload

    Returns:
        A tuple whose first element is the ResultStatus of the insert, and
        whose second element is the task id assigned to this task if the insert
        was successfull, otherwise EQ_ABORT.
    """
    try:
        global DB
        DB.execute("select nextval('emews_id_generator');")
        rs = DB.get()
        eq_task_id = rs[0]
        ts = datetime.now(timezone.utc).astimezone().isoformat()
        DB.insert("eq_tasks", ["eq_task_id", "eq_task_type", "json_out", "time_created"],
                  [eq_task_id, eq_type, Q(payload), Q(ts)])
        DB.insert("eq_exp_id_tasks", ["exp_id", "eq_task_id"],
                  [Q(exp_id), eq_task_id])
    except Exception as e:
        logger.error(f'insert_task error: {e}')
        logger.error(f'insert_task error {traceback.format_exc()}')
        return (ResultStatus.FAILURE, EQ_ABORT)

    return (ResultStatus.SUCCESS, eq_task_id)


def select_task_payload(eq_task_id: int) -> Tuple[ResultStatus, str]:
    """Selects the 'json_out' payload associated with the specified task id in
    the eq_tasks table, setting the start time of the task to
    the current time.

    Args:
        eq_task_id: the id of the task to get the json_out for

    Returns:
        A tuple containing the ResultStatus as its first element,
        and if successful the json_out payload
        for the specified task id as its second, otherwise the second element will
        be EQ_ABORT.
    """
    try:
        global DB
        DB.select("eq_tasks", "json_out", f'eq_task_id={eq_task_id}')
        rs = DB.get()
        ts = datetime.now(timezone.utc).astimezone().isoformat()
        DB.update("eq_tasks", ['time_start'], [Q(ts)], where=f'eq_task_id={eq_task_id}')
        result = rs[0]
    except Exception as e:
        logger.error(f'select_task_payload error: {e}')
        logger.error(f'select_task_payload error {traceback.format_exc()}')
        return (ResultStatus.FAILURE, EQ_ABORT)

    return (ResultStatus.SUCCESS, result)


def select_task_result(eq_task_id: int) -> Tuple[ResultStatus, str]:
    """Selects the result ('json_in') payload associated with the specified task id in
    the eq_tasks table.

    Args:
        eq_task_id: the id of the task to get the json_in for

    Returns:
        A tuple containing the ResultStatus, and if successful the result payload
        for the specified task id, otherwise EQ_ABORT.
    """
    try:
        global DB
        DB.select("eq_tasks", "json_in", f'eq_task_id={eq_task_id}')
        rs = DB.get()
        result = rs[0]
    except Exception as e:
        logger.error(f'select_task_result error: {e}')
        logger.error(f'select_task_result error {traceback.format_exc()}')
        return (ResultStatus.FAILURE, EQ_ABORT)

    return (ResultStatus.SUCCESS, result)


def update_task(eq_task_id: int, payload: str) -> ResultStatus:
    """Updates the specified task in the eq_tasks table with the specified
    result ('json_in') payload. This also updates the 'time_stop'
    to the time when the update occurred.

    Args:
        eq_task_id: the id of the task to update
        payload: the payload to update the task with
    Returns:
        ResultStatus.SUCCESS if the task was successfully updated, otherwise
        ResultStatus.FAILURE.
    """
    global DB
    ts = datetime.now(timezone.utc).astimezone().isoformat()
    try:
        DB.update("eq_tasks", ["json_in", 'time_stop'], [Q(payload), Q(ts)],
                  where=f'eq_task_id={eq_task_id}')
        return ResultStatus.SUCCESS
    except Exception as e:
        logger.error(f'update_task error: {e}')
        logger.error(f'update_task error {traceback.format_exc()}')
        return ResultStatus.FAILURE


def stop_worker_pool(eq_type: int) -> ResultStatus:
    """Stops any workers pools associated with the specified work type by
    pusing EQ_STOP into the queue.

    Args:
        eq_type: the work type for the pools to stop
    Returns:
        ResultStatus.SUCCESS if the stop message was successfully pushed, otherwise
        ResultStatus.FAILURE.
    """
    try:
        global DB
        DB.execute("select nextval('emews_id_generator');")
        rs = DB.get()
        eq_task_id = rs[0]
        DB.insert("eq_tasks", ["eq_task_id", "eq_task_type", "json_out"],
                  [eq_task_id, eq_type, Q("EQ_STOP")])
    except Exception as e:
        logger.error(f'stop_worker_pool error: {e}')
        logger.error(f'stop_worker_pool error {traceback.format_exc()}')
        return ResultStatus.FAILURE

    result_status = push_out_queue(eq_task_id, eq_type, priority=-1)
    return result_status


def query_task(eq_type: int, delay: float = 0.5, timeout: float = 2.0) -> Dict:
    """Queries for the highest priority task of the specified type.

    The query repeatedly polls for a task. The polling interval is specified by
    the delay such that the first interval is defined by the initial delay value
    which is increased exponentionally after the first poll. The polling will
    timeout after the amount of time specified by the timout value is has elapsed.

    Args:
        eq_type: the type of the task to query for
        delay: the initial polling delay value
        timeout: the duration after which the query will timeout

    Returns:
        A dictionary formatted message. If the query results in a
        status update, the dictionary will have the following format:
        {'type': 'status', 'payload': P} where P is one of 'EQ_STOP',
        'EQ_ABORT', or 'EQ_TIMEOUT'. If the query finds work to be done
        then the dictionary will be:  {'type': 'work', 'eq_task_id': eq_task_id,
        'payload': P} where P is the parameters for the work to be done.
    """
    status, result = pop_out_queue(eq_type, delay, timeout)
    logger.info(f'MSG: {status} {result}')
    if status == ResultStatus.SUCCESS:
        eq_task_id = result
        status, payload = select_task_payload(eq_task_id)
        if status == ResultStatus.SUCCESS:
            if payload == EQ_STOP:
                return {'type': 'status', 'payload': EQ_STOP}
            else:
                return {'type': 'work', 'eq_task_id': eq_task_id, 'payload': payload}
        else:
            return {'type': 'status', 'payload': payload}
    else:
        return {'type': 'status', 'payload': result}


def submit_task(exp_id: str, eq_type: int, payload: str, priority: int = 0) -> Tuple:
    """Submits work of the specified type and priority with the specified
    payload, returning the task id assigned to that task.

    Args:
        exp_id: the id of the experiment of which the work is part.
        eq_type: the type of work
        payload: the work payload
        priority: the priority of this work

    Returns:
        A tuple whose first element is the ResultStatus of the submission, and
        whose second element is the task id assigned to this task if the submission
        was successfull, otherwise -1.
    """
    status, eq_task_id = insert_task(exp_id, eq_type, payload)
    if status == ResultStatus.SUCCESS:
        status = push_out_queue(eq_task_id, eq_type, priority)
        if status == ResultStatus.SUCCESS:
            return (status, eq_task_id)

    return (ResultStatus.FAILURE, -1)


def report_task(eq_task_id: int, eq_type: int, result: str) -> ResultStatus:
    """Reports the result of the specified task of the specified type

    Args:
        eq_task_id: the id of the task whose results are being reported.
        eq_type: the type of the task whose results are being reported.
        result: the result of the task.
    Returns:
        ResultStatus.SUCCESS if the task was successfully reported, otherwise
        ResultStatus.FAILURE.
    """
    result_status = update_task(eq_task_id, result)
    if result_status == ResultStatus.SUCCESS:
        return push_in_queue(eq_task_id, eq_type)
    else:
        return result_status


def query_result(eq_task_id: int, delay: float = 0.5, timeout: float = 2.0) -> Tuple:
    """Queries for the result of the specified task.

    The query repeatedly polls for a result. The polling interval is specified by
    the delay such that the first interval is defined by the initial delay value
    which is increased exponentionally after the first poll. The polling will
    timeout after the amount of time specified by the timout value is has elapsed.

    Args:
        eq_task_id: the id of the task to query
        delay: the initial polling delay value
        timeout: the duration after which the query will timeout

    Returns:
        A tuple whose first element indicates the status of the query:
        ResultStatus.SUCCESS or ResultStatus.FAILURE, and whose second element
        is either the result of the task, or in the case of failure the reason
        for the failure (EQ_TIMEOUT, or EQ_ABORT)
    """
    msg = pop_in_queue(eq_task_id, delay, timeout)
    if msg[0] != ResultStatus.SUCCESS:
        return msg

    return select_task_result(eq_task_id)
