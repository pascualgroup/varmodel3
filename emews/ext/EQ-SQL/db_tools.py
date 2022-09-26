
# DB TOOLS PY
# These tools are independent of the application
#       (they do not know about COVID, Swift, Repast, ...)

import sys

assert sys.version[0] != '2', "This requires Python 3!"

from enum import IntEnum, unique
import logging
import os


def setup_db(log_level=logging.WARN, envs=False):
    """ Convenience function to use from Swift/T """
    if 'DB' not in globals():
        rank = os.getenv('PMIX_RANK')
        print('rank %s Connecting to DB...' % rank)
        global DB
        DB = workflow_sql(log_level=log_level, envs=envs)
    return DB


def setup_log(log_name, log_level, procname=""):
    logger = logging.getLogger(log_name)
    handlr = logging.StreamHandler()
    formtr = logging.Formatter("%(asctime)s " + procname +
                               " %(name)-9s %(message)s",
                               datefmt="%Y-%m-%d %H:%M:%S")
    handlr.setFormatter(formtr)
    logger.addHandler(handlr)
    logger.setLevel(log_level)
    return logger


@unique
class DB_Mode(IntEnum):
    NULL = 0  # undefined state
    OFF  = 1  # do not use the DB
    SOFT = 2  # try DB, run through errors
    ON   = 3  # insist on DB, errors are fatal


class workflow_sql:

    def __init__(self, host="127.0.0.1", port=5432,
                 user=os.environ['USER'],
                 mode=DB_Mode.ON,
                 dbname="EQ_SQL",
                 envs=False,
                 log_level=logging.WARN,
                 procname=""):
        """
        Sets up a wrapper around the SQL connection and cursor objects
        Also caches dicts that convert between names and ids for the
        features and studies tables
        envs: If True, self-configure based on the environment
        """
        self.conn   = None
        self.host   = host
        self.port   = port
        self.mode   = mode  # a DB_Mode
        self.dbname = dbname
        self.user = user
        if envs:
            self.configure_envs()
        self.autoclose = True
        self.procname = procname  # a unique process name
        self.logger = setup_log(__name__, log_level, self.procname)
        self.info("Initialized.")

    def configure_envs(self):
        def env_has(k):
            v = os.getenv(k)
            if v is None:
                return False
            if len(v.strip()) == 0:
                return False
            return True

        if env_has("DB_HOST"):
            self.host = os.getenv("DB_HOST")
        if env_has('DB_USER'):
            self.user = os.getenv('DB_USER')
        if env_has("DB_PORT"):
            try:
                port_string = os.getenv("DB_PORT")
                self.port = int(port_string)
            except ValueError as e:
                self.logger.fatal("DB_PORT is not an integer: " +
                                  "got: '%s'" % port_string)
                raise(e)
        if env_has("DB_MODE"):
            m = os.getenv("DB_MODE")
            try:
                self.mode = DB_Mode[m]
            except Exception:
                raise(ValueError("DB_MODE is unusable: got: '%s'" % m))
        if env_has("DB_NAME"):
            self.dbname = os.getenv("DB_NAME")

    def connect(self):
        import psycopg2
        if self.conn is None:
            if self.mode >= DB_Mode.SOFT:
                self.info("connect(): connecting to {} {} as {} ...".format(self.host, self.port, self.user))
                try:
                    self.conn = psycopg2.connect("dbname="+self.dbname,
                                                 host=self.host,
                                                 port=self.port,
                                                 user=self.user)
                except psycopg2.OperationalError as e:
                    self.info("connect(): could not connect!")
                    raise ConnectionException(e)
                self.info("connect(): connected.")
                self.debug("connect(): " + str(self.conn))
                self.cursor = self.conn.cursor()
                # self.debug("connect(): cursor: " + str(self.cursor))
            else:
                self.info("connect(): DB disabled.")
                self.conn = "DISABLED"
        else:
            if self.conn != "DISABLED":
                self.info("connect(): Already connected.")
        return "OK"

    def insert(self, table, names, values, echo=False):
        """ Do a SQL insert
            return rowid or -1 on SOFT error
        """
        if len(names) != len(values):
            raise ValueError("lengths of names, values must agree!")
        names_tpl  = sql_tuple(names)
        values_tpl = sql_tuple(values)
        cmd = "insert into {} {} values {};" \
            .format(table, names_tpl, values_tpl)
        try:
            self.execute(cmd, echo=echo)
            self.commit()
        except Exception as e:
            if self.mode == DB_Mode.SOFT:
                return -1
            self.info(e)
            raise(e)

    def update(self, table, names, values, where):
        """ Do a SQL update
            return rowid or -1 on SOFT error
        """
        if len(names) != len(values):
            raise ValueError("lengths of names, values must agree!")
        assign_list = []
        for n, v in zip(names, values):
            assign_list.append("%s=%s" % (n, str(v)))
        assigns = ", ".join(assign_list)
        cmd = "update {} set {} where {};" \
            .format(table, assigns, where)
        try:
            self.execute(cmd)
            self.commit()
        except Exception as e:
            if self.mode == DB_Mode.SOFT:
                return -1
            self.info(e)
            raise(e)

    def select(self, table, what, where=None):
        ''' Do a SQL select '''
        cmd = "select %s from %s" % (what, table)
        if where is not None:
            cmd += " where "
            cmd += where
        cmd += ";"
        self.execute(cmd)

    def execute(self, cmd, echo=False):
        self.info(cmd)
        if echo:
            print(cmd)
        try:
            self.cursor.execute(cmd)
            self.conn.commit()
        except Exception as e:
            print(str(e))  # Remove this line after debugging
            if self.mode == DB_Mode.SOFT:
                return
            print(str(e))
            raise(e)

    def executescript(self, cmds):
        try:
            self.cursor.executescript(cmds)
        except Exception as e:
            if self.mode == DB_Mode.SOFT:
                self.debug(e)
                return
            self.info(e)
            raise(e)

    def commit(self):
        """ Should not be called if mode==OFF """
        self.conn.commit()

    def get(self):
        """ Wrapper for fetchone() """
        result = self.cursor.fetchone()
        self.info("get(): %r" % (result is not None))
        return result

    def close(self):
        self.autoclose = False
        if self.mode == DB_Mode.OFF:
            return
        self.conn.close()
        self.conn = None

    def debug(self, message):
        if self.logger:
            self.logger.debug(message)

    def info(self, message):
        if self.logger:
            self.logger.info(message)

    def fatal(self, message):
        if self.logger:
            self.logger.fatal(message)
        else:
            print(message)

    def __del__(self):
        if self.mode == DB_Mode.OFF:
            return
        if not self.autoclose:
            return
        try:
            self.conn.commit()
            self.conn.close()
        except:  # noqa E722
            pass
        self.info("DB auto-closed.")


def workflow_sql_setup(sql, **kwargs):
    """ Setup SQL if not already done """
    if sql is not None:
        sql.logger.debug("workflow_sql_setup(): already setup.")
        return sql
    sql = workflow_sql(**kwargs)
    return sql


def Q(s):
    """ Quote the given string """
    return "'" + str(s) + "'"


def QL(L):
    """ Quote-List: Quote each list entry as a string """
    return map(q, L)


def QA(*args):
    """ Quote-Arguments: Quote each argument as a string,
        return list
    """
    return list(map(q, args))


def sql_tuple(L):
    """ Make the given list into a SQL-formatted tuple """
    L = list(map(str, L))
    result = ""
    result += "("
    result += ",".join(L)
    result += ")"
    return result


def sql_tuple_q(L):
    """ Make the given list into a Quoted SQL-formatted tuple """
    L = list(map(str, L))
    result = ""
    result += "("
    result += ",".join(qL(L))
    result += ")"
    return result


class ConnectionException(Exception):
    def __init__(self, cause):
        """ cause: another Exception """
        self.cause = cause
