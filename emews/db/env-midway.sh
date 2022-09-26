export DB_HOST=${DB_HOST:-$(hostname --long)} # thetalogin4
# Use python to find a free port
PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
export DB_PORT=${DB_PORT:-$PORT}
# The user name to use for the DB
export DB_USER=${DB_USER:-}
export DB_DATA=${DB_DATA:-/project2/pascualmm/ncollier}
export DB_MODE=ON

DB_ROOT=/project2/pascualmm/sfw/gcc-10.2.0/postgreql-14.2
export PATH=$DB_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$DB_ROOT/lib:$LD_LIBRARY_PATH
