
# ENV LOCAL
# Customize this for your local workstation

## Stan
# Installation locations
ROOT=$HOME/sfw
POSTGRES=$ROOT/postgresql-12.2
SWIFT=$ROOT/swift-t
EQR=$ROOT/EQ-R
R=$ROOT/R-3.6.0/lib64/R

# Environment setup
PATH=$SWIFT/stc/bin:$PATH
PATH=$SWIFT/turbine/bin:$PATH
PATH=$POSTGRES/bin:$PATH
export PYTHONPATH=$EMEWS_PROJECT_ROOT/db
export LD_LIBRARY_PATH=$R/lib:$POSTGRES/lib
MACHINE=""
