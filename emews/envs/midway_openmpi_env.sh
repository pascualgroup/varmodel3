# load the gcc used to compile swift-t and R
module load gcc/12.2.0
# load the mpi distribution used by swift-t
module load openmpi/5.0.2+gcc-12.2.0
# module load openmpi/4.1.1

source /project/jozik/ncollier/envs/varmodel-py3.9/bin/activate

SWIFT_T=/project/jozik/sfw/gcc-12.2.0/openmpi-5.0.2/swift-t-04092024
# SWIFT_T=/project/jozik/sfw/gcc-12.2.0/openmpi-4.1.1/swift-t-04092024
R=/project/ahotton/sfw/gcc-12.2.0/R-4.3.3
# Add swift-t and R to the PATH
export PATH=$SWIFT_T/stc/bin:$R/bin:$PATH
export LD_LIBRARY_PATH=$R/lib64/R/lib:$LD_LIBRARY_PATH