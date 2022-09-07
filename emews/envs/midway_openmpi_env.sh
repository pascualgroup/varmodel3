module load gcc/10.1.0
module load openmpi/3.0.0+gcc-10.1.0
module unload java
module load java/1.8
module load openmpi/3.0.0+gcc-10.1.0
module load python/anaconda-2021.05

# R setup
module load openblas/0.3.13 

export LD_LIBRARY_PATH=/software/R-4.1.0-el7-x86_64/lib64/R/lib:/software/R-4.1.0-el7-x86_64/lib64:$LD_LIBRARY_PATH
export R_LIBS_SITE=/project2/pascualmm/sfw/r-4.1.0-libs

export PATH=/project2/jozik/sfw/gcc-10.1.0/swift-t-08252022/stc/bin:$PATH