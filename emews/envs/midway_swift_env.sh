module load mpich/3.3+gcc-9.2.0    
module unload gcc
module load gcc/10.2.0
module unload java
module load R/4.2.0
module load python/anaconda-2021.05

export R_LIBS_SITE=/project2/pascualmm/sfw/r-4.2.0-libs
export PATH=/project2/jozik/sfw/gcc-10.2.0/swift-t-08102022/stc/bin:$PATH
# /project2/pascualmm/sfw/gcc-10.2.0/swift-t-072222/stc/bin:$PATH
