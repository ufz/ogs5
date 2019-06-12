module use /global/apps/modulefiles

module load foss/2018b
module load petsc-bilke/3.11.0-foss2018b

module load cmake
module load ninja/1.9.0

# Libraries
module load boost/1.55.0-4
module load doxygen/1.8.7-1_gcc_4.8.1
module load lapack/gcc/3.7.0-1

# Tools
module load numdiff/5.9_gcc_6.2.0_1
#module rm gmp/5.1.2-1
module load gmp-shared/6.1.2-1

export PATH=$PATH:/data/ogs/phreeqc/bin
