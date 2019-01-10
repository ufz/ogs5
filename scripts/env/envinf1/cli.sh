module use /global/apps/modulefiles

module load cmake
module load gcc/6.2.0-1
module load openmpi/gcc/1.8.8-1
module load ninja/1.8.2

# Libraries
module load boost/1.55.0-4
module load doxygen/1.8.7-1_gcc_4.8.1
module load lapack/gcc/3.7.0-1
module load petsc/3.7.6_maint_gcc6.2.0_openmpi_gcc_1.8.8-1

# Tools
module load openmpi/gcc/1.8.8-1
module load numdiff/5.8.1-1
module rm gmp/5.1.2-1
module load gmp-shared/6.1.2-1

export PATH=$PATH:/data/ogs/phreeqc/bin
