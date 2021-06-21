module use /global/apps/modulefiles

module load foss/2019b
module load petsc-bilke/3.11.0-foss2019b

module load git/2.23.0
module load git-lfs/2.7.1

module load CMake/3.15.3
module load ninja

# Libraries
module load boost/1.67.0-1
#module load doxygen/1.8.7-1_gcc_4.8.1
#module load lapack/gcc/3.7.0-1

# Tools
module load numdiff
#module rm gmp/5.1.2-1
#module load gmp-shared/6.1.2-1

export PATH=$PATH:/data/ogs/phreeqc/bin
