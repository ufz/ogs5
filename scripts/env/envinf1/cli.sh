module () { eval `/usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd sh $*`; }
export MODULEPATH=$MODULEPATH:/global/apps/modulefiles

module load cmake/3.1.3-1
module load gcc/6.2.0-1
module unload python
module load openmpi/gcc/1.8.8-1

# Libraries
module load boost/1.55.0-4
module load doxygen/1.8.7-1_gcc_4.8.1
module load lapack/gcc/3.7.0-1
module load petsc/3.7.6_maint_gcc6.2.0_openmpi_gcc_1.8.8-1

# Tools
module unload python
module load openmpi/gcc/1.8.8-1
module load numdiff/5.8.1-1
module rm gmp/5.1.2-1
module load gmp-shared/6.1.2-1

export PATH=$PATH:/data/ogs/phreeqc/bin
ln -s /opt/ogs/ogs5-libs Libs || :
