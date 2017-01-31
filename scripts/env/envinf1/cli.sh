module () { eval `/usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd sh $*`; }
export MODULEPATH=$MODULEPATH:/global/apps/modulefiles

module load cmake/3.1.3-1
module load gcc/4.8.1-3
module unload python
module load openmpi/gcc/1.8.4-2

# Libraries
module load boost/1.55.0-4
module load doxygen/1.8.7-1_gcc_4.8.1
module load lapack/3.5.0-1_gcc_4.8.1
module load petsc/3.5_maint_gcc_4.8.1-3_openmpi_gcc_1.8.2-1_gcc_4.8.1_CentOS6_envinf

# Tools
module unload python
module load openmpi/gcc/1.8.4-2

export PATH=$PATH:/data/ogs/phreeqc/bin
ln -s /opt/ogs/ogs5-libs Libs || :
