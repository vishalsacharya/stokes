module load valgrind
module unload petsc
module unload mpi.intel
module load mpi.intel/4.1
module load petsc/3.4
module unload mpi.intel
module load mpi.intel

make -j

module unload petsc
module unload mpi.intel
module load mpi.intel/4.1
module load petsc/3.5
module unload mpi.intel
module load mpi.intel
