#ifndef SPLIT_MPI_COMMUNICATOR_H
#define SPLIT_MPI_COMMUNICATOR_H

#include <mpi.h>

class SplitMPI_Communicator
{
public:
	static bool CreateCommunicator(MPI_Comm Comm, int np, int nb_ddc);
	static void Call_IPhreeqc(void);
};

extern MPI_Comm comm_DDC;
#endif
