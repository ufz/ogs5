
#include "SplitMPI_Communicator.h"

#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "makros.h"
#include "rf_react.h"

MPI_Comm comm_DDC;

bool SplitMPI_Communicator::CreateCommunicator(MPI_Comm comm_world, int np,
                                               int nb_ddc)
{
    int n_DDC;
    bool splitcomm;

    if ((nb_ddc > 0) && (nb_ddc < np))
    {  // if the number of total cores is larger than the number of DDCs is the
       // same, two new MPI groups will be
// generated will be generated
#ifdef OGS_FEM_IPQC
        splitcomm = true;
        n_DDC = nb_ddc;  // number of ddc

        int DDC_ranks[n_DDC];
        for (int k = 0; k < n_DDC; k++)
        {
            DDC_ranks[k] = k;
        }

        MPI_Comm comm_IPQC;
        MPI_Group group_base, group_DDC, group_IPQC;

        // define MPI group and communicator for DDC related processes WH
        MPI_Comm_group(comm_world, &group_base);
        MPI_Group_incl(group_base, n_DDC, DDC_ranks,
                       &group_DDC);  // define group flow and mass transport
        MPI_Comm_create(comm_world, group_DDC, &comm_DDC);

        // define MPI group and communicator for IPQC WH
        MPI_Group_difference(group_base, group_DDC, &group_IPQC);
        MPI_Comm_create(comm_world, group_IPQC, &comm_IPQC);

        int myrank_IPQC, mysize_IPQC;
        MPI_Group_size(group_DDC, &mysize);  // WH
        MPI_Group_rank(group_DDC, &myrank);  // WH
        MPI_Group_rank(group_IPQC, &myrank_IPQC);
        MPI_Group_size(group_IPQC, &mysize_IPQC);
        if (myrank_IPQC != MPI_UNDEFINED)  // WH
            std::cout << "After MPI_Init myrank_IPQC = " << myrank_IPQC << '\n';
        if (myrank != MPI_UNDEFINED)  // WH
            std::cout << "After MPI_Init myrank_DDC = " << myrank << '\n';

        if (myrank_IPQC !=
            MPI_UNDEFINED)  // ranks of group_IPQC will call to IPhreeqc
            Call_IPhreeqc();
#endif
    }
    else
    {  // if no -ddc is specified or the number of ddc is incorrect, make ddc =
       // np, no new MPI groups willnot be
        // generated;
        splitcomm = false;
        n_DDC = np;
        comm_DDC = comm_world;
        MPI_Comm_size(comm_DDC, &mysize);
        MPI_Comm_rank(comm_DDC, &myrank);
        std::cout << "After MPI_Init myrank_DDC = " << myrank << '\n';
    }

    return splitcomm;
}

#ifdef OGS_FEM_IPQC
void SplitMPI_Communicator::Call_IPhreeqc(void)
{
    int signal = 1;
    for (;;)
    {
        MPI_Status status;
        // signal: the length of the input character string from group_DDC
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        if (signal < 0)
            break;
        char string_DDC[signal + 1];
        MPI_Recv(string_DDC, signal + 1, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        std::stringstream input, output;
        input << string_DDC;
        // call to IPhreeqc
        Call_IPQC(&input, &output);
        // prepare the output string
        std::string tmp = output.str();
        char message[tmp.length() + 1];
        int strlength = tmp.length();
        for (std::size_t i = 0; i < tmp.length(); i++)
            message[i] = tmp[i];
        message[tmp.length()] = '\0';
        // send the output string back to the relevant rank of group_DDC
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Send(&strlength, 1, MPI_INT, status.MPI_SOURCE, rank,
                 MPI_COMM_WORLD);
        MPI_Send(message, tmp.length() + 1, MPI_CHAR, status.MPI_SOURCE, rank,
                 MPI_COMM_WORLD);
    }

    // terminate the MPI environment for ranks of group_IPQC
    MPI_Finalize();
    exit(0);
}
#endif
