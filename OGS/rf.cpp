/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: rf.c
 */
/* Aufgabe:
   ROCKFLOW-FEM - Hauptprogramm
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   06/1998     AH         Konfigurationsdatei
   08/1999     OK         RF-FEM Applikation
   10/1999     AH         Systemzeit

   last modified: OK 14.12.1999
 */
/**************************************************************************/

/**
 * the preprocessor directive RFW_FRACTURE is only useable until version 4.11 of OGS
 * */
#include "BuildInfo.h"

#include <iostream>

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
#include "par_ddc.h"
#include <mpi.h>
#endif
#ifdef LIS
#include "lis.h"
#include <omp.h>
#endif

/* Preprozessor-Definitionen */
#include "makros.h"
#include "display.h"
#include "memory.h"
#include "ogs_display.h"
#define TEST
/* Benutzte Module */
#include "break.h"
#include "timer.h"
// 16.12.2008. WW #include "rf_apl.h"
#include "FileTools.h"
#include "files0.h"
#ifdef SUPERCOMPUTER
// kg44 test for buffered outputh
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#endif
#include "problem.h"

/* Deklarationen */
int main(int argc, char* argv[]);
void ShowSwitches(void);
// LB,string FileName; //WW
// LB,string FilePath; //23.02.2009. WW
// ------  12.09.2007 WW:
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
double elapsed_time_mpi;
#include "SplitMPI_Communicator.h"
// ------
#endif

// Use PETSc. WW
#ifdef USE_PETSC
#include "petscksp.h"
#if (PETSC_VERSION_NUMBER > 3030)
#include "petsctime.h"
#endif
#endif

// declaration of defaultOutputPath
#include "rf_out_new.h"

/* Definitionen */

/**************************************************************************/
/* ROCKFLOW - Funktion: main
 */
/* Aufgabe:
   Hauptprogramm
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int argc: Anzahl der Kommandozeilenparameter (incl. Programmname)
   E char *argv[]: Zeiger auf Feld der argc Kommandozeilenparameter
 */
/* Ergebnis:
   Fehlerfreie Bearbeitung: Exit-Code 0
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   08/1999     OK         RF-FEM Applikation
 */
/**************************************************************************/
int main(int argc, char* argv[])
{
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
	int nb_ddc = 0; // number of cores for DDC related processes
#endif
#ifdef SUPERCOMPUTER
	// *********************************************************************
	// buffered output ... important for performance on cray
	// (unbuffered output is limited to 10 bytes per second)
	// georg.kosakowski@psi.ch 11.10.2007

	char buf[1024 * 1024];
	int bsize;

	bsize = 1024 * 1024; // question: what happens if buffer is full?
	// according to documentation the buffer is flushed when full.
	// If we have a lot of output, increasing buffer is usefull.
	if (bsize > 0)
		//        bufstd = malloc(bsize);
		setvbuf(stdout, buf, _IOFBF, bsize);
//**********************************************************************
#endif
/*---------- MPI Initialization ----------------------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
	printf("Before MPI_Init\n");
#if defined(USE_MPI_GEMS)
	int prov;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);
#else
	MPI_Init(&argc, &argv);
#endif
	MPI_Barrier(MPI_COMM_WORLD); // 12.09.2007 WW
	elapsed_time_mpi = -MPI_Wtime(); // 12.09.2007 WW
	bool splitcomm_flag;
	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	splitcomm_flag = SplitMPI_Communicator::CreateCommunicator(MPI_COMM_WORLD, np, nb_ddc);
	time_ele_paral = 0.0;
#endif
/*---------- MPI Initialization ----------------------------------*/

#ifdef USE_PETSC
	int rank, r_size;
	PetscLogDouble v1, v2;
	char help[] = "OGS with PETSc \n";
	// PetscInitialize(argc, argv, help);
	PetscInitialize(&argc, &argv, (char*)0, help);
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
	PetscTime(&v1);
#else
	PetscGetTime(&v1);
#endif
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &r_size);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nUse PETSc solver");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", r_size, rank);
#endif

/*---------- LIS solver -----------------------------------------*/
#ifdef LIS
	// Initialization of the lis solver.
	lis_initialize(&argc, &argv);
#endif

	std::vector<std::string> arg_strings;
	if (argc > 1)
	{
		for (int i = 1; i < argc; i++)
		{
			arg_strings.push_back(std::string(argv[i]));
		}
	}
	else
	{
#if defined(USE_MPI) // WW
		if (myrank == 0)
#endif
#ifdef USE_PETSC
		if (rank == 0)
#endif
		DisplayStartMsg();

		std::string s_buff;
		std::stringstream ss;
		getline(std::cin, s_buff);
		ss.str(s_buff);
		while (!ss.eof())
		{
			ss >> s_buff;
			arg_strings.push_back(s_buff);
		}
	}

	for (std::size_t i = 0; i < arg_strings.size(); i++)
	{
                const std::string anArg = arg_strings[i];
		if (anArg == "--help" || anArg == "-h")
		{
			std::cout << "Usage: ogs [MODEL_ROOT] [OPTIONS]\n"
			          << "Where OPTIONS are:\n"
			          << "  -h [--help]               print this message and exit\n"
			          << "  -b [--build-info]         print build info and exit\n"
			          << "  --output-directory DIR    put output files into DIR\n"
			          << "  --version                 print ogs version and exit"
			          << "\n";
			exit(0);
		}
		if (anArg == "--build-info" || anArg == "-b")
		{
			std::cout << "ogs version: " << BuildInfo::OGS_VERSION << "\n"
			          << "ogs date: " << BuildInfo::OGS_DATE << "\n";
			std::cout << "git commit info: " << BuildInfo::GIT_COMMIT_INFO << "\n";
			std::cout << "build timestamp: " << BuildInfo::BUILD_TIMESTAMP << "\n";
			exit(0);
		}
		if (anArg == "--version")
		{
			std::cout << BuildInfo::OGS_VERSION << "\n";
			exit(0);
		}
		if (anArg == "--output-directory")
		{
			if (i + 1 >= arg_strings.size())
			{
				std::cerr << "Error: Parameter " << anArg << " needs an additional argument" << std::endl;
				std::exit(EXIT_FAILURE);
			}
			std::string path = arg_strings[++i];

			if (!path.empty())
				defaultOutputPath = path;
			continue;
		}
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
		std::string decompositions;
		if (anArg == "--domain-decomposition" || anArg == "-ddc")
		{
			if (i + 1 >= arg_strings.size())
			{
				std::cerr << "Error: Parameter " << anArg << " needs an additional argument" << std::endl;
				std::exit(EXIT_FAILURE);
			}
			decompositions = std::string(arg_strings[++i]);
			nb_ddc = atoi(decompositions.c_str());
			continue;
		}
#endif
                else
                {
                    FileName = arg_strings[i];
                }
	} // end of parse argc loop

	if (argc > 1)
	{
#if defined(USE_MPI) // WW
		if (myrank == 0)
#endif
#ifdef USE_PETSC
		if (rank == 0)
#endif
		DisplayStartMsg();
        }

#ifdef TESTTIME
	TStartTimer(0);
#endif

	FilePath = pathDirname(FileName);

        // If no option is given, output files are placed in the same directory as the input files
	if (defaultOutputPath.empty())
		defaultOutputPath = FilePath;

	std::string solver_pkg_name = BuildInfo::SOLVER_PACKAGE_NAME;
	// No default linear solver package is in use.
	if (solver_pkg_name.find("Default") == std::string::npos)
	{
		std::cout << "\nWarning: " << solver_pkg_name << " other than the OGS default one is in use." << std::endl;
		std::cout << "         The solver setting may need to be adjusted for the solution accuracy!" << std::endl;
	}
    
	// LB Check if file exists
	std::string tmpFilename = FileName;
	tmpFilename.append(".pcs");
	if (!IsFileExisting(tmpFilename))
	{
		std::cout << " Error: Cannot find file " << FileName << "\n";
		return 1;
	}

	// ---------------------------WW
	Problem* aproblem = new Problem(FileName.data());
#ifdef USE_PETSC
	aproblem->setRankandSize(rank, r_size);
#endif
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
	aproblem->setRankandSize(myrank, mysize);

	if (myrank != MPI_UNDEFINED)
	{
#endif
		aproblem->Euler_TimeDiscretize();
		delete aproblem;
		aproblem = NULL;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_KRC)
	}

	// sending killing signals to ranks of group_IPQC, only when the group exists
	if (splitcomm_flag == true)
	{
#ifdef OGS_FEM_IPQC
		int signal = -1, rank_IPQC, mysize_IPQC = np - nb_ddc;
		for (int i = 0; i < mysize_IPQC; i++)
		{
			rank_IPQC = mysize + i;
			MPI_Send(&signal, 1, MPI_INT, rank_IPQC, 0, MPI_COMM_WORLD);
		}
#endif
	}

#endif

	if (ClockTimeVec.size() > 0)
		ClockTimeVec[0]->PrintTimes(); // CB time
	DestroyClockTime();
#ifdef TESTTIME
#if defined(USE_MPI)
	if (myrank == 0)
#endif
#if defined(USE_PETSC)
		if (rank == 0)
#endif
			std::cout << "Simulation time: " << TGetTimer(0) << "s"
			          << "\n";
#endif
/* Abspann ausgeben */
/*--------- MPI Finalize ------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_KRC)
	elapsed_time_mpi += MPI_Wtime(); // 12.09.2007 WW
	std::cout << "\n *** Total CPU time of parallel modeling: " << elapsed_time_mpi << "\n"; // WW
	// Count CPU time of post time loop WW
	MPI_Finalize();
#endif
/*--------- MPI Finalize ------------------*/
/*--------- LIS Finalize ------------------*/
#ifdef LIS
	lis_finalize();
#endif
	/*--------- LIS Finalize ------------------*/

#ifdef USE_PETSC
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
	PetscTime(&v2);
#else
	PetscGetTime(&v2);
#endif

	PetscPrintf(PETSC_COMM_WORLD, "\t\n>>Total elapsed time by using PETSC:%f s\n", v2 - v1);

	PetscFinalize();
#endif

	return 0;
}
