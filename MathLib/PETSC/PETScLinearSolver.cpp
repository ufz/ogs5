/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*!
   \brief Definition of member functions of class PETScLinearSolver

   10~11.2011. WW

*/
#include "PETScLinearSolver.h"

#include <iostream>

namespace petsc_group
{
PETScLinearSolver::PETScLinearSolver(const int size)
    : A(NULL), b(NULL), x(NULL), lsolver(NULL), prec(NULL), i_start(0), i_end(0), global_x(NULL)
{
	ltolerance = 1.e-10;
	m_size = size;
	time_elapsed = 0.0;
	d_nz = 10;
	o_nz = 10;
	nz = 10;
	m_size_loc = PETSC_DECIDE;
	mpi_size = 0;
	rank = 0;
}

PETScLinearSolver::~PETScLinearSolver()
{
	VecDestroy(&b);
	VecDestroy(&x);
	MatDestroy(&A);
	if (lsolver)
		KSPDestroy(&lsolver);
	// if(prec) PCDestroy(&prec);

	if (global_x)
		delete[] global_x;

	PetscPrintf(PETSC_COMM_WORLD, "\tNumber of Unknows: %d", m_size);
	PetscPrintf(PETSC_COMM_WORLD, "\n\tElapsed time in linear solver: %f s\n", time_elapsed);
}

void PETScLinearSolver::Init(const int* sparse_index)
{
	if (sparse_index)
	{
		d_nz = sparse_index[0];
		o_nz = sparse_index[1];
		nz = sparse_index[2];
		m_size_loc = sparse_index[3];
	}

	VectorCreate(m_size);
	MatrixCreate(m_size, m_size);

	global_x = new PetscScalar[m_size];
}

/*!
  \brief KSP and PC type

 KSPRICHARDSON "richardson"
 KSPCHEBYCHEV  "chebychev"
 KSPCG         "cg"
 KSPCGNE       "cgne"
 KSPNASH       "nash"
 KSPSTCG       "stcg"
 KSPGLTR       "gltr"
 KSPGMRES      "gmres"
 KSPFGMRES     "fgmres"
 KSPLGMRES     "lgmres"
 KSPDGMRES     "dgmres"
 KSPTCQMR      "tcqmr"
 KSPBCGS       "bcgs"
 KSPIBCGS        "ibcgs"
 KSPBCGSL        "bcgsl"
 KSPCGS        "cgs"
 KSPTFQMR      "tfqmr"
 KSPCR         "cr"
 KSPLSQR       "lsqr"
 KSPPREONLY    "preonly"
 KSPQCG        "qcg"
 KSPBICG       "bicg"
 KSPMINRES     "minres"
 KSPSYMMLQ     "symmlq"
 KSPLCD        "lcd"
 KSPPYTHON     "python"
 KSPBROYDEN    "broyden"
 KSPGCR        "gcr"
 KSPNGMRES     "ngmres"
 KSPSPECEST    "specest"

 PCNONE            "none"
 PCJACOBI          "jacobi"
 PCSOR             "sor"
 PCLU              "lu"
 PCSHELL           "shell"
 PCBJACOBI         "bjacobi"
 PCMG              "mg"
 PCEISENSTAT       "eisenstat"
 PCILU             "ilu"
 PCICC             "icc"
 PCASM             "asm"
 PCGASM            "gasm"
 PCKSP             "ksp"
 PCCOMPOSITE       "composite"
 PCREDUNDANT       "redundant"
 PCSPAI            "spai"
 PCNN              "nn"
 PCCHOLESKY        "cholesky"
 PCPBJACOBI        "pbjacobi"
 PCMAT             "mat"
 PCHYPRE           "hypre"
 PCPARMS           "parms"
 PCFIELDSPLIT      "fieldsplit"
 PCTFS             "tfs"
 PCML              "ml"
 PCPROMETHEUS      "prometheus"
 PCGALERKIN        "galerkin"
 PCEXOTIC          "exotic"
 PCHMPI            "hmpi"
 PCSUPPORTGRAPH    "supportgraph"
 PCASA             "asa"
 PCCP              "cp"
 PCBFBT            "bfbt"
 PCLSC             "lsc"
 PCPYTHON          "python"
 PCPFMG            "pfmg"
 PCSYSPFMG         "syspfmg"
 PCREDISTRIBUTE    "redistribute"
 PCSACUSP          "sacusp"
 PCSACUSPPOLY      "sacusppoly"
 PCBICGSTABCUSP    "bicgstabcusp"
 PCSVD             "svd"
 PCAINVCUSP        "ainvcusp"
 PCGAMG            "gamg"

*/
void PETScLinearSolver::Config(const PetscReal tol, const PetscInt maxits, const KSPType lsol, const PCType prec_type,
                               const std::string& prefix)
{
	ltolerance = tol;
	sol_type = lsol;
	pc_type = prec_type;

	KSPCreate(PETSC_COMM_WORLD, &lsolver);

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 4)
	KSPSetOperators(lsolver, A, A);
#else
	KSPSetOperators(lsolver, A, A, DIFFERENT_NONZERO_PATTERN);
#endif

	KSPSetType(lsolver, lsol);

	KSPGetPC(lsolver, &prec);
	PCSetType(prec, prec_type); //  PCJACOBI); //PCNONE);
	KSPSetTolerances(lsolver, ltolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxits);

	if (!prefix.empty())
	{
		KSPSetOptionsPrefix(lsolver, prefix.c_str());
		PCSetOptionsPrefix(prec, prefix.c_str());
	}

	KSPSetFromOptions(lsolver);
}

//-----------------------------------------------------------------
void PETScLinearSolver::VectorCreate(PetscInt m)
{
	// PetscErrorCode ierr;  // returned value from PETSc functions
	VecCreate(PETSC_COMM_WORLD, &b);
	////VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &b);
	// VecSetSizes(b, m_size_loc, m);
	VecSetSizes(b, PETSC_DECIDE, m);
	VecSetFromOptions(b);
	VecSetOption(b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
	VecSetUp(b); // kg44 for PETSC 3.3
	VecDuplicate(b, &x);
	VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

	VecGetLocalSize(x, &m_size_loc);

	// VecGetOwnershipRange(b, &i_start,&i_end);
}

void PETScLinearSolver::MatrixCreate(PetscInt m, PetscInt n)
{
	MatCreate(PETSC_COMM_WORLD, &A);
	// TEST  MatSetSizes(A, m_size_loc, PETSC_DECIDE, m, n);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, n);
	// MatSetSizes(A, m_size_loc, PETSC_DECIDE, m,  n);

	MatSetType(A, MATMPIAIJ);
	MatSetFromOptions(A);

	MatSeqAIJSetPreallocation(A, d_nz, PETSC_NULL);
	MatMPIAIJSetPreallocation(A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
	MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

	MatSetUp(A); // KG44 this seems to work with petsc 3.3 ..the commands below result in problems when assembling the
	// matrix with version 3.3

	MatGetOwnershipRange(A, &i_start, &i_end);
}

void PETScLinearSolver::getLocalRowColumnSizes(int* m, int* n)
{
	MatGetLocalSize(A, m, n);
}
void PETScLinearSolver::getOwnerRange(int* start_r, int* end_r)
{
	*start_r = i_start;
	*end_r = i_end;
}

void PETScLinearSolver::Solver()
{
// TEST
#ifdef TEST_MEM_PETSC
	PetscLogDouble mem1, mem2;
	PetscMemoryGetCurrentUsage(&mem1);
#endif

	/*
	//TEST
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject)x,"Solution");
	VecView(x, viewer);
	*/

	int its;
	PetscLogDouble v1, v2;
	KSPConvergedReason reason;

// #define PETSC34
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
	PetscTime(&v1);
#else
	PetscGetTime(&v1);
#endif

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 4)
	KSPSetOperators(lsolver, A, A);
#else
	KSPSetOperators(lsolver, A, A, DIFFERENT_NONZERO_PATTERN);
#endif

	KSPSolve(lsolver, b, x);

	KSPGetConvergedReason(lsolver, &reason); // CHKERRQ(ierr);
	if (reason == KSP_DIVERGED_INDEFINITE_PC)
	{
		PetscPrintf(PETSC_COMM_WORLD, "\nDivergence because of indefinite preconditioner;\n");
		PetscPrintf(PETSC_COMM_WORLD, "Run the executable again but with -pc_factor_shift_positive_definite option.\n");
	}
	else if (reason < 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "\nOther kind of divergence: this should not happen.\n");
	}
	else
	{
		const char* slv_type;
		const char* prc_type;
		KSPGetType(lsolver, &slv_type);
		PCGetType(prec, &prc_type);

		PetscPrintf(PETSC_COMM_WORLD, "\n================================================");
		PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner", slv_type, prc_type);
		KSPGetIterationNumber(lsolver, &its); // CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD, "\nConvergence in %d iterations.\n", (int)its);
		PetscPrintf(PETSC_COMM_WORLD, "\n================================================");
	}
	PetscPrintf(PETSC_COMM_WORLD, "\n");

// VecAssemblyBegin(x);
// VecAssemblyEnd(x);

// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
	PetscTime(&v2);
#else
	PetscGetTime(&v2);
#endif

	time_elapsed += v2 - v1;

#define nonTEST_OUT
#ifdef TEST_OUT
	// TEST
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_vector.txt", &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject)A, "Matrix");
	MatView(A, viewer);
	PetscObjectSetName((PetscObject)x, "Solution");
	VecView(x, viewer);
	PetscObjectSetName((PetscObject)b, "RHS");
	VecView(b, viewer);
/*
VecDestroy(&b);
VecDestroy(&x);
MatDestroy(&A);
if(lsolver) KSPDestroy(&lsolver);
if(global_x)
  delete []  global_x;

 PetscFinalize();
 exit(0);
 */
#endif

#ifdef TEST_MEM_PETSC
	// TEST
	PetscMemoryGetCurrentUsage(&mem2);
	PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2,
	            (int)(mem2 - mem1));
#endif
}

void PETScLinearSolver::AssembleRHS_PETSc()
{
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
}
void PETScLinearSolver::AssembleUnkowns_PETSc()
{
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
}
void PETScLinearSolver::AssembleMatrixPETSc(const MatAssemblyType type)
{
	MatAssemblyBegin(A, type);
	MatAssemblyEnd(A, type);
}

void PETScLinearSolver::UpdateSolutions(PetscScalar* u)
{
#ifdef TEST_MEM_PETSC
	// TEST
	PetscLogDouble mem1, mem2;
	PetscMemoryGetCurrentUsage(&mem1);
#endif

	PetscScalar* xp = NULL; // nullptr;
	VecGetArray(x, &xp);

	gatherLocalVectors(xp, u);

	// MPI_Barrier(PETSC_COMM_WORLD);

	VecRestoreArray(x, &xp);

// TEST
#ifdef TEST_MEM_PETSC
	PetscMemoryGetCurrentUsage(&mem2);
	PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2,
	            (int)(mem2 - mem1));
#endif
}

void PETScLinearSolver::MappingSolution()
{
	UpdateSolutions(global_x);
}

int PETScLinearSolver::GetLocalSolution(PetscScalar* x_l)
{
	PetscInt count;
	VecGetLocalSize(x, &count);

	VecGetArray(x, &x_l);

	return count;
}

int PETScLinearSolver::GetLocalRHS(PetscScalar* rhs_l)
{
	PetscInt count;
	VecGetLocalSize(b, &count);

	VecGetArray(b, &rhs_l);

	return count;
}

double* PETScLinearSolver::GetGlobalSolution() const
{
	return global_x;
}

/*!
  Get values of the specified elements from a global vector

  @param v_type - Indicator for vector: 0: x; 1: rhs
  @param ni 	- number of elements to get
  @param ix 	- indices where to get them from (in global 1d numbering)
*/
void PETScLinearSolver::GetVecValues(const int v_type, PetscInt ni, const PetscInt ix[], PetscScalar y[]) const
{
	if (v_type == 0)
		VecGetValues(x, ni, ix, y);
	else
		VecGetValues(b, ni, ix, y);
}

/*!
    Get norm of RHS
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i|
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormRHS(NormType nmtype)
{
	PetscReal norm = 0.;
	VecNorm(b, nmtype, &norm);
	return norm;
}
/*!
    Get norm of x
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i|
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormX(NormType nmtype)
{
	PetscReal norm = 0.;
	VecNorm(x, nmtype, &norm);
	return norm;
}

void PETScLinearSolver::RestoreLocalSolutionArray(PetscScalar* x_l)
{
	VecRestoreArray(x, &x_l);
}
void PETScLinearSolver::RestoreLocalRHSArray(PetscScalar* rhs_l)
{
	VecRestoreArray(b, &rhs_l);
}

void PETScLinearSolver::set_bVectorEntry(const int i, const double value)
{
	VecSetValues(b, 1, &i, &value, INSERT_VALUES);
}
void PETScLinearSolver::set_xVectorEntry(const int i, const double value)
{
	VecSetValues(x, 1, &i, &value, INSERT_VALUES);
}

void PETScLinearSolver::setArrayValues(int arr_idx, PetscInt ni, const PetscInt ix[], const PetscScalar y[],
                                       InsertMode iora)
{
	if (arr_idx == 0)
		VecSetValues(x, ni, ix, y, iora);
	else if (arr_idx == 1)
		VecSetValues(b, ni, ix, y, iora);
}

void PETScLinearSolver::add_bVectorEntry(const int i, const double value, InsertMode mode)
{
	VecSetValue(b, i, value, mode);
}
void PETScLinearSolver::add_xVectorEntry(const int i, const double value, InsertMode mode)
{
	VecSetValue(x, i, value, mode);
}

void PETScLinearSolver::Initialize()
{
	VecSet(b, 0.0);
	VecSet(x, 0.0);
	MatZeroEntries(A);
}

void PETScLinearSolver::addMatrixEntry(const int i, const int j, const double value)
{
	MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScLinearSolver::addMatrixEntries(const int m, const int idxm[], const int n, const int idxn[],
                                         const PetscScalar v[])
{
	MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
}

void PETScLinearSolver::zeroRows_in_Matrix(const int nrows, const PetscInt* rows)
{
	PetscScalar one = 1.0;
	// Each process indicates only rows it owns that are to be zeroed
	// MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
	if (nrows > 0)
		MatZeroRows(A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
	else
		MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScLinearSolver::gatherLocalVectors(PetscScalar local_array[], PetscScalar global_array[])
{
	// Collect vectors from processors.
	int size_rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

	// number of elements to be sent for each rank
	std::vector<PetscInt> i_cnt(size_rank);
	// offset in the receive vector of the data from each rank
	std::vector<PetscInt> i_disp(size_rank);

	MPI_Allgather(&m_size_loc, 1, MPI_INT, &i_cnt[0], 1, MPI_INT, PETSC_COMM_WORLD);

	// colloect local array
	PetscInt offset = 0;
	for (PetscInt i = 0; i < size_rank; i++)
	{
		i_disp[i] = offset;
		offset += i_cnt[i];
	}

	MPI_Allgatherv(local_array, m_size_loc, MPI_DOUBLE, global_array, &i_cnt[0], &i_disp[0], MPI_DOUBLE,
	               PETSC_COMM_WORLD);
}

void PETScLinearSolver::EQSV_Viewer(std::string file_name)
{
	PetscViewer viewer;
	std::string fname = file_name + "_eqs_dump.txt";
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

	AssembleRHS_PETSc();
	AssembleUnkowns_PETSc();
	AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);

	// PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
	PetscObjectSetName((PetscObject)A, "Stiffness_matrix");
	PetscObjectSetName((PetscObject)b, "RHS");
	PetscObjectSetName((PetscObject)x, "Solution");
	MatView(A, viewer);
	VecView(b, viewer);
	VecView(x, viewer);

//#define  EXIT_TEST
#ifdef EXIT_TEST
	VecDestroy(&b);
	VecDestroy(&x);
	MatDestroy(&A);
	if (lsolver)
		KSPDestroy(&lsolver);
	// if(prec) PCDestroy(&prec);
	if (global_x0)
		delete[] global_x0;
	if (global_x1)
		delete[] global_x1;
	PetscFinalize();
	exit(0);
#endif
}

} // end of namespace
