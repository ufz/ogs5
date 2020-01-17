/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*!
   \brief Declaration of class PETScLinearSolver

   11.2011. WW

*/
#ifndef PETSC_LSOLVER_INC
#define PETSC_LSOLVER_INC

#include <string>
#include <vector>
// kg needed for memcpy in petsc libs
#include <stdio.h>
#include <string.h>

#include "petscmat.h"
#include "petscksp.h"

#if (PETSC_VERSION_NUMBER > 3030)
#include "petsctime.h"
#endif

typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;

namespace petsc_group
{
class PETScLinearSolver
{
public:
    PETScLinearSolver(const int size);
    ~PETScLinearSolver();

    void Config(const PetscReal tol, const PetscInt maxits, const KSPType lsol,
                const PCType prec_type, const std::string& prefix = "");

    void Init(const int* sparse_index = NULL);

    void Solver();
    void AssembleRHS_PETSc();
    void AssembleUnkowns_PETSc();
    void AssembleMatrixPETSc(
        const MatAssemblyType type = MAT_FINAL_ASSEMBLY);  // MAT_FLUSH_ASSEMBLY

    void MappingSolution();

    int GetLocalSolution(PetscScalar* x_l);
    int GetLocalRHS(PetscScalar* rhs_l);
    double* GetGlobalSolution() const;
    void GetVecValues(const int v_type, PetscInt ni, const PetscInt ix[],
                      PetscScalar y[]) const;
    PetscReal GetVecNormRHS(NormType nmtype = NORM_2);
    PetscReal GetVecNormX(NormType nmtype = NORM_2);

    void RestoreLocalSolutionArray(PetscScalar* x_l);
    void RestoreLocalRHSArray(PetscScalar* rhs_l);
    void getLocalRowColumnSizes(int* m, int* n);
    void getOwnerRange(int* start_r, int* end_r);

    int Size() const { return m_size; }
    void set_xVectorEntry(const int i, const double value);
    void set_bVectorEntry(const int i, const double value);
    void setArrayValues(int arr_idx, PetscInt ni, const PetscInt ix[],
                        const PetscScalar y[], InsertMode iora = ADD_VALUES);

    void add_xVectorEntry(const int i, const double value, InsertMode mode);
    void add_bVectorEntry(const int i, const double value, InsertMode mode);
    void addMatrixEntry(const int i, const int j, const double value);
    void addMatrixEntries(const int m, const int idxm[], const int n,
                          const int idxn[], const PetscScalar v[]);

    void Initialize();

    void zeroRows_in_Matrix(const int nrow, const PetscInt* rows);
    void zeroMatrix() { MatZeroEntries(A); }
    void set_rank_size(const int m_rank, const int size)
    {
        mpi_size = size;
        rank = m_rank;
    }

    PetscInt getStartRow() const { return i_start; }
    PetscInt getEndRow() const { return i_end; }
    PetscInt getMPI_Size() const { return mpi_size; }
    PetscInt getMPI_Rank() const { return rank; }
    void EQSV_Viewer(std::string file_name);

private:
    PETSc_Mat A;
    PETSc_Vec b;
    PETSc_Vec x;
    KSP lsolver;
    PC prec;
    PetscInt i_start;
    PetscInt i_end;

    PetscScalar* global_x;

    // Slover and preconditioner names, only for log
    std::string sol_type;
    std::string pc_type;

    PetscLogDouble time_elapsed;

    PetscInt m_size;
    PetscInt m_size_loc;
    float ltolerance;
    // Number of nonzeros per row in DIAGONAL portion of
    // local submatrix (same value is used for all local rows)
    PetscInt d_nz;
    // Number of nonzeros per row in the OFF-DIAGONAL portion of
    // local submatrix (same value is used for all local rows).
    PetscInt o_nz;
    // Number of nonzeros per row (same for all rows)
    PetscInt nz;

    int mpi_size;
    int rank;

    void VectorCreate(PetscInt m);
    void MatrixCreate(PetscInt m, PetscInt n);

    /*!
         \brief  collect local vectors
         \param  local_array  local array
         \param  global_array global array
    */
    void gatherLocalVectors(PetscScalar local_array[],
                            PetscScalar global_array[]);

    void UpdateSolutions(PetscScalar* u);
};

// extern std::vector<PETScLinearSolver*> EQS_Vector;
}  // namespace petsc_group
#endif
