/*
 * GaussAlgorithm.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef GAUSSALGORITHM_H_
#define GAUSSALGORITHM_H_

#include <cstddef>
#include <cmath>

// BaseLib
#include "swap.h"

#include "Matrix.h"
#include "DenseDirectLinearSolver.h"
#include "TriangularSolve.h"

namespace MathLib
{
/**
 * This is a class for the direct solution of (dense) systems of
 * linear equations, \f$A x = b\f$. During the construction of
 * the object the matrix A is factorized in matrices L and U using
 * Gauss-Elimination with partial pivoting (rows are exchanged). In doing so
 * the entries of A change! The solution for a specific
 * right hand side is computed by the method execute().
 */
template <typename FLOAT_TYPE>
class GaussAlgorithm : public MathLib::DenseDirectLinearSolver
{
public:
	/**
	 * A direct solver for the (dense) linear system \f$A x = b\f$.
	 * @param A at the beginning the matrix A, at the end of the construction
	 * of the object the matrix contains the factor L (without the diagonal)
	 * in the strictly lower part and the factor U in the upper part.
	 * The diagonal entries of L are all 1.0 and are not explicitly stored.
	 * Attention: the given matrix will be destroyed!
	 * @return a object of type GaussAlgorithm
	 */
	GaussAlgorithm(Matrix<FLOAT_TYPE> &A) :
		_mat (A), _n(_mat.getNRows()), _perm (new size_t [_n])
	{
		size_t k, i, j, nr (_mat.getNRows()), nc(_mat.getNCols());

		for (k = 0; k < nc; k++) {
			// search pivot
			FLOAT_TYPE t = fabs(_mat(k, k));
			_perm[k] = k;
			for (i = k + 1; i < nr; i++) {
				if (fabs(_mat(i,k)) > t) {
					t = _mat(i,k);
					_perm[k] = i;
				}
			}

			// exchange rows
			if (_perm[k] != k) {
				for (j = 0; j < nc; j++) {
					BASELIB::swap (_mat(_perm[k],j), _mat(k,j));
				}
			}

			// eliminate
			for (i = k + 1; i < nr; i++) {
				const FLOAT_TYPE l (_mat(i,k) / _mat(k,k));
				for (j = k; j < nc; j++) {
					_mat(i,j) -= _mat(k,j) * l;
				}
				_mat(i,k) = l;
			}
		}
	}
	/**
	 * destructor, deletes the permutation
	 */
	~GaussAlgorithm()
	{
		delete [] _perm;
	}

	/**
	 * Method solves the linear system \f$A x = b\f$ (based on the LU factorization)
	 * using forward solve and backward solve
	 * @param b at the beginning the right hand side, at the end the solution
	 */
	void execute (FLOAT_TYPE* b) const
	{
		permuteRHS (b);
		forwardSolve<FLOAT_TYPE>(_mat, b); // L z = b, b will be overwritten by z
		backwardSolve<FLOAT_TYPE>(_mat, b); // U x = z, b (z) will be overwritten by x
	}

private:
	/**
	 * permute the right hand side vector according to the
	 * row permutations of the LU factorization
	 * @param b the entries of the vector b are permuted
	 */
	void permuteRHS (FLOAT_TYPE* b) const
	{
		for (size_t i = 0; i < _n; i++) {
			if (_perm[i] != i) {
				BASELIB::swap(b[i], b[_perm[i]]);
			}
		}
	}

	/**
	 * a reference to the matrix
	 */
	Matrix<FLOAT_TYPE>& _mat;

	/**
	 * the size of the matrix
	 */
	size_t _n;

	/**
	 * the permutation of the rows
	 */
	size_t* _perm;
};
} // end namespace MathLib

#endif /* GAUSSALGORITHM_H_ */
