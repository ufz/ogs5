/*
 * TriangularSolve.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef TRIANGULARSOLVE_H_
#define TRIANGULARSOLVE_H_

#include "Matrix.h"

namespace MathLib
{
/**
 * solves the \f$n \times n\f$ triangular linear system \f$L \cdot y = b\f$,
 * assumes \f$L_{ii} = 1.0\f$, \f$i=1,...,n\f$, \f$b\f$ is destroyed
 * @param L the lower triangular matrix
 * @param b at beginning the right hand side vector, at the end the solution vector
 */
template <typename FLOAT_TYPE>
void forwardSolve(const Matrix<FLOAT_TYPE>& L, FLOAT_TYPE* b)
{
	size_t m(L.getNRows());

	for (size_t r = 0; r < m; r++)
	{
		FLOAT_TYPE t(0.0);
		for (size_t c = 0; c < r; c++)
			t += L(r, c) * b[c];
		b[r] = b[r] - t;
	}
}

/**
 * solves the \f$n \times n\f$ triangular linear system \f$U \cdot x=z\f$,
 * \f$U\f$, where \f$U\f$ is a upper triangular matrix.
 * @param U upper triangular matrix
 * @param z at beginning the right hand side, at the end the solution
 */
template <typename FLOAT_TYPE>
void backwardSolve(const Matrix<FLOAT_TYPE>& U, FLOAT_TYPE* z)
{
	size_t m(U.getNRows()), n(U.getNCols());
	for (int r = m - 1; r >= 0; r--)
	{
		FLOAT_TYPE t(0.0);
		for (size_t c = r + 1; c < n; c++)
		{
			t += U(r, c) * z[c];
		}
		z[r] = (z[r] - t) / U(r, r);
	}
}

// backwardSolve mat * x = y, mat ... upper triangular matrix
/**
 * backward solve the system of linear equations \f$ U \cdot x = y\f$,
 * where \f$U\f$ is a upper triangular matrix
 * @param mat the upper triangular matrix
 * @param x the solution of the system of linear equations
 * @param b the right hand side
 */
template <typename FLOAT_TYPE>
void backwardSolve(Matrix<FLOAT_TYPE> const& mat, FLOAT_TYPE* x, FLOAT_TYPE* b)
{
	size_t n_cols(mat.getNCols());
	for (int r = (n_cols - 1); r >= 0; r--)
	{
		FLOAT_TYPE t(0.0);

		for (size_t c = r + 1; c < n_cols; c++)
			t += mat(r, c) * b[c];
		x[r] = (b[r] - t) / mat(r, r);
	}
}
} // end namespace MathLib

#endif /* TRIANGULARSOLVE_H_ */
