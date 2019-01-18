/*
 * TriangularSolve.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
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
 * @param b at beginning the right hand side vector, at the end the solution
 * @param dim the dimension of the matrix.
 * vector
 */
template <typename Matrix>
void forwardSolve(const Matrix& L, double* b, const std::size_t dim)
{
    for (size_t r = 0; r < dim; r++)
    {
        double t(0.0);
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
 * @param dim the dimension of the matrix.
 */
template <typename Matrix>
void backwardSolve(const Matrix& U, double* z, const std::size_t dim)
{
    const int m = static_cast<int>(dim);
    for (int r = m - 1; r >= 0; r--)
    {
        double t(0.0);
        for (size_t c = r + 1; c < dim; c++)
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
 * @param dim the dimension of the matrix.
 */
template <typename Matrix>
void backwardSolve(Matrix const& mat, double* x, double const* const b,
                   const std::size_t dim)
{
    const int n = static_cast<int>(dim);
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i];
    }

    backwardSolve(mat, x, dim);
}
}  // end namespace MathLib

#endif /* TRIANGULARSOLVE_H_ */
