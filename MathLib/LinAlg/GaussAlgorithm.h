/*
 * GaussAlgorithm.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
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
template <typename Matrix>
class GaussAlgorithm
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
    GaussAlgorithm(Matrix& A, const std::size_t dim)
        : _mat(A), _n(dim), _perm(std::vector<std::size_t>(dim, 0))
    {
    }

    /**
     * Method solves the linear system \f$A x = b\f$ (based on the LU
     * factorization) using forward solve and backward solve
     * @param b at the beginning the right hand side, at the end the solution
     */
    void execute(double* b)
    {
        const std::size_t nr(_n);
        const std::size_t nc(_n);

        for (std::size_t k = 0; k < nc; k++)
        {
            // search pivot
            double t = fabs(_mat(k, k));
            _perm[k] = k;
            for (std::size_t i = k + 1; i < nr; i++)
            {
                if (fabs(_mat(i, k)) > t)
                {
                    t = _mat(i, k);
                    _perm[k] = i;
                }
            }

            // exchange rows
            if (_perm[k] != k)
            {
                for (std::size_t j = 0; j < nc; j++)
                {
                    BASELIB::swap(_mat(_perm[k], j), _mat(k, j));
                }
            }

            // eliminate
            for (std::size_t i = k + 1; i < nr; i++)
            {
                const double l(_mat(i, k) / _mat(k, k));
                for (std::size_t j = k; j < nc; j++)
                {
                    _mat(i, j) -= _mat(k, j) * l;
                }
                _mat(i, k) = l;
            }
        }

        executeWithExistedElimination(b);
    }

    void executeWithExistedElimination(double* b)
    {
        permuteRHS(b);
        forwardSolve(_mat, b, _n);   // L z = b, b will be overwritten by z
        backwardSolve(_mat, b, _n);  // U x = z, b (z) will be overwritten by x
    }

private:
    /**
     * permute the right hand side vector according to the
     * row permutations of the LU factorization
     * @param b the entries of the vector b are permuted
     */
    void permuteRHS(double* b) const
    {
        for (size_t i = 0; i < _n; i++)
        {
            if (_perm[i] != i)
            {
                BASELIB::swap(b[i], b[_perm[i]]);
            }
        }
    }

    /**
     * a reference to the matrix
     */
    Matrix& _mat;

    /**
     * the size of the matrix
     */
    const size_t _n;

    /**
     * the permutation of the rows
     */

    std::vector<std::size_t> _perm;
};
}  // end namespace MathLib

#endif /* GAUSSALGORITHM_H_ */
