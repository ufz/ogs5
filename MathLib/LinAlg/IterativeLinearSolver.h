/*
 * IterativeLinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ITERATIVELINEARSOLVER_H_
#define ITERATIVELINEARSOLVER_H_

#include <LinearSolver.h>

namespace MathLib
{
class IterativeLinearSolver : public MathLib::LinearSolver
{
public:
	IterativeLinearSolver() {}
	virtual ~IterativeLinearSolver() {}
};
}

#endif /* ITERATIVELINEARSOLVER_H_ */
