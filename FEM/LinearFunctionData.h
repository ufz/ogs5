/*
 * LinearFunctionData.h
 *
 *  Created on: Sep 1, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef LINEARFUNCTIONDATA_H_
#define LINEARFUNCTIONDATA_H_

#include <fstream>

/*!
   \class LinearFunctionData
   \brief Define a linear function for IC, BC and ST

   WW 24.08.2011
 */
class LinearFunctionData
{
public:
	/*!
	   \brief Constrcutor of class LinearFunctionData

	   \param ifstream &ins: file
	   \param num_var: number of data
	   \param sub_domain: if sub_domain

	   WW 24.08.2011
	 */
	LinearFunctionData(std::ifstream& ins, int num_var = -1);
	~LinearFunctionData();

	double getValue(size_t dom_i, double x, double y, double z) const;
	double getValue(double x, double y, double z) const;
	size_t* getSubDomIndex() const { return _subdom_index; }
private:
	size_t _ndata;
	size_t* _subdom_index;
	// Coefficents for linear distribution function
	// f = a0+b0*x+c0*y+d0*z
	double *_a0, *_b0, *_c0, *_d0;
};

#endif /* LINEARFUNCTIONDATA_H_ */
