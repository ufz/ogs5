/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SWAP_H_
#define SWAP_H_

namespace BASELIB
{
/**
 * swap the content of arg0 and arg1
 */
template <class T>
void swap(T& arg0, T& arg1)
{
	T temp(arg0);
	arg0 = arg1;
	arg1 = temp;
}
} // end namespace BASELIB

#endif // SWAP_H_
