/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ISNAN_H
#define ISNAN_H

#include <math.h>
#ifdef _MSC_VER
#include <float.h>
#endif

namespace BASELIB
{
template <typename T>
inline bool isNAN(T v)
{
#if defined(isnan)
	return isnan(v);
#elif defined(_MSC_VER)
	return _isnan(v) != 0;
#else
	// Remark: the following line doesn't work if compiled with --fast-math flag
	return (v != v);
#endif
}
}

#endif // ISNAN_H
