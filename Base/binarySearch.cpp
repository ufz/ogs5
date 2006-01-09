/*
 * binarySearch.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: TF
 *
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "binarySearch.h"

size_t searchElement (double const& val, size_t beg, size_t end, const std::vector<double>& array)
{
	if (beg >= end)
		return std::numeric_limits<size_t>::max();
	size_t m ((end + beg) / 2);

	if (array[m] - val <= 0 && array[m + 1] - val > 0)
		return m;
	if (val < array[m])
		return searchElement (val, beg, m, array);
	return searchElement (val, m + 1, end, array);
}

size_t searchElement (double const& val, size_t beg, size_t end, const std::vector<double*>& array)
{
	if (beg >= end)
		return std::numeric_limits<size_t>::max();
	size_t m ((end + beg) / 2);

	if (*(array[m]) - val <= 0 && *(array[m + 1]) - val > 0)
		return m;
	if (val < *(array[m]))
		return searchElement (val, beg, m, array);
	return searchElement (val, m + 1, end, array);
}
