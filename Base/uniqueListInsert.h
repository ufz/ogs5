/*
 * uniqueListInsert.h
 *
 *  Created on: Feb 23, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef UNIQUELISTINSERT_H_
#define UNIQUELISTINSERT_H_

#include <list>

namespace BASELIB
{
void uniqueListInsert(std::list<size_t>& list, size_t element)
{
    // search element
    std::list<size_t>::const_iterator it;
    for (it = list.begin(); it != list.end(); it++)
        if (*it == element)
            return;
    // element not found -> insert
    list.push_back(element);
}
}  // end namespace BASELIB

#endif /* UNIQUELISTINSERT_H_ */
