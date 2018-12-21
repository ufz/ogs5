/*
 * printList.h
 *
 *  Created on: Feb 23, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PRINTLIST_H_
#define PRINTLIST_H_

// STL
#include <iostream>
#include <list>
#include <string>

namespace BASELIB
{
void printList(std::list<size_t> const& mylist, std::string const& title)
{
    std::cout << title << "\n";
    for (std::list<size_t>::const_iterator my_it(mylist.begin());
         my_it != mylist.end();
         my_it++)
        std::cout << *my_it << " ";
    std::cout << "\n";
}
}  // end namespace BASELIB

#endif /* PRINTLIST_H_ */
