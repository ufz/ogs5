/*
 * binarySearch.h
 *
 *  Created on: Jun 7, 2010
 *      Author: TF
 *
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// STL
#include <cstddef>
#include <limits>
#include <vector>

#ifndef BINARYSEARCH_H_
#define BINARYSEARCH_H_

/**
 * Binary search in a sorted vector of elements to get the
 * id of an element according its key.
 * @param key the key for the element
 * @param beg beginning index in the sorted vector of elements
 * @param end ending index in the sorted vector of elements
 * @param array the vector of elements
 * @return the id of the element in the vector or, if not found,
 * the value std::numeric_limits<size_t>::max()
 */
template <class T>
size_t searchElement(const T& key, size_t beg, size_t end,
                     const std::vector<T>& array)
{
    if (beg >= end)
        return std::numeric_limits<size_t>::max();
    size_t m((end + beg) / 2);

    if (key == array[m])
        return m;
    if (key < array[m])
        return searchElement(key, beg, m, array);
    return searchElement(key, m + 1, end, array);
}

template <class T>
size_t getLargestIndexSmallerThanElement(const T& key, size_t beg, size_t end,
                                         const std::vector<T>& array)
{
    if (beg >= end)
        return std::numeric_limits<size_t>::max();

    size_t m((end + beg) / 2);

    if (m + 1 == end)
    {
        if (key == array[m])
            return m;
        else
            std::numeric_limits<size_t>::max();
    }

    if (key >= array[m] && key < array[m + 1])
        return m;
    if (key < array[m])
        return getLargestIndexSmallerThanElement(key, beg, m, array);
    return getLargestIndexSmallerThanElement(key, m + 1, end, array);
}

size_t searchElement(double const& val, size_t beg, size_t end,
                     const std::vector<double>& array);
size_t searchElement(double const& val, size_t beg, size_t end,
                     const std::vector<double*>& array);

#endif /* BINARYSEARCH_H_ */
