/*
 * BruteForceClosestPair.h
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef BRUTEFORCECLOSESTPAIR_H_
#define BRUTEFORCECLOSESTPAIR_H_

#include "ClosestPair.h"

namespace GEOLIB
{
class BruteForceClosestPair : public ClosestPair
{
public:
    BruteForceClosestPair(std::vector<GEOLIB::Point*> const& pnts, size_t& id0,
                          size_t& id1);
};
}  // end namespace GEOLIB

#endif /* BRUTEFORCECLOSESTPAIR_H_ */
