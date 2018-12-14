/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 WW/OK Encapsulation from rf_ele_msh
   last modified
**************************************************************************/
#include "msh_core.h"

namespace MeshLib
{
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   03/2010 TF initialization in initialization list
**************************************************************************/
CCore::CCore(size_t id)
    : index(id), boundary_type('I'), mark(true), quadratic(false)
{
}
}  // namespace MeshLib
