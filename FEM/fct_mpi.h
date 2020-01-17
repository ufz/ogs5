/*
 * fct_mpi.h
 *
 *  Created on: Apr 19, 2013
 *      Author: NW
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FCT_MPI_H_
#define FCT_MPI_H_

#include <vector>
#include <string>

#include "../MSH/msh_mesh.h"

#include "SparseMatrixDOK.h"
#include "mathlib.h"

namespace FCT_MPI
{
#ifndef USE_PETSC

#define FCT_GLOB_ADDRESS(i) i

#else

#define FCT_GLOB_ADDRESS(i) m_msh->getNodes()[i]->GetEquationIndex()

typedef std::pair<long, long> Edge;  // internal node <-> ghost node

struct Neighbor
{
    int rank;
    std::vector<Edge> overlapping_edges;
    std::vector<long> overlapping_inner_nodes;
    std::vector<long> overlapping_ghost_nodes;

    Neighbor() : rank(0){};
};

struct CommunicationTable
{
    std::vector<Neighbor> vec_neighbors;
};

extern FCT_MPI::CommunicationTable ct;

void FCTCommRead(const std::string& file_base_name);

void gatherK(const CommunicationTable& ct,
             Math_Group::SparseMatrixDOK& globalK);

void gatherR(const CommunicationTable& ct, Math_Group::Vec& globalR_plus,
             Math_Group::Vec& globalR_min);

void computeD(const MeshLib::CFEMesh* m_msh,
              const Math_Group::SparseMatrixDOK& K,
              Math_Group::SparseMatrixDOK& d);

void debugADFlux(int myrank, MeshLib::CFEMesh* m_msh, size_t node_size,
                 Math_Group::SparseMatrixDOK::mat_t& fct_f);

#endif

}  // namespace FCT_MPI

#endif /* FCT_MPI_H_ */
