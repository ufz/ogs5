/*
 * fct_mpi.cpp
 *
 *  Created on: Apr 19, 2013
 *      Author: NW
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "fct_mpi.h"

#include <iostream>
#include <fstream>

#include <mpi.h>

#include "readNonBlankLineFromInputStream.h"

namespace FCT_MPI
{
// definition
FCT_MPI::CommunicationTable ct;

std::ios::pos_type ReadRank(std::ifstream* is, int myrank)
{
	std::string line_string;
	bool new_keyword = false;
	std::ios::pos_type position;

	std::string sub_string, strbuff;
	std::stringstream in;

	int temp = 0;
	size_t n_neighbors = 0;

	while (!new_keyword)
	{
		position = is->tellg();
		line_string = readNonBlankLineFromInputStream(*is);
		if (line_string.size() < 1)
			break;
		if (line_string.find("#") != std::string::npos)
		{
			new_keyword = true;
			break;
		}

		if (line_string.find("$MYRANK") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*is));
			in >> temp;
			in.clear();
			if (temp != myrank)
			{
				std::cout << "rank " << myrank << ": $MYRANK is wrong. Stop reading\n";
				break;
			}
		}
		else if (line_string.find("$NNEIGHBORS") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*is));
			in >> n_neighbors;
			in.clear();
		}
		else if (line_string.find("$NEIGHBOR") != std::string::npos)
		{
			Neighbor nei;
			in.str(readNonBlankLineFromInputStream(*is));
			in >> nei.rank;
			in.clear();
			in.str(readNonBlankLineFromInputStream(*is));
			size_t n_edges;
			in >> n_edges;
			in.clear();
			nei.overlapping_edges.resize(n_edges);
			int ie1, ie2;
			for (size_t i_e = 0; i_e < n_edges; i_e++)
			{
				in.str(readNonBlankLineFromInputStream(*is));
				in >> ie1 >> ie2;
				in.clear();
				nei.overlapping_edges[i_e] = Edge(ie1, ie2);
			}
			ct.vec_neighbors.push_back(nei);
		}
	}

	if (n_neighbors != ct.vec_neighbors.size())
	{
		std::cout << "rank " << myrank << ": Inconsistent NNEIGHBORS!\n";
	}

	for (size_t i = 0; i < ct.vec_neighbors.size(); i++)
	{
		std::set<long> overlapping_innder_nodes;
		std::set<long> overlapping_ghost_nodes;
		Neighbor& nei = ct.vec_neighbors[i];
		for (size_t i_e = 0; i_e < nei.overlapping_edges.size(); i_e++)
		{
			overlapping_innder_nodes.insert(nei.overlapping_edges[i_e].first);
			overlapping_ghost_nodes.insert(nei.overlapping_edges[i_e].second);
		}
		nei.overlapping_inner_nodes.assign(overlapping_innder_nodes.begin(), overlapping_innder_nodes.end());
		nei.overlapping_ghost_nodes.assign(overlapping_ghost_nodes.begin(), overlapping_ghost_nodes.end());
	}

	return position;
}

void FCTCommRead(const std::string& file_base_name)
{
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	std::string str_var = file_base_name + ".cct";
	std::ifstream bc_file(str_var.data(), std::ios::in);
	if (!bc_file.good())
	{
		return;
	}
	std::cout << "rank " << myrank << ": -> Read CCT ...\n";

	char line[MAX_ZEILE];
	std::string line_string;
	int i_table = 0;
	while (!bc_file.eof())
	{
		bc_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos)
		{
			std::cout << "done \n";
			break;
		}
		if (line_string.find("#COMMUNICATION_TABLE") != std::string::npos)
		{
			if (i_table == myrank)
			{
				std::ios::pos_type position = ReadRank(&bc_file, myrank);
				bc_file.seekg(position, std::ios::beg);
				break;
			}
			i_table++;
		} // keyword found
	} // eof

	bc_file.close();

	if (ct.vec_neighbors.size() > 0)
	{
		std::cout << "rank " << myrank << ": -> CCT data loaded\n";
	}
	else
	{
		std::cout << "rank " << myrank << ": -> No CCT data loaded\n";
	}
}

void gatherR(const CommunicationTable& ct, Math_Group::Vec& globalR_plus, Math_Group::Vec& globalR_min)
{
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// send R^+,R^- and get R^+,R^- to/from neighbors (i: inner nodes, j: ghost nodes)
	int n_neighbors = ct.vec_neighbors.size();
	std::vector<MPI_Status> vec_status(n_neighbors);
	std::vector<MPI_Request> vec_req(n_neighbors);
	std::vector<std::vector<double> > vec_sendbuf(n_neighbors);
	std::vector<std::vector<double> > vec_recvbuf(n_neighbors);
	// prepare send and receive buffer
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		vec_sendbuf[i].resize(neighbor.overlapping_inner_nodes.size() * 2);
		vec_recvbuf[i].resize(neighbor.overlapping_ghost_nodes.size() * 2);
		for (size_t j = 0; j < neighbor.overlapping_inner_nodes.size(); j++)
		{
			const long inner_nodeID = neighbor.overlapping_inner_nodes[j];
			vec_sendbuf[i][2 * j] = globalR_plus(inner_nodeID);
			vec_sendbuf[i][2 * j + 1] = globalR_min(inner_nodeID);
		}
	}
	// send
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		int neighbor_rank = neighbor.rank;
		MPI_Isend(&vec_sendbuf[i][0], vec_sendbuf[i].size(), MPI_DOUBLE, neighbor_rank, 0, MPI_COMM_WORLD, &vec_req[i]);
	}

	// receive
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		int neighbor_rank = neighbor.rank;
		MPI_Irecv(&vec_recvbuf[i][0], vec_recvbuf[i].size(), MPI_DOUBLE, neighbor_rank, 0, MPI_COMM_WORLD, &vec_req[i]);
	}

	// wait
	MPI_Waitall(n_neighbors, &vec_req[0], &vec_status[0]);

	// update R^+, R^-
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		size_t n_ghosts = neighbor.overlapping_ghost_nodes.size();
		for (size_t j = 0; j < n_ghosts; j++)
		{
			const long nodeID = neighbor.overlapping_ghost_nodes[j];
			globalR_plus(nodeID) = vec_recvbuf[i][2 * j];
			globalR_min(nodeID) = vec_recvbuf[i][2 * j + 1];
		}
	}
}

void gatherK(const CommunicationTable& ct, Math_Group::SparseMatrixDOK& globalK)
{
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// assume K is already assembled from local elements
	if (false)
	{
		int n_ranks;
		MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
		MPI_Status status;
		int done = 0;
		if (myrank > 0)
		{
			MPI_Recv(&done, 1, MPI_INT, myrank - 1, 0, MPI_COMM_WORLD, &status);
		}

		std::stringstream temp_in;
		globalK.Write(temp_in);

		if (myrank < n_ranks - 1)
		{
			MPI_Send(&done, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD);
		}
	}
	// send K_ij and get K_ji to/from neighbors (i: inner nodes, j: ghost nodes)
	int n_neighbors = ct.vec_neighbors.size();
	std::vector<MPI_Status> vec_status(n_neighbors);
	std::vector<MPI_Request> vec_req(n_neighbors);
	std::vector<std::vector<double> > vec_sendbuf(n_neighbors);
	std::vector<std::vector<double> > vec_recvbuf(n_neighbors);
	// prepare send and receive buffer
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		vec_sendbuf[i].resize(neighbor.overlapping_edges.size());
		vec_recvbuf[i].resize(neighbor.overlapping_edges.size());
		//		std::cout << "rank= " << myrank << ": neighbor rank " << neighbor.rank << "\n";
		for (size_t j = 0; j < neighbor.overlapping_edges.size(); j++)
		{
			const Edge& edge = neighbor.overlapping_edges[j];
			vec_sendbuf[i][j] = globalK(edge.first, edge.second);
			//			std::cout << "rank= " << myrank << ": neighbor rank " << neighbor.rank << "- edge (" <<
			// edge.first
			//<< "," << edge.second << ")\n";
		}
	}
	// send
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		int neighbor_rank = neighbor.rank;
		MPI_Isend(&vec_sendbuf[i][0], vec_sendbuf[i].size(), MPI_DOUBLE, neighbor_rank, 0, MPI_COMM_WORLD, &vec_req[i]);
		//		std::cout << "rank= " << myrank << ": sending " << vec_sendbuf[i].size() << " data to rank " <<
		// neighbor_rank << "\n";
	}

	// receive
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		int neighbor_rank = neighbor.rank;
		MPI_Irecv(&vec_recvbuf[i][0], vec_recvbuf[i].size(), MPI_DOUBLE, neighbor_rank, 0, MPI_COMM_WORLD, &vec_req[i]);
		//		std::cout << "rank= " << myrank << ": receiving " << vec_recvbuf[i].size() << " data from rank " <<
		// neighbor_rank << "\n";
	}

	// wait
	MPI_Waitall(n_neighbors, &vec_req[0], &vec_status[0]);
	//	std::cout << "rank= " << myrank << ": MPI_Waitall() done\n";

	// update K_ji
	for (int i = 0; i < n_neighbors; i++)
	{
		const Neighbor& neighbor = ct.vec_neighbors[i];
		size_t n_edges = neighbor.overlapping_edges.size();
		for (size_t j = 0; j < n_edges; j++)
		{
			const Edge& edge = neighbor.overlapping_edges[j];
			globalK(edge.second, edge.first) += vec_recvbuf[i][j];
		}
	}

	//	if (false) {
	//		int n_ranks;
	//		MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
	//		MPI_Status status;
	//		int done = 0;
	//		if (myrank>0) {
	//			MPI_Recv(&done, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, &status);
	//		}
	//
	//		for (int i=0; i<vec_recvbuf.size(); i++) {
	//			std::cout << "rank= " << myrank << ": neighbor rank " << i << " - ";
	//			for (size_t j=0; j<vec_recvbuf[i].size(); j++) {
	//				std::cout << vec_recvbuf[i][j] << " ";
	//			}
	//			std::cout << "\n";
	//		}
	//
	//		std::stringstream temp_in;
	//		globalK.Write(temp_in);
	//		std::cout << "rank= " << myrank << ": K is updated. K =\n" << temp_in.str();
	//
	//		if (myrank<n_ranks-1) {
	//			MPI_Send(&done, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD);
	//		}
	//		MPI_Barrier (MPI_COMM_WORLD);
	//	}
}

void computeD(const MeshLib::CFEMesh* m_msh, const Math_Group::SparseMatrixDOK& K, Math_Group::SparseMatrixDOK& d)
{
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	const Math_Group::SparseMatrixDOK::mat_t& fct_f = K.GetRawData();
	const Math_Group::SparseMatrixDOK::col_t* col;
	Math_Group::SparseMatrixDOK::mat_t::const_iterator ii;
	Math_Group::SparseMatrixDOK::col_t::const_iterator jj;

	const size_t node_size = m_msh->NodesInUsage();
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_glob = FCT_GLOB_ADDRESS(i);
		col = &fct_f[i_glob];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			const size_t j_glob = j; // FCT_GLOB_ADDRESS(j);
			if (i_glob > j_glob || i_glob == j_glob)
				continue; // symmetric part, off-diagonal

			// Get artificial diffusion operator D
			double K_ij = K(i_glob, j_glob);
			double K_ji = K(j_glob, i_glob);
			if (K_ij == 0.0 && K_ji == 0.0)
				continue;
			double d1 = GetFCTADiff(K_ij, K_ji);
			// double d0 = d1; //TODO should use AuxMatrix at the previous time step
			// if (list_bc_nodes.find(i)!=list_bc_nodes.end() || list_bc_nodes.find(j)!=list_bc_nodes.end()) {
			//  d1 = d0 = 0.0;
			//}
			d(i_glob, j_glob) += d1;
			d(j_glob, i_glob) += d1;
			d(i_glob, i_glob) += -d1;
			d(j_glob, j_glob) += -d1;
		}
	}
}

void debugADFlux(int myrank, MeshLib::CFEMesh* m_msh, size_t node_size, Math_Group::SparseMatrixDOK::mat_t& fct_f)
{
	using Math_Group::SparseMatrixDOK;

	SparseMatrixDOK::col_t* col;
	SparseMatrixDOK::mat_t::const_iterator ii;
	SparseMatrixDOK::col_t::const_iterator jj;
	int n_ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
	MPI_Status status;
	int done = 0;
	if (myrank > 0)
	{
		MPI_Recv(&done, 1, MPI_INT, myrank - 1, 0, MPI_COMM_WORLD, &status);
	}

	//        for (size_t i = 0; i < node_size; i++)
	//        {
	//            const size_t i_global = FCT_GLOB_ADDRESS(i);
	//            std::cout << "rank= " << myrank << ": node i=" << i_global << ", ML=" << (*ML)(i_global) << "\n";
	//        }
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_global = FCT_GLOB_ADDRESS(i);
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			const size_t j_global = FCT_GLOB_ADDRESS(j);
			std::cout << i_global << ", " << j_global << ": " << jj->second << "\n";
		}
	}
	std::cout << std::flush;

	if (myrank < n_ranks - 1)
	{
		MPI_Send(&done, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

} // end FCT_MPI
