/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//**************************************************************************
/*!
   \file msh_lib_ext.cpp

   Define the member function of class CFEMesh for reading
   subdomain mesh

   02/2012 WW/
   last modified
*/
//**************************************************************************
#include "msh_lib.h"

//#ifdef USET_PETSC
// kg44 included for memcpy needed by petsc
#include <stdio.h>
#include <string.h>
#include "petscksp.h"
//#elif USE_MPI
#include <mpi.h>
//#endif
#include <sstream>

#include "StringTools.h"

#if (PETSC_VERSION_NUMBER > 3030)
#include "petsctime.h"
#endif

using namespace std;
using namespace MeshLib;

void BuildNodeStruc(MeshNodes* anode, MPI_Datatype* MPI_Node_ptr);

void FEMRead(const string& file_base_name, vector<MeshLib::CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj, string* unique_name)
{
    int msize;
    int mrank;

    MPI_Comm_size(MPI_COMM_WORLD, &msize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mrank);

    PetscLogDouble v1, v2;
// #define PETSC34
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
    PetscTime(&v1);
#else
    PetscGetTime(&v1);
#endif

    //
    string s_msize = number2str(msize);

    string fname_new_base;
    MPI_File fh;
    int rc = 0;

    // Always try binary file first
    fname_new_base = file_base_name + "_partitioned_msh_cfg" + s_msize + ".bin";

    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc)
    {
        if (mrank == 0)
            printf("-->Reading ASCII mesh file ...");

        FEMRead_ASCII(msize, mrank, file_base_name, mesh_vec, geo_obj,
                      unique_name);
        // return;
    }
    else
    {
        MPI_File_close(&fh);
        if (mrank == 0)
            printf("-->Reading binary mesh file ...");

        FEMRead_BIN(msize, mrank, file_base_name, mesh_vec, geo_obj,
                    unique_name);
    }

// #define PETSC34
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER > 3030)
    PetscTime(&v2);
#else
    PetscGetTime(&v2);
#endif

    if (mrank == 0)
        printf("\t\n>>Total elapsed time in reading mesh:%f s\n", v2 - v1);

    MPI_Barrier(MPI_COMM_WORLD);

    // PetscFinalize();
    // exit(0);
}

void FEMRead_BIN(const int msize, const int mrank, const string& file_base_name,
                 vector<MeshLib::CFEMesh*>& mesh_vec,
                 GEOLIB::GEOObjects* geo_obj, string* unique_name)
{
    // 0 long size_sbd_nodes = 0;
    // 1 long size_sbd_nodes_l = 0;
    // 2 ong size_sbd_nodes_h = 0;
    // 3 long size_sbd_elems = 0;
    // 4 long size_g_elems = 0;
    const MyInt nheaders = 13;
    MyInt mesh_header[nheaders];

    MeshNodes* s_nodes = 0;
    MyInt* elem_info = 0;

    //
    string s_msize;
    stringstream ss;
    ss << msize;
    ss >> s_msize;
    ss.clear();

    MPI_File fh;
    string ftype = "native";
    int rc = 0;
    MPI_Offset offset_new;

    CFEMesh* mesh = new CFEMesh(geo_obj, unique_name);
    mesh_vec.push_back(mesh);

    //--------------------------------------------------------------------------
    // Header
    string fname_new_base =
        file_base_name + "_partitioned_msh_cfg" + s_msize + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc)
    {
        MPI_Finalize();
        if (mrank == 0)
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit(0);
    }

    //
    offset_new = 0;
    offset_new = mrank * nheaders * sizeof(MyInt);
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, &ftype[0],
                      MPI_INFO_NULL);
    MPI_File_read(fh, mesh_header, nheaders, MPI_LONG,
                  MPI_STATUS_IGNORE);  //_all
    MPI_File_close(&fh);

    //--------------------------------------------------------------------------
    // Node
    MPI_Datatype MPI_node;
    fname_new_base = file_base_name + "_partitioned_msh_nod" + s_msize + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc)
    {
        MPI_Finalize();
        if (mrank == 0)
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit(0);
    }

    s_nodes = (MeshNodes*)realloc(s_nodes, sizeof(MeshNodes) * mesh_header[0]);
    BuildNodeStruc(s_nodes, &MPI_node);

    offset_new = mesh_header[10];
    MPI_File_set_view(fh, offset_new, MPI_node, MPI_node, &ftype[0],
                      MPI_INFO_NULL);
    MPI_File_read(fh, s_nodes, mesh_header[0], MPI_node, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    mesh->setSubdomainNodes(mesh_header, s_nodes);
    free(s_nodes);
    s_nodes = NULL;
    MPI_Type_free(&MPI_node);

    //--------------------------------------------------------------------------
    // Element
    fname_new_base = file_base_name + "_partitioned_msh_ele" + s_msize + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc)
    {
        MPI_Finalize();
        if (mrank == 0)
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit(0);
    }

    MyInt size_elem_info = mesh_header[2] + mesh_header[8];
    elem_info = (MyInt*)realloc(elem_info, sizeof(MyInt) * size_elem_info);
    offset_new = mesh_header[11];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, &ftype[0],
                      MPI_INFO_NULL);
    MPI_File_read(fh, elem_info, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    mesh->setSubdomainElements(mesh_header, elem_info, true);

    free(elem_info);
    elem_info = NULL;

    //--------------------------------------------------------------------------
    // Ghost element
    fname_new_base =
        file_base_name + "_partitioned_msh_ele_g" + s_msize + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc)
    {
        MPI_Finalize();
        if (mrank == 0)
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit(0);
    }

    size_elem_info = mesh_header[3] + mesh_header[9];
    elem_info = (MyInt*)realloc(elem_info, sizeof(MyInt) * size_elem_info);
    offset_new = mesh_header[12];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, &ftype[0],
                      MPI_INFO_NULL);
    MPI_File_read(fh, elem_info, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    mesh->setSubdomainElements(mesh_header, elem_info, false);

    free(elem_info);
    elem_info = NULL;

    // MPI_Barrier (MPI_COMM_WORLD);
    mesh->ConstructGrid();
    mesh->FillTransformMatrix();
    // mesh->calMaximumConnectedNodes();
}

void FEMRead_ASCII(const int msize, const int mrank,
                   const string& file_base_name,
                   vector<MeshLib::CFEMesh*>& mesh_vec,
                   GEOLIB::GEOObjects* geo_obj, string* unique_name)
{
    // 0 long size_sbd_nodes = 0;
    // 1 long size_sbd_nodes_l = 0;
    // 2 ong size_sbd_nodes_h = 0;
    // 3 long size_sbd_elems = 0;
    // 4 long size_g_elems = 0;
    const int nheaders = 10;
    long mesh_header[nheaders];

    MeshNodes* s_nodes = (MeshNodes*)malloc(sizeof(MeshNodes));
    long* elem_info = (long*)malloc(1);

    string str_var = "";

    MPI_Datatype MPI_node;
    int tag[] = {0, 1, 2};
    // int ierr;
    // MPI_Request send_request, recv_request;
    MPI_Status status;

    //
    string s_msize;
    stringstream ss;
    ss << msize;
    ss >> s_msize;
    ss.clear();

    ifstream is_cfg;
    ifstream is_node;
    ifstream is_elem;
#ifdef MULTI_MESH_FILE
    stringstream ss(stringstream::in | stringstream::out);
#endif

    if (mrank == 0)
    {
        str_var = file_base_name + "_partitioned_cfg" + s_msize +
                  ".msh";  // "_partitioned.msh";
        is_cfg.open(str_var.c_str());
        if (!is_cfg.good())
        {
            std::cout << "Cannot open file " << str_var << " Quit now! "
                      << std::endl;
            exit(1);
        }

        str_var = file_base_name + "_partitioned_nodes_" + s_msize +
                  ".msh";  // "_partitioned.msh";
        is_node.open(str_var.c_str());
        if (!is_node.good())
        {
            std::cout << "Cannot open file " << str_var << " Quit now! "
                      << std::endl;
            exit(1);
        }

        str_var = file_base_name + "_partitioned_elems_" + s_msize +
                  ".msh";  // "_partitioned.msh";
        is_elem.open(str_var.c_str());
        if (!is_elem.good())
        {
            std::cout << "Cannot open file " << str_var << " Quit now! "
                      << std::endl;
            exit(1);
        }

        getline(is_cfg, str_var);
        int num_parts;
        is_cfg >> num_parts >> ws;
        if (num_parts != msize)
        {
            string str_m =
                "Sorry, I have to quit the simulation now because that "
                " the number of the requested computer cores "
                "is not identical to the number of subdomains.";
            cout << str_m << endl;
            MPI_Finalize();
            exit(1);
        }
    }

    CFEMesh* mesh = new CFEMesh(geo_obj, unique_name);
    mesh_vec.push_back(mesh);

    for (int i = 0; i < msize; i++)
    {
        if (mrank == 0)
        {
            cout << "-->Parallel reading the partitioned mesh: " << i << endl;

            for (int j = 0; j < nheaders; j++)
                is_cfg >> mesh_header[j];
            is_cfg >> ws;
        }
        MPI_Bcast(mesh_header, nheaders, MPI_LONG, 0, MPI_COMM_WORLD);

        // cout<<"\ncccccccccc "<<mesh_header[0]<<"    "<<mesh_header[1]
        //	<<"   " <<mesh_header[2]<<"  "<<mesh_header[3]<<endl;

        //-------------------------------------------------------------------------
        // Node
        s_nodes =
            (MeshNodes*)realloc(s_nodes, sizeof(MeshNodes) * mesh_header[0]);

        if (i > 0)
            BuildNodeStruc(s_nodes, &MPI_node);

        if (mrank == 0)
        {
            // Nodes
            for (long k = 0; k < mesh_header[0]; k++)
            {
                MeshNodes* anode = &s_nodes[k];
                is_node >> anode->index;

                is_node >> anode->x >> anode->y >> anode->z >> ws;
            }

            if (i == 0)
            {
                mesh->setSubdomainNodes(&mesh_header[0], s_nodes);
            }
            else
            {
                MPI_Send(s_nodes, mesh_header[0], MPI_node, i, tag[0],
                         MPI_COMM_WORLD);
            }
        }

        if (i > 0)
        {
            if (mrank == i)
            {
                MPI_Recv(s_nodes, mesh_header[0], MPI_node, 0, tag[0],
                         MPI_COMM_WORLD, &status);
                mesh->setSubdomainNodes(mesh_header, s_nodes);
            }
        }
        free(s_nodes);
        s_nodes = NULL;

        //-------------------------------------------------------------------------
        // Element
        const long size_elem_info = mesh_header[2] + mesh_header[8];
        elem_info = (long*)realloc(elem_info, sizeof(long) * size_elem_info);
        if (mrank == 0)
        {
            long counter = mesh_header[2];
            for (long j = 0; j < mesh_header[2]; j++)
            {
                elem_info[j] = counter;
                is_elem >> elem_info[counter];  // mat. idx
                counter++;
                is_elem >> elem_info[counter];  // type
                counter++;
                is_elem >> elem_info[counter];  // nnodes
                const int nn_e = elem_info[counter];
                counter++;
                for (int k = 0; k < nn_e; k++)
                {
                    is_elem >> elem_info[counter];
                    counter++;
                }
            }

            if (i == 0)
            {
                mesh->setSubdomainElements(mesh_header, elem_info, true);
            }
            else
            {
                MPI_Send(elem_info, size_elem_info, MPI_LONG, i, tag[1],
                         MPI_COMM_WORLD);
            }
        }
        if (i > 0)
        {
            if (mrank == i)
            {
                MPI_Recv(elem_info, size_elem_info, MPI_LONG, 0, tag[1],
                         MPI_COMM_WORLD, &status);
                mesh->setSubdomainElements(mesh_header, elem_info, true);
            }
        }

        // if(elem_info)
        //  {
        free(elem_info);
        elem_info = NULL;
        //  }

        //-------------------------------------------------------------------------
        // Ghost element
        const long size_elem_g_info = mesh_header[3] + mesh_header[9];
        elem_info = (long*)realloc(elem_info, sizeof(long) * size_elem_g_info);
        if (mrank == 0)
        {
            long counter = mesh_header[3];
            for (long j = 0; j < mesh_header[3]; j++)
            {
                elem_info[j] = counter;
                is_elem >> elem_info[counter];  // mat. idx
                counter++;
                is_elem >> elem_info[counter];  // type
                counter++;
                is_elem >> elem_info[counter];  // nnodes
                const int nn_e = elem_info[counter];
                counter++;
                for (int k = 0; k < nn_e; k++)
                {
                    is_elem >> elem_info[counter];
                    counter++;
                }
                is_elem >> elem_info[counter];
                //           const int nn_e_g =  elem_info[counter];
                counter++;
                // ghost nodes for linear element
                is_elem >> elem_info[counter];
                const int nn_e_g = elem_info[counter];
                counter++;
                for (int k = 0; k < nn_e_g; k++)
                {
                    is_elem >> elem_info[counter];
                    counter++;
                }
            }

            if (i == 0)
            {
                mesh->setSubdomainElements(mesh_header, elem_info, false);
            }
            else
            {
                MPI_Send(elem_info, size_elem_g_info, MPI_LONG, i, tag[2],
                         MPI_COMM_WORLD);
            }
        }
        if (i > 0)
        {
            if (mrank == i)
            {
                MPI_Recv(elem_info, size_elem_g_info, MPI_LONG, 0, tag[2],
                         MPI_COMM_WORLD, &status);
                mesh->setSubdomainElements(mesh_header, elem_info, false);
            }
        }

        free(elem_info);
        elem_info = NULL;
    }

    if (s_nodes)
    {
        free(s_nodes);
        s_nodes = NULL;
    }
    if (elem_info)
    {
        free(elem_info);
        elem_info = NULL;
    }

    if (mrank == 0)
    {
        is_cfg.close();
        is_node.close();
        is_elem.close();
    }

    MPI_Type_free(&MPI_node);

    MPI_Barrier(MPI_COMM_WORLD);
    mesh->ConstructGrid();
    mesh->FillTransformMatrix();
    // mesh->calMaximumConnectedNodes();
    MPI_Barrier(MPI_COMM_WORLD);

    // PetscFinalize();
    // exit(0);
}

namespace MeshLib
{
/*!
   Fill data for subdomain mesh

   @param header  : mesh header
   @param s_nodes : mesh nodes

   WW. 02~03.2012
*/
void CFEMesh::setSubdomainNodes(MyInt* header, const MeshNodes* s_nodes)
{
    int k;

    NodesNumber_Quadratic = header[0];
    NodesNumber_Linear = header[1];

    loc_NodesNumber_Linear = header[4];
    loc_NodesNumber_Quadratic = header[5];
    glb_NodesNumber_Linear = header[6];
    glb_NodesNumber_Quadratic = header[7];

    for (k = 0; k < header[0]; k++)
    {
        const MeshNodes* anode = &s_nodes[k];
        CNode* new_node = new CNode(k, anode->x, anode->y, anode->z);

        new_node->SetEquationIndex(anode->index);

        nod_vector.push_back(new_node);
    }
}

/*!
   Fill data for subdomain mesh

   @param header    : mesh header
   @param elem_info : element information
   @param inside :    indicator for inside domain
   WW. 02~03.2012
*/
void CFEMesh::setSubdomainElements(MyInt* header, const MyInt* elem_info,
                                   const bool inside)
{
    MyInt i;
    MyInt k;
    MyInt counter;
    int mat_idx;
    int e_type;
    int nnodes;

    MyInt ne = 0;
    if (inside)
        ne = header[2];
    else
        ne = header[3];

    // Element
    for (i = 0; i < ne; i++)
    {
        CElem* new_elem = new CElem(ele_vector.size());
        ele_vector.push_back(new_elem);

        counter = elem_info[i];

        mat_idx = elem_info[counter];
        counter++;
        e_type = elem_info[counter];
        counter++;
        nnodes = elem_info[counter];
        counter++;

        new_elem->nnodesHQ = nnodes;
        new_elem->nodes_index.resize(new_elem->nnodesHQ);
        for (k = 0; k < new_elem->nnodesHQ; k++)
        {
            new_elem->nodes_index[k] = elem_info[counter];
            counter++;
        }

        if (!inside)
        {
            const int nn_gl = elem_info[counter];
            counter++;
            const int nn_g = elem_info[counter];
            counter++;

            new_elem->g_index = new int[nn_g + 2];
            int* ele_gnidx = new_elem->g_index;
            ele_gnidx[0] = nn_gl;
            ele_gnidx[1] = nn_g;
            for (k = 2; k < nn_g + 2; k++)
            {
                ele_gnidx[k] = elem_info[counter];
                counter++;
            }
        }

        new_elem->patch_index = static_cast<size_t>(mat_idx);

        // element type
        switch (e_type)
        {
            case 1:
                new_elem->geo_type = MshElemType::LINE;
                new_elem->nnodes = 2;
                new_elem->ele_dim = 1;
                new_elem->nfaces = 2;
                new_elem->nedges = 0;
                break;
            case 2:
                new_elem->geo_type = MshElemType::QUAD;
                new_elem->nnodes = 4;
                new_elem->ele_dim = 2;
                new_elem->nfaces = 4;
                new_elem->nedges = 4;

                break;
            case 3:
                new_elem->geo_type = MshElemType::HEXAHEDRON;
                new_elem->nnodes = 8;
                new_elem->ele_dim = 3;
                new_elem->nfaces = 6;
                new_elem->nedges = 12;
                break;
            case 4:
                new_elem->geo_type = MshElemType::TRIANGLE;
                new_elem->nnodes = 3;
                new_elem->ele_dim = 2;
                new_elem->nfaces = 3;
                new_elem->nedges = 3;
                break;
            case 5:
                new_elem->geo_type = MshElemType::TETRAHEDRON;
                new_elem->nnodes = 4;
                new_elem->ele_dim = 3;
                new_elem->nfaces = 4;
                new_elem->nedges = 6;
                break;
            case 6:
                new_elem->geo_type = MshElemType::PRISM;
                new_elem->nnodes = 6;
                new_elem->ele_dim = 3;
                new_elem->nfaces = 5;
                new_elem->nedges = 9;
                break;
            case 7:  //:PYRAMID:
                new_elem->geo_type = MshElemType::PYRAMID;
                new_elem->nnodes = 5;
                new_elem->ele_dim = 3;
                new_elem->nfaces = 5;
                new_elem->nedges = 8;
                break;
            default:
                break;
        }

        new_elem->InitializeMembers();

        // new_elem->WriteIndex();
    }
}

/*!
    Configure the high order elements for parallel computing

    03.2012. WW
*/
void CFEMesh::ConfigHighOrderElements()
{
    int k;
    size_t i, kk;

    CEdge* edge = NULL;
    CElem* elem = NULL;

    int egde_line[] = {2};
    int egde_tri[] = {3, 4, 5};
    int egde_quad[] = {4, 5, 6, 7};
    int egde_tet[] = {4, 5, 6, 7, 8, 9};
    int egde_hex[] = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    int egde_pri[] = {6, 7, 8, 9, 10, 11, 12, 13, 14};

    int* middle_node = NULL;

    // TEST
    // int myrank;
    // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    /*
    for(i=0; i<edge_vector.size(); i++ )
    {
      edge_vector[i]->SetMark(false);
    }
    */

    const size_t e_size = ele_vector.size();

    for (i = 0; i < e_size; i++)
    {
        elem = ele_vector[i];

        // Nodes
        elem->nodes.resize(elem->nnodesHQ);
        for (k = 0; k < elem->nnodesHQ; k++)
            elem->nodes[k] = nod_vector[elem->nodes_index[k]];

        switch (elem->geo_type)
        {
            case MshElemType::LINE:
                middle_node = egde_line;
                break;
            case MshElemType::QUAD:
                middle_node = egde_quad;
                break;
            case MshElemType::HEXAHEDRON:
                middle_node = egde_hex;
                break;
            case MshElemType::TRIANGLE:
                middle_node = egde_tri;
                break;
            case MshElemType::TETRAHEDRON:
                middle_node = egde_tet;
                break;
            case MshElemType::PRISM:
                middle_node = egde_pri;
                break;
            case MshElemType::PYRAMID:
                break;
            default:
                break;
        }

        // Edges
        for (kk = 0; kk < elem->nedges; kk++)
        {
            edge = elem->edges[kk];

            //  if(edge->GetMark())
            //  continue;

            edge->SetNode(2, elem->nodes[middle_node[kk]]);

            // TEST
            // if(myrank == 0)
            //   edge->Write();
            //        edge->SetMark(true);
        }
    }
    /*
    for(i=0; i<edge_vector.size(); i++ )
    {
      edge_vector[i]->SetMark(false);
    }
    */

    CNode* node = NULL;
    for (size_t e = 0; e < e_size; e++)
    {
        elem = ele_vector[e];
        for (int i = elem->nnodes; i < elem->nnodesHQ; i++)
        {
            bool done = false;
            node = elem->GetNode(i);
            for (int k = 0; k < (int)node->getConnectedElementIDs().size(); k++)
                if (e == node->getConnectedElementIDs()[k])
                {
                    done = true;
                    break;
                }
            if (!done)
                node->getConnectedElementIDs().push_back(e);
        }
    }

    // For sparse matrix
    ConnectedNodes(true);
    ConnectedElements2Node(true);
}
//
int CFEMesh::calMaximumConnectedNodes()
{
    int max_connected_nodes = 0;
    for (size_t i = 0; i < nod_vector.size(); i++)
    {
        const int k = static_cast<int>(nod_vector[i]->getNumConnectedNodes());
        if (k > max_connected_nodes)
            max_connected_nodes = k;
    }

    int msize;
    // int mrank;
    int* i_cnt;
    int* i_disp;
    int* i_recv;

    MPI_Comm_size(MPI_COMM_WORLD, &msize);
    // MPI_Comm_rank(MPI_COMM_WORLD,&mrank);

    i_cnt = (int*)malloc(msize * sizeof(int));
    i_disp = (int*)malloc(msize * sizeof(int));
    i_recv = (int*)malloc(msize * sizeof(int));

    for (int i = 0; i < msize; i++)
    {
        i_cnt[i] = 1;
        i_disp[i] = i;
    }

    MPI_Allgatherv(&max_connected_nodes, 1, MPI_INT, i_recv, i_cnt, i_disp,
                   MPI_INT, MPI_COMM_WORLD);

    max_connected_nodes = 0;
    for (int i = 0; i < msize; i++)
    {
        if (i_recv[i] > max_connected_nodes)
            max_connected_nodes = i_recv[i];
    }
    free(i_cnt);
    free(i_recv);
    free(i_disp);

    return max_connected_nodes;
}

}  // namespace MeshLib

void BuildNodeStruc(MeshNodes* anode, MPI_Datatype* MPI_Node_ptr)
{
    MPI_Datatype my_comp_type[4];
    int nblocklen[4];
    MPI_Aint disp[4], base;
    int j;

    my_comp_type[0] = MPI_LONG;
    my_comp_type[1] = MPI_DOUBLE;
    my_comp_type[2] = MPI_DOUBLE;
    my_comp_type[3] = MPI_DOUBLE;
    nblocklen[0] = 1;
    nblocklen[1] = 1;
    nblocklen[2] = 1;
    nblocklen[3] = 1;

    // Compute displacement of struct MeshNodes

    /*
     disp[0] = 0;
     MPI_Get_address(&(anode.nnodes),&base);
     MPI_Get_address(&(anode.index), disp+1);
     disp[1] -= base;
     MPI_Get_address(&(anode.coord),disp+2);
     disp[2] -= base;

    for (j=0; j <3; j++)
    {

        cout<<"j=" <<j<<" "<<disp[j]<<endl;
    }
    */

    MPI_Get_address(anode, disp);
    MPI_Get_address(&(anode[0].x), disp + 1);
    MPI_Get_address(&(anode[0].y), disp + 2);
    MPI_Get_address(&(anode[0].z), disp + 3);
    base = disp[0];
    for (j = 0; j < 4; j++)
    {
        disp[j] -= base;

        //      cout<<"j=" <<j<<" "<<disp[j]<<endl;
    }

    // build datatype describing structure
    MPI_Type_create_struct(4, nblocklen, disp, my_comp_type, MPI_Node_ptr);
    MPI_Type_commit(MPI_Node_ptr);
}
