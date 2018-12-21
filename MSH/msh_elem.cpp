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

#include "mathlib.h"
#include <cmath>
#include <cstdlib>  //WW
#include <float.h>  //WW
// MSHLib
// WW#include "MSHEnums.h" // KR 2010/11/15
#include "msh_elem.h"

namespace MeshLib
{
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
CElem::CElem(size_t Index)
    : CCore(Index),
      normal_vector(NULL),
      /*geo_type(t), */ owner(NULL),
      ele_dim(1),
      nnodes(0),
      nnodesHQ(0),
      nodes(nnodes),
      nodes_index(nnodes),
      patch_index(0),
      transform_tensor(NULL)
{
    grid_adaptation = -1;
    volume = 0.0;
    face_index = -1;
    no_faces_on_surface = 0;
    gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
    normal_vector = NULL;
    area = 1.0;      // WW
    excavated = -1;  // WX

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
CElem::CElem() : CCore(0), normal_vector(NULL)
{
    selected = 0;
    // matgroup_view = 0;
    grid_adaptation = -1;
    nnodes = 0;
    nnodesHQ = 0;
    ele_dim = 1;  // Dimension of element
    patch_index = 0;
    //
    volume = 0.0;
    face_index = -1;
    no_faces_on_surface = 0;
    owner = NULL;
    nodes.resize(8);  // Nodes of face
    transform_tensor = NULL;
    gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
    area = 1.0;  // WW area = 1.0
    normal_vector = NULL;
    excavated = -1;  // WX

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
CElem::CElem(size_t Index, CElem* onwer, int Face)
    : CCore(Index), normal_vector(NULL), owner(onwer)
{
    int i, j, k, n, ne;
    int faceIndex_loc[10];
    int edgeIndex_loc[10] = {};
    no_faces_on_surface = 0;
    n = owner->GetElementFaceNodes(Face, faceIndex_loc);
    face_index = Face;
    gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
    transform_tensor = NULL;
    normal_vector = NULL;
    excavated = -1;  // WX
    //
    switch (owner->geo_type)
    {
        // case MshElemType::LINE:  // 1-D bar element //KR need not be
        // processed
        case MshElemType::QUAD:  // 2-D quadrilateral element
            this->setElementProperties(MshElemType::LINE);
            break;
        case MshElemType::HEXAHEDRON:  // 3-D hexahedral element
            this->setElementProperties(MshElemType::QUAD);
            break;
        case MshElemType::TRIANGLE:  // 2-D triagular element
            this->setElementProperties(MshElemType::LINE);
            break;
        case MshElemType::TETRAHEDRON:  // 3-D tetrahedral element
            this->setElementProperties(MshElemType::TRIANGLE);
            break;
        case MshElemType::PRISM:  // 3-D prismatic element
            if (Face < 2)         // top or bottom face of the prism
                this->setElementProperties(MshElemType::TRIANGLE);
            else  // side of the prism
                this->setElementProperties(MshElemType::QUAD);
            break;
        case MshElemType::PYRAMID:  // 3-D pyramid element
            if (Face < 1)           // bottom face
                this->setElementProperties(MshElemType::QUAD);
            else  // side faces
                this->setElementProperties(MshElemType::TRIANGLE);
            break;
        default:
            std::cerr << "CElem::CElem MshElemType not handled"
                      << "\n";
    }

    patch_index = owner->patch_index;
    quadratic = owner->quadratic;
    nodes_index.resize(n);
    nodes.resize(n);

    boundary_type = 'B';
    for (i = 0; i < n; i++)
    {
        nodes_index[i] = owner->nodes_index[faceIndex_loc[i]];
        nodes[i] = owner->nodes[faceIndex_loc[i]];
        // 18.02.2009. cf. changes in mapping & generation. WW
        if ((nodes[i]->GetBoundaryType() != '0') &&
            (nodes[i]->GetBoundaryType() != '1'))
            nodes[i]->SetBoundaryType('B');
    }
    // Face edges
    ne = owner->GetEdgesNumber();
    edges.resize(nnodes);
    edges_orientation.resize(nnodes);
    edges_orientation = 1;
    for (i = 0; i < nnodes; i++)
    {
        k = (i + 1) % nnodes;
        for (j = 0; j < ne; j++)
        {
            owner->GetLocalIndicesOfEdgeNodes(j, edgeIndex_loc);
            if ((faceIndex_loc[i] == edgeIndex_loc[0] &&
                 faceIndex_loc[k] == edgeIndex_loc[1]) ||
                (faceIndex_loc[i] == edgeIndex_loc[1] &&
                 faceIndex_loc[k] == edgeIndex_loc[0]))
            {
                edges[i] = owner->edges[j];
                if (faceIndex_loc[i] == edgeIndex_loc[1] &&
                    faceIndex_loc[k] == edgeIndex_loc[0])
                    edges_orientation[i] = -1;
                edges[i]->boundary_type = 'B';
                break;
            }
        }
    }

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW/OK Implementation
**************************************************************************/
CElem::CElem(size_t Index, CElem* m_ele_parent)
    : CCore(Index), normal_vector(NULL)
{
    //  static int faceIndex_loc[10];
    //  static int edgeIndex_loc[10];
    gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
    transform_tensor = NULL;
    this->setElementProperties(m_ele_parent->geo_type);

    patch_index = m_ele_parent->patch_index;
    quadratic = m_ele_parent->quadratic;

    boundary_type = 'I';

    // Initialize topological properties
    neighbors.resize(nfaces);
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i] = NULL;
    edges.resize(nedges);
    edges_orientation.resize(nedges);
    for (size_t i = 0; i < nedges; i++)
    {
        edges[i] = NULL;
        edges_orientation[i] = 1;
    }
    area = 1.0;  // WW

    excavated = -1;  // 12.08.2011. WW

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}

CElem::CElem(CElem const& elem)
    : CCore(elem.GetIndex()),
      mat_vector(elem.mat_vector),
      // matgroup_view (elem.matgroup_view),
      selected(elem.selected),
      normal_vector(new double[3]),
      representative_length(elem.representative_length),
      courant(elem.courant),
      neumann(elem.neumann),
      geo_type(elem.geo_type),
      owner(elem.owner),
      ele_dim(elem.ele_dim),
      nnodes(elem.nnodes),
      nnodesHQ(elem.nnodesHQ),
      nodes(elem.nodes),
      nedges(elem.nedges),
      edges(elem.edges),
      edges_orientation(elem.edges_orientation),
      nfaces(elem.nfaces),
      no_faces_on_surface(elem.no_faces_on_surface),
      face_index(elem.face_index),
      volume(elem.volume),
      grid_adaptation(elem.grid_adaptation),
      patch_index(elem.patch_index),
      area(elem.area),
      angle(elem.angle)
{
    for (size_t k(0); k < 3; k++)
    {
        normal_vector[k] = elem.normal_vector[k];
        gravity_center[k] = elem.gravity_center[k];
    }

    // copy nodes
    nodes.resize((int)((elem.nodes).Size()));
    for (size_t k(0); k < (elem.nodes).Size(); k++)
        nodes[k] =
            new CNode((elem.nodes[k])->GetIndex(), (elem.nodes[k])->getData());

        // copy edges
        //	edges = vec<CEdge*> ((elem.edges).Size());
        //	for (size_t k(0); k<(elem.edges).Size() ; k++) {
        //		edges[k] = new CNode ((elem.nodes[k])->GetIndex(),
        //(elem.nodes[k])->X(), (elem.nodes[k])->Y(), (elem.nodes[k])->Z());
        //	}

        //*angle = *(elem.angle);
        // copy transform tensor
        //	Matrix * transform_tensor;

        // copy neighbors
        //	vec<CElem*> neighbors;

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = elem.g_index;
#endif
}

CElem::CElem(MshElemType::type t, size_t node0, size_t node1, size_t node2,
             int mat)
    : CCore(0),
      normal_vector(NULL),
      geo_type(t),
      owner(NULL),
      ele_dim(2),
      nnodes(3),
      nnodesHQ(6),
      nodes(nnodes),
      nodes_index(nnodes),
      nedges(3),
      edges(nedges),
      edges_orientation(nedges),
      nfaces(3),
      patch_index(mat),
      transform_tensor(NULL),
      neighbors(nfaces)
{
    nodes_index[0] = node0;
    nodes_index[1] = node1;
    nodes_index[2] = node2;

    // Initialize topological properties
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i] = NULL;
    for (size_t i = 0; i < nedges; i++)
    {
        edges[i] = NULL;
        edges_orientation[i] = 1;
    }

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}

CElem::CElem(MshElemType::type t, size_t node0, size_t node1, size_t node2,
             size_t node3, int mat)
    : CCore(0),
      normal_vector(NULL),
      geo_type(t),
      owner(NULL),
      ele_dim(2),
      nnodes(4),
      nnodesHQ(9),
      nodes(nnodes),
      nodes_index(nnodes),
      nedges(4),
      edges(nedges),
      edges_orientation(nedges),
      nfaces(4),
      patch_index(mat),
      transform_tensor(NULL),
      neighbors(nfaces)
{
    nodes_index[0] = node0;
    nodes_index[1] = node1;
    nodes_index[2] = node2;
    nodes_index[3] = node3;

    // Initialize topological properties
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i] = NULL;
    for (size_t i = 0; i < nedges; i++)
    {
        edges[i] = NULL;
        edges_orientation[i] = 1;
    }
#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    g_index = NULL;
#endif
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
CElem::~CElem()
{
    // HACK LB those resize() do not free memory allocated with new
    // HACK LB Possible memory leak otherwise resize() not necessary
    nodes_index.resize(0);
    nodes.resize(0);
    edges.resize(0);
    neighbors.resize(0);
    mat_vector.resize(0);
    edges_orientation.resize(0);
    owner = NULL;
    if (transform_tensor)
        delete transform_tensor;
    transform_tensor = NULL;
    transform_tensor = NULL;
    if (normal_vector)
        delete[] normal_vector;
    normal_vector = NULL;

#if defined(USE_PETSC)  // || defined(using other parallel scheme). WW
    if (g_index)
        delete[] g_index;
    g_index = NULL;
#endif
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::FillTransformMatrix(const bool recompute_matrix)
{
    if (transform_tensor && recompute_matrix == false)
        return;

    double xx[3];
    double yy[3];
    double zz[3];
    transform_tensor = new Math_Group::Matrix(3, 3);
    if (geo_type == MshElemType::LINE)
    {
        // x"_vec
        double const* const pnt0(nodes[0]->getData());
        double const* const pnt1(nodes[1]->getData());
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        xx[2] = pnt1[2] - pnt0[2];
        NormalizeVector(xx, 3);
        // an arbitrary vector
        for (size_t i = 0; i < 3; i++)
            yy[i] = 0.0;
        // WW. 06.11.2007
        if (fabs(xx[0]) > 0.0 && fabs(xx[1]) + fabs(xx[2]) < DBL_MIN)
            yy[2] = 1.0;
        else if (fabs(xx[1]) > 0.0 && fabs(xx[0]) + fabs(xx[2]) < DBL_MIN)
            yy[0] = 1.0;
        else if (fabs(xx[2]) > 0.0 && fabs(xx[0]) + fabs(xx[1]) < DBL_MIN)
            yy[1] = 1.0;
        else
        {
            for (size_t i = 0; i < 3; i++)
                if (fabs(xx[i]) > 0.0)
                {
                    yy[i] = -xx[i];
                    break;
                }
        }
        // z"_vec
        CrossProduction(xx, yy, zz);
        NormalizeVector(zz, 3);
        // y"_vec
        CrossProduction(zz, xx, yy);
        NormalizeVector(yy, 3);
    }
    else if (geo_type == MshElemType::QUAD || geo_type == MshElemType::QUAD8 ||
             geo_type == MshElemType::TRIANGLE)
    {
        // x"_vec
        //			xx[0] = nodes[1]->X() - nodes[0]->X();
        //			xx[1] = nodes[1]->Y() - nodes[0]->Y();
        //			xx[2] = nodes[1]->Z() - nodes[0]->Z();
        double const* const pnt0(nodes[0]->getData());
        double const* const pnt1(nodes[1]->getData());
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        xx[2] = pnt1[2] - pnt0[2];
        NormalizeVector(xx, 3);
        // a vector on the plane
        //			yy[0] = nodes[2]->X() - nodes[1]->X();
        //			yy[1] = nodes[2]->Y() - nodes[1]->Y();
        //			yy[2] = nodes[2]->Z() - nodes[1]->Z();
        double const* const pnt2(nodes[2]->getData());
        yy[0] = pnt2[0] - pnt1[0];
        yy[1] = pnt2[1] - pnt1[1];
        yy[2] = pnt2[2] - pnt1[2];
        // z"_vec. off plane
        CrossProduction(xx, yy, zz);
        NormalizeVector(zz, 3);
        // y"_vec
        CrossProduction(zz, xx, yy);
        NormalizeVector(yy, 3);
    }
    for (size_t i = 0; i < 3; i++)
    {
        (*transform_tensor)(i, 0) = xx[i];
        (*transform_tensor)(i, 1) = yy[i];
        (*transform_tensor)(i, 2) = zz[i];
    }
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   02/2006 PCH Implementation
**************************************************************************/
double CElem::getTransformTensor(int idx)
{
    // WW
    int i = idx % 3;
    int j = idx / 3;
    return (*transform_tensor)(i, j);
    // return MatT[idx];
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::SetFace()
{
    nodes.resize(8);
    nodes_index.resize(8);
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::SetFace(CElem* onwer, const int Face)
{
    static int nodeIndex_loc[8];
    no_faces_on_surface = 0;
    owner = onwer;
    size_t n = owner->GetElementFaceNodes(Face, nodeIndex_loc);
    quadratic = owner->quadratic;
    face_index = Face;
    patch_index = owner->patch_index;
    switch (owner->geo_type)
    {
        // case MshElemType::LINE:  // 1-D bar element
        case MshElemType::QUAD:  // 2-D quadrilateral element
            this->setElementProperties(MshElemType::LINE);  // JOD 2014-11-10
            break;
        case MshElemType::HEXAHEDRON:  // 3-D hexahedral element
            this->setElementProperties(MshElemType::QUAD8);
            break;
        // case MshElemType::TRIANGLE:  // 2-D triagular element
        case MshElemType::TETRAHEDRON:  // 3-D tetrahedral element
            this->setElementProperties(MshElemType::TRIANGLE);
            break;
        case MshElemType::PRISM:
            if (Face < 2)
                this->setElementProperties(MshElemType::TRIANGLE);
            else
                this->setElementProperties(MshElemType::QUAD8);
            break;  // 3-D prismatic element
        case MshElemType::PYRAMID:
            if (Face < 1)
                this->setElementProperties(MshElemType::QUAD8);
            else
                this->setElementProperties(MshElemType::TRIANGLE);
            break;  // 3-D pyramid element
        default:
            std::cerr << "CElem::SetFace MshElemType not handled"
                      << "\n";
            break;
    }

    if (nodes.Size() < n)
        nodes.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        nodes[i] = owner->nodes[nodeIndex_loc[i]];
        nodes_index[i] = nodes[i]->GetIndex();
    }
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   11/2010 KR moved to MSHEnums.h
**************************************************************************/
std::string CElem::GetName() const
{
    return MshElemType2String(geo_type);
}

/**************************************************************************
   MSHLib-Method:
   Programing:
   06/2005 WW Implementation
   08/2005 WW/OK Extension for GMS/SOL files
   10/2008 OK FEFLOW
**************************************************************************/
void CElem::Read(std::istream& is, int fileType)
{
    // fileType=0: msh
    // fileType=1: rfi
    // fileType=2: gmsh
    // fileType=3: GMS
    // fileType=4: SOL
    // fileType=5: FLAC 3D. WW
    // fileType=5x: FLAC3D //MR
    // fileType=6: FEFLOW //OK
    // fileType=7: GMSH 2008 Version//TK
    int idummy(-1), et(-1);
    std::string name("");
    int gmsh_patch_index;  // TKOK

    //   is.ignore(numeric_limits<int>::max(), '\n');
    //----------------------------------------------------------------------
    // 1 Reading element type data
    switch (fileType)
    {
        //....................................................................
        case 0:  // msh
        {
            std::string buffer("");
            is >> index >> patch_index;
            is >> buffer;

            if (buffer.find("-1") != std::string::npos)
            {
                grid_adaptation = strtol(buffer.data(), NULL, 0);
                is >> name;
            }
            else
                name = buffer;

            geo_type = String2MshElemType(name);
            break;
        }
        //....................................................................
        case 1:  // rfi
            is >> index >> patch_index >> name;
            geo_type = String2MshElemType(name);
            break;
        //....................................................................
        case 2:  // gmsh
            is >> index >> et >> gmsh_patch_index >> idummy >> nnodes;
            patch_index = gmsh_patch_index - 1;  // OK
            switch (et)
            {
                case 1:
                    geo_type = MshElemType::LINE;
                    break;
                case 2:
                    geo_type = MshElemType::TRIANGLE;
                    break;
                case 3:
                    geo_type = MshElemType::QUAD;
                    break;
                case 4:
                    geo_type = MshElemType::TETRAHEDRON;
                    break;
                case 5:
                    geo_type = MshElemType::HEXAHEDRON;
                    break;
                case 6:
                    geo_type = MshElemType::PRISM;
                    break;
                case 7:
                    geo_type = MshElemType::PYRAMID;
                    break;
            }
            index--;
            break;
        case 7:  // GMSH 2008
            size_t nb_tags;

            is >> index >> et >> nb_tags >> idummy >> gmsh_patch_index;
            patch_index = gmsh_patch_index;
            for (size_t j = 2; j < nb_tags; j++)
                is >> idummy;
            switch (et)
            {
                case 1:
                    geo_type = MshElemType::LINE;
                    patch_index =
                        0;  // KR: can line elements have material ids?
                    nnodes = 2;
                    break;
                case 2:
                    geo_type = MshElemType::TRIANGLE;
                    nnodes = 3;
                    break;
                case 3:
                    geo_type = MshElemType::QUAD;
                    nnodes = 4;
                    break;
                case 4:
                    geo_type = MshElemType::TETRAHEDRON;
                    nnodes = 4;
                    break;
                case 5:
                    geo_type = MshElemType::HEXAHEDRON;
                    nnodes = 8;
                    break;
                case 6:
                    geo_type = MshElemType::PRISM;
                    nnodes = 6;
                    break;
                case 7:
                    geo_type = MshElemType::PYRAMID;
                    break;
                case 15:
                    geo_type = MshElemType::INVALID;
                    nnodes = 1;
                    break;
                default:
                    geo_type = MshElemType::INVALID;
                    break;
            }
            index--;
            break;
        //....................................................................
        case 3:  // GMS
            geo_type = MshElemType::TRIANGLE;
            break;
        //....................................................................
        case 4:  // gmsh
            geo_type = MshElemType::TRIANGLE;
            break;
        //....................................................................
        case 5:  // FLAC 3D. 14.01.2008 WW
            geo_type = MshElemType::HEXAHEDRON;
            break;
        case 56:                            // FLAC3D - pri (Wedge)        //MR
            geo_type = MshElemType::PRISM;  // MR
            fileType = 5;                   // MR
            break;                          // MR
    }

    // TF
    if (geo_type == MshElemType::INVALID)
    {
        // read rest of line
        std::string tmp;
        getline(is, tmp);
        return;
    }

    //----------------------------------------------------------------------
    // 2 Element configuration
    this->setElementProperties(geo_type);

    //----------------------------------------------------------------------
    // 3 Reading element node data
    switch (fileType)
    {
        case 0:  // msh
            for (int i = 0; i < nnodes; i++)
                is >> nodes_index[i];
            break;
        case 1:  // rfi
            for (int i = 0; i < nnodes; i++)
                is >> nodes_index[i];
            break;
        case 2:  // gmsh
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            break;
        case 7:  // GMSH 2008
            if (et != 15)
                for (int i = 0; i < nnodes; i++)
                {
                    is >> nodes_index[i];
                    nodes_index[i] -= 1;
                }
            else
            {
                // eat rest of line
                std::string dummy;
                is >> dummy;
            }
            break;
        case 3:  // GMS
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            break;
        case 4:  // SOL
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            is >> patch_index;
            break;
        case 5:  // FLAC 3D. 14.01.2008. WW
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            break;
        case 6:  // FEFLOLW
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            break;
        case 8:  // GMS_3DM
            is >> idummy;
            for (int i = 0; i < nnodes; i++)
            {
                is >> nodes_index[i];
                nodes_index[i] -= 1;
            }
            is >> patch_index;
            patch_index -= 1;
            break;
    }
    is >> std::ws;
    //----------------------------------------------------------------------
    // Initialize topological properties
    neighbors.resize(nfaces);
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i] = NULL;
    edges.resize(nedges);
    edges_orientation.resize(nedges);
    for (size_t i = 0; i < nedges; i++)
    {
        edges[i] = NULL;
        edges_orientation[i] = 1;
    }
}

void CElem::InitializeMembers()
{
    // Initialize topological properties
    neighbors.resize(nfaces);
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i] = NULL;
    edges.resize(nedges);
    edges_orientation.resize(nedges);
    for (size_t i = 0; i < nedges; i++)
    {
        edges[i] = NULL;
        edges_orientation[i] = 1;
    }
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::WriteIndex(std::ostream& os) const
{
    // Comment for GUI WW if(quadratic) nn = nnodesHQ;
    os << index << " " << patch_index << " " << GetName() << " ";
    for (int i = 0; i < nnodes - 1; i++)
        os << nodes_index[i] << " ";
    os << nodes_index[nnodes - 1] << "\n";
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::WriteIndex_TEC(std::ostream& os) const
{
    std::string deli = "  ";
    if (geo_type == MshElemType::LINE)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli
           << nodes_index[1] + 1 << deli << nodes_index[0] + 1;

    else if (geo_type == MshElemType::TRIANGLE)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli
           << nodes_index[2] + 1 << deli << nodes_index[0] + 1;
    else if (geo_type == MshElemType::PRISM)
        os << nodes_index[0] + 1 << deli << nodes_index[0] + 1 << deli
           << nodes_index[1] + 1 << deli << nodes_index[2] + 1 << deli
           << nodes_index[3] + 1 << deli << nodes_index[3] + 1 << deli
           << nodes_index[4] + 1 << deli << nodes_index[5] + 1 << deli;
    else if (geo_type == MshElemType::PYRAMID)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli
           << nodes_index[2] + 1 << deli << nodes_index[3] + 1 << deli
           << nodes_index[4] + 1 << deli << nodes_index[4] + 1 << deli
           << nodes_index[4] + 1 << deli << nodes_index[4] + 1 << deli;
    else

        for (int i = 0; i < nnodes; i++)
            os << nodes_index[i] + 1 << deli;
    os << '\n';  // GK44 for io performance do not flush buffer with endl...
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::WriteAll(std::ostream& os) const
{
    std::string deli = "  ";
    os << index << deli << patch_index << deli << GetName() << deli;
    os << "Index X Y Z: "
       << "\n";
    for (size_t i = 0; i < nodes.Size(); i++)
    {
        double const* const pnt_i(nodes[i]->getData());
        os << nodes_index[i] << deli << pnt_i[0] << deli << pnt_i[1] << deli
           << pnt_i[2] << "\n";
    }
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::WriteNeighbors(std::ostream& os) const
{
    os << "Neighbors of " << index << "\n";
    for (size_t i = 0; i < nfaces; i++)
        neighbors[i]->WriteAll(os);
    os << "End neighbors of " << index << "\n"
       << "\n";
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::MarkingAll(bool makop)
{
    this->mark = makop;
    int SizeV = nnodes;
    if (quadratic)
        SizeV = nnodesHQ;
    for (int i = 0; i < SizeV; i++)
        nodes[i]->SetMark(makop);

    size_t nedg = edges.Size();
    for (size_t i = 0; i < nedg; i++)
        edges[i]->SetMark(makop);
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::SetNodes(Math_Group::vec<CNode*>& ele_nodes, bool ReSize)
{
    int SizeV = nnodes;
    if (quadratic)
        SizeV = nnodesHQ;
    if (ReSize)
    {
        nodes.resize(SizeV);
#if !defined(USE_PETSC)  // && !defined (other parallel linear solver lib).
                         // //WW. 05.2013
        nodes_index.resize(SizeV);
#endif
    }
    for (int i = 0; i < SizeV; i++)
    {
        nodes[i] = ele_nodes[i];
#if !defined(USE_PETSC)  // && !defined (other parallel linear solver lib).
                         // //WW. 05.2013
        nodes_index[i] = nodes[i]->GetIndex();
#endif
    }
}

void CElem::setNodes(std::vector<CNode*> const& ele_nodes)
{
    if (nodes.Size() != ele_nodes.size())
    {
        nodes.resize(ele_nodes.size());
        nodes_index.resize(ele_nodes.size());
    }

    for (size_t k(0); k < ele_nodes.size(); k++)
    {
        nodes[k] = ele_nodes[k];
        nodes_index[k] = ele_nodes[k]->GetIndex();
    }
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CElem::GetLocalIndicesOfEdgeNodes(const int Edge, int* EdgeNodes)
{
    switch (geo_type)
    {
        case MshElemType::LINE:
            EdgeNodes[0] = 0;
            EdgeNodes[1] = 1;
            break;               // 1-D bar element
        case MshElemType::QUAD:  // 2-D quadrilateral element
            EdgeNodes[0] = Edge;
            EdgeNodes[1] = (Edge + 1) % 4;
            break;
        case MshElemType::HEXAHEDRON:  // 3-D hexahedral element
            if (Edge < 8)
            {
                EdgeNodes[0] = Edge;
                EdgeNodes[1] = (Edge + 1) % 4 + 4 * (int)(Edge / 4);
            }
            else
            {
                EdgeNodes[0] = Edge % 4;
                EdgeNodes[1] = Edge % 4 + 4;
            }
            break;
        case MshElemType::TRIANGLE:  // 2-D triagular element
            EdgeNodes[0] = Edge;
            EdgeNodes[1] = (Edge + 1) % 3;
            break;
        case MshElemType::TETRAHEDRON:  // 3-D tetrahedra
            if (Edge < 3)
            {
                EdgeNodes[0] = Edge;
                EdgeNodes[1] = (Edge + 1) % 3;
            }
            else
            {
                EdgeNodes[0] = 3;
                EdgeNodes[1] = (Edge + 1) % 3;
            }

            break;
        case MshElemType::PRISM:  // 3-D prismatic element
            if (Edge < 6)
            {
                EdgeNodes[0] = Edge;
                EdgeNodes[1] = (Edge + 1) % 3 + 3 * (int)(Edge / 3);
            }
            else
            {
                EdgeNodes[0] = Edge % 3;
                EdgeNodes[1] = Edge % 3 + 3;
            }
            break;
        case MshElemType::PYRAMID:  // 3-D pyramid element
            if (Edge < 4)
            {
                EdgeNodes[0] = Edge;
                EdgeNodes[1] = (Edge + 1) % 4;
            }
            else
            {
                EdgeNodes[0] = Edge % 4;
                EdgeNodes[1] = 4;
            }
            break;
        default:
            std::cerr << "CElem::GetLocalIndicesOfEdgeNodes() - MshElemType "
                         "not handled"
                      << "\n";
    }
}
/**************************************************************************
   GetElementFaceNodes
   Task: Get local indeces of an element face nodes
   Return: number of nodes of a face
   Programing:
   06/2004 WW
**************************************************************************/
int CElem::GetElementFaces1D(int* FaceNode)
{
    FaceNode[0] = 0;
    FaceNode[1] = 1;
    return 2;
}
/**************************************************************************
   GetElementFaceNodesTri
   Task: Get local indeces of a traingle element face nodes
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFacesTri(int Face, int* FaceNode)
{
    if (!quadratic)
    {
        FaceNode[0] = Face;
        FaceNode[1] = (Face + 1) % 3;
        return 2;
    }
    else
    {
        FaceNode[0] = Face;
        FaceNode[1] = (Face + 1) % 3;
        FaceNode[2] = Face + 3;
        return 3;
    }
}
/**************************************************************************
   GetElementFaceNodesQuad
   Task: Get local indeces of a quadralateral element face nodes
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFacesQuad(int Face, int* FaceNode)
{
    if (!quadratic)
    {
        FaceNode[0] = Face;
        FaceNode[1] = (Face + 1) % 4;
        return 2;
    }
    else
    {
        FaceNode[0] = Face;
        FaceNode[1] = (Face + 1) % 4;
        FaceNode[2] = Face + 4;
        return 3;
    }
}
/**************************************************************************
   GetElementFaceNodesHex
   Task: Get local indeces of a hexahedra element face nodes
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFacesHex(int Face, int* FaceNode)
{
    int nn = 4, k = 0;
    if (quadratic)
        nn = 8;
    switch (Face)
    {
        case 0:
            for (k = 0; k < 4; k++)
                FaceNode[k] = k;
            if (quadratic)
                for (k = 0; k < 4; k++)
                    FaceNode[k + 4] = k + 8;
            break;
        case 1:
            for (k = 0; k < 4; k++)
                FaceNode[k] = k + 4;
            if (quadratic)
                for (k = 0; k < 4; k++)
                    FaceNode[k + 4] = k + 12;
            break;
        case 2:
            FaceNode[0] = 0;
            FaceNode[1] = 4;
            FaceNode[2] = 5;
            FaceNode[3] = 1;
            if (quadratic)
            {
                FaceNode[4] = 16;
                FaceNode[5] = 12;
                FaceNode[6] = 17;
                FaceNode[7] = 8;
            }
            break;
        case 3:
            FaceNode[0] = 1;
            FaceNode[1] = 5;
            FaceNode[2] = 6;
            FaceNode[3] = 2;
            if (quadratic)
            {
                FaceNode[4] = 17;
                FaceNode[5] = 13;
                FaceNode[6] = 18;
                FaceNode[7] = 9;
            }

            break;
        case 4:
            FaceNode[0] = 2;
            FaceNode[1] = 6;
            FaceNode[2] = 7;
            FaceNode[3] = 3;
            if (quadratic)
            {
                FaceNode[4] = 18;
                FaceNode[5] = 14;
                FaceNode[6] = 19;
                FaceNode[7] = 10;
            }
            break;
        case 5:
            FaceNode[0] = 0;
            FaceNode[1] = 3;
            FaceNode[2] = 7;
            FaceNode[3] = 4;
            if (quadratic)
            {
                FaceNode[4] = 11;
                FaceNode[5] = 19;
                FaceNode[6] = 15;
                FaceNode[7] = 16;
            }
            break;
    }
    return nn;
}
/**************************************************************************
   GetElementFaceNodesTet
   Task: Get local indeces of a Tedrahedra element face nodes
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFacesTet(int Face, int* FaceNode)
{
    int nn = 3;
    if (quadratic)
        nn = 6;
    switch (Face)
    {
        case 0:
            FaceNode[0] = 1;
            FaceNode[1] = 2;
            FaceNode[2] = 3;
            if (quadratic)
            {
                FaceNode[3] = 5;
                FaceNode[4] = 8;
                FaceNode[5] = 7;
            }
            break;
        case 1:
            FaceNode[0] = 3;
            FaceNode[1] = 2;
            FaceNode[2] = 0;
            if (quadratic)
            {
                FaceNode[3] = 8;
                FaceNode[4] = 6;
                FaceNode[5] = 9;
            }
            break;
        case 2:
            FaceNode[0] = 1;
            FaceNode[1] = 3;
            FaceNode[2] = 0;
            if (quadratic)
            {
                FaceNode[3] = 7;
                FaceNode[4] = 9;
                FaceNode[5] = 4;
            }
            break;
        case 3:
            FaceNode[0] = 0;
            FaceNode[1] = 2;
            FaceNode[2] = 1;
            if (quadratic)
            {
                FaceNode[3] = 6;
                FaceNode[4] = 5;
                FaceNode[5] = 4;
            }
            break;
    }
    return nn;
}

/**************************************************************************
   GetElementFaceNodesPri
   Task: Get local indeces of a prismal element face nodes
   Augs.:
        int Face :  Local index of element face
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFacesPri(int Face, int* FaceNode)
{
    int nn = 3, k = 0;
    switch (Face)
    {
        case 0:
            nn = 3;
            for (k = 0; k < 3; k++)
                FaceNode[k] = k;
            if (quadratic)
            {
                for (k = 0; k < 3; k++)
                    FaceNode[k + 3] = k + 6;
                nn = 6;
            }
            break;
        case 1:
            for (k = 0; k < 3; k++)
                FaceNode[k] = k + 3;
            nn = 3;
            if (quadratic)
            {
                for (k = 0; k < 3; k++)
                    FaceNode[k + 3] = k + 9;
                nn = 6;
            }
            break;
        case 2:
            FaceNode[0] = 1;
            FaceNode[1] = 2;
            FaceNode[2] = 5;
            FaceNode[3] = 4;
            nn = 4;
            if (quadratic)
            {
                FaceNode[4] = 7;
                FaceNode[5] = 14;
                FaceNode[6] = 10;
                FaceNode[7] = 13;
                nn = 8;
            }
            break;
        case 3:
            FaceNode[0] = 5;
            FaceNode[1] = 2;
            FaceNode[2] = 0;
            FaceNode[3] = 3;
            nn = 4;
            if (quadratic)
            {
                FaceNode[4] = 14;
                FaceNode[5] = 8;
                FaceNode[6] = 12;
                FaceNode[7] = 11;
                nn = 8;
            }
            break;
        case 4:
            FaceNode[0] = 0;
            FaceNode[1] = 1;
            FaceNode[2] = 4;
            FaceNode[3] = 3;
            nn = 4;
            if (quadratic)
            {
                FaceNode[4] = 6;
                FaceNode[5] = 13;
                FaceNode[6] = 9;
                FaceNode[7] = 12;
                nn = 8;
            }
            break;
    }
    return nn;
}
/**************************************************************************
   GetElementFacesPyra
   Task: Get local indeces of a element face nodes
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   01/2010 NW
**************************************************************************/
int CElem::GetElementFacesPyramid(const int Face, int* FaceNode)
{
    int nn = 3, k = 0;
    switch (Face)
    {
        case 0:
            nn = 4;
            for (k = 0; k < nn; k++)
                FaceNode[k] = k;
            if (quadratic)
            {
                for (k = 0; k < nn; k++)
                    FaceNode[k + 4] = k + 5;
                nn = 8;
            }
            break;
        case 1:
            FaceNode[0] = 0;
            FaceNode[1] = 1;
            FaceNode[2] = 4;
            nn = 3;
            if (quadratic)
            {
                FaceNode[3] = 5;
                FaceNode[4] = 10;
                FaceNode[5] = 9;
                nn = 6;
            }
            break;
        case 2:
            FaceNode[0] = 1;
            FaceNode[1] = 2;
            FaceNode[2] = 4;
            nn = 3;
            if (quadratic)
            {
                FaceNode[3] = 6;
                FaceNode[4] = 11;
                FaceNode[5] = 10;
                nn = 6;
            }
            break;
        case 3:
            FaceNode[0] = 2;
            FaceNode[1] = 3;
            FaceNode[2] = 4;
            nn = 3;
            if (quadratic)
            {
                FaceNode[3] = 7;
                FaceNode[4] = 12;
                FaceNode[5] = 11;
                nn = 6;
            }
            break;
        case 4:
            FaceNode[0] = 3;
            FaceNode[1] = 0;
            FaceNode[2] = 4;
            nn = 3;
            if (quadratic)
            {
                FaceNode[3] = 8;
                FaceNode[4] = 9;
                FaceNode[5] = 12;
                nn = 6;
            }
            break;
    }
    return nn;
}
/**************************************************************************
   GetElementFaces
   Task: set element faces (Geometry)
   Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes

   Programing:
   09/2004 WW
**************************************************************************/
int CElem::GetElementFaceNodes(int Face, int* FacesNode)
{
    switch (geo_type)
    {
        case MshElemType::LINE:  // 1-D bar element
            return GetElementFaces1D(FacesNode);
        case MshElemType::QUAD:  // 2-D quadrilateral element
            return GetElementFacesQuad(Face, FacesNode);
        case MshElemType::HEXAHEDRON:  // 3-D hexahedral element
            return GetElementFacesHex(Face, FacesNode);
        case MshElemType::TRIANGLE:  // 2-D triagular element
            return GetElementFacesTri(Face, FacesNode);
        case MshElemType::TETRAHEDRON:  // 3-D tetrahedral element
            return GetElementFacesTet(Face, FacesNode);
        case MshElemType::PRISM:
            return GetElementFacesPri(Face, FacesNode);
        // 3-D prismatic element
        case MshElemType::PYRAMID:  // 3-D pyramid element
            return GetElementFacesPyramid(Face, FacesNode);
        default:
            std::cerr << "CElem::GetElementFaceNodes MshElemType not handled"
                      << "\n";
    }
    return 0;
}

/**************************************************************************
   MSHLib-Method:
   Task: Compute volume for elements with straight edges and surfaces
   Programing:
   06/2005 WW Implementation
   03/2011 KR cleaned up code
**************************************************************************/
void CElem::ComputeVolume()
{
    volume = calcVolume();

    if (this->geo_type == MshElemType::LINE)  // Line
        representative_length = volume;
    else if (this->geo_type == MshElemType::TRIANGLE)
        representative_length = sqrt(volume) * 4.0;
    else if (this->geo_type == MshElemType::QUAD)
        representative_length = sqrt(volume);  // kg44 reactivated
    else if (this->geo_type == MshElemType::TETRAHEDRON)
        representative_length = sqrt(volume) * 6.0;
    else if (this->geo_type == MshElemType::HEXAHEDRON)
        representative_length = pow(volume, 1. / 3.);
    else if (this->geo_type == MshElemType::PRISM)
        representative_length = pow(volume, 1. / 3.);
    else if (this->geo_type == MshElemType::PYRAMID)
        representative_length = 0.0;  // NW set zero because I don't know what
                                      // is the representative_length.
    else
        std::cerr << "Error in CElem::ComputeVolume() - MshElemType not found"
                  << "\n";
}

double CElem::calcVolume() const
{
    double elemVolume = 0.0;

    if (this->geo_type == MshElemType::LINE)  // Line
    {
        double const* const pnt0(nodes[0]->getData());
        double const* const pnt(nodes[nnodes - 1]->getData());
        double xDiff = pnt[0] - pnt0[0];
        double yDiff = pnt[1] - pnt0[1];
        double zDiff = pnt[2] - pnt0[2];
        elemVolume = sqrt(xDiff * xDiff + yDiff * yDiff +
                          zDiff * zDiff);  // CMCD kg44 reactivated
    }
    else if (this->geo_type == MshElemType::TRIANGLE)
        elemVolume = ComputeDetTri(nodes[0]->getData(), nodes[1]->getData(),
                                   nodes[2]->getData());  // kg44 reactivated
    else if (this->geo_type == MshElemType::QUAD)
        elemVolume = ComputeDetTri(nodes[0]->getData(), nodes[1]->getData(),
                                   nodes[2]->getData()) +
                     ComputeDetTri(nodes[2]->getData(), nodes[3]->getData(),
                                   nodes[0]->getData());
    else if (this->geo_type == MshElemType::TETRAHEDRON)
        elemVolume = ComputeDetTex(nodes[0]->getData(), nodes[1]->getData(),
                                   nodes[2]->getData(),
                                   nodes[3]->getData());  // kg44 reactivated
    else if (this->geo_type == MshElemType::HEXAHEDRON)
    {
        elemVolume = ComputeDetTex(nodes[4]->getData(), nodes[7]->getData(),
                                   nodes[5]->getData(), nodes[0]->getData());
        elemVolume += ComputeDetTex(nodes[5]->getData(), nodes[3]->getData(),
                                    nodes[1]->getData(), nodes[0]->getData());
        elemVolume += ComputeDetTex(nodes[5]->getData(), nodes[7]->getData(),
                                    nodes[3]->getData(), nodes[0]->getData());
        elemVolume += ComputeDetTex(nodes[5]->getData(), nodes[7]->getData(),
                                    nodes[6]->getData(), nodes[2]->getData());
        elemVolume += ComputeDetTex(nodes[1]->getData(), nodes[3]->getData(),
                                    nodes[5]->getData(), nodes[2]->getData());
        elemVolume += ComputeDetTex(nodes[3]->getData(), nodes[7]->getData(),
                                    nodes[5]->getData(), nodes[2]->getData());
        // kg44 reactivated
    }
    else if (this->geo_type == MshElemType::PRISM)
    {
        elemVolume = ComputeDetTex(nodes[0]->getData(), nodes[1]->getData(),
                                   nodes[2]->getData(), nodes[3]->getData());
        elemVolume += ComputeDetTex(nodes[1]->getData(), nodes[4]->getData(),
                                    nodes[2]->getData(), nodes[3]->getData());
        elemVolume += ComputeDetTex(nodes[2]->getData(), nodes[4]->getData(),
                                    nodes[5]->getData(), nodes[3]->getData());
        // kg44 reactivated ---------Here the direction of flow needs to be
        // taken into account, we need rep length in x,y,z direction
    }
    else if (this->geo_type == MshElemType::PYRAMID)
    {
        elemVolume = ComputeDetTex(nodes[0]->getData(), nodes[1]->getData(),
                                   nodes[2]->getData(), nodes[4]->getData());
        elemVolume += ComputeDetTex(nodes[2]->getData(), nodes[3]->getData(),
                                    nodes[0]->getData(), nodes[4]->getData());
    }
    else
        std::cerr << "Error in CElem::ComputeVolume() - MshElemType not found"
                  << "\n";

    return elemVolume;
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
**************************************************************************/
void CElem::FaceNormal(int index0, int index1, double* face)
{
    if (GetElementType() == MshElemType::TRIANGLE ||
        GetElementType() == MshElemType::QUAD)
    {
        double xx[3];
        double yy[3];
        double zz[3];

        //----plane normal----------------------------
        // transform_tensor = new Matrix(3,3);
        // face"_vec
        double const* pnt0(nodes[index0]->getData());
        double const* pnt1(nodes[index1]->getData());
        face[0] = pnt1[0] - pnt0[0];
        face[1] = pnt1[1] - pnt0[1];
        face[2] = pnt1[2] - pnt0[2];

        // x"_vec
        //         xx[0] = nodes[1]->X()-nodes[0]->X();
        //         xx[1] = nodes[1]->Y()-nodes[0]->Y();
        //         xx[2] = nodes[1]->Z()-nodes[0]->Z();
        pnt0 = nodes[0]->getData();
        pnt1 = nodes[1]->getData();
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        xx[2] = pnt1[2] - pnt0[2];
        // NormalizeVector(xx,3);
        // a vector on the plane
        //         yy[0] = nodes[2]->X()-nodes[1]->X();
        //         yy[1] = nodes[2]->Y()-nodes[1]->Y();
        //         yy[2] = nodes[2]->Z()-nodes[1]->Z();
        double const* pnt2(nodes[2]->getData());
        yy[0] = pnt2[0] - pnt1[0];
        yy[1] = pnt2[1] - pnt1[1];
        yy[2] = pnt2[2] - pnt1[2];
        // z"_vec

        CrossProduction(xx, yy, zz);
        NormalizeVector(zz, 3);
        // y"_vec
        CrossProduction(face, zz, yy);
        NormalizeVector(yy, 3);
        for (size_t i = 0; i < 3; i++)
            face[i] = yy[i];
    }
}

/**************************************************************************
   MSHLib-Method:
   06/2006 OK Implementation
**************************************************************************/
void CElem::SetNormalVector()
{
    if (!normal_vector)
        normal_vector = new double[3];                // WW
    if (this->GetElementType() == MshElemType::LINE)  // JOD 2014-11-10
    {
        double const* const p0(nodes[0]->getData());
        double const* const p1(nodes[1]->getData());
        double v1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        if (fabs(v1[2]) > 1.e-20)
        {
            double buffer = v1[1];
            v1[1] = v1[2];
            v1[2] = buffer;
        }
        const double v2[3] = {0, 0, 1};  // fluxes are on xy plane
        CrossProduction(v1, v2, normal_vector);
        NormalizeVector(normal_vector, 3);
        if (normal_vector[0] < 0)
            normal_vector[0] = -normal_vector[0];
    }
    if (this->GetElementType() == MshElemType::TRIANGLE)
    {
        double const* const p0(nodes[0]->getData());
        double const* const p1(nodes[1]->getData());
        double const* const p2(nodes[2]->getData());
        const double v1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const double v2[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        CrossProduction(v1, v2, normal_vector);
        NormalizeVector(normal_vector, 3);
    }
    if (this->GetElementType() ==
        MshElemType::QUAD)  // BG 04/2011: calculation of the normal vector of a
                            // quad element
    {
        double const* const p0(nodes[0]->getData());
        double const* const p1(nodes[1]->getData());
        //		double const* const p2 (nodes[2]->getData()); // TF unused
        // variable
        double const* const p3(nodes[3]->getData());
        const double v1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const double v2[3] = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
        // v1[0] = nodes[1]->X() - nodes[0]->X();
        // v1[1] = nodes[1]->Y() - nodes[0]->Y();
        // v1[2] = nodes[1]->Z() - nodes[0]->Z();
        // v2[0] = nodes[3]->X() - nodes[0]->X();
        // v2[1] = nodes[3]->Y() - nodes[0]->Y();
        // v2[2] = nodes[3]->Z() - nodes[0]->Z();
        CrossProduction(v1, v2, normal_vector);
        NormalizeVector(normal_vector, 3);
    }
}

// KR 2010/11/16
void CElem::setElementProperties(MshElemType::type t)
{
    switch (t)
    {
        case MshElemType::LINE:
            nnodes = 2;
            nnodesHQ = 3;
            ele_dim = 1;
            geo_type = MshElemType::LINE;
            nfaces = 2;
            nedges = 1;
            break;
        case MshElemType::QUAD8:
            nnodes = 4;
            nnodesHQ = 8;
            ele_dim = 2;
            geo_type = MshElemType::QUAD8;
            nfaces = 4;
            nedges = 4;
            break;
        case MshElemType::QUAD:
            nnodes = 4;
            nnodesHQ = 9;
            ele_dim = 2;
            geo_type = MshElemType::QUAD;
            nfaces = 4;
            nedges = 4;
            break;
        case MshElemType::HEXAHEDRON:
            nnodes = 8;
            nnodesHQ = 20;
            ele_dim = 3;
            nfaces = 6;
            nedges = 12;
            geo_type = MshElemType::HEXAHEDRON;
            break;
        case MshElemType::TRIANGLE:
            nnodes = 3;
            nnodesHQ = 6;
            ele_dim = 2;
            geo_type = MshElemType::TRIANGLE;
            nfaces = 3;
            nedges = 3;
            break;
        case MshElemType::TETRAHEDRON:
            nnodes = 4;
            nnodesHQ = 10;
            ele_dim = 3;
            geo_type = MshElemType::TETRAHEDRON;
            nfaces = 4;
            nedges = 6;
            break;
        case MshElemType::PRISM:
            nnodes = 6;
            nnodesHQ = 15;
            ele_dim = 3;
            geo_type = MshElemType::PRISM;
            nfaces = 5;
            nedges = 9;
            break;
        case MshElemType::PYRAMID:
            nnodes = 5;
            nnodesHQ = 13;
            ele_dim = 3;
            geo_type = MshElemType::PYRAMID;
            nfaces = 5;
            nedges = 8;
            break;
        default:
            std::cerr << "CElem::setElementProperties MshElemType not handled"
                      << "\n";
    }
    this->nodes_index.resize(quadratic ? nnodesHQ : nnodes);
}

// NW
double* CElem::ComputeGravityCenter()
{
    const size_t nnodes0 = this->nnodes;
    for (size_t i = 0; i < nnodes0; i++)  // Nodes
    {
        double const* const pnt(nodes[i]->getData());
        gravity_center[0] += pnt[0];
        gravity_center[1] += pnt[1];
        gravity_center[2] += pnt[2];
    }
    gravity_center[0] /= (double)nnodes0;
    gravity_center[1] /= (double)nnodes0;
    gravity_center[2] /= (double)nnodes0;
    return gravity_center;
}

/**************************************************************************
MSHLib-Method:
11/2014 JOD Implementation
**************************************************************************/
void CElem::DirectNormalVector()
{
    for (int i = 0; i < 3; i++)
        if (normal_vector[i] < 0)
            InvertNormalVector();
}

/**************************************************************************
MSHLib-Method:
11/2014 JOD Implementation
**************************************************************************/
void CElem::InvertNormalVector()
{
    normal_vector[0] = -normal_vector[0];
    normal_vector[1] = -normal_vector[1];
    normal_vector[2] = -normal_vector[2];
}

}  // namespace MeshLib
