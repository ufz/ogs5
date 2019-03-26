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
   08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
// C++
#include <string>
#include <vector>

// MSHLib
#include "msh_lib.h"

#include "display.h"

#include "memory.h"

// GEOLib
#include "files0.h"
#include "geo_lib.h"
// PCSLib
#include "mathlib.h"
#include "rf_mmp_new.h"  //OK411

extern int ReadRFIFile(std::string g_strFileNameBase);
#include "rf_pcs.h"

std::vector<MeshLib::CFEMesh*> fem_msh_vector;

#define FEM_FILE_EXTENSION ".msh"

using namespace Math_Group;

/**************************************************************************
   FEMLib-Method:
   03/2005 OK Implementation
   05/2005 TK modified
   05/2006 TK modified
**************************************************************************/
void FEMDeleteAll()
{
    for (int i = 0; i < (int)fem_msh_vector.size(); i++)
    {
        delete fem_msh_vector[i];
        fem_msh_vector[i] = NULL;
    }
    fem_msh_vector.clear();
}

#ifndef USE_PETSC  // && not defined(other parallel method with ddc)
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   08/2005 WW Topology construction and rfi compatible
   10/2005 OK BINARY
   08/2010 KR deleted binary mesh read
   03/2011 KR cleaned up code
   08/2011 WW Recovery multi-mesh
   09/2011 TF changed signature of function in order to read more than one mesh
**************************************************************************/
void FEMRead(const std::string& file_base_name, std::vector<CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj, std::string* unique_name)
{
    CFEMesh* mesh(NULL);
    std::string msh_file_name(file_base_name + FEM_FILE_EXTENSION);

    std::ifstream msh_file_ascii(msh_file_name.data(), std::ios::in);
    if (!msh_file_ascii.is_open())
        std::cout << "CFEMesh::FEMRead() - Could not open file...\n";

    Display::ScreenMessage("MSHRead:  ASCII file\n");

    std::string line_string("");
    getline(msh_file_ascii, line_string);

    bool more_mesh = false;                                 // 12.08.2011. WW
    if (line_string.find("#FEM_MSH") != std::string::npos)  // OGS mesh file
    {
        mesh = new CFEMesh(geo_obj, unique_name);
        more_mesh = mesh->Read(&msh_file_ascii);
        mesh_vec.push_back(mesh);  // TF

        // Multi-mesh 12.08.2011 WW
        if (more_mesh)
            while (!msh_file_ascii.eof())
            {
                // getline(msh_file_ascii, line_string);
                // if(line_string.find("#FEM_MSH")!=std::string::npos)
                mesh = new CFEMesh(geo_obj, unique_name);
                more_mesh = mesh->Read(&msh_file_ascii);
                mesh_vec.push_back(mesh);  // TF
                if (!more_mesh)
                    break;
            }
        //  if(line_string.find("#STOP")!=std::string::npos)
        //     break;
    }
    else  // RFI mesh file
    {
        msh_file_ascii.seekg(0L, std::ios::beg);
        mesh = new CFEMesh(geo_obj, unique_name);
        Read_RFI(msh_file_ascii, mesh);
        mesh_vec.push_back(mesh);  // 12.08.2011 WW
    }

    msh_file_ascii.close();
}
#endif

/**************************************************************************
   MSHLib-Method: Read rfi file ()
   Task:
   Programing:
   08/2005 WW Re-implememtation
**************************************************************************/
void Read_RFI(std::istream& msh_file, CFEMesh* m_msh)
{
    long id;
    long i = 0;
    int NumNodes = 0;
    int NumElements = 0;
    int End = 1;
    double x, y, z;
    std::string strbuffer;

    MeshLib::CNode* node = NULL;
    MeshLib::CElem* elem = NULL;
    //----------------------------------------------------------------------
    while (End)
    {
        getline(msh_file, strbuffer);  // The first line
        msh_file >> i >> NumNodes >> NumElements >> std::ws;
        //....................................................................
        // Node data
        for (i = 0; i < NumNodes; i++)
        {
            msh_file >> id >> x >> y >> z >> std::ws;
            node = new MeshLib::CNode(id, x, y, z);
            m_msh->nod_vector.push_back(node);
        }
        for (i = 0; i < NumElements; i++)
        {
            elem = new MeshLib::CElem(i);
            elem->Read(msh_file, 1);
            m_msh->ele_vector.push_back(elem);
        }
        End = 0;
    }
}

/**************************************************************************
   MSHLib-Method:
   02/2006 WW Implementation
**************************************************************************/
void CompleteMesh()
{
    for (int i = 0; i < (int)fem_msh_vector.size(); i++)
    {
        fem_msh_vector[i]->ConstructGrid();
        fem_msh_vector[i]->FillTransformMatrix();
    }
}

/**************************************************************************
   FEMLib-Method:
   Task: Master write functionn
   Programing:
   03/2005 OK Implementation
   10/2005 OK BINARY
   last modification:
   08/2010	KR binary case deleted
**************************************************************************/
void MSHWrite(std::string file_base_name)
{
    // File handling
    std::string fem_msh_file_name = file_base_name + FEM_FILE_EXTENSION;
    std::ofstream fem_msh_file(fem_msh_file_name.c_str(),
                               std::ios::trunc | std::ios::out);

    if (!fem_msh_file.good())
        return;

    for (size_t i = 0; i < fem_msh_vector.size(); i++)
    {
        FileIO::OGSMeshIO meshIO;
        meshIO.setPrecision(12);
        meshIO.setFormat(std::ios::scientific);
        meshIO.setMesh(fem_msh_vector[i]);
        fem_msh_file << meshIO.writeToString();
    }
    fem_msh_file << "#STOP";
    fem_msh_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modification:
**************************************************************************/
MeshLib::CFEMesh* FEMGet(const std::string& msh_name)
{
    size_t no_msh = fem_msh_vector.size();
    // If there is only one msh file available, use it for all process. WW
    if (no_msh == 1)
        return fem_msh_vector[0];  // WW
    for (size_t i = 0; i < no_msh; i++)
        if (fem_msh_vector[i]->pcs_name.compare(msh_name) == 0)
            return fem_msh_vector[i];
    return NULL;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
   11/2005 OK OO-ELE
**************************************************************************/
void MSHWriteTecplot()
{
    MshElemType::type ele_type = MshElemType::INVALID;
    long no_nodes;
    long no_elements;
    std::string delimiter(", ");
    long i;
    MeshLib::CElem* m_ele = NULL;
    vec<long> node_indeces(8);
    //----------------------------------------------------------------------
    // File handling
    std::string file_path = "MSH";
    //----------------------------------------------------------------------
    MeshLib::CFEMesh* m_msh = NULL;
    for (int j = 0; j < (int)fem_msh_vector.size(); j++)
    {
        m_msh = fem_msh_vector[j];
        no_nodes = (long)m_msh->nod_vector.size();
        no_elements = (long)m_msh->ele_vector.size();
        // Test ele_type
        if (no_elements > 0)
        {
            m_ele = m_msh->ele_vector[0];
            ele_type = m_ele->GetElementType();
        }
        // File handling
        std::string msh_file_name =
            file_path + "_" + m_msh->pcs_name + TEC_FILE_EXTENSION;
        std::fstream msh_file(msh_file_name.data(),
                              std::ios::trunc | std::ios::out);
        msh_file.setf(std::ios::scientific, std::ios::floatfield);
        msh_file.precision(12);
        if (!msh_file.good())
            return;
        msh_file.seekg(0L, std::ios::beg);
        msh_file << "VARIABLES = X,Y,Z"
                 << "\n";
        msh_file << "ZONE T = " << m_msh->pcs_name << delimiter
                 << "N = " << no_nodes << delimiter << "E = " << no_elements
                 << delimiter;
        msh_file << "F = FEPOINT" << delimiter;
        switch (ele_type)
        {
            //..................................................................
            case MshElemType::LINE:
                msh_file << "ET = QUADRILATERAL"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    msh_file << node_indeces[0] + 1 << " "
                             << node_indeces[1] + 1 << " "
                             << node_indeces[1] + 1 << " "
                             << node_indeces[0] + 1 << "\n";
                }
                break;
            //..................................................................
            case MshElemType::QUAD:
                msh_file << "ET = QUADRILATERAL"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    msh_file << node_indeces[0] + 1 << " "
                             << node_indeces[1] + 1 << " "
                             << node_indeces[2] + 1 << " "
                             << node_indeces[3] + 1 << "\n";
                }
                break;
            //..................................................................
            case MshElemType::HEXAHEDRON:
                msh_file << "ET = BRICK"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    msh_file
                        << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                        << " " << node_indeces[2] + 1 << " "
                        << node_indeces[3] + 1 << " " << node_indeces[4] + 1
                        << " " << node_indeces[5] + 1 << " "
                        << node_indeces[6] + 1 << " " << node_indeces[7] + 1
                        << "\n";
                }
                break;
            //..................................................................
            case MshElemType::TRIANGLE:
                msh_file << "ET = TRIANGLE"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    msh_file << node_indeces[0] + 1 << " "
                             << node_indeces[1] + 1 << " "
                             << node_indeces[2] + 1 << "\n";
                }
                break;
            //..................................................................
            case MshElemType::TETRAHEDRON:
                msh_file << "ET = TETRAHEDRON"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    msh_file << node_indeces[0] + 1 << " "
                             << node_indeces[1] + 1 << " "
                             << node_indeces[2] + 1 << " "
                             << node_indeces[3] + 1 << "\n";
                }
                break;
            //..................................................................
            case MshElemType::PRISM:
                msh_file << "ET = BRICK"
                         << "\n";
                for (i = 0; i < no_nodes; i++)
                {
                    double const* const pnt_i(m_msh->nod_vector[i]->getData());
                    msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2]
                             << "\n";
                }
                for (i = 0; i < no_elements; i++)
                {
                    m_ele = m_msh->ele_vector[i];
                    m_ele->GetNodeIndeces(node_indeces);
                    if (m_ele->GetElementType() == MshElemType::PRISM)
                        msh_file
                            << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                            << " " << node_indeces[2] + 1 << " "
                            << node_indeces[2] + 1 << " " << node_indeces[3] + 1
                            << " " << node_indeces[4] + 1 << " "
                            << node_indeces[5] + 1 << " " << node_indeces[5] + 1
                            << "\n";
                    if (m_ele->GetElementType() == MshElemType::HEXAHEDRON)
                        msh_file
                            << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                            << " " << node_indeces[2] + 1 << " "
                            << node_indeces[3] + 1 << " " << node_indeces[4] + 1
                            << " " << node_indeces[5] + 1 << " "
                            << node_indeces[6] + 1 << " " << node_indeces[7] + 1
                            << "\n";
                }
                break;
            default:
                std::cerr << "MSHWriteTecplot MshElemType not handled"
                          << "\n";
        }
    }
}

/**************************************************************************
   MSHLib-Method:
   12/2005 OK Implementation
   07/2007 OK PCS
**************************************************************************/
MeshLib::CFEMesh* MSHGet(const std::string& geo_name)
{
    MeshLib::CFEMesh* m_msh = NULL;
    for (int i = 0; i < (int)fem_msh_vector.size(); i++)
    {
        m_msh = fem_msh_vector[i];
        if (m_msh->geo_name.compare(geo_name) == 0)
            return m_msh;
        if (m_msh->pcs_name.compare(geo_name) == 0)
            return m_msh;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
MeshLib::CFEMesh* MSHGet(const std::string& pcs_type_name,
                         const std::string& geo_name)
{
    MeshLib::CFEMesh* m_msh = NULL;
    for (int i = 0; i < (int)fem_msh_vector.size(); i++)
    {
        m_msh = fem_msh_vector[i];
        if ((m_msh->pcs_name.compare(pcs_type_name) == 0) &&
            (m_msh->geo_name.compare(geo_name) == 0))
            return m_msh;
    }
    return NULL;
}

/**************************************************************************
   PCSLib-Method:
   12/2005 OK Implementation
**************************************************************************/
MeshLib::CFEMesh* MSHGetGEO(std::string geo_name)
{
    int no_msh = (int)fem_msh_vector.size();
    // If there is only one msh file available, use it for all process. WW
    if (no_msh == 1)
        return fem_msh_vector[0];  // WW
    //----------------------------------------------------------------------
    MeshLib::CFEMesh* m_msh = NULL;
    for (int i = 0; i < no_msh; i++)
    {
        m_msh = fem_msh_vector[i];
        if (m_msh->geo_name.compare(geo_name) == 0)
            return m_msh;
    }
    //----------------------------------------------------------------------
    return NULL;
}

/**************************************************************************
   MSHLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CompleteMesh(std::string pcs_name)
{
    bool succeed = false;
    for (int i = 0; i < (int)fem_msh_vector.size(); i++)
        if (fem_msh_vector[i]->pcs_name.compare(pcs_name) == 0)
        {
            fem_msh_vector[i]->ConstructGrid();
            fem_msh_vector[i]->FillTransformMatrix();
            succeed = true;
        }
    return succeed;
}

/**************************************************************************/
/* ROCKFLOW - Function: MSHGetNextNode
 */
/* Task:
   Find the next node to the starting node in user defined direction
 */
/* Parameter: (I: Input; R: Return; X: Both)
   I: startnode, direction
 */
/* Return:
   nextnode
 */
/* Programming:
   09/2002     MB         First Version
   08/2005     MB
 */
/**************************************************************************/
long MSHGetNextNode(long startnode, MeshLib::CFEMesh* m_msh)
{
    size_t NumberOfNodes(m_msh->nod_vector.size());
    long NumberOfNodesPerLayer =
        NumberOfNodes / (m_msh->getNumberOfMeshLayers() + 1);
    return startnode + NumberOfNodesPerLayer;
}

/**************************************************************************/
/* ROCKFLOW - Function: MSHSelectFreeSurfaceNodes
 */
/* Task:
   Selection of free surface nodes, i.e. Setting free surface node flag = 1
   for the uppermost row and free surface node flag = 2 for the lowermost
   row. (Lowermost row is not moving)
 */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
 */
/* Programming:
   03/2003     MB   First Version
   09/2004     MB   PCS
   08/2005     MB msh
 */
/**************************************************************************/
void MSHSelectFreeSurfaceNodes(MeshLib::CFEMesh* m_msh)
{
    // Number of nodes per node layer
    size_t NumberOfNodesPerLayer =
        m_msh->nod_vector.size() / (m_msh->getNumberOfMeshLayers() + 1);
    size_t no_unconfined_layer = 0;
    // create array with nodes in vertical column
    size_t* strang(new size_t[m_msh->getNumberOfMeshLayers() * 2]);
    for (size_t i = 0; i < NumberOfNodesPerLayer; i++)
    {
        if (m_msh->nod_vector[i]->free_surface == 4)
        {
            size_t nextnode = i;
            no_unconfined_layer = 0;
            for (size_t j = 0; j < m_msh->getNumberOfMeshLayers(); j++)
            {
                //				strang = (long*)
                // Realloc(strang,(j+1)*sizeof(long));
                strang[j] = nextnode;
                size_t startnode = nextnode;
                nextnode = MSHGetNextNode(startnode, m_msh);
                if (m_msh->nod_vector[nextnode]->free_surface == 4)
                {
                    strang[j + 1] = nextnode;
                    no_unconfined_layer++;
                }
                else
                    continue;
            }
        }  // endif free_surface==4

        // mark start of vertical column with 1 - end of column with 2
        // this is than used in MSHMoveNODUcFlow
        m_msh->nod_vector[strang[0]]->free_surface = 1;
        m_msh->nod_vector[strang[no_unconfined_layer]]->free_surface = 2;
    } /*endfor*/
    delete[] strang;
}

/**************************************************************************
   FEMLib-Method:
   Task: Searches mobile nodes and sets node->free_surface = 4
   Programing:
   09/2004 OK / MB Implementation
   05/2005 OK Bugfix
   07/2005 MB MMP keyword
   08/2005 MB m_pcs
**************************************************************************/
void MSHDefineMobile(CRFProcess* m_pcs)
{
    long* mobile_nodes = NULL;
    long no_mobile_nodes = -1;
    //----------------------------------------------------------------------
    // Define mobile MSH nodes
    //----------------------------------------------------------------------
    // MMP Groups
    if (mmp_vector.size() == 0)
        return;
    ////Schleife �ber alle Gruppen
    for (size_t i = 0; i < mmp_vector.size(); i++)
    {
        CMediumProperties const* const m_mat_mp(mmp_vector[i]);

        // WW    int test = m_pcs->m_msh->GetMaxElementDim();
        // m_pcs->m_msh->cross_section

        // if (m_mat_mp->unconfined_flow_group ==1 &&
        // m_pcs->m_msh->GetMaxElementDim() == 3){
        if ((m_mat_mp->unconfined_flow_group == 1 &&
             m_pcs->m_msh->GetMaxElementDim() == 3) ||
            m_pcs->m_msh->hasCrossSection())
        {
            // if (m_mat_mp->unconfined_flow_group ==1){
            // if (m_pcs->m_msh->cross_section){
            //....................................................................
            // DOMAIN
            //			if(m_mat_mp->geo_type_name.find("DOMAIN") !=
            // std::string::npos)
            if (m_mat_mp->getGeoType() == GEOLIB::GEODOMAIN)
            {
                // CGLDomain *m_domain = NULL;
                // m_domain = m_domain->Get(m_mat_mp->geo_name);
                // mobile_nodes = m_domain->GetPointsIn(&no_mobile_nodes);
                // ToDo einlesen von domains ????
                for (size_t nodes = 0; nodes < m_pcs->m_msh->nod_vector.size();
                     nodes++)
                {
                    mobile_nodes = (long*)Realloc(mobile_nodes,
                                                  sizeof(long) * (nodes + 1));
                    mobile_nodes[nodes] = nodes;
                }
                i = m_pcs->m_msh->nod_vector.size();
                no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
            }
        }  // end if unconfined flow group
    }      // end for mmp vector

    // Set mobile MSH nodes flag
    for (long i = 0; i < no_mobile_nodes; i++)
        m_pcs->m_msh->nod_vector[i]->free_surface = 4;

    if (no_mobile_nodes > 0)
    {
        m_pcs->mobile_nodes_flag = 1;
        MSHSelectFreeSurfaceNodes(m_pcs->m_msh);
    }
}

/**************************************************************************
   ROCKFLOW - Function: MSHGetNodesInColumn

   Task:
   Gets nodes of a column searching downward from startnode.

   Parameter: (I: Input; R: Return; X: Both)
           I: long node, int anz_zeilen

   Return:
  *long strang

   Programming:
   09/2002   MB   First Version
   08/2005   MB   m_msh
**************************************************************************/
long* MSHGetNodesInColumn(long nextnode, int anz_zeilen,
                          MeshLib::CFEMesh* m_msh)
{
    int i;
    long startnode;
    long* strang = NULL;

    for (i = 0; i < anz_zeilen + 1; i++)
    {
        strang = (long*)Realloc(strang, (i + 1) * sizeof(long));
        strang[i] = nextnode;
        startnode = nextnode;
        // nextnode = MSHGetNextNode (startnode, direction);
        nextnode = MSHGetNextNode(startnode, m_msh);
    }
    return strang;
}

/**************************************************************************/
/* ROCKFLOW - Function: MSHMoveNODUcFlow
 */
/* Task:
   Moves free surface nodes according to the pressure distribution
 */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
 */
/* Programming:
   09/2002     MB       First Version
   05/2003     MB       verallgemeinert f�r Prismen und Vierecke
   09/2004     MB       Methode vereinfacht
   09/2004     MB       PCS
   08/2005      MB       m_msh */
/**************************************************************************/
void MSHMoveNODUcFlow(CRFProcess* m_pcs)
{
    long nextnode = -1;
    long startnode;
    int anz_zeilen = 0;
    int i;
    double spanne_ges;
    double spanne_rel;
    long* strang = NULL;
    double head = 0.0;
    int xxflag;
    int nidy;
    // Number of nodes per node layer
    const size_t NumberOfNodesPerLayer(
        m_pcs->m_msh->nod_vector.size() /
        (m_pcs->m_msh->getNumberOfMeshLayers() + 1));
    double MinThickness = 1e-1;  // OKMB
    double z_bottom;             // OKMB

    for (size_t node = 0; node < NumberOfNodesPerLayer; node++)

        if (m_pcs->m_msh->nod_vector[node]->free_surface == 1)
        {
            /* Z�hlen der Zeilen (-> anz_zeilen) */
            anz_zeilen = 0;
            xxflag = 0;
            nextnode = node;
            do
            {
                startnode = nextnode;
                nextnode = MSHGetNextNode(startnode, m_pcs->m_msh);

                /* Test2: Geh�rt der n�chste Knoten zu unterer Reihe ==> Abbruch
                 */
                if (m_pcs->m_msh->nod_vector[nextnode]->free_surface == 2)
                    xxflag = 1;
                anz_zeilen++; /* Anzahl der beweglichen Zeilen (ohne die feste
                                 untere Zeile) */
            } while (xxflag != 1);
            /** Ende Z�hlen der Zeilen zwischen den oberen free surface node
             * etc... und den Unteren **/

            /* Die Knoten unterhalb eines Free Surface Knotens bilden einen
             * Strang */
            /* Die Knoten eines Stranges werden zwischengespeichert */
            strang = MSHGetNodesInColumn(node, anz_zeilen, m_pcs->m_msh);

            /* Die Knoten eines Stranges werden entsprechend der neuen
             * Druckverteilung  verformt */
            /* Standrohrspiegelh�he bestimmen */
            nidy = m_pcs->GetNodeValueIndex("HEAD") + 1;
            if (GetRFProcessDensityFlow()) /* mit Dichteunterschiede */
            {
                // OK_MOD     head = MODCalcHeadInColumn_MB(strang, anz_zeilen);
            }
            else /* ohne Dichteunterschiede */
                head = m_pcs->GetNodeValue(strang[0], nidy);

            /* nicht �ber surface elevation */
            CRFProcess* m_pcs_OLF = NULL;
            m_pcs_OLF = PCSGet("OVERLAND_FLOW");
            double SurfaceZ;

            if (m_pcs_OLF != NULL)
            {
                SurfaceZ =
                    m_pcs_OLF->m_msh->nod_vector[strang[0]]->getData()[2];
                if (head > SurfaceZ)
                    head = SurfaceZ;
            }

            /* Set minimum thickness */
            z_bottom =
                m_pcs->m_msh->nod_vector[strang[anz_zeilen]]->getData()[2];
            if (head - z_bottom < MinThickness)
                head = z_bottom + MinThickness;

            /* Berechnung der Differenz */
            spanne_ges = head - z_bottom;
            spanne_rel = spanne_ges / anz_zeilen;
            m_pcs->m_msh->nod_vector[strang[0]]->SetZ(head);

            if (spanne_ges != 0)
                /* Setzen der neuen Z-Werte entlang eines Stranges */
                for (i = 1; i < anz_zeilen;
                     i++) /* Schleife �ber Anzahl der Zeilen */
                    m_pcs->m_msh->nod_vector[strang[i]]->SetZ(head -
                                                              i * spanne_rel);

            strang = (long*)Free(strang);
        } /*endif index ==1 */
          /* end for Schleife �ber alle Knoten */
}

/**************************************************************************
   FEMLib-Method:
   Task: Searches mobile nodes and sets node->free_surface = 4
   Programing:
   09/2004 OK / MB Implementation
   05/2005 OK Bugfix
   07/2005 MB MMP keyword
   08/2005 MB m_pcs
   01/2006 OK LAYER
**************************************************************************/
void CFEMesh::DefineMobileNodes(CRFProcess* m_pcs)
{
    long* mobile_nodes = NULL;
    long no_mobile_nodes = -1;
    long i, j;
    //----------------------------------------------------------------------
    // Define mobile MSH nodes
    //----------------------------------------------------------------------
    //......................................................................
    // DOMAIN
    if (m_pcs->geo_type.find("DOMAIN") != std::string::npos)
    {
        for (i = 0; i < (long)nod_vector.size(); i++)
        {
            mobile_nodes = (long*)Realloc(mobile_nodes, sizeof(long) * (i + 1));
            mobile_nodes[i] = i;
        }
        no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
    }
    //......................................................................
    // LAYER
    if (m_pcs->geo_type.find("LAYER") != std::string::npos)
    {
        std::string m_string;
        long no_nodes_per_layer =
            (long)nod_vector.size() / (getNumberOfMeshLayers() + 1);
        int pos = 0;
        int layer_start = 0, layer_end = 0;
        if (m_pcs->geo_type_name.find("-") != std::string::npos)
        {
            pos = m_pcs->geo_type_name.find("-") != std::string::npos;
            m_string = m_pcs->geo_type_name.substr(0, pos);
            layer_start = strtol(m_string.c_str(), NULL, 0);
            m_string = m_pcs->geo_type_name.substr(pos + 1, std::string::npos);
            layer_end = strtol(m_string.c_str(), NULL, 0);
        }
        else
        {
            layer_start = strtol(m_pcs->geo_type_name.c_str(), NULL, 0);
            layer_end = layer_start;
        }
        int no_layers = layer_end - layer_start + 1;
        no_mobile_nodes = (no_layers + 1) * no_nodes_per_layer;
        mobile_nodes = new long[no_mobile_nodes];
        for (i = 0; i < no_layers + 1; i++)
            for (j = 0; j < no_nodes_per_layer; j++)
                mobile_nodes[i * no_nodes_per_layer + j] =
                    j + (layer_start - 1 + i) * no_nodes_per_layer;
    }
    //......................................................................
    //----------------------------------------------------------------------
    // Set mobile MSH nodes flag
    //----------------------------------------------------------------------
    for (i = 0; i < (long)nod_vector.size(); i++)
        nod_vector[i]->free_surface = -1;
    for (i = 0; i < no_mobile_nodes; i++)
        nod_vector[i]->free_surface = 4;
    // nod_vector[mobile_nodes[i]]->free_surface = 4;
    //----------------------------------------------------------------------
    if (no_mobile_nodes > 0)
    {
        m_pcs->mobile_nodes_flag = 1;
        MSHSelectFreeSurfaceNodes(this);
    }
    //----------------------------------------------------------------------
    delete[] mobile_nodes;
    mobile_nodes = NULL;
}

/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: Get element nodes in a material domain
   Programing:
   10/2004 WW Implementation
   06/2012 NW Made this function faster using std::set for large data set
**************************************************************************/
void GEOGetNodesInMaterialDomain(CFEMesh const* const msh, int MatIndex,
                                 std::vector<long>& Nodes, bool Order)
{
    Nodes.resize(0);
    std::set<long> set_nodes;

    std::vector<bool> node_status(msh->nod_vector.size());

    const size_t n_ele(msh->ele_vector.size());
    for (size_t e = 0; e < n_ele; e++)
    {
        MeshLib::CElem const* const elem(msh->ele_vector[e]);
        // if (elem->GetMark())              // Marked for use
        {
            const int nn(elem->GetNodesNumber(Order));
            if (elem->GetPatchIndex() == static_cast<size_t>(MatIndex))
            {
                for (int i = 0; i < nn; i++)
                {
                    const long node_id = elem->GetNodeIndex(i);
                    if(node_status[node_id])
                        continue;

                    set_nodes.insert(node_id);
                    node_status[node_id] = true;
                }
            }
        }  // if
    }      // For
    Nodes.assign(set_nodes.begin(), set_nodes.end());
}

/**************************************************************************
   MSHLib-Method:
   01/2006 OK Implementation
   06/2009 OK Bug fix
**************************************************************************/
int MSHSetMaxMMPGroups()
{
    int i;
    long j;
    CFEMesh* m_msh = NULL;
    //----------------------------------------------------------------------
    size_t msh_max_mmp_groups;
    for (i = 0; i < (int)fem_msh_vector.size(); i++)
    {
        m_msh = fem_msh_vector[i];
        m_msh->max_mmp_groups = 0;
        msh_max_mmp_groups = 0;
        for (j = 0; j < (long)m_msh->ele_vector.size(); j++)
            if ((m_msh->ele_vector[j]->GetPatchIndex() + 1) >
                msh_max_mmp_groups)
                msh_max_mmp_groups++;
        m_msh->max_mmp_groups = msh_max_mmp_groups;
    }
    //----------------------------------------------------------------------
    size_t g_msh_max_mmp_groups = 0;
    for (i = 0; i < (int)fem_msh_vector.size(); i++)
        if (m_msh->max_mmp_groups > g_msh_max_mmp_groups)
            g_msh_max_mmp_groups++;
    //----------------------------------------------------------------------
    return g_msh_max_mmp_groups;
}

/**************************************************************************
   MSHLib-Method:
   12/2014 NW Implement
**************************************************************************/
size_t MSHGetMaxPatchIndex(const CFEMesh* m_msh)
{
    size_t max_mat_id = 0;
    for (size_t j = 0; j < m_msh->ele_vector.size(); j++)
        max_mat_id =
            std::max(max_mat_id, m_msh->ele_vector[j]->GetPatchIndex());
    return max_mat_id;
}

/**************************************************************************
   MSHLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool MSHTestMATGroups()
{
    int g_max_mmp_groups = MSHSetMaxMMPGroups();
    if (g_max_mmp_groups > (int)mmp_vector.size())
    {
        std::cout << "Error: not enough MMP data";
        return false;  // abort();
    }
    return true;
}
