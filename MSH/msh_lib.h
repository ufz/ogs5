/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#ifndef msh_lib_INC
#define msh_lib_INC

#include "msh_mesh.h"

extern std::vector<MeshLib::CFEMesh*> fem_msh_vector;
extern MeshLib::CFEMesh* FEMGet(const std::string& msh_name);

/**
 * reads a mesh. The following formats are possible:
 * <ol>
 * <li>native OGS meshes</li>
 * <li>TetGen</li>
 * <li>GMSH</li>
 * </ol>
 * @param mesh_fname file name of the mesh
 * @param mesh_vec a vector, the new mesh will be put in this vector
 * @param geo_obj object, that manages the geometric entities
 * @param unique_name the name of geometric data
 */
void FEMRead(const std::string& mesh_fname,
             std::vector<MeshLib::CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj = NULL,
             std::string* unique_name = NULL);
#if defined(USE_PETSC)  // || defined(using other parallel scheme)
void FEMRead_ASCII(const int msize, const int mrank,
                   const std::string& file_base_name,
                   std::vector<MeshLib::CFEMesh*>& mesh_vec,
                   GEOLIB::GEOObjects* geo_obj = NULL,
                   std::string* unique_name = NULL);
void FEMRead_BIN(const int msize, const int mrank,
                 const std::string& file_base_name,
                 std::vector<MeshLib::CFEMesh*>& mesh_vec,
                 GEOLIB::GEOObjects* geo_obj = NULL,
                 std::string* unique_name = NULL);
#endif

extern void MSHWrite(std::string);
extern void CompleteMesh();  // WW
extern void FEMDeleteAll();
void Read_RFI(std::istream& msh_file, MeshLib::CFEMesh* m_msh);
extern int MSHSetMaxMMPGroups();  // OK
extern bool MSHTestMATGroups();   // OK
extern size_t MSHGetMaxPatchIndex(const MeshLib::CFEMesh* m_msh);
extern void GEOGetNodesInMaterialDomain(MeshLib::CFEMesh const* const, int,
                                        std::vector<long>&, bool);

extern void MSHDefineMobile(CRFProcess*);   // OK411
extern void MSHMoveNODUcFlow(CRFProcess*);  // OK411
extern MeshLib::CFEMesh* MSHGet(const std::string& mat_type_name);
extern void DATWriteParticleFile(int);  // PCH

#endif
