/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
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
extern MeshLib::CFEMesh* FEMGet(const std::string &msh_name);
//OK
extern void MSHCreateNOD2ELERelations(MeshLib::CFEMesh*);

//MeshLib::CFEMesh* FEMRead(const std::string& fname,
//                          GEOLIB::GEOObjects* geo_obj = NULL,
//                          std::string* unique_name = NULL);

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
#if defined(USE_PETSC) // || defined(using other parallel scheme)
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
extern void CompleteMesh();                       //WW
extern bool CompleteMesh(std::string);            //OK
extern void FEMDeleteAll();
//KR extern void MSHCalcMinMaxMidCoordinates();        //OK
//KR extern double msh_x_min,msh_x_max;                //OK
//KR extern double msh_y_min,msh_y_max;                //OK
//KR extern double msh_z_min,msh_z_max;                //OK
//KR extern double msh_x_mid,msh_y_mid,msh_z_mid;      //OK
// Might be removed
void Read_RFI(std::istream& msh_file, MeshLib::CFEMesh* m_msh);
extern void GMSH2MSH(const char*, MeshLib::CFEMesh*);
extern int MSHSetMaxMMPGroups();                  //OK
extern bool MSHTestMATGroups();                   //OK
extern size_t MSHGetMaxPatchIndex(const MeshLib::CFEMesh* m_msh);
extern void GEOGetNodesInMaterialDomain(MeshLib::CFEMesh const*const, int, std::vector<long>&, bool);

extern void MSHDefineMobile(CRFProcess*);         //OK411
extern void MSHMoveNODUcFlow (CRFProcess*);       //OK411
extern MeshLib::CFEMesh* MSHGet(const std::string &mat_type_name);
extern void DATWriteParticleFile(int);            // PCH
extern long* GetPointsIn(Surface*,long*);         //OK411
extern long* MSHGetNodesClose(long*,CGLPolyline*); //OK411
                                                  //OK411
#ifdef ObsoleteGUI //WW 03.2012
extern void MSHWriteVOL2TEC(std::string);         //OK
extern void MSHLayerWriteTecplot();               //OK
#endif
#ifdef Obsolete //WW 03.2012
extern void MSHAssignMATGroup2Elements(std::string);
extern void MSHCreateQuadsFromPLY(CGLPolyline*,int);
//OK411 extern void MSHCreatePrismsFromTriangles();
extern void MSHCreateNodes();
extern void MSHDeleteDoubleElements(int);
extern long MSHMarkDoubleElementsType(int);
extern void RFIWriteTecplot();
extern void MSHWriteTecplot();
extern void MSHAssignMATGroup2LineElements();
extern void MSHAssignMATGroup2TrisElements(std::string);
extern void MSHAssignMATGroup2QuadElements();
extern void MSHAssignMATGroup2TetsElements();
extern void MSHAssignMATGroup2PrisElements();
extern void MSHAssignMATGroup2PrisElementsNew();
extern void MSH2MATPris();
extern void MSHAssignMATGroup2HexsElements();
extern void MSHDestroy();
extern void MSHDelete(std::string);
extern void DATWriteRFIFile(const char* file_name);
//KR extern bool msh_file_binary; //OK
//TK
extern void Mesh_Single_Surface(std::string surface_name, const char* file_name_const_char);
//extern void Select_Nodes_Elements_by_TINFile(const char *file_name_const_char);
extern void Clear_Selected_Nodes_Elements();
extern void GMSH2TIN(const char* file_name_const_char);

//OK
extern MeshLib::CFEMesh* MSHGet(const std::string &mat_type_name);
//OK
extern MeshLib::CFEMesh* MSHGet(const std::string &pcs_type_name,const std::string &mat_type_name);
extern MeshLib::CFEMesh* MSHGetGEO(std::string);           //OK

//extern bool IsPointInSurface(Surface*,CGLPoint*); //OK411

extern void SetRFIPointsClose(CGLLine*);          //OK411
                                                  //OK411
extern void MSHGetNodesClose(std::vector<long>&,CGLPoint*);
#endif //#ifdef Obsolete //WW 03.2012

#endif
