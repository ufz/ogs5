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
   08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
// C++
#include <string>
#include <vector>

#include "memory.h"

// GEOLib
#include "files0.h"
#include "geo_lib.h"
// MSHLib
#include "msh_lib.h"
// PCSLib
#include "mathlib.h"
#include "rf_mmp_new.h" //OK411

// WW extern void RFConfigRenumber(void);
#ifndef NEW_EQS // WW. 07.11.2008
extern void ConfigRenumberProperties(void);
#endif
extern int ReadRFIFile(std::string g_strFileNameBase);
#include "rf_pcs.h"

std::vector<MeshLib::CFEMesh*> fem_msh_vector;

#define FEM_FILE_EXTENSION ".msh"

// KR double msh_x_min,msh_x_max;                       //OK
// KR double msh_y_min,msh_y_max;                       //OK
// KR double msh_z_min,msh_z_max;                       //OK
// KR double msh_x_mid,msh_y_mid,msh_z_mid;             //OK

//#define MSH_SIZE 1e5

using namespace Math_Group;

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
**************************************************************************/
void MSHDelete(std::string m_msh_name)
{
	CFEMesh* m_fem_msh = NULL;
	size_t fem_msh_vector_size = fem_msh_vector.size();
	for (size_t i = 0; i < fem_msh_vector_size; i++)
	{
		m_fem_msh = fem_msh_vector[i];
		if (m_fem_msh->pcs_name.compare(m_msh_name) == 0)
		{
			delete m_fem_msh;
			fem_msh_vector.erase((fem_msh_vector.begin() + i));
		}
	}
}

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

#ifndef USE_PETSC // && not defined(other parallel method with ddc)
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
void FEMRead(const std::string& file_base_name, std::vector<CFEMesh*>& mesh_vec, GEOLIB::GEOObjects* geo_obj,
             std::string* unique_name)
{
	CFEMesh* mesh(NULL);
	std::string msh_file_name(file_base_name + FEM_FILE_EXTENSION);

	std::ifstream msh_file_ascii(msh_file_name.data(), std::ios::in);
	if (!msh_file_ascii.is_open())
		std::cout << "CFEMesh::FEMRead() - Could not open file...\n";

	std::cout << "MSHRead:  ASCII file"
	          << "\n";
	std::string line_string("");
	getline(msh_file_ascii, line_string);

	bool more_mesh = false; // 12.08.2011. WW
	if (line_string.find("#FEM_MSH") != std::string::npos) // OGS mesh file
	{
		mesh = new CFEMesh(geo_obj, unique_name);
		more_mesh = mesh->Read(&msh_file_ascii);
		mesh_vec.push_back(mesh); // TF

		// Multi-mesh 12.08.2011 WW
		if (more_mesh)
			while (!msh_file_ascii.eof())
			{
				// getline(msh_file_ascii, line_string);
				// if(line_string.find("#FEM_MSH")!=std::string::npos)
				mesh = new CFEMesh(geo_obj, unique_name);
				more_mesh = mesh->Read(&msh_file_ascii);
				mesh_vec.push_back(mesh); // TF
				if (!more_mesh)
					break;
			}
		//  if(line_string.find("#STOP")!=std::string::npos)
		//     break;
	}
	else // RFI mesh file
	{
		msh_file_ascii.seekg(0L, std::ios::beg);
		mesh = new CFEMesh(geo_obj, unique_name);
		Read_RFI(msh_file_ascii, mesh);
		mesh_vec.push_back(mesh); // 12.08.2011 WW
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
		getline(msh_file, strbuffer); // The first line
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
	std::ofstream fem_msh_file(fem_msh_file_name.c_str(), std::ios::trunc | std::ios::out);

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
		return fem_msh_vector[0]; // WW
	for (size_t i = 0; i < no_msh; i++)
		if (fem_msh_vector[i]->pcs_name.compare(msh_name) == 0)
			return fem_msh_vector[i];
	return NULL;
}

/* KR method not used

 **************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
 **************************************************************************
   void MSHCalcMinMaxMidCoordinates()
   {
   double m_dXMin1 = 1.e+19;
   double m_dXMax1 = -1.e+19;
   double m_dYMin1 = 1.e+19;
   double m_dYMax1 = -1.e+19;
   double m_dZMin1 = 1.e+19;
   double m_dZMax1 = -1.e+19;
   double value;
   CFEMesh* m_msh = NULL;
   //----------------------------------------------------------------------
   for(int j=0;j<(int)fem_msh_vector.size();j++)
   {
      m_msh = fem_msh_vector[j];
      for(long i=0;i<(long)m_msh->nod_vector.size();i++)
      {
         value = m_msh->nod_vector[i]->X();
         if(value<m_dXMin1) m_dXMin1 = value;
         if(value>m_dXMax1) m_dXMax1 = value;
         value = m_msh->nod_vector[i]->Y();
         if(value<m_dYMin1) m_dYMin1 = value;
         if(value>m_dYMax1) m_dYMax1 = value;
         value = m_msh->nod_vector[i]->Z();
         if(value<m_dZMin1) m_dZMin1 = value;
         if(value>m_dZMax1) m_dZMax1 = value;
         //..................................................................
         // Shrink a bit
         msh_x_min = m_dXMin1 - 0.05*(m_dXMax1-m_dXMin1);
         msh_x_max = m_dXMax1 + 0.05*(m_dXMax1-m_dXMin1);
         msh_y_min = m_dYMin1 - 0.05*(m_dYMax1-m_dYMin1);
         msh_y_max = m_dYMax1 + 0.05*(m_dYMax1-m_dYMin1);
         msh_z_min = m_dZMin1 - 0.05*(m_dZMax1-m_dZMin1);
         msh_z_max = m_dZMax1 + 0.05*(m_dZMax1-m_dZMin1);
      }
   }
   //----------------------------------------------------------------------
   msh_x_mid = 0.5*(msh_x_min+msh_x_max);
   msh_y_mid = 0.5*(msh_y_min+msh_y_max);
   msh_z_mid = 0.5*(msh_z_min+msh_z_max);
   //----------------------------------------------------------------------
   }
 */

#ifdef ObsoleteGUI // WW 03.2012
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   04/2004 OK Implementation
   01/2005 OK File handling
   09/2005 OK MSH ToDo
   last modification:
**************************************************************************/
void MSHWriteVOL2TEC(std::string m_msh_name)
{
	long i, j;
	CGLVolume* m_vol = NULL;
	std::vector<CGLVolume*>::const_iterator p_vol;
	std::string name("VOLUMES");
	std::vector<Surface*>::const_iterator p_sfc;
	double x, y, z;
	//  CGLPoint m_point; // TF
	std::ios::pos_type position;
	int vol_number = -1;
	Surface* m_sfc = NULL;
	//--------------------------------------------------------------------
	MeshLib::CFEMesh* m_msh(FEMGet(m_msh_name));
	if (!m_msh)
		return;
	long no_nodes = (long)m_msh->nod_vector.size();
	long ep_layer = (long)m_msh->ele_vector.size() / m_msh->getNumberOfMeshLayers();
	//--------------------------------------------------------------------
	// File handling
	std::string tec_path;
	CGSProject* m_gsp = GSPGetMember("gli");
	if (m_gsp)
		tec_path = m_gsp->path;
	//======================================================================
	p_vol = volume_vector.begin();
	while (p_vol != volume_vector.end())
	{
		m_vol = *p_vol;
		if (m_vol->layer == 0) // OK
		{
			p_vol++;
			continue;
		}
		p_sfc = m_vol->surface_vector.begin();
		m_sfc = *p_sfc;
		if (!m_sfc)
			return;
		//--------------------------------------------------------------------
		long jb = (m_vol->layer - 1) * ep_layer;
		long je = jb + ep_layer;
		vol_number++;
		//--------------------------------------------------------------------
		// file handling
		std::string vol_file_name = tec_path + "VOL_" + m_vol->name + TEC_FILE_EXTENSION;
		std::fstream vol_file(vol_file_name.data(), std::ios::trunc | std::ios::out);
		vol_file.setf(std::ios::scientific, std::ios::floatfield);
		vol_file.precision(12);
		if (!vol_file.good())
			return;
		vol_file.seekg(0L, std::ios::beg);
		//--------------------------------------------------------------------
		vol_file << "VARIABLES = X,Y,Z,VOL"
		         << "\n";
		//--------------------------------------------------------------------
		long no_mat_elements = 0;
		MeshLib::CElem* m_ele = NULL;
		vec<long> node_indeces(6);
		for (i = jb; i < je; i++)
		{
			m_ele = m_msh->ele_vector[i];
			if (m_ele->GetElementType() == MshElemType::PRISM)
			{
				m_ele->GetNodeIndeces(node_indeces);
				// nodes = m_msh->ele_vector[i]->nodes;
				x = 0.0;
				y = 0.0;
				z = 0.0;
				for (j = 0; j < 6; j++)
				{
					double const* const pnt(m_msh->nod_vector[node_indeces[j]]->getData());
					x += pnt[0];
					y += pnt[1];
					z += pnt[2];
				}
				x /= double(6);
				y /= double(6);
				z /= double(6);
				// TF m_sfc->PointInSurface(&m_point) returns always false
				//        m_point.x = x;
				//        m_point.y = y;
				//        m_point.z = z;
				//        if(m_sfc->PointInSurface(&m_point)){
				//          no_mat_elements++;
				//        }
			}
		}
		//--------------------------------------------------------------------
		position = vol_file.tellg();
		vol_file << "ZONE T = " << m_vol->name << ", "
		         << "N = " << no_nodes << ", "
		         << "E = " << no_mat_elements << ", "
		         << "F = FEPOINT"
		         << ", "
		         << "ET = BRICK"
		         << "\n";
		for (i = 0; i < no_nodes; i++)
		{
			double const* const pnt_i(m_msh->nod_vector[i]->getData());
			vol_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << " " << vol_number << "\n";
		}
		for (long i = jb; i < je; i++)
		{
			m_ele = m_msh->ele_vector[i];
			if (m_ele->GetElementType() == MshElemType::PRISM)
			{
				m_ele->GetNodeIndeces(node_indeces);
				x = 0.0;
				y = 0.0;
				z = 0.0;
				for (j = 0; j < 6; j++)
				{
					double const* const pnt_j(m_msh->nod_vector[node_indeces[j]]->getData());
					x += pnt_j[0];
					y += pnt_j[1];
					z += pnt_j[2];
				}
				x /= double(6);
				y /= double(6);
				z /= double(6);
				// TF m_sfc->PointInSurface(&m_point) returns always false
				//        m_point.x = x;
				//        m_point.y = y;
				//        m_point.z = z;
				//        if(m_sfc->PointInSurface(&m_point)){
				//          vol_file
				//            << node_indeces[0]+1 << " " << node_indeces[0]+1 << " " << node_indeces[1]+1 << " " <<
				//            node_indeces[2]+1 << " "
				//            << node_indeces[3]+1 << " " << node_indeces[3]+1 << " " << node_indeces[4]+1 << " " <<
				//            node_indeces[5]+1 << "\n";
				//        }
			}
		}
		++p_vol;
		//======================================================================
	}
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
   11/2005 OK OO-ELE
**************************************************************************/
void MSHLayerWriteTecplot()
{
	MshElemType::type ele_type = MshElemType::INVALID;
	// WW long no_nodes;
	long no_elements;
	std::string delimiter(", ");
	MeshLib::CElem* m_ele = NULL;
	vec<long> node_indeces(8);
	std::string no_layer_str;
	char no_layer_char[3];
	//----------------------------------------------------------------------
	// File handling
	std::string file_path = "MSH";
	CGSProject* m_gsp = NULL;
	m_gsp = GSPGetMember("msh");
	if (m_gsp)
		file_path = m_gsp->path;
	//----------------------------------------------------------------------
	MeshLib::CFEMesh* m_msh = NULL;
	for (int j = 0; j < (int)fem_msh_vector.size(); j++)
	{
		m_msh = fem_msh_vector[j];
		for (size_t k = 0; k < m_msh->getNumberOfMeshLayers(); k++)
		{
			sprintf(no_layer_char, "%lu", static_cast<long unsigned>(k) + 1);
			no_layer_str = no_layer_char;
			// WW no_nodes = (long) m_msh->nod_vector.size() / (m_msh->getNumberOfMeshLayers() + 1);
			no_elements = (long)m_msh->ele_vector.size() / m_msh->getNumberOfMeshLayers();
			// Test ele_type
			if (no_elements > 0)
			{
				m_ele = m_msh->ele_vector[0];
				ele_type = m_ele->GetElementType();
			}
			// File handling
			std::string msh_file_name
			    = file_path + "MSH_LAYER" + no_layer_str + "_" + m_msh->pcs_name + TEC_FILE_EXTENSION;
			std::fstream msh_file(msh_file_name.data(), std::ios::trunc | std::ios::out);
			msh_file.setf(std::ios::scientific, std::ios::floatfield);
			msh_file.precision(12);
			if (!msh_file.good())
				return;
			msh_file.seekg(0L, std::ios::beg);
			msh_file << "VARIABLES = X,Y,Z"
			         << "\n";
			msh_file << "ZONE T = " << m_msh->pcs_name << delimiter << "N = " << (long)m_msh->nod_vector.size()
			         << delimiter << "E = " << no_elements << delimiter;
			msh_file << "F = FEPOINT" << delimiter;
			switch (ele_type)
			{
				//..................................................................
				case MshElemType::LINE:
					msh_file << "ET = QUADRILATERAL"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[1] + 1
						         << " " << node_indeces[0] + 1 << "\n";
					}
					break;
				//..................................................................
				case MshElemType::QUAD:
					msh_file << "ET = QUADRILATERAL"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[3] + 1 << "\n";
					}
					break;
				//..................................................................
				case MshElemType::HEXAHEDRON:
					msh_file << "ET = BRICK"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[3] + 1 << " " << node_indeces[4] + 1 << " "
						         << node_indeces[5] + 1 << " " << node_indeces[6] + 1 << " " << node_indeces[7] + 1
						         << "\n";
					}
					break;
				//..................................................................
				case MshElemType::TRIANGLE:
					msh_file << "ET = TRIANGLE"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << "\n";
					}
					break;
				//..................................................................
				case MshElemType::TETRAHEDRON:
					msh_file << "ET = TETRAHEDRON"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[3] + 1 << "\n";
					}
					break;
				//..................................................................
				case MshElemType::PRISM:
					msh_file << "ET = BRICK"
					         << "\n";
					for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
					{
						double const* const pnt_i(m_msh->nod_vector[i]->getData());
						msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
					}
					for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
					{
						m_ele = m_msh->ele_vector[i];
						m_ele->GetNodeIndeces(node_indeces);
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[2] + 1 << " " << node_indeces[3] + 1 << " "
						         << node_indeces[4] + 1 << " " << node_indeces[5] + 1 << " " << node_indeces[5] + 1
						         << "\n";
					}
					break;
				default:
					std::cerr << "MSHLayerWriteTecplot MshElemType not handled"
					          << "\n";
			}
		} // layer
	}
}
#endif //#ifdef ObsoleteGUI //WW 03.2012

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
		std::string msh_file_name = file_path + "_" + m_msh->pcs_name + TEC_FILE_EXTENSION;
		std::fstream msh_file(msh_file_name.data(), std::ios::trunc | std::ios::out);
		msh_file.setf(std::ios::scientific, std::ios::floatfield);
		msh_file.precision(12);
		if (!msh_file.good())
			return;
		msh_file.seekg(0L, std::ios::beg);
		msh_file << "VARIABLES = X,Y,Z"
		         << "\n";
		msh_file << "ZONE T = " << m_msh->pcs_name << delimiter << "N = " << no_nodes << delimiter
		         << "E = " << no_elements << delimiter;
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
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[1] + 1 << " "
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
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1 << " "
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
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1 << " "
					         << node_indeces[3] + 1 << " " << node_indeces[4] + 1 << " " << node_indeces[5] + 1 << " "
					         << node_indeces[6] + 1 << " " << node_indeces[7] + 1 << "\n";
				}
				break;
			//..................................................................
			case MshElemType::TRIANGLE:
				msh_file << "ET = TRIANGLE"
				         << "\n";
				for (i = 0; i < no_nodes; i++)
				{
					double const* const pnt_i(m_msh->nod_vector[i]->getData());
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1 << "\n";
				}
				break;
			//..................................................................
			case MshElemType::TETRAHEDRON:
				msh_file << "ET = TETRAHEDRON"
				         << "\n";
				for (i = 0; i < no_nodes; i++)
				{
					double const* const pnt_i(m_msh->nod_vector[i]->getData());
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1 << " "
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
					msh_file << pnt_i[0] << " " << pnt_i[1] << " " << pnt_i[2] << "\n";
				}
				for (i = 0; i < no_elements; i++)
				{
					m_ele = m_msh->ele_vector[i];
					m_ele->GetNodeIndeces(node_indeces);
					if (m_ele->GetElementType() == MshElemType::PRISM)
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[2] + 1 << " " << node_indeces[3] + 1 << " "
						         << node_indeces[4] + 1 << " " << node_indeces[5] + 1 << " " << node_indeces[5] + 1
						         << "\n";
					if (m_ele->GetElementType() == MshElemType::HEXAHEDRON)
						msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1 << " " << node_indeces[2] + 1
						         << " " << node_indeces[3] + 1 << " " << node_indeces[4] + 1 << " "
						         << node_indeces[5] + 1 << " " << node_indeces[6] + 1 << " " << node_indeces[7] + 1
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
MeshLib::CFEMesh* MSHGet(const std::string& pcs_type_name, const std::string& geo_name)
{
	MeshLib::CFEMesh* m_msh = NULL;
	for (int i = 0; i < (int)fem_msh_vector.size(); i++)
	{
		m_msh = fem_msh_vector[i];
		if ((m_msh->pcs_name.compare(pcs_type_name) == 0) && (m_msh->geo_name.compare(geo_name) == 0))
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
		return fem_msh_vector[0]; // WW
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
	long NumberOfNodesPerLayer = NumberOfNodes / (m_msh->getNumberOfMeshLayers() + 1);
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
	size_t NumberOfNodesPerLayer = m_msh->nod_vector.size() / (m_msh->getNumberOfMeshLayers() + 1);
	size_t no_unconfined_layer = 0;
	// create array with nodes in vertical column
	size_t* strang(new size_t[m_msh->getNumberOfMeshLayers()]);
	for (size_t i = 0; i < NumberOfNodesPerLayer; i++)
	{
		if (m_msh->nod_vector[i]->free_surface == 4)
		{
			size_t nextnode = i;
			no_unconfined_layer = 0;
			for (size_t j = 0; j < m_msh->getNumberOfMeshLayers(); j++)
			{
				//				strang = (long*) Realloc(strang,(j+1)*sizeof(long));
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
		} // endif free_surface==4

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

		// if (m_mat_mp->unconfined_flow_group ==1 && m_pcs->m_msh->GetMaxElementDim() == 3){
		if ((m_mat_mp->unconfined_flow_group == 1 && m_pcs->m_msh->GetMaxElementDim() == 3)
		    || m_pcs->m_msh->hasCrossSection())
		{
			// if (m_mat_mp->unconfined_flow_group ==1){
			// if (m_pcs->m_msh->cross_section){
			//....................................................................
			// DOMAIN
			//			if(m_mat_mp->geo_type_name.find("DOMAIN") != std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::GEODOMAIN)
			{
				// CGLDomain *m_domain = NULL;
				// m_domain = m_domain->Get(m_mat_mp->geo_name);
				// mobile_nodes = m_domain->GetPointsIn(&no_mobile_nodes);
				// ToDo einlesen von domains ????
				for (size_t nodes = 0; nodes < m_pcs->m_msh->nod_vector.size(); nodes++)
				{
					mobile_nodes = (long*)Realloc(mobile_nodes, sizeof(long) * (nodes + 1));
					mobile_nodes[nodes] = nodes;
				}
				i = m_pcs->m_msh->nod_vector.size();
				no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
			}

			// SURFACE
			//			if(m_mat_mp->geo_type_name.find("SURFACE") != std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::SURFACE)
			{
				Surface* m_surface = NULL;
				// CC
				m_surface = GEOGetSFCByName(m_mat_mp->geo_name);
				// CC
				mobile_nodes = GetPointsIn(m_surface, &no_mobile_nodes);
			}

			// VOLUME
			//			if(m_mat_mp->geo_type_name.find("VOLUME") != std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::VOLUME)
			{
				// WW CGLVolume *m_volume = NULL;
				// CC 10/05
				// WW  m_volume = GEOGetVOL(m_mat_mp->geo_name);
				// ToDo TK
				// OK411 mobile_nodes =GetPointsInVolume(m_volume,&no_mobile_nodes);//CC 10/05
			}
		} // end if unconfined flow group
	} // end for mmp vector

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
long* MSHGetNodesInColumn(long nextnode, int anz_zeilen, MeshLib::CFEMesh* m_msh)
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
   08/2005      MB       m_msh                                                                   */
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
	const size_t NumberOfNodesPerLayer(m_pcs->m_msh->nod_vector.size() / (m_pcs->m_msh->getNumberOfMeshLayers() + 1));
	double MinThickness = 1e-1; // OKMB
	double z_bottom; // OKMB

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

				/* Test2: Geh�rt der n�chste Knoten zu unterer Reihe ==> Abbruch */
				if (m_pcs->m_msh->nod_vector[nextnode]->free_surface == 2)
					xxflag = 1;
				anz_zeilen++; /* Anzahl der beweglichen Zeilen (ohne die feste untere Zeile) */
			} while (xxflag != 1);
			/** Ende Z�hlen der Zeilen zwischen den oberen free surface node etc... und den Unteren **/

			/* Die Knoten unterhalb eines Free Surface Knotens bilden einen Strang */
			/* Die Knoten eines Stranges werden zwischengespeichert */
			strang = MSHGetNodesInColumn(node, anz_zeilen, m_pcs->m_msh);

			/* Die Knoten eines Stranges werden entsprechend der neuen Druckverteilung  verformt */
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
				SurfaceZ = m_pcs_OLF->m_msh->nod_vector[strang[0]]->getData()[2];
				if (head > SurfaceZ)
					head = SurfaceZ;
			}

			/* Set minimum thickness */
			z_bottom = m_pcs->m_msh->nod_vector[strang[anz_zeilen]]->getData()[2];
			if (head - z_bottom < MinThickness)
				head = z_bottom + MinThickness;

			/* Berechnung der Differenz */
			spanne_ges = head - z_bottom;
			spanne_rel = spanne_ges / anz_zeilen;
			m_pcs->m_msh->nod_vector[strang[0]]->SetZ(head);

			if (spanne_ges != 0)
				/* Setzen der neuen Z-Werte entlang eines Stranges */
				for (i = 1; i < anz_zeilen; i++) /* Schleife �ber Anzahl der Zeilen */
					m_pcs->m_msh->nod_vector[strang[i]]->SetZ(head - i * spanne_rel);

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
		long no_nodes_per_layer = (long)nod_vector.size() / (getNumberOfMeshLayers() + 1);
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
				mobile_nodes[i * no_nodes_per_layer + j] = j + (layer_start - 1 + i) * no_nodes_per_layer;
	}
	//......................................................................
	// SURFACE
	if (m_pcs->geo_type.find("SURFACE") != std::string::npos)
	{
		Surface* m_sfc = NULL;
		// CC
		m_sfc = GEOGetSFCByName(m_pcs->geo_type_name);
		if (m_sfc)
			// CC
			mobile_nodes = GetPointsIn(m_sfc, &no_mobile_nodes);
		else
			std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data"
			          << "\n";
	}
	//......................................................................
	// VOLUME
	/*OK411
	   if(m_pcs->geo_type.find("VOLUME")!=std::string::npos)
	   {
	    CGLVolume *m_vol = NULL;
	    m_vol = GEOGetVOL(m_pcs->geo_type_name);//CC 10/05
	    if(m_vol)
	      mobile_nodes = GetPointsInVolume(m_vol,&no_mobile_nodes);//CC 10/05
	    else
	      std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data" << "\n";
	   }
	 */
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
   MSHLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   ToDo evtl. vector<CGLPoint>
   08/2005 CC Modification: CGLPoint* e_pnt - Move from GeoLib to MshLib
**************************************************************************/
void MSHGetNodesClose(std::vector<long>& msh_point_vector, CGLPoint* e_pnt)
{
	e_pnt = e_pnt;
	msh_point_vector.size();
	/*OK411
	   long i;
	   CGLPoint m_pnt;
	   // Node loop
	   for (i=0;i<NodeListSize();i++) {
	    if (GetNode(i)==NULL) continue;
	    m_pnt.x = GetNodeX(i);
	    m_pnt.y = GetNodeY(i);
	    m_pnt.z = GetNodeZ(i);
	    if(e_pnt->PointDis(&m_pnt)<=(e_pnt->epsilon+MKleinsteZahl))
	      msh_point_vector.push_back(i);
	   }
	 */
}

/*************************************************************************
   ROCKFLOW - Function: MSHGetNodesClose
   Task: Searching grid points which are close to a polyline
   Programming:
   10/2002 OK Encapsulated from ExecuteSourceSinkMethod11 (CT)
   01/2003 OK Test
   last modified: 20.01.2003 OK
   08/2005 CC Modification Move from GeoLib to MSHLib
 **************************************************************************/
long* MSHGetNodesClose(long* number_of_nodes, CGLPolyline* m_ply)
{
	long* nodes_all = NULL;
	m_ply = m_ply;
	number_of_nodes = number_of_nodes;
	/*OK411
	   long j,k,l;
	   double pt1[3],line1[3],line2[3],pt0[3];
	   double mult_eps = 1.0;
	   double dist1p,dist2p,*length,laenge;
	   long anz_relevant = 0;
	   typedef struct {
	     long knoten;
	     long abschnitt;
	     double laenge;
	   } INFO;
	   INFO *relevant=NULL;
	   int weiter;
	   double w1,w2;
	   long knoten_help;
	   double laenge_help;
	   double gesamte_laenge = 0.;
	   long polyline_point_vector_size;

	   m_ply->sbuffer.clear();
	   m_ply->ibuffer.clear();

	   if (m_ply) {

	   length = (double*) Malloc(sizeof(double) *(long)m_ply->point_vector.size());

	   pt0[0] = m_ply->point_vector[0]->x;
	   pt0[1] = m_ply->point_vector[0]->y;
	   pt0[2] = m_ply->point_vector[0]->z;

	   polyline_point_vector_size =(long)m_ply->point_vector.size();
	   for (k=0;k<polyline_point_vector_size-1;k++) {
	   line1[0] = m_ply->point_vector[k]->x;
	   line1[1] = m_ply->point_vector[k]->y;
	   line1[2] = m_ply->point_vector[k]->z;
	   line2[0] = m_ply->point_vector[k+1]->x;
	   line2[1] = m_ply->point_vector[k+1]->y;
	   line2[2] = m_ply->point_vector[k+1]->z;
	   length[k] = MCalcDistancePointToPoint(line2, line1);
	   gesamte_laenge += length[k];
	   }

	   // Wiederholen bis zumindest ein Knoten gefunden wurde
	   while(anz_relevant==0) {

	   for (j=0;j<NodeListSize();j++) {
	   if (GetNode(j)==NULL) continue;

	   polyline_point_vector_size =(long)m_ply->point_vector.size();
	   for (k=0;k<polyline_point_vector_size-1;k++) {

	   pt1[0] = GetNodeX(j);
	   pt1[1] = GetNodeY(j);
	   pt1[2] = GetNodeZ(j);

	   line1[0] = m_ply->point_vector[k]->x;
	   line1[1] = m_ply->point_vector[k]->y;
	   line1[2] = m_ply->point_vector[k]->z;
	   line2[0] = m_ply->point_vector[k+1]->x;
	   line2[1] = m_ply->point_vector[k+1]->y;
	   line2[2] = m_ply->point_vector[k+1]->z;

	   if ( MCalcDistancePointToLine(pt1,line1,line2) <= mult_eps*m_ply->epsilon ) {
	   MCalcProjectionOfPointOnLine(pt1,line1,line2,pt1);
	   dist1p = MCalcDistancePointToPoint(line1, pt1);
	   dist2p = MCalcDistancePointToPoint(line2, pt1);
	   if ((dist1p+dist2p-length[k]) <=  mult_eps*m_ply->epsilon ) {

	   // For boundary conditions. WW
	   m_ply->sbuffer.push_back(dist1p);
	   m_ply->ibuffer.push_back(k);
	   // ---------------------------

	   anz_relevant++;
	   nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
	   relevant = (INFO *) Realloc(relevant, sizeof(INFO) * anz_relevant);
	   nodes_all[anz_relevant-1] = j;
	   laenge = 0.;
	   for (l=0; l < k; l++)
	   laenge += length[l];
	   relevant[anz_relevant-1].knoten = j;
	   relevant[anz_relevant-1].laenge = laenge + dist1p;
	   k =(long)m_ply->point_vector.size();
	   }
	   }
	   }
	   }
	   if(anz_relevant==0) mult_eps *=2.;
	   }

	   if (mult_eps > 1.)
	   std::cout << "!!! Epsilon increased in sources!" << "\n";

	   do {
	   weiter = 0;
	   for (k=0;k<anz_relevant-1;k++) {
	   w1=relevant[k].laenge;
	   w2=relevant[k+1].laenge;
	   if (w1>w2) { // Die Eintraege vertauschen
	   knoten_help = relevant[k].knoten;
	   laenge_help = relevant[k].laenge;
	   relevant[k].knoten = relevant[k+1].knoten;
	   relevant[k].laenge = relevant[k+1].laenge;
	   relevant[k+1].knoten = knoten_help;
	   relevant[k+1].laenge = laenge_help;
	   weiter=1;
	   }
	   }
	   } while (weiter);

	   relevant = (INFO*) Free(relevant);
	   *number_of_nodes = anz_relevant;
	   }
	 */
	return nodes_all;
}

/**************************************************************************
   GeoLib-Method: GetPointsIn
   Task:
   Programing:
   01/2004 OK Implementation
   08/2005 CC Modification Move from Geolib to Mshlib
**************************************************************************/
long* GetPointsIn(Surface* m_sfc, long* number_of_nodes)
{
	long* nodes = NULL;
	number_of_nodes = number_of_nodes;
	m_sfc = m_sfc;
	/*OK411
	   long i;
	   double *xp=NULL,*yp=NULL,*zp=NULL;
	   long anz_relevant = 0;
	   CGLPoint m_pnt;
	   // Inside polygon
	   if(!m_sfc->polygon_point_vector.empty()) {
	    xp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    yp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    zp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    long polygon_point_vector_length = (long)m_sfc->polygon_point_vector.size();
	   for(i=0;i<polygon_point_vector_length;i++) {
	   xp[i] = m_sfc->polygon_point_vector[i]->x;
	   yp[i] = m_sfc->polygon_point_vector[i]->y;
	   zp[i] = m_sfc->polygon_point_vector[i]->z;
	   }

	   //-----------------------------------------------------------------
	   for(i=0;i<NodeListSize();i++) {
	   if (GetNode(i)==NULL) continue;
	   m_pnt.x = GetNodeX(i);
	   m_pnt.y = GetNodeY(i);
	   m_pnt.z = GetNodeZ(i);
	   if(m_pnt.IsInsidePolygonPlain(
	   xp,yp,zp,\
	   (long)m_sfc->polygon_point_vector.size())) {
	   anz_relevant++;
	   nodes = (long *) Realloc(nodes,sizeof(long)*anz_relevant);
	   nodes[anz_relevant-1] = i;
	   }
	   }
	   }
	   // Destructions
	   // nodes extern
	   xp = (double*) Free(xp);
	   yp = (double*) Free(yp);
	   zp = (double*) Free(zp);
	   //
	   *number_of_nodes = anz_relevant;
	 */
	return nodes;
}

/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: Get element nodes in a material domain
   Programing:
   10/2004 WW Implementation
**************************************************************************/
void GEOGetNodesInMaterialDomain(const int MatIndex, std::vector<long>& Nodes)
{
	(void)MatIndex;
	(void)Nodes;
	/*OK411
	   MatIndex;
	   Nodes.size();
	   long index, *element_nodes;
	   int i, j, Size, nn, order = 2;
	   const int L_Nodes = GetLowOrderNodeNumber();
	   bool exist;
	   if(L_Nodes==NodeListSize()) order = 1;
	   if(L_Nodes==0) order = 1;

	   Nodes.resize(0);
	   nn = 0;
	   for (index=0;index<ElListSize();index++)
	   {
	   if (ElGetElement(index)!=NULL)
	   {  // Eelement exist
	   if (ElGetElementActiveState(index))
	   {  // Element active
	   if(order==1) nn = NumbersOfElementNode(index);
	   if(order==2) nn = NumbersOfElementNodeHQ(index);
	   if(ElGetElementGroupNumber(index)==MatIndex)
	   {
	   Size = (int)Nodes.size();
	   element_nodes = ElGetElementNodes(index);
	   for(i=0; i<nn; i++)
	   {
	   exist = false;
	   for(j=0; j<Size; j++)
	   {
	   if(element_nodes[i]==Nodes[j])
	   {
	   exist = true;
	   break;
	   }
	   }
	   if(!exist) Nodes.push_back(element_nodes[i]);
	   }
	   }
	   }
	   }
	   }
	 */
}

/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: Get element nodes in a material domain
   Programing:
   10/2004 WW Implementation
   06/2012 NW Made this function faster using std::set for large data set
**************************************************************************/
void GEOGetNodesInMaterialDomain(CFEMesh const* const msh, int MatIndex, std::vector<long>& Nodes, bool Order)
{
	Nodes.resize(0);
	std::set<long> set_nodes;
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
					set_nodes.insert(elem->GetNodeIndex(i));
				}
			}
		} // if
	} // For
	Nodes.assign(set_nodes.begin(), set_nodes.end());
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation based on algorithm originally by CT
   09/2005 CC Move from geolib to mshlib
**************************************************************************/
void SetRFIPointsClose(CGLLine* m_lin)
{
	m_lin = m_lin;
	/*OK411
	   long j;
	   double pt1[3],pt2[3],line1[3],line2[3];
	   double mult_eps = 1.0;
	   double dist1p,dist2p,length;
	   long anz_relevant;
	   double *dist;
	   typedef struct {
	     long knoten;
	     long abschnitt;
	     double laenge;
	   } INFO;
	   INFO *relevant=NULL;
	   long knoten_help;
	   double laenge_help;
	   double w1,w2;
	   int weiter;
	   //----------------------------------------------------------------------
	   // Tests
	   if(!ELEListExists()) {
	   return;
	   }
	   if(m_lin->point1<0)
	   return;
	   if(m_lin->point2<0)
	   return;
	   //----------------------------------------------------------------------
	   // Initializations
	   anz_relevant = 0;
	   m_lin->msh_nodes = NULL;
	   CGLPoint *m_point=NULL;
	   m_point = GEOGetPointById(m_lin->point1);//CC
	   line1[0] = m_point->x;
	   line1[1] = m_point->y;
	   line1[2] = m_point->z;
	   m_point = GEOGetPointById(m_lin->point2);//CC
	   line2[0] = m_point->x;
	   line2[1] = m_point->y;
	   line2[2] = m_point->z;
	   length = MCalcDistancePointToPoint(line2,line1);
	   //----------------------------------------------------------------------
	   // Repeat untill at least one node is found
	   while(anz_relevant==0) {
	   // NOD list loop
	   for (j=0;j<NodeListSize();j++) {
	   if (GetNode(j)==NULL) continue;
	   pt1[0] = GetNodeX(j);
	   pt1[1] = GetNodeY(j);
	   pt1[2] = GetNodeZ(j);
	   // Is MSH point near to line
	   if ( MCalcDistancePointToLine(pt1,line1,line2) <= mult_eps*m_lin->epsilon ) {
	   // Calc projection of pt1 to line and use this in the following
	   MCalcProjectionOfPointOnLine(pt1,line1,line2,pt1);
	   // Abstand des Punktes zum ersten Punkt des Polygonabschnitts
	   dist1p = MCalcDistancePointToPoint(line1, pt1);
	   // Abstand des Punktes zum zweiten Punkt des Polygonabschnitts
	   dist2p = MCalcDistancePointToPoint(line2, pt1);
	   // Ist der Knoten innerhalb des Intervalls?
	   if ((dist1p+dist2p-length) <=  mult_eps*m_lin->epsilon ) {
	   anz_relevant++;
	   // Feld anpassen
	   m_lin->msh_nodes = (long *) Realloc(m_lin->msh_nodes,sizeof(long)*anz_relevant);
	   relevant = (INFO *) Realloc(relevant, sizeof(INFO)*anz_relevant);
	   // Ablegen von Knotennummer und Position
	   m_lin->msh_nodes[anz_relevant-1] = j;
	   }
	   } // endif
	   } // Ende Schleife ueber Knoten
	   if(anz_relevant==0) mult_eps *=2.;
	   } // Ende Schleife Wiederholungen
	   if (mult_eps > 1.)
	   printf("Warning: Epsilon increased in CGLLine::SetRFIPointsClose.");
	   m_lin->no_msh_nodes = anz_relevant;
	   //----------------------------------------------------------------------
	   // Sort MSH nodes, beginning from first line point
	   //......................................................................
	   // Calc distances from first line point
	   m_point = GEOGetPointById(m_lin->point1);//CC
	   pt1[0] = m_point->x;
	   pt1[1] = m_point->y;
	   pt1[2] = m_point->z;
	   dist = (double*) Malloc(sizeof(double)*m_lin->no_msh_nodes);
	   for(j=0;j<m_lin->no_msh_nodes;j++) {
	   pt2[0] = GetNodeX(m_lin->msh_nodes[j]);
	   pt2[1] = GetNodeY(m_lin->msh_nodes[j]);
	   pt2[2] = GetNodeZ(m_lin->msh_nodes[j]);
	   dist[j] = MCalcDistancePointToPoint(pt1,pt2);
	   }
	   //......................................................................
	   // Sorting by pair switching
	   do {
	   weiter = 0;
	   for (j=0;j<m_lin->no_msh_nodes-1;j++) {
	   w1=dist[j];
	   w2=dist[j+1];
	   if (w1>w2) { // Die Eintraege vertauschen
	   knoten_help = m_lin->msh_nodes[j];
	   laenge_help = dist[j];
	   m_lin->msh_nodes[j] = m_lin->msh_nodes[j+1];
	   dist[j] = dist[j+1];
	   m_lin->msh_nodes[j+1] = knoten_help;
	   dist[j+1] = laenge_help;
	   weiter=1;
	   }
	   }
	   } while (weiter);
	   //----------------------------------------------------------------------
	   // Destructions
	   dist = (double*) Free(dist);
	   relevant = (INFO*) Free(relevant);
	 */
}

/**************************************************************************
   GeoSys-GUI Function
   Programing:
   01/2004 OK Implementation
   09/2005 CC Modification No MFC function
**************************************************************************/
// bool IsPointInSurface(Surface* m_fsc, CGLPoint *m_point)
//{
//  bool ok = false;
//  m_point = m_point;
//  if(!m_fsc)
//    return false; //OK
//  return ok;
//}

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
			if ((m_msh->ele_vector[j]->GetPatchIndex() + 1) > msh_max_mmp_groups)
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
		max_mat_id = std::max(max_mat_id, m_msh->ele_vector[j]->GetPatchIndex());
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
		return false; // abort();
	}
	return true;
}
