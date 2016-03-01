/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "vtk.h"
#include <fstream>
#if defined(WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
#include <sstream>
#endif
#include "Output.h"
#include "fem_ele_std.h"                          // for element velocity
#include "makros.h"
#include "rf_mmp_new.h"
#include "FileTools.h"

using namespace std;

const std::string INDEX_STR = "  ";
const std::string velocity_name[3][4] =
{
	{
		"VELOCITY_X1", "VELOCITY_Y1", "VELOCITY_Z1", "NODAL_VELOCITY1"
	}
	,
	{
		"VELOCITY_X2", "VELOCITY_Y2", "VELOCITY_Z2", "NODAL_VELOCITY2"
	}
	,
	{
		"VELOCITY1_X", "VELOCITY1_Y", "VELOCITY1_Z", "GL_NODAL_VELOCITY1"
	}
};

//#################################################################################################
// Functions for Paraview Data File (PVD)

bool CVTK::InitializePVD(const string &file_base_name, const string &pcs_type_name, bool binary)
{
	//PVD
	vec_dataset.clear();
	pvd_file_name = pathJoin(defaultOutputPath, pathBasename(file_base_name));
	// pvd_file_name = file_base_name;
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
		pvd_file_name += mrank_str;
#endif
	if(pcs_type_name.size() > 0)          // PCS
		pvd_file_name += "_" + pcs_type_name;
	pvd_file_name += ".pvd";

	/* // Make the following lines as comments by WW
	//VTK
	int ibs = (int)file_base_name.find_last_of("\\");
	int is = (int)file_base_name.find_last_of("/");
	//OK411
	if (ibs != (int)string::npos  || is != (int)string::npos)
	{
		int ibegin = ibs;
		if (is > ibs)
			ibegin = is;
		ibegin += 1;
		this->pvd_vtk_file_name_base = file_base_name.substr(ibegin);
		this->pvd_vtk_file_path_base = file_base_name.substr(0, ibegin);
	}
	else
    {
		this->pvd_vtk_file_name_base = file_base_name;
		this->pvd_vtk_file_path_base = "";
	}
	if (pcs_type_name.size() > 0)        // PCS
		this->pvd_vtk_file_name_base += "_" + pcs_type_name;
	*/

	//
	pvd_vtk_file_path_base = defaultOutputPath;
	pvd_vtk_file_name_base = pathBasename(file_base_name) + "_" + pcs_type_name; //WW
	useBinary = binary;

	return true;
}

bool CVTK::WriteHeaderOfPVD(std::fstream &fin)
{
	fin << "<?xml version=\"1.0\"?>" << "\n";
	fin <<
	"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
	    << "\n";
	fin << INDEX_STR << "<Collection>" << "\n";
	return true;
}

bool CVTK::WriteEndOfPVD(std::fstream &fin)
{
	fin << INDEX_STR << "</Collection>" << "\n";
	fin << "</VTKFile>" << "\n";
	return true;
}

bool CVTK::WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile)
{
	fin.setf(ios::scientific,std::ios::floatfield);
	fin.precision(12);
	fin << INDEX_STR << INDEX_STR << "<DataSet timestep=\"" << timestep <<
	"\" group=\"\" part=\"0\" file=\"" << vtkfile << "\"/>" << "\n";
	return true;
}

bool CVTK::UpdatePVD(const string &pvdfile, const vector<VTK_Info> &vec_vtk)
{
	fstream fin(pvdfile.data(), ios::out);
	if (!fin.good())
		return false;

	CVTK::WriteHeaderOfPVD(fin);
	for (int i = 0; i < (int)vec_vtk.size(); i++)
		CVTK::WriteDatasetOfPVD(fin, vec_vtk[i].timestep, vec_vtk[i].vtk_file);
	CVTK::WriteEndOfPVD(fin);

	fin.close();

	return true;
}

bool CVTK::CreateDirOfPVD(const string &pvdfile)
{
	string pvd_dir_path = pvdfile + ".d";
#if defined(WIN32)
	if (mkdir(pvd_dir_path.c_str()) == -1)
	{
#else
	if (mkdir(pvd_dir_path.c_str(), 0777) == -1)
	{
#endif
		//error
		cout << "***ERROR: Fail to create a PVD directory: " << pvd_dir_path << "\n";
		return false;
	}
	return true;
}

//#################################################################################################
// Functions for VTU files (XML UnstructuredGrid file)

unsigned char CVTK::GetVTKCellType(const MshElemType::type ele_type)
{
	unsigned char cell_type = 0;

	switch(ele_type)
	{
	case MshElemType::LINE:               // vtk_line=3
		cell_type = 3;
		break;
	case MshElemType::QUAD:               // quadrilateral=9
		cell_type = 9;
		break;
	case MshElemType::HEXAHEDRON:         // hexahedron=12
		cell_type = 12;
		break;
	case MshElemType::TRIANGLE:           // triangle=5
		cell_type = 5;
		break;
	case MshElemType::TETRAHEDRON:        // tetrahedron=10
		cell_type = 10;
		break;
	case MshElemType::PRISM:              // wedge=13
		cell_type = 13;
		break;
	case MshElemType::PYRAMID:              // pyramid=14
		cell_type = 14;
		break;
	default:
		std::cerr << "***ERROR: NO CORRESPONDING VTK CELL TYPE FOUND. (ELEMENT TYPE=" <<
		ele_type << ")" << "\n";
	}
	return cell_type;
}

void CVTK::InitializeVTU()
{
	//if (this->useBinary) {
	//======================================================================
	//# Set machine dependent stuff
	//Data type
	if (sizeof(unsigned char) == 1)
		type_UChar = CVTK::UInt8;
	else if (sizeof(unsigned char) == 2)
		type_UChar = CVTK::UInt16;
	if (sizeof(int) == 4)
		type_Int = CVTK::Int32;
	else if (sizeof(int) == 8)
		type_Int = CVTK::Int64;
	if (sizeof(unsigned int) == 4)
		type_UInt = CVTK::UInt32;
	else if (sizeof(unsigned int) == 8)
		type_UInt = CVTK::UInt64;
	if (sizeof(long) == 4)
		type_Long = CVTK::Int32;
	else if (sizeof(long) == 8)
		type_Long = CVTK::Int64;
	if (sizeof(double) == 4)
		type_Double = CVTK::Float32;
	else if (sizeof(double) == 8)
		type_Double = CVTK::Float64;
	//
	SIZE_OF_BLOCK_LENGTH_TAG = sizeof(unsigned int);
	//Endian(byte order)
	isLittleEndian = IsLittleEndian();
	//}

	this->isInitialized = true;
}

bool CVTK::WriteDataArrayHeader(std::fstream &fin,
                                VTK_XML_DATA_TYPE data_type,
                                const std::string &str_name,
                                int nr_components,
                                const std::string &str_format,
                                long offset)
{
	std::string str_data_type;
	switch (data_type)
	{
	case CVTK::Int8: str_data_type = "Int8";
		break;
	case CVTK::UInt8: str_data_type = "UInt8";
		break;
	case CVTK::Int16: str_data_type = "Int16";
		break;
	case CVTK::UInt16: str_data_type = "UInt16";
		break;
	case CVTK::Int32: str_data_type = "Int32";
		break;
	case CVTK::UInt32: str_data_type = "UInt32";
		break;
	case CVTK::Int64: str_data_type = "Int64";
		break;
	case CVTK::UInt64: str_data_type = "UInt64";
		break;
	case CVTK::Float32: str_data_type = "Float32";
		break;
	case CVTK::Float64: str_data_type = "Float64";
		break;
	}
	fin << "        <DataArray type=\"" << str_data_type << "\"";
	if (str_name != "")
		fin << " Name=\"" << str_name << "\"";
	if (nr_components > 1)
		fin << " NumberOfComponents=\"" << nr_components << "\"";
	fin << " format=\"" << str_format << "\"";
	if (useBinary)
		fin << " offset=\"" << offset << "\" /";
	fin << ">" << "\n";

	return true;
}

bool CVTK::WriteDataArrayFooter(std::fstream &fin)
{
	if (!this->useBinary)
		fin << "        </DataArray>" << "\n";

	return true;
}

bool CVTK::WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                    COutput* out,
                                    const int time_step_number)
{
	if (!this->isInitialized)
		this->InitializeVTU();

	//-------------------------------------------------------------------------
	//# Setup file stream
	//-------------------------------------------------------------------------
	std::fstream fin;
	if (this->useBinary)
		fin.open(vtkfile.data(), std::ios::out | std::ios::binary);
	else
		fin.open(vtkfile.data(), std::ios::out);

	if (!fin.good())
	{
		std::cout << "***Warning: Cannot open the output file, " << vtkfile << "\n";
		return false;
	}

	if (!this->useBinary)
	{
		fin.setf(std::ios::scientific,std::ios::floatfield);
		fin.precision(12);
	}

	//-------------------------------------------------------------------------
	//# Output
	//-------------------------------------------------------------------------
	//  CFEMesh *msh = out->GetMSH();
	CFEMesh* msh = out->getMesh();
	long offset = 0;

	string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";
	bool data_out = !useBinary;

	//# Header
	fin << "<?xml version=\"1.0\"?>" << "\n";
	fin << "<!-- Time step: " << time_step_number << " | Time: " << out->getTime() << " -->" <<
	"\n";
	fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
	if (!this->useBinary || isLittleEndian)
		fin << " byte_order=\"LittleEndian\"";
	else
		fin << " byte_order=\"BigEndian\"";
	fin << ">" << "\n";
	//  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << "\n";

	//# Unstructured Grid information
	fin << "  <UnstructuredGrid>" << "\n";
	fin << "    <Piece NumberOfPoints=\"" << msh->GetNodesNumber(false) <<
	"\" NumberOfCells=\"" << msh->ele_vector.size() << "\">" << "\n";
	//....................................................................
	// Nodes
	//OK411 CNode *nod = NULL;
	fin << "      <Points>" << "\n";
	WriteDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
	WriteMeshNodes(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	fin << "      </Points>" << "\n";
	//....................................................................
	// Elements
	//OK411 CElem * ele = NULL;
	fin << "      <Cells>" << "\n";
	//connectivity
	WriteDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
	long sum_ele_components = 0;
	WriteMeshElementConnectivity(fin, data_out, msh, offset, sum_ele_components);
	WriteDataArrayFooter(fin);
	//offset
	WriteDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
	WriteMeshElementOffset(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	//type
	WriteDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
	WriteMeshElementType(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	fin << "      </Cells>" << "\n";
	//....................................................................
	// Nodal values
	if (out->_nod_value_vector.size() > 0)
		fin << "      <PointData Scalars=\"" << out->_nod_value_vector[0] << "\">" << "\n";
	else
		fin << "      <PointData Scalars=\"scalars\">" << "\n";
	WriteNodalValue(fin, data_out, out, msh, offset);
	fin << "      </PointData>" << "\n";

	//======================================================================
	//....................................................................
	// Element values
	fin << "      <CellData>" << "\n";
	WriteElementValue(fin, data_out, out, msh, offset);
	fin << "      </CellData>" << "\n";
	fin << "    </Piece>" << "\n";
	fin << "  </UnstructuredGrid>" << "\n";

	//======================================================================
	// Raw data (for binary mode)
	if (useBinary)
	{
		fin << "  <AppendedData encoding=\"raw\">" << "\n";
		fin << "    _";

		//Node
		this->WriteMeshNodes(fin, true, msh, offset);
		//Element
		//conncectivity
		this->WriteMeshElementConnectivity(fin, true, msh, offset, sum_ele_components);
		//offset
		this->WriteMeshElementOffset(fin, true, msh, offset);
		//type
		this->WriteMeshElementType(fin, true, msh, offset);
		// Nodal values
		this->WriteNodalValue(fin, true, out, msh, offset);
		// Elemental values
		this->WriteElementValue(fin, true, out, msh, offset);

		fin << "\n";
		fin << "  </AppendedData>" << "\n";
	}

	fin << "</VTKFile>" << "\n";
	fin.close();

	return true;
}

bool CVTK::IsLittleEndian()
{
	int x = 0x00000001;
	if (*(char*)&x)
		return true;              //am little
	else
		return false;             //am big
}

template <typename T> void CVTK::write_value_binary(std::fstream &fin, T val)
{
	fin.write((const char*)&val, sizeof(T));
}

bool CVTK::WriteMeshNodes(std::fstream &fin, bool output_data, CFEMesh* msh, long &offset)
{
	const size_t n_msh_nodes = msh->GetNodesNumber(false);
	if (output_data)
	{
		if (!useBinary)
			for (size_t i = 0; i < n_msh_nodes; i++)
			{
				double const* const pnt (msh->nod_vector[i]->getData());
				fin << "          " << pnt[0] << " " << pnt[1] << " " << pnt[2] <<
				"\n";
			}
		else
		{
			//OK411
			write_value_binary<unsigned int>(fin, sizeof(double) * 3 * n_msh_nodes);
			for (size_t i = 0; i < n_msh_nodes; i++)
			{
				double const* const pnt (msh->nod_vector[i]->getData());
				write_value_binary(fin, pnt[0]);
				write_value_binary(fin, pnt[1]);
				write_value_binary(fin, pnt[2]);
			}
		}
	}
	else if (useBinary)
		//OK411
		offset += n_msh_nodes * sizeof(double) * 3 + SIZE_OF_BLOCK_LENGTH_TAG;

	return true;
}

bool CVTK::WriteMeshElementConnectivity(std::fstream &fin,
                                        bool output_data,
                                        CFEMesh* msh,
                                        long &offset,
                                        long &sum_ele_components)
{
	if (output_data)
	{
		MeshLib::CElem* ele = NULL;
		if (!useBinary)
			for (long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				fin << "          ";
				for (size_t j = 0; j < ele->GetNodesNumber(false); j++)
					fin << ele->GetNodeIndex(j) << " ";
				fin << "\n";
			}
		else
		{
			write_value_binary<unsigned int>(fin, sizeof(long) * sum_ele_components);
			for (long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				for (size_t j = 0; j < msh->ele_vector[i]->GetNodesNumber(false);
				     j++)
					write_value_binary<long>(fin, ele->GetNodeIndex(j));
			}
		}
	}
	else if (useBinary)
	{
		sum_ele_components = 0;
		for (size_t i = 0; i < msh->ele_vector.size(); i++)
			sum_ele_components += msh->ele_vector[i]->GetNodesNumber(false);
		offset += sum_ele_components * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
	}

	return true;
}

bool CVTK::WriteMeshElementOffset(std::fstream &fin, bool output_data, CFEMesh* msh, long &offset)
{
	if (output_data)
	{
		MeshLib::CElem* ele = NULL;

		if (!useBinary)
		{
			fin << "          ";
			long ele_offset = 0;
			for (long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				ele_offset += ele->GetNodesNumber(false);
				fin << ele_offset << " ";
			}
			fin << "\n";
		}
		else
		{
			//OK411
			write_value_binary<unsigned int>(fin,
			                                 sizeof(long) * (long)msh->ele_vector.size());
			long ele_offset = 0;
			for (long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				ele_offset += ele->GetNodesNumber(false);
				write_value_binary(fin, ele_offset);
			}
		}
	}
	else if (useBinary)
		//OK411
		offset += (long)msh->ele_vector.size() * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;

	return true;
}

bool CVTK::WriteMeshElementType(std::fstream &fin, bool output_data, CFEMesh* msh, long &offset)
{
	if (output_data)
	{
		MeshLib::CElem* ele = NULL;
		if (!useBinary)
		{
			fin << "          ";
			for(long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				fin << (int)this->GetVTKCellType(ele->GetElementType()) << " ";
			}
			fin << "\n";
		}
		else
		{
			//OK411
			write_value_binary<unsigned int>(
			        fin,
			        sizeof(unsigned char) *
			        (long)msh->ele_vector.size());
			for(long i = 0; i < (long)msh->ele_vector.size(); i++)
			{
				ele = msh->ele_vector[i];
				write_value_binary(fin, this->GetVTKCellType(ele->GetElementType()));
			}
		}
	}
	else if (useBinary)
		//OK411
		offset += (long)msh->ele_vector.size() * sizeof(unsigned char) +
		          SIZE_OF_BLOCK_LENGTH_TAG;

	return true;
}

bool CVTK::WriteNodalValue(std::fstream &fin,
                           bool output_data,
                           COutput* out,
                           CFEMesh* msh,
                           long &offset)
{
	CRFProcess* m_pcs = NULL;
	std::vector<int> NodeIndex(out->_nod_value_vector.size());
	//  if (out->m_pcs == NULL && out->pcs_type_name.compare("NO_PCS")!=0)
	//      out->m_pcs = PCSGet(out->pcs_type_name);
	if (out->m_pcs == NULL && out->getProcessType() == FiniteElement::NO_PCS)
		out->m_pcs = PCSGet(out->getProcessType());
	if (out->m_pcs != NULL)
		m_pcs = out->m_pcs;

 	string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";

    bool isXZplane = (msh->GetCoordinateFlag()==22);
    bool is3D = (msh->GetCoordinateFlag() / 10 == 3);
	bool outNodeVelocity = false;
    bool outNodeDisplacement = false;
    bool outNodePrincipleStressDirections = false;

	//Nodal values
	for (int i = 0; i < (int) out->_nod_value_vector.size(); i++)
	{
        const string &internal_val_name = out->_nod_value_vector[i];
        const string &external_val_name = out->_alias_nod_value_vector[i];
		//is velocity
		if (internal_val_name.find("VELOCITY") != string::npos)
		{
			outNodeVelocity = true;
			continue;
		}
        if (internal_val_name.find("DISPLACEMENT") != string::npos)
        {
            outNodeDisplacement = true;
            continue;
        }

        if (internal_val_name.find("NORM_STRESS") != string::npos)
        {
                            outNodePrincipleStressDirections = true;
                            continue;
        }


		//    if (out->m_pcs == NULL || out->pcs_type_name.compare("NO_PCS")==0)
		if (out->m_pcs == NULL || out->getProcessType() == FiniteElement::NO_PCS)
			m_pcs = PCSGet(internal_val_name, true);
		if (!m_pcs)
			continue;

		NodeIndex[i] = m_pcs->GetNodeValueIndex(internal_val_name,true); //JT: Latest
		if (NodeIndex[i] < 0)
			continue;


		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, type_Double, external_val_name, 0, str_format, offset);

		if (output_data)
		{
			/* JT: Just get the latest index. No need for all this extra looping.
			for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
				if (internal_val_name.compare(m_pcs->pcs_primary_function_name[j]) == 0) {
					NodeIndex[i]++; //current step
					break;
				}
			*/
			if (!useBinary) {
				fin << "          ";
			} else {
				write_value_binary<unsigned int> (fin, sizeof(double)* msh->GetNodesNumber(false));
			}
            for (size_t j = 0; j < msh->GetNodesNumber(false); j++) {
                double v = m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), NodeIndex[i]);
                if (!useBinary) {
                    fin << v << " ";
                } else {
                    write_value_binary(fin, v);
                }
            }
            if (!useBinary) {
                fin << "\n";
            }
		}
		else
			offset += msh->GetNodesNumber(false) * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;

		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

	// Nodal velocities
	if (outNodeVelocity)
	{
		unsigned int velocity_id = 0;
		for (int i = 0; i < (int) out->_nod_value_vector.size(); i++)
		{
            const string &internal_val_name = out->_nod_value_vector[i];
//            const string &external_val_name = out->_alias_nod_value_vector[i];
			if (internal_val_name.find("VELOCITY_X1") != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 0;
			}
			else if (internal_val_name.find("VELOCITY_X2")
			         != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 1;
			}
			else if (internal_val_name.find("VELOCITY1_X")
			         != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 2;
			}
			else
				continue;
			if (!m_pcs)
				continue;

			if (!useBinary || !output_data)
				WriteDataArrayHeader(fin,
				                     this->type_Double,
				                     velocity_name[velocity_id][3],
				                     3,
				                     str_format,
				                     offset);
			if (output_data)
			{
				int ix, iy, iz;
				ix = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][0],true); // JT: Fix. Need latest value.
				iy = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][1],true);
				iz = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][2],true);
				if (!useBinary)
				{
					fin << "          ";
					for (size_t j = 0l; j < msh->GetNodesNumber(false); j++)
					{
						fin << m_pcs->GetNodeValue(
						        msh->nod_vector[j]->GetIndex(), ix) << " ";
                        if (m_pcs->getProcessType()==FiniteElement::FLUID_MOMENTUM || !isXZplane) {
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iy) << " ";
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iz) << " ";
                        } else {
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iz) << " ";
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iy) << " ";
                        }
					}
					fin << "\n";
				}
				else
				{
					write_value_binary<unsigned int> (fin, sizeof(double)
					                                  * msh->GetNodesNumber(
					                                          false) * 3);
					for (size_t j = 0l; j < msh->GetNodesNumber(false); j++)
					{
						write_value_binary(fin, m_pcs->GetNodeValue(
						                           msh->nod_vector[j]->
						                           GetIndex(), ix));
                        if (!isXZplane) {
                            write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iy));
                            write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iz));
                        } else {
                            write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iz));
                            write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iy));
                        }
					}
				}
			}
			else
				offset += msh->GetNodesNumber(false) * 3 * sizeof( double) +
				          SIZE_OF_BLOCK_LENGTH_TAG;

			if (!useBinary || !output_data)
				WriteDataArrayFooter(fin);
		}
	}

    //Displacement
    if (outNodeDisplacement)
    {
//        unsigned int disp_id = 0;
        for (int i = 0; i < (int)out->_nod_value_vector.size(); i++)
        {
            const string &internal_val_name = out->_nod_value_vector[i];
//            const string &external_val_name = out->_alias_nod_value_vector[i];
            if (internal_val_name.find("DISPLACEMENT_X1") != string::npos)
            {
                if (out->m_pcs == NULL)
                    m_pcs = PCSGet(internal_val_name,true);
//                disp_id = 0;
            }
            else
                continue;
            if(!m_pcs)
                continue;

            if (!useBinary || !output_data)
                WriteDataArrayHeader(fin, this->type_Double, "DISPLACEMENT", 3, str_format, offset);
            if (output_data)
            {
                int var_id[3] = {};
                var_id[0] = m_pcs->GetNodeValueIndex("DISPLACEMENT_X1");
                var_id[1] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
                var_id[2] = -1;
                if (is3D) {
                    var_id[2] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
                } else if (isXZplane) {
                    var_id[1] = -1;
                    var_id[2] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
                }
                //
                if (!useBinary) {
                    fin << "          ";
                } else {
                    write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                }
                double u[3] = {};
                for(size_t j = 0l; j < msh->GetNodesNumber(false); j++)
                {
                    for (size_t k=0; k<3; k++) {
                        if (var_id[k]<0)
                            u[k] = .0;
                        else
                            u[k] =  m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), var_id[k]);
                    }

                    if (!useBinary) {
                        for (size_t k=0; k<3; k++)
                            fin << u[k] << " ";
                    } else {
                        for (size_t k=0; k<3; k++)
                            write_value_binary(fin, u[k]);
                    }
                }
                if (!useBinary) {
                    fin << "\n";
                } else {
                    write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                }
            }
            else
                offset += msh->GetNodesNumber(false) * 3 * sizeof(double) +
                SIZE_OF_BLOCK_LENGTH_TAG;
            if (!useBinary || !output_data)
                WriteDataArrayFooter(fin);
        }
    }


    if (outNodePrincipleStressDirections)
        {
            for (size_t i = 0; i < out->_nod_value_vector.size(); i++)
            {
                const string &internal_val_name = out->_nod_value_vector[i];
    //            const string &external_val_name = out->_alias_nod_value_vector[i];
                if (internal_val_name.find("NORM_STRESS_1_X") != string::npos)
                {
                    if (out->m_pcs == NULL)
                        m_pcs = PCSGet(internal_val_name,true);
    //                disp_id = 0;
                }
                else
                    continue;
                if(!m_pcs)
                    continue;

                if (!useBinary || !output_data)
                    WriteDataArrayHeader(fin, this->type_Double, "NORMSTRESS_1", 3, str_format, offset);
                if (output_data)
                {
                    int var_id[3] = {};
                    var_id[0] = m_pcs->GetNodeValueIndex("NORM_STRESS_1_X");
                    var_id[1] = m_pcs->GetNodeValueIndex("NORM_STRESS_1_Y");
                    var_id[2] = m_pcs->GetNodeValueIndex("NORM_STRESS_1_Z");
                    //
                    if (!useBinary) {
                        fin << "          ";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                    double u[3] = {};
                    for(size_t j = 0l; j < msh->GetNodesNumber(false); j++)
                    {
                        for (size_t k=0; k<3; k++) {
                            if (var_id[k]<0)
                                u[k] = .0;
                            else
                                u[k] =  m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), var_id[k]);
                        }

                        if (!useBinary) {
                            for (size_t k=0; k<3; k++)
                                fin << u[k] << " ";
                        } else {
                            for (size_t k=0; k<3; k++)
                                write_value_binary(fin, u[k]);
                        }
                    }
                    if (!useBinary) {
                        fin << "\n";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                }
                else
                    offset += msh->GetNodesNumber(false) * 3 * sizeof(double) +
                    SIZE_OF_BLOCK_LENGTH_TAG;
                if (!useBinary || !output_data)
                    WriteDataArrayFooter(fin);
            }
            for (size_t i = 0; i < out->_nod_value_vector.size(); i++)
            {
                const string &internal_val_name = out->_nod_value_vector[i];
    //            const string &external_val_name = out->_alias_nod_value_vector[i];
                if (internal_val_name.find("NORM_STRESS_2_X") != string::npos)
                {
                    if (out->m_pcs == NULL)
                        m_pcs = PCSGet(internal_val_name,true);
    //                disp_id = 0;
                }
                else
                    continue;
                if(!m_pcs)
                    continue;

                if (!useBinary || !output_data)
                    WriteDataArrayHeader(fin, this->type_Double, "NORMSTRESS_2", 3, str_format, offset);
                if (output_data)
                {
                    int var_id[3] = {};
                    var_id[0] = m_pcs->GetNodeValueIndex("NORM_STRESS_2_X");
                    var_id[1] = m_pcs->GetNodeValueIndex("NORM_STRESS_2_Y");
                    var_id[2] = m_pcs->GetNodeValueIndex("NORM_STRESS_2_Z");
                    //
                    if (!useBinary) {
                        fin << "          ";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                    double u[3] = {};
                    for(size_t j = 0l; j < msh->GetNodesNumber(false); j++)
                    {
                        for (size_t k=0; k<3; k++) {
                            if (var_id[k]<0)
                                u[k] = .0;
                            else
                                u[k] =  m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), var_id[k]);
                        }

                        if (!useBinary) {
                            for (size_t k=0; k<3; k++)
                                fin << u[k] << " ";
                        } else {
                            for (size_t k=0; k<3; k++)
                                write_value_binary(fin, u[k]);
                        }
                    }
                    if (!useBinary) {
                        fin << "\n";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                }
                else
                    offset += msh->GetNodesNumber(false) * 3 * sizeof(double) +
                    SIZE_OF_BLOCK_LENGTH_TAG;
                if (!useBinary || !output_data)
                    WriteDataArrayFooter(fin);
            }
            for (size_t i = 0; i < out->_nod_value_vector.size(); i++)
            {
                const string &internal_val_name = out->_nod_value_vector[i];
    //            const string &external_val_name = out->_alias_nod_value_vector[i];
                if (internal_val_name.find("NORM_STRESS_3_X") != string::npos)
                {
                    if (out->m_pcs == NULL)
                        m_pcs = PCSGet(internal_val_name,true);
    //                disp_id = 0;
                }
                else
                    continue;
                if(!m_pcs)
                    continue;

                if (!useBinary || !output_data)
                    WriteDataArrayHeader(fin, this->type_Double, "NORMSTRESS_3", 3, str_format, offset);
                if (output_data)
                {
                    int var_id[3] = {};
                    var_id[0] = m_pcs->GetNodeValueIndex("NORM_STRESS_3_X");
                    var_id[1] = m_pcs->GetNodeValueIndex("NORM_STRESS_3_Y");
                    var_id[2] = m_pcs->GetNodeValueIndex("NORM_STRESS_3_Z");
                    //
                    if (!useBinary) {
                        fin << "          ";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                    double u[3] = {};
                    for(size_t j = 0l; j < msh->GetNodesNumber(false); j++)
                    {
                        for (size_t k=0; k<3; k++) {
                            if (var_id[k]<0)
                                u[k] = .0;
                            else
                                u[k] =  m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), var_id[k]);
                        }

                        if (!useBinary) {
                            for (size_t k=0; k<3; k++)
                                fin << u[k] << " ";
                        } else {
                            for (size_t k=0; k<3; k++)
                                write_value_binary(fin, u[k]);
                        }
                    }
                    if (!useBinary) {
                        fin << "\n";
                    } else {
                        write_value_binary<unsigned int>(fin, sizeof(double)*msh->GetNodesNumber(false)*3);
                    }
                }
                else
                    offset += msh->GetNodesNumber(false) * 3 * sizeof(double) +
                    SIZE_OF_BLOCK_LENGTH_TAG;
                if (!useBinary || !output_data)
                    WriteDataArrayFooter(fin);
            }




        }
	//MFP
      for ( size_t i=0; i<out->mfp_value_vector.size(); i++) {

          const std::string &mfp_name = out->mfp_value_vector[i];
          if (!useBinary || !output_data)
          {
              WriteDataArrayHeader(fin, type_Double, mfp_name, 0, str_format, offset);
          }

          if (output_data)
          {
              if (!useBinary)
              {
                  fin << "          ";
                  for (size_t j = 0; j < msh->GetNodesNumber(false); j++)
                  {
                      const double v = MFPGetNodeValue(msh->nod_vector[j]->GetIndex(), mfp_name, atoi(&mfp_name[mfp_name.size() - 1]) - 1);
                      fin << v << " ";
                  }
                  fin << endl;
              }
              else
              {
                  write_value_binary<unsigned int> (fin, sizeof(double)
                      * msh->GetNodesNumber(false));
                  for (size_t j = 0; j < msh->GetNodesNumber(false); j++)
                  {
                      const double v = MFPGetNodeValue(msh->nod_vector[j]->GetIndex(), mfp_name, atoi(&mfp_name[mfp_name.size() - 1]) - 1);
                      write_value_binary(fin, v);
                  }
              }
          }
          else
          {
              offset += msh->GetNodesNumber(false) * sizeof(double)
                  + SIZE_OF_BLOCK_LENGTH_TAG;
          }

          if (!useBinary || !output_data)
          {
              WriteDataArrayFooter(fin);
          }
      }
    return true;
}

bool CVTK::WriteElementValue(std::fstream &fin,
                             bool output_data,
                             COutput* out,
                             CFEMesh* msh,
                             long &offset)
{
	std::vector<int> ele_value_index_vector(out->getElementValueVector().size());
	if (ele_value_index_vector.size() > 0) // GetELEValuesIndexVector() should check this!
		out->GetELEValuesIndexVector(ele_value_index_vector);
	CRFProcess* m_pcs = NULL;
	MeshLib::CElem* ele = NULL;

    bool isXZplane = (msh->GetCoordinateFlag()==22);

	string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";

	bool outEleVelocity = false;

	//Element values
	for (int i = 0; i < (int) ele_value_index_vector.size(); i++)
	{
		if (out->getElementValueVector()[i].find("VELOCITY") != string::npos)
		{
			outEleVelocity = true;
			continue;
		}
		m_pcs = out->GetPCS_ELE(out->getElementValueVector()[i]);

		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Double,
			                     out->getElementValueVector()[i], 0, str_format, offset);
		if (output_data)
		{
			if (!useBinary)
			{
				fin << "          ";
				for (long j = 0; j < (long) msh->ele_vector.size(); j++)
					if (ele_value_index_vector[i]>=0)
						fin << m_pcs->GetElementValue(j, ele_value_index_vector[i])
						    << " ";
				fin << "\n";
			}
			else
			{
				write_value_binary<unsigned int> (fin, sizeof(double)
				                                  * (long) msh->ele_vector.size()); //OK411
				for (long j = 0; j < (long) msh->ele_vector.size(); j++)
					if (ele_value_index_vector[i]>=0)
						write_value_binary(fin,
						                   m_pcs->GetElementValue(j,ele_value_index_vector[i]));
			}
		}
		else
			offset += (long) msh->ele_vector.size() * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;  //OK411
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

	//Element velocity
	if (outEleVelocity)
	{
		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Double, "ELEMENT_VELOCITY", 3,
			                     str_format, offset);
		if (output_data)
		{
			if (!useBinary)
			{
				fin << "          ";
				static double ele_vel[3] = { 0.0, 0.0, 0.0 };
				for (long i = 0; i < (long) msh->ele_vector.size(); i++)
				{
					ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
					fin << ele_vel[0] << " ";
                    if (!isXZplane) {
                        fin << ele_vel[1] << " ";
                        fin << ele_vel[2] << " ";
                    } else {
                        fin << ele_vel[2] << " ";
                        fin << ele_vel[1] << " ";
                    }
				}
				fin << "\n";
			}
			else
			{
				static double ele_vel[3] = { 0.0, 0.0, 0.0 };
				write_value_binary<unsigned int> (fin, sizeof(double) * 3 *
				                                  (long )msh->ele_vector.size()); //OK411
				for(long i = 0; i < (long)msh->ele_vector.size(); i++)
				{
					ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
                    write_value_binary(fin, ele_vel[0]);
                    if (!isXZplane) {
                        write_value_binary(fin, ele_vel[1]);
                        write_value_binary(fin, ele_vel[2]);
                    } else {
                        write_value_binary(fin, ele_vel[2]);
                        write_value_binary(fin, ele_vel[1]);
                    }
				}
			}
		}
		else
			//OK411
			offset += (long)msh->ele_vector.size() * sizeof(double) * 3 +
			          SIZE_OF_BLOCK_LENGTH_TAG;
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);

		//    if(out->pcs_type_name.compare("FLUID_MOMENTUM")==0)
		if(out->getProcessType () == FiniteElement::FLUID_MOMENTUM)
		{
			if (!useBinary || !output_data)
				WriteDataArrayHeader(fin,
				                     this->type_Double,
				                     "GLOBAL_VELOCITY",
				                     3,
				                     str_format,
				                     offset);
			if (output_data)
			{
				CRFProcess* pch_pcs = PCSGet("FLUID_MOMENTUM");
				if (!this->useBinary)
				{
					fin << "          ";
					for(long i = 0; i < (long)msh->ele_vector.size(); i++)
					{
						fin << pch_pcs->GetElementValue(
						        i,
						        pch_pcs->
						        GetElementValueIndex("VELOCITY1_X") +
						        1) << " ";
                        if (!isXZplane) {
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Y") +
                                1) << " ";
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Z") +
                                1) << " ";
                        } else {
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Z") +
                                1) << " ";
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Y") +
                                1) << " ";
                        }
					}
					fin << "\n";
				}
				else
				{
					//OK411
					write_value_binary<unsigned int>(fin, sizeof(double) * 3 *
					                                 (long)msh->ele_vector.size());
					for(long i = 0; i < (long)msh->ele_vector.size(); i++)
					{
						write_value_binary(fin,
						                   pch_pcs->GetElementValue(i,
						                                            pch_pcs
						                                            ->
						                                            GetElementValueIndex(
						                                                    "VELOCITY1_X")
						                                            + 1));
                        if (!isXZplane) {
                            write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Y")
                                + 1));
                            write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Z")
                                + 1));
                        } else {
                            write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Z")
                                + 1));
                        }
					}
				}
			}
			else
				//OK411
				offset += (long)msh->ele_vector.size() * sizeof(double) * 3 +
				          SIZE_OF_BLOCK_LENGTH_TAG;
			if (!useBinary || !output_data)
				WriteDataArrayFooter(fin);
		}
	}
	//Material information
	if(mmp_vector.size() > 1)
	{
		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
		if (output_data)
		{
			if (!this->useBinary)
			{
				fin << "          ";
				for(long i = 0; i < (long)msh->ele_vector.size(); i++)
				{
					ele = msh->ele_vector[i];
					fin << ele->GetPatchIndex() << " ";
				}
				fin << "\n";
			}
			else
			{
				//OK411
				write_value_binary<unsigned int>(
				        fin,
				        sizeof(int) *
				        (long)msh->ele_vector.size());
				for (long i = 0; i < (long)msh->ele_vector.size(); i++)
					write_value_binary(fin, msh->ele_vector[i]->GetPatchIndex());
			}
		}
		else
			//OK411
			offset += (long)msh->ele_vector.size() * sizeof(int) +
			          SIZE_OF_BLOCK_LENGTH_TAG;
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}
	return true;
}
