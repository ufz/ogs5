/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   mHMPreprocessor.cpp
 *  Created on May 31, 2018, 1:39 PM
 */

#include "mHMPreprocessor.h"

#include <limits>

#include "FEM/fem_ele.h"
#include "FEM/mathlib.h"
#include "FEM/ShapeFunctionPool.h"

namespace MeshLib
{
/*!
   Find element nodes on the top surface of a mesh domain
   07.06.2010
   By WW
 */
void mHMPreprocessor::MarkInterface_mHM_Hydro_3D()
{
#ifdef output_top_z
	/// For output z coordinate of all nodes on the top surface
	/// 13.08.2010. WW
	vector<bool> node_mark(NodesNumber_Linear);
	for (std::size_t i = 0; i < NodesNumber_Linear; i++)
		node_mark[i] = false;
#endif

	const double tol = sqrt(std::numeric_limits<double>::epsilon()); // 1.e-5;
	for (std::size_t i = 0; i < face_vector.size(); i++)
	{
		CElem* elem = face_vector[i];
		CElem const* const own_elem = elem->GetOwner();

		//// In element
		//		// compute center of mass
		double cent[] = {0., 0., 0.};
		for (std::size_t k = 0; k < own_elem->GetNodesNumber(false); k++)
		{
			//				node = own_elem->nodes[k];
			double const* const pnt_k(own_elem->GetNode(k)->getData());
			cent[0] += pnt_k[0];
			cent[1] += pnt_k[1];
			cent[2] += pnt_k[2];
		}
		for (int k = 0; k < 3; k++)
			cent[k] /= (double)own_elem->GetNodesNumber(false);

		//			node = elem->nodes[0];
		double const* const pnt_0(elem->GetNode(0)->getData());
		cent[0] -= pnt_0[0];
		cent[1] -= pnt_0[1];
		cent[2] -= pnt_0[2];
		NormalizeVector(cent, 3);
		elem->ComputeVolume();
		elem->FillTransformMatrix();
		/// Compute the normal to this surface element
		Math_Group::Matrix const* const transform_tensor = elem->getTransformTensor();
		double fac = cent[0] * (*transform_tensor)(0, 2) + cent[1] * (*transform_tensor)(1, 2)
		                   + cent[2] * (*transform_tensor)(2, 2);
		if (fac > 0.0)
			fac = -1.0;
		else
			fac = 1.0;
		//////

		/// If n.z>0
		if ((*transform_tensor)(2, 2) * fac > tol)
		{
			elem->SetMark(true);
			for (int k = 0; k < 3; k++)
				(*transform_tensor)(k, 2) *= fac;

#ifdef output_top_z
			for (std::size_t k = 0; k < elem->GetNodesNumber(false); k++)
				node_mark[elem->GetNodeIndex(k)] = true;
#endif
		}
		else if ((*transform_tensor)(2, 2) * fac < -tol)
		{
			elem->SetMark(false);
			for (int k = 0; k < 3; k++)
				(*transform_tensor)(k, 2) *= fac;

#ifdef output_top_z
			for (k = 0; k < elem->GetNodesNumber(quad); k++)
				node_mark[elem->GetNodeIndex(k)] = bottom;
#endif
		}
		else
			elem->SetMark(false);
	}

#ifdef output_top_z
	string ccc = FileName + "_top_head.asc";
	ofstream ofile_asci(ccc.c_str(), ios::trunc);

	for (std::size_t i = 0; i < (long)nod_vector.size(); i++)
	{
		node = nod_vector[i];
		if (!node_mark[i])
			continue;
		ofile_asci << node->GetIndex() << " " << node->Z() << "\n";
	}
	ofile_asci << "#STOP"
	           << "\n";
#endif
}

/*!
   \brief Transform GIS shapfile stored precitation data into finite element node values

   Assume the precipitation data is stored in GIS shapfiles. This funtion read the data
   and then performs the numerical integration for each  GIS shapfile.

   06/2010  WW

 */
void mHMPreprocessor::mHM2NeumannBC(const std::string output_path)
{
	double ratio, step;

	std::string aline;
	std::stringstream ss;

	std::string fname = *_geo_name + ".pcp";

	std::ifstream ins(fname.c_str());
	if (!ins.good())
	{
		std::cout << "Can not open file " << fname << "\n";
		return;
	}

	ConstructGrid();

	MarkInterface_mHM_Hydro_3D();

	std::string key, uname, ofname;
	// char stro[1024];

	getline(ins, aline);
	ss.str(aline);
	ss >> key >> uname;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> key >> ratio;
	ss.clear();

	step = 0.;

	std::string infiltration_files;
	infiltration_files = *_geo_name + ".ifl";
	std::ofstream infil(infiltration_files.c_str(), std::ios::trunc);

	std::basic_string<char>::size_type indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = _geo_name->find_last_of('\\');
	indexChLinux = _geo_name->find_last_of('/');
	//
	std::string file_path;
	if (indexChWin != std::string::npos)
		file_path = _geo_name->substr(0, indexChWin) + "\\";
	else if (indexChLinux != std::string::npos)
		file_path = _geo_name->substr(0, indexChLinux) + "/";

	while (!ins.eof())
	{
		getline(ins, aline);
		ss.str(aline);
		ss >> key;
		ss.clear();

		if (key.size() == 0) // An empty line
			continue;

		if (key.find("#STOP") != std::string::npos)
			break;

		// sprintf(stro, "%f",step);
		// ofname = stro;
		ofname = file_path + key + ".bin";
		infil << step << " " << key + ".bin"
		      << "\n";

		key = file_path + key;
		Precipitation2NeumannBC(key, ofname, ratio);

		step += 1.0;
	}
	infil << "#STOP"
	      << "\n";
	infil.close();
}

/*!
   Compute int {f} a dA on top surface.

   WW. 29.11.2010
 */
void mHMPreprocessor::TopSurfaceIntegration()
{
	ConstructGrid();
	MarkInterface_mHM_Hydro_3D();

	std::vector<double> val;
	val.resize(NodesNumber_Linear);
	for (std::size_t i = 0; i < nod_vector.size(); i++)
	{
		nod_vector[i]->SetMark(false);
		val[i] = 0.0;
	}

	std::string ofname = *_geo_name + "_top_surface_Neumann_BC.txt";
	std::ofstream ofile_asci(ofname.c_str(), std::ios::trunc);
	ofile_asci.setf(std::ios::scientific, std::ios::floatfield);
	ofile_asci.precision(14);

	// Compute shape functions
	// Check element types of meshes
	std::vector<MshElemType::type> elem_types;
	elem_types.reserve(MshElemType::NUM_ELEM_TYPES);
	for (std::size_t i = 0; i < static_cast<std::size_t>(MshElemType::NUM_ELEM_TYPES); i++)
	{
		elem_types.push_back(MshElemType::INVALID);
	}
	elem_types[static_cast<int>(MshElemType::QUAD) - 1] = MshElemType::QUAD;
	elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] = MshElemType::TRIANGLE;

	FiniteElement::CElement* fem = new FiniteElement::CElement(GetCoordinateFlag());
	FiniteElement::ShapeFunctionPool* line_shapefunction_pool
	    = new FiniteElement::ShapeFunctionPool(elem_types, *fem, 3);
	fem->setShapeFunctionPool(line_shapefunction_pool, line_shapefunction_pool);

	double node_val[8];
	for (std::size_t  i = 0; i < face_vector.size(); i++)
	{
		CElem* elem = face_vector[i];
		if (!elem->GetMark())
			continue;

		for (std::size_t  k = 0; k < elem->GetNodesNumber(false); k++)
			node_val[k] = 1.0;

		elem->ComputeVolume();
		fem->setOrder(getOrder() + 1);
		fem->ConfigElement(elem, true);
		fem->FaceIntegration(node_val);
		for (std::size_t  k = 0; k < elem->GetNodesNumber(false); k++)
		{
			CNode* node = elem->GetNode(k);
			node->SetMark(true);
			val[node->GetIndex()] += node_val[k];
		}
	}

	for (std::size_t  i = 0; i < nod_vector.size(); i++)
	{
		CNode* node = nod_vector[i];
		if (!node->GetMark())
			continue;
		ofile_asci << node->GetIndex() << " " << val[i] << "\n";
	}
	ofile_asci << "#STOP "
	           << "\n";

	ofile_asci.close();
	delete fem;
	fem = NULL;
	delete line_shapefunction_pool;
	line_shapefunction_pool = NULL;
	val.clear();
}

} // end of namespace MeshLib
