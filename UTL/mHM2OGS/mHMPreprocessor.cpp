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

#include "Base/FileTools.h"

#include "FEM/fem_ele.h"
#include "FEM/mathlib.h"
#include "FEM/ShapeFunctionPool.h"

namespace MeshLib
{
struct RasterDataGIS
{
	std::size_t nrows;  /// Rows of GIS cells.
	std::size_t ncols;  /// Columns of GIS cells.
	/// (x_0, y_0): coordinate of the left down corner
	double x0, y0;
	double csize;            /// cell size.
	double ndata_v;          /// invalid data.
	std::vector<double> zz;  /// of cells.
};

mHMPreprocessor::~mHMPreprocessor()
{
	if (_fem)
		delete _fem;
}

void mHMPreprocessor::transform_mHMData(const std::string& output_path)
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

	mark_mHM_Interface3D();
	// Compute shape functions
	// Check element types of meshes
	std::vector<MshElemType::type> elem_types;
	elem_types.reserve(MshElemType::NUM_ELEM_TYPES);
	for (std::size_t i = 0;
	     i < static_cast<std::size_t>(MshElemType::NUM_ELEM_TYPES);
	     i++)
	{
		elem_types.push_back(MshElemType::INVALID);
	}
	elem_types[static_cast<int>(MshElemType::QUAD) - 1] = MshElemType::QUAD;
	elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] =
	    MshElemType::TRIANGLE;

	_fem = new FiniteElement::CElement(GetCoordinateFlag());
	FiniteElement::ShapeFunctionPool* line_shapefunction_pool =
	    new FiniteElement::ShapeFunctionPool(elem_types, *_fem, 3);
	_fem->setShapeFunctionPool(line_shapefunction_pool,
	                           line_shapefunction_pool);

	std::string key, uname, ofname;

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

	std::string file_path = pathDirname(*_geo_name);
	if (!output_path.empty())
		file_path = output_path;
	while (!ins.eof())
	{
		getline(ins, aline);
		ss.str(aline);
		ss >> key;
		ss.clear();

		if (key.size() == 0)  // An empty line
			continue;

		if (key.find("#STOP") != std::string::npos)
			break;

		ofname = file_path + key + ".bin";
		infil << step << " " << key + ".bin"
		      << "\n";

		key = file_path + key;
		transfromSingle_mHMdataToNodalFlux(key, ofname, ratio);

		step += 1.0;
	}
	infil << "#STOP"
	      << "\n";
	infil.close();
}

//---------------------------------------------------------------------------
/*!
   \brief Read GIS shapfile

   \param fname The file name.

   03/2010 WW
 */
void ReadShapeFile(std::string const& fname, RasterDataGIS& raster_data)
{
	std::ifstream ins(fname.c_str());
	if (!ins.good())
	{
		std::cout << "Can not find file "
		          << "\n";
		return;
	}

	std::string aline;
	std::stringstream ss;
	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.ncols;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.nrows;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.x0;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.y0;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.csize;
	ss.clear();

	getline(ins, aline);
	ss.str(aline);
	ss >> aline >> raster_data.ndata_v;
	ss.clear();

	raster_data.zz.clear();
	raster_data.zz.resize(raster_data.nrows * raster_data.ncols);
	for (std::size_t i = 0; i < raster_data.zz.size(); i++)
		ins >> raster_data.zz[i];
	ins.close();
}

void mHMPreprocessor::transfromSingle_mHMdataToNodalFlux(
    std::string const& fname, std::string const& ofname, double ratio)
{
	std::vector<double> val;
	val.resize(NodesNumber_Linear);
	const std::size_t nod_vector_size(nod_vector.size());
	for (std::size_t i = 0; i < nod_vector_size; i++)
	{
		nod_vector[i]->SetMark(false);
		val[i] = 0.0;
	}

	//
	RasterDataGIS raster_data;
	ReadShapeFile(fname, raster_data);

	//
	std::ofstream ofile_bin(ofname.c_str(), std::ios::trunc | std::ios::binary);
	ofile_bin.setf(std::ios::scientific, std::ios::floatfield);
	ofile_bin.precision(14);

	double node_val[8];
	for (std::size_t i = 0; i < face_vector.size(); i++)
	{
		CElem* elem = face_vector[i];
		if (!elem->GetMark())
			continue;

		for (std::size_t k = 0; k < elem->GetNodesNumber(false); k++)
			node_val[k] = 0.0;

		for (std::size_t k = 0; k < elem->GetNodesNumber(false); k++)
		{
			CNode const* const node = elem->GetNode(k);
			double const* const pnt_k(node->getData());

			long nx = static_cast<long>((pnt_k[0] - raster_data.x0) /
			                            raster_data.csize);
			long ny = static_cast<long>((pnt_k[1] - raster_data.y0) /
			                            raster_data.csize);
			ny = raster_data.nrows - ny;
			if (ny < 0)
				ny = 0;
			if (ny > static_cast<long>(raster_data.nrows))
				ny = raster_data.nrows;

			if (nx * raster_data.csize + raster_data.x0 >= pnt_k[0])
				nx -= 1;
			if (ny * raster_data.csize + raster_data.y0 >= pnt_k[1])
				ny -= 1;
			if (nx >= static_cast<long>(raster_data.ncols) - 1)
				nx = raster_data.ncols - 2;
			if (ny >= static_cast<long>(raster_data.nrows) - 1)
				ny = raster_data.nrows - 2;
			if (nx < 0)
				nx = 0;
			if (ny < 0)
				ny = 0;

			node_val[k] = raster_data.zz[raster_data.ncols * ny + nx];
			if (fabs(node_val[k] - raster_data.ndata_v) <
			    std::numeric_limits<double>::min())
				node_val[k] = 0.;
		}

		elem->ComputeVolume();
		_fem->setOrder(getOrder() + 1);
		_fem->ConfigElement(elem);
		_fem->FaceIntegration(node_val);
		for (std::size_t k = 0; k < elem->GetNodesNumber(false); k++)
		{
			CNode* node = elem->GetNode(k);
			node->SetMark(true);
			val[node->GetIndex()] += node_val[k];
		}
	}

	std::size_t counter = 0;
	for (size_t i = 0; i < nod_vector_size; i++)
	{
		if (!nod_vector[i]->GetMark())
			continue;
		counter++;
		val[i] *= ratio;
	}
	ofile_bin.write((char*)(&counter), sizeof(counter));

	for (std::size_t i = 0; i < nod_vector_size; i++)
	{
		CNode* node = nod_vector[i];
		if (!node->GetMark())
			continue;
		std::size_t nx = node->GetIndex();
		ofile_bin.write((char*)(&nx), sizeof(nx));
		ofile_bin.write((char*)(&val[i]), sizeof(val[i]));
	}

	ofile_bin.close();
	val.clear();
}

}  // end of namespace MeshLib
