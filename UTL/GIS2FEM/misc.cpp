/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "misc.h"

#include <iomanip>

#include "msh_mesh.h"

using namespace std;

void writeOGSMesh(const MeshLib::CFEMesh &mesh, const std::string file_base_name)
{
    const string fname = file_base_name + "_mesh_converted_from_GIS_shape_data.msh";
	ofstream os(fname.c_str(), ios::trunc);
    os.setf(ios::scientific, ios::floatfield);
    setw(14);
    os.precision(10);

    os << "#FEM_MSH \n $PCS_TYPE\n NO_PCS\n $NODES" <<endl;
	const size_t nn = mesh.GetNodesNumber(false);
	os << nn <<endl;

	for(size_t i=0; i<nn; i++)
	{
        MeshLib::CNode *n_ptr = mesh.nod_vector[i];
		const double *xyz_ptr = n_ptr->getData();
		os<< n_ptr->GetIndex() <<" "<<xyz_ptr[0] <<" "<<xyz_ptr[1] <<" "<<xyz_ptr[2] <<"\n";
	}

    os << "$ELEMENTS" <<endl;
    const size_t ne = mesh.ele_vector.size();

	for(size_t i=0; i<ne; i++)
	{
       MeshLib::CElem *e_ptr = mesh.ele_vector[i];
	   e_ptr->WriteIndex(os);
    }

	os << "#STOP" <<endl;

	os.close();
}
