/**
 * \file LegacyVtkInterface.h
 * 05/04/2011 LB Initial implementation
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef LEGACYVTKINTERFACE_H
#define LEGACYVTKINTERFACE_H

#include <fstream>
#include <string>
#include <vector>

#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#include "petscksp.h"
#include "petscmat.h"
typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;
#endif

namespace MeshLib
{
class CFEMesh;
}
class CRFProcess;
class ProcessInfo;

/// @brief Writes a legacy ascii vtk file of a mesh.
// TODO decouple from COutput
class LegacyVtkInterface
{
public:
	LegacyVtkInterface(MeshLib::CFEMesh* mesh,
	                   std::vector<std::string>
	                       pointArrayNames,
	                   std::vector<std::string>
	                       cellArrayNames,
	                   std::vector<std::string>
	                       materialPropertyArrayNames,
	                   std::string meshTypeName,
	                   ProcessInfo* processInfo);
	virtual ~LegacyVtkInterface();

	void WriteDataVTK(int number, double simulation_time, std::string baseFilename) const;
#if defined(USE_PETSC)
	void WriteDataVTKPETSC(int number, double simulation_time, std::string baseFilename) const;
#endif
	double RoundDoubleVTK(double MyZahl);

protected:
	void WriteVTKHeader(std::fstream&, int, double) const;
	void WriteVTKPointData(std::fstream&) const;
	void WriteVTKCellData(std::fstream&) const;
	void WriteVTKDataArrays(std::fstream&) const;
	void WriteELEVelocity(std::fstream& vtk_file) const;
#if defined(USE_PETSC)
	void WriteVTKPointDataPETSC(PetscViewer) const;
	void WriteVTKCellDataPETSC(PetscViewer) const;
	void WriteVTKDataArraysPETSC(PetscViewer) const;
#endif

	void printScalarArray(std::string arrayName, std::fstream& vtk_file) const;

	// Copied from COutput
	CRFProcess* GetPCS_ELE(const std::string& var_name) const;

	MeshLib::CFEMesh* _mesh;
	std::string _processType;
	std::vector<std::string> _pointArrayNames;
	std::vector<std::string> _cellArrayNames;
	std::vector<std::string> _materialPropertyArrayNames;
	std::string _meshTypeName;
	ProcessInfo* _processInfo;
};

#endif // LEGACYVTKINTERFACE_H
