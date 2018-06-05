/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

////////////////////////////////////////////////////////////////
//   Filename: rf_fluid_momentum.h
//
//   written by PCH, 05/2005
////////////////////////////////////////////////////////////////

// C++ STL
#include "rf_pcs.h"
#include <fstream>
#include <list>
#include <string>
#include <vector>
//#include "rf_vel_new.h"
#include "fem_ele.h"
#include "fem_ele_std.h"
#ifndef NEW_EQS // WW. 06.11.2008
#include "matrix_routines.h"
#endif
#include "mathlib.h"

class PlaneSet
{
public:
	int eleIndex;
	double ratio; // ratio: contribution of velocity to this plane.

	double V[3];
	double norm[3];
	double Eele[3]; // The vector from the crossroad to the center of the connected element
	// that lies in one of the connected planes.

	// Constructor
	PlaneSet(void);

	PlaneSet& operator=(const PlaneSet& B)
	{
		eleIndex = B.eleIndex;
		ratio = B.ratio;
		for (int i = 0; i < 3; ++i)
		{
			V[i] = B.V[i];
			norm[i] = B.norm[i];
		}

		return *this;
	}
};

class CrossRoad
{
public:
	int numOfThePlanes;
	int Index; // This can be node or edge index
	// depending on crossroad or edge

	PlaneSet* plane;

	// Constructor and destructor
	CrossRoad(void);
	~CrossRoad(void);

	void CreatePlaneSet(const int index);

	// Some operator overloading
	CrossRoad& operator=(const CrossRoad& B)
	{
		Index = B.Index;
		numOfThePlanes = B.numOfThePlanes;
		plane = B.plane;

		return *this;
	}
};

class CFluidMomentum : public CRFProcess
{
public:
	CFluidMomentum(void);
	~CFluidMomentum(void);

	int RWPTSwitch;

	std::vector<CrossRoad*> crossroads;
	std::vector<CrossRoad*> joints;

	void Create(void);
	virtual double Execute(int loop_process_number);
	void SolveDarcyVelocityOnNode();
	void ConstructFractureNetworkTopology();
	void SolveForEdgeVelocity(void);

protected:
	FiniteElement::CFiniteElementStd* fem;

private:
	CRFProcess* m_pcs;
};

extern void FMRead(std::string pcs_name = "");
// WWextern void DATWriteHETFile(const char* file_name);
