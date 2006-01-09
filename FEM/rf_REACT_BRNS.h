/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef RF_REACT_BRNS_H
#define RF_REACT_BRNS_H

#include "problem.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"
#include "rfmat_cp.h"
#include <time.h>

#ifdef BRNS
#ifdef WIN32
#pragma once
extern "C" __declspec( dllimport )  void brnsIsAlive();
extern "C" __declspec( dllimport )  void invokebrns(double* theCurArray,
                                                    double* thePreArray,
                                                    double* outputArray,
                                                    int* sizeOfArray,
                                                    double* time_step,
                                                    int* boundary_flag,
                                                    int* returnValue,
                                                    double* pos_x,
                                                    double* pos_y,
                                                    double* pos_z,
                                                    double* porosity,
                                                    double* waterSaturation,
                                                    double* parameterVector);
#endif
#endif

class REACT_BRNS
{
public:
	REACT_BRNS(void);
	~REACT_BRNS(void);

	// flag, whether BRNS is initialized; 1-true; 0-false;
	bool init_flag;
	int flowflag;
	bool flag_update_porosity;

	// Data structure storing the Concentration;
	double* cur_ts_Conc;
	double* pre_ts_Conc;
	double* m_dC_Chem_delta;
	int* boundary_flag;

	// storing porosity
	// its dimension is the number of nodes
	double* m_porosity_Node;
	double* m_porosity_Elem;

	// array of return values of BRNS;
	int* rt_BRNS;

#ifdef GCC
	void* hDll;
//	void* hDll, * hDll_1, * hDll_2;
	typedef void (*LPFNDLLFUNC)(double*, double*, double*, int*, double*, int*, int*, double*,
	                            double*, double*, double*, double*, double*);
	LPFNDLLFUNC invokebrns;
#endif

#ifdef USE_MPI_BRNS
	// buffer for the MPI implementation
	double* pre_ts_Conc_buf;
	int* rt_BRNS_buf;
#endif

	// pointer to the PCS Class;
	CRFProcess* m_pcs, * this_pcs, * m_flow_pcs;
	Problem* myProblem;

	// pointer to MFP class
	CFluidProperties* m_FluidProp;

	// pointer to the MCP Class;
	CompProperties* m_cp;

	// Media property
	CMediumProperties* MediaProp;

	// number of Components passed to BRNS;
	int num_Comp;

	// number of nodes;
	long nNodes;
	long nElems;

	// This is just a test run of BRNS dll;
	void TestRUN(void);

	// The Run function
	void RUN(double time_step);

	// Get the number of Nodes in GeoSys;
	long GetNodesNumber(void);
	long GetElemNumber(void);

	// Get flow flag;
	int GetFlowType_MT(void);

	// Get the number of Components in GeoSys;
	int GetCompsNumber(void);

	// Initialize Data Structure;
	void InitBRNS(Problem* myProblem);

	// Data transfer btw GeoSys and BRNS;
	void GSRF2Buffer( long i );
	void Buffer2GSRF( long i );

	// BC node checking
	int IsThisPointBCIfYesStoreValue(int index, CRFProcess* m_pcs, double& value);

	// porosity setting function
	int SetPorosityValue_MT ( long ele_Index,  double m_porosity_Elem, int i_timestep );
	void GetFluidProperty_MT ( void );
	void ConvPorosityNodeValue2Elem( int i_timestep );

private:
	// For measuring the time spent in BRNS calls
	double timeSpentInBrnsCoupling;
	clock_t startTime;

#ifdef USE_MPI_BRNS
	// MPI Buffer Value Manipulation
	void GetBRNSResult_MPI(void);
	void CleanMPIBuffer(void);
#endif
};
#endif
