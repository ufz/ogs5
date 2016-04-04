/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Object: Process PCS_Deformation
   Task:
   Programing:
   07/2003 WW Implementation
**************************************************************************/
#ifndef pcs_dm_INC
#define pcs_dm_INC

#include <vector>

#include "rf_pcs.h"

// Strong discontinuity
extern bool Localizing; // for tracing localization
typedef struct
{
	int ElementIndex;
	int NumInterFace; // Number of intersection faces
	// Local indeces of intersection faces (3D)
	int* InterFace;
} DisElement;
extern std::vector<DisElement*> LastElement; // Last discontinuity element correponding to SeedElement
extern std::vector<long> ElementOnPath; // Element on the discontinuity path

namespace FiniteElement
{
class CFiniteElementVec;
}
using FiniteElement::CFiniteElementVec;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
class CPARDomain;
#endif

namespace process
{
enum InitDataReadWriteType
{
	none,
	read_write,
	read_all_binary,
	write_all_binary,
	read_all_asci,
	write_all_asci,
	read_stress_binary,
	write_stress_binary,
	read_displacement,
	write_displacement,
	read_pressure,
	write_pressure
};

// Elasto-plastic Deformation
class CRFProcessDeformation : public CRFProcess
{
public:
	CRFProcessDeformation();
	virtual ~CRFProcessDeformation();

	void Initialization();
	void InitialNodeValueHpcs(); // WX:08.2011
	void CalIniTotalStress(); // WX:04.2013

	// Assemble system equation
	void GlobalAssembly();
	void GlobalAssembly_DM();

	// overloaded
	double Execute(int loop_process_number);

	// Aux. Memory
	double* GetAuxArray() const { return ARRAY; }
	void ScalingNodeForce(const double SFactor);
	void InitGauss();
	//
	void SetInitialGuess_EQS_VEC();
	void UpdateIterativeStep(const double damp, const int u_type);
	void InitializeNewtonSteps(const bool ini_excav = false);
	double NormOfUpdatedNewton();
	void StoreLastSolution(const int ty = 0);
	void RecoverSolution(const int ty = 0);
	double NormOfDisp();
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS
	double NormOfUnkonwn_orRHS(bool isUnknowns = true);
#endif
	// Stress
	// For partitioned HM coupled scheme
	void ResetCouplingStep();
	void ResetTimeStep();
	//
	void UpdateStress();
	void UpdateInitialStress(bool ZeroInitialS);
	void Extropolation_GaussValue();

	// Excavation computation
	void ReleaseLoadingByExcavation();
	void CreateInitialState4Excavation();

	// Dynamic
	bool CalcBC_or_SecondaryVariable_Dynamics(bool BC = false);
	// Calculate scaling factor for load increment
	double CaclMaxiumLoadRatio();

	// Write stresses
	void WriteGaussPointStress(const bool last_step = false);
	void ReadGaussPointStress();
	void ReadElementStress();

	// Access members
	CFiniteElementVec* GetFEM_Assembler() const { return fem_dm; }
	// WX:07.2011
	void PostExcavation();
	// WX:10.2011
	void UpdateIniStateValue();

private:
	CFiniteElementVec* fem_dm;
	void InitialMBuffer();
	double* ARRAY;

	int counter;
	double InitialNorm;
	double InitialNormU;
	double InitialNormU0;

	InitDataReadWriteType idata_type;

	bool _has_initial_stress_data;

	//
	double error_k0;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	// Domain decompisition
	void DomainAssembly(CPARDomain* m_dom);
#endif

	// For strong discontinuity approach
	void Trace_Discontinuity();
	long MarkBifurcatedNeighbor(const int PathIndex);
};
} // end namespace

extern void CalStressInvariants(const long Node_Inex, double* StressInv);
// For visualization
extern void CalMaxiumStressInvariants(double* StressInv);
extern double LoadFactor;
extern double Tolerance_global_Newton;
extern double Tolerance_Local_Newton;
extern int enhanced_strain_dm;
extern int number_of_load_steps;
extern int problem_dimension_dm;
extern int PreLoad;
extern bool GravityForce;
#endif
