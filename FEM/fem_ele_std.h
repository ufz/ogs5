/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   Class element declaration
   class for finite element.
   Designed and programmed by WW/OK, 12/2004
   modified by TF, 10/2010
*/

#ifndef fem_std_INC
#define fem_std_INC

//#include "FEMEnums.h"

#include "fem_ele.h"
#include "matrix_class.h"

// Problems
#include "rf_mfp_new.h"
//#include "rf_msp_new.h"
#include "rf_out_new.h"                           //OK

#include "Eigen/Eigen"

//-----------------------------------------------------
// Process type
//L: Liquid flow
//U: Unconfined flow
//G: Groundwater flow
//T: Two-phase flow
//C: Componental flow
//H: heat transport
//M: Mass transport
//O: Overland flow
//R: Richards flow
//F: Fluid momentum
//A: Gas flow
//N: Thermal nonequilibrium
enum EnumProcessType { EPT_LIQUID_FLOW, EPT_UNCONFINED_FLOW, EPT_GROUNDWATER_FLOW, EPT_TWOPHASE_FLOW, EPT_COMPONENTAL_FLOW,
                       EPT_HEAT_TRANSPORT, EPT_MASS_TRANSPORT, EPT_OVERLAND_FLOW, EPT_RICHARDS_FLOW, EPT_FLUID_MOMENTUM,
                       EPT_GAS_FLOW, EPT_MULTIPHASE_FLOW, EPT_PSGLOBAL, EPT_MULTI_COMPONENTIAL_FLOW,
                       EPT_THERMAL_NONEQUILIBRIUM, EPT_TES };
//-----------------------------------------------------

namespace process
{class CRFProcessDeformation;
}
// using SolidProp::CSolidProperties; // evil in header!
namespace SolidProp
{ class CSolidProperties;
}
// Predeclared classes  01/07, WW
class CMediumProperties;
class CFluidProperties;

class CRFProcess;
namespace FiniteElement
{
using Math_Group::SymMatrix;
using Math_Group::Matrix;
using Math_Group::DiagonalMatrix;
using Math_Group::Vec;
using process::CRFProcessDeformation;
using ::CRFProcess;

class CFiniteElementStd : public CElement
{
public:
	CFiniteElementStd(CRFProcess* Pcs, const int C_Sys_Flad, const int order = 1);
	~CFiniteElementStd();

	// Set material data
	void SetMaterial(const int phase = 0);
	// Set memory for local matrices
	void SetMemory();
	// Set variable  YD
	void SetVariable();

	// Set coupling information
	void ConfigureCoupling(CRFProcess* pcs, const int* Shift, bool dyn = false);

	// Element claculation
	// 1. Mass matrix
	void CalcMass();
	void CalcMass2();
	void CalcMassMCF();                   //AKS/NB
	void CalcMassTNEQ();                      //AKS/NB
	void CalcMassTES();                       //AKS/NB
	void CalcMassPSGLOBAL();              // PCH
	// 2. Lumped mass matrix
	void CalcLumpedMass();
	void CalcLumpedMass2();
	void CalcLumpedMassTES();
	void CalcLumpedMassMCF();  //AKS
	void CalcLumpedMassPSGLOBAL();        // PCH
	// 3. Laplace matrix
	void CalcLaplace();
	void CalcLaplaceMCF();//AKS
	// 4. Gravity term
	void CalcGravity();
	// 5. Strain coupling matrix
	void CalcStrainCoupling(int phase = 0);
	// 6. Thermal coupling
	void CalcRHS_by_ThermalDiffusion();
	// 7. Advection matrix
	void CalcAdvection();
	void CalcAdvectionMCF();
	void CalcAdvectionTNEQ();
	void CalcAdvectionTES();
	// 8. Storage matrix
	void CalcStorage();
	// 9. Content matrix
	void CalcContent();
	void CalcContentTNEQ(); //NW
	void CalcContentTES(); //NW
	//
	void CalcSatution();                  //WW
	//
#ifdef E_NORM
	//25.08.2008. WW
	void CalcEnergyNorm(double &err_norm0, double &err_normn);
	//25.09.2008. WW
	void CalcEnergyNorm_Dual(double &err_norm0, double &err_normn);
	//
#endif
	void CalcNodeMatParatemer();          //WW
	// Assembly
	void Assembly();
	void Assembly(int option, int dimension); // PCH for Fluid Momentum
	void Cal_Velocity();

	void CalcSolidDensityRate();                     //HS thermal storage application, calculate rho_s
	void Cal_VelocityMCF();//AKS
	void Cal_Velocity_2();                //CB this is to provide velocity only at the element center of gravity
	void Cal_GP_Velocity_FM(int* i_ind); //SB 4900 interpolate node velocities to Gauss point velocities
	                                      //BG
	std::string Cal_GP_Velocity_ECLIPSE(std::string tempstring,
	                                    bool output_average,
	                                    int phase_index,
	                                    std::string phase);
	//BG coupling to DuMux
	std::string Cal_GP_Velocity_DuMux(int* i_ind, CRFProcess* m_pcs, int phase_index);
	// BG, 04/2012: Provides the average element velocity over all gauss points
	double Get_Element_Velocity(int Index, CRFProcess* m_pcs, int phase_index, int dimension);
	// necessary for using precalculated density and viscosity BG, 11/2010
	double InterpolatePropertyToGausspoint(int GPIndex, CRFProcess* m_pcs, int Variableindex);
	//
	//OK
	void AssembleParabolicEquationRHSVector();

	// CVFEM functions for overland flow   JOD
	// to move
	void GetOverlandBasisFunctionMatrix_Line();
	// to move
	void GetOverlandBasisFunctionMatrix_Quad();
	void CalcOverlandCoefficients(double* head, double* axx, double* ayy, double* ast);
	void CalcOverlandCoefficientsLine(double* head, double* axx, double* ast );
	void CalcOverlandCoefficientsQuad(double* head, double* axx, double* ayy, double* ast );
	void CalcOverlandCoefficientsTri(double* head, double* axx, double* ayy, double* ast );
	void CalcOverlandNLTERMS(double* H, double* HaaOld, double* swval, double* swold);
	void CalcOverlandNLTERMSRills(double* H, double* HaaOld, double* swval, double* swold);
	void CalcOverlandNLTERMSChannel(double* H, double* HaaOld, double* swval, double* swold);
	void CalcOverlandCKWR(double* head, double* ckwr, int* iups);
	void CalcOverlandCKWRatNodes(int i, int j, double* head, double* ckwr, int* iups);
	void CalcOverlandResidual(double* head,
	                          double* swval,
	                          double* swold,
	                          double ast,
	                          double* residuall,
	                          double** amat);
	double CalcOverlandJacobiNodes(int i,
	                               int j,
	                               double* depth,
	                               double* depth_keep,
	                               double akrw,
	                               double axx,
	                               double ayy,
	                               double** amatij,
	                               double* sumjac);
	void CalcOverlandUpwindedCoefficients(double** amat, double* ckwr, double axx, double ayy);
	//
	//CB added by CB: 090507
	void UpwindAlphaMass(double* alpha);
	//CB added by CB: 090507
	void UpwindUnitCoord(int p, int point, int ind);
	int UpwindElement(int option, int phase); // PCH
	//CB added by CB: 090507
	void UpwindSummandMass(const int gp,
	                       int& gp_r,
	                       int& gp_s,
	                       int& gp_t,
	                       double* alpha,
	                       double* summand);
	//NW
	double CalcSUPGCoefficient(double* vel,int ip);
	//NW
	void CalcSUPGWeightingFunction(double* vel, int ip, double &tau, double* v_dN);
	//NW
	double CalcSUPGEffectiveElemenetLength(double* vel);
	// Gauss value
	void ExtropolateGauss(CRFProcess* m_pcs, const int idof);
	// Extrapolate reaction rates on TNEQ flow
	void ExtrapolateGauss_ReactRate_TNEQ_TES(CRFProcess *m_pcs);
	void UpdateSolidDensity(size_t elem_idx);       // HS
	// CB _ctx_ CB_merge_0513
	//void Set_ctx_(long ele_index, double val, int gaussp, int i_dim);
	//double Get_ctx_(long ele_index, int gaussp, int i_dim);

private:
	bool newton_raphson;                  //24.05.2007 WW
	long index;
	int dof_index;                        //24.02.2007 WW
	// Column index in the node value table
	int idx0, idx1, idxS, idxSn0, idxSn1, idx3, idxMCF[12];
		 int idx_t2_0, idx_t2_1, idx_x0, idx_x1;
	int idxp0,idxp1, idxp20, idxp21, idxt0, idxt1;
	int phase;
	int comp;                             // Component
	int LocalShift;                       // For RHS
	// Danymic
	int* idx_vel_disp, idx_pres;
	// Velocity
	int* idx_vel;                         //WW
	// Material properties
	double* mat;
	double MassMatrixElements[36];
	double AdvectionMatrixElements[36];
	double ContentMatrixElements[36];
	double LaplaceMatrixElements[36][9];
	double* eqs_rhs;                      //For DDC WW
	bool heat_phase_change;

	//     /**
	//      * process type, \sa enum ProcessType
	//      */
	//     ProcessType _pcs_type; // TF

	bool dynamic;
	CRFProcess* mfp_pcs;
	SolidProp::CSolidProperties* SolidProp;
	CFluidProperties* FluidProp;
	CFluidProperties* GasProp;
	CMediumProperties* MediaProp;
	CMediumProperties* MediaProp1;        // Matrix for the dual model. YD/WW
	SolidProp::CSolidProperties* SolidProp1;         // Matrix for the dual model. YD/WW
	CRFProcess* pcs;
	::CRFProcess* cpl_pcs;                // Pointer to coupled process. WW
	process::CRFProcessDeformation* dm_pcs;
	bool flag_cpl_pcs;                    //OK
	//-------------------------------------------------------
	// Auxillarary matrices
	Matrix* StiffMatrix;
	Matrix* AuxMatrix;
	Matrix* AuxMatrix1;
	// Gravity matrix;
	//25.2.2007.WW  SymMatrix *GravityMatrix;
	// Gauss point value. Buffers. // Some changes. 27.2.2007 WW
	double TG, TG0, PG, PG0, PG2,PG20, drho_gw_dT;
	double Sw, rhow, poro, dSdp;
	double rho_gw, rho_ga, rho_g, p_gw, M_g, tort, Xw, eos_arg[5], heat_capacity, heat_conductivity, viscosity;


	//
	double* edlluse;                      // WW edlluse[16]
	double* edttuse;                      // WW edlluse[16]

	// Local matrices
	Matrix* Mass;                         //CB symMatrix *Mass; // unsymmetric in case of upwinding
	Matrix* Mass2;
	Matrix* Laplace;
	Matrix* Advection;                    //SB4200
	Matrix* Storage;                      //SB4200
	Matrix* Content;                      //SB4209
	Matrix* StrainCoupling;
	Vec* RHS;
	DiagonalMatrix* FCT_MassL;            //NW
	//-------------------------------------------------------
	void SetHighOrderNodes();             // 25.2.2007 WW
	// Primary as water head
	bool HEAD_Flag;
	//
	void Config();
	//
	double CalCoefMass();
	// 25.2.2007 WW
	double CalCoefMass2(int dof_index);
	double CalCoefMasstneq(int dof_index);
	// 03.3.2009 PCH
	double CalCoefMassPSGLOBAL(int dof_index);
	void CalCoefLaplace(bool Gravity, int ip = 0);
	// 10 2008 PCH
	void CalCoefLaplaceMultiphase(int phase, int ip = 0);
	void CalCoefLaplace2(bool Gravity, int dof_index);
	void CalCoefLaplaceTNEQ(const int dof_index);
	void CalCoefLaplaceTES(const int dof_index);
	void CalCoefLaplacePSGLOBAL(bool Gravity, int dof_index);
	double CalCoefAdvection();        //SB4200 OK/CMCD
	//AKS/NB
	double CalCoefAdvectionTNEQ(const int dof_index);
	double CalCoefAdvectionTES(const int dof_index);

	double CalCoefMassTNEQ(const int dof_index);
	double CalCoefMassTES(const int dof_index);
	void CalCoefMassMCF();
	void CalCoefAdvectionMCF();
	void CalCoefLaplaceMCF(int ip);
	void CalcContentMCF();
	void CalCoefContentMCF();

	double CalCoefStorage();       //SB4200
	double CalCoefContent();
	double CalCoefContentTNEQ(const int dof_index); //NW
	double CalCoefContentTES(const int dof_index);
	double CalCoefStrainCouping(const int phase = 0);

	double  CalcCoefDualTransfer();
	// 27.2.2007 WW
	double CalCoef_RHS_T_MPhase(int dof_index);
	double CalCoef_RHS_TNEQ(const int dof_index);
	double CalCoef_RHS_TES(const int dof_index);
	// 27.2.2007 WW
	double CalCoef_RHS_M_MPhase(int dof_index);
	double CalCoef_RHS_PSGLOBAL(int dof_index);
	//  NB
	double CalCoef_RHS_T_PSGlobal(int dof_index);
	// 03.2007 PCH
	void CalCoef_RHS_Pc(int dof_index);
	//AKS
	double CalCoef_RHS_AIR_FLOW(int dof_index);
	//AKS
	double CalCoef_RHS_HEAT_TRANSPORT(int dof_index);
	//AKS
	double CalCoef_RHS_HEAT_TRANSPORT2(int dof_index);
	void CalNodalEnthalpy();

	void ComputeAdditionalJacobi_H2(); //WW
	void ComputeAdditionalJacobi_Richards(); //WW

	//-----------------------------------------------------
	// Process type
	//L: Liquid flow
	//U: Unconfined flow
	//G: Groundwater flow
	//T: Two-phase flow
	//C: Componental flow
	//H: heat transport
	//M: Mass transport
	//O: Liquid flow
	//R: Richards flow
	//A: Gas flow
	//F: Fluid Momentum
	//N: Thermal nonequilibrium
	EnumProcessType PcsType;
	//-----------------------------------------------------
	// Local Assembly
	// Assembly of parabolic equation
	void AssembleParabolicEquation();     //OK4104
	void AssembleMixedHyperbolicParabolicEquation();
	void AssembleParabolicEquationNewton();
	// JOD
	void AssembleParabolicEquationNewtonJacobian(double** jacob,
	                                             double* Haa,
	                                             double* HaaOld,
	                                             double axx,
	                                             double ayy,
	                                             double** amat,
	                                             double ast,
	                                             double* swold,
	                                             double* residuall,
	                                             int* iups);
	void Assemble_strainCPL(const int phase = 0);     // Assembly of strain coupling
	void Assemble_strainCPL_Matrix(const double fac, const int phase = 0);

	void AssembleMassMatrix(int option);  // PCH
	// Assembly of RHS by Darcy's gravity term
	void Assemble_Gravity();
		void Assemble_GravityMCF();//AKS
	void Assemble_Gravity_Multiphase();
	// Assembly of RHS by temperature for m-phase flow 27.2.2007 WW
	void Assemble_RHS_T_MPhaseFlow();
	// Assembly of RHS by deformation. 27.2.2007 WW
	void Assemble_RHS_M();
	void Assemble_RHS_Pc();               // 03.2009 PCH
	void Assemble_RHS_AIR_FLOW();         //AKS
	void Assemble_RHS_HEAT_TRANSPORT();   //AKS
	void Assemble_RHS_TNEQ();      //AKS
	void Assemble_RHS_TES();      //AKS
	void Assemble_RHS_HEAT_TRANSPORT2();  //AKS
	void Assemble_RHS_T_PSGlobal();       // Assembly of RHS by temperature for PSGlobal
	void AssembleRHS(int dimension);      // PCH
	void Assemble_RHS_LIQUIDFLOW();       //NW
	void Assemble_DualTransfer();
	bool check_matrices;                  //OK4104
	void AssembleRHSVector();             //OK
	void AssembleCapillaryEffect();       // PCH
	                                      // PCH for debugging
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	void add2GlobalMatrixII();
#else
	void add2GlobalMatrixII(const int block_cols = 2);            //WW. 06.2011
#endif
	void PrintTheSetOfElementMatrices(std::string mark);
	// Friend classes, 01/07, WW
	friend class ::CMediumProperties;
	friend class SolidProp::CSolidProperties;
	friend class ::CFluidProperties;
	// Friend functions. WW
	friend double ::MFPCalcFluidsHeatCapacity(CFiniteElementStd * assem);

	// Auxillarary vectors for node values
	// Vector of local node values, e.g. pressure, temperature.
	// Assume maximium element nodes is 20
	//double OldMatrix[64]; // For grid adapting
	double NodalValue[12][40];
	double* NodalVal;
	double* NodalVal0;                    //?? NodalValueSaturation, NodalValueTemperature; ...
	double* NodalVal1;
	double* NodalVal2;
	double* NodalVal3;
	double* NodalVal4;
	double *NodalVal5;
	double* NodalValC;
	double* NodalValC1;
	double* NodalVal_Sat;
	double* NodalVal_SatNW;
	double* NodalVal_p2;
	double* NodalVal_p20;                 //AKS
	double *NodalVal_t0;                     // for TEMPERATURE1
	double *NodalVal_t1;                     //AKS
	double *NodalVal_t2_0;                   // FOR TEMPERATURE2 previous time step
	double *NodalVal_t2_1;                   // for TEMPERATURE2 current time step
	double *NodalVal_X0;                     // for CONCENTRATION previous time step
	double *NodalVal_X1;                     // for CONCENTRATION current time step
	//
	double* weight_func;                  //NW
	void CalcFEM_FCT();                   //NW
	//
	friend class ::CRFProcess;
};

// Vector for storing element values WW
class ElementValue
{
public:
	ElementValue(CRFProcess* m_pcs, CElem* ele);
	~ElementValue();
	void getIPvalue_vec(const int IP, double* vec);
	//SB 09/2010
	void getIPvalue_vec_phase(const int IP, int phase, double* vec);
	void GetEleVelocity(double* vec);
	Matrix Velocity;

	// HS Thermal Storage parameters---------------
	// Array of parameters on each Gauss point
	double *rho_s_prev, *rho_s_curr;
	double *q_R;
	// End of Thermal Storage parameters---------------
#ifdef USE_TRANSPORT_FLUX
	Matrix TransportFlux;  // Fick or Fourier law  with dispersion      JOD 2014-11-10
#endif
private:
	// Friend class
	friend class ::CRFProcess;
	friend class FiniteElement::CFiniteElementStd;
	friend class ::COutput;               //OK
	// Process
	CRFProcess* pcs;
	// Data
	Matrix Velocity_g;
	// CB _ctx_ CB_merge_0513
	//Matrix _ctx_Gauss;
};


// controls if gas mass form is used in matrix elements of TNEQ process
#define GAS_MASS_FORM
#ifdef GAS_MASS_FORM
const bool GasMassForm = true;
#else
const bool GasMassForm = false;
#endif
#undef GAS_MASS_FORM

}                                                 // end namespace

/*------------------------------------------------------------------
   Finite element calculation for standard PDE.
   12.12.2004 WW
   ------------------------------------------------------------------*/
#endif
