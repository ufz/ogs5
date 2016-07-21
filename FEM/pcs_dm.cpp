/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pcs_dm.h"

#include "makros.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>

#include "StringTools.h"

#include "FEMEnums.h"
#include "mathlib.h"
//#include "femlib.h"
// Element
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
// BC_Dynamic
#include "rf_bc_new.h"
#include "rf_pcs.h" //OK_MOD"
#include "tools.h"
//
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
#include "matrix_routines.h"
#endif
#include "fem_ele_vec.h"
#include "rf_msp_new.h"
#include "rf_tim_new.h"
// Excavation
#include "rf_out_new.h"
#include "rf_st_new.h"
// GEOLib
#include "geo_sfc.h"
// MSHLib
#include "msh_elem.h"
// IC
#include "rf_ic_new.h"

#include "rf_node.h"

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
#include "PETSC/PETScLinearSolver.h"
#endif

using namespace std;

// Solver
#if defined(NEW_EQS)
#include "equation_class.h"
#endif

double LoadFactor = 1.0;
double Tolerance_global_Newton = 0.0;
double Tolerance_Local_Newton = 0.0;
int enhanced_strain_dm = 0;
int number_of_load_steps = 1;
int problem_dimension_dm = 0;
int PreLoad = 0;
bool GravityForce = true;

bool Localizing = false; // for tracing localization
// Last discontinuity element correponding to SeedElement
vector<DisElement*> LastElement(0);
vector<long> ElementOnPath(0); // Element on the discontinuity path

using namespace std;
using FiniteElement::CFiniteElementVec;
using FiniteElement::CFiniteElementStd;
using FiniteElement::ElementValue_DM;
using SolidProp::CSolidProperties;
using Math_Group::Matrix;

namespace process
{
CRFProcessDeformation::CRFProcessDeformation()
    : CRFProcess(), fem_dm(NULL), ARRAY(NULL), counter(0), InitialNorm(0.0), idata_type(none),
      _has_initial_stress_data(false), error_k0(1.0e10)

{
}

CRFProcessDeformation::~CRFProcessDeformation()
{
	const bool last_step = true;
	WriteGaussPointStress(last_step);
	if (type == 41 && (idata_type == write_all_binary || idata_type == read_write))
	{
		// mono-deformation-liquid
		WriteSolution();
	}

	if (ARRAY)
		delete[] ARRAY;
	if (fem_dm)
		delete fem_dm;

	fem_dm = NULL;
	ARRAY = NULL;
	// Release memory for element variables
	MeshLib::CElem* elem = NULL;
	for (std::size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			delete ele_value_dm[i];
			ele_value_dm[i] = NULL;
		}
	}
	if (enhanced_strain_dm > 0)
	{
		while (ele_value_dm.size() > 0)
			ele_value_dm.pop_back();
		for (std::size_t i = 0; i < LastElement.size(); i++)
		{
			DisElement* disEle = LastElement[i];
			delete disEle->InterFace;
			delete disEle;
			disEle = NULL;
		}
		while (LastElement.size() > 0)
			LastElement.pop_back();
	}
}

/*************************************************************************
   Task: Initilization for deformation process
   Programming:
   05/2003 OK/WW Implementation
   08/2003 WW   Some changes for monolithic scheme
   08/2004 WW   Changes based on PCSCreateMProcess(obsolete)
   last modified: WW
 **************************************************************************/
void CRFProcessDeformation::Initialization()
{
	//-- NW 25.10.2011
	// this section has to be executed at latest before calling InitGauss()
	// Control for reading and writing solution
	if (reload == 1)
		idata_type = write_all_binary;
	if (reload == 2)
		idata_type = read_all_binary;
	if (reload == 3)
		idata_type = read_write;

	// Local assembliers
	// An instaniate of CFiniteElementVec
	int i, Axisymm = 1; // ani-axisymmetry
	//
	if (m_msh->isAxisymmetry())
		Axisymm = -1; // Axisymmetry is true
	fem_dm = new CFiniteElementVec(this, Axisymm * m_msh->GetCoordinateFlag());
	fem_dm->SetGaussPointNumber(m_num->ele_gauss_points);
	//
	// Monolithic scheme
	if (type / 10 == 4)
		fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
	//
	pcs_number_deformation = pcs_number;
	//
	if (m_num)
	{
		Tolerance_Local_Newton = m_num->nls_plasticity_local_tolerance;
		Tolerance_global_Newton = m_num->nls_error_tolerance[0];
	}

	// Initialize material transform tensor for tansverse isotropic elasticity
	// UJG/WW. 25.11.2009
	for (i = 0; i < (int)msp_vector.size(); i++)
		msp_vector[i]->CalculateTransformMatrixFromNormalVector(problem_dimension_dm);

	if (!msp_vector.size())
	{
		std::cout << "***ERROR: MSP data not found!"
		          << "\n";
		return;
	}
	InitialMBuffer();
	////////////////////////////////////
	// WX:08.2011 initialise node value of h_pcs
	if (Neglect_H_ini == 2)
		InitialNodeValueHpcs();
	if (Neglect_H_ini == 1)
		CalIniTotalStress();

#ifdef DECOVALEX
	// DECOVALEX test
	size_t i;
	int idv0 = 0, idv1 = 0;
	CRFProcess* h_pcs = NULL;
	h_pcs = fem_dm->h_pcs;
	if (h_pcs->type == 14) // Richards
	{
		idv0 = h_pcs->GetNodeValueIndex("PRESSURE_I");
		idv1 = h_pcs->GetNodeValueIndex("PRESSURE1");
		for (i = 0; i < m_msh->GetNodesNumber(false); i++)
			h_pcs->SetNodeValue(i, idv0, h_pcs->GetNodeValue(i, idv1));
	}
#endif
	///////////////////////////
	if (fem_dm->dynamic)
		CalcBC_or_SecondaryVariable_Dynamics();

	// TEST
	//   De_ActivateElement(false);
}

/*************************************************************************
WX:08.2011 initialise node value of h pcs
*************************************************************************/
void CRFProcessDeformation::InitialNodeValueHpcs()
{
	CRFProcess* tmp_h_pcs = NULL;
	if (fem_dm->h_pcs == NULL)
		return;

	tmp_h_pcs = fem_dm->h_pcs;
	int h_pcs_type = tmp_h_pcs->type;
	int idv_p_ini, idv_p1_ini, idv_p2_ini, idv_p_1, idv_p1_1, idv_p2_1; // idv_sw_ini, idv_sw_1;
	size_t i;
	if (h_pcs_type == 1 || h_pcs_type == 41) // Liquide Flow
	{
		idv_p_ini = tmp_h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
		idv_p_1 = tmp_h_pcs->GetNodeValueIndex("PRESSURE1");
		for (i = 0; i < m_msh->GetNodesNumber(false); i++)
			tmp_h_pcs->SetNodeValue(i, idv_p_ini, tmp_h_pcs->GetNodeValue(i, idv_p_1));
	}
	else if (h_pcs_type == 14) // Richards Flow
	{
		idv_p_ini = tmp_h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
		// idv_sw_ini = tmp_h_pcs->GetNodeValueIndex("SATURATION1_Ini");
		idv_p_1 = tmp_h_pcs->GetNodeValueIndex("PRESSURE1");
		// idv_sw_1 = tmp_h_pcs->GetNodeValueIndex("SATURATION1");
		for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		{
			tmp_h_pcs->SetNodeValue(i, idv_p_ini, tmp_h_pcs->GetNodeValue(i, idv_p_1));
			// tmp_h_pcs->SetNodeValue(i,idv_sw_ini,tmp_h_pcs->GetNodeValue(i,idv_sw_1));
		}
	}
	else if (h_pcs_type == 1212 || h_pcs_type == 42) // Multi Phase Flwo
	{
		idv_p1_ini = tmp_h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
		idv_p2_ini = tmp_h_pcs->GetNodeValueIndex("PRESSURE2_Ini");
		// idv_sw_ini = tmp_h_pcs->GetNodeValueIndex("SATURSTION1_Ini");
		idv_p1_1 = tmp_h_pcs->GetNodeValueIndex("PRESSURE1");
		idv_p2_1 = tmp_h_pcs->GetNodeValueIndex("PRESSURE2");
		// idv_sw_1 = tmp_h_pcs->GetNodeValueIndex("SATURATION1");
		for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		{
			tmp_h_pcs->SetNodeValue(i, idv_p1_ini, tmp_h_pcs->GetNodeValue(i, idv_p1_1));
			tmp_h_pcs->SetNodeValue(i, idv_p2_ini, tmp_h_pcs->GetNodeValue(i, idv_p2_1));
			// tmp_h_pcs->SetNodeValue(i,idv_sw_ini,tmp_h_pcs->GetNodeValue(i,idv_sw_1));
		}
	}
	return;
}
/*************************************************************************
WX:04.2013 calculate initial total stress for neglect h initial effect
*************************************************************************/
void CRFProcessDeformation::CalIniTotalStress()
{
	if (fem_dm->h_pcs == NULL)
		return;
	CRFProcess* tmp_h_pcs = NULL;
	tmp_h_pcs = fem_dm->h_pcs;
	ElementValue_DM* eleV_DM = NULL;
	CSolidProperties* SMat = NULL;
	MeshLib::CElem* elem = NULL;
	int h_pcs_type = tmp_h_pcs->type;
	int idx_p1, idx_p2, idx_s;
	double pw = 0.;
	idx_p1 = tmp_h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
	idx_p2 = tmp_h_pcs->GetNodeValueIndex("PRESSURE2") + 1;
	idx_s = tmp_h_pcs->GetNodeValueIndex("SATURATION1") + 1;
	for (std::size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		pw = 0.;
		elem = m_msh->ele_vector[i];
		const int nnodes = elem->GetNodesNumber(false);
		const int MatGroup = elem->GetPatchIndex();
		SMat = msp_vector[MatGroup];
		for (int j = 0; j < nnodes; j++)
		{
			if (h_pcs_type == 1 || h_pcs_type == 41) // liquid flow
			{
				pw += tmp_h_pcs->GetNodeValue(j, idx_p1);
			}
			if (h_pcs_type == 14) // richards flow
			{
				if (SMat->bishop_model == 1)
					pw += tmp_h_pcs->GetNodeValue(j, idx_p1) * SMat->bishop_model_value;
				else if (SMat->bishop_model == 2)
					pw += tmp_h_pcs->GetNodeValue(j, idx_p1)
					      * pow(tmp_h_pcs->GetNodeValue(j, idx_s), SMat->bishop_model_value);
				else if (SMat->bishop_model == 3)
				{
					tmp_h_pcs->GetNodeValue(j, idx_s) >= SMat->bishop_model_value
					    ? pw += tmp_h_pcs->GetNodeValue(j, idx_p1)
					    : pw += 0.0;
				}
				else
					pw += tmp_h_pcs->GetNodeValue(j, idx_p1) * tmp_h_pcs->GetNodeValue(j, idx_s);
			}
			if (h_pcs_type == 1212 || h_pcs_type == 42) // multi phase flow  pg-sw*pc
			{
				if (SMat->bishop_model == 1)
					pw += tmp_h_pcs->GetNodeValue(j, idx_p2)
					      - SMat->bishop_model_value * tmp_h_pcs->GetNodeValue(j, idx_p1);
				else if (SMat->bishop_model == 2)
					pw += tmp_h_pcs->GetNodeValue(j, idx_p2)
					      - pow(tmp_h_pcs->GetNodeValue(j, idx_s), SMat->bishop_model_value)
					            * tmp_h_pcs->GetNodeValue(j, idx_p1);
				else if (SMat->bishop_model == 3)
				{
					tmp_h_pcs->GetNodeValue(j, idx_p1) >= SMat->bishop_model_value
					    ? pw += tmp_h_pcs->GetNodeValue(j, idx_p2) - tmp_h_pcs->GetNodeValue(j, idx_p1)
					    : pw += tmp_h_pcs->GetNodeValue(j, idx_p2);
				}
				else
					pw += tmp_h_pcs->GetNodeValue(j, idx_p2)
					      - tmp_h_pcs->GetNodeValue(j, idx_s) * tmp_h_pcs->GetNodeValue(j, idx_p1);
			}
		}
		pw /= nnodes; // average node value, could be also interp. value

		if (elem->GetMark()) // Marked for use, it is necessary here
		{
			// elem->SetOrder(true);
			// fem_dm->ConfigElement(elem);
			eleV_DM = ele_value_dm[i];
			const int NGS = fem_dm->GetNumGaussPoints();

			for (int gp = 0; gp < NGS; gp++)
			{
				for (int j = 0; j < 3; j++)
				{
					if (SMat->biot_const < 0 && pw < 0) // if biot is negative value, negative pw is not considered
						(*eleV_DM->Stress0)(j, gp) = (*eleV_DM->Stress0)(j, gp);
					else
						(*eleV_DM->Stress0)(j, gp) = (*eleV_DM->Stress0)(j, gp) - SMat->biot_const * pw;
				}
				// for postExcavation Stress0 is already updated with Stress, as well as UpdateIniStateValue()
			}
		}
	} // end loop elements
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitialMBuffer
   Task:  Initialize the temporarily used variables
   Programming:
   12/2003 WW
 **************************************************************************/
void CRFProcessDeformation::InitialMBuffer()
{
	if (!msp_vector.size())
	{
		cout << "No .msp file.   "
		     << "\n";
		abort();
	}

	size_t bufferSize(0);
	bool HM_Stagered = false;
	if (GetObjType() == 4)
	{
		bufferSize = GetPrimaryVNumber() * m_msh->GetNodesNumber(true);
		if (H_Process)
			HM_Stagered = true;
	}
	else if (GetObjType() == 41)
		bufferSize = (GetPrimaryVNumber() - 1) * m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false);
	else if (GetObjType() == 42)
		bufferSize = (GetPrimaryVNumber() - 2) * m_msh->GetNodesNumber(true) + 2 * m_msh->GetNodesNumber(false);

	// Allocate memory for  temporal array
	if (m_num->nls_method != 2)
		ARRAY = new double[bufferSize];

	// Allocate memory for element variables
	MeshLib::CElem* elem = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		//       if (elem->GetMark()) // Marked for use
		//       {
		ElementValue_DM* ele_val = new ElementValue_DM(elem, m_num->ele_gauss_points, HM_Stagered);
		ele_value_dm.push_back(ele_val);
		//       }
	}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  Solve plastic deformation by generalized Newton-Raphson method
   Programming:
   02/2003 OK Implementation
   05/2003 WW Polymorphism function by OK
   last modified: 23.05.2003
 **************************************************************************/
double CRFProcessDeformation::Execute(int loop_process_number)
{
#if defined(USE_MPI) || defined(USE_PETSC)
	if (myrank == 1)
	{
#endif
		std::cout << "\n      ================================================"
		          << "\n";
		std::cout << "    ->Process " << loop_process_number << ": " << convertProcessTypeToString(getProcessType())
		          << "\n";
		std::cout << "      ================================================"
		          << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
	}
#endif

	clock_t dm_time;

	//  const int MaxLoadsteps=10; //20;
	//  const double LoadAmplifier =2.0;
	//  double MaxLoadRatio=1.0, maxStress=0.0;
	//  double LoadFactor0 = 0.0;
	//  double minLoadRatio = 0.0000001;
	//  const int defaultSteps=100;

	double damping = 1.0;
	double NormU, Norm = 0.0, Error1, Error = 0.0;
	double ErrorU1, ErrorU = 0.0;

	// const int defaultSteps=100;
	int MaxIteration = m_num->nls_max_iterations;
	int elasticity = 0;
	// int monolithic=0;
	//
	string delim = " | ";
	//----------------------------------------------------------
	dm_time = -clock();
	//
	m_msh->SwitchOnQuadraticNodes(true);
	//
	// TEST if(num_type_name.find("EXCAVATION")!=0)
	if (hasAnyProcessDeactivatedSubdomains || NumDeactivated_SubDomains > 0
	    || num_type_name.find("EXCAVATION") != string::npos)
		// if(NumDeactivated_SubDomains>0||num_type_name.find("EXCAVATION")!=string::npos)
		CheckMarkedElement();
	// MarkNodesForGlobalAssembly();
	if (ExcavMaterialGroup > -1)
		CheckExcavedElement(); // WX:07.2011

	counter++; // Times of this method  to be called
	LoadFactor = 1.0;

	// For pure elesticity
	if (pcs_deformation <= 100 && !fem_dm->dynamic)
	{
		elasticity = 1;
		MaxIteration = 1;
	}

	// For monolithic scheme
	if (type / 10 == 4) // Modified at 05.07.2010 WW

		//        monolithic=1;
		number_of_load_steps = 1;
// system matrix
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
	eqs_new->Initialize();
#elif defined(NEW_EQS) // WW
//
#if defined(USE_MPI)
	CPARDomain* dom = dom_vector[myrank];
	long global_eqs_dim = pcs_number_of_primary_nvals * m_msh->GetNodesNumber(true);
	dom->ConfigEQS(m_num, global_eqs_dim, true);
#else
	eqs_new->ConfigNumerics(m_num); // 27.11.2007 WW
#endif
//
#else
	SetZeroLinearSolver(eqs);
#endif

	/*

	   // dom->test();

	   //TEST_MPI
	   string test = "rank";
	   char stro[1028];
	   sprintf(stro, "%d",myrank);
	   string test1 = test+(string)stro+"Assemble.txt";
	   ofstream Dum(test1.c_str(), ios::out); // WW
	   dom->eqsH->Write(Dum);
	   Dum.close();
	   //  MPI_Finalize();
	   exit(1);

	 */

	// JT//if(CouplingIterations == 0 && m_num->nls_method != 2)
	if (this->first_coupling_iteration && m_num->nls_method != 2)
		StoreLastSolution(); // u_n-->temp
	//  Reset stress for each coupling step when partitioned scheme is applied to HM
	if (H_Process && (type / 10 != 4))
		ResetCouplingStep();
	//
	// Compute the maxium ratio of load increment and
	//   predict the number of load steps
	// ---------------------------------------------------------------
	// Compute the ratio of the current load to initial yield load
	// ---------------------------------------------------------------
	number_of_load_steps = 1;
	/*
	   if(!elasticity&&!fem_dm->dynamic&&(type!=41))
	   {
	   InitializeNewtonSteps(0); // w=du=0
	   number_of_load_steps = 1;

	   if(counter==1) // The first time this method is called.
	   {
	       // This is may needed by pure mechacial process
	       //  Auto increment loading  1
	       //-----------  Predict the number of load increment -------------------
	   //  Prepared for the special usage. i.e, for test purpose.
	   //  Activate this segement if neccessary as well as Auto increment loading  2
	   //---------------------------------------------------------------------
	   // InitializeNewtonSteps(1); // u=0
	   DisplayMsgLn("\nEvaluate load ratio: ");
	   PreLoad = 1;
	   GlobalAssembly();
	   PreLoad = 0;
	   //		 {MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);  abort();}

	   ExecuteLinearSolver();

	   UpdateIterativeStep(1.0, 0); // w = w+dw
	   MaxLoadRatio=CaclMaxiumLoadRatio();

	   if(MaxLoadRatio<=1.0&&MaxLoadRatio>=0.0) number_of_load_steps=1;
	   else if(MaxLoadRatio<0.0)
	   {
	   if( fabs(maxStress)> MKleinsteZahl)
	   {
	   int Incre=(int)fabs(MaxLoadRatio/maxStress);
	   if(Incre)
	   number_of_load_steps=defaultSteps;
	   }
	   number_of_load_steps=defaultSteps;
	   }
	   else
	   {
	   number_of_load_steps=(int)(LoadAmplifier*MaxLoadRatio);
	   if(number_of_load_steps>=MaxLoadsteps) number_of_load_steps=MaxLoadsteps;
	   }
	   cout<<"\n***Load ratio: "<< MaxLoadRatio<<"\n";

	   }
	   }
	 */
	// ---------------------------------------------------------------
	// Load steps
	// ---------------------------------------------------------------
	/*
	   //TEST
	   if(Cam_Clay)
	   {
	    if(counter==1) number_of_load_steps = MaxLoadsteps;
	    else number_of_load_steps = 1;
	     if(fluid) number_of_load_steps = 10;
	   }
	 */
	for (int l = 1; l <= number_of_load_steps; l++)
	{
		// This is may needed by pure mechacial process
		//  Auto increment loading  2
		//-----------  Predict the load increment ration-------------------
		//     Prepared for the special usage. i.e, for test purpose.
		//     Agitate this segement if neccessary as well as Auto increment loading  1
		//---------------------------------------------------------------------
		/*
		   if(elasticity!=1&&number_of_load_steps>1)
		   {
		    // Predictor the size of the load increment, only perform at the first calling
		     minLoadRatio = (double)l/(double)number_of_load_steps;

		   // Caculate load ratio \gamma_k
		     if(l==1)
		     {
		       if(MaxLoadRatio>1.0)
		      {
		   LoadFactor=1.0/MaxLoadRatio;
		   if(LoadFactor<minLoadRatio) LoadFactor=minLoadRatio;
		   LoadFactor0=LoadFactor;
		   }
		   else LoadFactor = (double)l/(double)number_of_load_steps;
		   }
		   else
		   {
		   if(MaxLoadRatio>1.0)
		   LoadFactor+=(1.0-LoadFactor0)/(number_of_load_steps-1);
		   else LoadFactor = (double)l/(double)number_of_load_steps;
		   }
		   // End of Caculate load ratio \gamma_k

		   }
		   //TEST     else if(fluid)   LoadFactor = (double)l/(double)number_of_load_steps;
		 */

		//
		// Initialize inremental displacement: w=0
		InitializeNewtonSteps();
		//
		// Begin Newton-Raphson steps
		if (elasticity != 1)
		{
			// ite_steps = 0;
			Error = 1.0e+8;
			ErrorU = 1.0e+8;
			Norm = 1.0e+8;
			NormU = 1.0e+8;

#if defined(USE_MPI) || defined(USE_PETSC)
			if (myrank == 0)
			{
#endif
				// Screan printing:
				std::cout << "      Starting loading step " << l << "/" << number_of_load_steps
				          << ".  Load factor: " << LoadFactor << "\n";
				std::cout << "      ------------------------------------------------"
				          << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
			}
#endif
		}
		ite_steps = 0;
		while (ite_steps < MaxIteration)
		{
			ite_steps++;
// Refresh solver
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
//	   InitializeRHS_with_u0();
#elif defined(NEW_EQS) // WW
#ifndef USE_MPI
			eqs_new->Initialize(); // 27.11.2007 WW
#endif
#ifdef JFNK_H2M
			/// If JFNK method (1.09.2010. WW):
			if (m_num->nls_method == 2)
			{
				Jacobian_Multi_Vector_JFNK();
				eqs_new->setPCS(this);
			}
#endif
#else // ifdef NEW_EQS
			SetZeroLinearSolver(eqs);
#endif

// Assemble and solve system equation
/*
   #ifdef MFC
        CString m_str;
        m_str.Format("Time step: t=%e sec, %s, Load step: %i, NR-Iteration: %i, Calculate element matrices",\
                      aktuelle_zeit,pcs_type_name.c_str(),l,ite_steps);
        pWin->SendMessage(WM_SETMESSAGESTRING,0,(LPARAM)(LPCSTR)m_str);
   #endif
 */
#if defined(USE_MPI) || defined(USE_PETSC) // WW
			if (myrank == 0)
#endif
				std::cout << "      Assembling equation system..."
				          << "\n";

#if defined(USE_MPI) || defined(USE_PETSC) // WW
			clock_t cpu_time = 0; // WW
			cpu_time = -clock();
#endif
			if (m_num->nls_method != 2) // Not JFNK method. 05.08.2010. WW
				GlobalAssembly();
#if defined(USE_MPI) || defined(USE_PETSC) // WW
			cpu_time += clock();
			cpu_time_assembly += cpu_time;
#endif
			//
			if (type != 41)
#if defined(USE_PETSC) // WW
				InitializeRHS_with_u0();
#else
				SetInitialGuess_EQS_VEC();
#endif
#ifdef USE_MPI // WW
// No initial guess for deformation.
// for(long ll=0; ll<eqs->dim; ll++)
//  eqs->x[ll] = 0.0;
//
#endif

#if defined(USE_MPI) || defined(USE_PETSC) // WW
			if (myrank == 0)
#endif
				std::cout << "      Calling linear solver..."
				          << "\n";
/// Linear solver
#if defined(USE_PETSC) //|| defined(other parallel libs)//03~04.3012. WW
			eqs_new->Solver();
			eqs_new->MappingSolution();
			if (!elasticity)
				Norm = eqs_new->GetVecNormRHS();

#elif defined(NEW_EQS) // WW
//
#if defined(USE_MPI)
			// 21.12.2007
			dom->eqsH->Solver(eqs_new->x, global_eqs_dim);
#else
#if defined(LIS) || defined(MKL)
			eqs_new->Solver(this->m_num); // NW
#else
			eqs_new->Solver(); // 27.11.2007
#endif
#endif
#else // ifdef NEW_EQS
			ExecuteLinearSolver();
#endif
//
// Norm of b (RHS in eqs)
#ifdef USE_MPI
			if (!elasticity)
				Norm = dom->eqsH->NormRHS();
#else
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif defined(NEW_EQS)
			if (!elasticity)
				Norm = eqs_new->NormRHS();
#else
			if (!elasticity)
				Norm = NormOfUnkonwn_orRHS(false);
#endif
#endif
			if (!elasticity)
			{
				// Check the convergence
				Error1 = Error;
				ErrorU1 = ErrorU;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
				NormU = eqs_new->GetVecNormX();
#elif defined(NEW_EQS)
				NormU = eqs_new->NormX();
#else
				NormU = NormOfUnkonwn_orRHS();
#endif

				// JT//if(ite_steps == 1 && CouplingIterations == 0)
				if (ite_steps == 1 && this->first_coupling_iteration)
				{
					InitialNorm = Norm;
					InitialNormU0 = NormU;
					if (counter == 1)
						InitialNormU = NormU;
				}

				Error = Norm / InitialNorm;
				ErrorU = NormU / InitialNormU0;
//				if (Norm < Tolerance_global_Newton && Error > Norm)
//					Error = Norm;
//				//           if(Norm<TolNorm)  Error = 0.01*Tolerance_global_Newton;
//				if ((NormU / InitialNormU) <= Tolerance_global_Newton)
//					Error = NormU / InitialNormU;

				// Compute damping for Newton-Raphson step
				damping = 1.0;
				//           if(Error/Error1>1.0e-1) damping=0.5;
//				if (Error / Error1 > 1.0e-1 || ErrorU / ErrorU1 > 1.0e-1)
//					damping = 0.5;
//				if (ErrorU < Error)
//					Error = ErrorU;
#if defined(NEW_EQS) && defined(JFNK_H2M)
				/// If JFNK, get w from the buffer
				if (m_num->nls_method == 2)
					// damping = LineSearch();
					Recovery_du_JFNK();

#endif
				// JT: Store the process and coupling errors
				pcs_num_dof_errors = 1;
				if (ite_steps == 1)
				{
					pcs_absolute_error[0] = NormU;
					pcs_relative_error[0] = pcs_absolute_error[0] / Tolerance_global_Newton;
					cpl_max_relative_error = pcs_relative_error[0];
					cpl_num_dof_errors = 1;
				}
				else
				{
					pcs_absolute_error[0] = Error;
					pcs_relative_error[0] = Error / Tolerance_global_Newton;
				}
//
#if defined(USE_MPI) || defined(USE_PETSC)
				if (myrank == 0)
				{
#endif
				//Screan printing:
				std::cout<<"      -->End of Newton-Raphson iteration: "<<ite_steps<<"/"<< MaxIteration <<"\n";
				cout.width(8);
 				cout.precision(2);
				cout.setf(ios::scientific);
				cout<<"         NR-Error"<<"  "<<"RHS Norm 0"<<"  "<<"RHS Norm  "<<"  "<<"Unknowns Norm"<<"  "<<"Damping"<<"\n";
				cout<<"         "<<Error<<"  "<<InitialNorm<<"  "<<Norm<<"   "<<NormU<<"   "<<"   "<<damping<<"\n";
				std::cout <<"      ------------------------------------------------"<<"\n";
#if defined(USE_MPI) || defined(USE_PETSC)
				}
#endif
//				if (Error > 100.0 && ite_steps > 1)
//				{
//					printf("\n  Attention: Newton-Raphson step is diverged. Programme halt!\n");
//					exit(1);
//				}
//				if (InitialNorm < 10 * Tolerance_global_Newton)
//					break;
//				if (Norm < 0.001 * InitialNorm)
//					break;
				//Test on absolute and relative norms with same tolerance
				//TODO: Move to input file control
				//if((Error <= Tolerance_global_Newton && ErrorU <= Tolerance_global_Newton) || (Norm <= Tolerance_global_Newton && NormU <= Tolerance_global_Newton))
				if (NormU <= Tolerance_global_Newton) //This tests on displacement norm only
				{
					if (ite_steps == 1) // WX:05.2012
					{
						UpdateIterativeStep(damping, 0);
						break;
					}
					else
						break;
				}
			}
			// w = w+dw for Newton-Raphson
			UpdateIterativeStep(damping, 0); // w = w+dw
		} // Newton-Raphson iteration

		// Update stresses
		UpdateStress();
		if (fem_dm->dynamic)
			CalcBC_or_SecondaryVariable_Dynamics();

		// Update displacements, u=u+w for the Newton-Raphson
		// u1 = u0 for pure elasticity
		UpdateIterativeStep(1.0, 1);
	}
	// Load step
	//
	// For coupling control
	std::cout << "      Deformation process converged."
	          << "\n";
	Error = 0.0;
	if (type / 10 != 4) // Partitioned scheme
	{
#ifdef USE_PETSC
		NormU = eqs_new->GetVecNormX();
#else
		for (size_t n = 0; n < m_msh->GetNodesNumber(true); n++)
			for (int l = 0; l < pcs_number_of_primary_nvals; l++)
			{
				NormU = GetNodeValue(n, fem_dm->Idx_dm1[l]);
				Error += NormU * NormU;
			}
		NormU = Error;
		Error = sqrt(NormU);
#endif
	}

	// Determine the discontinuity surface if enhanced strain methods is on.

	if (enhanced_strain_dm > 0)
		Trace_Discontinuity();
	//
	dm_time += clock();
#if defined(USE_MPI) || defined(USE_PETSC) // WW
	if (myrank == 0)
	{
#endif
		std::cout <<"      CPU time elapsed in deformation: " << (double)dm_time / CLOCKS_PER_SEC<<"s"<<"\n";
		std::cout <<"      ------------------------------------------------"<<"\n";
#if defined( USE_MPI) || defined( USE_PETSC)                             //WW
	}
#endif
	// Recovery the old solution.  Temp --> u_n	for flow proccess
	if (m_num->nls_method != 2)
		RecoverSolution();
//
#ifdef NEW_EQS // WW
#if defined(USE_MPI)
	dom->eqsH->Clean();
#else
	// Also allocate temporary memory for linear solver. WW
	eqs_new->Clean();
#endif
#endif
	//
	// JT//if(CouplingIterations > 0)
	if (this->first_coupling_iteration)
		Error = fabs(Error - error_k0) / error_k0;
	error_k0 = Error;
	//

	/*  WX: Move this part into PostExcavation() in PostLoop()
	//----------------------------------------------------------------------
	//Excavation. .. .12.2009. WW
	//----------------------------------------------------------------------
	std::vector<int> deact_dom;
	for(l = 0; l < (long)msp_vector.size(); l++)
	    if(msp_vector[l]->excavated)
	        deact_dom.push_back(l);
	if(ExcavMaterialGroup >= 0 && PCS_ExcavState < 0) //WX:01.2010.update pcs excav state
	{
	    for(l = 0; l < (long)m_msh->ele_vector.size(); l++)
	        if((m_msh->ele_vector[l]->GetExcavState() > 0) &&
	           !(m_msh->ele_vector[l]->GetMark()))
	        {
	            PCS_ExcavState = 1;
	            break;
	        }
	}
	//WX:01.2011 modified for coupled excavation
	if(deact_dom.size() > 0 || PCS_ExcavState > 0)
	{
	    //	  MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);  //abort();}

	    //
	    // 07.04.2010 WW
	    size_t i;
	    bool done;
	    MeshLib::CElem* elem = NULL;
	    MeshLib::CNode* node = NULL;
	    ElementValue_DM* eleV_DM = NULL;
	    for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
	    {
	        eleV_DM = ele_value_dm[l];
	        (*eleV_DM->Stress0) =  (*eleV_DM->Stress);

	        elem = m_msh->ele_vector[l];
	        done = false;
	        for(i = 0; i < deact_dom.size(); i++)
	            if(elem->GetPatchIndex() == static_cast<size_t>(deact_dom[i]))
	            {
	                elem->MarkingAll(false);
	                done = true;
	                break;
	            }
	        if(ExcavMaterialGroup >= 0) //WX
	            if(elem->GetExcavState() >= 0)
	            {
	                elem->MarkingAll(false);
	                done = true;
	            }
	        if(done)
	            continue;
	        else
	            elem->MarkingAll(true);
	    }

	    size_t mesh_node_vector_size (m_msh->nod_vector.size());
	    for (size_t l = 0; l < mesh_node_vector_size; l++)
	        while(m_msh->nod_vector[l]->getConnectedElementIDs().size())
	            m_msh->nod_vector[l]->getConnectedElementIDs().pop_back();

	    size_t mesh_ele_vector_size (m_msh->ele_vector.size());
	    //WX:07.2011 error fixed
	    for (size_t l = 0; l < mesh_ele_vector_size; l++)
	    {
	        elem = m_msh->ele_vector[l];
	        if(!elem->GetMark())
	            continue;
	        for(i = 0; i < elem->GetNodesNumber(m_msh->getOrder()); i++)
	        {
	            done = false;
	            node = elem->GetNode(i);
	            for(size_t j = 0; j < node->getConnectedElementIDs().size(); j++)
	                if(l == node->getConnectedElementIDs()[j])
	                {
	                    done = true;
	                    break;
	                }
	            if(!done)
	                node->getConnectedElementIDs().push_back(l);
	        }
	    }                         //
	}
	*/

	return Error;
}

/**************************************************************************
   ROCKFLOW - Funktion: InitializeStress

   Aufgabe:
   Initilize all Gausss values and others

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   01/2003  WW  Erste Version
   09/2007  WW  Parallelize the released load for the excavation modeling
   letzte Aenderung:
**************************************************************************/
void CRFProcessDeformation::InitGauss(void)
{
	const int LenMat = 7;
	size_t i;
	int j, k, gp, NGS, MatGroup, n_dom;
	int PModel = 1;
	//  double z=0.0;
	double xyz[3];
	static double Strs[6];
	ElementValue_DM* eleV_DM = NULL;
	CSolidProperties* SMat = NULL;
	CInitialCondition* m_ic = NULL;
	std::vector<CInitialCondition*> stress_ic(6);

	// double M_cam = 0.0;
	double pc0 = 0.0;
	double OCR = 1.0;
	n_dom = k = 0;

	int Idx_Strain[9];

	int NS = 4;
	Idx_Strain[0] = GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = GetNodeValueIndex("STRAIN_XY");

	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Strain[4] = GetNodeValueIndex("STRAIN_XZ");
		Idx_Strain[5] = GetNodeValueIndex("STRAIN_YZ");
	}
	Idx_Strain[NS] = GetNodeValueIndex("STRAIN_PLS");

	for (j = 0; j < NS; j++)
		stress_ic[j] = NULL;
	for (j = 0; j < (long)ic_vector.size(); j++)
	{
		m_ic = ic_vector[j];
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XX)
			stress_ic[0] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_YY)
			stress_ic[1] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_ZZ)
			stress_ic[2] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XY)
			stress_ic[3] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XZ)
			stress_ic[4] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_YZ)
			stress_ic[5] = m_ic;
	}
	int ccounter = 0;
	for (j = 0; j < NS; j++)
		if (stress_ic[j])
			ccounter++;

	// Initial stresses are given in .ic file, reload file is therefore disabled,
	if (ccounter > 0)
	{
		_has_initial_stress_data = true;
	}

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		for (j = 0; j < NS + 1; j++)
			SetNodeValue(i, Idx_Strain[j], 0.0);
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			MatGroup = elem->GetPatchIndex();
			SMat = msp_vector[MatGroup];
			elem->SetOrder(true);
			eleV_DM = ele_value_dm[i];
			*(eleV_DM->Stress0) = 0.0;
			*(eleV_DM->Stress) = 0.0;
			PModel = SMat->Plasticity_type;

			for (j = 3; j < fem_dm->ns; j++)
				Strs[j] = 0.0;

			if (PModel == 2)
				*(eleV_DM->xi) = 0.0;

			if (PModel == 3)
			{
				// WW M_cam = (*SMat->data_Plasticity)(0);
				pc0 = (*SMat->data_Plasticity)(3); // The initial preconsolidation pressure
				// Void ratio
				*(eleV_DM->e_i) = (*SMat->data_Plasticity)(4);
				OCR = (*SMat->data_Plasticity)(5); // Over consolidation ratio
				for (j = 0; j < 3; j++)
					Strs[j] = (*SMat->data_Plasticity)(6 + j);

				/*
				   g_s = GetSolidDensity(i);
				   if(g_s<=0.0)
				   {
				   printf("\n !!! Input error. Gravity density should not be less than zero with Cam-Clay model\n  ");
				   abort();
				   }

				   if(EleType== TriEle) // Triangle
				   nh = 6;
				   // Set soil profile. Cam-Clay. Step 2
				   for (j = 0; j < nh; j++)
				   h_node[j]=GetNodeY(element_nodes[j]); //Note: for 3D, should be Z
				 */
			}

			fem_dm->ConfigElement(elem);
			fem_dm->setOrder(2);
			fem_dm->SetIntegrationPointNumber(elem->GetElementType());
			NGS = fem_dm->GetNumGaussPoints();
			//

			for (gp = 0; gp < NGS; gp++)
			{
				if (ccounter > 0)
				{
					fem_dm->getShapefunctValues(gp, 2);
					fem_dm->RealCoordinates(xyz);
					for (j = 0; j < NS; j++)
					{
						m_ic = stress_ic[j];
						if (!m_ic)
							continue;
						n_dom = m_ic->GetNumDom();
						for (k = 0; k < n_dom; k++)
						{
							if (MatGroup != m_ic->GetDomain(k))
								continue;
							(*eleV_DM->Stress)(j, gp) = m_ic->getLinearFunction()->getValue(k, xyz[0], xyz[1], xyz[2]);
							(*eleV_DM->Stress0)(j, gp) = (*eleV_DM->Stress)(j, gp);
						}
					}
				}
				else
				{
					switch (PModel)
					{
						case 2: // Weimar's model
							// Initial stress_xx, yy,zz
							for (j = 0; j < 3; j++)
								(*eleV_DM->Stress)(j, gp) = (*SMat->data_Plasticity)(20 + j);
							break;
						case 3: // Cam-Clay
							for (j = 0; j < 3; j++)
								(*eleV_DM->Stress0)(j, gp) = Strs[j];
							(*eleV_DM->Stress) = (*eleV_DM->Stress0);
							break;
					}
				}
				if (eleV_DM->Stress_j)
					(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
				//
				switch (PModel)
				{
					case 2: // Weimar's model
						for (j = 0; j < LenMat; j++)
							(*eleV_DM->MatP)(j, gp) = (*SMat->data_Plasticity)(j);
						break;
					case 3: // Cam-Clay
						pc0 *= OCR; /// TEST
						(*eleV_DM->prep0)(gp) = pc0;
						break;
				}
				//
			}
// Initial condition by LBNL
////////////////////////////////////////////////////////
//#define  EXCAVATION
#ifdef EXCAVATION
			int gp_r, gp_s, gp_t;
			double z = 0.0;
			double xyz[3];
			fem_dm->getShapeFunctionPtr(elem->GetElementType(), 1);

			for (gp = 0; gp < NGS; gp++)
			{
				fem_dm->GetGaussData(gp, gp_r, gp_s, gp_t);
				fem_dm->getShapefunctValues(gp, 2);
				fem_dm->RealCoordinates(xyz);
				/*
				   //THM2
				   z = 250.0-xyz[1];
				   (*eleV_DM->Stress)(1, gp) = -2360*9.81*z;
				   (*eleV_DM->Stress)(2, gp) = 0.5*(*eleV_DM->Stress)(1, gp);
				   (*eleV_DM->Stress)(0, gp) = 0.6*(*eleV_DM->Stress)(1, gp);
				 */

				// THM1
				z = 500 - xyz[2]; // 3D xyz[1]; //2D
				(*eleV_DM->Stress)(2, gp) = -(0.02 * z + 0.6) * 1.0e6;
				(*eleV_DM->Stress)(1, gp) = -2700 * 9.81 * z;
				(*eleV_DM->Stress)(0, gp) = -(0.055 * z + 4.6) * 1.0e6;

				if (eleV_DM->Stress_j)
					(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
			}
#endif
			////////////////////////////////////////////////////////
			elem->SetOrder(false);
		}
	}
	// Reload the stress results of the previous simulation
	if (idata_type == read_all_binary || idata_type == read_write)
	{
		ReadGaussPointStress();
		if (type == 41) // mono-deformation-liquid
			ReadSolution();
	}
	// For excavation simulation. Moved here on 05.09.2007 WW
	if (num_type_name.find("EXCAVATION") != 0)
		Extropolation_GaussValue();
	//
}
/*************************************************************************
   ROCKFLOW - Function: Calculations of initial stress and released load
   Programming:
   09/2007 WW
 **************************************************************************/
void CRFProcessDeformation::CreateInitialState4Excavation()
{
	size_t i;
	int j;
	int Idx_Strain[9];
	int NS = 4;
	if (num_type_name.find("EXCAVATION") != 0)
		return;
	//
	Idx_Strain[0] = GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = GetNodeValueIndex("STRAIN_XY");

	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Strain[4] = GetNodeValueIndex("STRAIN_XZ");
		Idx_Strain[5] = GetNodeValueIndex("STRAIN_YZ");
	}
	Idx_Strain[NS] = GetNodeValueIndex("STRAIN_PLS");
	// For excavation simulation. Moved here on 05.09.2007 WW
	if (!_has_initial_stress_data)
	{
		GravityForce = true;
		cout << "\n ***Excavation simulation: 1. Establish initial stress profile..."
		     << "\n";
		counter = 0;
		Execute(0);
	}
	else
		UpdateInitialStress(true); // s0 = 0
	//
	Extropolation_GaussValue();
	//
	cout << "\n ***Excavation simulation: 2. Excavating..."
	     << "\n";
	counter = 0;
	InitializeNewtonSteps(true);
	GravityForce = false;
	//
	ReleaseLoadingByExcavation();
	// GravityForce = true;
	UpdateInitialStress(false); // s-->s0
	m_msh->ConnectedElements2Node();
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		for (j = 0; j < NS + 1; j++)
			SetNodeValue(i, Idx_Strain[j], 0.0);
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	if (dom_vector.size() > 0)
	{
		bc_node_value_in_dom.clear();
		bc_local_index_in_dom.clear();
		rank_bc_node_value_in_dom.clear();
		st_node_value_in_dom.clear();
		st_local_index_in_dom.clear();
		rank_st_node_value_in_dom.clear();
		CountDoms2Nodes(this);
		SetBoundaryConditionSubDomain();
	}
#endif
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitializeStress_EachCouplingStep()
   Programming:
   12/2005 WW
 **************************************************************************/
void CRFProcessDeformation::ResetCouplingStep()
{
	long i, e;
	int j;
	long number_of_nodes;
	long shift = 0;
	ElementValue_DM* eleV_DM = NULL;
	for (e = 0; e < (long)m_msh->ele_vector.size(); e++)
		if (m_msh->ele_vector[e]->GetMark())
		{
			eleV_DM = ele_value_dm[e];
			eleV_DM->ResetStress(true);
		}
	shift = 0;
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			SetNodeValue(j, p_var_index[i], ARRAY[shift + j]);
		shift += number_of_nodes;
	}
}
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitializeStress_EachCouplingStep()
   Programming:
   12/2005 WW
 **************************************************************************/
void CRFProcessDeformation::ResetTimeStep()
{
	long e;
	ElementValue_DM* eleV_DM = NULL;
	for (e = 0; e < (long)m_msh->ele_vector.size(); e++)
		if (m_msh->ele_vector[e]->GetMark())
		{
			eleV_DM = ele_value_dm[e];
			eleV_DM->ResetStress(false);
		}
}

/*************************************************************************
   ROCKFLOW - Funktion: TransferNodeValuesToVectorLinearSolver

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Programmaenderungen:
   02/2000   OK   aus TransferNodeValuesToVectorLinearSolver abgeleitet
   07/2005   WW  aus  TransferNodeValuesToVectorLinearSolver(OK)
   11/2010   WW  Modification for H2M
*************************************************************************/
void CRFProcessDeformation::SetInitialGuess_EQS_VEC()
{
	int i;
	long j, v_idx = 0;
	long number_of_nodes;
	long shift = 0;
	double* eqs_x = NULL;
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
	eqs_x = eqs_new->x;
#else
	eqs_x = eqs->x;
#endif
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		v_idx = p_var_index[i];
		if (i < problem_dimension_dm)
		{
			v_idx--;
			for (j = 0; j < number_of_nodes; j++)
				eqs_x[shift + j] = GetNodeValue(j, v_idx);
		}
		else
			for (j = 0; j < number_of_nodes; j++)
				eqs_x[shift + j] = 0.;
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: UpdateIterativeStep(LINEAR_SOLVER * ls, const int Type)

   Aufgabe:
   Update solution in Newton-Raphson procedure

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver
   const double damp : damping for Newton-Raphson method
   const int type    : 0,  update w=w+dw
                       1,  update u=u+w

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
void CRFProcessDeformation::UpdateIterativeStep(const double damp, const int u_type)
{
	int i, j;
	long shift = 0;
	long number_of_nodes;
	int ColIndex = 0;
	double* eqs_x = NULL;

#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
	eqs_x = eqs_new->GetGlobalSolution();
#elif defined(NEW_EQS)
	eqs_x = eqs_new->x;
#else
	eqs_x = eqs->x;
#endif

	if (type == 41 && fem_dm->dynamic)
	{
#if defined(USE_PETSC) // || defined (other parallel solver lib). 06.2013 WW
		// const long size_q = num_nodes_p_var[0];
		const long size_l = num_nodes_p_var[pcs_number_of_primary_nvals - 1];
#endif
		for (i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			number_of_nodes = num_nodes_p_var[i];
			//
			if (u_type == 0)
			{
				ColIndex = p_var_index[i] - 1;
				for (j = 0; j < number_of_nodes; j++)
				{
#if defined(USE_PETSC) // || defined (other parallel solver lib). 06.2013 WW
					long eqs_r = 0;
					const long ja = m_msh->Eqs2Global_NodeIndex[j];

					if (j < size_l)
						eqs_r = pcs_number_of_primary_nvals * ja + i;
					else
					{
						if (i > problem_dimension_dm)
							continue;

						eqs_r = pcs_number_of_primary_nvals * size_l + problem_dimension_dm * (ja - size_l) + i;
					}

					const long eqs_row = eqs_r;
#else
					const long eqs_row = j + shift;
#endif
					SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[eqs_row] * damp);
				}
				shift += number_of_nodes;
			}
			else
			{
				ColIndex = p_var_index[i];
				for (j = 0; j < number_of_nodes; j++)
				{
					SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + GetNodeValue(j, ColIndex - 1));
				}
			}
		}
		return;
	}

	//
	for (i = 0; i < problem_dimension_dm; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		//
		ColIndex = p_var_index[i] - 1;
		///  Update Newton step: w = w+dw
		if (u_type == 0)
		{
			for (j = 0; j < number_of_nodes; j++)
			{
#if defined(USE_PETSC) // || defined (other parallel solver lib). 06.2013 WW
				const long eqs_row = problem_dimension_dm * m_msh->Eqs2Global_NodeIndex[j] + i;
#else
				const long eqs_row = j + shift;
#endif
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[eqs_row] * damp);
			}
			shift += number_of_nodes;
		}
		else
		{
			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, ColIndex + 1, GetNodeValue(j, ColIndex + 1) + GetNodeValue(j, ColIndex));
		}
	}

#if defined(USE_PETSC) // || defined (other parallel solver lib). 06.2013 WW
	// The node wise storage has already included the following stuff.
	return;
#else
	// if(type == 42&&m_num->nls_method>0)         //H2M, Newton-Raphson. 06.09.2010. WW
	if (type / 10 == 4) // H2M, HM. 28.09.2011. WW
	{
		/// $p_{n+1}=p_{n+1}+\Delta p$ is already performed when type = 0
		if (u_type == 1)
			return;

		for (i = problem_dimension_dm; i < pcs_number_of_primary_nvals; i++)
		{
			number_of_nodes = num_nodes_p_var[i];
			//
			ColIndex = p_var_index[i];

			for (j = 0; j < number_of_nodes; j++)
			{
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[j + shift] * damp);
			}
			shift += number_of_nodes;
		}
	}
#endif
}

/**************************************************************************
   ROCKFLOW - Funktion: InitializeNewtonSteps(LINEAR_SOLVER * ls)

   Aufgabe:
   Initialize the incremental unknows in Newton-Raphson procedure

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:
   LINEAR_SOLVER * ls: linear solver
   const int type    : 0,  update w=0 (u0=0)
                       1,  update u=0 (u1=0)

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
   06/2007   WW   Rewrite
**************************************************************************/
void CRFProcessDeformation::InitializeNewtonSteps(const bool ini_excav)
{
	long i, j;
	long number_of_nodes;
	int col0, Col = 0, start, end;
	//
	//
	start = 0;
	end = pcs_number_of_primary_nvals;
	//

	/// u_0 = 0
	if (type == 42) // H2M
		end = problem_dimension_dm;

	/// Dynamic: plus p_0 = 0
	if (type == 41 && !fem_dm->dynamic)
	{
		// p_1 = 0
		for (i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			Col = p_var_index[i];
			col0 = Col - 1;
			number_of_nodes = num_nodes_p_var[i];
			if (i < problem_dimension_dm)
				for (j = 0; j < number_of_nodes; j++)
					// SetNodeValue(j, Col, 0.0);
					SetNodeValue(j, col0, 0.0);

			else
			{
				if (m_num->nls_method > 0) // If newton. 29.09.2011. WW
					continue;

				for (j = 0; j < number_of_nodes; j++)
					SetNodeValue(j, Col, 0.0);
			}
		}
	}
	else // non HM monolithic
	{
		for (i = start; i < end; i++)
		{
			Col = p_var_index[i] - 1;
			number_of_nodes = num_nodes_p_var[i];
			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, Col, 0.0);

			if (fem_dm->dynamic)
				continue;
		}
	}
	/// Excavation: plus u_1 = 0;
	if (ini_excav)
		// p_1 = 0
		for (i = 0; i < problem_dimension_dm; i++)
		{
			Col = p_var_index[i];
			number_of_nodes = num_nodes_p_var[i];
			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, Col, 0.0);
		}
}

/**************************************************************************
   ROCKFLOW - Funktion: NormOfUpdatedNewton

   Aufgabe:
   Compute the norm of Newton increment

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   12/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
double CRFProcessDeformation::NormOfUpdatedNewton()
{
	int i, j;
	long number_of_nodes;
	double NormW = 0.0;
	double val;
	int Colshift = 1;
	//
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
		{
			val = GetNodeValue(j, p_var_index[i] - Colshift);
			NormW += val * val;
		}
	}
	return sqrt(NormW);
}

/**************************************************************************
   ROCKFLOW - Funktion: StoreDisplacement

   Aufgabe:
   Copy the displacement of the previous time interval to a vector
   temporarily

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/

void CRFProcessDeformation::StoreLastSolution(const int ty)
{
	int i, j;
	long number_of_nodes;
	long shift = 0;

	// Displacement
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			ARRAY[shift + j] = GetNodeValue(j, p_var_index[i] - ty);
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: RetrieveDisplacement(LINEAR_SOLVER * ls)

   Aufgabe:
   Retrive the displacement from the temporary array
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
void CRFProcessDeformation::RecoverSolution(const int ty)
{
	int i, j, idx;
	long number_of_nodes;
	int Colshift = 1;
	long shift = 0;
	double tem = 0.0;

	int start, end;

	start = 0;
	end = pcs_number_of_primary_nvals;

	// If monolithic scheme for p-u coupling,  p_i-->p_0 only
	if (pcs_deformation % 11 == 0 && ty > 0)
	{
		start = problem_dimension_dm;
		for (i = 0; i < start; i++)
			shift += num_nodes_p_var[i];

		// TODO: end = problem_dimension_dm;
	}
	for (i = start; i < end; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		idx = p_var_index[i] - Colshift;
		for (j = 0; j < number_of_nodes; j++)
		{
			if (ty < 2)
			{
				if (ty == 1)
					tem = GetNodeValue(j, idx);
				SetNodeValue(j, idx, ARRAY[shift + j]);
				if (ty == 1)
					ARRAY[shift + j] = tem;
			}
			else if (ty == 2)
			{
				tem = ARRAY[shift + j];
				ARRAY[shift + j] = GetNodeValue(j, idx);
				SetNodeValue(j, idx, tem);
			}
		}
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: NormOfDisp

   Aufgabe:
   Compute the norm of  u_{n+1}
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
double CRFProcessDeformation::NormOfDisp()
{
	int i, j;
	long number_of_nodes;
	double Norm1 = 0.0;
	//
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			Norm1 += GetNodeValue(j, p_var_index[i]) * GetNodeValue(j, p_var_index[i]);
	}
	return Norm1;
}

/**************************************************************************
   ROCKFLOW - Funktion: NormOfUnkonwn

   Aufgabe:
   Compute the norm of unkowns of a linear equation

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   07/2011   WW

**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC)
double CRFProcessDeformation::NormOfUnkonwn_orRHS(bool isUnknowns)
{
	int i, j;
	long number_of_nodes;
	long v_shift = 0;
	double NormW = 0.0;
	double val;

#ifdef G_DEBUG
	if (!eqs)
	{
		printf(" \n Warning: solver not defined, exit from loop_ww.cc");
		exit(1);
	}
#endif

	double* vec = NULL;
	if (isUnknowns)
		vec = eqs->x;
	else
		vec = eqs->b;

	int end = pcs_number_of_primary_nvals;
	if (fem_dm->dynamic)
		end = problem_dimension_dm;

	for (i = 0; i < end; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
		{
			val = vec[v_shift + j];
			NormW += val * val;
		}

		v_shift += number_of_nodes;
	}
	return sqrt(NormW);
}
#endif
/**************************************************************************
   ROCKFLOW - Funktion: MaxiumLoadRatio

   Aufgabe:
   Calculate the muxium effective stress, Smax, of all Gauss points.
   (For 2-D 9 nodes element only up to now). Then compute the maxium ration
   by:

   Smax/Y0(initial yield stress)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   01/2003  WW  Erste Version

   letzte Aenderung:

**************************************************************************/
//#define Modified_B_matrix
double CRFProcessDeformation::CaclMaxiumLoadRatio(void)
{
	double* dstrain;

	double S0 = 0.0, p0 = 0.0;
	double MaxS = 0.000001;
	double EffS = 0.0;

	ElementValue_DM* eleV_DM = NULL;
	CSolidProperties* SMat = NULL;

	// Weimar's model
	Matrix* Mat = NULL;
	double II = 0.0;
	double III = 0.0;

	double PRatio = 0.0;
	const double MaxR = 20.0;

	// gp_t = 0;

	for (std::size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		MeshLib::CElem* elem = m_msh->ele_vector[i];
		if (elem->GetMark())      // Marked for use
		{
			fem_dm->ConfigElement(elem);
			fem_dm->SetMaterial();
			eleV_DM = ele_value_dm[i];
			SMat = fem_dm->smat;
			SMat->axisymmetry = m_msh->isAxisymmetry();
			const int PModel = SMat->Plasticity_type;

			//
			switch (PModel)
			{
				case 1:
#ifdef RFW_FRACTURE
					SMat->Calculate_Lame_Constant(elem);
#endif
#ifndef RFW_FRACTURE
					SMat->Calculate_Lame_Constant();
#endif
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
					SMat->CalulateCoefficent_DP();
					S0 = MSqrt2Over3 * SMat->BetaN * SMat->Y0;
					break;
				case 2:
#ifdef RFW_FRACTURE
					SMat->Calculate_Lame_Constant(elem);
#endif
#ifndef RFW_FRACTURE
					SMat->Calculate_Lame_Constant();
#endif
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
					Mat = eleV_DM->MatP;
					break;
				case 3:
					Mat = SMat->data_Plasticity;
					S0 = (*Mat)(3);
					break;
			}
			const int NGS = fem_dm->GetNumGaussPoints();
			//
			for (int gp = 0; gp < NGS; gp++)
			{
				if (!(    elem->GetElementType() == MshElemType::TRIANGLE
				       || elem->GetElementType() == MshElemType::QUAD) )
				{
					std::cerr << "CRFProcessDeformation::CaclMaxiumLoadRatio MshElemType not handled" << std::endl;
				}
				fem_dm->getGradShapefunctValues(gp, 2);
				fem_dm->ComputeStrain(gp);

				dstrain = fem_dm->GetStrain();

				if (PModel == 3) // Cam-Clay
				{
					p0 = ((*eleV_DM->Stress)(0, gp) + (*eleV_DM->Stress)(1, gp) + (*eleV_DM->Stress)(2, gp)) / 3.0;
					// Swelling index: (*SMat->data_Plasticity)(2)
					if (fabs(p0) < MKleinsteZahl)
						// The initial preconsolidation pressure
						p0 = (*SMat->data_Plasticity)(3);

					SMat->K = (1.0 + (*eleV_DM->e_i)(gp)) * fabs(p0) / (*SMat->data_Plasticity)(2);
					SMat->G = 1.5 * SMat->K * (1 - 2.0 * SMat->PoissonRatio) / (1 + SMat->PoissonRatio);
					SMat->Lambda = SMat->K - 2.0 * SMat->G / 3.0;
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
				}

				// Stress of the previous time step
				for (int j = 0; j < fem_dm->ns; j++)
					fem_dm->dstress[j] = (*eleV_DM->Stress)(j, gp);

				// Compute try stress, stress incremental:
				fem_dm->De->multi(dstrain, fem_dm->dstress);

				p0 = DeviatoricStress(fem_dm->dstress) / 3.0;

				switch (PModel)
				{
					case 1: // Drucker-Prager model
						EffS = sqrt(TensorMutiplication2(fem_dm->dstress, fem_dm->dstress, fem_dm->Dim()))
						       + 3.0 * SMat->Al * p0;

						if (EffS > S0 && EffS > MaxS && fabs(S0) > MKleinsteZahl)
						{
							MaxS = EffS;
							PRatio = MaxS / S0;
						}
						break;

					case 2: // Single yield surface
						// Compute try stress, stress incremental:
						II = TensorMutiplication2(fem_dm->dstress, fem_dm->dstress, fem_dm->Dim());
						III = TensorMutiplication3(fem_dm->dstress, fem_dm->dstress, fem_dm->dstress, fem_dm->Dim());
						p0 *= 3.0;
						EffS
						    = sqrt(II * pow(1.0 + (*Mat)(5) * III / pow(II, 1.5), (*Mat)(6)) + 0.5 * (*Mat)(0) * p0 * p0
						           + (*Mat)(2) * (*Mat)(2) * p0 * p0 * p0 * p0) + (*Mat)(1) * p0
						      + (*Mat)(3) * p0* p0;

						if (EffS > (*Mat)(4))
						{
							if ((*Mat)(4) > 0.0)
							{
								if (EffS > MaxS)
									MaxS = EffS;
								PRatio = MaxS / (*Mat)(4);
								if (PRatio > MaxR)
									PRatio = MaxR;
							}
							else
								PRatio = EffS;
						}
						break;

					case 3: // Cam-Clay
						II = 1.5 * TensorMutiplication2(fem_dm->dstress, fem_dm->dstress, fem_dm->Dim());
						if (S0 > 0.0)
						{
							EffS = II / (p0 * (*Mat)(0) * (*Mat)(0)) + p0;
							if (EffS > S0)
								PRatio = EffS / S0;
							else
								PRatio = 1.0;
						}
						else
							PRatio = 1.0;
						break;
				}
			}
		}
	}

	return PRatio;
}

/**************************************************************************
   ROCKFLOW - Funktion: Extropolation_GaussValue

   Aufgabe:
   Calculate the stresses of element nodes using the values at Gauss points.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2002  WW  Erste Version
   07/2003  WW  Extroplolation in quadraitc triangle element is added
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::Extropolation_GaussValue()
{
	int k, NS;
	long i = 0;
	int Idx_Stress[7];
	const long LowOrderNodes = m_msh->GetNodesNumber(false);
	MeshLib::CElem* elem = NULL;

	// Clean nodal stresses
	NS = 4;
	Idx_Stress[0] = GetNodeValueIndex("STRESS_XX");
	if (m_msh->isAxisymmetry()) // indices are swapped if axisymmetry!
	{
		Idx_Stress[2] = GetNodeValueIndex("STRESS_YY");
		Idx_Stress[1] = GetNodeValueIndex("STRESS_ZZ");
	}
	else
	{
		Idx_Stress[1] = GetNodeValueIndex("STRESS_YY");
		Idx_Stress[2] = GetNodeValueIndex("STRESS_ZZ");
	}
	Idx_Stress[3] = GetNodeValueIndex("STRESS_XY");
	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Stress[4] = GetNodeValueIndex("STRESS_XZ");
		Idx_Stress[5] = GetNodeValueIndex("STRESS_YZ");
	}
	Idx_Stress[NS] = GetNodeValueIndex("STRAIN_PLS");
	NS++;

	// NB, TN
	int stressPrincipleIndices[3];

	stressPrincipleIndices[0] = GetNodeValueIndex("STRESS_1");
	stressPrincipleIndices[1] = GetNodeValueIndex("STRESS_2");
	stressPrincipleIndices[2] = GetNodeValueIndex("STRESS_3");

	int PrinStressDirectionIndices[9];

	PrinStressDirectionIndices[0] = GetNodeValueIndex("NORM_STRESS_1_X");
	PrinStressDirectionIndices[1] = GetNodeValueIndex("NORM_STRESS_1_Y");
	PrinStressDirectionIndices[2] = GetNodeValueIndex("NORM_STRESS_1_Z");
	PrinStressDirectionIndices[3] = GetNodeValueIndex("NORM_STRESS_2_X");
	PrinStressDirectionIndices[4] = GetNodeValueIndex("NORM_STRESS_2_Y");
	PrinStressDirectionIndices[5] = GetNodeValueIndex("NORM_STRESS_2_Z");
	PrinStressDirectionIndices[6] = GetNodeValueIndex("NORM_STRESS_3_X");
	PrinStressDirectionIndices[7] = GetNodeValueIndex("NORM_STRESS_3_Y");
	PrinStressDirectionIndices[8] = GetNodeValueIndex("NORM_STRESS_3_Z");

	for (i = 0; i < LowOrderNodes; i++)
	{
		for (k = 0; k < NS; k++)
			SetNodeValue(i, Idx_Stress[k], 0.0);
		for (k = 0; k < 3; k++)
			SetNodeValue(i, stressPrincipleIndices[k], 0.0);
	}

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			fem_dm->setElement(elem);
			fem_dm->setOrder(2);
			fem_dm->SetIntegrationPointNumber(elem->GetElementType());

			fem_dm->SetMaterial();
			//         eval_DM = ele_value_dm[i];
			// TEST        (*eval_DM->Stress) += (*eval_DM->Stress0);
			fem_dm->ExtropolateGuassStress();
			// TEST        if(!update)
			//           (*eval_DM->Stress) -= (*eval_DM->Stress0);
			// calculation of principal stresses for postprocessing on nodes

			for (int i = 0; i < fem_dm->nnodes; i++)
			{
				// MeshLib::CNode *node;
				int node_index = fem_dm->nodes[i];
				// node = m_msh->nod_vector[node_index];
				double stress[6];
				stress[0] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[0]); // sigma_xx
				stress[1] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[1]); // sigma_yy
				stress[2] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[2]); // sigma_zz
				stress[3] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[3]); // sigma_xy
				if (problem_dimension_dm == 3)
				{
					stress[4] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[4]); // sigma_xz
					stress[5] = fem_dm->pcs->GetNodeValue(node_index, Idx_Stress[5]); // sigma_yz
				}
				else
				{
					stress[4] = 0; // sigma_xz
					stress[5] = 0; // sigma_yz
				}
				double prin_str[3];
				double prin_dir[9];
				fem_dm->smat->CalPrinStrDir(stress, prin_str, prin_dir, 3);
				// transpose rotation tensor for principal directions
				for (size_t i = 0; i < 3; i++)
					fem_dm->pcs->SetNodeValue(node_index, stressPrincipleIndices[i], prin_str[i]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[0], prin_dir[0]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[1], prin_dir[3]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[2], prin_dir[6]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[3], prin_dir[1]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[4], prin_dir[4]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[5], prin_dir[7]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[6], prin_dir[2]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[7], prin_dir[5]);
				fem_dm->pcs->SetNodeValue(node_index, PrinStressDirectionIndices[8], prin_dir[8]);
			}
		}
	}
}

/*--------------------------------------------------------------------------
   Trace discontinuity path. Belong to Geometry
   --------------------------------------------------------------------------*/
void CRFProcessDeformation::Trace_Discontinuity()
{
	long k, l;
	int i, nn, Size, bFaces, bFacesCounter, intP;
	int b_node_counter;
	int locEleFound, numf, nPathNodes;
	bool thisLoop; //, neighborSeed;

	// Element value
	ElementValue_DM* eleV_DM;

	static double xn[20], yn[20], zn[20];
	static double xa[3], xb[3], n[3], ts[3];
	double v1, v2;

	int FNodes0[8];
	std::vector<long> SeedElement;

	locEleFound = 0;
	nPathNodes = 2; // 2D element
	bFaces = 0;

	intP = 0;

	// Check all element for bifurcation
	MeshLib::CElem* elem = NULL;
	for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
	{
		elem = m_msh->ele_vector[l];
		if (elem->GetMark()) // Marked for use

			if (fem_dm->LocalAssembly_CheckLocalization(elem))
			{
				locEleFound++;
				// If this is first bifurcated element, call them as seeds
				if (!Localizing)
					SeedElement.push_back(l);
			}
	}

	if (locEleFound > 0 && !Localizing) // Bifurcation inception
	{
		// TEST
		// mesh1   de =23;
		// mesh2_iregular de = 76
		// mesh coarst de = 5;
		// mesh quad de=23
		// crack tri de=0

		/*
		   SeedElement.clear();
		   int de = 39; //64; //itri  //39; //Quad //72; //rtri crack
		   SeedElement.push_back(de);
		 */

		// TEST

		// Determine the seed element
		Size = (long)SeedElement.size();
		for (l = 0; l < Size; l++)
		{
			k = SeedElement[l];
			elem = m_msh->ele_vector[k];

			numf = elem->GetFacesNumber();
			eleV_DM = ele_value_dm[k];

			// If seed element are neighbor. Choose one
			/*//TEST
			   neighborSeed = false;
			   for(m=0; m<Size; m++)
			   {
			   if(m==l) continue;
			   for(i=0; i<numf; i++)
			   {
			       if(neighbor[i]==SeedElement[m])
			       {
			          neighborSeed = true;
			          break;
			   }
			   }
			   }
			   if(neighborSeed)
			   {
			   delete eleV_DM->orientation;
			   eleV_DM->orientation = NULL;
			   continue;
			   }*/
			//

			nn = elem->GetNodesNumber(true);
			for (i = 0; i < nn; i++)
			{
				// Coordinates of all element nodes
				//               xn[i] = elem->nodes[i]->X();
				//               yn[i] = elem->nodes[i]->Y();
				//               zn[i] = elem->nodes[i]->Z();
				double const* const coords(elem->nodes[i]->getData());
				xn[i] = coords[0];
				yn[i] = coords[1];
				zn[i] = coords[2];
			}
			// Elements which have one only boundary face are chosen as seed element
			bFaces = -1;
			bFacesCounter = 0;
			for (i = 0; i < numf; i++)
				if (elem->neighbors[i]->GetDimension() != elem->GetDimension())
				{
					bFaces = i;
					bFacesCounter++;
				}

			// Elements which have only one boundary face or one boundary node are chosen as seed element
			if (bFacesCounter != 1)
			{
				//
				b_node_counter = 0;
				for (i = 0; i < elem->GetVertexNumber(); i++)
				{
					bFaces = i;
					b_node_counter++;
				}
				if (b_node_counter != 1)
				{
					eleV_DM->Localized = false;
					delete eleV_DM->orientation;
					eleV_DM->orientation = NULL;
					continue;
				}
			}

			fem_dm->ConfigElement(elem);
			// 2D
			elem->GetElementFaceNodes(bFaces, FNodes0);
			if (elem->GetElementType() == MshElemType::QUAD || elem->GetElementType() == MshElemType::TRIANGLE)
				nPathNodes = 2;
			// Locate memory for points on the path of this element
			eleV_DM->NodesOnPath = new Matrix(3, nPathNodes);
			*eleV_DM->NodesOnPath = 0.0;
			if (nPathNodes == 2) // 2D
			{
				// Departure point
				(*eleV_DM->NodesOnPath)(0, 0) = 0.5 * (xn[FNodes0[0]] + xn[FNodes0[1]]);
				(*eleV_DM->NodesOnPath)(1, 0) = 0.5 * (yn[FNodes0[0]] + yn[FNodes0[1]]);
				(*eleV_DM->NodesOnPath)(2, 0) = 0.5 * (zn[FNodes0[0]] + zn[FNodes0[1]]);

				xa[0] = (*eleV_DM->NodesOnPath)(0, 0);
				xa[1] = (*eleV_DM->NodesOnPath)(1, 0);
				xa[2] = (*eleV_DM->NodesOnPath)(2, 0);

				// Check oreintation again.
				ts[0] = xn[FNodes0[1]] - xn[FNodes0[0]];
				ts[1] = yn[FNodes0[1]] - yn[FNodes0[0]];
				// ts[2] = zn[FNodes0[1]]-zn[FNodes0[0]];
				v1 = sqrt(ts[0] * ts[0] + ts[1] * ts[1]);
				ts[0] /= v1;
				ts[1] /= v1;

				n[0] = cos(eleV_DM->orientation[0]);
				n[1] = sin(eleV_DM->orientation[0]);
				v1 = n[0] * ts[0] + n[1] * ts[1];
				n[0] = cos(eleV_DM->orientation[1]);
				n[1] = sin(eleV_DM->orientation[1]);
				v2 = n[0] * ts[0] + n[1] * ts[1];
				if (fabs(v2) > fabs(v1))
				{
					v1 = eleV_DM->orientation[0];
					eleV_DM->orientation[0] = eleV_DM->orientation[1];
					eleV_DM->orientation[1] = v1;
				}

				intP = fem_dm->IntersectionPoint(bFaces, xa, xb);

				(*eleV_DM->NodesOnPath)(0, 1) = xb[0];
				(*eleV_DM->NodesOnPath)(1, 1) = xb[1];
				(*eleV_DM->NodesOnPath)(2, 1) = xb[2];
			}

			// Last element to this seed
			DisElement* disEle = new DisElement;
			disEle->NumInterFace = 1;
			disEle->ElementIndex = k;
			disEle->InterFace = new int[1];
			disEle->InterFace[0] = intP;
			LastElement.push_back(disEle);
			ElementOnPath.push_back(k);
		}

		Localizing = true;
	}

	// Seek path from the last bifurcated element of corrsponding seeds.
	if (Localizing)
	{
		Size = (long)LastElement.size();
		for (i = 0; i < Size; i++)
		{
			thisLoop = true;
			while (thisLoop)
			{
				nn = MarkBifurcatedNeighbor(i);
				if (nn < 0)
					break;
				ElementOnPath.push_back(nn);
			}
		}

		Size = (long)ElementOnPath.size();
		for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
		{
			elem = m_msh->ele_vector[l];
			if (elem->GetMark()) // Marked for use
			{
				eleV_DM = ele_value_dm[l];
				if (eleV_DM->Localized)
				{
					thisLoop = false;
					for (k = 0; k < Size; k++)
						if (l == ElementOnPath[k])
						{
							thisLoop = true;
							break;
						}
					if (!thisLoop)
					{
						eleV_DM->Localized = false;
						delete eleV_DM->orientation;
						eleV_DM->orientation = NULL;
					}
				}
			}
		}
	}
}

// WW
long CRFProcessDeformation::MarkBifurcatedNeighbor(const int PathIndex)
{
	int j;
	int f1, f2, nb, numf1;
	long index, Extended;
	bool adjacent, onPath;
	ElementValue_DM *eleV_DM, *eleV_DM1;
	DisElement* disEle;
	static double n1[2], n2[2], xA[3], xB[3];
	// WW static int Face_node[8];                    // Only 2D
	MeshLib::CElem* elem;
	MeshLib::CElem* elem1;

	double pd1, pd2;

	Extended = -1;
	// 2D only
	disEle = LastElement[PathIndex];
	index = disEle->ElementIndex;
	f1 = disEle->InterFace[0];

	elem = m_msh->ele_vector[index];
	eleV_DM = ele_value_dm[index];

	// numf = elem->GetFacesNumber();

	n1[0] = cos(eleV_DM->orientation[0]);
	n1[1] = sin(eleV_DM->orientation[0]);

	xA[0] = (*eleV_DM->NodesOnPath)(0, 1);
	xA[1] = (*eleV_DM->NodesOnPath)(1, 1);
	xA[2] = (*eleV_DM->NodesOnPath)(2, 1);

	// nfnode = elem->GetElementFaceNodes(f1, Face_node);

	// Check discintinuity path goes to which neighbor
	elem1 = elem->neighbors[f1];
	// Boundary reached
	if (elem1->GetDimension() != elem->GetDimension())
		return -1;
	nb = elem1->GetIndex();
	// Check if the element is already in discontinuity line/surface
	onPath = false;
	for (j = 0; j < (int)ElementOnPath.size(); j++)
		if (nb == ElementOnPath[j])
		{
			onPath = true;
			break;
		}

	// If has neighbor and it is not on the discontinuity surface.
	if (!onPath) // Has neighbor
		if (ele_value_dm[nb]->Localized)
		{
			// TEST OUT
			// cout <<" element on track  " <<nb<<"\n";

			adjacent = false;
			numf1 = elem1->GetFacesNumber();
			eleV_DM1 = ele_value_dm[nb];
			fem_dm->ConfigElement(elem1);
			// Search faces of neighbor's neighbors
			for (j = 0; j < numf1; j++)
			{
				// Neighbor is a face on surface
				if (elem1->neighbors[j]->GetDimension() != elem1->GetDimension())
					continue;
				if ((size_t)index != elem1->neighbors[j]->GetIndex())
					continue;
				{
					adjacent = true;
					Extended = nb;
					// Choose a smooth direction
					n2[0] = cos(eleV_DM1->orientation[0]);
					n2[1] = sin(eleV_DM1->orientation[0]);
					pd1 = n1[0] * n2[0] + n1[1] * n2[1];
					n2[0] = cos(eleV_DM1->orientation[1]);
					n2[1] = sin(eleV_DM1->orientation[1]);
					pd2 = n1[0] * n2[0] + n1[1] * n2[1];
					if (pd2 > pd1) // Always use the first entry of orientation
					{
						// Swap the values
						pd1 = eleV_DM1->orientation[1];
						eleV_DM1->orientation[1] = eleV_DM1->orientation[0];
						eleV_DM1->orientation[0] = pd1;
					}
					eleV_DM1->NodesOnPath = new Matrix(3, 2);
					*eleV_DM1->NodesOnPath = 0.0;
					// Get another intersection point
					f2 = fem_dm->IntersectionPoint(j, xA, xB);
					(*eleV_DM1->NodesOnPath)(0, 0) = xA[0];
					(*eleV_DM1->NodesOnPath)(1, 0) = xA[1];
					(*eleV_DM1->NodesOnPath)(2, 0) = xA[2];
					(*eleV_DM1->NodesOnPath)(0, 1) = xB[0];
					(*eleV_DM1->NodesOnPath)(1, 1) = xB[1];
					(*eleV_DM1->NodesOnPath)(2, 1) = xB[2];

					// Last element
					disEle->ElementIndex = nb;
					disEle->NumInterFace = 1;
					disEle->InterFace[0] = f2;
					break;
				}
			}

			// If not on the discontinuity surface
			// Release memory
			if (!adjacent)
			{
				delete eleV_DM1->orientation;
				eleV_DM1->orientation = NULL;
				eleV_DM1->Localized = false;
			}
		}

	// If true, discontinuity extended to its neighbor
	return Extended;
}

/**************************************************************************
   FEMLib-Method:
   Task: Assembly in the sense of sub-domains
   Programing:
   04/2006 WW
**************************************************************************/
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
void CRFProcessDeformation::DomainAssembly(CPARDomain* m_dom)
{
	long i;
	MeshLib::CElem* elem = NULL;
#ifdef NEW_EQS
	m_dom->InitialEQS(this);
#else
	SetLinearSolver(m_dom->eqs);
	SetZeroLinearSolver(m_dom->eqs);
#endif
	for (i = 0; i < (long)m_dom->elements.size(); i++)
	{
		elem = m_msh->ele_vector[m_dom->elements[i]];
		if (elem->GetMark()) // Marked for use
		{
			elem->SetOrder(true);
			// WW
			fem_dm->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
			fem_dm->ConfigElement(elem);
			fem_dm->m_dom = m_dom;
			fem_dm->LocalAssembly(0);
		}
	}
	if (type == 41) // p-u monolithic scheme
	{
		if (!fem_dm->dynamic)
			RecoverSolution(1); // p_i-->p_0
		// 2.
		// Assemble pressure eqs
		for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
		{
			elem = m_msh->ele_vector[m_dom->elements[i]];
			if (elem->GetMark()) // Marked for use
			{
				elem->SetOrder(false);
				// WW
				fem->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
				fem->ConfigElement(elem);
				fem->m_dom = m_dom;
				fem->Assembly();
			}
		}
		if (!fem_dm->dynamic)
			RecoverSolution(2); // p_i-->p_0
	}

	/*

	   //TEST
	   string test = "rank";
	   char stro[1028];
	   sprintf(stro, "%d",myrank);
	   string test1 = test+(string)stro+"dom_eqs.txt";

	   ofstream Dum(test1.c_str(), ios::out);
	   m_dom->eqsH->Write(Dum);
	   Dum.close();
	   exit(1);

	 */
}
#endif
/**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices and RHS for each element
   Programing:
   02/2005 WW
**************************************************************************/
void CRFProcessDeformation::GlobalAssembly()
{
#if !defined(USE_PETSC) // && !defined(other parallel libs)//10.3012. WW
#ifdef USE_MPI
	if (dom_vector.size() > 0)
	{
		std::cout << "      Domain Decomposition " << myrank << '\n';

		CPARDomain* m_dom = NULL;
		m_dom = dom_vector[myrank];
		DomainAssembly(m_dom);

		/*
		   //TEST
		   string test = "rank";
		   char stro[64];
		   sprintf(stro, "%d",myrank);
		   string test1 = test+(string)stro+"dom_eqs.txt";

		   ofstream Dum(test1.c_str(), ios::out);
		   m_dom->eqsH->Write(Dum);
		   Dum.close();
		   exit(1);
		 */

		// Apply Neumann BC
		IncorporateSourceTerms(myrank);
		// Apply Dirchlete bounday condition
		IncorporateBoundaryConditions(myrank);
		//....................................................................

		// Assemble global system
		// DDCAssembleGlobalMatrix();
		// MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);
	}
#else // ifdef USE_MPI
	//----------------------------------------------------------------------
	// DDC
	if (dom_vector.size() > 0)
	{
		cout << "      Domain Decomposition" << '\n';
		CPARDomain* m_dom = NULL;
		for (int j = 0; j < (int)dom_vector.size(); j++)
		{
			m_dom = dom_vector[j];
			DomainAssembly(m_dom);
			// Apply Neumann BC
			IncorporateSourceTerms(j);
			// Apply Dirchlete bounday condition
			IncorporateBoundaryConditions(j);
		}
		//....................................................................
		// Assemble global system
		DDCAssembleGlobalMatrix();

		//      ofstream Dum("rf_pcs.txt", ios::out); // WW
		//  eqs_new->Write(Dum);
		//  Dum.close();

		// MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); abort();
	}
#endif
	//----------------------------------------------------------------------
	// STD
	else
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//10.3012. WW
	{
		GlobalAssembly_DM();

		if (type / 10 == 4) // p-u monolithic scheme

			// if(!fem_dm->dynamic)   ///
			//  RecoverSolution(1);  // p_i-->p_0
			// 2.
			// Assemble pressure eqs
			// Changes for OpenMP
			GlobalAssembly_std(true);
		// if(!fem_dm->dynamic)
		//   RecoverSolution(2);  // p_i-->p_0

		//----------------------------------------------------------------------
		//
		// {        		MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); // abort();}

		// DumpEqs("rf_pcs1.txt");
		/*
		   ofstream Dum("rf_pcs_omp.txt", ios::out); // WW
		   eqs_new->Write(Dum);
		   Dum.close();
		 */

		// Apply Neumann BC
		IncorporateSourceTerms();

#if defined(USE_PETSC) // || defined(other parallel libs)//05.3013.
		eqs_new->AssembleRHS_PETSc();
		eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);
#endif

		// DumpEqs("rf_pcs2.txt");

		// {			MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); // abort();}

		/// If not JFNK or if JFNK but the Newton step is greater than one. 11.11.2010. WW
		if (!(m_num->nls_method == 2 && ite_steps == 1))
		{
			// Apply Dirchlete bounday condition
			if (!fem_dm->dynamic)
				IncorporateBoundaryConditions();
			else
				CalcBC_or_SecondaryVariable_Dynamics(true);
		}
		//  {  		MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x);  //abort();}
		//

#define atest_dump
#ifdef test_dump
		string fname = FileName + "rf_pcs_omp.txt";
		ofstream Dum1(fname.c_str(), ios::out); // WW
		eqs_new->Write(Dum1);
		Dum1.close(); //   abort();
#endif

#define atest_bin_dump
#ifdef test_bin_dump // WW
		string fname = FileName + ".eqiation_binary.bin";

		ofstream Dum1(fname.data(), ios::out | ios::binary | ios::trunc);
		if (Dum1.good())
			eqs_new->Write_BIN(Dum1);
		Dum1.close();
#endif
		//
	}
}

/*!  \brief Assembe matrix and vectors
      for deformation process

      24.11.2010. WW
 */
void CRFProcessDeformation::GlobalAssembly_DM()
{
	long i;
	MeshLib::CElem* elem = NULL;
	/// If JFNK method. 10.08.2010. WW
	//   if(m_num->nls_method==2&&ite_steps==1)
	//      IncorporateBoundaryConditions();

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark()) // Marked for use
			continue;

		elem->SetOrder(true);
		fem_dm->ConfigElement(elem);
		fem_dm->LocalAssembly(0);
	}
}

/**************************************************************************
   FEMLib-Method:
Task: post process for excavation
Programing:
07/2011 WX
**************************************************************************/
void CRFProcessDeformation::PostExcavation()
{
	std::vector<int> deact_dom;
	// PCS_ExcavState = -1;
	bool now_Excav = false;
	MeshLib::CElem* elem = NULL;
	MeshLib::CNode* node = NULL;
	for (size_t l = 0; l < msp_vector.size(); l++)
	{
		if (msp_vector[l]->excavated)
			deact_dom.push_back(l);
	}
	// if(ExcavMaterialGroup>=0&&PCS_ExcavState<0)	//WX:01.2010.update pcs excav state
	if (ExcavMaterialGroup >= 0)
	{
		for (size_t l = 0; l < m_msh->ele_vector.size(); l++)
		{
			// if((m_msh->ele_vector[l]->GetExcavState()>0)&&!(m_msh->ele_vector[l]->GetMark()))//WX:07.2011 HM excav
			if (ExcavMaterialGroup == static_cast<int>(m_msh->ele_vector[l]->GetPatchIndex()))
			{
				if ((m_msh->ele_vector[l]->GetExcavState() > -1) && m_msh->ele_vector[l]->GetMark())
				{
					if (m_msh->ele_vector[l]->GetExcavState() == 1)
						m_msh->ele_vector[l]->SetExcavState(0); // 1=now, 0=past
					PCS_ExcavState = 1; // not necessary
					now_Excav = true; // new elems are excavated at this time step
					// break;
				}
			}
		}
	}

	if (deact_dom.size() > 0 || now_Excav) // WX:01.2011 modified for coupled excavation
	{
		//	 		  MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);  //abort();}
		// 07.04.2010 WW
		bool done;
		ElementValue_DM* eleV_DM = NULL;
		if ((fem_dm->Flow_Type == 0 || fem_dm->Flow_Type == 2 || fem_dm->Flow_Type == 1) && Neglect_H_ini == 2)
		{
			int tmp_type, idx_p1_ini, idx_p2_ini, idx_p1_1, idx_p2_1; //, idx_sw_1;//idx_sw_ini,
			CRFProcess* h_pcs = NULL;
			h_pcs = fem_dm->h_pcs;
			tmp_type = h_pcs->type;
			idx_p1_ini = h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
			idx_p1_1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
			if (tmp_type == 1212)
			{
				idx_p2_ini = h_pcs->GetNodeValueIndex("PRESSURE2_Ini");
				// idx_sw_ini = h_pcs->GetNodeValueIndex("SATURATION1_Ini");
				idx_p2_1 = h_pcs->GetNodeValueIndex("PRESSURE2") + 1;
				// idx_sw_1 = h_pcs->GetNodeValueIndex("SATURATION1")+1;
			}
			for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			{
				// if(tmp_type==1)//LiquideFlow, RichardsFlow, MultiPhaseFlow
				h_pcs->SetNodeValue(i, idx_p1_ini, h_pcs->GetNodeValue(i, idx_p1_1));
				if (tmp_type == 1212)
				{
					// h_pcs->SetNodeValue(i, idx_p1_ini, h_pcs->GetNodeValue(i, idx_p1_1));
					h_pcs->SetNodeValue(i, idx_p2_ini, h_pcs->GetNodeValue(i, idx_p2_1));
					// h_pcs->SetNodeValue(i, idx_sw_ini, h_pcs->GetNodeValue(i, idx_sw_1));
				}
			}
		}
		for (size_t l = 0; l < m_msh->ele_vector.size(); l++)
		{
			eleV_DM = ele_value_dm[l];
			//
			(*eleV_DM->Stress0) = (*eleV_DM->Stress); // if not assembly the excavated eles in next time step
			//
			elem = m_msh->ele_vector[l];
			done = false;
			for (size_t i = 0; i < deact_dom.size(); i++)
			{
				if (elem->GetPatchIndex() == static_cast<std::size_t>(deact_dom[i]))
				{
					elem->MarkingAll(false);
					*(eleV_DM->Stress) = 0.;
					*(eleV_DM->Stress0) = 0.;
					if (eleV_DM->Stress_j)
						(*eleV_DM->Stress_j) = 0.0;
					if (eleV_DM->Stress_i)
						(*eleV_DM->Stress_i) = 0.0;
					done = true;
					break;
				}
			}
			if (ExcavMaterialGroup >= 0) // WX
			{
				if (elem->GetExcavState() >= 0)
				{
					elem->MarkingAll(false);
					// if update stress0
					(*eleV_DM->Stress) = 0.;
					(*eleV_DM->Stress0) = 0.;
					if (eleV_DM->Stress_j)
						(*eleV_DM->Stress_j) = 0.0;
					if (eleV_DM->Stress_i)
						(*eleV_DM->Stress_i) = 0.0;
					//
					done = true;
				}
			}
			if (hasAnyProcessDeactivatedSubdomains) // WX:11.2012 if there is deactivated subdomain when excavated
			{
				for (int i = 0; i < NumDeactivated_SubDomains; i++)
				{
					if (elem->GetPatchIndex() == static_cast<std::size_t>(Deactivated_SubDomain[i]))
					{
						elem->MarkingAll(false);
						done = true;
						break;
					}
				}
			}
			if (done)
				continue;
			else
				elem->MarkingAll(true);
		}

		size_t mesh_node_vector_size(m_msh->nod_vector.size());
		for (size_t l = 0; l < mesh_node_vector_size; l++)
		{
			while (m_msh->nod_vector[l]->getConnectedElementIDs().size())
				m_msh->nod_vector[l]->getConnectedElementIDs().pop_back();
		}
		// for (size_t l = 0; l < mesh_node_vector_size; l++)
		size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t l = 0; l < mesh_ele_vector_size; l++) // WX: 07.2011
		{
			elem = m_msh->ele_vector[l];
			if (!elem->GetMark())
				continue;
			// for(size_t i=0; i<elem->GetNodesNumber(m_msh->getOrder()); i++)
			for (size_t i = 0; i < elem->GetNodesNumber(1); i++) // WX:10.2011 change for one way coup. M->H
			{
				done = false;
				node = elem->GetNode(i);
				for (size_t j = 0; j < node->getConnectedElementIDs().size(); j++)
				{
					if (l == node->getConnectedElementIDs()[j])
					{
						done = true;
						break;
					}
				}
				if (!done)
					node->getConnectedElementIDs().push_back(l);
			}
		}
		if (Neglect_H_ini == 1) // WX:04.2013
			CalIniTotalStress();
	}

	// WX:10.2011 strain update for excavated node, the excavated node Strain = 0
	int Idx_Strain[6];
	Idx_Strain[0] = GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = GetNodeValueIndex("STRAIN_XY");
	if (problem_dimension_dm == 3)
	{
		Idx_Strain[4] = GetNodeValueIndex("STRAIN_XZ");
		Idx_Strain[5] = GetNodeValueIndex("STRAIN_YZ");
	}
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		node = m_msh->nod_vector[i];
		if (node->getConnectedElementIDs().size() == 0)
			for (int ii = 0; ii < (2 * problem_dimension_dm); ii++)
				fem_dm->pcs->SetNodeValue(i, Idx_Strain[ii], 0);
	}
	return;
}
/**************************************************************************
FEMLib-Method:
Task: update initial stress. if (Neglect_H_ini == 2) also update pw,pg,pc ini
Programing:
07/2011 WX
**************************************************************************/
void CRFProcessDeformation::UpdateIniStateValue()
{
	for (size_t l = 0; l < m_msh->ele_vector.size(); l++)
	{
		ElementValue_DM* eleV_DM = NULL;
		eleV_DM = ele_value_dm[l];
		//
		(*eleV_DM->Stress0) = (*eleV_DM->Stress);
	}
	if (Neglect_H_ini == 1)
		CalIniTotalStress();
	if ((fem_dm->Flow_Type == 0 || fem_dm->Flow_Type == 1 || fem_dm->Flow_Type == 2) && Neglect_H_ini == 2)
	{
		int tmp_type;
		CRFProcess* h_pcs = NULL;
		h_pcs = fem_dm->h_pcs;
		tmp_type = fem_dm->Flow_Type;
		int idx_p1_ini = h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
		int idx_p1_1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		if (tmp_type == 0 || tmp_type == 1) // LiquideFlow,RichardsFlow
		{
			for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			{
				h_pcs->SetNodeValue(i, idx_p1_ini, h_pcs->GetNodeValue(i, idx_p1_1));
			}
		}
		else if (tmp_type == 2)
		{
			int idx_p2_ini = h_pcs->GetNodeValueIndex("PRESSURE2_Ini");
			// int idx_sw_ini = h_pcs->GetNodeValueIndex("SATURATION1_Ini");
			int idx_p2_1 = h_pcs->GetNodeValueIndex("PRESSURE2") + 1;
			// int idx_sw_1 = h_pcs->GetNodeValueIndex("SATURATION1")+1;
			for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			{
				h_pcs->SetNodeValue(i, idx_p1_ini, h_pcs->GetNodeValue(i, idx_p1_1));
				h_pcs->SetNodeValue(i, idx_p2_ini, h_pcs->GetNodeValue(i, idx_p2_1));
				// h_pcs->SetNodeValue(i, idx_sw_ini, h_pcs->GetNodeValue(i, idx_sw_1));
			}
		}
	}
	return;
}

/**************************************************************************
   FEMLib-Method:
   Task: Update stresses and straines at each Gauss points
   Argument:
   Programing:
   02/2005 WW
   06/2005 WW  Parallelization
**************************************************************************/
void CRFProcessDeformation::UpdateStress()
{
	long i;
	MeshLib::CElem* elem = NULL;
	/*
	   long j, irank;
	   j = 0;
	   if(dom_vector.size()>0)
	   {
	   CPARDomain* m_dom = NULL;
	   #ifdef USE_MPI
	      irank = myrank;
	   #else
	   for(int j=0;j<(int)dom_vector.size();j++)
	   {
	   irank = j;
	   #endif
	   m_dom = dom_vector[irank];
	   for (i = 0; i < (long)m_dom->elements.size(); i++)
	   {
	   elem = m_msh->ele_vector[m_dom->elements[i]];
	   if (elem->GetMark()) // Marked for use
	   {
	   elem->SetOrder(true);
	   fem_dm->SetElementNodesDomain(m_dom->element_nodes_dom[i]); //WW
	   fem_dm->ConfigElement(elem);
	   fem_dm->m_dom = m_dom;
	   fem_dm->LocalAssembly(1);
	   }
	   }
	   #ifndef USE_MPI
	   }
	   #endif
	   }
	   else
	   {
	 */
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			elem->SetOrder(true);
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
			fem_dm->m_dom = NULL;
#endif
			fem_dm->ConfigElement(elem);
			fem_dm->LocalAssembly(1);
		}
	}
	//}
}

/**************************************************************************
   ROCKFLOW - Funktion: WriteGaussPointStress()

   Aufgabe:
   Write Gauss point stresses to a file

   Programmaenderungen:
   03/2005  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::WriteGaussPointStress(const bool last_step)
{
	if (!(idata_type == write_all_binary || idata_type == read_write))
		return;

	if ((aktueller_zeitschritt % nwrite_restart) > 0 && (!last_step))
		return;

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
	const std::string StressFileName = FileName + "_" + number2str(myrank) + ".sts";
#else
	const std::string StressFileName = FileName + ".sts";
#endif

	fstream file_stress(StressFileName.data(), ios::binary | ios::out | ios::trunc);
	ElementValue_DM* eleV_DM = NULL;

	std::size_t ActiveElements = 0;
	MeshLib::CElem* elem = NULL;
	for (std::size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (ExcavMaterialGroup > -1) // WX:if excavation write all eles
			ActiveElements++;
		else
		{
			if (elem->GetMark()) // Marked for use
				ActiveElements++;
		}
	}
	file_stress.write((char*)(&ActiveElements), sizeof(ActiveElements));
	for (std::size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark() || ExcavMaterialGroup > -1) // Marked for use//WX:if excavation write all eles
		{
			eleV_DM = ele_value_dm[i];
			file_stress.write((char*)(&i), sizeof(i));
			//          *eleV_DM->Stress_i += *eleV_DM->Stress0;
			// TEST           *eleV_DM->Stress0 = 0.0;
			eleV_DM->Write_BIN(file_stress, last_step);
		}
	}
	//
	file_stress.close();
}
/**************************************************************************
   ROCKFLOW - Funktion: ReadGaussPointStress()

   Aufgabe:
   Read Gauss point stresses

   Programmaenderungen:
   03/2005  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReadGaussPointStress()
{
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
	const std::string StressFileName = FileName + "_" + number2str(myrank) + ".sts";
#else
	const std::string StressFileName = FileName + ".sts";
#endif

	fstream file_stress(StressFileName.data(), ios::binary | ios::in);
	ElementValue_DM* eleV_DM = NULL;
	//
	std::size_t ActiveElements;
	file_stress.read((char*)(&ActiveElements), sizeof(ActiveElements));
	for (std::size_t i = 0; i < ActiveElements; i++)
	{
		std::size_t index;
		file_stress.read((char*)(&index), sizeof(index));
		eleV_DM = ele_value_dm[index];
		eleV_DM->Read_BIN(file_stress);
		(*eleV_DM->Stress0) = (*eleV_DM->Stress);
		if (eleV_DM->Stress_j)
			(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
	}
	//
	file_stress.close();
}

/**************************************************************************
   ROCKFLOW - Funktion: ReadGaussPointStress()

   Aufgabe:
   Read element-wise stress data

   Programmaenderungen:
   10/2011  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReadElementStress()
{
	long i, index, ActiveElements;
	string StressFileName = FileName + ".ele_stress.asc";
	fstream file_stress(StressFileName.data());
	ElementValue_DM* eleV_DM = NULL;
	//
	file_stress >> ActiveElements;
	for (i = 0; i < ActiveElements; i++)
	{
		file_stress >> index;
		eleV_DM = ele_value_dm[index];
		eleV_DM->ReadElementStressASCI(file_stress);
		(*eleV_DM->Stress) = (*eleV_DM->Stress0);
		if (eleV_DM->Stress_j)
			(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
	}
	//
	file_stress.close();
}

/**************************************************************************
   ROCKFLOW - Funktion: ReleaseLoadingByExcavation()

   Aufgabe:
   Compute the nodal forces produced by excavated body

   Programmaenderungen:
   04/2005  WW  Erste Version
   09/2007  WW  Set as a boundary condition
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReleaseLoadingByExcavation()
{
	long i, actElements;
	int j, k, l, SizeSt, SizeSubD;
	ElementValue_DM* ele_val = NULL;

	std::vector<int> ExcavDomainIndex;
	std::vector<long> NodesOnCaveSurface;

	CSourceTerm* m_st = NULL;
	SizeSt = (int)st_vector.size();
	bool exist = false;
	double* eqs_b = NULL;

#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
	eqs_b = eqs_new->b;
#else
	eqs_b = eqs->b;
#endif

	for (k = 0; k < SizeSt; k++)
	{
		m_st = st_vector[k];
		if (m_st->getProcessPrimaryVariable() == FiniteElement::EXCAVATION)
		{
			// ---- 16.01.2009 WW
			exist = false;

			for (j = k + 1; j < SizeSt; j++)
				if (m_st->getSubDomainIndex() == st_vector[j]->getSubDomainIndex())
				{
					exist = true;
					break;
				}
			if (!exist)
				ExcavDomainIndex.push_back(m_st->getSubDomainIndex());
		}
	}
	SizeSubD = (int)ExcavDomainIndex.size();
	if (SizeSubD == 0)
		return; // 05.09.2007 WW
	exist = false; // 16.02
	// 1. De-active host domain to be exvacated
	actElements = 0;
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		elem->SetMark(false);
		for (k = 0; k < SizeSubD; k++)
			if (elem->GetPatchIndex() == static_cast<size_t>(ExcavDomainIndex[k]))
				elem->SetMark(true);
		if (elem->GetMark())
			actElements++;
	}
	if (actElements == 0)
	{
		cout << "No element specified for excavation. Please check data in .st file "
		     << "\n";
		abort();
	}
// 2. Compute the released node loading

#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 06.11.2008, 04.2012
	SetLinearSolver(eqs);
	SetZeroLinearSolver(eqs);
#endif
	for (i = 0; i < 4; i++) // In case the domain decomposition is employed
		fem_dm->NodeShift[i] = Shift[i];
	//
	PreLoad = 11;
	LoadFactor = 1.0;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			fem_dm->ConfigElement(elem);
			fem_dm->LocalAssembly(0);
			ele_val = ele_value_dm[i];
			// Clear stresses in excavated domain
			(*ele_val->Stress0) = 0.0;
			(*ele_val->Stress) = 0.0;
			if (ele_val->Stress_j)
				(*ele_val->Stress_j) = 0.0;
		}
	}

	// 3 --------------------------------------------------------
	// Store the released loads to source term buffer
	long number_of_nodes;
	CNodeValue* m_node_value = NULL;
	std::vector<long> nodes_vector(0);

	number_of_nodes = 0;
	RecordNodeVSize((long)st_node_value.size());

	// TEST
	st_node_value.clear();
	//

	for (k = 0; k < SizeSt; k++)
	{
		// Get nodes on cave surface
		m_st = st_vector[k];
		if (m_st->getProcessPrimaryVariable() != FiniteElement::EXCAVATION)
			continue;
		if (m_st->getGeoType() == GEOLIB::POLYLINE)
		{
			CGLPolyline* m_polyline(GEOGetPLYByName(m_st->getGeoName()));

			// reset the min edge length of mesh
			double mesh_min_edge_length(m_msh->getMinEdgeLength());
			m_msh->setMinEdgeLength(m_polyline->epsilon);

			if (m_st->getGeoObj())
			{
				m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(m_st->getGeoObj()), nodes_vector);
				// reset min edge length of mesh
				m_msh->setMinEdgeLength(mesh_min_edge_length);
			}
			m_msh->setMinEdgeLength(mesh_min_edge_length);
		}
		if (m_st->getGeoType() == GEOLIB::SURFACE)
		{
			// CC 10/05
			Surface* m_surface = GEOGetSFCByName(m_st->getGeoName());
			//			 07/2010 TF ToDo: to do away with the global vector surface_vector
			//			                  fetch the geometry from CFEMesh
			//			Surface *m_surface (surface_vector[m_st->getGeoObjIdx()]);
			if (m_surface)
			{
				if (m_surface->type == 100)
					m_msh->GetNodesOnCylindricalSurface(m_surface, nodes_vector);
				else
					m_msh->GetNODOnSFC_PLY(m_surface, nodes_vector);
			}
		}
		// Set released node forces from eqs->b;
		number_of_nodes = (int)nodes_vector.size();
		for (j = 0; j < problem_dimension_dm; j++)
			for (i = 0; i < number_of_nodes; i++)
			{
				m_node_value = new CNodeValue();
				m_node_value->msh_node_number = nodes_vector[i] + Shift[j];
				m_node_value->geo_node_number = nodes_vector[i];
				m_node_value->node_value = -eqs_b[m_node_value->geo_node_number + Shift[j]];
				m_node_value->CurveIndex = m_st->CurveIndex;
				// Each node only take once
				exist = false;
				for (l = 0; l < (int)st_node_value.size(); l++)
					if (st_node_value[l]->msh_node_number == m_node_value->msh_node_number)
					{
						exist = true;
						break;
					}
				if (!exist)
					st_node_value.push_back(m_node_value);
			}
	}
	//
	// Deactivate the subdomains to be excavated
	if (Deactivated_SubDomain)
		delete[] Deactivated_SubDomain;
	Deactivated_SubDomain = new int[SizeSubD];
	NumDeactivated_SubDomains = SizeSubD;
	for (j = 0; j < SizeSubD; j++)
		Deactivated_SubDomain[j] = ExcavDomainIndex[j];

	// Activate the host domain for excavtion analysis
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark())
			elem->SetMark(true);
	}
	PreLoad = 1;
	// TEST OUTPUT
	//   {MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);  abort();}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::UpdateInitialStress()
   Task:  Compute number of element neighbors to a node
   Dim : Default=2
   Programming:
   12/2003 WW
 **************************************************************************/
void CRFProcessDeformation::UpdateInitialStress(bool ZeroInitialS)
{
	long i;
	ElementValue_DM* eval_DM;

	// Over all elements
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			eval_DM = ele_value_dm[i];
			if (ZeroInitialS)
				(*eval_DM->Stress0) = 0.0;
			else
				(*eval_DM->Stress0) = (*eval_DM->Stress);
		}
	}
}
/**************************************************************************
   GEOSYS - Funktion: CalcBC_or_SecondaryVariable_Dynamics(bool BC);
   Programmaenderungen:
   05/2005  WW  Erste Version
   letzte Aenderung:
   09/2011 TF substituted pow by fastpow
**************************************************************************/
bool CRFProcessDeformation::CalcBC_or_SecondaryVariable_Dynamics(bool BC)
{
	const char* function_name[7];
	size_t i;
	long j;
	double v, bc_value, time_fac = 1.0;

	std::vector<int> bc_type;
	long bc_msh_node;
	long bc_eqs_index;
	int interp_method = 0;
	int curve, valid = 0;
	int idx_disp[3], idx_vel[3], idx_acc[3], idx_acc0[3];
	int idx_pre, idx_dpre, idx_dpre0;
	int nv, k;

	size_t Size = m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false);
	CBoundaryCondition* m_bc = NULL;
	bc_type.resize(Size);

	v = 0.0;
	// 0: not given
	// 1, 2, 3: x,y, or z is given
	for (size_t i = 0; i < Size; i++)
		bc_type[i] = 0;

	idx_dpre0 = GetNodeValueIndex("PRESSURE_RATE1");
	idx_dpre = idx_dpre0 + 1;
	idx_pre = GetNodeValueIndex("PRESSURE1");

	nv = 2 * problem_dimension_dm + 1;
	if (m_msh->GetCoordinateFlag() / 10 == 2) // 2D
	{
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
		idx_disp[0] = GetNodeValueIndex("DISPLACEMENT_X1");
		idx_disp[1] = GetNodeValueIndex("DISPLACEMENT_Y1");
		idx_vel[0] = GetNodeValueIndex("VELOCITY_DM_X");
		idx_vel[1] = GetNodeValueIndex("VELOCITY_DM_Y");
		idx_acc0[0] = GetNodeValueIndex("ACCELERATION_X1");
		idx_acc0[1] = GetNodeValueIndex("ACCELERATION_Y1");
		idx_acc[0] = idx_acc0[0] + 1;
		idx_acc[1] = idx_acc0[1] + 1;
	}
	else if (m_msh->GetCoordinateFlag() / 10 == 3) // 3D
	{
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
		idx_disp[0] = GetNodeValueIndex("DISPLACEMENT_X1");
		idx_disp[1] = GetNodeValueIndex("DISPLACEMENT_Y1");
		idx_disp[2] = GetNodeValueIndex("DISPLACEMENT_Z1");
		idx_vel[0] = GetNodeValueIndex("VELOCITY_DM_X");
		idx_vel[1] = GetNodeValueIndex("VELOCITY_DM_Y");
		idx_vel[2] = GetNodeValueIndex("VELOCITY_DM_Z");
		idx_acc0[0] = GetNodeValueIndex("ACCELERATION_X1");
		idx_acc0[1] = GetNodeValueIndex("ACCELERATION_Y1");
		idx_acc0[2] = GetNodeValueIndex("ACCELERATION_Z1");
		for (k = 0; k < 3; k++)
			idx_acc[k] = idx_acc0[k] + 1;
	}

	//
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* m_bc_node = bc_node_value[i];
		m_bc = bc_node[i];
		for (j = 0; j < nv; j++)
			if (convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).compare(function_name[j]) == 0)
				break;
		if (j == nv)
		{
			cout << "No such primary variable found in CalcBC_or_SecondaryVariable_Dynamics."
			     << "\n";
			abort();
		}
		bc_msh_node = m_bc_node->geo_node_number;
		if (!m_msh->nod_vector[bc_msh_node]->GetMark())
			continue;
		if (bc_msh_node >= 0)
		{
			curve = m_bc_node->CurveIndex;
			if (curve > 0)
			{
				time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
				if (!valid)
					continue;
			}
			else
				time_fac = 1.0;
			bc_value = time_fac * m_bc_node->node_value;
			bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();

			if (BC)
			{
				if (j < problem_dimension_dm) // da
				{
					bc_eqs_index += Shift[j];
// da = v = 0.0;
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS) // WW
					eqs_new->SetKnownX_i(bc_eqs_index, 0.);
#else
					MXRandbed(bc_eqs_index, 0.0, eqs->b);
#endif
				}
				else if (j == nv - 1) // P
				{
					bc_eqs_index += Shift[problem_dimension_dm];
// da = v = 0.0;
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS) // WW
					eqs_new->SetKnownX_i(bc_eqs_index, 0.);
#else
					MXRandbed(bc_eqs_index, 0.0, eqs->b);
#endif
				}
			}
			else
			{
				// Bit operator
				if (!(bc_type[bc_eqs_index] & (int)MathLib::fastpow(2, j)))
					bc_type[bc_eqs_index] += (int)MathLib::fastpow(2, j);
				if (j < problem_dimension_dm) // Disp
				{
					SetNodeValue(bc_eqs_index, idx_disp[j], bc_value);
					SetNodeValue(bc_eqs_index, idx_vel[j], 0.0);
					SetNodeValue(bc_eqs_index, idx_acc[j], 0.0);
				}
				// Vel
				else if (j >= problem_dimension_dm && j < nv - 1)
				{
					v = GetNodeValue(bc_eqs_index, idx_disp[j]);
					v += bc_value * dt
					     + 0.5 * dt * dt
					           * (ARRAY[bc_eqs_index + Shift[j]]
					              + m_num->GetDynamicDamping_beta2() * GetNodeValue(bc_eqs_index, idx_acc0[j]));
					SetNodeValue(bc_eqs_index, idx_disp[j], v);
					SetNodeValue(bc_eqs_index, idx_vel[j], bc_value);
				}
				else if (j == nv - 1) // Vel
				{ // p
					SetNodeValue(bc_eqs_index, idx_pre, bc_value);
					SetNodeValue(bc_eqs_index, idx_dpre, 0.0);
				}
			}
		}
	}
	if (BC)
		return BC;

	// BC
	for (i = 0; i < m_msh->GetNodesNumber(true); i++)
		for (k = 0; k < problem_dimension_dm; k++)
		{
			// If boundary
			if (bc_type[i] & (int)MathLib::fastpow(2, k))
				continue; // u
			//
			v = GetNodeValue(i, idx_disp[k]);
			v += GetNodeValue(i, idx_vel[k]) * dt
			     + 0.5 * dt * dt
			           * (ARRAY[i + Shift[k]] + m_num->GetDynamicDamping_beta2() * GetNodeValue(i, idx_acc0[k]));
			SetNodeValue(i, idx_disp[k], v);
			if (bc_type[i] & (int)MathLib::fastpow(2, k + problem_dimension_dm))
				continue;
			// v
			v = GetNodeValue(i, idx_vel[k]);
			v += dt * ARRAY[i + Shift[k]] + m_num->GetDynamicDamping_beta1() * dt * GetNodeValue(i, idx_acc0[k]);
			SetNodeValue(i, idx_vel[k], v);
		}

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		if (bc_type[i] & (int)MathLib::fastpow(2, (nv - 1)))
			continue;
		v = GetNodeValue(i, idx_pre);
		v += ARRAY[i + Shift[problem_dimension_dm]] * dt
		     + m_num->GetDynamicDamping_beta1() * dt * GetNodeValue(i, idx_dpre0);
		SetNodeValue(i, idx_pre, v);
	}

	return BC;
}

bool CRFProcessDeformation::isDynamic() const
{
	return fem_dm->dynamic;
}

} // end namespace
