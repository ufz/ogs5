/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_std.h"

#include "rf_mmp_new.h"
#include "rf_msp_new.h"

#include "rfmat_cp.h"

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
#include "PETSC/PETScLinearSolver.h"
#else
#ifndef NEW_EQS // WW. 06.11.2008
// Sytem matrix
#include "matrix_routines.h"
#endif
#endif

#ifdef NEW_EQS
#include "equation_class.h"
using Math_Group::CSparseMatrix;
#endif

#include "PhysicalConstant.h"

namespace
{
static inline double ipol(double const* const a, double const* const b, const double theta,
                          FiniteElement::CElement const* const obj)
{
	return (1.0 - theta) * obj->interpolate(a) + theta * obj->interpolate(b);
}
}

namespace FiniteElement
{
/*************************************************************************
    Programming:
    07/2013 TN
**************************************************************************/
void CFiniteElementStd::CalcMassTES()
{
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		// Get local coordinates and weights
		// Compute Jacobian matrix and its determinate
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		getShapefunctValues(gp, 1);       // Linear interpolation function

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				// Material
				const double coeff = CalCoefMassTES(in * pcs->dof + jn);

				const double mat_fac = fkt * coeff;
// Calculate mass matrix

#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
				const int jn_offset = jn * nnodes;
				for (int i = 0; i < act_nodes; i++)
				{
					const int ia = local_idx[i];
					const int ib = i + in * nnodes;
					for (int j = 0; j < nnodes; j++)
					{
						(*Mass2)(ib, j + jn_offset) += mat_fac * shapefct[ia] * shapefct[j];
					}
				}
#else
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
						(*Mass2)(i + in * nnodes, j + jn * nnodes) += mat_fac * shapefct[i] * shapefct[j];
				}
#endif
			}
		}
	}

	// std::cout << __FUNCTION__ << ":" << __LINE__ << ":\n"
	//           << (*Mass2) << std::endl;
}

/*************************************************************************
    Programming:
    07/2013 TN
**************************************************************************/
void CFiniteElementStd::CalcLumpedMassTES()
{
	int gp_r, gp_s, gp_t;
	const int nDF = pcs->dof;
	double vol = 0.0;

	// Volume
	if (axisymmetry)
	{ // This calculation should be done in CompleteMesh.
		// However, in order not to destroy the concise of the code,
		// it is put here. Anyway it is computational cheap. WW
		vol = 0.0;
		for (gp = 0; gp < nGaussPoints; gp++)
			//  Get local coordinates and weights
			//  Compute Jacobian matrix and its determinate
			vol += GetGaussData(gp, gp_r, gp_s, gp_t);
	}
	else
		vol = MeshElement->GetVolume(); //* MeshElement->area;

	// Initialize
	(*Mass2) = 0.0;
	// Center of the reference element
	getShapeFunctionCentroid();               // Linear interpolation function

	for(int in = 0; in < nDF; in++)
	{
		const int ish = in * nnodes;
		for (int jn = 0; jn < nDF; jn++)
		{
			const int jsh = jn * nnodes;
			double factor = CalCoefMassTES(in * nDF + jn);

			//			pcs->timebuffer = factor; // Tim Control "Neumann"
			factor *= vol;
			for (int i = 0; i < nnodes; i++)
			{
				(*Mass2)(i + ish, i + jsh) = shapefct[i] * factor;
			}
		}
	}
}

/**************************************************************************
    FEMLib-Method:
    Task: Calculate material coefficient for mass matrix
    Implementaion:
    03/2011 AKS /  NB
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoefMassTES(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	double const* const X0 = NodalVal_X0;
	double const* const X1 = NodalVal_X1;

	double& p = eos_arg[0];
	double& T = eos_arg[1];
	double& X = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	const int Index = MeshElement->GetIndex();
	const ElementValue* gp_ele = ele_gp_value[Index];
	poro = MediaProp->Porosity(Index, theta);

	double val = 0.0;

	switch (dof_index)
	{
		case 0: // M_pp
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);
			val = poro / p * FluidProp->Density(eos_arg);
			break;

		case 1: // M_pT
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);
			val = -poro / T * FluidProp->Density(eos_arg);
			break;

		case 2: // M_px
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			const double M0 = cp_vec[0]->molar_mass; // inert
			const double M1 = cp_vec[1]->molar_mass; // reactive

			double dxn_dxm = M0 * M1; // 0 is inert, 1 is reactive
			dxn_dxm /= (M0 * X + M1 * (1.0 - X)) * (M0 * X + M1 * (1.0 - X));

			val = (M1 - M0) * p / (PhysicalConstant::IdealGasConstant * T) * dxn_dxm * poro;
			break;
		}

		case 3: // M_Tp
			val = -poro;
			break;

		case 4: // M_TT
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			const double rhoSR = gp_ele->rho_s_curr[gp];
			const double rhoGR = FluidProp->Density(eos_arg);
			const double cpG = FluidProp->SpecificHeatCapacity(eos_arg);
			const double cpS = SolidProp->Heat_Capacity(rhoSR);

			val = poro * rhoGR * cpG + (1.0 - poro) * rhoSR * cpS;
			break;
		}

		//    case 5:
		//    case 6:
		//    case 7:
		//        break;

		case 8: // M_xx
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);
			val = poro * FluidProp->Density(eos_arg);
			break;
	}

	return val;
}

/**************************************************************************
    FEMLib-Method:
    Task: Calculate material coefficient for Laplacian matrix
    Implementaion:
    03/2011 AKS /  NB
    07/2013 TN
    04/2015 CL
**************************************************************************/
void CFiniteElementStd::CalCoefLaplaceTES(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	double const* const X0 = NodalVal_X0;
	double const* const X1 = NodalVal_X1;

	double& p = eos_arg[0];
	double& T = eos_arg[1];
	double& X = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	const int Index = MeshElement->GetIndex();

	switch (dof_index)
	{
		case 0:
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			double* tensor = MediaProp->PermeabilityTensor(Index);
			double k_rel = 1.0;
			if (MediaProp->flowlinearity_model > 0)
			{
				k_rel = MediaProp->NonlinearFlowFunction(Index, gp, theta, this);
			}

			double val = FluidProp->Density(eos_arg) * k_rel / FluidProp->Viscosity(eos_arg);

			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = tensor[i] * val;
			}
			break;
		}

		//    case 1:
		//    case 2:
		//    case 3:
		//        break;

		case 4:
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			// TODO [CL]: only diagonal neeeded, and only one array needed
			double fluid_heat_conductivity_tensor[9] = {0.};
			double solid_heat_conductivity_tensor[9] = {0.};

			poro = MediaProp->Porosity(Index, theta);
			const double lamf = FluidProp->HeatConductivity(eos_arg);
			double lams = SolidProp->Heat_Conductivity();

			if (SolidProp->getSolidReactiveSystem() == FiniteElement::Z13XBF)
			{
				ElementValue* gp_ele = ele_gp_value[Index];
				double C = gp_ele->rho_s_curr[gp] / SolidProp->lower_solid_density_limit - 1.;
				const double lambda_ads = 0.7; // TODO [CL] Find relation for this
				lams += C * SolidProp->lower_solid_density_limit / pcs->m_conversion_rate->get_adsorbate_density(T)
				        * (lambda_ads - lamf);
			}

			for (size_t i = 0; i < dim; i++)
			{
				fluid_heat_conductivity_tensor[i * dim + i] = poro * lamf;
				solid_heat_conductivity_tensor[i * dim + i] = (1.0 - poro) * lams;
			}

			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = fluid_heat_conductivity_tensor[i] + solid_heat_conductivity_tensor[i];
			}
			break;
		}

		//    case 5:
		//    case 6:
		//    case 7:
		//        break;

		case 8:
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			double diffusion_tensor[9] = {0.};

			poro = MediaProp->Porosity(Index, theta);
			tort = MediaProp->TortuosityFunction(Index, unit, theta);
			const double diffusion_coefficient_component
			    = cp_vec[1]->CalcDiffusionCoefficientCP(Index, theta, pcs); // coefficient of reactive (2nd) component
			const double rhoGR = FluidProp->Density(eos_arg);

			for (size_t i = 0; i < dim; i++)
			{
				// TODO [CL] poro?
				diffusion_tensor[i * dim + i] = tort * poro * rhoGR * diffusion_coefficient_component;
			}

			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = diffusion_tensor[i]; // TN
			}
			break;
		}

		default:
			std::fill_n(mat, dim * dim, 0);
	}
}

/***************************************************************************
    GeoSys - Funktion:
         CFiniteElementStd:: CalcAdvection
    Aufgabe:  Calculate the advection matrix

    Programming:
    01/2005   WW
    02/2005   OK GEO factor
    09/2005   SB - adapted to advection
    03/2007   WW - Fluid advection with multiphase flow
    05/2008   WW - General densty for multiphase flow
    01/2010   NW - SUPG
    07/2013 adapted for TES
**************************************************************************/
void CFiniteElementStd::CalcAdvectionTES()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;
	ElementValue* gp_ele = ele_gp_value[Index];

	for (gp = 0; gp < nGaussPoints; gp++)
	{
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		getShapefunctValues(gp, 1);
		getGradShapefunctValues(gp, 1);

		// Velocity
		// TODO [CL] vel includes porosity? cf. \tilde w
		double vel[] = {gp_ele->Velocity(0, gp), gp_ele->Velocity(1, gp), gp_ele->Velocity(2, gp)};

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				const double coeff = CalCoefAdvectionTES(in * pcs->dof + jn);
				const double mat_fac = fkt * coeff;

#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
				const int jn_offset = jn * nnodes;
				for (int i = 0; i < act_nodes; i++)
				{
					const int ia = local_idx[i];
					const int ib = i + in * nnodes;
					for (int j = 0; j < nnodes; j++)
					{
						for (size_t k = 0; k < dim; k++)
						{
							(*Advection)(ib, j + jn_offset)
							    += mat_fac * shapefct[ia] * vel[k] * dshapefct[k * nnodes + j];
						}
					}
				}
#else
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
					{
						for (size_t k = 0; k < dim; k++)
						{
							(*Advection)(i + in * nnodes, j + jn * nnodes)
							    += mat_fac * shapefct[i] * vel[k] * dshapefct[k * nnodes + j];
						}
					}
				}
#endif
			}
		}
	}

	// std::cout << __FUNCTION__ << ":" << __LINE__ << ":\n"
	//           << (*Advection) << std::endl;
}

/**************************************************************************
   FEMLib-Method:
    Task: Calculate material coefficient for advection matrix
    Programing:
    01/2005 WW/OK Implementation
    03/2005 WW Heat transport
    07/2005 WW Change for geometry element object
    09/2005 SB
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoefAdvectionTES(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	double const* const X0 = NodalVal_X0;
	double const* const X1 = NodalVal_X1;

	double& p = eos_arg[0];
	double& T = eos_arg[1];
	double& X = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	double val = 0.0;

	switch (dof_index)
	{
		//    case 0:
		//    case 1:
		//    case 2:
		//    case 3:
		//        break;

		case 4:
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			val = FluidProp->Density(eos_arg) * FluidProp->SpecificHeatCapacity(eos_arg);
			break;

		//    case 5:
		//    case 6:
		//    case 7:
		//        break;

		case 8:
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			val = FluidProp->Density(eos_arg);
			break;
	}

	return val;
}

/***************************************************************************
    FEMLib-Method:
    Task: Assembly of ContentMatrix for
    TES
    Implementaion:
    03/2011 AKS
    07/2013 TN
**************************************************************************/
void CFiniteElementStd::CalcContentTES()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;

	for (gp = 0; gp < nGaussPoints; gp++)
	{
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		getShapefunctValues(gp, 1);

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				const double coeff = CalCoefContentTES(in * pcs->dof + jn);
				double mat_fac = fkt * coeff;

#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
				const int jn_offset = jn * nnodes;
				for (int i = 0; i < act_nodes; i++)
				{
					const int ia = local_idx[i];
					const int ib = i + in * nnodes;
					for (int j = 0; j < nnodes; j++)
					{
						(*Content)(ib, j + jn_offset) += mat_fac * shapefct[ia] * shapefct[j];
					}
				}
#else
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
					{
						(*Content)(i + in * nnodes, j + jn * nnodes) += mat_fac * shapefct[i] * shapefct[j];
					}
				}
#endif
			}
		}
	}

	// std::cout << __FUNCTION__ << ":" << __LINE__ << ":\n"
	//           << (*Content) << std::endl;
}

/**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoefContentTES(const int dof_index)
{
	const int Index = MeshElement->GetIndex();
	const ElementValue* gp_ele = ele_gp_value[Index];

	const double theta = pcs->m_num->ls_theta;

	double val = 0.0;

	switch (dof_index)
	{
		// gas flow
		//    case 0:
		//    case 1:
		//    case 2:
		//    case 3:
		//    case 4:
		//    case 5:
		//    case 6:
		//    case 7:
		//        break;
		case 8: // x x
			val = (MediaProp->Porosity(Index, theta) - 1.0) * gp_ele->q_R[gp];
			break;
	}

	return val;
}

/***************************************************************************
    GeoSys - Funktion:
    Assemble_RHS_TES:
    11/2011   AKS
    07/2013 TN
**************************************************************************/
void CFiniteElementStd::Assemble_RHS_TES()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;

	for (int i = 0; i < pcs->dof * nnodes; i++)
		NodalVal[i] = 0.0;

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		// Compute geometry
		getShapefunctValues(gp, 1);

		for (int ii = 0; ii < pcs->dof; ii++)
		{
			const double fac = CalCoef_RHS_TES(ii);

			for (int i = 0; i < nnodes; i++)
				NodalVal[i + ii * nnodes] += fac * fkt * shapefct[i];
		}
	}

	for (int ii = 0; ii < pcs->dof; ii++)
	{
		const int ii_sh = ii * nnodes;
		// std::cout << ii << " " << i_sh << " " << ii_sh << "\n";
		for (int i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC) // && !defined(other parallel libs)//07~07.2014. TN
			const long i_sh = NodeShift[ii];
			eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i + ii_sh];
#else
			(*RHS)[i + LocalShift + ii_sh] += NodalVal[i + ii_sh];
#endif
		}
	}
}

/**************************************************************************
    FEMLib-Method:
    Task: Calculate  coefficient of temperature induced RHS of multi-phase flow
    Programing:
    02/2007 WW Implementation
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_TES(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	double const* const X0 = NodalVal_X0;
	double const* const X1 = NodalVal_X1;

	double& pg = eos_arg[0];
	double& Tg = eos_arg[1];
	double& Xw = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	const int Index = MeshElement->GetIndex();
	poro = MediaProp->Porosity(Index, theta);
	const ElementValue* gp_ele = ele_gp_value[Index];
	const double q_r = gp_ele->q_R[gp]; // reaction rate

	double val = 0.0;

	// NodalVal0 and NodalVal1 is the pressure on previous and current time step.
	//
	switch (dof_index)
	{
		case 0:
			val = (poro - 1.0) * q_r;
			break;

		case 1:
		{
			pg = ipol(p0, p1, theta, this);
			Tg = ipol(T0, T1, theta, this);
			Xw = ipol(X0, X1, theta, this);

			val = FluidProp->Density(eos_arg) * poro * FluidProp->specific_heat_source;

			double H_vap(0.);
			if (SolidProp->getSolidReactiveSystem() == FiniteElement::Z13XBF)
			{
				const double mole_frac = pcs->m_conversion_rate->get_mole_fraction(Xw);
				H_vap = pcs->m_conversion_rate->get_enthalpy(Tg, pg * mole_frac);
			} else if (SolidProp->getSolidReactiveSystem() != FiniteElement::INERT){
				// sign convention:
				// defined negative for exothermic composition reaction but equ. written as:
				// AB + \Delta H <--> A + B
				H_vap = - SolidProp->reaction_enthalpy;
				//enthalpy correction
				const double rhoSR = gp_ele->rho_s_curr[gp];
				const double dcp_drhoSR((((*SolidProp->data_Capacity)(1)*SolidProp->upper_solid_density_limit -
								(*SolidProp->data_Capacity)(0)*SolidProp->lower_solid_density_limit)/
							(SolidProp->upper_solid_density_limit-SolidProp->lower_solid_density_limit) -
							(*SolidProp->data_Capacity)(0)) * SolidProp->lower_solid_density_limit/(rhoSR*rhoSR));

				const double cpS = SolidProp->Heat_Capacity(rhoSR);
				const double cpG = FluidProp->SpecificHeatCapacity(eos_arg);
				H_vap -= (cpS - cpG + rhoSR * dcp_drhoSR)*(Tg - 573.15);//TODO: Move IC to input file

			}
			val += (1.0-poro) * q_r * H_vap;
			val += gp_ele->rho_s_curr[gp] * (1.0-poro) * SolidProp->specific_heat_source;
		}
		break;

		case 2:
			val = (poro - 1.0) * q_r;
			if (Xw < 0.0)
				val += 100.;
			break;
	}

	return val;
}
}
