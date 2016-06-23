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
void CFiniteElementStd::CalcMassTNEQ()
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
		ComputeShapefct(1); // Linear interpolation function

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				// Material
				double mat_fac = fkt * CalCoefMassTNEQ(in * pcs->dof + jn);
				// Calculate mass matrix
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
						(*Mass2)(i + in * nnodes, j + jn * nnodes) += mat_fac * shapefct[i] * shapefct[j];
				}
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
double CFiniteElementStd::CalCoefMassTNEQ(const int dof_index)
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
			if (GasMassForm)
			{
				p = ipol(p0, p1, theta, this);
				T = ipol(T0, T1, theta, this);
				X = ipol(X0, X1, theta, this);
				val = poro / p * FluidProp->Density(eos_arg);
			}
			else
			{
				p = ipol(p0, p1, theta, this);
				val = poro / p;
			}
			break;

		case 1: // M_pT
			if (GasMassForm)
			{
				p = ipol(p0, p1, theta, this);
				T = ipol(T0, T1, theta, this);
				X = ipol(X0, X1, theta, this);
				val = -poro / T * FluidProp->Density(eos_arg);
			}
			else
			{
				T = ipol(T0, T1, theta, this);
				val = -poro / T;
			}
			break;

		//    case 2:
		//        val = 0.0;
		//        break;

		case 3: // M_px
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			const double M0 = cp_vec[0]->molar_mass; // inert
			const double M1 = cp_vec[1]->molar_mass; // reactive

			double dxn_dxm = M0 * M1; // 0 is inert, 1 is reactive
			dxn_dxm /= (M0 * X + M1 * (1.0 - X)) * (M0 * X + M1 * (1.0 - X));

			if (GasMassForm)
			{
				val = (M1 - M0) * p / (PhysicalConstant::IdealGasConstant * T) * dxn_dxm * poro;
			}
			else
			{
				val = (M1 - M0) * p / (PhysicalConstant::IdealGasConstant * T) * dxn_dxm * poro
				      / FluidProp->Density(eos_arg);
			}
			break;
		}

		case 4: // M_Tp
			if (FluidProp->beta_T == 0.0)
			{
				// TN: beta_T read as 0 from input file. This leads to neglection of this term.
				val = -poro;
			}
			else
			{
				T = ipol(T0, T1, theta, this);
				val = -poro * FluidProp->beta_T * T;
			}
			break;

		case 5: // M_TT
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			val = poro * FluidProp->Density(eos_arg) * FluidProp->SpecificHeatCapacity(eos_arg);

			break;

		//    case 6:
		//    case 7:
		//    case 8:
		//    case 9:
		//        val = 0.0;
		//        break;

		case 10: // M_TT^S
		{
			const double rho_s = gp_ele->rho_s_curr[gp];
			val = (1 - poro) * rho_s * SolidProp->Heat_Capacity(rho_s); // SolidProp->Density()
			break;
		}

		//    case 11:
		//    case 12:
		//    case 13:
		//    case 14:
		//        val = 0.0;
		//        break;

		case 15: // M_xx
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
  Pressure-Temperature Coupled global approach
  Implementaion:
    03/2011 AKS /  NB
    07/2013 TN
    04/2015 CL
**************************************************************************/
void CFiniteElementStd::CalCoefLaplaceTNEQ(const int dof_index)
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

	double fluid_heat_conductivity_tensor[9];
	double solid_heat_conductivity_tensor[9];
	double diffusion_tensor[9];

	for (size_t i = 0; i < dim * dim; i++)
	{
		fluid_heat_conductivity_tensor[i] = 0.0;
		solid_heat_conductivity_tensor[i] = 0.0;
		diffusion_tensor[i] = 0.0;
		mat[i] = 0.0;
	}

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

			double val;
			if (GasMassForm)
			{
				val = FluidProp->Density(eos_arg) * k_rel / FluidProp->Viscosity(eos_arg);
			}
			else
			{
				val = k_rel / FluidProp->Viscosity(eos_arg);
			}

			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = tensor[i] * val;
			}
			break;
		}

		//    case 1:
		//    case 2:
		//    case 3:
		//    case 4:
		//        break;

		case 5:
		{
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			poro = MediaProp->Porosity(Index, theta);
			const double lamf = FluidProp->HeatConductivity(eos_arg);

			for (size_t i = 0; i < dim; i++)
			{
				fluid_heat_conductivity_tensor[i * dim + i] = poro * lamf;
			}
			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = fluid_heat_conductivity_tensor[i];
			}
			break;
		}

		//    case 6:
		//    case 7:
		//    case 8:
		//    case 9:
		//        break;

		case 10:
		{
			poro = MediaProp->Porosity(Index, theta);
			const double lams = SolidProp->Heat_Conductivity();

			for (size_t i = 0; i < dim; i++)
			{
				solid_heat_conductivity_tensor[i * dim + i] = (1 - poro) * lams;
			}
			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = solid_heat_conductivity_tensor[i];
			}
			break;
		}

		//    case 11:
		//    case 12:
		//    case 13:
		//    case 14:
		//        break;

		case 15:
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			poro = MediaProp->Porosity(Index, theta);
			tort = MediaProp->TortuosityFunction(Index, unit, theta);
			const double diffusion_coefficient_component
			    = cp_vec[1]->CalcDiffusionCoefficientCP(Index, theta, pcs); // coefficient of reactive (2nd) component
			const double rhoGR = FluidProp->Density(eos_arg);

			for (size_t i = 0; i < dim; i++)
			{
				diffusion_tensor[i * dim + i] = tort * poro * rhoGR * diffusion_coefficient_component;
			}

			for (size_t i = 0; i < dim * dim; i++)
			{
				mat[i] = diffusion_tensor[i]; // TN
			}
			break;
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
    07/2013 adapted for TNEQ
**************************************************************************/
void CFiniteElementStd::CalcAdvectionTNEQ()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;
	ElementValue* gp_ele = ele_gp_value[Index];
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		ComputeGradShapefct(1);
		ComputeShapefct(1);

		// Velocity
		double vel[] = {gp_ele->Velocity(0, gp), gp_ele->Velocity(1, gp), gp_ele->Velocity(2, gp)};

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				double mat_fac = fkt * CalCoefAdvectionTNEQ(in * pcs->dof + jn);
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
			}
		}
	}
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
double CFiniteElementStd::CalCoefAdvectionTNEQ(const int dof_index)
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

	if (!GasMassForm)
	{
		switch (dof_index)
		{
			case 0:
				val = 1.0 / p;
				break;

			case 1:
				T = ipol(T0, T1, theta, this);
				// T = (1-pcs->m_num->ls_theta)*interpolate(NodalVal_t0) +
				// pcs->m_num->ls_theta*interpolate(NodalVal_t1); //NW include theta
				val = -1.0 / T;

				break;

			// case 2:

			case 3:
				poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
				val = 0.0;
				break;
		}
	}

	switch (dof_index)
	{
		case 4:
			if (FluidProp->beta_T == 0.0)
			{
				// TODO [CL] check logic
				if (MediaProp->getFrictionPhase() == FiniteElement::FLUID)
				{
					val = 0.0;
				}
				else
				{
					val = -1.0;
				}
			}
			else
			{
				T = ipol(T0, T1, theta, this);
				val = -FluidProp->beta_T * T;

				if (MediaProp->getFrictionPhase() == FiniteElement::FLUID)
				{
					val += 1.0;
				}
			}
			break;

		case 5:
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = ipol(X0, X1, theta, this);

			val = FluidProp->Density(eos_arg) * FluidProp->SpecificHeatCapacity(eos_arg);
			break;

		// case 6:
		// case 7:
		// case 8:
		// case 9:
		// case 10:
		// case 11:
		// case 12:
		// case 13:
		// case 14:
		//    val = 0.0;
		//    break;

		case 15:
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
   TNEQ
   Implementaion:
   03/2011 AKS
**************************************************************************/
void CFiniteElementStd::CalcContentTNEQ()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		ComputeShapefct(1);

		for (int in = 0; in < pcs->dof; in++)
		{
			for (int jn = 0; jn < pcs->dof; jn++)
			{
				double mat_fac = fkt * CalCoefContentTNEQ(in * pcs->dof + jn);
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
					{
						(*Content)(i + in * nnodes, j + jn * nnodes) += mat_fac * shapefct[i] * shapefct[j];
					}
				}
			}
		}
	}
}

/**************************************************************************
  FEMLib-Method:
  Programing:
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoefContentTNEQ(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	// double const * const X0 = NodalVal_X0;
	// double const * const X1 = NodalVal_X1;

	double& p = eos_arg[0];
	double& T = eos_arg[1];
	double& X = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	int Index = MeshElement->GetIndex();
	ElementValue* gp_ele = ele_gp_value[Index];

	double val = 0.0;

	switch (dof_index)
	{
		// gas flow
		// case 0:
		// case 1:
		// case 2:
		// case 3:
		//     val = 0.0;
		//     break;

		// heat in gas
		// case 4:
		//     val = 0.0;//-phi_g^-1 * grad phi_g
		//     break;

		case 5: // T T
			val = MediaProp->HeatTransferCoefficient(Index, theta, this);
			break;
		case 6: // T T^S
			val = -MediaProp->HeatTransferCoefficient(Index, theta, this);
			break;

		// case 7: // T x
		// case 8: // T^S p
		//     break;

		case 9: // T^ST
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = 1.0; // TN - only reactive component for specific heat capacity here!
			val = -MediaProp->HeatTransferCoefficient(Index, theta, this);
			val -= (1.0 - MediaProp->Porosity(Index, theta)) * gp_ele->q_R[gp]
			       * FluidProp->SpecificHeatCapacity(eos_arg); // TN
			break;

		case 10: // T^S T^S
			val = MediaProp->HeatTransferCoefficient(Index, theta, this);
			p = ipol(p0, p1, theta, this);
			T = ipol(T0, T1, theta, this);
			X = 1.0; // TN - only reactive component for specific heat capacity here!
			val += (1.0 - MediaProp->Porosity(Index, theta)) * gp_ele->q_R[gp]
			       * FluidProp->SpecificHeatCapacity(eos_arg); // TN
			break;

		//    case 11:
		//    case 12:
		//    case 13:
		//    case 14:
		//        val = 0.0;
		//        break;

		case 15: // x x
			val = (MediaProp->Porosity(Index, theta) - 1.0) * gp_ele->q_R[gp];
			break;
	}

	return val;
}

/***************************************************************************
   GeoSys - Funktion:
   Assemble_RHS_TNEQ:
   11/2011   AKS
**************************************************************************/
void CFiniteElementStd::Assemble_RHS_TNEQ()
{
	int gp_r = 0, gp_s = 0, gp_t = 0;
	for (int i = 0; i < pcs->dof * nnodes; i++)
		NodalVal[i] = 0.0;

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		const double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		// Compute geometry
		ComputeShapefct(1);
		ComputeGradShapefct(1); // this is needed for CalCoef_RHS_TNEQ() !!

		for (int ii = 0; ii < pcs->dof; ii++)
		{
			const double fac = CalCoef_RHS_TNEQ(ii);
			for (int i = 0; i < nnodes; i++)
				NodalVal[i + ii * nnodes] += fac * fkt * shapefct[i];
		}
	}

	for (int ii = 0; ii < pcs->dof; ii++)
	{
		const long i_sh = NodeShift[ii];
		const int ii_sh = ii * nnodes;
		for (int i = 0; i < nnodes; i++)
		{
			eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i + ii_sh];
			(*RHS)[i + LocalShift + ii_sh] += NodalVal[i + ii_sh];
		}
	}
}

/**************************************************************************
    FEMLib-Method:
    Task: Calculate coefficient of RHS

    Precondition: The gradient of the shape function must have been computed
                  beforehand, because it is used for certain components.

    Programing:
    07/2013 TN
    04/2015 CL
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_TNEQ(const int dof_index)
{
	double const* const p0 = NodalVal0;
	double const* const p1 = NodalVal1;
	double const* const T0 = NodalVal_t0;
	double const* const T1 = NodalVal_t1;
	double const* const X0 = NodalVal_X0;
	double const* const X1 = NodalVal_X1;
	double const* const T0s = NodalVal_t2_0;
	double const* const T1s = NodalVal_t2_1;

	double& pg = eos_arg[0];
	double& Tg = eos_arg[1];
	double& Xw = eos_arg[2];

	const double theta = pcs->m_num->ls_theta;

	const int Index = MeshElement->GetIndex();
	poro = MediaProp->Porosity(Index, theta);
	const ElementValue* gp_ele = ele_gp_value[Index];
	const double q_r = gp_ele->q_R[gp];

	double val = 0.0;

	// NodalVal0 and NodalVal1 is the pressure on previous and current time step.
	//
	switch (dof_index)
	{
		case 0:
			val = (poro - 1.0) * q_r;
			break;

		case 1:
			pg = ipol(p0, p1, theta, this);
			Tg = ipol(T0, T1, theta, this);
			Xw = ipol(X0, X1, theta, this);

			val = FluidProp->Density(eos_arg) * poro * FluidProp->specific_heat_source;
			break;

		case 2:
		{
			pg = ipol(p0, p1, theta, this);
			Tg = ipol(T0, T1, theta, this);
			const double Ts = ipol(T0s, T1s, theta, this);
			Xw = ipol(X0, X1, theta, this);
			// end of adding eos_arg-----------------
			double H_vap = - SolidProp->reaction_enthalpy; //sign convention: defined negative for exothermic composition reaction but equ. written as: AB + \Delta H <--> A + B
			//Correction for temperature-dependent enthalpy in case of variable cpS
			const double rhoSR = gp_ele->rho_s_curr[gp];
			const double dcp_drhoSR((((*SolidProp->data_Capacity)(1)*SolidProp->upper_solid_density_limit -
							(*SolidProp->data_Capacity)(0)*SolidProp->lower_solid_density_limit)/
						(SolidProp->upper_solid_density_limit-SolidProp->lower_solid_density_limit) -
						(*SolidProp->data_Capacity)(0)) * SolidProp->lower_solid_density_limit/(rhoSR*rhoSR));

			const double cpS = SolidProp->Heat_Capacity(rhoSR);
			const double cpG = FluidProp->SpecificHeatCapacity(eos_arg);
			H_vap -= (cpS - cpG + rhoSR * dcp_drhoSR)*(Ts - 573.15);//TODO: Move IC to input file

			val = (1.0 - poro) * q_r * H_vap;
			val += gp_ele->rho_s_curr[gp] * (1.0 - poro) * SolidProp->specific_heat_source;

			if (MediaProp->getFrictionPhase() == FiniteElement::SOLID)
			{
				// HS, implementing the friction term here.

				const double* tensor = MediaProp->PermeabilityTensor(Index); // pointer to permeability tensor;
				double k_rel = 1.0; // relative permeability;

				// HS, added to get viscosity------------
				pg = ipol(p0, p1, theta, this);
				Tg = ipol(T0, T1, theta, this);
				Xw = ipol(X0, X1, theta, this);
				// end of adding eos_arg-----------------

				if (MediaProp->flowlinearity_model > 0)
				{
					k_rel = MediaProp->NonlinearFlowFunction(Index, gp, theta, this);
				}

				double* grad_pg = new double[dim]; // gradient of gas pressure.

				for (size_t i = 0; i < dim; i++) // loop over all dimensions
				{
					grad_pg[i] = 0.0; // clear to zero;
					for (int j = 0; j < nnodes; j++) // loop over all connecting nodes
					{
						double pg_tmp
						    = (1.0 - theta) * NodalVal0[j] + theta * NodalVal1[j]; // tmp value of gas pressure;
						int index_tmp = i * nnodes + j;
						grad_pg[i] += dshapefct[index_tmp] * pg_tmp;
					}
				}

				double friction_term = 0.0; // friction = poro \cdot grad_pg \cdot v_gas;
				// double * vel_Darcy = new double[dim]; // velocity term;

				for (size_t i = 0; i < dim; i++)
				{
					double vel_Darcy = 0.0;
					for (size_t j = 0; j < dim; j++)
					{
						size_t index_tmp = i * dim + j;
						vel_Darcy += tensor[index_tmp] * grad_pg[i];
					}
					friction_term += vel_Darcy * grad_pg[i] * k_rel / FluidProp->Viscosity(eos_arg);
				}

				val += friction_term;

				delete[] grad_pg; // clean tmp memory
				// delete[] vel_Darcy; // clean the memory of velocity darcy
			}
			break;
		}

		case 3:
			val = (poro - 1.0) * q_r;
			break;
	}

	return val;
}
}
