/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*! \file extend of file fem_ele_std1.cpp
   The definitions of members of class CFiniteElementStd.
 */

// C++ STL
#include <cfloat>
// Method
#include "fem_ele_std.h"
#include "mathlib.h"
// Problems
#include "rf_mmp_new.h"

#include "pcs_dm.h"
#include "rfmat_cp.h"

#ifndef NEW_EQS
#include "matrix_routines.h"
#endif

// Solver
#ifdef NEW_EQS
#include "equation_class.h"
using Math_Group::CSparseMatrix;
#endif

#include "pcs_dm.h" // displacement coupled

using namespace std;
namespace FiniteElement
{
/*!
   \brief Compute the additional term of Jacobian

       for the Newton-Raphson method for the p-p scheme.

     07.2011. WW

 */
void CFiniteElementStd::ComputeAdditionalJacobi_H2()
{
	int l; //, m;
	// int dm_shift = problem_dimension_dm;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt;
	double perturb = sqrt(DBL_EPSILON);
	double *tensor, *p2, *p2_0;
	double S1, vsc1, vsc2;
	double dkdp1, dkdp2;
	// double phi_dP_dt, d_ds_dp, phi_dP_dt_g;

	double relax = pcs->m_num->nls_relaxation;

	double gradPw[3], gradPg[3];
	double f_buff;
	double vw[3], vg[3];
	const double g_constant = 9.81;

	double dens_arg[3];
	dens_arg[1] = 293.15;

	p2_0 = NodalVal0 + nnodes;
	p2 = NodalVal1 + nnodes;

	//*StiffMatrix = 0.;

	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = relax * GetGaussData(gp, gp_r, gp_s, gp_t);
		ComputeShapefct(1); // Linear interpolation function
		ComputeGradShapefct(1); // Linear interpolation function

		// poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
		tensor = MediaProp->PermeabilityTensor(Index);
		PG = interpolate(NodalVal1);
		PG2 = interpolate(p2);
		Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
		S1 = Sw + perturb; // MediaProp->SaturationCapillaryPressureFunction(PG+perturb);

		dens_arg[0] = PG;
		rhow = FluidProp->Density(dens_arg);
		dens_arg[0] = PG2;
		rho_ga = GasProp->Density(dens_arg);
		vsc1 = FluidProp->Viscosity();
		vsc2 = GasProp->Viscosity();

		// dSdp = MediaProp->SaturationPressureDependency(Sw);
		dSdp = (MediaProp->SaturationCapillaryPressureFunction(PG + perturb) - Sw) / perturb;

		// Velocity
		for (size_t i = 0; i < dim; i++)
		{
			gradPw[i] = 0.0;
			gradPg[i] = 0.;
			for (int j = 0; j < nnodes; j++)
			{
				gradPw[i] += (p2[j] - NodalVal1[j]) * dshapefct[i * nnodes + j];
				gradPg[i] += p2[j] * dshapefct[i * nnodes + j];
			}
		}

		if ((coordinate_system) % 10 == 2)
		{
			gradPw[dim - 1] += g_constant * rhow;
			gradPg[dim - 1] += g_constant * rho_ga;
		}

		dkdp1 = dSdp
		        * (MediaProp->PermeabilitySaturationFunction(S1, 0) - MediaProp->PermeabilitySaturationFunction(Sw, 0))
		        / perturb;
		dkdp2 = dSdp
		        * (MediaProp->PermeabilitySaturationFunction(S1, 1) - MediaProp->PermeabilitySaturationFunction(Sw, 1))
		        / perturb;

		for (size_t i = 0; i < dim && i < 3; i++)
		{
			vw[i] = 0.0;
			vg[i] = 0.;
			const size_t ish = i * dim;
			for (size_t j = 0; j < dim; j++)
			{
				vw[i] += tensor[ish + j] * gradPw[j];
				vg[i] += tensor[ish + j] * gradPg[j];
			}
			vw[i] *= dkdp1 * time_unit_factor / vsc1;
			vg[i] *= dkdp2 * rho_ga * time_unit_factor / (vsc2 * rhow);
		}

		/// For the Laplace
		for (int i = 0; i < nnodes; i++)
		{
			l = i + nnodes;
			for (int j = 0; j < nnodes; j++)
				// m = j+nnodes;
				for (size_t k = 0; k < dim; k++)
				{
					f_buff = fkt * dshapefct[k * nnodes + i] * shapefct[j];
					(*StiffMatrix)(i, j) += f_buff * vw[k];
					(*StiffMatrix)(l, j) += f_buff * vg[k];
				}
		}

#define Take_Deformation_to_Jacobian
#ifdef Take_Deformation_to_Jacobian

		// d(dS/dp)/dS  may not be available.
		// Instead dS = dS(p)/dp * dp
		/*
		   // Mass related
		   /// d(dS/dp)/dS
		   d_ds_dp = ( MediaProp->SaturationPressureDependency(S1) - dSdp)/perturb;

		   phi_dP_dt = 0.;
		   phi_dP_dt_g = 0.;
		   for(i=0; i<nnodes; i++)
		   {
		   phi_dP_dt += (NodalVal1[i] - NodalVal0[i]) * shapefct[i];
		   phi_dP_dt_g += (p2[i] - p2_0[i]) * shapefct[i];
		   }
		   phi_dP_dt *= poro*dSdp*d_ds_dp/dt;
		   phi_dP_dt_g *= poro*dSdp*d_ds_dp/dt;
		 */

		double ddens_g_dt;
		dens_arg[0] = PG2;
		rho_ga = GasProp->Density(dens_arg);
		dens_arg[0] = interpolate(p2_0);
		ddens_g_dt = -poro * dSdp * (rho_ga - GasProp->Density(dens_arg)) / dt;
		ddens_g_dt /= rhow;

		if (dm_pcs)
		{
			// setOrder(2);
			// GetGaussData(gp, gp_r, gp_s, gp_t);
			// ComputeGradShapefct(2);
			// setOrder(1);

			/// if deformation is coupled
			vw[0] = 0.; // Here for dSdp*grad u/dt
			for (int i = 0; i < nnodes; i++)
			//            for (i=0;i<nnodesHQ;i++)
			{
				vw[0] += NodalVal2[i] * dshapefct[i] + NodalVal3[i] * dshapefct[i + nnodes];
				//               vw[0]  += NodalVal2[i]*dshapefctHQ[i]+NodalVal3[i]*dshapefctHQ[i+nnodesHQ];
				if (dim == 3) // 3D.
					//                  vw[0]  +=   NodalVal4[i]*dshapefctHQ[2*nnodesHQ+i];
					vw[0] += NodalVal4[i] * dshapefct[2 * nnodes + i];
			}
			vw[0] *= dSdp / dt;
		}
		else
			vw[0] = 0.;

		for (int i = 0; i < nnodes; i++)
		{
			l = i + nnodes;
			for (int j = 0; j < nnodes; j++)
			{
				// m = j+nnodes;
				f_buff = fkt * shapefct[i] * shapefct[j];
				//(*StiffMatrix)(i,j) += f_buff*(phi_dP_dt+vw[0]);
				//(*StiffMatrix)(l,j) -= rho_ga*f_buff*(phi_dP_dt_g+vw[0])/rhow;

				(*StiffMatrix)(l, j) += f_buff * ddens_g_dt;

				(*StiffMatrix)(i, j) += f_buff * vw[0];
				(*StiffMatrix)(l, j) -= rho_ga * f_buff * vw[0] / rhow;
			}
		}
#endif
	} // loop gauss points

	// add2GlobalMatrixII(1);

	// StiffMatrix->Write();
}

//------------------------------------------------------------------
/*!
   \brief Compute the additional term of Jacobian

       for the Newton-Raphson method for the p-p scheme.

     07.2011. WW

 */
void CFiniteElementStd::ComputeAdditionalJacobi_Richards()
{
	// int dm_shift = problem_dimension_dm;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt; //, mat_fac;
	double perturb = sqrt(DBL_EPSILON);
	double* tensor;
	double S1, vsc1;
	double dkdp1;
	// double phi_dP_dt, d_ds_dp, phi_dP_dt_g;

	double relax = pcs->m_num->nls_relaxation;

	double gradPw[3];
	double vw[3];
	const double g_constant = 9.81;

	//*StiffMatrix = 0.;

	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = relax * GetGaussData(gp, gp_r, gp_s, gp_t);
		ComputeShapefct(1); // Linear interpolation function
		ComputeGradShapefct(1); // Linear interpolation function

		tensor = MediaProp->PermeabilityTensor(Index);
		PG = -interpolate(NodalVal1);
		Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
		S1 = Sw + perturb; // MediaProp->SaturationCapillaryPressureFunction(PG+perturb,0);

		vsc1 = FluidProp->Viscosity();

		// dSdp = MediaProp->SaturationPressureDependency(Sw);
		dSdp = (MediaProp->SaturationCapillaryPressureFunction(PG + perturb) - Sw) / perturb;

		// Velocity
		for (size_t i = 0; i < dim; i++)
		{
			gradPw[i] = 0.0;
			for (int j = 0; j < nnodes; j++)
				gradPw[i] += NodalVal1[j] * dshapefct[i * nnodes + j];
		}

		if ((coordinate_system) % 10 == 2)
			gradPw[dim - 1] += g_constant * rhow;

		dkdp1 = dSdp
		        * (MediaProp->PermeabilitySaturationFunction(S1, 0) - MediaProp->PermeabilitySaturationFunction(Sw, 0))
		        / perturb;

		for (size_t i = 0; i < dim && i < 3; i++)
		{
			vw[i] = 0.0;
			for (size_t j = 0; j < dim; j++)
				vw[i] += tensor[i * dim + j] * gradPw[j];

			vw[i] *= dkdp1 * time_unit_factor / vsc1;
		}

		/// For the Laplace
		for (int i = 0; i < nnodes; i++)
			for (int j = 0; j < nnodes; j++)
				// m = j+nnodes;
				for (size_t k = 0; k < dim; k++)
					(*StiffMatrix)(i, j) += fkt * dshapefct[k * nnodes + i] * shapefct[j] * vw[k];

#define Take_Deformation_to_Jacobian
#ifdef Take_Deformation_to_Jacobian

		if (dm_pcs)
		{
			// setOrder(2);
			// GetGaussData(gp, gp_r, gp_s, gp_t);
			// ComputeGradShapefct(2);
			// setOrder(1);

			/// if deformation is coupled
			vw[0] = 0.; // Here for dSdp*grad u/dt
			for (int i = 0; i < nnodes; i++)
			//            for (i=0;i<nnodesHQ;i++)
			{
				vw[0] += NodalVal2[i] * dshapefct[i] + NodalVal3[i] * dshapefct[i + nnodes];
				//               vw[0]  += NodalVal2[i]*dshapefctHQ[i]+NodalVal3[i]*dshapefctHQ[i+nnodesHQ];
				if (dim == 3) // 3D.
					//                  vw[0]  +=   NodalVal4[i]*dshapefctHQ[2*nnodesHQ+i];
					vw[0] += NodalVal4[i] * dshapefct[2 * nnodes + i];
			}
			vw[0] *= dSdp / dt;
		}
		else
			vw[0] = 0.;

		for (int i = 0; i < nnodes; i++)
			for (int j = 0; j < nnodes; j++)
				(*StiffMatrix)(i, j) += fkt * shapefct[i] * shapefct[j] * vw[0];

#endif
	} // loop gauss points

	// add2GlobalMatrixII(1);

	// StiffMatrix->Write();
}
} // namespace
