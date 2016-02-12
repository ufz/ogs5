/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   Class: RandomWalk
   Task: Random Walk - an alternative for FDM or FEM of transport equation
   Programing:
   07/2005 PCH Implementation
**************************************************************************/

#include "rf_random_walk.h"

#include "FileTools.h"
#include "Output.h"
#include "matrix_class.h"
#include "rf_fluid_momentum.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"

// C++ STL
#include <sstream>
#include <cfloat>
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Math_Group;

#define PCT_FILE_EXTENSION ".pct"
#define SWAP(x, y) \
	{              \
		double t;  \
		t = x;     \
		x = y;     \
		y = t;     \
	} // WW data type is change to double
//#define CountParticleNumber                            //YS: to count the number of particles leave the domain

/**************************************************************************
   Class: RandomWalk
   Task: constructor
   Programing:
   07/2005 PCH Implementation
   last modification:
**************************************************************************/
RandomWalk::RandomWalk(int srand_seed)
{
	m_pcs = NULL;
	fem = NULL;

	// This is going to be reset by user input.
	// Further, used for allocating dynamic memory.
	numOfParticles = 0;
	leavingParticles = 0; // YS: to count the number of particles in the outflow
	UniformOrNormal = 1; // Uniform random number generation
	RWPTMode = 0; // Initialized to be homogeneous media
	PURERWPT = 0;
	CurrentTime = 0.0;
	FDMIndexSwitch = 0;
	GridOption = 0;
	ChanceOfIrreversed = NULL; // YS: judgement for decay

	// To produce a different pseudo-random series each time your program is run.
	if (srand_seed == 0)
		srand((int)time(0));
	// To produce same pseudo-random series each time your program is run.
	else if (srand_seed == 1)
		srand(1);

	// These are the allowable outputs (input as options to file <file_base_name>.out
	rwpt_out_strings.push_back("PARTICLES"); // output particle locations
	// output particles as elemental concentration
	rwpt_out_strings.push_back("PARTICLE_CONCENTRATION");
	// FM_TEST
	flow_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* pcs_e = pcs_vector[i];
		if (pcs_e->getProcessType() == FiniteElement::GROUNDWATER_FLOW
		    || pcs_e->getProcessType() == FiniteElement::LIQUID_FLOW
		    || pcs_e->getProcessType() == FiniteElement::RICHARDS_FLOW
		    || pcs_e->getProcessType() == FiniteElement::MULTI_PHASE_FLOW
		    || pcs_e->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
		{
			flow_pcs = pcs_e;
			break;
		}
	}
}

// Constructor
Particle::Particle(void)
{
	// Position Vector
	x = y = z = 0.0;

	// Velocity Vector
	Vx = Vy = Vz = 0.0;
	K = 0.0;

	dVxdx = dVydy = dVzdz = 0.0;

	for (int i = 0; i < 9; ++i)
		D[i] = 0.0;

	// Time
	t = 0.0;
	identity = 0;
	on_boundary = 0; // YS
}

/**************************************************************************
   Class: RandomWalk
   Task: destructor
   Programing:
   07/2005 PCH Implementation
   last modification:
**************************************************************************/
RandomWalk::~RandomWalk(void)
{
	//    if(X) delete [] X;
	if (ChanceOfIrreversed)
		delete[] ChanceOfIrreversed;
	ChanceOfIrreversed = NULL;
}

/**************************************************************************
   Class: RandomWalk
   Task: Create a random number from N(0,1) distribution
      Marsaglia Algorithm adopted...
      Refer Math 4255 Assignment
      NOTE:
      To make a use of this function, srand() should called beforehand.
      Otherwise, the random number keeps continuing the same sequence.
   Programing:
   08/2005 PCH Implementation
   last modification:
**************************************************************************/
double RandomWalk::Marsaglia(void)
{
	int whichOne = 0;
	double u1 = 0.0, u2 = 0.0;
	double v1 = 0.0, v2 = 0.0;
	double s = 0.0;

	do
	{
		// Create two random numbers which are uniform between 0 to 1.
		u1 = (double)(1.0 * rand() / (RAND_MAX + 1.0));
		u2 = (double)(1.0 * rand() / (RAND_MAX + 1.0));

		v1 = 2. * u1 - 1.0;
		v2 = 2. * u2 - 1.0;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.0 || s == 0.0); // To fit log definition

	double fac = sqrt(-2.0 * log(s) / s);

	// This will create either 0 or 1.
	whichOne = (int)(1.0 * rand() / (RAND_MAX + 1.0));

	if (whichOne == 0)
		return v1 * fac;
	else
		return v2 * fac;
}

double RandomWalk::randomMinusOneToOne(void)
{
	return (double)(2.0 * rand() / (RAND_MAX + 1.0) - 1.0);
}

double RandomWalk::randomZeroToOne(void)
{
	return (double)(1.0 * rand() / (RAND_MAX + 1.0));
}

CFEMesh* RandomWalk::selectMeshForFluidMomentumProcess()
{
	size_t n_pcs(pcs_vector.size());
	FiniteElement::ProcessType pcs_type;
	CFEMesh* msh(NULL);

	for (size_t i = 0; i < n_pcs; ++i)
	{
		pcs_type = pcs_vector[i]->getProcessType();
		// Select the mesh whose process name has the mesh for Fluid_Momentum
		if (pcs_type == FiniteElement::RICHARDS_FLOW)
			msh = FEMGet("RICHARDS_FLOW");
		else if (pcs_type == FiniteElement::LIQUID_FLOW)
			msh = FEMGet("LIQUID_FLOW");
		else if (pcs_type == FiniteElement::GROUNDWATER_FLOW)
			msh = FEMGet("GROUNDWATER_FLOW");
	}

	if (!msh) // FM_TEST
		msh = FEMGet("FLUID_MOMENTUM");
	return msh;
}

/**************************************************************************
   Class: RandomWalk
   Task: This function interpolates velocity in reference space
   Programing:
   01/2007 PCH Implementation
   10/2010 TF changed access to process type
**************************************************************************/
void RandomWalk::InterpolateVelocity(Particle* A)
{
	// Get the mesh first
	// TF
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());

	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos) TF
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->getProcessType () == LIQUID_FLOW)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->getProcessType () == GROUNDWATER_FLOW)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	// Mount the element fromthe first particle from particles initially
	MeshLib::CElem* theEle = m_msh->ele_vector[A->elementIndex];
	// OK411 double tolerance = 1e-8;

	// If a quad element,
	int nnode = theEle->GetEdgesNumber();
	m_pcs = PCSGet("FLUID_MOMENTUM");

	// FM_TEST
	int idx, idy, idz;
	if (m_pcs)
	{
		idx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
		idy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
		idz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;
	}
	else if (flow_pcs)
	{
		idx = flow_pcs->GetNodeValueIndex("VELOCITY_X1");
		idy = flow_pcs->GetNodeValueIndex("VELOCITY_Y1");
		idz = flow_pcs->GetNodeValueIndex("VELOCITY_Z1");

		if (m_msh->GetCoordinateFlag() / 10 == 1)
		{
			if (m_msh->GetCoordinateFlag() == 11)
			{
				int ibuff = idy;
				idy = idx;
				idx = ibuff;
			}
			if (m_msh->GetCoordinateFlag() == 12)
			{
				int ibuff = idz;
				idz = idx;
				idx = ibuff;
			}
		}

		m_pcs = flow_pcs;
	}

	/*
	// Let's solve pore velocity.
	// It is simple because Sw stuff automatically handles in Richards Flow.
	// Thus, I only divide Darcy velocity by porosity only to get pore velocity.
	CMediumProperties *MediaProp = mmp_vector[theEle->GetPatchIndex()];
	double porosity = 0.0;
	if(MediaProp->porosity > 10-6)
	    porosity = MediaProp->porosity;             // This is for simple one.
	else
	    // This will get you porosity.
	    porosity = MediaProp->porosity_model_values[0];
	*/

	// I guess for Dual Porosity stuff,
	// this code should be revisited.

	if (nnode == 4)
	{
		// Get physical coordinates of four corner points
		double x[4], y[4] /*, z[4]*/;
		// double vx[4], vy[4], vz[4];
		for (int i = 0; i < nnode; ++i)
		{
			double const* const pnt(theEle->GetNode(i)->getData());
			x[i] = pnt[0];
			y[i] = pnt[1];
			// z[i] = pnt[2];

			// vx[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idx)/porosity; //FM_TEST
			// vy[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idy)/porosity;
			// vz[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idz)/porosity;
		}

		// solve for Jm at xm = (xhat,yhat)=(1/2,1/2) <- RT0
		double Jm;
		Jm = 0.5 * (x[0] * y[1] - y[0] * x[1] + x[1] * y[2] - y[1] * x[2] + y[0] * x[3] - x[0] * y[3] + x[2] * y[3]
		            - y[2] * x[3]);

		// Mount the edges of the element
		vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
		theEle->GetEdges(theEdgesOfThisElement);
		// Mount the nodes of the edge
		vec<MeshLib::CNode*> theNodesOfThisEdge(3);

		// Solve for proper index in reference space
		double const* Ecenter = theEle->GetGravityCenter();
		double u[4]; // Edge flux in reference space
		for (int i = 0; i < 4; ++i)
			u[i] = 0.0;

		for (int i = 0; i < 4; ++i)
		{
			// Get the nodes of the edge i
			theEdgesOfThisElement[i]->GetNodes(theNodesOfThisEdge);

			double node1[3] = {(theNodesOfThisEdge[0]->getData())[0], (theNodesOfThisEdge[0]->getData())[1], 0.0};
			double node2[3] = {(theNodesOfThisEdge[1]->getData())[0], (theNodesOfThisEdge[1]->getData())[1], 0.0};
			// TF         node1[0] = theNodesOfThisEdge[0]->X(); node1[1] = theNodesOfThisEdge[0]->Y();
			// TF        node2[0] = theNodesOfThisEdge[1]->X(); node2[1] = theNodesOfThisEdge[1]->Y();
			// Get the referece position of these two ending points of the edge
			IsoparametricMappingQuadfromPtoR(A->elementIndex, node1);
			IsoparametricMappingQuadfromPtoR(A->elementIndex, node2);

			double dx, dy;
			dx = theEdgesOfThisElement[i]->GetVelocity(0);
			dy = theEdgesOfThisElement[i]->GetVelocity(1);

			double edge2mid[3], EdgeMidPoint[3];
			theEdgesOfThisElement[i]->GetEdgeMidPoint(EdgeMidPoint);
			for (int j = 0; j < 3; ++j)
				edge2mid[j] = Ecenter[j] - EdgeMidPoint[j];

			double angle
			    = acos((dx * edge2mid[0] + dy * edge2mid[1])
			           / (sqrt(dx * dx + dy * dy) * sqrt(edge2mid[0] * edge2mid[0] + edge2mid[1] * edge2mid[1])));

			// EA and EB: x1hat = x2hat and x1hat*x2hat > 0
			if (fabs(node1[0] - node2[0]) < 0.5 && (node1[0] * node2[0]) > 0.0)
			{
				// EA: x1hat < 0
				if (node1[0] < 0.0)
				{
					if (angle > 3.141592 / 2.0)
						u[0] = -sqrt(dx * dx + dy * dy) / Jm;
					else
						u[0] = sqrt(dx * dx + dy * dy) / Jm;
				}
				else // EB
				{
					if (angle > 3.141592 / 2.0)
						u[1] = sqrt(dx * dx + dy * dy) / Jm;
					else
						u[1] = -sqrt(dx * dx + dy * dy) / Jm;
				}
			}
			// Else then, EC and ED
			else
			{
				// EC: y1hat < 0
				if (node1[1] < 0.0)
				{
					if (angle > 3.141592 / 2.0)
						u[2] = -sqrt(dx * dx + dy * dy) / Jm;
					else
						u[2] = sqrt(dx * dx + dy * dy) / Jm;
				}
				else // ED
				{
					if (angle > 3.141592 / 2.0)
						u[3] = sqrt(dx * dx + dy * dy) / Jm;
					else
						u[3] = -sqrt(dx * dx + dy * dy) / Jm;
				}
			}
		}

		// Let's do Pollock's method here.
		// Solve for the reference position of this particle;
		double R[3];
		R[0] = A->x;
		R[1] = A->y;
		R[2] = A->z;
		IsoparametricMappingQuadfromPtoR(A->elementIndex, R);

		// Find the correct entry face from here.
		double Gx = u[1] - u[0];
		double Gy = u[3] - u[2];
		double uT = Gx * R[0] + u[0];
		double vT = Gy * R[1] + u[2];

		double t[4], bigTime = 1e10, tmin;

		t[0] = TA(Gx, uT, u[0], R[0]);
		if (t[0] < 0.0)
			t[0] = bigTime;
		t[1] = TB(Gx, uT, u[0], u[1], R[1]);
		if (t[1] < 0.0)
			t[1] = bigTime;
		t[2] = TA(Gy, vT, u[2], R[1]);
		if (t[2] < 0.0)
			t[2] = bigTime;
		t[3] = TB(Gy, vT, u[2], u[3], R[1]);
		if (t[3] < 0.0)
			t[3] = bigTime;

		int idx = -10;
		tmin = Tmin(t, &idx);

		R[0] = 1.0 / Gx * (exp(-Gx * tmin) * uT - u[0]);
		R[1] = 1.0 / Gy * (exp(-Gy * tmin) * vT - u[2]);
		R[2] = 0.0;

		// Solve exit point in physical space back.
		IsoparametricMappingQuadfromRtoP(A->elementIndex, R);

		A->x = R[0];
		A->y = R[1];
		A->z = R[2];
	}
}

/**************************************************************************
   Class: RandomWalk
   Task: This function traces exit point in this element.
     The element should be regular.
     Not really working.
   Programing:
   02/2007 PCH Implementation
   03/2010 JT modified cell search algorithm
**************************************************************************/
void RandomWalk::TracePathlineInThisElement(Particle* A)
{
	// Get the mesh first
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	// Mount the element fromthe first particle from particles initially
	MeshLib::CElem* theEle = m_msh->ele_vector[A->elementIndex];
	// double tolerance = 1e-8;

	// If a quad element,
	int nnode = theEle->GetEdgesNumber();
	m_pcs = PCSGet("FLUID_MOMENTUM");
	if (!m_pcs) // FM_TEST
		m_pcs = flow_pcs;

	// Let's solve pore velocity.
	// It is simple because Sw stuff automatically handles in Richards Flow.
	// Thus, I only divide Darcy velocity by porosity only to get pore velocity.
	// WW CMediumProperties *MediaProp = mmp_vector[theEle->GetPatchIndex()];
	/*   //WW
	   double porosity = 0.0;
	   if(MediaProp->porosity > 10-6)
	   porosity = MediaProp->porosity;             // This is for simple one.
	   else
	                                               // This will get you porosity.
	      porosity = MediaProp->porosity_model_values[0];
	 */
	// I guess for Dual Porocity stuff,
	// this code should be revisited.

	if (nnode == 4)
	{
		// Mount the edges of the element
		vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
		theEle->GetEdges(theEdgesOfThisElement);
		// Mount the nodes of the edge
		vec<MeshLib::CNode*> theNodesOfThisEdge(3);

		// Let's do Pollock's method here.
		// Solve for the reference position of this particle;
		double R[3];
		R[0] = A->x;
		R[1] = A->y;
		R[2] = A->z;

		double x0, x1, y0, y1;
		// Get the nodes of the edge i
		theEdgesOfThisElement[0]->GetNodes(theNodesOfThisEdge);
		//	x0 = theNodesOfThisEdge[0]->X();
		x0 = R[0];
		double const* const pnt(theNodesOfThisEdge[1]->getData());
		x1 = pnt[0];
		theEdgesOfThisElement[1]->GetNodes(theNodesOfThisEdge);
		//	y0 = theNodesOfThisEdge[0]->Y();
		y0 = R[1];
		y1 = pnt[1];

		double Fx0, Fx1, Fy0, Fy1; // Flux at each edge in physical space
		Fx0 = theEdgesOfThisElement[3]->GetVelocity(0);
		Fx1 = theEdgesOfThisElement[1]->GetVelocity(0);
		Fy0 = theEdgesOfThisElement[0]->GetVelocity(1);
		Fy1 = theEdgesOfThisElement[2]->GetVelocity(1);

		double ax, bx, cx, ay, by, cy;
		ax = Fx1 - Fx0;
		bx = Fx0 * x1 - x0 * x1;
		cx = x1 - x0;
		ay = Fy1 - Fy0;
		by = Fy0 * y1 - y0 * y1;
		cy = y1 - y0;

		// Solve the T
		double tolerance = 1e-8;
		double tx, ty;
		if (fabs(ax) > tolerance)
			tx = cx / ax * log((ax * x1 + bx) / (ax * x0 + bx));
		else
			abort();

		if (fabs(ay) > tolerance)
			ty = cy / ay * log((ay * y1 + by) / (ay * y0 + by));
		else
			abort();

		double tmin;
		if (tx > ty)
			tmin = ty;
		else
			tmin = tx;

		if (cx > tolerance)
			R[0] = (x0 + bx / ax) * exp(ax / cx * tmin) - bx / ax;

		if (cy > tolerance)
			R[1] = (y0 + by / ay) * exp(ay / cy * tmin) - by / ay;

		R[2] = pnt_z_min;
		A->x = R[0];
		A->y = R[1];
		A->z = R[2];
		// update the element index // JT 2010;
		long i, j, k;
		i = (int)((A->x - pnt_x_min) / dx);
		j = (int)((A->y - pnt_y_min) / dy);
		k = (int)((A->z - pnt_z_min) / dz);
		long iFDM = k * (nx * ny) + j * nx + i;

		if ((size_t)iFDM < indexFDM.size())
			A->elementIndex = indexFDM[iFDM].eleIndex;
		else
			A->elementIndex = -10; // Outside of the domain
	}
}

double RandomWalk::Tmin(double* a, int* idx)
{
	double Tmin = 1e20;

	for (int i = 0; i < 4; ++i)
		if (Tmin > a[i])
		{
			Tmin = a[i];
			*idx = i;
		}

	return Tmin;
}

double RandomWalk::TA(double Gx, double uT, double uA, double xTA)
{
	double tolerance = 1e-10, bigTime = 1e10;

	if (fabs(Gx) > tolerance) // If Gx is not zero
	{
		if (fabs(uA) > tolerance)
		{
			if (uA * uT > 0.0)
			{
				if (1.0 / Gx * log(uT / uA) < 0.0)
					return bigTime;
				else
					return 1.0 / Gx * log(uT / uA);
			}
			else
				return bigTime;
		}
		else
			return bigTime; // Assume T is infinite.
	}
	else
	{
		if (fabs(uA) > tolerance) // uA is not zero
		{
			double xA = -1.0; // In my reference space
			if ((xTA - xA) * uA < 0.0)
				return bigTime;
			else
				return xTA / uA;
		}
		else
			return bigTime;
	}
}

double RandomWalk::TB(double Gx, double uT, double uA, double uB, double xTB)
{
	double tolerance = 1e-10, bigTime = 1e10;

	if (fabs(Gx) > tolerance) // If Gx is not zero
	{
		if (fabs(uB) > tolerance)
		{
			if (uB * uT > 0.0)
			{
				if (1.0 / Gx * log(uT / uB) < 0.0)
					return bigTime;
				else
					return 1.0 / Gx * log(uT / uB);
			}
			else
				return bigTime;
		}
		else
			return bigTime; // Assume T is infinite.
	}
	else
	{
		if (fabs(uB) > tolerance) // uA is not zero
		{
			double xB = 1.0; // In my reference space
			if ((xTB - xB) * uB < 0.0)
				return bigTime;
			else
				return (xTB - 1.0) / uA;
		}
		else
			return bigTime;
	}
}

/**************************************************************************
   Class: RandomWalk
   Task: This function interpolates velocity of the particle
     based on the inverse distance method
     RWPT-IM This function should OK with the real plane.
   Programing:
   10/2005 PCH Implementation
   02/2006 PCH The function is updated to solve for velocity in the element
         that has a joint or crossroads.
   05/2006 PCH This one gets hydraulic conductivity as well.
   last modification:
**************************************************************************/
void RandomWalk::InterpolateVelocityOfTheParticleByInverseDistance(Particle* A)
{
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();

	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];

	// Let's get the hydraulic conductivity first.
	CMediumProperties* MediaProp = mmp_vector[m_ele->GetPatchIndex()];
	// OK411 int phase = 0;
	CFluidProperties* FluidProp = mfp_vector[0];
	double* kTensor = MediaProp->PermeabilityTensor(A->elementIndex);
	double k = kTensor[0];

	A->K = k * FluidProp->Density() * 9.81 / FluidProp->Viscosity();

	// Get the number of nodes
	int nnodes = m_ele->GetVertexNumber();
	// Allocate the memory accordingly
	Particle* vertex = NULL;
	vertex = new Particle[nnodes];
	double* d = NULL;
	d = new double[nnodes];
	double SumOfdInverse = 0.0;

	// Get the cooridinate of the nodes in the element
	for (int i = 0; i < nnodes; ++i)
	{
		double const* const coords(m_ele->GetNode(i)->getData());
		vertex[i].x = coords[0];
		vertex[i].y = coords[1];
		vertex[i].z = coords[2];

		// Compute the each distance
		double x = vertex[i].x - A->x;
		double y = vertex[i].y - A->y;
		double z = vertex[i].z - A->z;
		d[i] = sqrt(x * x + y * y + z * z);
		SumOfdInverse += 1.0 / d[i];
	}

	// Let's get the weight of each node
	double* w = NULL;
	w = new double[nnodes];
	// Initialize the velocity
	A->Vx = A->Vy = A->Vz = 0.0;

	// FM_TEST
	m_pcs = PCSGet("FLUID_MOMENTUM");
	int idx = -1, idy = -1, idz = -1;
	if (m_pcs)
	{
		idx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
		idy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
		idz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;
	}
	else if (flow_pcs)
	{
		idx = flow_pcs->GetNodeValueIndex("VELOCITY_X1");
		idy = flow_pcs->GetNodeValueIndex("VELOCITY_Y1");
		idz = flow_pcs->GetNodeValueIndex("VELOCITY_Z1");
		m_pcs = flow_pcs;

		if (m_msh->GetCoordinateFlag() / 10 == 1)
		{
			if (m_msh->GetCoordinateFlag() == 11)
			{
				int ibuff = idy;
				idy = idx;
				idx = ibuff;
			}
			if (m_msh->GetCoordinateFlag() == 12)
			{
				int ibuff = idz;
				idz = idx;
				idx = ibuff;
			}
		}
	}
	for (int i = 0; i < nnodes; ++i)
	{
		w[i] = 1.0 / (d[i] * SumOfdInverse);

		double vx = 0.0, vy = 0.0, vz = 0.0;
		// If this node is crossroad,
		if (m_msh->nod_vector[m_ele->GetNodeIndex(i)]->crossroad)
		{
			// Get the velocity contributed in this element
			CrossRoad* crossroad = NULL;
			for (int j = 0; j < (int)(m_msh->fm_pcs->crossroads.size()); ++j)
				if ((size_t)m_msh->fm_pcs->crossroads[j]->Index
				    == m_msh->nod_vector[m_ele->GetNodeIndex(i)]->GetIndex())
					crossroad = m_msh->fm_pcs->crossroads[j];

			if (crossroad)
			{
			}
			else // Failed to find the crossroad although it is a crossroad
				abort();

			// Find the velocity of the crossroad associated with the connected planes.
			for (int k = 0; k < crossroad->numOfThePlanes; ++k)
			{
				// I am going to check the normal vector of the element and the connected plane.
				double tolerance = 1e-10;
				double E[3], P[3];
				for (int p = 0; p < 3; ++p)
				{
					E[p] = m_ele->getTransformTensor(6 + p);
					P[p] = crossroad->plane[k].norm[p];
				}

				double same
				    = (E[0] - P[0]) * (E[0] - P[0]) + (E[1] - P[1]) * (E[1] - P[1]) + (E[2] - P[2]) * (E[2] - P[2]);

				if (same < tolerance)
				{
					vx = crossroad->plane[k].V[0];
					vy = crossroad->plane[k].V[1];
					vz = crossroad->plane[k].V[2];
				}
			}
		}
		else
		{
			vx = m_pcs->GetNodeValue(m_ele->GetNodeIndex(i), idx); // FM_TEST
			vy = m_pcs->GetNodeValue(m_ele->GetNodeIndex(i), idy);
			vz = m_pcs->GetNodeValue(m_ele->GetNodeIndex(i), idz);

			// Let's solve pore velocity.
			// It is simple because Sw stuff automatically handles in Richards Flow.
			// Thus, I only divide Darcy velocity by porosity only to get pore velocity.
			CMediumProperties* MediaProp = mmp_vector[m_ele->GetPatchIndex()];
			double porosity = 0.0;
			// This will get you porosity.
			porosity = MediaProp->porosity_model_values[0];
			// I guess for Dual Porocity stuff,
			// this code should be revisited.

			vx /= porosity;
			vy /= porosity;
			vz /= porosity;
		}

		A->Vx += w[i] * vx;
		A->Vy += w[i] * vy;
		A->Vz += w[i] * vz;
	}

	// Release the temperary memory in this function
	delete[] vertex;
	delete[] d;
	delete[] w;
}

/**************************************************************************
   Class: RandomWalk
   Task: This function interpolates velocity of the particle
     based on the bilinear method
     The function requires FDM-like grid.
     option 0: FDM method
     option 1: Transformation for irregular grids
   Programing:
   01/2007 PCH Implementation
   02/2007 PCH Modification
   last modification:
**************************************************************************/
void RandomWalk::InterpolateVelocityOfTheParticleByBilinear(int option, Particle* A)
{
	// Get the element that the particle belongs
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	// Let's allocate some memory for miniFEM
	CFEMesh* m_mini(new CFEMesh(m_msh->getGEOObjects(), m_msh->getProjectName()));
	MeshLib::CElem* miniEle = new MeshLib::CElem();

	if (option == 0) // FDM method
	{
		int eleIndex = IndexOfTheElementThatThisParticleBelong(0, A);

		if (eleIndex == -5)
			// This is a temperary measure.
			// Gotta be written better later on.
			eleIndex = 0;

		// Element is outside of the domain or in the sink
		// Do not do anything. If not, then proceed.
		if (eleIndex != -10)
		{
			A->elementIndex = eleIndex;
			MeshLib::CElem* m_ele = m_msh->ele_vector[eleIndex];

			int nnode = m_ele->GetEdgesNumber();

			// Let's get the hydraulic conductivity first.
			CMediumProperties* MediaProp = mmp_vector[m_ele->GetPatchIndex()];
			CFluidProperties* FluidProp = mfp_vector[0];
			double* kTensor = MediaProp->PermeabilityTensor(eleIndex);
			double k = kTensor[0];

			A->K = k * FluidProp->Density() * 9.81 / FluidProp->Viscosity();

			// Get porosity
			// I guess for Dual Porocity stuff this code should be revisited.
			double porosity = 0.0;
			if (MediaProp->porosity > 10 - 6)
				porosity = MediaProp->porosity; // This is for simple one.
			else
				// This will get you porosity.
				porosity = MediaProp->porosity_model_values[0];

			// Get the number of nodes
			// int nnodes = m_ele->GetVertexNumber();

			// Mount the edges of the element
			vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
			m_ele->GetEdges(theEdgesOfThisElement);

			double E1[3], E2[3], E3[3], E4[3]; // mid points of each edge
			theEdgesOfThisElement[0]->GetEdgeMidPoint(E1);
			theEdgesOfThisElement[1]->GetEdgeMidPoint(E2);
			theEdgesOfThisElement[2]->GetEdgeMidPoint(E3);
			theEdgesOfThisElement[3]->GetEdgeMidPoint(E4);

			// ux = a+b(x-x0); uy = c + d(y-y0)

			double vx1, vx2, vy1, vy2, a, b, x0, y0;

			vx2 = theEdgesOfThisElement[1]->GetVelocity(0);
			vx1 = theEdgesOfThisElement[3]->GetVelocity(0);
			vy2 = theEdgesOfThisElement[0]->GetVelocity(1);
			vy1 = theEdgesOfThisElement[2]->GetVelocity(1);

			x0 = E4[0];
			y0 = E1[1];
			a = (vx2 - vx1) / dx;
			b = (vy1 - vy2) / dy;

			// Let's solve pore velocity.
			// It is simple because Sw stuff automatically handles in Richards Flow.
			// Thus, I only divide Darcy velocity by porosity only to get pore velocity.
			A->Vx = (vx1 + a * (A->x - x0)) / porosity;
			A->Vy = (vy2 + b * (A->y - y0)) / porosity;
			A->Vz = 0.0;
		}
		// else
		//	;	// Think later for out of boundary particles.
	}
	else if (option == 1) // miniFEM Way
	{
		m_pcs = PCSGet("FLUID_MOMENTUM");
		// FM_TEST
		int idx, idy, idz;
		if (m_pcs)
		{
			idx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
			idy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
			idz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;
		}
		else if (flow_pcs)
		{
			idx = flow_pcs->GetNodeValueIndex("VELOCITY_X1");
			idy = flow_pcs->GetNodeValueIndex("VELOCITY_Y1");
			idz = flow_pcs->GetNodeValueIndex("VELOCITY_Z1");
			m_pcs = flow_pcs;

			if (m_msh->GetCoordinateFlag() / 10 == 1)
			{
				if (m_msh->GetCoordinateFlag() == 11)
				{
					int ibuff = idy;
					idy = idx;
					idx = ibuff;
				}
				if (m_msh->GetCoordinateFlag() == 12)
				{
					int ibuff = idz;
					idz = idx;
					idx = ibuff;
				}
			}
		}
		MeshLib::CElem* theEle = m_msh->ele_vector[A->elementIndex];
		int eleIndex = A->elementIndex;
		// Set the pointer that leads to the nodes of element

		// If Element is outside of the domain or in the sink or on the edge or at the corner vertix
		// Do not do anything. If not, then proceed.
		// Check if the particle is on the edge or node.
		int IsOnTheEdge = IsParticleOnTheEdge(A);
		int nnode = theEle->GetEdgesNumber();
		if (eleIndex != -10 && IsOnTheEdge == 0)
		{
			// Let's solve pore velocity.
			// It is simple because Sw stuff automatically handles in Richards Flow.
			// Thus, I only divide Darcy velocity by porosity only to get pore velocity.

			/*
			CMediumProperties *MediaProp = mmp_vector[theEle->GetPatchIndex()];
			double porosity = 0.0;
			if(MediaProp->porosity > 10-6)
			    porosity = MediaProp->porosity;       // This is for simple one.
			else
			    // This will get you porosity.
			    porosity = MediaProp->porosity_model_values[0];
			*/

			// I guess for Dual Porocity stuff,
			// this code should be revisited.

			// Get the number of nodes

			// Mount the edges of the element
			vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
			theEle->GetEdges(theEdgesOfThisElement);

			// Get physical coordinates of four corner points
			/*
			double x[8], y[8], z[8];
			double vx[8], vy[8], vz[8];
			for(int i=0; i<nnode; ++i)
			{
			    MeshLib::CNode* theNode = NULL;
			    theNode = theEle->GetNode(i);
			    x[i] = theNode->X();
			    y[i] = theNode->Y();
			    z[i] = theNode->Z();

			    vx[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idx)/porosity; //FM_TEST
			    vy[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idy)/porosity;
			    vz[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idz)/porosity;
			}
			*/

			// Mount the nodes of the edge
			vec<MeshLib::CNode*> theNodesOfThisEdge(3);

			// MiniFEM for 2D elements starts here
			// Find all the nodes for miniFEM first
			if (nnode == 4) // Quad
			{
				MeshLib::CNode* theNode = NULL;
				// Node 0 for miniFEM for quad
				theNode = theEle->GetNode(0);
				m_mini->nod_vector.push_back(theNode);
				// Node 1 For edge 0
				GetNodeOfMiniFEMforTheEdge(theNode, theEdgesOfThisElement[0], A);
				m_mini->nod_vector.push_back(theNode);
				// Node 2 for miniFEM for quad
				theNode = theEle->GetNode(1);
				m_mini->nod_vector.push_back(theNode);
				// Node 3 For edge 1
				GetNodeOfMiniFEMforTheEdge(theNode, theEdgesOfThisElement[1], A);
				m_mini->nod_vector.push_back(theNode);
				// Node 4 for miniFEM for quad
				theNode = theEle->GetNode(2);
				m_mini->nod_vector.push_back(theNode);
				// Node 5 For edge 2
				GetNodeOfMiniFEMforTheEdge(theNode, theEdgesOfThisElement[2], A);
				m_mini->nod_vector.push_back(theNode);
				// Node 6 for miniFEM for quad
				theNode = theEle->GetNode(3);
				m_mini->nod_vector.push_back(theNode);
				// Node 7 For edge 3
				GetNodeOfMiniFEMforTheEdge(theNode, theEdgesOfThisElement[3], A);
				m_mini->nod_vector.push_back(theNode);
				// Node 8 Finally this particle
				theNode->SetX(A->x);
				theNode->SetY(A->y);
				theNode->SetZ(A->z);
				m_mini->nod_vector.push_back(theNode);

				// Now elements for miniFEM
				// ele 0 for miniFEM
				miniEle->SetNodeIndex(0, m_mini->nod_vector[0]->GetIndex());
				miniEle->SetNodeIndex(1, m_mini->nod_vector[1]->GetIndex());
				miniEle->SetNodeIndex(2, m_mini->nod_vector[8]->GetIndex());
				miniEle->SetNodeIndex(3, m_mini->nod_vector[7]->GetIndex());
				m_mini->ele_vector.push_back(miniEle);
				// ele 1 for miniFEM
				miniEle->SetNodeIndex(0, m_mini->nod_vector[1]->GetIndex());
				miniEle->SetNodeIndex(1, m_mini->nod_vector[2]->GetIndex());
				miniEle->SetNodeIndex(2, m_mini->nod_vector[3]->GetIndex());
				miniEle->SetNodeIndex(3, m_mini->nod_vector[8]->GetIndex());
				m_mini->ele_vector.push_back(miniEle);
				// ele 2 for miniFEM
				miniEle->SetNodeIndex(0, m_mini->nod_vector[8]->GetIndex());
				miniEle->SetNodeIndex(1, m_mini->nod_vector[3]->GetIndex());
				miniEle->SetNodeIndex(2, m_mini->nod_vector[4]->GetIndex());
				miniEle->SetNodeIndex(3, m_mini->nod_vector[5]->GetIndex());
				m_mini->ele_vector.push_back(miniEle);
				// ele 3 for miniFEM
				miniEle->SetNodeIndex(0, m_mini->nod_vector[7]->GetIndex());
				miniEle->SetNodeIndex(1, m_mini->nod_vector[8]->GetIndex());
				miniEle->SetNodeIndex(2, m_mini->nod_vector[5]->GetIndex());
				miniEle->SetNodeIndex(3, m_mini->nod_vector[6]->GetIndex());
				m_mini->ele_vector.push_back(miniEle);

				// Assemble
				fem = new CFiniteElementStd(m_pcs, m_msh->GetCoordinateFlag());
				// I am going to create global matrix here
				MeshLib::CElem* elem = NULL;
				for (int d = 0; d < 2; ++d)
				{
/* Initializations */
/* System matrix */
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// Todo
#elif defined(NEW_EQS) // WW
					m_pcs->EQSInitialize();
#else
					SetZeroLinearSolver(m_pcs->getEQSPointer());
#endif

					for (int i = 0; i < (long)m_msh->ele_vector.size(); i++)
					{
						elem = m_mini->ele_vector[i];
						fem->ConfigElement(elem, m_pcs->m_num->ele_gauss_points);
						// Assembly gotta be written different way
						fem->Assembly(0, d);
					}

					m_pcs->IncorporateBoundaryConditions(-1, d);

// Solve for velocity
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// Todo
#elif defined(NEW_EQS)

#if defined(LIS)
					double* x;
					int size = m_msh->nod_vector.size();
					x = new double[size];
					m_pcs->EQSSolver(x); // an option added to tell FLUID_MOMENTUM for sparse matrix system.
					cout << "Solver passed in FLUID_MOMENTUM."
					     << "\n";
#endif
#else
					m_pcs->ExecuteLinearSolver(m_pcs->getEQSPointer());
#endif
				}
			}
			else if (nnode == 3)
			{
			}
			// else
			//	;	// This shouldn't happen here. There are only tri or quad ele's in 2D
		}
		else if (IsOnTheEdge != 0) // The particle is on the edge.
		{
			// We have update edge index for this particle performing
			// IsParticleOnTheEdge() function
			// Thus, we only need to interpolate particle velocity along this edge.
			// Even for this, I will use 1d Galerkin method.
		}
		// else
		//	;	// Think later for out of boundary particles.
	}
	else // Real Space and reference space way
	{
		m_pcs = PCSGet("FLUID_MOMENTUM");
		// FM_TEST
		int idx, idy, idz;
		if (m_pcs)
		{
			idx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
			idy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
			idz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;
		}
		else if (flow_pcs)
		{
			idx = flow_pcs->GetNodeValueIndex("VELOCITY_X1");
			idy = flow_pcs->GetNodeValueIndex("VELOCITY_Y1");
			idz = flow_pcs->GetNodeValueIndex("VELOCITY_Z1");
			m_pcs = flow_pcs;
			if (m_msh->GetCoordinateFlag() / 10 == 1)
			{
				if (m_msh->GetCoordinateFlag() == 11)
				{
					int ibuff = idy;
					idy = idx;
					idx = ibuff;
				}
				if (m_msh->GetCoordinateFlag() == 12)
				{
					int ibuff = idz;
					idz = idx;
					idx = ibuff;
				}
			}
		}
		MeshLib::CElem* theEle = m_msh->ele_vector[A->elementIndex];
		int eleIndex = A->elementIndex;
		// Set the pointer that leads to the nodes of element
		// OK411 CNode* node = NULL;

		// Element is outside of the domain or in the sink
		// Do not do anything. If not, then proceed.
		if (eleIndex != -10)
		{
			int nnode = theEle->GetEdgesNumber();

			// Let's solve pore velocity.
			// It is simple because Sw stuff automatically handles in Richards Flow.
			// Thus, I only divide Darcy velocity by porosity only to get pore velocity.

			/*
			CMediumProperties *MediaProp = mmp_vector[theEle->GetPatchIndex()];
			double porosity = 0.0;
			if(MediaProp->porosity > 10-6)
			    porosity = MediaProp->porosity;       // This is for simple one.
			else
			    // This will get you porosity.
			    porosity = MediaProp->porosity_model_values[0];
			*/

			// I guess for Dual Porocity stuff,
			// this code should be revisited.

			// Get the number of nodes
			// OK411 int nnodes = theEle->GetVertexNumber();

			// Mount the edges of the element
			vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
			theEle->GetEdges(theEdgesOfThisElement);

			// Get physical coordinates of four corner points
			double x[4], y[4] /*, z[4]*/;
			// double vx[4], vy[4], vz[4];
			for (int i = 0; i < nnode; ++i)
			{
				double const* const pnt(theEle->GetNode(i)->getData());
				x[i] = pnt[0];
				y[i] = pnt[1];
				/*
				z[i] = pnt[2];

				vx[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idx)/porosity; //FM_TEST
				vy[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idy)/porosity;
				vz[i] = m_pcs->GetNodeValue(theEle->GetNodeIndex(i),idz)/porosity;
				*/
			}
			// solve for Jm at xm = (xhat,yhat)=(1/2,1/2) <- RT0
			double Jm;
			Jm = 0.5 * (x[0] * y[1] - y[0] * x[1] + x[1] * y[2] - y[1] * x[2] + y[0] * x[3] - x[0] * y[3] + x[2] * y[3]
			            - y[2] * x[3]);

			// Mount the nodes of the edge
			vec<MeshLib::CNode*> theNodesOfThisEdge(3);

			// Solve for proper index in reference space
			double const* Ecenter = theEle->GetGravityCenter();
			double u[4]; // Edge flux in reference space
			for (int i = 0; i < 4; ++i)
				u[i] = 0.0;

			for (int i = 0; i < 4; ++i)
			{
				// Get the nodes of the edge i
				theEdgesOfThisElement[i]->GetNodes(theNodesOfThisEdge);

				double node1[3] = {theNodesOfThisEdge[0]->getData()[0], theNodesOfThisEdge[0]->getData()[1], 0.0};
				double node2[3] = {theNodesOfThisEdge[1]->getData()[0], theNodesOfThisEdge[1]->getData()[1], 0.0};
				// TF           node1[0] = theNodesOfThisEdge[0]->X(); node1[1] = theNodesOfThisEdge[0]->Y();
				// TF           node2[0] = theNodesOfThisEdge[1]->X(); node2[1] = theNodesOfThisEdge[1]->Y();
				// Get the referece position of these two ending points of the edge
				IsoparametricMappingQuadfromPtoR(A->elementIndex, node1);
				IsoparametricMappingQuadfromPtoR(A->elementIndex, node2);

				double dx, dy;
				dx = theEdgesOfThisElement[i]->GetVelocity(0);
				dy = theEdgesOfThisElement[i]->GetVelocity(1);

				double edge2mid[3], EdgeMidPoint[3];
				theEdgesOfThisElement[i]->GetEdgeMidPoint(EdgeMidPoint);
				for (int j = 0; j < 3; ++j)
					edge2mid[j] = Ecenter[j] - EdgeMidPoint[j];

				double angle
				    = acos((dx * edge2mid[0] + dy * edge2mid[1])
				           / (sqrt(dx * dx + dy * dy) * sqrt(edge2mid[0] * edge2mid[0] + edge2mid[1] * edge2mid[1])));

				// EA and EB: x1hat = x2hat and x1hat*x2hat > 0
				if (fabs(node1[0] - node2[0]) < 0.5 && (node1[0] * node2[0]) > 0.0)
				{
					// EA: x1hat < 0
					if (node1[0] < 0.0)
					{
						if (angle > 3.141592 / 2.0)
							u[0] = -sqrt(dx * dx + dy * dy) / Jm;
						else
							u[0] = sqrt(dx * dx + dy * dy) / Jm;
					}
					else // EB
					{
						if (angle > 3.141592 / 2.0)
							u[1] = sqrt(dx * dx + dy * dy) / Jm;
						else
							u[1] = -sqrt(dx * dx + dy * dy) / Jm;
					}
				}
				// Else then, EC and ED
				else
				{
					// EC: y1hat < 0
					if (node1[1] < 0.0)
					{
						if (angle > 3.141592 / 2.0)
							u[2] = -sqrt(dx * dx + dy * dy) / Jm;
						else
							u[2] = sqrt(dx * dx + dy * dy) / Jm;
					}
					else // ED
					{
						if (angle > 3.141592 / 2.0)
							u[3] = sqrt(dx * dx + dy * dy) / Jm;
						else
							u[3] = -sqrt(dx * dx + dy * dy) / Jm;
					}
				}
			}

			double E1[3], E2[3], E3[3], E4[3]; // mid points of each edge
			theEdgesOfThisElement[0]->GetEdgeMidPoint(E1);
			theEdgesOfThisElement[1]->GetEdgeMidPoint(E2);
			theEdgesOfThisElement[2]->GetEdgeMidPoint(E3);
			theEdgesOfThisElement[3]->GetEdgeMidPoint(E4);
			IsoparametricMappingQuadfromPtoR(A->elementIndex, E1);
			IsoparametricMappingQuadfromPtoR(A->elementIndex, E2);
			IsoparametricMappingQuadfromPtoR(A->elementIndex, E3);
			IsoparametricMappingQuadfromPtoR(A->elementIndex, E4);

			// ux = a+b(x-x0); uy = c + d(y-y0)
			double a, b, c, d, x0, y0;

			x0 = E4[0];
			y0 = E1[1];
			a = u[0];
			b = (u[1] - u[0]) / 2.0; // set dx for reference space unit
			c = u[2];
			d = (u[3] - u[2]) / 2.0;

			double R[3];
			R[0] = A->x;
			R[1] = A->y;
			R[2] = 0.0;
			IsoparametricMappingQuadfromPtoR(A->elementIndex, R);
			A->Vx = (a + b * (R[0] - x0));
			A->Vy = (c + d * (R[1] - y0));
			A->Vz = 0.0;
		}
		// else
		//	;	// Think later for out of boundary particles.
	}

	delete miniEle;
	delete m_mini;
}

void RandomWalk::GetNodeOfMiniFEMforTheEdge(MeshLib::CNode* theNode, MeshLib::CEdge* theEdge, Particle* A)
{
	// Just for two nodes of the edge
	double x[2] = {theEdge->GetNode(0)->getData()[0], theEdge->GetNode(0)->getData()[1]};
	double y[2] = {theEdge->GetNode(1)->getData()[0], theEdge->GetNode(1)->getData()[1]};

	// Two rules. The angle is perpendicular. So the dot product is zero
	// The particle is on this edge. Line equation
	// a = x1-x2, b=y1-y2, aa=ax0+by0, bb=bx1-ay1
	double a, b, aa, bb;
	a = x[0] - x[1];
	b = y[0] - y[1];
	aa = a * A->x + b * A->y;
	bb = b * x[0] - a * y[0];
	double det = -a * a - b * b;
	double X, Y;
	X = -(a * aa + b * bb) / det;
	Y = (-b * aa + a * bb) / det;

	theNode->SetX(X);
	theNode->SetY(Y);
	theNode->SetZ(0.0);
}

/**************************************************************************
   Class: RandomWalk
   Task: The function returns 0 if particle is not on any of edge of the element
     1 if particle is on any of edge of the element
   Programing:
   02/2007 PCH Implementation
   last modification:
**************************************************************************/
int RandomWalk::IsParticleOnTheEdge(Particle* A)
{
	double tolerance = 1e-6;
	MeshLib::CElem* theEle = m_msh->ele_vector[A->elementIndex];
	int nnode = theEle->GetVertexNumber();

	// If element is triangle or quadrilateral.
	if (nnode == 3 || nnode == 4)
	{
		// Mount the edges of the element
		vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
		theEle->GetEdges(theEdgesOfThisElement);

		//      double x1[3], x2[3];                        //OK411 , v[3];
		vec<MeshLib::CNode*> theNodesOfThisEdge(3);

		for (int i = 0; i < nnode; ++i)
		{
			theEdgesOfThisElement[i]->GetNodes(theNodesOfThisEdge);
			double const* const x1(theNodesOfThisEdge[0]->getData());
			//         x1[0]=pnt1[0];
			//         x1[1]=pnt1[1];
			//         x1[2]=pnt1[2];
			double const* const x2(theNodesOfThisEdge[1]->getData());
			//         x2[0]=pnt2[0];
			//         x2[1]=pnt2[1];
			//         x2[2]=pnt2[2];

			double x2x1[3], x1x0[3];
			x2x1[0] = x2[0] - x1[0];
			x2x1[1] = x2[1] - x1[1];
			x2x1[2] = x2[2] - x1[2];
			x1x0[0] = x1[0] - A->x;
			x1x0[1] = x1[1] - A->y;
			x1x0[2] = x1[2] - A->z;
			double x2x1square;
			x2x1square = (x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1])
			             + (x2[2] - x1[2]) * (x2[2] - x1[2]);

			double cvec[3];
			CrossProduction(x2x1, x1x0, cvec);

			double distance;
			distance = sqrt((cvec[0] * cvec[0] + cvec[1] * cvec[1] + cvec[2] * cvec[2]) / x2x1square);

			if (distance < tolerance) // The particle is on this edge
			{
				// Update the edge that the particle is on now.
				A->edgeIndex = theEdgesOfThisElement[i]->GetIndex();

				return 1; // Yes, it is on one of the edges in the element
			}
		}
	}
	else
		return 0; // For 3D element, later....

	return 0; // No, The particle is not any edge in the element.
}

/**************************************************************************
   Class: RandomWalk
   Task: This function interpolates location of the particle for a given dt
     based on the bilinear method
     The function requires FDM-like grid.
   Programing:
   02/2007 PCH Implementation
   last modification:
**************************************************************************/
double* RandomWalk::InterpolateLocationOfTheParticleByBilinear(Particle* A, double dt)
{
	double* x = new double[3];

	// Get the element that the particle belongs
	m_msh = selectMeshForFluidMomentumProcess();
	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	int eleIndex = IndexOfTheElementThatThisParticleBelong(0, A);
	if (eleIndex == -5)
		// This is a temperary measure.
		// Gotta be written better later on.
		eleIndex = 0;

	// Element is outside of the domain or in the sink
	// Do not do anything. If not, then proceed.
	if (eleIndex != -10)
	{
		A->elementIndex = eleIndex;
		MeshLib::CElem* m_ele = m_msh->ele_vector[eleIndex];

		int nnode = m_ele->GetEdgesNumber();

		// Let's get the hydraulic conductivity first.
		CMediumProperties* MediaProp = mmp_vector[m_ele->GetPatchIndex()];
		// OK411 int phase = 0;
		CFluidProperties* FluidProp = mfp_vector[0];
		double* kTensor = MediaProp->PermeabilityTensor(eleIndex);
		double k = kTensor[0];

		A->K = k * FluidProp->Density() * 9.81 / FluidProp->Viscosity();

		// Get the number of nodes
		// OK411 int nnodes = m_ele->GetVertexNumber();

		// Mount the edges of the element
		vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
		m_ele->GetEdges(theEdgesOfThisElement);

		double E1[3], E2[3], E3[3], E4[3]; // mid points of each edge
		theEdgesOfThisElement[0]->GetEdgeMidPoint(E1);
		theEdgesOfThisElement[1]->GetEdgeMidPoint(E2);
		theEdgesOfThisElement[2]->GetEdgeMidPoint(E3);
		theEdgesOfThisElement[3]->GetEdgeMidPoint(E4);

		// ux = a+b(x-x0); uy = c + d(y-y0)
		double a, b, c, d, x0, y0;
		x0 = E4[0];
		y0 = E1[1];
		a = theEdgesOfThisElement[3]->GetVelocity(0);
		b = (theEdgesOfThisElement[1]->GetVelocity(0) - theEdgesOfThisElement[3]->GetVelocity(0)) / dx;
		c = theEdgesOfThisElement[0]->GetVelocity(1);
		d = (theEdgesOfThisElement[2]->GetVelocity(1) - theEdgesOfThisElement[0]->GetVelocity(1)) / dy;

		A->x = (a + b * (A->x - x0)) / b * exp(b * dt) - a / b;
		A->y = (c + d * (A->y - y0)) / d * exp(d * dt) - c / d;
		A->z = 0.0;
	}
	// else
	//	;	// Think later for out of boundary particles.

	return x;
}

/**************************************************************************
   Class: RandomWalk
   Programing: Locate element location of particle, for FDM method
   PCH Implementation
   03/2010 JT modified cellular search algorithm
**************************************************************************/
int RandomWalk::IndexOfTheElementThatThisParticleBelong(int option, Particle* A)
{
	int index = -10;

	if (option == 0)
	{
		long i, j, k, iFDM;
		double x, y, z;
		x = A->x;
		y = A->y;
		z = A->z;
		i = j = k = 0;

		// Set off the domain first
		if (xrw_range > 1.e-12) // only if non-negligible range in this direction
		{
			if (x > pnt_x_max || x < pnt_x_min)
				return index;
			i = (long int)floor((x - pnt_x_min) / dx);
		}
		if (yrw_range > 1.e-12)
		{
			if (y > pnt_y_max || y < pnt_y_min)
				return index;
			j = (long int)floor((y - pnt_y_min) / dy);
		}
		if (zrw_range > 1.e-12)
		{
			if (z > pnt_z_max || x < pnt_z_min)
				return index;
			k = (long int)floor((z - pnt_z_min) / dz);
		}

		iFDM = k * (nx * ny) + j * nx + i;
		index = indexFDM[iFDM].eleIndex;
		/*
		   if(index == 322 || index == 323 || index == 342 || index == 343)	// Sink condition
		   return -10;
		   else
		 */
		return index;
	}
	else
		return GetTheElementOfTheParticleFromNeighbor(A);
}

/**************************************************************************
   Class: RandomWalk
   Task: The function solves two intersections along x or y or z axis.
     2: The function returns two intersections
     1: The function returns one intersection
    -1: The function failed
   axis = 0: a line parallel to the x axis
   axis = 1: a line parallel to the y axis
   axis = 2: a line parallel to the z axis
   Programing:
   11/2005 PCH Implementation
   02/2006 PCH Improvement for RWPT in Fracture networks.
   last modification:
**************************************************************************/
int RandomWalk::SolveForTwoIntersectionsInTheElement(Particle* A, double* P1, double* P2, int axis)
{
	// Get the element that the particle belongs
	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];
	// Set the pointer that leads to the nodes of element
	MeshLib::CNode* node = NULL;

	// Get the number of nodes
	int nnodes = m_ele->GetVertexNumber();
	// Allocate the memory accordingly
	Particle* vertex = NULL;
	vertex = new Particle[nnodes]();
	int R = 0, L = 0;

	// Set the size of displacement
	double disp = 1e4; // This should be bigger the largest element size.

	// RWPT-IM
	// Get the cooridinate of the nodes in the element
	for (int i = 0; i < nnodes; ++i)
	{
		node = m_ele->GetNode(i);
		double X[3] = {node->getData()[0], node->getData()[1], node->getData()[2]};
		ToTheXYPlane(m_ele, X);
		vertex[i].x = X[0];
		vertex[i].y = X[1];
		vertex[i].z = X[2];
	}

	// Solve for the line equation
	for (int i = 0; i < nnodes; ++i)
	{
		double p1[3], p2[3], p3[3], p4[3];
		// Need coordinate transform here.
		p1[0] = vertex[i % nnodes].x;
		p1[1] = vertex[i % nnodes].y;
		p1[2] = vertex[i % nnodes].z;
		p2[0] = vertex[(i + 1) % nnodes].x;
		p2[1] = vertex[(i + 1) % nnodes].y;
		p2[2] = vertex[(i + 1) % nnodes].z;
		// RWPT-IM
		double X[3];
		X[0] = A->x;
		X[1] = A->y;
		X[2] = A->z;
		ToTheXYPlane(m_ele, X);
		for (int p = 0; p < 3; ++p)
			p3[p] = p4[p] = X[p];

		for (int j = 0; j < 2; ++j)
		{
			// See if there is an intersection in this line.
			// if a line is set to be parallel to x axis on the right
			if (axis == 0 && j == 0)
				p4[0] = X[0] + disp;
			// if a line is set to be parallel to y axis on the right,
			else if (axis == 1 && j == 0)
				p4[1] = X[1] + disp;
			// if a line is set to be parallel to z axis on the right,
			else if (axis == 2 && j == 0)
				p4[2] = X[2] + disp;
			// if a line is set to be parallel to x axis on the left
			else if (axis == 0 && j == 1)
				p4[0] = X[0] - disp;
			// if a line is set to be parallel to y axis on the left,
			else if (axis == 1 && j == 1)
				p4[1] = X[1] - disp;
			// if a line is set to be parallel to z axis on the left,
			else if (axis == 2 && j == 1)
				p4[2] = X[2] - disp;
			else
			{
				printf("Axis type in searching the intersection failed. Wrong axis type.\n");
				abort();
			}

			double x = 0.0, y = 0.0, ra = 0.0, rb = 0.0;

			int status
			    = G_intersect_line_segments(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], &ra, &rb, &x, &y);
			// RWPT-IM P1 and P2 are already on the XY plane.
			if (status == 1 && j == 0)
			{
				P1[0] = x;
				P1[1] = y;
				P1[2] = 0.0;
				// Transform back the coordinates.
				ToTheRealPlane(m_ele, P1);

				++R;
			}
			else if (status == 1 && j == 1)
			{
				P2[0] = x;
				P2[1] = y;
				P2[2] = 0.0;
				// Transform back the coordinates.
				ToTheRealPlane(m_ele, P2);

				++L;
			}
		}
	}

	// Free the memory for this function
	delete[] vertex;

	if (R + L == 2)
		return 2; // The job succeeded
	else if (R + L == 1)
		return 1;
	else
		return -1; // The job failed.
}

/**************************************************************************
   Class: RandomWalk
   Task: The function solves three displacement by derivatives of
     dispersion tensor.
     1: The function succeeded
    -1: The function failed
   Programing:
   11/2005 PCH Implementation
   02/2006 PCH Improved for the RWPT method in Fracture Networks.
   last modification:
**************************************************************************/
int RandomWalk::SolveForDisplacementByDerivativeOfDispersion(Particle* A, double* dD)
{
	double TensorOfdD[9];

	// Solve for the derivative of velocity first
	// statusForDeivativeOfVelocity is never further used down the code.
	// WW int statusForDeivativeOfVelocity = -10;
	// WW statusForDeivativeOfVelocity = SolveForDerivativeOfVelocity(A);

	// Solve for the tensor of dispersion derivatives
	// Extract the dispersivities from the group that the particle belongs
	// To extract dispersivities from material properties
	// This should be checked if the dispersivity gets correctly.
	CMediumProperties* m_mat_mp = NULL;
	double alphaL = 0.0, alphaT = 0.0;
	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];
	int group = m_ele->GetPatchIndex();
	m_mat_mp = mmp_vector[group];
	alphaL = m_mat_mp->mass_dispersion_longitudinal;
	alphaT = m_mat_mp->mass_dispersion_transverse;

	// RWPT - IM
	// This thing should be done on the XY plane too.
	double V[3];
	V[0] = A->Vx;
	V[1] = A->Vy;
	V[2] = A->Vz;
	ToTheXYPlane(A->elementIndex, V);

	double U = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);

	TensorOfdD[0] = V[0] * A->dVxdx * (alphaL * (2.0 / U - V[0] * V[0] / (U * U * U))
	                                   - alphaT * (V[1] * V[1] + V[2] * V[2]) / (U * U * U));
	TensorOfdD[1] = (alphaL - alphaT) * (A->dVydy * V[0] / U - V[0] * V[1] * V[1] / (U * U * U) * A->dVydy);
	TensorOfdD[2] = (alphaL - alphaT) * (A->dVzdz * V[0] / U - V[0] * V[2] * V[2] / (U * U * U) * A->dVzdz);
	TensorOfdD[3] = (alphaL - alphaT) * (A->dVxdx * V[1] / U - V[1] * V[0] * V[0] / (U * U * U) * A->dVxdx);
	TensorOfdD[4] = V[1] * A->dVydy * (alphaL * (2.0 / U - V[1] * V[1] / (U * U * U))
	                                   - alphaT * (V[0] * V[0] + V[2] * V[2]) / (U * U * U));
	TensorOfdD[5] = (alphaL - alphaT) * (A->dVzdz * V[1] / U - V[1] * V[2] * V[2] / (U * U * U) * A->dVzdz);
	TensorOfdD[6] = (alphaL - alphaT) * (A->dVxdx * V[2] / U - V[2] * V[0] * V[0] / (U * U * U) * A->dVxdx);
	TensorOfdD[7] = (alphaL - alphaT) * (A->dVydy * V[2] / U - V[2] * V[1] * V[1] / (U * U * U) * A->dVydy);
	TensorOfdD[8] = V[2] * A->dVzdz * (alphaL * (2.0 / U - V[2] * V[2] / (U * U * U))
	                                   - alphaT * (V[0] * V[0] + V[1] * V[1]) / (U * U * U));

	// Solve the three displacement by the tensor of dispersion derivative.
	dD[0] = TensorOfdD[0] + TensorOfdD[1] + TensorOfdD[2];
	dD[1] = TensorOfdD[3] + TensorOfdD[4] + TensorOfdD[5];
	dD[2] = TensorOfdD[6] + TensorOfdD[7] + TensorOfdD[8];

	return 1;
}

/**************************************************************************
   Class: RandomWalk
   Task: The function solves three main derivative of velocity. The rest of
     the components is assumed to be zero.
     1: The function succeeded
    -1: The function failed
   Programing:
   11/2005 PCH Implementation
   02/2006 PCH Improvement for fracture networks.
   last modification:
**************************************************************************/
int RandomWalk::SolveForDerivativeOfVelocity(Particle* A)
{
	int status = -10; // Set to be meaningliss in the beginning
	m_msh = fem_msh_vector[0];
	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];

	// If not 1D,
	if (m_ele->GetDimension() == 3)
	{
		/// 1. Get coordinates of the particles, x, y, z
		/// 2. Find the element where the particle is.
		/// 3. Compute Gradient of velocity

		/// 4. Convert real coordinates to local coordinates
		/// Gradient of shape  function: CElement::dshapefct
		/// Gradient of velocity:
		/// =V_i*dshapefct[i];
		return -1; // not supported yet
	}
	else if (m_ele->GetDimension() == 2)
	{
		// intersections for x and y axis
		double x1[3], x2[3], y1[3], y2[3]; // I don't put the intersections for z direction for now.

		// RWPT-IM x1 and x2 are the intersection coordinates on the XY plane.
		// But the position of Particle A is on the realy plane.
		// Get the two intersecitions parallel to x axis
		status = SolveForTwoIntersectionsInTheElement(A, x1, x2, 0);
		// RWPT-IM After SolveForTwoIntersectionsInTheElement,
		// All the coordinates are on the real plane.
		// Check if the function succeeded.
		if (status == -1)
			//		printf("Solving two intersections parallel to x axis failed\n");
			return -1; // Failed
		// Solve for the velocity for two intersections
		Particle XR, XL;
		// RWPT-IM
		XR = XL = *A;
		// Again, the real plane coordinates.
		XR.x = x1[0];
		XR.y = x1[1];
		XR.z = x1[2];
		XL.x = x2[0];
		XL.y = x2[1];
		XL.z = x2[2];

		// Interpolating velocity by the real coordinates should be no problem.
		if (PURERWPT != 2)
		{
			InterpolateVelocityOfTheParticleByInverseDistance(&XR);
			InterpolateVelocityOfTheParticleByInverseDistance(&XL);
		}
		else
		{
			InterpolateVelocityOfTheParticleByBilinear(GridOption, &XR);
			InterpolateVelocityOfTheParticleByBilinear(GridOption, &XL);
		}

		// Solve for dVxdx
		double x = XR.x - XL.x;
		double y = XR.y - XL.y;
		double z = XR.z - XL.z;
		double dx = sqrt(x * x + y * y + z * z); // The distance does not make any difference.
		// RWPT-IM
		// Let me think if velocity should projected to the connected plane or treated in true 3D.
		// Yes. Velocity should be on the XY plane
		double Vx[3];
		Vx[0] = XR.Vx - XL.Vx;
		Vx[1] = XR.Vy - XL.Vy;
		Vx[2] = XR.Vz - XL.Vz;
		ToTheXYPlane(A->elementIndex, Vx);
		A->dVxdx = Vx[0] / dx; // A->dVxdx = (XR.Vx - XL.Vx) / dx;

		// RWPT-IM Just the same thing one more time.
		// Get the two intersecitions parallel to y axis
		status = SolveForTwoIntersectionsInTheElement(A, y1, y2, 1);
		if (status == -1)
			//		printf("Solving two intersections parallel to y axis failed\n");
			return -1; // Failed

		// Solve for the velocity for two intersections
		Particle YR, YL;
		YR = YL = *A;
		YR.x = y1[0];
		YR.y = y1[1];
		YR.z = y1[2];
		YL.x = y2[0];
		YL.y = y2[1];
		YL.z = y2[2];
		if (PURERWPT != 2)
		{
			InterpolateVelocityOfTheParticleByInverseDistance(&YR);
			InterpolateVelocityOfTheParticleByInverseDistance(&YL);
		}
		else
		{
			InterpolateVelocityOfTheParticleByBilinear(GridOption, &XR);
			InterpolateVelocityOfTheParticleByBilinear(GridOption, &XL);
		}

		// Solve for dVydy
		x = YR.x - YL.x;
		y = YR.y - YL.y;
		z = YR.z - YL.z;
		double dy = sqrt(x * x + y * y + z * z);
		double Vy[3];
		Vy[0] = YR.Vx - YL.Vx;
		Vy[1] = YR.Vy - YL.Vy;
		Vy[2] = YR.Vz - YL.Vz;
		ToTheXYPlane(A->elementIndex, Vy);
		A->dVydy = Vy[1] / dy; // A->dVydy = (YR.Vy - YL.Vy) / dy;

		// Just set dVzdz to be zero for now
		A->dVzdz = 0.0;

		// Return 1 for success
		return 1;
	}
	else // 1D line element
	{
		// Solve the length of the element
		double length = m_ele->GetVolume();

		m_pcs = PCSGet("FLUID_MOMENTUM");

		// FM_TEST
		int idx = -1, idy = -1, idz = -1;
		if (m_pcs)
		{
			idx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
			idy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
			idz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;
		}
		else if (flow_pcs)
		{
			idx = flow_pcs->GetNodeValueIndex("VELOCITY_X1");
			idy = flow_pcs->GetNodeValueIndex("VELOCITY_Y1");
			idz = flow_pcs->GetNodeValueIndex("VELOCITY_Z1");
			m_pcs = flow_pcs;
			if (m_msh->GetCoordinateFlag() / 10 == 1)
			{
				if (m_msh->GetCoordinateFlag() == 11)
				{
					int ibuff = idy;
					idy = idx;
					idx = ibuff;
				}
				if (m_msh->GetCoordinateFlag() == 12)
				{
					int ibuff = idz;
					idz = idx;
					idx = ibuff;
				}
			}
		}

		double v1[3], v2[3];
		v1[0] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idx); // FM_TEST
		v1[1] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idy);
		v1[2] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idz);
		v2[0] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idx);
		v2[1] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idy);
		v2[2] = m_pcs->GetNodeValue(m_ele->GetNodeIndex(0), idz);

		int coordinateflag = m_msh->GetCoordinateFlag();
		if (coordinateflag == 10) // x only
		{
			A->dVxdx = (v2[0] - v1[0]) / length;
			A->dVydy = A->dVzdz = 0.0;
		}
		else if (coordinateflag == 11) // y only
		{
			A->dVydy = (v2[1] - v1[1]) / length;
			A->dVxdx = A->dVzdz = 0.0;
		}
		else if (coordinateflag == 12) // z only
		{
			A->dVzdz = (v2[2] - v1[2]) / length;
			A->dVydy = A->dVxdx = 0.0;
		}
		else // Something Wrong.
			abort();

		return 1;
	}
}

void RandomWalk::CopyParticleCoordToArray(Particle* A, double* x1buff, double* x2buff, double* x3buff, double* x4buff)
{
	x1buff[0] = A[0].x;
	x1buff[1] = A[0].y;
	x1buff[2] = A[0].z;
	x2buff[0] = A[1].x;
	x2buff[1] = A[1].y;
	x2buff[2] = A[1].z;
	x3buff[0] = A[2].x;
	x3buff[1] = A[2].y;
	x3buff[2] = A[2].z;
	x4buff[0] = A[3].x;
	x4buff[1] = A[3].y;
	x4buff[2] = A[3].z;
}

/**************************************************************************
   MSHLib-Method:
   Task:Compute the volume of the object
   Programing:
   09/2005 PCH Implementation
**************************************************************************/
double RandomWalk::ComputeVolume(Particle* A, MeshLib::CElem* m_ele)
{
	//   double x1buff[3];
	//   double x2buff[3];
	//   double x3buff[3];
	//   double x4buff[3];
	double volume = 0.0;
	double* PieceOfVolume = NULL;

	double A2buff[3];

	A2buff[0] = A->x;
	A2buff[1] = A->y;
	A2buff[2] = A->z;

	// If this is not a line element, get three verteces.
	//   if(m_ele->GetElementType()!=MshElemType::LINE)
	//   {
	//      node = m_ele->GetNode(0);
	//      x1buff[0] = node->X();
	//      x1buff[1] = node->Y();
	//      x1buff[2] = node->Z();
	//
	//      node = m_ele->GetNode(1);
	//      x2buff[0] = node->X();
	//      x2buff[1] = node->Y();
	//      x2buff[2] = node->Z();
	//
	//      node = m_ele->GetNode(2);
	//      x3buff[0] = node->X();
	//      x3buff[1] = node->Y();
	//      x3buff[2] = node->Z();
	//   }

	// LINES = 1
	if (m_ele->GetElementType() == MshElemType::LINE)
	{
		PieceOfVolume = new double[2]();
		for (int i = 0; i < 2; ++i)
		{
			double const* const pnt(m_ele->GetNode(i)->getData());
			double x2buff[3] = {pnt[0] - A2buff[0], pnt[1] - A2buff[1], pnt[2] - A2buff[2]};
			PieceOfVolume[i] = sqrt(x2buff[0] * x2buff[0] + x2buff[1] * x2buff[1] + x2buff[2] * x2buff[2]);
			volume += PieceOfVolume[i];
		}
	}
	// RECTANGLES = 2
	if (m_ele->GetElementType() == MshElemType::QUAD)
	{
		PieceOfVolume = new double[4]();

		//      node = m_ele->GetNode(3);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();
		double const* const x1buff(m_ele->GetNode(0)->getData());
		double const* const x2buff(m_ele->GetNode(1)->getData());
		double const* const x3buff(m_ele->GetNode(2)->getData());
		double const* const x4buff(m_ele->GetNode(3)->getData());

		PieceOfVolume[0] = ComputeDetTri(x1buff, x2buff, A2buff);
		PieceOfVolume[1] = ComputeDetTri(x2buff, x3buff, A2buff);
		PieceOfVolume[2] = ComputeDetTri(x3buff, x4buff, A2buff);
		PieceOfVolume[3] = ComputeDetTri(x4buff, x1buff, A2buff);

		for (int i = 0; i < 4; ++i)
			volume += PieceOfVolume[i];
	}
	// HEXAHEDRA = 3
	if (m_ele->GetElementType() == MshElemType::HEXAHEDRON)
	{
		PieceOfVolume = new double[12]();

		// 2,1,4,3 face
		//      node = m_ele->GetNode(1);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(0);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(3);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff(m_ele->GetNode(1)->getData());
		double const* const x2buff(m_ele->GetNode(0)->getData());
		double const* const x3buff(m_ele->GetNode(3)->getData());
		double const* const x4buff(m_ele->GetNode(2)->getData());

		PieceOfVolume[0] = ComputeDetTex(A2buff, x1buff, x2buff, x4buff);
		PieceOfVolume[1] = ComputeDetTex(A2buff, x2buff, x3buff, x4buff);

		// 5,6,7,8 face
		//      node = m_ele->GetNode(4);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(6);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(7);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff1(m_ele->GetNode(4)->getData());
		double const* const x2buff1(m_ele->GetNode(5)->getData());
		double const* const x3buff1(m_ele->GetNode(6)->getData());
		double const* const x4buff1(m_ele->GetNode(7)->getData());

		PieceOfVolume[2] = ComputeDetTex(A2buff, x1buff1, x2buff1, x4buff1);
		PieceOfVolume[3] = ComputeDetTex(A2buff, x2buff1, x3buff1, x4buff1);

		// 1,5,8,4 face
		//      node = m_ele->GetNode(0);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(4);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(7);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(3);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff2(m_ele->GetNode(0)->getData());
		double const* const x2buff2(m_ele->GetNode(4)->getData());
		double const* const x3buff2(m_ele->GetNode(7)->getData());
		double const* const x4buff2(m_ele->GetNode(3)->getData());

		PieceOfVolume[4] = ComputeDetTex(A2buff, x1buff2, x2buff2, x4buff2);
		PieceOfVolume[5] = ComputeDetTex(A2buff, x2buff2, x3buff2, x4buff2);

		// 8,7,3,4 face
		//      node = m_ele->GetNode(7);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(6);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(3);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff3(m_ele->GetNode(7)->getData());
		double const* const x2buff3(m_ele->GetNode(6)->getData());
		double const* const x3buff3(m_ele->GetNode(2)->getData());
		double const* const x4buff3(m_ele->GetNode(3)->getData());

		PieceOfVolume[6] = ComputeDetTex(A2buff, x1buff3, x2buff3, x4buff3);
		PieceOfVolume[7] = ComputeDetTex(A2buff, x2buff3, x3buff3, x4buff3);

		// 2,3,7,6 face
		//      node = m_ele->GetNode(1);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(6);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff4(m_ele->GetNode(1)->getData());
		double const* const x2buff4(m_ele->GetNode(2)->getData());
		double const* const x3buff4(m_ele->GetNode(6)->getData());
		double const* const x4buff4(m_ele->GetNode(5)->getData());

		PieceOfVolume[8] = ComputeDetTex(A2buff, x1buff4, x2buff4, x4buff4);
		PieceOfVolume[9] = ComputeDetTex(A2buff, x2buff4, x3buff4, x4buff4);

		// 1,2,6,5 face
		//      node = m_ele->GetNode(0);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(1);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(4);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff5(m_ele->GetNode(0)->getData());
		double const* const x2buff5(m_ele->GetNode(1)->getData());
		double const* const x3buff5(m_ele->GetNode(5)->getData());
		double const* const x4buff5(m_ele->GetNode(4)->getData());

		PieceOfVolume[10] = ComputeDetTex(A2buff, x1buff5, x2buff5, x4buff5);
		PieceOfVolume[11] = ComputeDetTex(A2buff, x2buff5, x3buff5, x4buff5);

		for (int i = 0; i < 12; ++i)
			volume += PieceOfVolume[i];
	}
	// TRIANGLES = 4
	if (m_ele->GetElementType() == MshElemType::TRIANGLE)
	{
		PieceOfVolume = new double[3]();

		double const* const x1buff(m_ele->GetNode(0)->getData());
		double const* const x2buff(m_ele->GetNode(1)->getData());
		double const* const x3buff(m_ele->GetNode(2)->getData());

		PieceOfVolume[0] = ComputeDetTri(x1buff, x2buff, A2buff);
		PieceOfVolume[1] = ComputeDetTri(x2buff, x3buff, A2buff);
		PieceOfVolume[2] = ComputeDetTri(x3buff, x1buff, A2buff);

		for (int i = 0; i < 3; ++i)
			volume += PieceOfVolume[i];
	}
	// TETRAHEDRAS = 5
	if (m_ele->GetElementType() == MshElemType::TETRAHEDRON)
	{
		PieceOfVolume = new double[4]();

		//      node = m_ele->GetNode(3);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff(m_ele->GetNode(0)->getData());
		double const* const x2buff(m_ele->GetNode(1)->getData());
		double const* const x3buff(m_ele->GetNode(2)->getData());
		double const* const x4buff(m_ele->GetNode(3)->getData());

		PieceOfVolume[0] = ComputeDetTex(A2buff, x1buff, x2buff, x3buff);
		PieceOfVolume[1] = ComputeDetTex(A2buff, x1buff, x3buff, x4buff);
		PieceOfVolume[2] = ComputeDetTex(A2buff, x1buff, x4buff, x2buff);
		PieceOfVolume[3] = ComputeDetTex(A2buff, x2buff, x3buff, x4buff);

		for (int i = 0; i < 4; ++i)
			volume += PieceOfVolume[i];
	}
	// PRISMS = 6
	if (m_ele->GetElementType() == MshElemType::PRISM)
	{
		PieceOfVolume = new double[8]();

		// 2,1,3 face
		//      node = m_ele->GetNode(1);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(0);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();

		double const* const x1buff(m_ele->GetNode(1)->getData());
		double const* const x2buff(m_ele->GetNode(0)->getData());
		double const* const x3buff(m_ele->GetNode(2)->getData());

		PieceOfVolume[0] = ComputeDetTex(A2buff, x1buff, x2buff, x3buff);

		// 4,5,6 face
		//      node = m_ele->GetNode(3);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(4);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();

		double const* const x1buff1(m_ele->GetNode(3)->getData());
		double const* const x2buff1(m_ele->GetNode(4)->getData());
		double const* const x3buff1(m_ele->GetNode(5)->getData());

		PieceOfVolume[1] = ComputeDetTex(A2buff, x1buff1, x2buff1, x3buff1);

		// 1,4,6,3 face
		//      node = m_ele->GetNode(0);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(3);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff2(m_ele->GetNode(0)->getData());
		double const* const x2buff2(m_ele->GetNode(3)->getData());
		double const* const x3buff2(m_ele->GetNode(5)->getData());
		double const* const x4buff2(m_ele->GetNode(2)->getData());

		PieceOfVolume[2] = ComputeDetTex(A2buff, x1buff2, x2buff2, x4buff2);
		PieceOfVolume[3] = ComputeDetTex(A2buff, x2buff2, x3buff2, x4buff2);

		// 2,5,4,1 face
		//      node = m_ele->GetNode(1);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(4);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(3);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(0);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff3(m_ele->GetNode(1)->getData());
		double const* const x2buff3(m_ele->GetNode(4)->getData());
		double const* const x3buff3(m_ele->GetNode(3)->getData());
		double const* const x4buff3(m_ele->GetNode(0)->getData());

		PieceOfVolume[4] = ComputeDetTex(A2buff, x1buff3, x2buff3, x4buff3);
		PieceOfVolume[5] = ComputeDetTex(A2buff, x2buff3, x3buff3, x4buff3);

		// 5,2,3,6 face
		//      node = m_ele->GetNode(4);
		//      x1buff[0] = node->X();
		//      x1buff[1] = node->Y();
		//      x1buff[2] = node->Z();
		//      node = m_ele->GetNode(1);
		//      x2buff[0] = node->X();
		//      x2buff[1] = node->Y();
		//      x2buff[2] = node->Z();
		//      node = m_ele->GetNode(2);
		//      x3buff[0] = node->X();
		//      x3buff[1] = node->Y();
		//      x3buff[2] = node->Z();
		//      node = m_ele->GetNode(5);
		//      x4buff[0] = node->X();
		//      x4buff[1] = node->Y();
		//      x4buff[2] = node->Z();

		double const* const x1buff4(m_ele->GetNode(4)->getData());
		double const* const x2buff4(m_ele->GetNode(1)->getData());
		double const* const x3buff4(m_ele->GetNode(2)->getData());
		double const* const x4buff4(m_ele->GetNode(5)->getData());

		PieceOfVolume[6] = ComputeDetTex(A2buff, x1buff4, x2buff4, x4buff4);
		PieceOfVolume[7] = ComputeDetTex(A2buff, x2buff4, x3buff4, x4buff4);

		for (int i = 0; i < 8; ++i)
			volume += PieceOfVolume[i];
	}

	// Release the memory
	delete[] PieceOfVolume;

	return volume;
}

/**************************************************************************
   MSHLib-Method:
   Task:Compute the volume of the object via the particle inside of the object
   Programing:
   09/2005 PCH Implementation
**************************************************************************/
double RandomWalk::ComputeVolume(Particle* A, Particle* element, MeshLib::CElem* m_ele)
{
	double x1buff[3];
	double x2buff[3];
	double x3buff[3];
	double x4buff[3];
	double volume = 0.0;
	double* PieceOfVolume = NULL;

	double A2buff[3];

	A2buff[0] = A->x;
	A2buff[1] = A->y;
	A2buff[2] = A->z;

	x1buff[0] = element[0].x;
	x1buff[1] = element[0].y;
	x1buff[2] = element[0].z;
	x2buff[0] = element[1].x;
	x2buff[1] = element[1].y;
	x2buff[2] = element[1].z;
	x3buff[0] = element[2].x;
	x3buff[1] = element[2].y;
	x3buff[2] = element[2].z;

	// TRIANGLES = 4, RECTANGLE = 2
	int eleType = m_ele->GetElementType();
	if (eleType == MshElemType::TRIANGLE || eleType == MshElemType::QUAD)
	{
		PieceOfVolume = new double[3]();

		PieceOfVolume[0] = ComputeDetTri(x1buff, x2buff, A2buff);
		PieceOfVolume[1] = ComputeDetTri(x2buff, x3buff, A2buff);
		PieceOfVolume[2] = ComputeDetTri(x3buff, x1buff, A2buff);

		for (int i = 0; i < 3; ++i)
			volume += PieceOfVolume[i];
	}
	// TETRAHEDRAS = 5, HEXAHEDRA = 3, PRISM = 6
	else if (eleType == MshElemType::TETRAHEDRON || eleType == MshElemType::HEXAHEDRON || eleType == MshElemType::PRISM)
	{
		PieceOfVolume = new double[4]();

		x4buff[0] = element[3].x;
		x4buff[1] = element[3].y;
		x4buff[2] = element[3].z;

		PieceOfVolume[0] = ComputeDetTex(A2buff, x1buff, x2buff, x3buff);
		PieceOfVolume[1] = ComputeDetTex(A2buff, x1buff, x3buff, x4buff);
		PieceOfVolume[2] = ComputeDetTex(A2buff, x1buff, x4buff, x2buff);
		PieceOfVolume[3] = ComputeDetTex(A2buff, x2buff, x3buff, x4buff);

		for (int i = 0; i < 4; ++i)
			volume += PieceOfVolume[i];
	}
	else
		abort();

	// Release the memory
	delete[] PieceOfVolume;

	return volume;
}

/**************************************************************************
   MSHLib-Method:
   Task:The function advances the set of particles by advection
    and dispersion
   Programing:
   10/2005 PCH Implementation
**************************************************************************/
void RandomWalk::AdvanceBySplitTime(double dt, int numOfSplit)
{
	double subdt = dt / (double)numOfSplit;
	double ctime = 0.0; // JT 05.2010: for continuous source compatible with SplitTimes.
	leavingParticles = 0;
	for (int i = 0; i < numOfSplit; ++i)
	{
		AdvanceToNextTimeStep(subdt, ctime);
		ctime += subdt;
	}
}

/**************************************************************************
   Task:The function trace streamline or pathlines for a given number of
      particles
   Programing:
   01/2007 PCH Implementation
**************************************************************************/
void RandomWalk::TraceStreamline(void)
{
	// OK411 double tolerance = 1e-18;
	// Loop over all the particles
	for (int i = 0; i < numOfParticles; ++i)
		if (X[i].Now.elementIndex != -10)
		{
			// Let's record the current to the past
			//			InterpolateVelocity(&(X[i].Now));
			TracePathlineInThisElement(&(X[i].Now));

			// Record path
			if (i < 50)
				RecordPath(i, &(X[i].Now));
		}

#ifdef _FEMPCHDEBUG_
	// PCH Let's monitor what's going on in the FEM
	// This messagebox is for debugging the primary variables at every time step.
	// Should combine with the picking...
	CWnd* pWnd = NULL;
	pWnd->MessageBox("Split second!!!", "Debug help", MB_ICONINFORMATION);
#endif
}

/**************************************************************************
   Task:The function advances the set of particles by advection
    and dispersion
   Programing:
   10/2005 PCH Implementation
   05/2009 PCH mobiility of particle now defined in .mcp via components
   06/2009 PCH Case specific apps for RWPT defined by rwpt_app
**************************************************************************/
void RandomWalk::AdvanceToNextTimeStep(double dt, double ctime)
{
	double tolerance = 1e-18;
	double tol = 1e-10;
	int TimeMobility;
	// 05.2010 JT
	double exceedtime = aktuelle_zeit + MKleinsteZahl + ctime;
// Loop over all the particles
#ifdef _OPENMP
	std::cout << "[RandomWalk::AdvanceToNextTimeStep] OpenMP parallelized loop ... " << std::flush;
#else
	std::cout << "[RandomWalk::AdvanceToNextTimeStep] loop ... " << std::flush;
#endif
	clock_t start, stop;
	start = clock();

	m_pcs = PCSGet("FLUID_MOMENTUM");
	if (!m_pcs && flow_pcs)
	{
		m_pcs = flow_pcs;
	}

	if (!m_pcs)
	{
		std::cerr << "No flow process defined\n";
		exit(1);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < numOfParticles; ++i)
	{
		TimeMobility
		    = 0; // JT 2010, using this for now. Setting identity = 1 causes simulation failure... not sure why??
		// X[i].Now.identity=1;
		if ((X[i].Now.StartingTime < exceedtime) || fabs(X[i].Now.StartingTime - exceedtime) < tol)
			TimeMobility = 1;
		// X[i].Now.identity=0;

		// components defined in .mcp should be syncronized with identity of particles.
		CompProperties* m_cp = cp_vec[0];

		// If mobile, do transport.
		if (m_cp->mobil && TimeMobility > 0)
		{
			Particle Y; // the displaced particle
			int Astatus = 100; // Set to be meaningless

			if (X[i].Now.elementIndex != -10 && X[i].Now.identity != 2)
			{
				// Let's record the current to the past
				if (PURERWPT != 2)
					InterpolateVelocityOfTheParticleByInverseDistance(&(X[i].Now));
				else
					InterpolateVelocityOfTheParticleByBilinear(GridOption, &(X[i].Now));

				// If the mode is for heterogeneous media
				if (RWPTMode == 0 || RWPTMode == 1 || RWPTMode == 3)
					SolveForDerivativeOfVelocity(&(X[i].Now));
				if (RWPTMode < 2 || RWPTMode > 3) // 0 or 1 for advection and dispersion cases.
					SolveDispersionCoefficient(&(X[i].Now));

				// Initialize the reference and past particles
				Y = X[i].Past = X[i].Now;

				// Initialize
				Y.t = dt;

				// Record path
				// if(i<50) // JT :: no longer needed... just use .particle outpus and Paraview
				// RecordPath(i, &Y);

				do
				{
					// Let's update the info of Particle Y.
					if (PURERWPT != 2)
						InterpolateVelocityOfTheParticleByInverseDistance(&Y);
					else
						InterpolateVelocityOfTheParticleByBilinear(GridOption, &Y);

					if (RWPTMode < 4)
						SolveForDerivativeOfVelocity(&Y);
					if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
						SolveDispersionCoefficient(&Y);

					if (Astatus == -1)
						Y.t = dt;
					// YS: absorbed and filtered particles don't move, but absorbed particles can be desorpted
					if (X[i].Now.identity == 0)
						Astatus = SolveForNextPosition(&(X[i].Now), &Y);

					//#ifdef CountParticleNumber
					if ((m_pcs->rwpt_count == 1) && (Y.identity != 2))
					{
						if ((fabs(Y.x - 0.1) < 0.001 && fabs(Y.y - 0.0) < 0.001 && fabs(Y.z - 0.0) < 0.001))
						{
							leavingParticles++;
							Y.elementIndex = -10; // YS: out of the domain
						}
					}
					//#endif

					// Just get the element index after this movement
					// if not Homogeneous aquifer
					if (RWPTMode % 2 == 1)
					{
						if (PURERWPT != 2)
							Y.elementIndex = IndexOfTheElementThatThisParticleBelong(1, &Y);
						//			Y.elementIndex = GetTheElementOfTheParticle(&(X[i].Now), &Y);
						else
							Y.elementIndex = IndexOfTheElementThatThisParticleBelong(0, &Y);
					}

					// We let the particle go outside of the domain
					if (Y.elementIndex == -10)
					{
						// Before letting this particle outside of the domain,
						// record the particle postion that includes the element index = -10.
						X[i].Now = Y;
						break;
					}

					// 07.2011 YS
					if (Y.on_boundary != 0)
					{
						// record particle position
						X[i].Now = Y;
						break;
					}
					// The result of the function is unknown error.
					if (Astatus == -2)
					{
						printf("Astatus = %d\n", Astatus);
						abort();
					}
					// Particle goes outside of the domain.
					// Thus, do it again.
					else if (Astatus == -1)
						Y = X[i].Now;
					// Right on track. Keep going.
					else
					{
						// If particle stays in the element, Y.t = dt.
						if (Y.t < tolerance)
							Y.t = dt;

						// Update the current info
						// where the advected particle is in either the element or
						// the neighboring element.
						X[i].Now = Y;
					}

					// Keep looping if the time is not spent all or
					// if the particle is outside of the domain or
					// if the function fails
				} while (Y.t < dt);

				// Record path
				// if(i<50) // JT :: no longer needed
				// RecordPath(i, &Y);
				// Update the correct dt
				X[i].Now.t = X[i].Past.t + dt;
			}
		}

		// Now ODE parts in RWPT
		if (m_pcs->rwpt_app == 1) // Is the application is Cell Dynamics?
		{
			/*
			         // Do mobile-Immobile by switching the identity of particles
			         double ChanceOfMobile = randomZeroToOne();
			         double ChanceOfImmobile = randomZeroToOne();

			         // Rate coefficients definition here
			         // For the use of the existing parameters, I borrow Freundlich non-isotherm parameter set,
			         // which is isotherm_model = 2
			         int numOfComps = (int)cp_vec.size(); //OK411
			         double* Kon = NULL; double* Koff = NULL;
			         Kon = new double [numOfComps];
			   Koff = new double [numOfComps];
			   double FeqSum=0.0;
			   for(int i=0; i< numOfComps; ++i)
			   {
			   Kon[i] = cp_vec[i]->isotherm_model_values[0];
			   Koff[i] = cp_vec[i]->isotherm_model_values[1];
			   FeqSum+=Kon[i]/Koff[i];
			   }
			   if(X[i].Now.identity == 0)	// Among mobile particles,
			   {
			   double Feq=1./(1.+FeqSum);

			   if( ChanceOfMobile < (1.0-exp(-Koff[0]*dt))*Kon[0]/Koff[0] )
			   {
			   if(numOfComps > 1)
			   X[i].Now.identity = 1;	// Make it immobile
			   else
			   ;
			   }
			   else if( ChanceOfMobile < ( (1.0-exp(-Koff[0]*dt))*Kon[0]/Koff[0] +
			   (1.0-exp(-Koff[1]*dt))*Kon[1]/Koff[1]) )
			   {
			   if(numOfComps > 2)	// If the number of components is bigger than 2,
			   X[i].Now.identity = 2;	// Make another kind of immobile
			   else
			   ;
			   }
			   else
			   ;	// Leave it as mobile
			   }
			   else if(X[i].Now.identity == 1 && numOfComps > 1)	// Among immobile particles,
			   {
			   if( ChanceOfImmobile < (1.0-exp(-Koff[0]*dt)) )
			   X[i].Now.identity = 0;	// Make it mobile
			   else
			   ;	// Leave it as mobile
			   }
			   else if(X[i].Now.identity == 2 && numOfComps > 2)
			   {
			   if( ChanceOfImmobile < (1.0-exp(-Koff[1]*dt)) )
			   X[i].Now.identity = 0;	// Make it mobile
			   else
			   ;	// Leave it as mobile
			   }
			   else
			   {
			   cout<< "Only Identity 0 and 1 are covered. There are more than 2 identities detected" << "\n";
			   abort();
			   }

			   // Release memory
			   delete [] Kon;	delete [] Koff;
			 */
		}
		else if (m_pcs->rwpt_app == 2) // Is the application Cryptosporidium oocysts?
		{
			if (X[i].Now.elementIndex != -10 && X[i].Now.identity != 2)
			{
				// Do sorption-desorption by switching the identity of particles
				double ChanceOfSorbtion = randomZeroToOne();
				// Two-Rate Model: N/N0=Ae^(-k1t)+(1-A)e^(-k2t)
				if (m_cp->isotherm_model == 5)
				{
					double A = m_cp->isotherm_model_values[0];
					double k1 = m_cp->isotherm_model_values[1];
					double k2 = m_cp->isotherm_model_values[2];
					double FractionRemainingOnMedia = Two_rateModel(A, k1, k2, X[i].Now.t);

					if (ChanceOfSorbtion < FractionRemainingOnMedia)
						X[i].Now.identity = 1; // Temporarily absorbed
					else
						X[i].Now.identity = 0;
				}

				// Irreversible Reactions - Oocysts filtered for good or some chemicals decayed
				// C/C0 = exp(-kt)
				// For oocyst irreversible filtration, k = vp * lambda
				double lamda = m_cp->decay_model_values[0];
				double k = X[i].Now.Vx * lamda;
				double Irreversed = exp(-k * X[i].Now.t);

				if (ChanceOfIrreversed[i] > Irreversed)
					X[i].Now.identity = 2; // Permanently filtered
			}
		}
		else // For some other applications with different kinetics
		{
		}
	}
	stop = clock();
	std::cout << " took " << (stop - start) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
}

/**************************************************************************
   MSHLib-Method:
   Task:The function records pathlines of 100 particles or less
   Programing:
   01/2007 PCH Implementation
**************************************************************************/
void RandomWalk::RecordPath(int no, Particle* P)
{
	Position p;

	p.p[0] = P->x;
	p.p[1] = P->y;
	p.p[2] = P->z;

	pathline[no].path.push_back(p);
}

/**************************************************************************
   MSHLib-Method:
   Task: Determines when to generate output files, and introduces call to data output
   Programing:
   03/2010 JT
**************************************************************************/
void RandomWalk::RandomWalkOutput(double dbl_time, int current_time_step)
{
	bool outputornot;
	CurrentTime = dbl_time;
	std::string current_name;

	for (size_t i = 0; i < rwpt_out_strings.size(); i++)
	{
		current_name = rwpt_out_strings[i];
		COutput* out = OUTGetRWPT(current_name);
		if (!out)
			continue;

		outputornot = false;
		out->setTime(CurrentTime);
		size_t no_times = out->getRWPTTimeVector().size();

		if (no_times == 0 && (out->getNSteps() > 0) && (current_time_step % out->getNSteps() == 0))
			outputornot = true;
		if (current_time_step < 2)
			outputornot = true;

		if (outputornot)
			if (current_name.compare("PARTICLES") == 0)
				DATWriteParticleFile(current_time_step);
		// else if(current_name.compare("PARTICLE_CONCENTRATION")==0)
		// DATWriteParticleConcFile(current_time_step); // routine not yet configured
		{
			for (size_t j = 0; j < no_times; j++)
				if (CurrentTime >= out->getRWPTTimeVector()[j])
				{
					if (current_name.compare("PARTICLES") == 0)
						DATWriteParticleFile(current_time_step);

					// else if(current_name.compare("PARTICLE_CONCENTRATION")==0)
					// DATWriteParticleConcFile(current_time_step); // routine not yet configured (see commented text
					// below)

					out->getRWPTTimeVector().erase(out->getRWPTTimeVector().begin() + j);
					break;
				}
		}
	}

	//  OUTPUT PARTICLES AS ELEMENTAL CONCENTRATION
	//  ---------------------------------------------------
	/*
	   m_out = OUTGetRWPT("PARTICLE_CONCENTRATION");
	   outputornot=false;
	   m_out->time = CurrentTime;
	   no_times = (int)m_out->time_vector.size();
	   if(current_time_step%m_out->nSteps==0)
	   outputornot = true;
	   if(current_time_step<2)
	   outputornot = true;

	   if(outputornot)
	   DATWriteParticleFile(current_time_step);
	   else
	   {
	   for(j=0;j<no_times;j++)
	   {
	   if(CurrentTime>=m_out->time_vector[j])
	   DATWriteParticleFile(current_time_step);
	   }
	   }

	   double UnitConcentration = 0.0;

	   // Here's definition for unit concentration
	   UnitConcentration = (double)numOfParticles / (double) m_msh->ele_vector.size();

	   for(int i=0; i< (int)m_msh->ele_vector.size(); ++i)
	   {
	   // Now count the number of particles in each element.
	   int CountInThisElement = 0;
	   for(int j=0; j< numOfParticles; ++j)
	   {
	   if(X[j].Now.elementIndex == i)
	   ++CountInThisElement;
	   }

	   // Store the normalized concentration
	   double NormConcentrationOfTheElement = (double)CountInThisElement / UnitConcentration;
	   //		double NormConcentrationOfTheElement = (double)CountInThisElement/numOfParticles;
	   SetElementValue(i, GetElementValueIndex("CONCENTRATION0")+1, NormConcentrationOfTheElement);
	   }
	 */

	/*
	   if( ((int)(X[0].Now.t*1000))== 5)
	   {
	   char now[100];
	   sprintf(now, "%f", X[0].Now.t);
	   ConcPTFile(now);
	   }
	 */
}

void RandomWalk::ConcPTFile(const char* file_name)
{
	FILE* pct_file = NULL;
	char pct_file_name[MAX_ZEILE];

	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	sprintf(pct_file_name, "%s.conc", file_name);
	pct_file = fopen(pct_file_name, "w+t");

	// Make a grid
	int gridDensity = 0;
	// Search Max and Min of each axis
	double MaxX = -1e6, MinX = 1e6;
	for (size_t i = 0; i < m_msh->nod_vector.size(); ++i)
	{
		// TF	   MeshLib::CNode* thisNode = m_msh->nod_vector[i]
		double x_coord(m_msh->nod_vector[i]->getData()[0]);
		if (MaxX < x_coord)
			MaxX = x_coord;
		if (MinX > x_coord)
			MinX = x_coord;
	}

	gridDensity = (int)(MaxX - MinX);

	fprintf(pct_file, "VARIABLES = X,Y,Z,CONCENTRATION0\n");
	fprintf(pct_file, "ZONE T=\"%fs\", I=%d, F=POINT, C=BLACK\n", X[0].Now.t, gridDensity);

	for (int i = 0; i < gridDensity; ++i)
	{
		int count = 0;
		double seg_start = 0.0, seg_end = 0.0;
		for (int j = 0; j < numOfParticles; ++j)
		{
			seg_start = MinX + i;
			seg_end = MinX + i + 1.0;

			if ((X[j].Now.x >= seg_start) && (X[j].Now.x < seg_end))
				++count;
		}
		//	fprintf(pct_file, "%f 0.0 0.0 %f\n", (seg_start+seg_end)/2.0, count / numOfParticles);
		fprintf(pct_file, "%f 0.0 0.0 %f\n", (seg_start + seg_end) / 2.0, (double)count);
	}

	// Let's close it, now
	fclose(pct_file);
}

/**************************************************************************
   Task:Solve fraction remaining on media
     N/N0=Ae^(-k1t)+(1-A)e^(-k2t)
   Programing:
   09/2007 PCH Implementation
**************************************************************************/
double RandomWalk::Two_rateModel(double A, double k1, double k2, double t)
{
	return A * exp(-k1 * t) + (1.0 - A) * exp(-k2 * t);
}

/**************************************************************************
   MSHLib-Method:
   Task:Give the beat boys and free my soul. I wanna get lost in your rock &
     roll (random displace) and DRIFT AWAY. Return three component of
     random drift at the particle position.
   Programing:
   10/2005 PCH Implementation
**************************************************************************/
void RandomWalk::RandomlyDriftAway(Particle* A, double dt, double* delta, int type)
{
	// OK411
	type = type;
	dt = dt;

	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];

	// Let's generate three random components N(0,1) and use it to compute deltaOfX
	double Z[3];

	// Here I tell the dimension for the element that contains the particle A
	int ele_dim = m_ele->GetDimension();
	if (ele_dim == 1)
	{
		if (UniformOrNormal == 1)
		{
			Z[0] = randomMinusOneToOne();
			delta[0] = sqrt(6.0 * A->D[0] * dt) * Z[0];
			delta[1] = 0.0;
			delta[2] = 0.0;
		}
		else
		{
			Z[0] = Marsaglia();
			delta[0] = sqrt(2.0 * A->D[0] * dt) * Z[0];
			delta[1] = 0.0;
			delta[2] = 0.0;
		}
	}
	else if (ele_dim == 2)
	{
		if (UniformOrNormal == 1)
		{
			Z[0] = randomMinusOneToOne();
			Z[1] = randomMinusOneToOne();
			delta[0] = sqrt(6.0 * A->D[0] * dt) * Z[0] + sqrt(6.0 * A->D[1] * dt) * Z[1];
			delta[1] = sqrt(6.0 * A->D[3] * dt) * Z[0] + sqrt(6.0 * A->D[4] * dt) * Z[1];
			delta[2] = 0.0;
		}
		else
		{
			Z[0] = Marsaglia();
			Z[1] = Marsaglia();
			delta[0] = sqrt(2.0 * A->D[0] * dt) * Z[0] + sqrt(2.0 * A->D[1] * dt) * Z[1];
			delta[1] = sqrt(2.0 * A->D[3] * dt) * Z[0] + sqrt(2.0 * A->D[4] * dt) * Z[1];
			delta[2] = 0.0;
		}
	}
	else if (ele_dim == 3)
	{
		if (UniformOrNormal == 1)
		{
			Z[0] = randomMinusOneToOne();
			Z[1] = randomMinusOneToOne();
			Z[2] = randomMinusOneToOne();
			delta[0]
			    = sqrt(6.0 * A->D[0] * dt) * Z[0] + sqrt(6.0 * A->D[1] * dt) * Z[1] + sqrt(6.0 * A->D[2] * dt) * Z[2];
			delta[1]
			    = sqrt(6.0 * A->D[3] * dt) * Z[0] + sqrt(6.0 * A->D[4] * dt) * Z[1] + sqrt(6.0 * A->D[5] * dt) * Z[2];
			delta[2]
			    = sqrt(6.0 * A->D[6] * dt) * Z[0] + sqrt(6.0 * A->D[7] * dt) * Z[1] + sqrt(6.0 * A->D[8] * dt) * Z[2];
		}
		else
		{
			Z[0] = Marsaglia();
			Z[1] = Marsaglia();
			Z[2] = Marsaglia();
			delta[0]
			    = sqrt(2.0 * A->D[0] * dt) * Z[0] + sqrt(2.0 * A->D[1] * dt) * Z[1] + sqrt(2.0 * A->D[2] * dt) * Z[2];
			delta[1]
			    = sqrt(2.0 * A->D[3] * dt) * Z[0] + sqrt(2.0 * A->D[4] * dt) * Z[1] + sqrt(2.0 * A->D[5] * dt) * Z[2];
			delta[2]
			    = sqrt(2.0 * A->D[6] * dt) * Z[0] + sqrt(2.0 * A->D[7] * dt) * Z[1] + sqrt(2.0 * A->D[8] * dt) * Z[2];
		}
	}
}

/**************************************************************************
   MSHLib-Method:
   Task:The function solves random displacement by random number generation.

   Programing:
   12/2005 PCH Implementation
**************************************************************************/
int RandomWalk::RandomWalkDrift(double* Z, int dim)
{
	if (dim == 1) // Generate the faster one.
	{
		Z[0] = randomMinusOneToOne();
		Z[1] = Z[2] = 0.0;

		return 1;
	}
	else if (dim == 2) // Generate the normal distribution one
	{
		Z[0] = randomMinusOneToOne();
		Z[1] = randomMinusOneToOne();
		Z[2] = 0.0;

		return 1;
	}
	else if (dim == 3)
	{
		Z[0] = randomMinusOneToOne();
		Z[1] = randomMinusOneToOne();
		Z[2] = randomMinusOneToOne();

		return 1;
	}
	else
	{
		printf("Something wrong in generation random drift\n");
		abort();
	}

	return -1; // Failed.
}

/**************************************************************************
   Task: SolveDispersionCoefficient(Particle* A)
   Programing: This function solves velocity tensor from the velocity of
         particle.
   10/2005 PCH
   02/2006 PCH	Extension to cover 2D elements in 3D
   05/2009 PCH mobility of components now defined in .mcp. However,
         the particle idensity should be syncronized with components
         .mcp files.
**************************************************************************/
void RandomWalk::SolveDispersionCoefficient(Particle* A)
{
	// To extract dispersivities from material properties
	CMediumProperties* m_mat_mp = NULL;
	double alphaL = 0.0, alphaT = 0.0;
	double V[3];
	double tolerance = 1e-18;

	// Extract the dispersivities from the group that the particle belongs
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* m_ele = m_msh->ele_vector[A->elementIndex];
	int group = m_ele->GetPatchIndex();
	m_mat_mp = mmp_vector[group];
	alphaL = m_mat_mp->mass_dispersion_longitudinal;
	alphaT = m_mat_mp->mass_dispersion_transverse;

	// Let's solve pore velocity.
	// It is simple because Sw stuff automatically handles in Richards Flow.
	// Thus, I only divide Darcy velocity by porosity only to get pore velocity.
	CMediumProperties* MediaProp = mmp_vector[m_ele->GetPatchIndex()];
	double porosity = 0.0;
	if (MediaProp->porosity > 10 - 6)
		porosity = MediaProp->porosity; // This is for simple one.
	else
		// This will get you porosity.
		porosity = MediaProp->porosity_model_values[0];
	// I guess for Dual Porocity stuff can also be handled here.
	double molecular_diffusion_value = 0.0;
	// components defined in .mcp should be syncronized with identity of particles.
	CompProperties* m_cp = cp_vec[0];
	double g[3] = {0., 0., 0.};
	double theta = 1.0; // I'll just set it to be unity for moment.
	molecular_diffusion_value = m_cp->CalcDiffusionCoefficientCP(A->elementIndex, 1.0, m_pcs)
	                            * MediaProp->TortuosityFunction(A->elementIndex, g, theta);
	molecular_diffusion_value /= porosity; // This should be divided by porosity in this RWPT method.

	// Just solve for the magnitude of the velocity to compute the dispersion tensor
	V[0] = A->Vx;
	V[1] = A->Vy;
	V[2] = A->Vz;

	// RWPT-IM
	// Let's transform this velocity to be on the xy plane
	// Some nice if condition to tell the need for transform will be nice. Later....
	if (m_ele->GetDimension() < 3 && (m_ele->GetDimension() != m_msh->GetMaxElementDim()))
		ToTheXYPlane(m_ele, V);

	double Vmagnitude = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);

	// Compute the dispersion tensor at the particle location
	// If the magnitude of velocity is not zero.
	if (Vmagnitude > tolerance)
	{
		// Dxx
		A->D[0]
		    = (alphaT * (V[1] * V[1] + V[2] * V[2]) + alphaL * V[0] * V[0]) / Vmagnitude + molecular_diffusion_value;
		// Dxy = Dyz
		A->D[1] = A->D[3] = (alphaL - alphaT) * V[0] * V[1] / Vmagnitude;
		// Dxz = Dzx
		A->D[2] = A->D[6] = (alphaL - alphaT) * V[0] * V[2] / Vmagnitude;
		// Dyy
		A->D[4]
		    = (alphaT * (V[0] * V[0] + V[2] * V[2]) + alphaL * V[1] * V[1]) / Vmagnitude + molecular_diffusion_value;
		// Dyz = Dzy
		A->D[5] = A->D[7] = (alphaL - alphaT) * V[1] * V[2] / Vmagnitude;
		// Dzz
		A->D[8]
		    = (alphaT * (V[0] * V[0] + V[1] * V[1]) + alphaL * V[2] * V[2]) / Vmagnitude + molecular_diffusion_value;
	}
	else
	{
		A->D[0] = molecular_diffusion_value; // Dxx
		A->D[1] = A->D[3] = 0.0; // Dxy = Dyz
		A->D[2] = A->D[6] = 0.0; // Dxz = Dzx
		A->D[4] = molecular_diffusion_value; // Dyy
		A->D[5] = A->D[7] = 0.0; // Dyz = Dzy
		A->D[8] = molecular_diffusion_value; // Dzz
	}
}
//------------------------------------------------------------------------------
/*
Task: Check whether the particle is out of the domain and get the pierced point
Programing:
11/2011 YS/WW Implementation
*/
void RandomWalk::CheckBoundary2D(Particle* A, Particle* B)
{
	// firstly, check if the particle is inside of the domain
	int eleindex = SearchElementFromNeighbor(B);

	if (eleindex > -1)
	{
		B->elementIndex = eleindex;
		B->on_boundary = 0;
		return;
	}

	// If the particle is out of domain, find the intersection point between the path line and
	// the boundary.
	// But check the old position first, in case it is already on the boundary.

	if (A->on_boundary)
	{
		// drag B back to the A position
		B->x = A->x;
		B->y = A->y;
		B->elementIndex = A->elementIndex;
		B->on_boundary = A->on_boundary;
		B->identity = A->identity;
	}
	else
	{
		double p1[2], p2[2], p3[2], p4[2];
		double xx, yy, ua, ub;

		p3[0] = A->x;
		p3[1] = A->y;
		p4[0] = B->x;
		p4[1] = B->y;

		MeshLib::CElem* n_ele;
		MeshLib::CNode* a_node;

		for (size_t j = 0; j < m_msh->face_vector.size(); j++)
		{
			n_ele = m_msh->face_vector[j];
			a_node = n_ele->GetNode(0);
			p1[0] = a_node->X();
			p1[1] = a_node->Y();
			a_node = n_ele->GetNode(1);
			p2[0] = a_node->X();
			p2[1] = a_node->Y();

			xx = yy = ua = ub = 0.;
			if (G_intersect_line_segments(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], &ua, &ub, &xx, &yy)
			    == 1)
			{
				if (ua >= 0.0 && ua <= 1.0)
				{
					B->x = xx;
					B->y = yy;
					B->elementIndex = n_ele->GetOwner()->GetIndex();
					B->on_boundary = 1;
					break;
				}
			}
		}
	}
}
//-------------------------------------------------------------------------------
/*
Task: Check whether the particle is out of the domain and get the pierced point

Implementation:
interesection between straight line and surface
straight line is described as L(lambda)=p+lambda*(q-p)
surface is described as E(miu,ita)=r+miu*(s-r)+ita*(t-r)
with the abbrevations us=s-r,vs=t-r,ws=p-q,bs=p-r
to get the intersection is to solve this equation in matrix vector form
  (us0  vs0  ws0) (miu   )   (bs0)
  (us1  vs1  ws1).(ita   ) = (bs1)
  (us2  vs2  ws2) (lambda)   (bs2)
the intersection point is in the surface inside the triangle spanned by the three points
if and only if 0<miu,ita,lambda<1 and miu+ita<=1

Programing: 11/2011 YS/WW
*/
void RandomWalk::CheckBoundary3D(Particle* A, Particle* B)
{
	// Only for tetrahedral and prismatic mesh for now
	// firstly, check if the particle is inside of the domain
	int eleindex = SearchElementFromNeighbor(B);

	if (eleindex > -1)
	{
		B->elementIndex = eleindex;
		B->on_boundary = 0;
		return;
	}

	// The new position is in none of the elements, i.e., outside of the domain.
	// Find the intersection point between the path line and the boundary.
	// But check the old position first, in case it is already on the boundary.

	if (A->on_boundary)
	{
		// drag B back to the A position
		B->x = A->x;
		B->y = A->y;
		B->z = A->z;
		B->elementIndex = A->elementIndex;
		B->on_boundary = A->on_boundary;
		B->identity = A->identity;
	}
	else
	{
		int kk, itr, ntr;
		double x1[3], x2[3], x3[3], x4[3];
		MeshLib::CNode* a_node;
		MeshLib::CElem* n_ele;
		double us[3], vs[3], ws[3], bs[3];
		double f_buff;
		bool s_cond;

		x4[0] = A->x; // point p on the path line
		x4[1] = A->y;
		x4[2] = A->z;

		double x_new[3];
		x_new[0] = B->x; // point q on the path line
		x_new[1] = B->y;
		x_new[2] = B->z;

		// search all the element facets which are on the boundary of the domain
		for (int j = 0; j < (long)m_msh->face_vector.size(); j++)
		{
			n_ele = m_msh->face_vector[j];

			ntr = 1;
			// if the element facet is quadrilateral, divide it into two triangles
			if (n_ele->GetElementType() == MshElemType::QUAD)
				ntr = 2;

			for (itr = 0; itr < ntr; itr++)
			{
				if (itr == 0)
				{
					a_node = n_ele->GetNode(0); // point r on the boundary plane
					x1[0] = a_node->X();
					x1[1] = a_node->Y();
					x1[2] = a_node->Z();
					a_node = n_ele->GetNode(1); // point s on the boundary plane
					x2[0] = a_node->X();
					x2[1] = a_node->Y();
					x2[2] = a_node->Z();
					a_node = n_ele->GetNode(2); // point t on the boundary plane
					x3[0] = a_node->X();
					x3[1] = a_node->Y();
					x3[2] = a_node->Z();
				}
				else
				{
					a_node = n_ele->GetNode(0); // point r on the boundary plane
					x1[0] = a_node->X();
					x1[1] = a_node->Y();
					x1[2] = a_node->Z();
					a_node = n_ele->GetNode(2); // point s on the boundary plane
					x2[0] = a_node->X();
					x2[1] = a_node->Y();
					x2[2] = a_node->Z();
					a_node = n_ele->GetNode(3); // point t on the boundary plane
					x3[0] = a_node->X();
					x3[1] = a_node->Y();
					x3[2] = a_node->Z();
				}

				for (kk = 0; kk < 3; kk++)
				{
					us[kk] = x2[kk] - x1[kk]; // s-r
					vs[kk] = x3[kk] - x1[kk]; // t-r
					ws[kk] = x4[kk] - x_new[kk]; // p-q
					bs[kk] = x4[kk] - x1[kk]; // p-r
				}

				if (fabs(us[0]) < DBL_MIN)
				{
					if (fabs(us[1]) > DBL_MIN)
					{
						SWAP(us[0], us[1]);
						SWAP(vs[0], vs[1]);
						SWAP(ws[0], ws[1]);
						SWAP(bs[0], bs[1]);
					}
					else if (fabs(us[2]) > DBL_MIN)
					{
						SWAP(us[0], us[2]);
						SWAP(vs[0], vs[2]);
						SWAP(ws[0], ws[2]);
						SWAP(bs[0], bs[2]);
					}
					else
						continue;
				}

				// Solve the mini equation
				f_buff = -us[1] / us[0];
				x1[0] = vs[1] + vs[0] * f_buff; // a11
				x1[1] = ws[1] + ws[0] * f_buff; // a12
				x1[2] = bs[1] + bs[0] * f_buff; // b1

				f_buff = -us[2] / us[0];
				x2[0] = vs[2] + vs[0] * f_buff; // a21
				x2[1] = ws[2] + ws[0] * f_buff; // a22
				x2[2] = bs[2] + bs[0] * f_buff; // b2

				f_buff = x1[0] * x2[1] - x1[1] * x2[0];

				if (fabs(f_buff) < DBL_MIN)
					continue;

				if (fabs(x1[0]) < DBL_MIN)
				{
					// lambda
					x3[2] = x1[2] / x1[1];
					// ita
					x3[1] = (x2[2] - x2[1] * x3[2]) / x2[0];
				}

				else
				{
					// lambda
					x3[2] = (x1[0] * x2[2] - x1[2] * x2[0]) / f_buff;
					// ita
					x3[1] = (x1[2] - x1[1] * x3[2]) / x1[0];
				}

				// mu
				x3[0] = (bs[0] - vs[0] * x3[1] - ws[0] * x3[2]) / us[0];

				// Whether 0<=lambda, ita, mu<=1
				s_cond = true;
				for (kk = 0; kk < 3; kk++)
				{
					if (x3[kk] < 0.0 || x3[kk] > 1.0)
					{
						s_cond = false;
						break;
					}
				}
				if (!s_cond)
					continue;

				if (x3[0] + x3[1] > 1.0)
					continue;

				// Assign the intersection point to B
				// B_X = p + \lambda * (q -p)
				B->x = x4[0] + x3[2] * (x_new[0] - x4[0]);
				B->y = x4[1] + x3[2] * (x_new[1] - x4[1]);
				B->z = x4[2] + x3[2] * (x_new[2] - x4[2]);

				B->elementIndex = n_ele->GetOwner()->GetIndex();
				B->on_boundary = 1;
				break;
			}
		}
	}
}

/**************************************************************************
Task: Get the index of the element that contains the particle
      by searching the neighbors
Programing: 03/2012 YS Implementation
 **************************************************************************/
using namespace MeshLib;
int RandomWalk::SearchElementFromNeighbor(Particle* A)
{
	int index;
	double xyz[3];
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	CElem* theElement = m_msh->ele_vector[A->elementIndex];

	// Let's check this element first.
	index = CheckElementIndex(A);
	if (index != -1)
		return index;

	// First neighbor's search around the main element
	for (size_t i = 0; i < theElement->GetFacesNumber(); ++i)
	{
		CElem* theNeighbor = theElement->GetNeighbor(i);

		if (theNeighbor->GetDimension() == theElement->GetDimension())
		{
			A->elementIndex = theNeighbor->GetIndex();

			index = CheckElementIndex(A);
			if (index != -1)
				return index;

			// Second, search the neighbor's neighbor
			for (size_t j = 0; j < theNeighbor->GetFacesNumber(); ++j)
			{
				CElem* theNNeighbor = theNeighbor->GetNeighbor(j);

				if (theNNeighbor->GetDimension() == theNeighbor->GetDimension())
				{
					A->elementIndex = theNNeighbor->GetIndex();

					index = CheckElementIndex(A);
					if (index != -1)
						return index;

					// Third, search the neighbor's neighbor's neighbor
					for (size_t k = 0; k < theNNeighbor->GetFacesNumber(); ++k)
					{
						CElem* theNNNeighbor = theNNeighbor->GetNeighbor(k);

						if (theNNNeighbor->GetDimension() == theNNeighbor->GetDimension())
						{
							A->elementIndex = theNNNeighbor->GetIndex();

							index = CheckElementIndex(A);
							if (index != -1)
								return index;

							// Try one more time, the neighbor's neighbor's neighbor's neighbor
							for (size_t l = 0; l < theNNNeighbor->GetFacesNumber(); ++l)
							{
								CElem* theNNNNeighbor = theNNNeighbor->GetNeighbor(l);

								if (theNNNNeighbor->GetDimension() == theNNNeighbor->GetDimension())
								{
									A->elementIndex = theNNNNeighbor->GetIndex();

									index = CheckElementIndex(A);
									if (index != -1)
										return index;
								}
							}
						}
					}
				}
			}
		}
	}

	// The particle is not in the neighboring elements. Have to search all.

	xyz[0] = A->x;
	xyz[1] = A->y;
	xyz[2] = A->z;
	index = m_msh->FindElementByPoint(xyz);
	if (index != -1)
		return index;

	// The particle is in none of the elements.
	return -1;
}
//---------------------------------------------------------------------------
/*
  Task: find out if the particle position (*xyz) belongs to current element (indicated by A->elementIndex)
  Programing:
  11/2011 YS/WW Implementation
*/
int RandomWalk::CheckElementIndex(Particle* A)
{
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();

	double xyz[3], x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], x7[3], x8[3];
	double a, a_sub[12];
	CNode* a_node;
	CElem* n_ele;
	double tol = 1e-9;

	xyz[0] = A->x;
	xyz[1] = A->y;
	xyz[2] = A->z;

	n_ele = m_msh->ele_vector[A->elementIndex];
	a = n_ele->GetVolume();

	a_node = n_ele->GetNode(0);
	x1[0] = a_node->X();
	x1[1] = a_node->Y();
	x1[2] = a_node->Z();
	a_node = n_ele->GetNode(1);
	x2[0] = a_node->X();
	x2[1] = a_node->Y();
	x2[2] = a_node->Z();
	if (n_ele->GetElementType() != MshElemType::LINE)
	{
		a_node = n_ele->GetNode(2);
		x3[0] = a_node->X();
		x3[1] = a_node->Y();
		x3[2] = a_node->Z();
	}

	if (n_ele->GetElementType() == MshElemType::QUAD || n_ele->GetElementType() == MshElemType::TETRAHEDRON)
	{
		a_node = n_ele->GetNode(3);
		x4[0] = a_node->X();
		x4[1] = a_node->Y();
		x4[2] = a_node->Z();
	}

	if (n_ele->GetElementType() == MshElemType::PYRAMID)
	{
		a_node = n_ele->GetNode(3);
		x4[0] = a_node->X();
		x4[1] = a_node->Y();
		x4[2] = a_node->Z();
		a_node = n_ele->GetNode(4);
		x5[0] = a_node->X();
		x5[1] = a_node->Y();
		x5[2] = a_node->Z();
	}

	if (n_ele->GetElementType() == MshElemType::PRISM)
	{
		a_node = n_ele->GetNode(3);
		x4[0] = a_node->X();
		x4[1] = a_node->Y();
		x4[2] = a_node->Z();
		a_node = n_ele->GetNode(4);
		x5[0] = a_node->X();
		x5[1] = a_node->Y();
		x5[2] = a_node->Z();
		a_node = n_ele->GetNode(5);
		x6[0] = a_node->X();
		x6[1] = a_node->Y();
		x6[2] = a_node->Z();
	}

	if (n_ele->GetElementType() == MshElemType::HEXAHEDRON)
	{
		a_node = n_ele->GetNode(3);
		x4[0] = a_node->X();
		x4[1] = a_node->Y();
		x4[2] = a_node->Z();
		a_node = n_ele->GetNode(4);
		x5[0] = a_node->X();
		x5[1] = a_node->Y();
		x5[2] = a_node->Z();
		a_node = n_ele->GetNode(5);
		x6[0] = a_node->X();
		x6[1] = a_node->Y();
		x6[2] = a_node->Z();
		a_node = n_ele->GetNode(6);
		x7[0] = a_node->X();
		x7[1] = a_node->Y();
		x7[2] = a_node->Z();
		a_node = n_ele->GetNode(7);
		x8[0] = a_node->X();
		x8[1] = a_node->Y();
		x8[2] = a_node->Z();
	}

	switch (n_ele->GetElementType())
	{
		case MshElemType::LINE:
			double d1, d2, d3;
			d1 = 0.;
			d2 = 0.;
			d3 = 0.;
			for (int kk = 0; kk < 3; kk++)
			{
				d1 += (x1[kk] - xyz[kk]) * (x1[kk] - xyz[kk]);
				d2 += (x2[kk] - xyz[kk]) * (x2[kk] - xyz[kk]);
				d3 += (x1[kk] - x2[kk]) * (x1[kk] - x2[kk]);
			}
			d1 = sqrt(d1);
			d2 = sqrt(d2);
			d3 = sqrt(d3);

			if (fabs((d1 + d2 - d3) / d3) < tol)
			{
				return A->elementIndex;
			}
			break;
		case MshElemType::TRIANGLE:
			a_sub[0] = ComputeDetTri(x1, x2, xyz);
			a_sub[1] = ComputeDetTri(x2, x3, xyz);
			a_sub[2] = ComputeDetTri(x3, x1, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] - a) / a) < tol)
			{
				return A->elementIndex;
			}
			break;

		case MshElemType::QUAD:
			a_sub[0] = ComputeDetTri(x1, x2, xyz);
			a_sub[1] = ComputeDetTri(x2, x3, xyz);
			a_sub[2] = ComputeDetTri(x3, x4, xyz);
			a_sub[3] = ComputeDetTri(x4, x1, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] - a) / a) < tol)
			{
				return A->elementIndex;
			}
			break;

		case MshElemType::TETRAHEDRON:
			a_sub[0] = ComputeDetTex(x2, x4, x3, xyz);
			a_sub[1] = ComputeDetTex(x1, x3, x4, xyz);
			a_sub[2] = ComputeDetTex(x2, x1, x4, xyz);
			a_sub[3] = ComputeDetTex(x2, x3, x1, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] - a) / a) < tol)
			{
				return A->elementIndex;
			}
			break;

		case MshElemType::PYRAMID:
			a_sub[0] = ComputeDetTex(x1, x2, x4, xyz);
			a_sub[1] = ComputeDetTex(x2, x3, x4, xyz);
			a_sub[2] = ComputeDetTex(x1, x5, x2, xyz);
			a_sub[3] = ComputeDetTex(x2, x5, x4, xyz);
			a_sub[4] = ComputeDetTex(x3, x5, x4, xyz);
			a_sub[5] = ComputeDetTex(x4, x5, x1, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] - a) / a) < tol)
			{
				return A->elementIndex;
			}
			break;

		case MshElemType::PRISM:
			a_sub[0] = ComputeDetTex(x1, x2, x3, xyz);
			a_sub[1] = ComputeDetTex(x4, x6, x5, xyz);
			a_sub[2] = ComputeDetTex(x1, x4, x2, xyz);
			a_sub[3] = ComputeDetTex(x2, x4, x5, xyz);
			a_sub[4] = ComputeDetTex(x2, x5, x3, xyz);
			a_sub[5] = ComputeDetTex(x3, x5, x6, xyz);
			a_sub[6] = ComputeDetTex(x3, x6, x1, xyz);
			a_sub[7] = ComputeDetTex(x1, x6, x4, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] + a_sub[6] + a_sub[7] - a) / a)
			    < tol)
			{
				return A->elementIndex;
			}
			break;

		case MshElemType::HEXAHEDRON:
			a_sub[0] = ComputeDetTex(x1, x2, x4, xyz);
			a_sub[1] = ComputeDetTex(x4, x2, x3, xyz);
			a_sub[2] = ComputeDetTex(x2, x6, x3, xyz);
			a_sub[3] = ComputeDetTex(x3, x6, x7, xyz);
			a_sub[4] = ComputeDetTex(x3, x7, x4, xyz);
			a_sub[5] = ComputeDetTex(x4, x7, x8, xyz);
			a_sub[6] = ComputeDetTex(x4, x8, x1, xyz);
			a_sub[7] = ComputeDetTex(x1, x8, x5, xyz);
			a_sub[8] = ComputeDetTex(x1, x5, x2, xyz);
			a_sub[9] = ComputeDetTex(x2, x5, x6, xyz);
			a_sub[10] = ComputeDetTex(x5, x8, x6, xyz);
			a_sub[11] = ComputeDetTex(x6, x8, x7, xyz);

			if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] + a_sub[6] + a_sub[7] + a_sub[8]
			          + a_sub[9]
			          + a_sub[10]
			          + a_sub[11]
			          - a)
			         / a)
			    < tol)
			{
				return A->elementIndex;
			}
			break;

		default:
			// do nothing with other elements
			break;
	}
	// the particle is not in this element
	return -1;
}

/**************************************************************************
   MSHLib-Method:
   Task:This function solves the next info of the particle by advection only
    1: The particle moved to the neighbor element
    0: Particle stays in the same element.
   -1: Particle displaced outside of the domain
   -2: The function failed
   Programing:
   09/2005 PCH Implementation
   03/2006 PCH Upgraded as one.
**************************************************************************/
int RandomWalk::SolveForNextPosition(Particle* A, Particle* B)
{
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* theElement = m_msh->ele_vector[B->elementIndex];

	// Getting the number of the edges in the element that Particle P belongs
	int nEdges = theElement->GetEdgesNumber();

	// TF - to avoid an abort while processing line elements
	// Is the handling for the mesh element type LINE not correct implemented?
	if (theElement->GetElementType() == MshElemType::LINE)
	{
		nEdges = 0;
	}

	int countNoIntersection = 0;

	// The estimated position advected for the given B->t
	double dD[3];
	dD[0] = dD[1] = dD[2] = 0.0;
	double Z[3];
	Z[0] = Z[1] = Z[2] = 0.0;
	int ele_dim = theElement->GetDimension();

	// Initialize some variables.
	double dtt = 0.0, dt1 = 0.0, dt2 = 0.0, d1 = 0.0, d = 0.0;
	double tolerance = 1e-6;
	double timeSplit = 100; // Important: This timeSplit is a bit sensitive.

	if (ele_dim < 3)
	{
		// Loop over the edges
		for (int i = 0; i < nEdges; ++i)
		{
			// Get the edges of the element
			vec<MeshLib::CEdge*> theEdges(nEdges);
			theElement->GetEdges(theEdges);

			// Get the nodes of the edge
			vec<MeshLib::CNode*> theNodes(3);
			theEdges[i]->GetNodes(theNodes);

			double p1[3], p2[3], p3[3], p4[3];
			// RWPT - IM
			// Two points in the edge
			double X1[3] = {theNodes[0]->getData()[0], theNodes[0]->getData()[1], theNodes[0]->getData()[2]};
			double X2[3] = {theNodes[1]->getData()[0], theNodes[1]->getData()[1], theNodes[1]->getData()[2]};
			//         X1[0] = theNodes[0]->X(); X1[1] = theNodes[0]->Y(); X1[2] = theNodes[0]->Z();
			//         X2[0] = theNodes[1]->X(); X2[1] = theNodes[1]->Y(); X2[2] = theNodes[1]->Z();
			ToTheXYPlane(theElement, X1);
			ToTheXYPlane(theElement, X2);
			for (int j = 0; j < 3; ++j)
			{
				p1[j] = X1[j];
				p2[j] = X2[j];
			}
			// The starting point displaced by pure advection
			p3[0] = B->x;
			p3[1] = B->y;
			p3[2] = B->z;
			ToTheXYPlane(theElement, p3);
			p3[2] = theElement->GetAngle(2);

			int dDStatus = 1;

			// If the mode is for heterogeneous media
			if (RWPTMode % 2 == 1)
				// This currently only return TRUE (1)
				dDStatus = SolveForDisplacementByDerivativeOfDispersion(A, dD);

			// Let's get the local vector for particle velocity
			double V[3];

			// Create random drift according to the element dimension
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
				RandomWalkDrift(Z, ele_dim);
			if (dDStatus == 1)
			{
				if (ele_dim == 2)
				{
					// RWPT - IM
					// This should be done carefully. Velocity should be transformed to be on the XY plane.
					// dD[] should be fine because it is handled in the SolveForDisplacementByDerivativeOfDispersion
					// function.
					// Z[] should also be fine. Just randome nubmers.
					// D[] Yes, this should be fine too. It is handled in the SolveDispersionCoefficient function.
					// OK. Just velocity left.
					// In fact, V[2] gotta be zero.
					V[0] = B->Vx;
					V[1] = B->Vy;
					V[2] = B->Vz;
					ToTheXYPlane(B->elementIndex, V);
					double dsp[3];
					GetDisplacement(B, Z, V, dD, B->t, dsp);

					// Fix for translation
					p4[0] = p3[0] + dsp[0];
					p4[1] = p3[1] + dsp[1];
					p4[2] = theElement->GetAngle(2);
				}
				else
				{
					printf("Other dimensions are not implemented yet.\n");
					abort();
				}
			}
			else
			{
				// This should never be the case by now. Later on, maybe.
				printf("SolveForDisplacementByDerivativeOfDispersion failed\n");
				abort();
			}

			// Initialize the values for getting the intersection
			dtt = B->t;
			double x = 0.0, y = 0.0, ra = 0.0, rb = 0.0;

			int status
			    = G_intersect_line_segments(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], &ra, &rb, &x, &y);

			// If intersection is a sinle point
			if (status == 1)
			{
				// Compute the time left over.
				double I[3];
				// Fix for translation
				I[0] = x;
				I[1] = y;
				I[2] = theElement->GetAngle(2);
				d1 = SolveDistanceBetweenTwoPoints(p3, I);
				d = SolveDistanceBetweenTwoPoints(p3, p4);
				dt1 = dtt * d1 / d;
				dt2 = dtt - dt1;

				// dt2 should be positive
				if (dt2 < 0.0)
					/*
					   printf("The program aborts because dt2 < 0.0\n");
					   abort();
					 */
					// I and P4 are almost identical.
					++countNoIntersection;
				else
				{
					double dsp[3];
					dsp[0] = dsp[1] = dsp[2] = 0.0;
					if (d1 > tolerance)
					{
						// Update the record.
						B->t = dt2;
						// Adjust the position for the obtained dt1.
						// But, keep in mind in this displacement there is no advective displament.
						double Vzero[3];
						Vzero[0] = Vzero[1] = Vzero[2] = 0.0;
						GetDisplacement(B, Z, Vzero, dD, dt1, dsp);

						double IC[3];
						// Fix for translation
						IC[0] = x + dsp[0];
						IC[1] = y + dsp[1];
						IC[2] = theElement->GetAngle(2);
						// Let's convert these XY plance coordinates to the real plane coordinates.
						ToTheRealPlane(B->elementIndex, IC);
						B->x = IC[0];
						B->y = IC[1];
						B->z = IC[2];
						CheckBoundary2D(A, B);
						return 1; // The element index switched to the neighbor element
					}
					else // It finds the wrong intersection.
					{
						// Just advance a little
						dt1 = dtt / timeSplit;
						dt2 = dtt - dt1;
						B->t = dt2;
						double IC[3];
						GetDisplacement(B, Z, V, dD, dt1, dsp);

						IC[0] = x + dsp[0];
						IC[1] = y + dsp[1];
						IC[2] = theElement->GetAngle(2);
						// Let's convert these XY plance coordinates to the real plane coordinates.
						ToTheRealPlane(B->elementIndex, IC);
						B->x = IC[0];
						B->y = IC[1];
						B->z = IC[2];
						CheckBoundary2D(A, B);
						return 1;
					}
				}
			}
			// It couldn't reach to the edge
			else if (status == 0)
				++countNoIntersection;
			// If two segments are parallel
			else if (status == -1)
				++countNoIntersection;
			// keep going.
			// If two segments are colinear
			else if (status == 2)
			{
				printf("The program aborts because two segments are colinear.\n");
				++countNoIntersection;
				// keep going.
			}
			else
			{
				printf("The program aborts because status of intersection search is not 1 or 0 or -1.\n");
				abort();
			}
		}
		// Check if the time left advances the particle within this element
		if (countNoIntersection == nEdges)
		{
			if (ele_dim == 2)
			{
				double V[3];
				// In fact, V[2] gotta be zero.
				V[0] = B->Vx;
				V[1] = B->Vy;
				V[2] = B->Vz;
				ToTheXYPlane(B->elementIndex, V); // V XY planed.

				double dsp[3];
				dsp[0] = dsp[1] = dsp[2] = 0.0;
				GetDisplacement(B, Z, V, dD, B->t, dsp);

				// Assigning the next postion of the particle. The index of element in this if condition
				// should be one of the connected planes randomly chosen.
				// Now just solve the real plane coordinates for the particle at the next position.
				double P[3];
				P[0] = B->x;
				P[1] = B->y;
				P[2] = B->z;
				ToTheXYPlane(B->elementIndex, P);
				P[0] += dsp[0];
				P[1] += dsp[1];
				P[2] = theElement->GetAngle(2);
				ToTheRealPlane(B->elementIndex, P);
				B->x = P[0];
				B->y = P[1];
				B->z = P[2];
				CheckBoundary2D(A, B);
			}
			else if (ele_dim == 1)
			{
				// Create random numbers according to dimension
				RandomWalkDrift(Z, ele_dim);
				double V[3];
				V[0] = B->Vx;
				V[1] = B->Vy;
				V[2] = B->Vz;

				double dsp[3];
				dsp[0] = dsp[1] = dsp[2] = 0.0;
				GetDisplacement(B, Z, V, dD, B->t, dsp);
				B->x += dsp[0];
				B->y += dsp[1];
				B->z += dsp[2];

				// OK411 double dt = B->t;
			}
			else
			{
				printf("Other dimensions are not implemented yet.\n");
				abort();
			}
			B->t = 0.0;

			return 0; // Particle stays in the same element.
		}
		else
			return -2; // The function failed
	}
	else // For 3D elements
	{
		// Currently 3D elements only work for dispersion in Homogeneous.
		// WW int dDStatus = 1;
		// If the mode is for heterogeneous media
		if (RWPTMode % 2 == 1)
			// This currently only return TRUE (1)
			// WW dDStatus = SolveForDisplacementByDerivativeOfDispersion(A, dD);
			SolveForDisplacementByDerivativeOfDispersion(A, dD);

		// Let's get the local vector for particle velocity
		double V[3];
		V[0] = B->Vx;
		V[1] = B->Vy;
		V[2] = B->Vz;

		// Create random drift according to the element dimension
		if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
		{
			RandomWalkDrift(Z, ele_dim);
			double dsp[3];
			GetDisplacement(B, Z, V, dD, B->t, dsp);
			B->x += dsp[0];
			B->y += dsp[1];
			B->z += dsp[2];
		}

		// Initialize the values for getting the intersection
		B->t = 0.0;
		CheckBoundary3D(A, B);
		return 1;
	}
	return -2;
}

/**************************************************************************
   MSHLib-Method:
   Task:The function solves four different types of displacement.
   Programing:
   03/2006 PCH Implementation
**************************************************************************/
void RandomWalk::GetDisplacement(Particle* B, double* Z, double* V, double* dD, double time, double* dsp)
{
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* theElement = m_msh->ele_vector[B->elementIndex];
	int ele_dim = theElement->GetDimension();
	double Dxx = 0.0, Dxy = 0.0, Dxz = 0.0, Dyx = 0.0, Dyy = 0.0, Dyz = 0.0, Dzx = 0.0, Dzy = 0.0, Dzz = 0.0;

	if (ele_dim == 2) // If 2D,
	{
		// If the mode is for heterogeneous media
		if (RWPTMode % 2 == 1)
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
			{
				Dxx = sqrt(6.0 * B->D[0] * time) * Z[0];
				Dxy = sqrt(6.0 * fabs(B->D[1]) * time) * Z[1];
				Dyx = sqrt(6.0 * fabs(B->D[3]) * time) * Z[0];
				Dyy = sqrt(6.0 * B->D[4] * time) * Z[1];

				dsp[0] = V[0] * time + dD[0] * time + Dxx + Dxy;
				dsp[1] = V[1] * time + dD[1] * time + Dyx + Dyy;
			}
			else // advection only
			{
				// Do nothing for dipsersive transport.
				dsp[0] = V[0] * time + dD[0] * time;
				dsp[1] = V[1] * time + dD[1] * time;
			}
		}
		else // Homogeneous case
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
			{
				// Homo and hetero in this case are the same.
				Dxx = sqrt(6.0 * B->D[0] * time) * Z[0];
				Dxy = sqrt(6.0 * fabs(B->D[1]) * time) * Z[1];
				Dyx = sqrt(6.0 * fabs(B->D[3]) * time) * Z[0];
				Dyy = sqrt(6.0 * B->D[4] * time) * Z[1];

				dsp[0] = V[0] * time + Dxx + Dxy;
				dsp[1] = V[1] * time + Dyx + Dyy;
			}
			else // advection only
			{
				dsp[0] = V[0] * time;
				dsp[1] = V[1] * time;
			}
		}
		// Fix for translation
		dsp[2] = theElement->GetAngle(2);
	}
	else if (ele_dim == 1) // If 1D,
	{
		double VV = 0.0, DD = 0.0, Dsp = 0.0;
		int coordinateflag = m_msh->GetCoordinateFlag();
		if (coordinateflag == 10) // x only
		{
			VV = V[0];
			DD = sqrt(6.0 * B->D[0] * time) * Z[0];
		}
		else if (coordinateflag == 11) // y only
		{
			VV = V[1];
			DD = sqrt(6.0 * B->D[4] * time) * Z[0];
		}
		else if (coordinateflag == 12) // z only
		{
			VV = V[2];
			DD = sqrt(6.0 * B->D[8] * time) * Z[0];
		}
		else // Something Wrong.
			abort();

		// If the mode is for heterogeneous media
		if (RWPTMode % 2 == 1)
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on

				Dsp = VV * time + dD[0] * time + DD;
			else // advection only

				// Do nothing for dipsersive transport.
				Dsp = VV * time + dD[0] * time;
		}
		else // Homogeneous case
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on

				// Homo and hetero in this case are the same.
				Dsp = VV * time + DD;
			else // advection only

				Dsp = VV * time;
		}

		if (coordinateflag == 10) // x only
			dsp[0] = Dsp;
		else if (coordinateflag == 11) // y only
			dsp[1] = Dsp;
		else if (coordinateflag == 12) // z only
			dsp[2] = Dsp;
		else // Something Wrong.
			abort();
	}
	else // 3D elements
	{
		// If the mode is for heterogeneous media
		if (RWPTMode % 2 == 1)
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
			{
				Dxx = sqrt(6.0 * B->D[0] * time) * Z[0];
				Dxy = sqrt(6.0 * fabs(B->D[1]) * time) * Z[1];
				Dxz = sqrt(6.0 * fabs(B->D[2]) * time) * Z[2];
				Dyx = sqrt(6.0 * fabs(B->D[3]) * time) * Z[0];
				Dyy = sqrt(6.0 * B->D[4] * time) * Z[1];
				Dyz = sqrt(6.0 * fabs(B->D[5]) * time) * Z[2];
				Dzx = sqrt(6.0 * fabs(B->D[6]) * time) * Z[0];
				Dzy = sqrt(6.0 * fabs(B->D[5]) * time) * Z[1];
				Dzz = sqrt(6.0 * B->D[8] * time) * Z[2];

				dsp[0] = V[0] * time + dD[0] * time + Dxx + Dxy + Dxz;
				dsp[1] = V[1] * time + dD[1] * time + Dyx + Dyy + Dxz;
				dsp[2] = V[2] * time + dD[2] * time + Dzx + Dzy + Dzz;
			}
			else // advection only
			{
				// Do nothing for dipsersive transport.
				dsp[0] = V[0] * time + dD[0] * time;
				dsp[1] = V[1] * time + dD[1] * time;
				dsp[2] = V[2] * time + dD[2] * time;
			}
		}
		else // Homogeneous case
		{
			if (RWPTMode < 2 || RWPTMode > 3) // whenever dispersion is on
			{
				// Homo and hetero in this case are the same.
				Dxx = sqrt(6.0 * B->D[0] * time) * Z[0];
				Dxy = sqrt(6.0 * fabs(B->D[1]) * time) * Z[1];
				Dxz = sqrt(6.0 * fabs(B->D[2]) * time) * Z[2];
				Dyx = sqrt(6.0 * fabs(B->D[3]) * time) * Z[0];
				Dyy = sqrt(6.0 * B->D[4] * time) * Z[1];
				Dyz = sqrt(6.0 * fabs(B->D[5]) * time) * Z[2];
				Dzx = sqrt(6.0 * fabs(B->D[6]) * time) * Z[0];
				Dzy = sqrt(6.0 * fabs(B->D[5]) * time) * Z[1];
				Dzz = sqrt(6.0 * B->D[8] * time) * Z[2];

				dsp[0] = V[0] * time + Dxx + Dxy + Dxz;
				dsp[1] = V[1] * time + Dyx + Dyy + Dyz;
				dsp[2] = V[2] * time + Dzx + Dzy + Dzz;
			}
			else // advection only
			{
				dsp[0] = V[0] * time;
				dsp[1] = V[1] * time;
				dsp[2] = V[2] * time;
			}
		}
	}
}

/**************************************************************************
   MSHLib-Method:
   Task:This function returns the index of the element that contains
    the particle by comparing the previous position with current position
   Programing:
   01/2007 PCH Implementation
**************************************************************************/
int RandomWalk::GetTheElementOfTheParticle(Particle* Pold, Particle* Pnew)
{
	int index = -10;

	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
//	for(int i=0; i< (int)pcs_vector.size(); ++i)
//	{
//		m_pcs = pcs_vector[i];
//
//		// Select the mesh whose process name has the mesh for Fluid_Momentum
//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
//			m_msh = FEMGet("RICHARDS_FLOW");
//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
//			m_msh = FEMGet("LIQUID_FLOW");
//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
//			m_msh = FEMGet("GROUNDWATER_FLOW");
//	}

#ifdef ALLOW_PARTICLES_GO_OUTSIDE
	if (Pold->elementIndex != -10)
	{
#endif

		MeshLib::CElem* theElement = m_msh->ele_vector[Pold->elementIndex];
		// Let's check this element first.
		index = IsTheParticleInThisElement(Pold);
		if (index != -1)
			return index;

		// First find the edge of the previous element
		// Mount the edges of the element
		int nnode = theElement->GetEdgesNumber();
		vec<MeshLib::CEdge*> theEdgesOfThisElement(nnode);
		theElement->GetEdges(theEdgesOfThisElement);
		// Mount the nodes of the edge
		vec<MeshLib::CNode*> theNodesOfThisEdge(3);

		for (int i = 0; i < nnode; ++i)
		{
			// Get the nodes of the edge
			theEdgesOfThisElement[i]->GetNodes(theNodesOfThisEdge);

			double p1[3], p2[3], p3[3], p4[3];
			// RWPT - IM
			// Two points in the edge
			double X1[3] = {theNodesOfThisEdge[0]->getData()[0], theNodesOfThisEdge[0]->getData()[1],
			                theNodesOfThisEdge[0]->getData()[2]};
			double X2[3] = {theNodesOfThisEdge[1]->getData()[0], theNodesOfThisEdge[1]->getData()[1],
			                theNodesOfThisEdge[1]->getData()[2]};
//         double X1[3], X2[3];
//         X1[0] = theNodesOfThisEdge[0]->X(); X1[1] = theNodesOfThisEdge[0]->Y(); X1[2] = theNodesOfThisEdge[0]->Z();
//         X2[0] = theNodesOfThisEdge[1]->X(); X2[1] = theNodesOfThisEdge[1]->Y(); X2[2] = theNodesOfThisEdge[1]->Z();
#ifdef TWODINTHREED
			ToTheXYPlane(theElement, X1);
			ToTheXYPlane(theElement, X2);
#endif
			for (int j = 0; j < 3; ++j)
			{
				p1[j] = X1[j];
				p2[j] = X2[j];
			}

			// The starting point which is the previous position
			p3[0] = Pold->x;
			p3[1] = Pold->y;
			p3[2] = Pold->z;
#ifdef TWODINTHREED
			// RWPT - IM
			ToTheXYPlane(theElement, p3);
#endif
			// The ending point which is the current position
			p4[0] = Pnew->x;
			p4[1] = Pnew->y;
			p4[2] = Pnew->z;

			double x = 0.0, y = 0.0, ra = 0.0, rb = 0.0;

			int status
			    = G_intersect_line_segments(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], &ra, &rb, &x, &y);

			if (status == 0) // Not intersect but extension intersects
				;
			else if (status == -1) // Parallel
				;
			else if (status == 1) // single intersection - Means the current position is outside of this element

				// Do further implementation here.
				index = theEdgesOfThisElement[i]->connected_elements[1];
			else if (status == 2) // Overlap just do nothing
				index = Pold->elementIndex; // This should indicate the particle in this element, then.
			else
			{
				printf("Not making any sense.\n");
				abort();
			}
		}

#ifdef ALLOW_PARTICLES_GO_OUTSIDE
	}
#endif

	return index;
}

/**************************************************************************
   MSHLib-Method:
   Task:This function returns the index of the element that contains
    the particle from neighboring elements only.
   Programing:
   12/2005 PCH Implementation
**************************************************************************/
int RandomWalk::GetTheElementOfTheParticleFromNeighbor(Particle* A)
{
	int index = -10;

	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
//	for(int i=0; i< (int)pcs_vector.size(); ++i)
//	{
//		m_pcs = pcs_vector[i];
//
//		// Select the mesh whose process name has the mesh for Fluid_Momentum
//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
//			m_msh = FEMGet("RICHARDS_FLOW");
//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
//			m_msh = FEMGet("LIQUID_FLOW");
//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
//			m_msh = FEMGet("GROUNDWATER_FLOW");
//	}

#ifdef ALLOW_PARTICLES_GO_OUTSIDE
	if (A->elementIndex != -10)
	{
#endif

		MeshLib::CElem* theElement = m_msh->ele_vector[A->elementIndex];
		// Let's check this element first.
		index = IsTheParticleInThisElement(A);
		if (index != -1)
		{
			A->elementIndex = index;
			return index;
		}

		// First meighbor's search around the main element
		for (size_t i = 0; i < theElement->GetFacesNumber(); ++i)
		{
			MeshLib::CElem* thisNeighbor = theElement->GetNeighbor(i);

			// If element type has more dimension than line.
			if (thisNeighbor->GetElementType() != MshElemType::LINE)
			{
				// If the particle belongs to this element
				A->elementIndex = thisNeighbor->GetIndex();
				index = IsTheParticleInThisElement(A);
				if (index != -1)
					return index;

				// Second, search the neighbor's neighbor
				for (size_t j = 0; j < thisNeighbor->GetFacesNumber(); ++j)
				{
					MeshLib::CElem* theNeighborsNeighbor = thisNeighbor->GetNeighbor(j);

					if (theNeighborsNeighbor->GetElementType() != MshElemType::LINE)
					{
						// If the particle belongs to this element
						A->elementIndex = theNeighborsNeighbor->GetIndex();
						index = IsTheParticleInThisElement(A);
						if (index != -1)
							return index;

						// Third, search the neighbor's neighbor's neighbor
						for (size_t k = 0; k < theNeighborsNeighbor->GetFacesNumber(); ++k)
						{
							MeshLib::CElem* theNeighborsNeighborsNeighbor = theNeighborsNeighbor->GetNeighbor(k);

							if (theNeighborsNeighborsNeighbor->GetElementType() != MshElemType::LINE)
							{
								// If the particle belongs to this element
								A->elementIndex = theNeighborsNeighborsNeighbor->GetIndex();
								index = IsTheParticleInThisElement(A);
								if (index != -1)
									return index;
							}
						}
					}
				}
			}
		}

		// If the code pases the following loop, it means I am not lucky in this neighbor search.
		int numberOfElements = (int)m_msh->ele_vector.size();
		for (int i = 0; i < numberOfElements; ++i)
		{
			MeshLib::CElem* thisElement = m_msh->ele_vector[i];

			if (thisElement->GetElementType() != MshElemType::LINE)
			{
				// If the particle belongs to this element
				A->elementIndex = thisElement->GetIndex();
				index = IsTheParticleInThisElement(A);
				if (index != -1)
					return index;
			}
		}

		// The search failed
		if (index == -1)
		{
			index = -10;
			printf("Searching the index from the neighbor failed\n");
			printf("The particle should be outside of the domain.\n");
		}

#ifdef ALLOW_PARTICLES_GO_OUTSIDE
	}
#endif
	return index;
}

/**************************************************************************
   MSHLib-Method:
   Task:This function returns the index of the element that contains
    the particle if the particle exists in the element.
    Or return -1 if the particle is not in the element.
   Programing:
   12/2005 PCH Implementation
   02/2006 PCH Improved for the RWPT method in Fracture Networks.
   02/2006 PCH The ray method implemented based on the proven theory.
**************************************************************************/
int RandomWalk::IsTheParticleInThisElement(Particle* A)
{
	// Mount the proper mesh
	m_msh = selectMeshForFluidMomentumProcess();
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* theElement = m_msh->ele_vector[A->elementIndex];

	int ele_dim = theElement->GetDimension();

	if (ele_dim == 2)
	{
		// Getting the number of the edges in the element that Particle P belongs
		int nEdges = theElement->GetEdgesNumber();
		int countOfInterception = 0;
		// WW int parallel = 0;
		// Loop over the edges
		for (int i = 0; i < nEdges; ++i)
		{
			// Get the edges of the element
			vec<MeshLib::CEdge*> theEdges(nEdges);
			theElement->GetEdges(theEdges);

			// Get the nodes of the edge
			vec<MeshLib::CNode*> theNodes(3);
			theEdges[i]->GetNodes(theNodes);

			double p1[3], p2[3], p3[3], p4[3];
			// RWPT - IM
			// Two points in the edge

			double X1[3] = {theNodes[0]->getData()[0], theNodes[0]->getData()[1], theNodes[0]->getData()[2]};
			double X2[3] = {theNodes[1]->getData()[0], theNodes[1]->getData()[1], theNodes[1]->getData()[2]};

			//         double X1[3], X2[3];
			//         X1[0] = theNodes[0]->X(); X1[1] = theNodes[0]->Y(); X1[2] = theNodes[0]->Z();
			//         X2[0] = theNodes[1]->X(); X2[1] = theNodes[1]->Y(); X2[2] = theNodes[1]->Z();
			ToTheXYPlane(theElement, X1);
			ToTheXYPlane(theElement, X2);
			for (int j = 0; j < 3; ++j)
			{
				p1[j] = X1[j];
				p2[j] = X2[j];
			}

			// The starting point which is the particle position
			p3[0] = A->x;
			p3[1] = A->y;
			p3[2] = A->z;
			// RWPT - IM
			ToTheXYPlane(theElement, p3);
			// Make p4 very long in x direction in the XY plane
			double big = 1e3;
			p4[0] = p3[0] + big;
			p4[1] = p3[1];
			p4[2] = p3[2];

			double x = 0.0, y = 0.0, ra = 0.0, rb = 0.0;

			int status
			    = G_intersect_line_segments(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], &ra, &rb, &x, &y);

			if (status == 0) // Not intersect but extension intersects
				;
			// WW else if(status == -1)                    // Parallel
			// WW   parallel = 1;
			else if (status == 1)
				++countOfInterception; // single intersection
			else if (status == 2) // Overlap just do nothing
				; // This should indicate the particle in this element, then.
			else
			{
				printf("Not making any sense.\n");
				abort();
			}
		}
		// Check if this particle is inside of the element
		// If the number of interceptions is odd,
		// then, it is inside of this element.
		if (countOfInterception % 2 == 1)
			return A->elementIndex;
		// if the number is even,
		// then, it is outside
		else
			return -1; // This element does not have the particle.
	}
	else if (ele_dim == 3)
	{
		// double volumePism = theElement->GetVolume();//  ComputeVolume();
		// theElement->nodes
	}
	else if (ele_dim == 1)
	{
		// Since this is 1D, I'll do exhaustive search.
		// Checking the coordinateflag for proper solution.
		int coordinateflag = m_msh->GetCoordinateFlag();
		// OK411??? long
		for (int i = 0; i < (int)m_msh->ele_vector.size(); ++i)
		{
			double const* const pnt1(m_msh->nod_vector[m_msh->ele_vector[i]->GetNodeIndex(0)]->getData());
			double const* const pnt2(m_msh->nod_vector[m_msh->ele_vector[i]->GetNodeIndex(1)]->getData());
			if (coordinateflag == 10) // x only
			{
				if ((A->x >= pnt1[0] && A->x < pnt2[0]) || (A->x >= pnt2[0] && A->x < pnt1[0]))
					return i;
			}
			else if (coordinateflag == 11) // y only
			{
				if ((A->y >= pnt1[1] && A->y < pnt2[1]) || (A->y >= pnt2[1] && A->y < pnt1[1]))
					return i;
			}
			else if (coordinateflag == 12) // z only
			{
				if ((A->z >= pnt1[2] && A->z < pnt2[2]) || (A->z >= pnt2[2] && A->z < pnt1[2]))
					return i;
			}
			else
				return -1; // Something is wrong.
		}

		return -10; // The particle is outside of domain
	}

	return -1;
}

/**************************************************************************
   MSHLib-Method:
   Task:This function solves the distance between two points
   Programing:
   09/2005 PCH Implementation
**************************************************************************/
double RandomWalk::SolveDistanceBetweenTwoPoints(double* p1, double* p2)
{
	double x = p2[0] - p1[0];
	double y = p2[1] - p1[1];
	double z = p2[2] - p1[2];

	return sqrt(x * x + y * y + z * z);
}

/**************************************************************
* find interesection between two lines defined by points on the lines
* line segment A is (ax1,ay1) to (ax2,ay2)
* line segment B is (bx1,by1) to (bx2,by2)
* returns
*   -1 segment A and B do not intersect (parallel without overlap)
*    0 segment A and B do not intersect but extensions do intersect
*    1 intersection is a single point
*    2 intersection is a line segment (colinear with overlap)
* x,y intersection point
* ra - ratio that the intersection divides A
* rb - ratio that the intersection divides B
*
*                              B2
*                              /
*                             /
*   r=p/(p+q) : A1---p-------*--q------A2
*                           /
*                          /
*                         B1
*
**************************************************************/

/**************************************************************
 *
 * A point P which lies on line defined by points A1=(x1,y1) and A2=(x2,y2)
 * is given by the equation r * (x2,y2) + (1-r) * (x1,y1).
 * if r is between 0 and 1, p lies between A1 and A2.
 *
 * Suppose points on line (A1, A2) has equation
 *     (x,y) = ra * (ax2,ay2) + (1-ra) * (ax1,ay1)
 * or for x and y separately
 *     x = ra * ax2 - ra * ax1 + ax1
 *     y = ra * ay2 - ra * ay1 + ay1
 * and the points on line (B1, B2) are represented by
 *     (x,y) = rb * (bx2,by2) + (1-rb) * (bx1,by1)
 * or for x and y separately
 *     x = rb * bx2 - rb * bx1 + bx1
 *     y = rb * by2 - rb * by1 + by1
 *
 * when the lines intersect, the point (x,y) has to
 * satisfy a system of 2 equations:
 *     ra * ax2 - ra * ax1 + ax1 = rb * bx2 - rb * bx1 + bx1
 *     ra * ay2 - ra * ay1 + ay1 = rb * by2 - rb * by1 + by1
 *
 * or
 *
 *     (ax2 - ax1) * ra - (bx2 - bx1) * rb = bx1 - ax1
 *     (ay2 - ay1) * ra - (by2 - by1) * rb = by1 - ay1
 *
 * by Cramer's method, one can solve this by computing 3
 * determinants of matrices:
 *
 *    M  = (ax2-ax1)  (bx1-bx2)
 *         (ay2-ay1)  (by1-by2)
 *
 *    M1 = (bx1-ax1)  (bx1-bx2)
 *         (by1-ay1)  (by1-by2)
 *
 *    M2 = (ax2-ax1)  (bx1-ax1)
 *         (ay2-ay1)  (by1-ay1)
 *
 * Which are exactly the determinants D, D2, D1 below:
 *
 *   D  ((ax2-ax1)*(by1-by2) - (ay2-ay1)*(bx1-bx2))
 *
 *   D1 ((bx1-ax1)*(by1-by2) - (by1-ay1)*(bx1-bx2))
 *
 *   D2 ((ax2-ax1)*(by1-ay1) - (ay2-ay1)*(bx1-ax1))
 ***********************************************************************/

int RandomWalk::G_intersect_line_segments(double ax1, double ay1, double ax2, double ay2, double bx1, double by1,
                                          double bx2, double by2, double* ra, double* rb, double* x, double* y)
{
	double D = ((ax2 - ax1) * (by1 - by2) - (ay2 - ay1) * (bx1 - bx2));
	double D1 = ((bx1 - ax1) * (by1 - by2) - (by1 - ay1) * (bx1 - bx2));
	double D2 = ((ax2 - ax1) * (by1 - ay1) - (ay2 - ay1) * (bx1 - ax1));

	double d;
	d = D;

	if (d) /* lines are not parallel */
	{
		*ra = D1 / d;
		*rb = D2 / d;

		*x = ax1 + (*ra) * (ax2 - ax1);
		*y = ay1 + (*ra) * (ay2 - ay1);
		return *ra >= 0.0 && *ra <= 1.0 && *rb >= 0.0 && *rb <= 1.0;
	}

	if (D1 || D2)
		return -1; /* lines are parallel, not colinear */

	if (ax1 > ax2)
	{
		SWAP(ax1, ax2)
	}
	if (bx1 > bx2)
	{
		SWAP(bx1, bx2)
	}
	if (ax1 > bx2)
		return -1;
	if (ax2 < bx1)
		return -1;

	/* there is overlap */
	if (ax1 == bx2)
	{
		*x = ax1;
		*y = ay1;
		return 1; /* at endpoints only */
	}
	if (ax2 == bx1)
	{
		*x = ax2;
		*y = ay2;
		return 1; /* at endpoints only */
	}

	return 2; /* colinear with overlap on an interval, not just a single point*/
}

/**************************************************************
* Task: find interesection between two lines defined
*		by points on the lines
*
* returns
*   -1 parallel
*    0 do not intersect so in this element
*    1 intersection is a single point and stored in pi
*    2 on the plane
**************************************************************/

int RandomWalk::G_intersect_line_segments_3D(double* pl1, double* pl2, double* pp1, double* pp2, double* pp3,
                                             double* pi)
{
	// Solve for norm of the plane by performing cross product to solve for the plane
	double p2p1[3], p3p1[3], normOfThePlane[3];

	CrossProduction(pp2, pp1, p2p1);
	CrossProduction(pp3, pp1, p3p1);
	CrossProduction(p2p1, p3p1, normOfThePlane);

	double a, b, c, d;
	a = normOfThePlane[0];
	b = normOfThePlane[1];
	c = normOfThePlane[2];
	d = -(a * pp1[0] + b * pp1[1] + c * pp1[2]); // Refer pp1 point for solving d of the plane

	// Now solution of intersection
	double u, denominator; // The line eqn: P=P1+u(P2-P1) where P1 and P2 is on the line
	denominator = a * (pl1[0] - pl2[0]) + b * (pl1[1] - pl2[1]) + c * (pl1[2] - pl2[2]);

	// Check if denominator is zero or not
	if (fabs(denominator) < 10e-8) // If this is zero
	{
		// The line is either parallel or on the plane.
		// Check if the line is on the plane first by sustituting P1 into the plane eqn.
		double onThePlane = a * pl1[0] + b * pl1[1] + c * pl1[2] + d;

		if (fabs(onThePlane) < 10e-8) // If the line is on the plane
			return -1;
		else // If the line is parallel
			return 2;
	}
	else // The line segment either is in the element or has an intersectional point.
	{
		u = (a * pl1[0] + b * pl1[1] + c * pl1[2] + d) / denominator;

		if (u > 0 && u < 1) // The line is in this element
			return 0;
		else // The line intersects the plane
		{
			// Solve for the intersection
			pi[0] = pl1[0] + u * (pl2[0] - pl1[0]);
			pi[1] = pl1[1] + u * (pl2[1] - pl1[1]);
			pi[2] = pl1[2] + u * (pl2[2] - pl1[2]);

			return 1;
		}
	}
}

/**************************************************************************
   Task: The function solves reference coordinates of quad elements for a given physical
      coordinates. For bilinear elements such as quadrilateral elements,
      there is an analytical solution of the inverse bilinear
      transformation

   01/2007 PCH
**************************************************************************/
void RandomWalk::IsoparametricMappingQuadfromPtoR(int index, double* R)
{
	// Get the mesh first
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	// Mount the element fromthe first particle from particles initially
	MeshLib::CElem* theEle = m_msh->ele_vector[index];
	double tolerance = 1e-8;

	// If a quad element,
	int nnode = theEle->GetEdgesNumber();

	if (nnode == 4)
	{
		// Get physical coordinates of four corner points
		double x[4], y[4];
		for (int i = 0; i < nnode; ++i)
		{
			double const* const coords(theEle->GetNode(i)->getData());
			x[i] = coords[0];
			y[i] = coords[1];
		}

		// Some coeff's for convenience
		double ax, ay, bx, by, cx, cy, dx, dy;
		double X = R[0];
		double Y = R[1];
		/*
		   ax = x[0]-x[1]+x[2]-x[3]; bx = x[0]-x[1];
		   cx = x[0]-x[3]; dx = X - x[0];
		   ay = y[0]-y[1]+y[2]-y[3]; by = y[0]-y[1];
		   cy = y[0]-y[3]; dy = Y - y[0];
		 */

		ax = -0.25 * (x[0] - x[1] + x[2] - x[3]);
		bx = -0.25 * (x[0] - x[1] - x[2] + x[3]);
		cx = -0.25 * (x[0] + x[1] - x[2] - x[3]);
		dx = -X + 0.25 * (x[0] + x[1] + x[2] + x[3]);
		ay = -0.25 * (y[0] - y[1] + y[2] - y[3]);
		by = -0.25 * (y[0] - y[1] - y[2] + y[3]);
		cy = -0.25 * (y[0] + y[1] - y[2] - y[3]);
		dy = -Y + 0.25 * (y[0] + y[1] + y[2] + y[3]);

		// Cases for solution
		if (fabs(ax) > tolerance) // ax is not zero CASE 1
		{
			double s, r, t;

			s = ay / ax * bx - by;
			r = ay / ax * cx - cy;
			t = ay / ax * dx - dy;

			if (fabs(s) > tolerance) // s is not zero
			{
				if (fabs(r) > tolerance) // r is not zero
				{
					double A, B, C;
					A = -r * ax;
					B = r * bx - t * ax - s * cx;
					C = t * bx - s * dx;

					double d, y1, y2;
					d = B * B - 4.0 * A * C;
					if (d >= 0.0) // Only for real roots
					{
						y1 = (-B + sqrt(d)) / (2.0 * A);
						y2 = (-B - sqrt(d)) / (2.0 * A);
						// Get the right y
						//		if(y1 >= 0.0 && y1 <= 1.0)
						if (y1 >= -1.0 && y1 <= 1.0)
						{
							R[1] = y1;
							R[0] = -r / s * y1 - t / s;
							R[2] = 0.0; // For now for 2D elemenets
						}
						//		else if(y2 >= 0.0 && y2 <= 1.0)
						else if (y2 >= -1.0 && y2 <= 1.0)
						{
							R[1] = y2;
							R[0] = -r / s * y2 - t / s;
							R[2] = 0.0; // For now for 2D elemenets
						}
						else
						{
							printf("Failed to solve reference position for the particle\n");
							abort(); // Failed find the solution.
						}
					}
					else // Case 2
					{
						// y undefined. This should not happen.
					}
				}
				else
				{
				}
			}
			else // Case 3
			{
				// x undefined. This should not happen
			}
		}
		else // if ax = 0 CASE 2
		{
			if (fabs(cx) > tolerance) // cx is not zero
			{
				// If ay is not zero and bx is not zero
				if ((fabs(ay) > tolerance) && (fabs(bx) > tolerance))
				{
					double A, B, C;
					A = -ay * bx;
					B = bx * cy - by * cx - ay * dx;
					C = cy * dx - cx * dy;

					double d, x1, x2;

					// If ay = 0 or bx = 0
					if (fabs(ay) < tolerance || fabs(bx) < tolerance)
					{
						R[0] = (dy * cx - cy * dx) / (cy * bx - dy * cx - ay * dx);
						R[1] = (-bx * R[0] - dx) / cx;
						R[2] = 0.0;
					}
					else // If ay is not zero and bx is not zero,
					{
						d = B * B - 4.0 * A * C;
						if (d >= 0.0) // Only for real roots
						{
							x1 = (-B + sqrt(d)) / (2.0 * A);
							x2 = (-B - sqrt(d)) / (2.0 * A);

							// Get the right y
							if (x1 >= -1.0 && x1 <= 1.0)
							//	if(x1 >= 0.0 && x1 <= 1.0)
							{
								R[0] = x1;
								R[1] = (-bx * R[0] - dx) / cx;
								R[2] = 0.0; // For now for 2D elemenets
							}
							//		else if(x2 >= 0.0 && x2 <= 1.0)
							else if (x2 >= -1.0 && x2 <= 1.0)
							{
								R[0] = x2;
								R[1] = (-bx * R[0] - dx) / cx;
								R[2] = 0.0; // For now for 2D elemenets
							}
							else
							{
								printf("Failed to solve reference position for the particle\n");
								abort(); // Failed find the solution.
							}
						}
					}
				}
			}
			else // If cx = 0, CASE 3
			{
				R[0] = -dx / bx;
				R[1] = (by * dx - dy * dx) / (ay * dx + cy * bx);
				R[2] = 0.0;
			}
		}
		/*
		   // Convert it to -1 to 1
		   R[0] = 2*R[0]-1.0;
		   R[1] = 2*R[1]-1.0;
		 */
	}
	// else	// the element is not quad.
	//{
	//}
}

/**************************************************************************
   Task: The function solves physical coordinates of quad elements for a given reference
      coordinates. For bilinear elements such as quadrilateral elements,
      there is an analytical solution of the bilinear
      transformation

   01/2007 PCH
**************************************************************************/
void RandomWalk::IsoparametricMappingQuadfromRtoP(int index, double* P)
{
	// Get the mesh first
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	// Mount the element fromthe first particle from particles initially
	MeshLib::CElem* theEle = m_msh->ele_vector[index];
	// OK411 double tolerance = 1e-8;

	// If a quad element,
	int nnode = theEle->GetEdgesNumber();

	if (nnode == 4)
	{
		// Get physical coordinates of four corner points
		double x[4], y[4];
		for (int i = 0; i < nnode; ++i)
		{
			double const* const coords(theEle->GetNode(i)->getData());
			x[i] = coords[0];
			y[i] = coords[1];
		}

		double phat[3];
		phat[0] = P[0];
		phat[1] = P[1];
		phat[2] = P[2];

		P[0] = 0.25 * (x[0] * (1.0 - phat[0]) * (1.0 - phat[1]) + x[1] * (1.0 + phat[0]) * (1.0 - phat[1])
		               + x[2] * (1.0 + phat[0]) * (1.0 + phat[1])
		               + x[3] * (1.0 - phat[0]) * (1.0 + phat[1]));
		P[1] = 0.25 * (y[0] * (1.0 - phat[0]) * (1.0 - phat[1]) + y[1] * (1.0 + phat[0]) * (1.0 - phat[1])
		               + y[2] * (1.0 + phat[0]) * (1.0 + phat[1])
		               + y[3] * (1.0 - phat[0]) * (1.0 + phat[1]));
		P[2] = 0.0;
	}
	// else	// the element is not quad.
	//{
	//}
}

/**************************************************************************
   FEMLib-Method:
   Task: DoJointEffectOfElementInitially(void)
   Programing: This function does make a choice for each particle
         that lies on a crossroad or a joint. The contribution is
         determined by Fluid Momentum. Roulette Wheel Selection (RWE)
         determines which brach the particle continue to travel.
   02/2006 PCH
**************************************************************************/
void RandomWalk::DoJointEffectOfElementInitially(void)
{
	// Get the mesh first
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	// Looping all over the particles to have a choice which plane to go.
	// Because all of the particles are on the joint initially.
	for (int p = 0; p < m_msh->PT->numOfParticles; ++p)
	{
		// Mount the element fromthe first particle from particles initially
		int eleIdx = m_msh->PT->X[p].Now.elementIndex;
		MeshLib::CElem* theEle = m_msh->ele_vector[eleIdx];
		// Let's get the number of edges in the element and mount them
		int numOfEdgeIntheElement = theEle->GetEdgesNumber();
		vec<MeshLib::CEdge*> theEdges(numOfEdgeIntheElement);
		theEle->GetEdges(theEdges);
		MeshLib::CEdge* theJoint = NULL;

		// Now, 1. find the joint out of theses edges
		for (int i = 0; i < numOfEdgeIntheElement; ++i)
			// Is this a joint?
			if (theEdges[i]->GetJoint() == 1)
				theJoint = theEdges[i];

		// 2. Get multiple planes out of the joint
		// Now we need one of crossroad from the joint
		// Get the nodes of the edge
		vec<MeshLib::CNode*> theNodes(3);
		theJoint->GetNodes(theNodes);
		// I will use the first node of the joint as a crossroad
		MeshLib::CNode* crossnode = theNodes[0];
		// Let's mount the crossroad class
		CrossRoad* crossroad = NULL;
		for (int i = 0; i < (int)(m_msh->fm_pcs->crossroads.size()); ++i)
			if ((size_t)m_msh->fm_pcs->crossroads[i]->Index == crossnode->GetIndex())
				crossroad = m_msh->fm_pcs->crossroads[i];
		// Let's get the contribution of each connected plane.
		double chances[100]; // I just set 100 as a maximum number of
		// connected planes.
		for (int i = 0; i < crossroad->numOfThePlanes; ++i)
			chances[i] = crossroad->plane[i].ratio;

		// 3. Roulette Wheel Selection
		int whichWay = RouletteWheelSelection(chances, crossroad->numOfThePlanes);
		m_msh->PT->X[p].Now.elementIndex = m_msh->PT->X[p].Past.elementIndex = crossroad->plane[whichWay].eleIndex;
	}
}

/**************************************************************************
   Task: ToTheXYPlane(MeshLib::CElem* E, double* X)
   Programing: This function rotate-transforms the vector to be on the xy plane
   02/2006 PCH
**************************************************************************/
void RandomWalk::ToTheXYPlane(MeshLib::CElem* E, double* X)
{
	if (E->GetDimension() == m_msh->GetMaxElementDim())
		return;
	double x[3], xx[3];

	// Get the norm of the element plane and do some initialization
	for (int k = 0; k < 3; ++k)
		x[k] = xx[k] = 0.0;

	double alpha = E->GetAngle(0);
	double beta = E->GetAngle(1);
	// Let's rotate the original Enorm to the BB coordinate system
	// along the y axis
	x[0] = cos(alpha) * X[0] + sin(alpha) * X[2];
	x[1] = X[1];
	x[2] = -sin(alpha) * X[0] + cos(alpha) * X[2];
	// Let's rotate the BB coordinate system to the BBB coordinate system
	// along the x axis
	xx[0] = x[0];
	xx[1] = cos(beta) * x[1] - sin(beta) * x[2];
	xx[2] = sin(beta) * x[1] + cos(beta) * x[2];

	for (int i = 0; i < 3; ++i)
		X[i] = xx[i];
	// Do translation along z'' axis.
	//	X[2] -= E->GetAngle(2);
}

void RandomWalk::ToTheXYPlane(int idx, double* X)
{
	CFEMesh* m_msh = NULL;
	if (fem_msh_vector.size() == 0)
		return; // OK
	m_msh = selectMeshForFluidMomentumProcess();
	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* E = m_msh->ele_vector[idx];

	if (E->GetDimension() == m_msh->GetMaxElementDim())
		return;
	double x[3], xx[3];

	// Get the norm of the element plane and do some initialization
	for (int k = 0; k < 3; ++k)
		x[k] = xx[k] = 0.0;

	double alpha = E->GetAngle(0);
	double beta = E->GetAngle(1);
	// Let's rotate the original Enorm to the BB coordinate system
	// along the y axis
	x[0] = cos(alpha) * X[0] + sin(alpha) * X[2];
	x[1] = X[1];
	x[2] = -sin(alpha) * X[0] + cos(alpha) * X[2];
	// Let's rotate the BB coordinate system to the BBB coordinate system
	// along the x axis
	xx[0] = x[0];
	xx[1] = cos(beta) * x[1] - sin(beta) * x[2];
	xx[2] = sin(beta) * x[1] + cos(beta) * x[2];

	for (int i = 0; i < 3; ++i)
		X[i] = xx[i];
	// Do translation along z'' axis.
	//	X[2] -= E->GetAngle(2);
}

/**************************************************************************
   Task: ToTheRealPlane(MeshLib::CElem* E, double* X)
   Programing: This function transform the vector on the xy plane to the
         original plane of the element in 3D.
   02/2006 PCH
**************************************************************************/
void RandomWalk::ToTheRealPlane(MeshLib::CElem* E, double* X)
{
	if (E->GetDimension() == m_msh->GetMaxElementDim())
		return;
	double x[3], xx[3];

	// Get the norm of the element plane and do some initialization
	for (int k = 0; k < 3; ++k)
		x[k] = xx[k] = 0.0;

	double alpha = E->GetAngle(0);
	double beta = E->GetAngle(1);
	// Let's rotate the original Enorm to the BB coordinate system
	// along the y axis
	x[0] = cos(alpha) * X[0] - sin(alpha) * X[2];
	x[1] = X[1];
	x[2] = sin(alpha) * X[0] + cos(alpha) * X[2];
	// Let's rotate the BB coordinate system to the BBB coordinate system
	// along the x axis
	xx[0] = x[0];
	xx[1] = cos(beta) * x[1] + sin(beta) * x[2];
	xx[2] = -sin(beta) * x[1] + cos(beta) * x[2];

	for (int i = 0; i < 3; ++i)
		X[i] = xx[i];
	// Let's translate back to z axis.
	//	X[2] += E->GetAngle(2);
}

void RandomWalk::ToTheRealPlane(int idx, double* X)
{
	CFEMesh* m_msh = NULL;
	if (fem_msh_vector.size() == 0)
		return;
	m_msh = selectMeshForFluidMomentumProcess();
	//	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	MeshLib::CElem* E = m_msh->ele_vector[idx];

	if (E->GetDimension() == m_msh->GetMaxElementDim())
		return;
	double x[3], xx[3];

	// Get the norm of the element plane and do some initialization
	for (int k = 0; k < 3; ++k)
		x[k] = xx[k] = 0.0;

	double alpha = E->GetAngle(0);
	double beta = E->GetAngle(1);
	// Let's rotate the original Enorm to the BB coordinate system
	// along the y axis
	x[0] = cos(alpha) * X[0] - sin(alpha) * X[2];
	x[1] = X[1];
	x[2] = sin(alpha) * X[0] + cos(alpha) * X[2];
	// Let's rotate the BB coordinate system to the BBB coordinate system
	// along the x axis
	xx[0] = x[0];
	xx[1] = cos(beta) * x[1] + sin(beta) * x[2];
	xx[2] = -sin(beta) * x[1] + cos(beta) * x[2];

	for (int i = 0; i < 3; ++i)
		X[i] = xx[i];
	// Let's translate back to z axis.
	//	X[2] += E->GetAngle(2);
}

/**************************************************************************
   Task: SolveAnglesOfTheElment(MeshLib::CElem* E)
   Programing: This function solves two angles for rotation transformation
   02/2006 PCH
**************************************************************************/
void RandomWalk::SolveAnglesOfTheElment(MeshLib::CElem* E)
{
	CFEMesh* m_msh = NULL;
	if (fem_msh_vector.size() == 0)
		return;
	m_msh = selectMeshForFluidMomentumProcess();
	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}

	double tolerance = 1e-20, Enorm[3];
	// Allocate angle memory dynamically.
	// E->AllocateMeomoryforAngle();

	// Get the norm of the element plane and do some initialization
	int coordinate_system = m_msh->GetCoordinateFlag();
	// If the coordinate system is xz plane or xyz, solve the angles.
	if (coordinate_system != 32 && coordinate_system != 22)
	{
		Enorm[0] = 0.0;
		Enorm[1] = 0.0;
		Enorm[2] = 1.0;
	}
	else
		for (int k = 0; k < 3; ++k)
			Enorm[k] = E->getTransformTensor(k + 6);

	// solve for two angles for two rotation transformation.
	// Solving alpha that will be used for rotation along y axis.
	double alpha = 0.0;

	if (Enorm[0] * Enorm[0] + Enorm[2] * Enorm[2] < tolerance)
		; // have alpha to be zero. No need to rotate.
	else
		alpha = acos(Enorm[2] / sqrt(Enorm[0] * Enorm[0] + Enorm[2] * Enorm[2]));
	// The following if condition is required because
	// the acos function is not distintive in the case that Enorm[0]'s of
	// the two planes are opposite each other.
	if (Enorm[0] < 0.0)
		E->SetAngle(0, alpha);
	else
		E->SetAngle(0, alpha + 2.0 * (PI - alpha));

	// Solving beta that will be used for rotation along x' axis
	double beta = 0.0, BB[3], TranZ;
	// Let's rotate the original Enorm to this coordinate system.
	BB[0] = cos(E->GetAngle(0)) * Enorm[0] + sin(E->GetAngle(0)) * Enorm[2];
	BB[1] = Enorm[1];
	BB[2] = -sin(E->GetAngle(0)) * Enorm[0] + cos(E->GetAngle(0)) * Enorm[2];
	if (BB[2] > tolerance)
		beta = atan(BB[1] / BB[2]);
	else // if BB[2] is zero
		beta = 0.5 * PI;

	E->SetAngle(1, beta);

	// Solve for the translation.
	// I'll use the center of the element for this translation.
	double const* center = E->GetGravityCenter();
	double x[3], xx[3];
	// Get the norm of the element plane and do some initialization
	for (int k = 0; k < 3; ++k)
		x[k] = xx[k] = 0.0;
	// Let's rotate the original Enorm to the BB coordinate system
	// along the y axis
	x[0] = cos(E->GetAngle(0)) * center[0] + sin(E->GetAngle(0)) * center[2];
	x[1] = center[1];
	x[2] = -sin(E->GetAngle(0)) * center[0] + cos(E->GetAngle(0)) * center[2];
	// Let's rotate the BB coordinate system to the BBB coordinate system
	// along the x axis
	xx[0] = x[0];
	xx[1] = cos(E->GetAngle(1)) * x[1] - sin(E->GetAngle(1)) * x[2];
	xx[2] = sin(E->GetAngle(1)) * x[1] + cos(E->GetAngle(1)) * x[2];
	TranZ = xx[2];
	E->SetAngle(2, TranZ);
}

/**************************************************************************
   FEMLib-Method:
   Task: RouletteWheelSelection(double *chances, int numOfCases)
   Programing: This function makes a choice by RWS that is based on
         velocity contribution on each of the connected planes.
   02/2006 PCH
**************************************************************************/
int RandomWalk::RouletteWheelSelection(double* chances, int numOfCases)
{
	int whichOne = -1000; // Set it meaningless
	double* roulette;
	roulette = new double[numOfCases]();

	MakeRoulette(chances, roulette, numOfCases);
	whichOne = Select(roulette, numOfCases);

	delete[] roulette;

	return whichOne;
}

/**************************************************************************
   FEMLib-Method:
   Task: MakeRoulette(double* fit, double* roulette, int numOfCases)
   Programing: This function makes a roulette according to chances (fit)
   02/2006 PCH
**************************************************************************/
void RandomWalk::MakeRoulette(double* fit, double* roulette, int numOfCases)
{
	double* pi = NULL;
	double* fitProbability = NULL;
	double fitTotal = 0.0, ProbTotal = 0.0;

	// Create memory for these two arrays dynamically
	pi = new double[numOfCases]();
	fitProbability = new double[numOfCases]();

	for (int i = 0; i < numOfCases; ++i)
	{
		// Function modification can be done here.
		pi[i] = 1. / exp(-fit[i]);
		fitTotal += pi[i];
	}

	// Making Roulette
	for (int i = 0; i < numOfCases; ++i)
	{
		fitProbability[i] = pi[i] / fitTotal;
		ProbTotal += fitProbability[i];
		roulette[i] = ProbTotal;
	}

	delete[] pi;
	delete[] fitProbability;
}

/**************************************************************************
   FEMLib-Method:
   Task: Select(double* roulette)
   Programing: This function select a choice out of the roulette
   02/2006 PCH
**************************************************************************/
int RandomWalk::Select(double* roulette, int numOfCases)
{
	double probability;

	probability = randomZeroToOne();
	for (int i = 0; i < numOfCases; ++i)
		if (probability < roulette[i])
			return i;
	return 0;
}

/**************************************************************************
   FEMLib-Method:
   Task: ReadInVelocityFieldOnNodes(string file_base_name)
   Programming: This function gets velocity fields from a separate file.
         A COMPLETE bypass of the FEM.
   05/2006 PCH
**************************************************************************/
int RandomWalk::ReadInVelocityFieldOnNodes(string file_base_name)
{
	// Something must be done later on here.
	CFEMesh* m_msh(selectMeshForFluidMomentumProcess());
	// Mount the proper mesh
	//	for(int i=0; i< (int)pcs_vector.size(); ++i)
	//	{
	//		m_pcs = pcs_vector[i];
	//
	//		// Select the mesh whose process name has the mesh for Fluid_Momentum
	//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
	//			m_msh = FEMGet("RICHARDS_FLOW");
	//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
	//			m_msh = FEMGet("LIQUID_FLOW");
	//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
	//			m_msh = FEMGet("GROUNDWATER_FLOW");
	//	}
	CRFProcess* m_pcs = PCSGet("FLUID_MOMENTUM");

	// File handling
	string vel_file_name;
	ios::pos_type position;
	vel_file_name = file_base_name + ".vel";

	ifstream vel_file(vel_file_name.data(), ios::in);

	int End = 1;
	string strbuffer;

	while (End)
	{
		// OK411??? long
		for (int i = 0; i < (int)m_msh->nod_vector.size(); ++i)
		{
			double v[3];
			for (int p = 0; p < 3; ++p)
				v[p] = 0.0;
			vel_file >> v[0] >> v[1] >> v[2] >> ws;

			// Let's assign the velocity
			for (int j = 0; j < 3; ++j)
			{
				int nidx1 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[j]) + 1;
				m_pcs->SetNodeValue(m_msh->Eqs2Global_NodeIndex[i], nidx1, v[j]);
			}
		}
		End = 0;
	}

	return 1;
}

/**************************************************************************
   FEMLib-Method:
   Programing: This function buildS FDM index for quads
   01/2007 PCH
   Modified:: 03/2010 JT ... re-designed and expanded to 3-D
**************************************************************************/
void RandomWalk::buildFDMIndex(void)
{
	double xmax, ymax, zmax;
	double xmin, ymin, zmin;
	long i, j, k, iel, ic, jc, kc, nels;
	// WW long ne, nels;
	int index;
	neFDM = -1;

	// get mesh
	CFEMesh* m_msh(NULL);
	for (index = 0; index < (int)pcs_vector.size(); index++)
	{
		m_pcs = pcs_vector[index];
		//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos){
		if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
		{
			m_msh = FEMGet("RICHARDS_FLOW");
			break;
		}
		//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos){
		else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
		{
			m_msh = FEMGet("LIQUID_FLOW");
			break;
		}
		//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos){
		else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		{
			m_msh = FEMGet("GROUNDWATER_FLOW");
			break;
		}
	}
	selectMeshForFluidMomentumProcess();

	// total geometry size
	GEOCalcPointMinMaxCoordinates();
	xrw_range = pnt_x_max - pnt_x_min;
	yrw_range = pnt_y_max - pnt_y_min;
	zrw_range = pnt_z_max - pnt_z_min;

	// element size
	xmax = ymax = zmax = -1e+12;
	xmin = ymin = zmin = 1e+12;
	for (index = 0; index < 4; index++) // 4 nodes because FDM method only for equa-sized quads
	{
		//      x = m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(index)]->X();
		//      y = m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(index)]->Y();
		//      z = m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(index)]->Z();
		double const* const xyz(m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(index)]->getData());
		if (xyz[0] > xmax)
			xmax = xyz[0];
		if (xyz[0] < xmin)
			xmin = xyz[0];
		if (xyz[1] > ymax)
			ymax = xyz[1];
		if (xyz[1] < ymin)
			ymin = xyz[1];
		if (xyz[2] > zmax)
			zmax = xyz[2];
		if (xyz[2] < zmin)
			zmin = xyz[2];
	}
	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;
	if (dx < 1.e-12)
		dx = 1.e-12;
	if (dy < 1.e-12)
		dy = 1.e-12;
	if (dz < 1.e-12)
		dz = 1.e-12;
	nx = (int)floor(xrw_range / dx) + 1;
	ny = (int)floor(yrw_range / dy) + 1;
	nz = (int)floor(zrw_range / dz) + 1;
	// WW ne = nx*ny*nz;
	nels = m_msh->ele_vector.size();

	for (k = 0; k < nz; k++) // loop over the dummy element set
	{
		for (j = 0; j < ny; j++)
			for (i = 0; i < nx; i++)
			{
				FDMIndex one; // store dummy index by default, initialize class vector
				one.i = i;
				one.j = j;
				one.k = k;
				one.eleIndex = -5; // eleIndex -5 is dummy index different from -10
				for (iel = 0; iel < nels; iel++) // loop over mesh elements, assign them to dummy elements
				{
					double const* center = m_msh->ele_vector[iel]->GetGravityCenter();
					ic = (int)floor((center[0] - pnt_x_min) / dx);
					jc = (int)floor((center[1] - pnt_y_min) / dy);
					kc = (int)floor((center[2] - pnt_z_min) / dz);
					if (ic != i || jc != j || kc != k)
						continue;
					one.eleIndex = iel;
					break;
				}
				indexFDM.push_back(one);
				neFDM += 1;
			}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: ReadInVelocityFieldOnNodes(string file_base_name)
   Programing: This function build FDM index for quads
   01/2007 PCH
**************************************************************************
   void RandomWalk::buildFDMIndexOLDDDDDDDDDDDDDDDDDD(void)
   {
   double* Cx; double* Cy;

   Cx= new double[gli_points_vector.size()];
   Cy= new double[gli_points_vector.size()];

   // Solve for four corners from .gli
   for(int i=0; i<(int)gli_points_vector.size(); ++i) //OK411??? long
   {
   Cx[i] = gli_points_vector[i]->x;
   Cy[i] = gli_points_vector[i]->y;
   }

   // Solve for min and max of x and y respectively
   double xmax=-1e+12, ymax=-1e+12;
   double xmin=1e+12, ymin=1e+12;
   double minX;
   double minY;
   for(int i=0;i<(int)gli_points_vector.size(); ++i)  //OK411??? long
   {
   if(Cx[i] > xmax)
   xmax = Cx[i];
   if(Cy[i] > ymax)
   ymax = Cy[i];
   if(Cx[i] < xmin)
   xmin = Cx[i];
   if(Cy[i] < ymin)
   ymin = Cy[i];

   }
   minX = xmin; minY=ymin;
   XT = xmax - xmin; YT = ymax - ymin;

   // solve for dx and dy
   CFEMesh* m_msh = NULL;
   // Mount the proper mesh
   for(int i=0; i< (int)pcs_vector.size(); ++i)
   {
   m_pcs = pcs_vector[i];

   // Select the mesh whose process name has the mesh for Fluid_Momentum
   if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
   m_msh = FEMGet("RICHARDS_FLOW");
   else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
   m_msh = FEMGet("LIQUID_FLOW");
   else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
   m_msh = FEMGet("GROUNDWATER_FLOW");
   else;
   }
   // Solve for four corners from element #0.
   for(int i=0; i<4; ++i)
   {
   Cx[i] = m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(i)]->X();
   Cy[i] = m_msh->nod_vector[m_msh->ele_vector[0]->GetNodeIndex(i)]->Y();
   }

   // Use the same variables again
   xmax=-1e+12, ymax=-1e+12;
   xmin=1e+12, ymin=1e+12;
   for(int i=0;i<4; ++i)
   {
   if(Cx[i] > xmax)
   xmax = Cx[i];
   if(Cy[i] > ymax)
   ymax = Cy[i];
   if(Cx[i] < xmin)
   xmin = Cx[i];
   if(Cy[i] < ymin)
   ymin = Cy[i];

   }
   dx = xmax - xmin; dy = ymax - ymin;

   nx = (int)((XT+1e-12)/dx);
   ny = (int)((YT+1e-12)/dy);
   // We can allocate memory for indexFDM now.

   for(double j=0; j<YT; j += dx)
   {
   double seg_startY = 0.0, seg_endY = 0.0;
   seg_startY = minY+j;
   seg_endY = minY+j+dy;
   for(double i=0; i< XT; i += dy)
   {
   double seg_startX = 0.0, seg_endX = 0.0;
   seg_startX = minX+i;
   seg_endX = minX+i+dx;

   // Store the dummy index by default
   FDMIndex one;
   // eleIndex -5 is dummy index different from -10
   one.i = (int)(i/dx); one.j = (int)(j/dy); one.eleIndex = -5;

   for(int k=0; k<(int)m_msh->ele_vector.size(); ++k) //OK411??? long
   {
   double* center = m_msh->ele_vector[k]->GetGravityCenter();
   if( (center[0]>= seg_startX)	&&      (center[0]<= seg_endX)	&&
   (center[1]>= seg_startY)	&&      (center[1]<= seg_endY))
   {
   one.eleIndex = k;
   }
   }
   indexFDM.push_back(one);
   }
   }

   delete [] Cx; delete [] Cy;
   }

**************************************************************************/

/**************************************************************************
   FEMLib-Method:
   Task: Random Walk read function
   Programming:
   09/2005 PCH Destruct before read
**************************************************************************/
void PCTRead(string file_base_name)
{
	CFEMesh* m_msh = NULL;

	if (fem_msh_vector.size() == 0)
		return; // OK

	double xyz[3];
	// Mount the proper mesh
	m_msh = fem_msh_vector[0]; // Something must be done later on here.

	// File handling
	string pct_file_name;
	ios::pos_type position;
	pct_file_name = file_base_name + PCT_FILE_EXTENSION;

	ifstream pct_file(pct_file_name.data(), ios::in);
	if (!pct_file)
		return;

	std::stringstream ss;
	int srand_seed;
	string s_flag;
	getline(pct_file, s_flag);
	ss.str(s_flag);
	ss >> srand_seed;
	ss.clear();

	int End = 1;
	string strbuffer;
	RandomWalk* RW = NULL;
	m_msh->PT = new RandomWalk(srand_seed); // PCH
	RW = m_msh->PT;

	// Create pathline
	RandomWalk::Pathline path;
	while (End)
	{
		// Later on from this line, I can put which mesh I am dealing with.
		// pct_file >> RW->numOfParticles >> ws;

		getline(pct_file, s_flag);
		ss.str(s_flag);
		ss >> RW->numOfParticles;
		ss.clear();
		Trace one;

		RW->ChanceOfIrreversed = new double[RW->numOfParticles]; // YS

		int counter = 0;
		for (int i = 0; i < RW->numOfParticles; ++i)
		{
			// Assign the number to the particle
			int idx = 0, identity = 0;
			double x = 0.0, y = 0.0, z = 0.0, Starting = 0, vx = 0.0, vy = 0.0, vz = 0.0, K = 0.0;

			pct_file >> idx >> x >> y >> z >> identity >> Starting >> vx >> vy >> vz >> K >> ws;

			xyz[0] = x;
			xyz[1] = y;
			xyz[2] = z;

			idx = m_msh->FindElementByPoint(xyz);

			if (idx == -1)
				continue;

			one.Past.elementIndex = one.Now.elementIndex = idx;
			one.Past.x = one.Now.x = x;
			one.Past.y = one.Now.y = y;
			one.Past.z = one.Now.z = z;
			// JT 2010
			one.Past.StartingTime = one.Now.StartingTime = Starting;
			one.Past.identity = one.Now.identity = identity;
			one.Past.Vx = one.Now.Vx = vx;
			one.Past.Vy = one.Now.Vy = vy;
			one.Past.Vz = one.Now.Vz = vz;
			one.Past.K = one.Now.K = K;

			RW->X.push_back(one);

			RW->ChanceOfIrreversed[counter] = RW->randomZeroToOne();
			counter++;

			// Creat pathline
			if (i < 50)
				RW->pathline.push_back(path);
		}
		RW->numOfParticles = counter;
		End = 0;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: DATWriteFile
   Task: Write PCT file
   Programing:
   09/2005   PCH   Implementation
**************************************************************************/
void DATWriteParticleFile(int current_time_step)
{
	CFEMesh* m_msh = NULL;
	RandomWalk* RW = NULL;

	// Gather the momentum mesh
	size_t pcs_vector_size(pcs_vector.size());
	for (size_t i = 0; i < pcs_vector_size; ++i)
	{
		//		m_pcs = pcs_vector[i];
		const FiniteElement::ProcessType pcs_type(pcs_vector[i]->getProcessType());
		//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos){
		if (pcs_type == FiniteElement::RICHARDS_FLOW)
		{
			m_msh = FEMGet("RICHARDS_FLOW");
			break;
		}
		//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos){
		else if (pcs_type == FiniteElement::LIQUID_FLOW)
		{
			m_msh = FEMGet("LIQUID_FLOW");
			break;
		}
		else if (pcs_type == FiniteElement::GROUNDWATER_FLOW)
		{
			m_msh = FEMGet("GROUNDWATER_FLOW");
			break;
		}
	}

	RW = m_msh->PT;
	int np = RW->numOfParticles;

	// file naming
	char now[10];
	sprintf(now, "%i", current_time_step);
	string nowstr = now;

	string vtk_file_name =  pathJoin(defaultOutputPath, FileName);
	vtk_file_name += "_RWPT_" + nowstr + ".particles.vtk";
	fstream vtk_file (vtk_file_name.data(),ios::out);
	vtk_file.setf(ios::scientific,ios::floatfield);
	vtk_file.precision(12);
	if (!vtk_file.good())
		return;
	vtk_file.seekg(0L, ios::beg);

	// Write Header
	vtk_file << "# vtk DataFile Version 3.6.2"
	         << "\n";
	vtk_file << "Particle file: OpenGeoSys->Paraview. Current time (s) = " << RW->CurrentTime << "\n";
	vtk_file << "ASCII"
	         << "\n";
	vtk_file << "\n";
	vtk_file << "DATASET POLYDATA"
	         << "\n"; // KR vtk_file << "DATASET PARTICLES"  << "\n";
	vtk_file << "POINTS " << RW->numOfParticles << " double"
	         << "\n";

	// Write particle locations
	for (int i = 0; i < np; ++i)
		vtk_file << RW->X[i].Now.x << " " << RW->X[i].Now.y << " " << RW->X[i].Now.z << "\n";

	// KR add "vertices" block to create a correct VTK file
	vtk_file << "VERTICES " << np << " " << (2 * np) << "\n";
	for (int i = 0; i < np; ++i)
		vtk_file << 1 << " " << i << "\n";

	// Write particle identities
	vtk_file << "\n";
	vtk_file << "POINT_DATA " << RW->numOfParticles << "\n";
	vtk_file << "SCALARS identity float 1"
	         << "\n";
	vtk_file << "LOOKUP_TABLE default"
	         << "\n";
	for (int i = 0; i < np; ++i)
		vtk_file << RW->X[i].Now.identity << "\n";

	// Write particle on_boundary or not
	vtk_file << endl;
	vtk_file << "SCALARS on_boundary float 1" << endl;
	vtk_file << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < np; ++i)
		vtk_file << RW->X[i].Now.on_boundary << endl;

	// Write particle vectors
	vtk_file << endl;
	vtk_file << "VECTORS velocity float" << endl;
	for (int i = 0; i < np; ++i)
		vtk_file << RW->X[i].Now.Vx << " " << RW->X[i].Now.Vy << " " << RW->X[i].Now.Vz << endl;

	// Let's close it, now
	vtk_file.close();
}
