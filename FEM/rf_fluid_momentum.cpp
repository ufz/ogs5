/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// Some change

/**************************************************************************
   FEMLib-Object:
   Task: MediumProperties
   Programing:
   05/2005 PCH Implementation
**************************************************************************/

#include <iostream>
using namespace std;

#include "matrix_class.h"
#include "rf_fluid_momentum.h"
#include "rf_random_walk.h"
using namespace Math_Group;

/**************************************************************************
   FEMLib-Method: ThreeComponet
   Task: constructor
   Programing:
   02/2006 PCH Implementation
   last modification:
**************************************************************************/
PlaneSet::PlaneSet(void)
{
	eleIndex = 0;
	ratio = 0.0;
	for (int i = 0; i < 3; ++i)
		V[i] = norm[i] = 0.0;
}

/**************************************************************************
   FEMLib-Method: CrossRoad
   Task: constructor
   Programing:
   01/2006 PCH Implementation
   last modification:
**************************************************************************/
CrossRoad::CrossRoad(void)
{
	Index = 0;
	numOfThePlanes = 0;
	plane = NULL;
}

CrossRoad::~CrossRoad(void)
{
	if (plane)
		delete[] plane;
}

/**************************************************************************
   FEMLib-Method: CreatePlaneSet
   Task: Create a set of PlaneSet
   Programing:
   02/2006 PCH Implementation
   last modification:
**************************************************************************/
void CrossRoad::CreatePlaneSet(const int index)
{
	plane = new PlaneSet[index]();
}

/**************************************************************************
   FEMLib-Method: CFluidMomentum
   Task: constructor
   Programing:
   05/2005 PCH Implementation
   last modification:
**************************************************************************/
CFluidMomentum::CFluidMomentum(void)
{
	m_pcs = NULL;
	RWPTSwitch = 0; // Set to be no
}

/**************************************************************************
   FEMLib-Method: CFluidMomentum
   Task: destructor
   Programing:
   05/2005 PCH Implementation
   last modification:
**************************************************************************/
CFluidMomentum::~CFluidMomentum(void)
{
}

/**************************************************************************
   FEMLib-Method: double Execute()
   Task: compute the Darcy velocity on node
   Programing:
   05/2005 PCH Implementation
   last modification:
   01/2010 TF changed access to process type
**************************************************************************/
double CFluidMomentum::Execute(int loop_process_number)
{
	double pcs_error = 0.0;

	int no_processes = (int)pcs_vector.size();

	Create();

	// JT: Printout message. Give zero error to coupling scheme.
	std::cout << "\n      ================================================"
	          << "\n";
	std::cout << "    ->Process " << loop_process_number << ": "
	          << "FLUID_MOMENTUM"
	          << "\n";
	std::cout << "      ================================================"
	          << "\n";
	cout << "      Solving for nodal velocity..."
	     << "\n";
	cpl_max_relative_error = pcs_error;

	bool isFlow = false;
	CRFProcess* a_pcs = NULL;
	// CRFProcess *f_pcs = NULL;
	for (int k = 0; k < no_processes; k++)
	{
		a_pcs = pcs_vector[k];
		if (!a_pcs)
			continue;
		if (a_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW
		    || a_pcs->getProcessType() == FiniteElement::LIQUID_FLOW
		    || a_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW
		    || a_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW
		    || a_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
		{
			isFlow = true;
			break;
		}
	}
	for (int i = 0; i < no_processes; ++i)
	{
		m_pcs = pcs_vector[i];

		// Select the mesh whose process name has the mesh for Fluid_Momentum
		//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos) TF
		if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
			m_msh = FEMGet("RICHARDS_FLOW");
		//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos) TF
		else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
			m_msh = FEMGet("LIQUID_FLOW");
		//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos) TF
		else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			m_msh = FEMGet("GROUNDWATER_FLOW");

		//		if(m_pcs->pcs_type_name.find("FLUID_MOMENTUM")!=string::npos) TF
		if (m_pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM)
		{
			if (isFlow)
				SolveDarcyVelocityOnNode();
			else
			{
				m_msh = FEMGet("FLUID_MOMENTUM");
				string vel_file = FileName + ".vel";
				ifstream ins(vel_file.c_str());
				double vx, vy, vz;

				int nidx = m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1;
				int nidy = m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1;
				int nidz = m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1;

				for (size_t i = 0; i < m_pcs->m_msh->nod_vector.size(); i++)
				{
					ins >> vx >> vy >> vz >> ws;
					m_pcs->SetNodeValue(i, nidx, vx);
					m_pcs->SetNodeValue(i, nidy, vy);
					m_pcs->SetNodeValue(i, nidz, vz);
				}
			}
		}
	}

	// Just one time execution. Needs improvement later on.
	m_pcs = PCSGet("RANDOM_WALK");
	if (m_pcs && RWPTSwitch == 0)
	{
		if (m_msh->GetCoordinateFlag() != 32)
			ConstructFractureNetworkTopology();
		RWPTSwitch = 1;
	}

	return pcs_error;
}

/**************************************************************************
   FEMLib-Method: SolveDarcyVelocityOnNode(CRFProcess*m_pcs)
   Task: compute the Darcy velocity on node
   Programing:
   05/2005 PCH Implementation
   last modification:
**************************************************************************/
void CFluidMomentum::SolveDarcyVelocityOnNode()
{
#if !defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	int nidx1 = 0;
	long i;
	MeshLib::CElem* elem = NULL;

	CheckMarkedElement();
	fem = new CFiniteElementStd(m_pcs, m_msh->GetCoordinateFlag());

	// Checking the coordinateflag for proper solution.
	int coordinateflag = m_msh->GetCoordinateFlag();
	int dimension = 0;
	int axis = 0;
	if (coordinateflag == 10)
	{
		dimension = 1;
		axis = 0; // x only
	}
	else if (coordinateflag == 11)
	{
		dimension = 1;
		axis = 1; // y only
	}
	else if (coordinateflag == 12)
	{
		dimension = 1;
		axis = 2; // z only
	}
	else if (coordinateflag == 21)
	{
		dimension = 2;
		axis = 1; // x, y only
	}
	else if (coordinateflag == 22)
	{
		dimension = 2;
		axis = 2; // x, z only
	}
	else if (coordinateflag == 32)
	{
		dimension = 3;
		axis = 3; // x, y, z only
	}

	// Loop over three dimension to solve three velocity components
	for (int phase = 0; phase < GetRFProcessNumPhases(); phase++)
	{
		for (int d = 0; d < dimension; ++d)
		{
/* Initializations */
/* System matrix */
#if defined(NEW_EQS) // WW
			m_pcs->EQSInitialize();
#else
			SetLinearSolverType(m_pcs->getEQSPointer(), m_num); // NW
			SetZeroLinearSolver(m_pcs->getEQSPointer());
#endif

			for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
			{
				elem = m_msh->ele_vector[i];
				if (elem->GetMark()) // Marked for use
				{
					fem->ConfigElement(elem, m_num->ele_gauss_points);
					fem->Assembly(0, d);
				}
			}

			//		MXDumpGLS("rf_pcs.txt",1,m_pcs->eqs->b,m_pcs->eqs->x); //abort();
			m_pcs->IncorporateBoundaryConditions(-1, d);
// Solve for velocity
#if defined(NEW_EQS)

			double* x;
			int size = (int)m_msh->nod_vector.size(); // OK411??? long
			x = new double[size];
#if defined(LIS)
			m_pcs->EQSSolver(x); // an option added to tell FLUID_MOMENTUM for sparse matrix system.
			cout << "Solver passed in FLUID_MOMENTUM."
			     << "\n";
#endif
#else
			ExecuteLinearSolver(m_pcs->getEQSPointer());
#endif

			/* Store solution vector in model node values table */
			if (dimension == 1)
				nidx1 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[axis]) + 1;
			else if (dimension == 2)
			{
				if (axis == 1) // x,y only
					nidx1 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[(axis - d + 1) % dimension]) + 1;
				else if (axis == 2) // x,z only
					nidx1 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[(axis - d + 1) % 3]) + 1;
				else
					abort(); // Just stop something's wrong.
			}
			else if (dimension == 3)
				nidx1 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[d]) + 1;
			else
				abort(); // Just stop something's wrong.

#if defined(NEW_EQS)
			for (int j = 0; j < size; j++)
				m_pcs->SetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1, x[j]);

			delete[] x;
#else
			LINEAR_SOLVER* eqs = m_pcs->getEQSPointer();
			for (int j = 0; j < eqs->dim; j++)
				m_pcs->SetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1, eqs->x[j]);
#endif
		}

		/*
		   if(m_msh->GetCoordinateFlag() == 32)
		   ; // do nothing
		   else
		   {
		   ConstructFractureNetworkTopology();
		   for(i = 0; i < (long)m_msh->nod_vector.size(); i++)
		   {
		      // Let's get the norm of the first connected element plane.
		      double norm[3];
		      if(m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
		   {
		   norm[0] = 0.0; norm[1] = 0.0; norm[2] = 1.0;
		   }
		   else
		   {
		   // I assume that all the element stay on the same plane.
		   // So, I use element No.1 as the reference element or plane
		   for(int j=0; j<3; ++j)
		   norm[j] = m_msh->ele_vector[0]->getTransformTensor(j+6);
		   }

		   // Do some proper projection of velocity computed from Fluid Momentum.
		   // Get the fluid velocity for this node
		   double V[3];
		   V[0] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_X")+1);
		   V[1] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Y")+1);
		   V[2] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Z")+1);

		   // Let's solve the projected velocity on the element plane
		   // by  Vp = norm X (V X norm) assuming norm is a unit vector
		   double VxNorm[3], Vp[3];
		   CrossProduction(V,norm,VxNorm);
		   CrossProduction(norm,VxNorm, Vp);

		   // Store the projected velocity back to the node velocity
		   m_pcs->SetNodeValue(i,m_pcs->GetNodeValueIndex("VELOCITY1_X")+1,Vp[0]);
		   m_pcs->SetNodeValue(i,m_pcs->GetNodeValueIndex("VELOCITY1_Y")+1,Vp[1]);
		   m_pcs->SetNodeValue(i,m_pcs->GetNodeValueIndex("VELOCITY1_Z")+1,Vp[2]);
		   }
		   }
		 */
		// Obtain the edge velocity
		SolveForEdgeVelocity();

		// Obtain element-based velocity
		for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
		{
			elem = m_msh->ele_vector[i];

			double vx = 0.0, vy = 0.0, vz = 0.0;
			int numOfNodeInElement = elem->GetVertexNumber();

			for (int j = 0; j < numOfNodeInElement; ++j)
			{
				vx += m_pcs->GetNodeValue(elem->GetNodeIndex(j), m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
				vy += m_pcs->GetNodeValue(elem->GetNodeIndex(j), m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
				vz += m_pcs->GetNodeValue(elem->GetNodeIndex(j), m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);
			}
			vx /= (double)numOfNodeInElement;
			vy /= (double)numOfNodeInElement;
			vz /= (double)numOfNodeInElement;

			/*
			         switch(phase)
			         {
			            case 0:
			 */

			m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY1_X") + 1, vx);
			m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY1_Y") + 1, vy);
			m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY1_Z") + 1, vz);

			/*
			               break;
			            case 1:
			               m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY2_X")+1, vx);
			                    m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY2_Y")+1, vy);
			                    m_pcs->SetElementValue(i, m_pcs->GetElementValueIndex("VELOCITY2_Z")+1, vz);
			               break;
			            default:
			               cout << "Error in VELCalcElementVelocity: invalid phase number" << "\n";
			         }
			 */
		}
	}

	// Release memroy
	delete fem;
#endif
}

/**************************************************************************
   FEMLib-Method: void Create()
   Task: This only creates NUM nothing else
   Programing:
   06/2005 PCH Implementation
   last modification:
**************************************************************************/
void CFluidMomentum::Create()
{
	// NUM_NEW
	int no_numerics = (int)num_vector.size();

	CNumerics* m_num_tmp = NULL;

	for (int i = 0; i < no_numerics; i++)
	{
		m_num_tmp = num_vector[i];

		if (m_num_tmp->pcs_type_name.compare("FLUID_MOMENTUM") == 0)
			m_num = m_num_tmp;
	}
	if (!m_num)
		cout << "Warning in CRFProcess::Create() - no numerical properties"
		     << "\n";
}

/**************************************************************************
   FEMLib-Method: void ConstructFractureNetworkTopology()
   Task: The function constructs the joint edges and crossroad nodes.
   Programing:
   01/2006 PCH Implementation
   last modification:
   10/2010 TF changed access to process type
**************************************************************************/
void CFluidMomentum::ConstructFractureNetworkTopology()
{
	// Mount the process and the mesh
	CFEMesh* m_msh = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); ++i)
	{
		m_pcs = pcs_vector[i];

		// Select the mesh whose process name has the mesh for Fluid_Momentum
		//		if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos) TF
		if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
			m_msh = FEMGet("RICHARDS_FLOW");
		//		else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos) TF
		else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
			m_msh = FEMGet("LIQUID_FLOW");
		//		else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos) TF
		else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			m_msh = FEMGet("GROUNDWATER_FLOW");
	}
	m_pcs = PCSGet("FLUID_MOMENTUM");
	if (!m_msh)
		m_msh = m_pcs->m_msh;
	// Something must be done later on here.
	double tolerance = 1e-12;

	// Checking the node is a crossroad starts here
	// Loop over all the nodes
	for (int i = 0; i < (int)m_msh->nod_vector.size(); ++i)
	{
		MeshLib::CNode* thisNode = m_msh->nod_vector[i];
		int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();

		// Let's get the norm of the first connected element plane.
		double norm[3];
		int index = thisNode->getConnectedElementIDs()[0];
		// Let's store the index of the reference element
		// to the connected_planes of thisNode
		thisNode->connected_planes.push_back(index);
		if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
		{
			norm[0] = 0.0;
			norm[1] = 0.0;
			norm[2] = 1.0;
		}
		else
			for (int j = 0; j < 3; ++j)
				norm[j] = m_msh->ele_vector[index]->getTransformTensor(j + 6);

		// Let's compare this norm with other norms of the connected elements
		for (int j = 1; j < NumOfNeighborElements; ++j)
		{
			double normOther[3];
			// Let's get the element one by one.
			int indexOther = thisNode->getConnectedElementIDs()[j];
			if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
			{
				normOther[0] = 0.0;
				normOther[1] = 0.0;
				normOther[2] = 1.0;
			}
			else
				for (int k = 0; k < 3; ++k)
					normOther[k] = m_msh->ele_vector[indexOther]->getTransformTensor(k + 6);

			// Check two norms are same.
			if (fabs(norm[0] - normOther[0]) < tolerance && fabs(norm[1] - normOther[1]) < tolerance
			    && fabs(norm[2] - normOther[2]) < tolerance)
				; // Two elements stay on the same plane
			else
			{
				thisNode->crossroad = true;
				// I am going to store all the element indeces which are different from
				// the reference element. Then, I will get rid of the duplicate elements
				// of the same plane.
				int indexOther = -1;
				indexOther = thisNode->getConnectedElementIDs()[j];
				thisNode->connected_planes.push_back(indexOther);
			}
		}

		// Get rid of duplicates of the elements that have the same norm for the potential crossroads
		// This works great with h_frac example in RWPT
		if (thisNode->crossroad)
		{
			// Let's get rid of duplicates of the elements on the same plane
			int numOfPlanesAtCrossroad = (int)thisNode->connected_planes.size();
			for (int j = 0; j < numOfPlanesAtCrossroad; ++j)
			{
				double normOther[3];
				// Let's get the element one by one.
				int indexOther = thisNode->connected_planes[j];
				if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
				{
					normOther[0] = 0.0;
					normOther[1] = 0.0;
					normOther[2] = 1.0;
				}
				else
					for (int k = 0; k < 3; ++k)
						normOther[k] = m_msh->ele_vector[indexOther]->getTransformTensor(k + 6);

				for (int l = j + 1; l < numOfPlanesAtCrossroad; ++l)
				{
					double normAnother[3];
					// Let's get the element one by one.
					int indexAnother = thisNode->connected_planes[l];
					if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
					{
						normAnother[0] = 0.0;
						normAnother[1] = 0.0;
						normAnother[2] = 1.0;
					}
					else
						for (int k = 0; k < 3; ++k)
							normAnother[k] = m_msh->ele_vector[indexAnother]->getTransformTensor(k + 6);

					// Check two norms are same.
					// If two norms of the elemenets are same,
					if (fabs(normOther[0] - normAnother[0]) < tolerance
					    && fabs(normOther[1] - normAnother[1]) < tolerance
					    && fabs(normOther[2] - normAnother[2]) < tolerance)
					{
						// Two elements stay on the same plane
						for (int m = l; m < (numOfPlanesAtCrossroad - 1); ++m)
							thisNode->connected_planes[m] = thisNode->connected_planes[m + 1];

						// Erase the element of the vector and adjust the number of the planes at crossroad
						thisNode->connected_planes.erase(thisNode->connected_planes.begin() + numOfPlanesAtCrossroad
						                                 - 1);
						numOfPlanesAtCrossroad = (int)thisNode->connected_planes.size();
						--l; // Very important. Huh.
					}
				}
			}
		}
		// Now we got the number of connected planes on nodes.
		// Getting rid of duplicates ends here.

		// Do some proper projection of velocity computed from Fluid Momentum.
		// Get the fluid velocity for this node
		double V[3];
		V[0] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
		V[1] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
		V[2] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);
		if (thisNode->crossroad == false)
		{
			// Let's solve the projected velocity on the element plane
			// by  Vp = norm X (V X norm) assuming norm is a unit vector
			double VxNorm[3], Vp[3];
			CrossProduction(V, norm, VxNorm);
			CrossProduction(norm, VxNorm, Vp);

			// Store the projected velocity back to the node velocity
			m_pcs->SetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1, Vp[0]);
			m_pcs->SetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1, Vp[1]);
			m_pcs->SetNodeValue(i, m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1, Vp[2]);
		}
		else
		{
			// For the velocity of the node that has more than one connected planes.
			double Vmag = 0.0;
			// Mount cross
			CrossRoad* thisCross;
			thisCross = new CrossRoad();
			thisCross->numOfThePlanes = (int)thisNode->connected_planes.size();
			thisCross->CreatePlaneSet(thisCross->numOfThePlanes);
			// Loop over the number of connected planes
			for (int j = 0; j < thisCross->numOfThePlanes; ++j)
			{
				// Some local variables within this else
				double norm[3], VxNorm[3], Vp[3];
				// Solve for the norm of this plane.
				int index = thisNode->connected_planes[j];
				if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
				{
					thisCross->plane[j].norm[0] = norm[0] = 0.0;
					thisCross->plane[j].norm[1] = norm[1] = 0.0;
					thisCross->plane[j].norm[2] = norm[2] = 1.0;
				}
				else
					for (int k = 0; k < 3; ++k)
						thisCross->plane[j].norm[k] = norm[k] = m_msh->ele_vector[index]->getTransformTensor(k + 6);

				// Store the position vector defined below
				double const* CenterOfEle(m_msh->ele_vector[index]->GetGravityCenter());
				double const* const coords(thisNode->getData());
				thisCross->plane[j].Eele[0] = CenterOfEle[0] - coords[0]; // thisNode->X();
				thisCross->plane[j].Eele[1] = CenterOfEle[1] - coords[1]; // thisNode->Y();
				thisCross->plane[j].Eele[2] = CenterOfEle[2] - coords[2]; // thisNode->Z();

				// Solve the velocity contribution for this plane.
				CrossProduction(V, norm, VxNorm);
				CrossProduction(norm, VxNorm, Vp);
				for (int k = 0; k < 3; ++k)
					thisCross->plane[j].V[k] = Vp[k];

				// For ratio and Vmag
				thisCross->plane[j].ratio = sqrt(Vp[0] * Vp[0] + Vp[1] * Vp[1] + Vp[2] * Vp[2]);
				Vmag += thisCross->plane[j].ratio;
				// Update the eleIndex for this plane
				thisCross->plane[j].eleIndex = thisNode->connected_planes[j];
			}
			// Let's sort the contribution of each plane and the index of this node
			for (int j = 0; j < thisCross->numOfThePlanes; ++j)
				thisCross->plane[j].ratio = thisCross->plane[j].ratio / Vmag;
			thisCross->Index = i;
			// velocity for each connected planes ends here.

			// Extract the real crossroads from all the potential crossroads.
			// First, we will solve vectors from each potential crossroad to one of the center
			// of the connected elements.
			// Now some dot product.
			for (int j = 0; j < thisCross->numOfThePlanes; ++j)
			{
				double angle = 0.0;
				double a[3], b[3];
				for (int p = 0; p < 3; ++p)
				{
					a[p] = thisCross->plane[j].Eele[p];
					b[p] = thisCross->plane[j].V[p];
				}
				// Let's normalize these two vectors first.
				NormalizeVector(a, 3);
				NormalizeVector(b, 3);
				double dotProduct = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

				angle = acos(dotProduct);

				// If this angle is bigger than Pi/2 (90 degree),
				// then this crossroad is not a realone.
				if (angle > PI / 2.0)
					thisNode->crossroad = false;
			}
			// Extraction ends here.

			// Now add this crossroad to the vector of all crossroads in the domain
			if (thisNode->crossroad)
				crossroads.push_back(thisCross);
		}
	}
	// Checking the node is a crossroad ends here

	// Checking the edge is a joint starts here
	// Loop over all the edges
	for (int i = 0; i < (int)m_msh->edge_vector.size(); ++i)
	{
		// Mount the nodes of the edge
		vec<MeshLib::CNode*> theNodesOfThisEdge(3);
		m_msh->edge_vector[i]->GetNodes(theNodesOfThisEdge);

		// Do some proper projection of velocity computed from Fluid Momentum.
		// Get the fluid velocity for this edge
		double V[3], V0[3], V1[3];
		V0[0] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
		V0[1] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
		V0[2] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);
		V1[0] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
		V1[1] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
		V1[2] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);
		V[0] = (V0[0] + V1[0]) / 2.0;
		V[1] = (V0[1] + V1[1]) / 2.0;
		V[2] = (V0[2] + V1[2]) / 2.0;

		if (theNodesOfThisEdge[0]->crossroad && theNodesOfThisEdge[1]->crossroad)
		{
			m_msh->edge_vector[i]->SetJoint(1);

			// Constructing the connected elements of an edge starts here
			int numOfCEfromNode0 = (int)theNodesOfThisEdge[0]->getConnectedElementIDs().size();
			int numOfCEfromNode1 = (int)theNodesOfThisEdge[1]->getConnectedElementIDs().size();

			for (int j = 0; j < numOfCEfromNode0; ++j)
				for (int k = 0; k < numOfCEfromNode1; ++k)
				{
					int indexOfCEfromNode0 = theNodesOfThisEdge[0]->getConnectedElementIDs()[j];
					int indexOfCEfromNode1 = theNodesOfThisEdge[1]->getConnectedElementIDs()[k];
					if (indexOfCEfromNode0 == indexOfCEfromNode1)
						m_msh->edge_vector[i]->connected_elements.push_back(indexOfCEfromNode0);
				}
			// Constructing the connected elements of an edge ends here

			// Construct Joint vectors
			// For the velocity of the edge that has more than one connected planes.
			double Vmag = 0.0;
			CrossRoad* thisJoint;
			thisJoint = new CrossRoad();
			thisJoint->numOfThePlanes = (int)m_msh->edge_vector[i]->connected_elements.size();
			thisJoint->CreatePlaneSet(thisJoint->numOfThePlanes);

			// Loop over the number of connected planes
			for (int j = 0; j < thisJoint->numOfThePlanes; ++j)
			{
				// Some local variables within this else
				double norm[3], VxNorm[3], Vp[3];
				// Solve for the norm of this plane.
				int index = m_msh->edge_vector[i]->connected_elements[j];
				if (m_msh->GetCoordinateFlag() != 32 && m_msh->GetCoordinateFlag() != 22)
				{
					norm[0] = 0.0;
					norm[1] = 0.0;
					norm[2] = 1.0;
				}
				else
					for (int k = 0; k < 3; ++k)
						thisJoint->plane[j].norm[k] = norm[k] = m_msh->ele_vector[index]->getTransformTensor(k + 6);

				// Store the position vector defined below
				double const* CenterOfEle(m_msh->ele_vector[index]->GetGravityCenter());
				// I am using the center position of the joint
				//            thisJoint->plane[j].Eele[0] = CenterOfEle[0] -
				//            (theNodesOfThisEdge[0]->X()+theNodesOfThisEdge[1]->X())/2.0;
				//            thisJoint->plane[j].Eele[1] = CenterOfEle[1] -
				//            (theNodesOfThisEdge[0]->Y()+theNodesOfThisEdge[1]->Y())/2.0;
				//            thisJoint->plane[j].Eele[2] = CenterOfEle[2] -
				//            (theNodesOfThisEdge[0]->Z()+theNodesOfThisEdge[1]->Z())/2.0;
				double const* const pnt0(theNodesOfThisEdge[0]->getData());
				double const* const pnt1(theNodesOfThisEdge[1]->getData());
				thisJoint->plane[j].Eele[0] = CenterOfEle[0] - (pnt0[0] + pnt1[0]) / 2.0;
				thisJoint->plane[j].Eele[1] = CenterOfEle[1] - (pnt0[1] + pnt1[1]) / 2.0;
				thisJoint->plane[j].Eele[2] = CenterOfEle[2] - (pnt0[2] + pnt1[2]) / 2.0;

				// Solve the velocity contribution for this plane.
				CrossProduction(V, norm, VxNorm);
				CrossProduction(norm, VxNorm, Vp);
				for (int k = 0; k < 3; ++k)
					thisJoint->plane[j].V[k] = Vp[k];

				// For ratio and Vmag
				thisJoint->plane[j].ratio = sqrt(Vp[0] * Vp[0] + Vp[1] * Vp[1] + Vp[2] * Vp[2]);
				Vmag += thisJoint->plane[j].ratio;
				// Update the eleIndex for this plane
				thisJoint->plane[j].eleIndex = m_msh->edge_vector[i]->connected_elements[j];
			}
			// Let's sort the contribution of each plane and the index of this node
			for (int j = 0; j < thisJoint->numOfThePlanes; ++j)
				thisJoint->plane[j].ratio = thisJoint->plane[j].ratio / Vmag;
			thisJoint->Index = i;
			// velocity for each connected planes ends here.

			// Extract the real joints from all the potential joints.
			// First, we will solve vectors from each potential joint to one of the center
			// of the connected elements.
			// Now some dot product.
			for (int j = 0; j < thisJoint->numOfThePlanes; ++j)
			{
				double angle = 0.0;
				double a[3], b[3];
				for (int p = 0; p < 3; ++p)
				{
					a[p] = thisJoint->plane[j].Eele[p];
					b[p] = thisJoint->plane[j].V[p];
				}
				// Let's normalize these two vectors first.
				NormalizeVector(a, 3);
				NormalizeVector(b, 3);
				double dotProduct = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

				angle = acos(dotProduct);

				// If this angle is bigger than Pi/2 (90 degree),
				// then this crossroad is not a realone.
				if (angle > PI / 2.0)
					m_msh->edge_vector[i]->SetJoint(0);
			}
			// Extraction ends here.

			// Now add this crossroad to the vector of all crossroads in the domain
			if (m_msh->edge_vector[i]->GetJoint() == 1)
				joints.push_back(thisJoint);
		}
	}
	// Checking the edge is a joint ends here

	// Compute the angles of the element for rotation.
	// This process can be embedded into CElem class when in need.
	for (int i = 0; i < (int)m_msh->ele_vector.size(); ++i)
		m_msh->PT->SolveAnglesOfTheElment(m_msh->ele_vector[i]);
}

/**************************************************************************
   Task: The function solves constant edge velocity in 2D
   01/2007 PCH Implementation
   last modification:
**************************************************************************/
void CFluidMomentum::SolveForEdgeVelocity(void)
{
	// Mount the process and the mesh
	m_pcs = PCSGet("FLUID_MOMENTUM");
	CFEMesh* m_msh = fem_msh_vector[0]; // Something must be done later on here.
	double tolerance = 1e-12;

	// Checking the edge is a joint starts here
	// Loop over all the edges
	// Mount the nodes of the edge
	vec<MeshLib::CNode*> theNodesOfThisEdge(3);
	for (int i = 0; i < (int)m_msh->edge_vector.size(); ++i)
	{
		m_msh->edge_vector[i]->GetNodes(theNodesOfThisEdge);

		double const* const coords1(theNodesOfThisEdge[1]->getData());
		double const* const coords0(theNodesOfThisEdge[0]->getData());
		// Norma vector of the edge
		double VectorOfEdge[3] = {coords1[0] - coords0[0], coords1[1] - coords0[1], 0.0};

		// Note MagOfVector is never zero
		double MagOfVector = sqrt(VectorOfEdge[0] * VectorOfEdge[0] + VectorOfEdge[1] * VectorOfEdge[1]);
		// Now VectorOfEdge is unit vector
		VectorOfEdge[0] = VectorOfEdge[0] / MagOfVector;
		VectorOfEdge[1] = VectorOfEdge[1] / MagOfVector;
		VectorOfEdge[2] = 0.0;
		// TF Why not simply set to one at this points, i.e. MagOfVector = 1.0?
		MagOfVector = sqrt(VectorOfEdge[0] * VectorOfEdge[0] + VectorOfEdge[1] * VectorOfEdge[1]);

		// Geting velocity on two ending nodes of the edge
		double V[3], V0[3], V1[3], Ve0[3], Ve1[3];

		V0[0] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
		V0[1] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
		V0[2] = m_pcs->GetNodeValue(theNodesOfThisEdge[0]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);
		V1[0] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_X") + 1);
		V1[1] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Y") + 1);
		V1[2] = m_pcs->GetNodeValue(theNodesOfThisEdge[1]->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY1_Z") + 1);

		double cos, sin, MagOfV0, MagOfV1;
		MagOfV0 = sqrt(V0[0] * V0[0] + V0[1] * V0[1]);
		MagOfV1 = sqrt(V1[0] * V1[0] + V1[1] * V1[1]);

		// There are two unit normal vectors for edge
		// such as (-y,x) or (y,-x) for (x,y) edge vector
		double UNVOfE[2];
		// I will only use (-y,x) for rightys
		UNVOfE[0] = -VectorOfEdge[1];
		UNVOfE[1] = VectorOfEdge[0];

		// Check if denominator is zero
		if (MagOfV0 < tolerance)
		{
			Ve0[0] = 0.0;
			Ve0[1] = 0.0;
			Ve0[2] = 0.0;
		}
		else
		{
			// For Ve0 first
			cos = (VectorOfEdge[0] * V0[0] + VectorOfEdge[1] * V0[1]) / (MagOfVector * MagOfV0);
			sin = sqrt(1.0 - cos * cos);

			double magOfNormalV = sin * MagOfV0;

			// Therefore V0 which is normal to the edge
			Ve0[0] = magOfNormalV * UNVOfE[0];
			Ve0[1] = magOfNormalV * UNVOfE[1];
			Ve0[2] = 0.0;
		}
		// Now for the other node of the edge
		if (MagOfV1 < tolerance)
		{
			Ve1[0] = 0.0;
			Ve1[1] = 0.0;
			Ve1[2] = 0.0;
		}
		else
		{
			// For Ve1 second
			cos = (VectorOfEdge[0] * V1[0] + VectorOfEdge[1] * V1[1]) / (MagOfVector * MagOfV1);
			sin = sqrt(1.0 - cos * cos);

			double magOfNormalV = sin * MagOfV1;

			// Therefore V0 which is normal to the edge
			Ve1[0] = magOfNormalV * UNVOfE[0];
			Ve1[1] = magOfNormalV * UNVOfE[1];
			Ve1[2] = 0.0;
		}

		// Solve for normal velocity vector for the edge at the midpoint.
		V[0] = (Ve0[0] + Ve1[0]) / 2.0;
		V[1] = (Ve0[1] + Ve1[1]) / 2.0;
		V[2] = 0.0;
		double Vc[3];
		Vc[0] = (V0[0] + V1[0]) / 2.0;
		Vc[1] = (V0[1] + V1[1]) / 2.0;
		Vc[2] = 0.0;

		// Checking the angle between node averaged velocity and flux averaged velocity along the edge
		double angle = acos((V[0] * Vc[0] + V[1] * Vc[1]) / (sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2])
		                                                     * sqrt(Vc[0] * Vc[0] + Vc[1] * Vc[1] + Vc[2] * Vc[2])));
		if (angle > 3.141592 / 2.0)
		{
			V[0] = -V[0];
			V[1] = -V[1];
		}

		// Assign this velocity to the edge.
		m_msh->edge_vector[i]->AllocateMeomoryforV();
		m_msh->edge_vector[i]->SetVelocity(V);
	}
}

void FMRead(string file_base_name)
{
	(void)file_base_name;
	// Fluid_Momentum memory allocation is moved here. by PCH
	CRFProcess* m_pcs = PCSGet("FLUID_MOMENTUM");
	// WW  if(!m_pcs)
	if (m_pcs) // WW
	{
		CFEMesh* m_msh = fem_msh_vector[0]; // Something must be done later on here.
		m_msh->fm_pcs = new CFluidMomentum();
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: DATWriteFile
   Task: Write PCT file
   Programing:
   10/2005   PCH   Implementation
**************************************************************************/
/* //WW
void DATWriteHETFile(const char* file_name)
{
    FILE* tet_file = NULL;
    char tet_file_name[MAX_ZEILE];
    CFEMesh* m_msh = NULL;
    m_msh = fem_msh_vector[0];            // Something must be done later on here.
    MeshLib::CElem* elem = NULL;

    sprintf(tet_file_name,"%s.%s",file_name,"tet");
    tet_file = fopen(tet_file_name,"w+t");
    // Obtain element-based velocity
    for (int i = 0; i < (long)m_msh->ele_vector.size(); i++)
    {
        elem = m_msh->ele_vector[i];
        double const* center(elem->GetGravityCenter());

        fprintf(tet_file, "%17.12e %17.12e %17.12e %17.12e\n",
                center[0], center[1], center[2], elem->mat_vector(0) * 1e7);
    }

    // Let's close it, now
    fclose(tet_file);
}
*/
