/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_REACT_BRNS.h"
#include "rf_pcs.h"
#include "rfmat_cp.h"
#include <iostream>

#ifdef BRNS

#ifdef BRNS_OMP
#include <omp.h>
#define NUM_THREADS 2
#endif

#ifdef USE_MPI_BRNS
//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END
#include "mpi.h" //Parallel Computing Support
#include "par_ddc.h"
// HS 07.01.2008: Comment the following 2 lines on LiClus.
// int size;
// int myrank;
#endif

#ifdef GCC
#include <dlfcn.h>
#endif

#include "display.h"

using namespace std;

REACT_BRNS::REACT_BRNS(void)
{
	// number of Components;
	num_Comp = 0;

	// set initialized flag to false by default;
	init_flag = false;

	// set timer to zero
	timeSpentInBrnsCoupling = 0.0;
	cout << "Debugging #1." << endl;
#ifdef GCC
	hDll = dlopen("./brns.so", RTLD_NOW);
	if (hDll == NULL)
		std::cout << "***error: failed to load a library ./brns.so" << std::endl;
	cout << "Debugging #2." << endl;
	invokebrns = (LPFNDLLFUNC)dlsym(hDll, "invokebrns_");
	if (invokebrns == NULL)
		std::cout << "***error: failed to find a symbol \"invokebrns_\"" << std::endl;
#endif
	cout << "Debugging #3." << endl;
}

REACT_BRNS::~REACT_BRNS(void)
{
	// reclaim the memory;
	delete[] rt_BRNS;
	delete[] cur_ts_Conc;
	delete[] pre_ts_Conc;
	delete[] m_dC_Chem_delta;
	delete[] boundary_flag;
	delete[] m_porosity_Node;
	delete[] m_porosity_Elem;

#ifdef GCC
	dlclose(hDll);
#endif

#ifdef USE_MPI_BRNS
	delete[] rt_BRNS_buf;
	delete[] pre_ts_Conc_buf;
#endif

#ifdef TESTTIME
	cout << "Total time spent in BRNS Coupling, including Solver: " << timeSpentInBrnsCoupling << "s" << endl;
#endif
}

void REACT_BRNS::TestRUN(void)
{
	cout << "Hello World!" << endl;
	cout << "Trying to call BrnsDLL ..." << endl;
	// brnsIsAlive();
	cout << "Hm, looks good, I guess." << endl;

	double* myArray = NULL;
	int sizeOfArray;
	int returnValue = -111;

	sizeOfArray = 6;
	myArray = new double[sizeOfArray];

	myArray[0] = 1.2;
	myArray[1] = 3.4;
	myArray[2] = 5.6;
	myArray[3] = 7.8;
	myArray[4] = 9.10;
	myArray[5] = 11.12;

	cout << "Initial Concentrations:" << endl;
	for (int i = 0; i < sizeOfArray; i++)
		cout << i << ": " << myArray[i] << endl;

	cout << "Calling invokebrns() ..." << endl;

	// invokebrns(myArray, sizeOfArray, &returnValue);

	cout << "Returning from invokebrns(), return code " << returnValue << endl;
	cout << "New Concentrationvector:" << endl;

	for (int i = 0; i < sizeOfArray; i++)
		cout << i << ": " << myArray[i] << endl;

	delete[] myArray;
	myArray = NULL;
}

long REACT_BRNS::GetNodesNumber(void)
{
	long number;
	number = 0;
	//------------read number of nodes--------------
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0) {
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			number = (long)m_pcs->m_msh->GetNodesNumber(false);
			return number;
		}
	}
	//------------end of reading number of nodes----
	return number;
}

long REACT_BRNS::GetElemNumber(void)
{
	long number;
	//------------read number of elems--------------
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//		if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			number = (long)m_pcs->m_msh->ele_vector.size();
			return number;
		}
	}
	//------------end of reading number of nodes----
	return 0;
}

int REACT_BRNS::GetCompsNumber(void)
{
	long number;
	number = 0;

	number = (int)cp_vec.size();

	return number;
}

void REACT_BRNS::InitBRNS(Problem* myProblem)
{
	num_Comp = GetCompsNumber();
	nNodes = GetNodesNumber();
	nElems = GetElemNumber();
	flowflag = GetFlowType_MT();

	this->myProblem = myProblem;

	// initialize the array;
	cur_ts_Conc = new double[num_Comp * nNodes];
	pre_ts_Conc = new double[num_Comp * nNodes];
	m_dC_Chem_delta = new double[num_Comp * nNodes];
	m_porosity_Node = new double[nNodes];
	m_porosity_Elem = new double[nElems];

	rt_BRNS = new int[nNodes];
	boundary_flag = new int[num_Comp * nNodes];

#ifdef USE_MPI_BRNS
	pre_ts_Conc_buf = new double[num_Comp * nNodes];
	rt_BRNS_buf = new int[nNodes];
#endif

	// Giving initial values
	for (int i = 0; i < num_Comp * nNodes; i++)
	{
		cur_ts_Conc[i] = 0.0;
		pre_ts_Conc[i] = 0.0;
		m_dC_Chem_delta[i] = 0.0;
	}
	for (long k = 0; k < nNodes; k++)
	{
		rt_BRNS[k] = 0;
		m_porosity_Node[k] = 0.0;
	}

	if (num_Comp > 0 && nNodes > 0)
		init_flag = true;
	else
	{
#ifdef MFC
		AfxMessageBox("!!! Node number and Components Number must be bigger than zero!");
#endif
		DisplayErrorMsg("!!! Node number and Components Number must be bigger than zero!");
		abort();
	}

	// Marking species with fixed concentrations (in boundary nodes)

	double BCValue = 0.0;
	int p, k;
	for (long i = 0; i < nNodes; i++)
		for (k = 0; k < num_Comp; k++)
		{
			this_pcs = NULL;
			m_cp = cp_vec[k];

			// Get the pointer to the proper PCS.
			this_pcs = PCSGet("MASS_TRANSPORT", m_cp->compname);
			if (this_pcs)
				for (p = 0; p < this_pcs->pcs_number_of_primary_nvals; ++p)
				{
					// Let's print BC and ST values
					CBoundaryConditionsGroup* m_bc_group = NULL;

					const std::string tmp_pcs_type(convertProcessTypeToString(this_pcs->getProcessType()));
					m_bc_group = BCGetGroup(tmp_pcs_type, this_pcs->pcs_primary_function_name[p]);

					// BC printing
					if (IsThisPointBCIfYesStoreValue(i, this_pcs, BCValue))
					{
						// If this node is on the fixed boudnary for this component
						boundary_flag[i * num_Comp + k] = 1;
						cout << "Node " << i << ", Comp " << k << ",Value " << BCValue << endl;
					}
					else
						// If this node is NOT on the fixed boudnary for this component
						boundary_flag[i * num_Comp + k] = 0;
				} // end of for
			// end of if (this_pcs)
			else // not getting the pointer to the proper PCS.
			{
				DisplayErrorMsg("!!! In InitBRNS, can not find corresponding PCS!");
				abort();
			}
		} // end of for components
	// fill in the porosity values.
	for (size_t i = 0; i < m_flow_pcs->m_msh->ele_vector.size(); i++)
	{
		int idx_patch = m_flow_pcs->m_msh->ele_vector[i]->GetPatchIndex();
		MediaProp = mmp_vector[idx_patch];
		m_porosity_Elem[i] = MediaProp->Porosity(i, 0.0);
	}
	// averaging to Nodes
	int idx_elem;
	double elem_volume;
	double total_volume;
	for (long i = 0; i < nNodes; i++)
	{
		if (m_flow_pcs->m_msh->nod_vector[i]->getConnectedElementIDs().size() > 0)
		{
			total_volume = 0.0;
			for (int j = 0; j < (int)m_flow_pcs->m_msh->nod_vector[i]->getConnectedElementIDs().size(); j++)
			{
				idx_elem = m_flow_pcs->m_msh->nod_vector[i]->getConnectedElementIDs()[j];
				elem_volume = m_flow_pcs->m_msh->ele_vector[idx_elem]->GetVolume();
				m_porosity_Node[i] += m_porosity_Elem[idx_elem] * elem_volume;
				total_volume += elem_volume;
			}
			m_porosity_Node[i] = m_porosity_Node[i] / total_volume;
		}
		else
			m_porosity_Node[i] = 0.0;
	}
}

void REACT_BRNS::GSRF2Buffer(long i)
{
	long i_times_num_Comp_plus_k = i * num_Comp;
	// for this node, loop over all the chemical components and update the values
	for (int k = 0; k < num_Comp; k++, i_times_num_Comp_plus_k++)
	{
		this_pcs = NULL;
		m_cp = cp_vec[k];

		// Get the pointer to the proper PCS.
		this_pcs = PCSGet("MASS_TRANSPORT", m_cp->compname);
		if (this_pcs)
		{
			// Set the Concentration of this component
			cur_ts_Conc[i_times_num_Comp_plus_k]
			    = this_pcs->GetNodeValue(i, this_pcs->GetNodeValueIndex(this_pcs->pcs_primary_function_name[0]) + 1);
			pre_ts_Conc[i_times_num_Comp_plus_k]
			    = this_pcs->GetNodeValue(i, this_pcs->GetNodeValueIndex(this_pcs->pcs_primary_function_name[0]) + 0);
		}
		else // not getting the pointer to the proper PCS.
		{
#ifdef MFC
			AfxMessageBox("!!! In Data transfer for BRNS, can not find corresponding PCS!");
#endif
			DisplayErrorMsg("!!! In Data transfer for BRNS, can not find corresponding PCS!");
			abort();
		}
	}
}

void REACT_BRNS::Buffer2GSRF(long i)
{
	long i_times_num_Comp_plus_k = i * num_Comp;
	// for this node, loop over all the chemical components and update the values
	for (int k = 0; k < num_Comp; k++, i_times_num_Comp_plus_k++)
	{
		this_pcs = NULL;
		m_cp = cp_vec[k];

		// Get the pointer to the proper PCS.
		this_pcs = PCSGet("MASS_TRANSPORT", m_cp->compname);
		if (this_pcs)
			// Set the Concentration of this component at current time step;
			this_pcs->SetNodeValue(
			    i,
			    this_pcs->GetNodeValueIndex(this_pcs->pcs_primary_function_name[0]) + 1 /*1-current time step*/,
			    pre_ts_Conc[i_times_num_Comp_plus_k]);
		else // not getting the pointer to the proper PCS.
		{
#ifdef MFC
			AfxMessageBox("!!! In Data transfer for BRNS, can not find corresponding PCS!");
#endif
			DisplayErrorMsg("!!! In Data transfer for BRNS, can not find corresponding PCS!");
			abort();
		}
	}
}

void REACT_BRNS::RUN(double time_step)
{
	long i;
	// this is the last call to the chemical solver and
	if (myProblem->getEndTime() < myProblem->getCurrentTime() + time_step)
	{
		rt_BRNS[0] = -1; // reaction rates should be printed by BRNS.dll for the whole domain
		for (i = 1; i < nNodes; i++) // -1: delete a potentially existing rate file
			rt_BRNS[i] = -2;
	}
#ifdef USE_MPI_BRNS
	if (myrank == 1)
		startTime = clock(); // in MPI: make statistics only when running on one processor
#else
	startTime = clock();
#endif

	// Loop over all nodes to transfer data
	for (i = 0; i < nNodes; i++)
		// Get Conc Data from GSRF;
		GSRF2Buffer(i);

	double pos_x;
	double pos_y;
	double pos_z;

	// double porosity;
	double waterSaturation;

//------------end of reading number of nodes----

// Run BRNS
#ifdef BRNS_OMP
	// ADD CONSIDERATION OF porosity and waterSaturation!!!
	// do it with OpenMP
	int num_Comp_temp;
	double* cur_ts_Conc_temp;
	double* pre_ts_Conc_temp;
	long nNodes_temp;
	int* boundary_flag_temp;
	num_Comp_temp = num_Comp;
	cur_ts_Conc_temp = cur_ts_Conc;
	pre_ts_Conc_temp = pre_ts_Conc;
	boundary_flag_temp = boundary_flag;
	nNodes_temp = nNodes;
	omp_set_num_threads((int)NUM_THREADS);
	cout << "Max. thread num is:" << omp_get_max_threads() << endl;

#pragma omp parallel for shared(                                                                            \
    nNodes_temp, num_Comp_temp, time_step, cur_ts_Conc_temp, pre_ts_Conc_temp, boundary_flag_temp) private( \
        i, pos_x, pos_y, pos_z)
	for (i = 0; i < nNodes_temp; i++)
	{
		pos_x = m_pcs->m_msh->nod_vector[i]->X();
		pos_y = m_pcs->m_msh->nod_vector[i]->Y();
		pos_z = m_pcs->m_msh->nod_vector[i]->Z();
		//------------------------
		// Prepare porosity and waterSaturation values (weighted sums?), or move to Buffer<->GSFR methods

		// porosity = 0.0;
		waterSaturation = 0.0;

		//------------------------

		// Run BRNS;
		invokebrns(&(cur_ts_Conc_temp[i * num_Comp_temp]),
		           &(pre_ts_Conc_temp[i * num_Comp_temp]),
		           &(pre_ts_Conc_temp[i * num_Comp_temp]),
		           &num_Comp_temp,
		           &time_step,
		           &(boundary_flag_temp[i * num_Comp_temp]),
		           &(rt_BRNS[i]),
		           &pos_x,
		           &pos_y,
		           &pos_z,
		           REACT_BRNS::m_porosity_Node + i,
		           &waterSaturation,
		           NULL);
		int num = omp_get_thread_num();
		cout << "#" << num << "thread reporting. My i is: " << i /*cur_ts_Conc_temp[i*num_Comp_temp]*/ << endl;
		// cout << endl << "Number of threads " << omp_get_num_threads ();
	}
#else // ifdef BRNS_OMP
#ifdef USE_MPI_BRNS
	MPI_Bcast(&nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&num_Comp, 1, MPI_INT, 0, MPI_COMM_WORLD);
	for (i = myrank; i < nNodes; i += mysize)
#else
	for (i = 0; i < nNodes; i++)
#endif
	{
		// Check out if this node is on the boundary - A fixed boundary condition
		pos_x = m_pcs->m_msh->nod_vector[i]->getData()[0];
		pos_y = m_pcs->m_msh->nod_vector[i]->getData()[1];
		pos_z = m_pcs->m_msh->nod_vector[i]->getData()[2];

		//------------------------
		// Prepare porosity and waterSaturation values (weighted sums?), or move to Buffer<->GSRF methods

		// porosity = 0.0;
		waterSaturation = 0.0;

//------------------------

// Run BRNS;
//	 solverTime=0.0;
//	 if (pos_x== 0 && pos_y==0) { porosity=8; pos_z=99.9; waterSaturation=3.3; solverTime=7.0; cout << "solverTime
// passed to BRNS: " << solverTime << endl; }

#ifdef USE_MPI_BRNS
		invokebrns(&(cur_ts_Conc[i * num_Comp]),
		           &(pre_ts_Conc[i * num_Comp]),
		           &(pre_ts_Conc_buf[i * num_Comp]),
		           &num_Comp,
		           &time_step,
		           &(boundary_flag[i * num_Comp]),
		           &(rt_BRNS[i]),
		           &pos_x,
		           &pos_y,
		           &pos_z,
		           REACT_BRNS::m_porosity_Node + i,
		           &waterSaturation,
		           NULL);
#else // ifdef USE_MPI_BRNS
		invokebrns(&(cur_ts_Conc[i * num_Comp]),
		           &(pre_ts_Conc[i * num_Comp]),
		           &(pre_ts_Conc[i * num_Comp]),
		           &num_Comp,
		           &time_step,
		           &(boundary_flag[i * num_Comp]),
		           &(rt_BRNS[i]),
		           &pos_x,
		           &pos_y,
		           &pos_z,
		           REACT_BRNS::m_porosity_Node + i,
		           &waterSaturation,
		           NULL);
#endif
	}
#endif
#ifdef USE_MPI_BRNS
	GetBRNSResult_MPI();
	CleanMPIBuffer();
#endif

	// calculate dC
	for (i = 0; i < nNodes * num_Comp; i++)
		m_dC_Chem_delta[i] = cur_ts_Conc[i] - pre_ts_Conc[i];

	// Loop over all nodes to retrieve data
	for (i = 0; i < nNodes; i++)
	{
		// Set data back to GSRF
		Buffer2GSRF(i);
		if (rt_BRNS[i] == 2 || rt_BRNS[i] == 3)
		{
			cout << "In Node " << i << ": BRNS exceeded max newton iterations, aborting." << endl;
			exit(10);
		}
#ifndef BRNS_NO_LOG // should be eventually changed to #ifdef BRNS_LOG
		switch (rt_BRNS[i])
		{
			case 1:
				cout << "In Node " << i << ": BRNS calculated negative concentration!" << endl;
				break;
			case 2:
				cout << "In Node " << i << ": BRNS exceeded max newton iterations!" << endl;
				break;
			case 3:
				cout << "In Node " << i << ": BRNS calculated negative concentration, "
				     << "and exceeded max newton iterations!" << endl;
				break;
				//			default:
		}
		if (rt_BRNS[i] == 1 || rt_BRNS[i] == 3)
		{
			cout << "Still negative concentrations after chemical simulation for:";
			for (int j = 0; j < num_Comp; j++)
				if (pre_ts_Conc[i * num_Comp + j] < 0.0)
					cout << " Species " << j + 1 << ": " << pre_ts_Conc[i * num_Comp + j];
			cout << endl;
		}
#endif
	}
	// do not forget to write the porosity value into the Node data structure here.
	// 1 stands for new time step.
	ConvPorosityNodeValue2Elem(0);
	ConvPorosityNodeValue2Elem(1);

#ifdef USE_MPI_BRNS
	if (myrank == 1)
		timeSpentInBrnsCoupling += (double)(clock() - startTime) / CLOCKS_PER_SEC;

#else
	timeSpentInBrnsCoupling += (double)(clock() - startTime) / CLOCKS_PER_SEC;

//			cout << "solverTime=" << solverTime << "total "<< timeSpentInBrnsSolver << endl;
#endif

#ifdef BRNS_LOG
#ifdef USE_MPI_BRNS
	if (myrank == 1)
		cout << "Total time spent in BRNS Coupling and Solver: " << timeSpentInBrnsCoupling << "s, "
		     << timeSpentInBrnsSolver << "s " << endl;

#else
	cout << "Total time spent in BRNS Coupling and Solver: " << timeSpentInBrnsCoupling << "s, "
	     << timeSpentInBrnsSolver << "s " << endl;
#endif
#endif
}

int REACT_BRNS::IsThisPointBCIfYesStoreValue(int index, CRFProcess* m_pcs, double& value)
{
	for (int p = 0; p < (int)m_pcs->bc_node_value.size(); ++p)
		if (index == m_pcs->bc_node_value[p]->msh_node_number)
		{
			value = m_pcs->bc_node_value[p]->node_value;
			return 1; // Yes, found it.
		}

	return 0;
}

// i_timestep 0: old timestep 1: new timestep
int REACT_BRNS::SetPorosityValue_MT(long ele_Index, double m_porosity_Elem, int i_timestep)
{
	int idx;
	double old_porosity;

	idx = -1;
	old_porosity = 0.0;

	CRFProcess* m_pcs = NULL;

	if (flowflag > 0)
	{
		GetFluidProperty_MT();
		switch (flowflag)
		{
			case 1:
				m_pcs = PCSGet("GROUNDWATER_FLOW");
				idx = m_pcs->GetElementValueIndex("POROSITY");

				// set new porosity;
				m_pcs->SetElementValue(ele_Index, idx + i_timestep, m_porosity_Elem);
				break;
			case 2:
				m_pcs = PCSGet("LIQUID_FLOW");
				idx = m_pcs->GetElementValueIndex("POROSITY");
				// always write into the new step
				m_pcs->SetElementValue(ele_Index, idx + 1, m_porosity_Elem);
				break;
			case 3:

				m_pcs = PCSGet("RICHARDS_FLOW");
				idx = m_pcs->GetElementValueIndex("POROSITY");
				// always write into the new step
				if (m_pcs->GetElementValue(ele_Index, idx + 1) == 0)
					m_pcs->SetElementValue(ele_Index, idx + 1, m_porosity_Elem);
				else if (m_porosity_Elem != m_pcs->GetElementValue(ele_Index, idx + 1))
				{
					cout << "Richards flow not supported - quitting program" << endl;
					exit(1);
				}
				break;
				m_pcs->SetElementValue(ele_Index, idx + 1, m_porosity_Elem);
				break;
			case 4: // kg44: do we have to update POROSITY_IL and POROSITY_SW?
				m_pcs = PCSGet("TWO_PHASE_FLOW");
				idx = m_pcs->GetElementValueIndex("POROSITY");
				// always write into the new step
				m_pcs->SetElementValue(ele_Index, idx + 1, m_porosity_Elem);
				break;
			default:
#ifdef USE_MPI_BRNS
				if (myrank == 0 /*should be set to root*/)
#endif
					DisplayErrorMsg("Error: Not implemented for the flow in BRNS case!!!");
				break;
		}
	}
	return 1;
}

void REACT_BRNS::GetFluidProperty_MT(void)
{
	m_FluidProp = MFPGet("LIQUID");
}

int REACT_BRNS::GetFlowType_MT(void)
{
	// flow type
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//		if ( m_pcs->pcs_type_name.compare ( "GROUNDWATER_FLOW" ) ==0 ) {
		if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		{
			m_flow_pcs = m_pcs;
			return 1;
		}
		//		else if ( m_pcs->pcs_type_name.compare ( "LIQUID_FLOW" ) ==0 ) {
		else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
		{
			m_flow_pcs = m_pcs;
			return 2;
		}
		//		else if ( m_pcs->pcs_type_name.compare ( "RICHARDS_FLOW" ) ==0 ) {
		else if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
		{
			m_flow_pcs = m_pcs;
			return 3;
		}
		//		else if ( m_pcs->pcs_type_name.compare ( "TWO_PHASE_FLOW" ) ==0 ) {
		else if (m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
		{
			m_flow_pcs = m_pcs;
			return 4;
		}
	}
	return 0;
}

// i_timestep: 0: old timestep 1: new timestep
void REACT_BRNS::ConvPorosityNodeValue2Elem(int i_timestep)
{
	long i, idx_Node;
	size_t j, number_of_nodes;
	double pormin = 2.0, pormax = 0.0;
	MeshLib::CNode* m_Node;
	MeshLib::CElem* m_Elem;

	for (i = 0; i < nElems; i++)
	{
		m_Elem = m_pcs->m_msh->ele_vector[i];

		// first set the parameters to zero;
		m_porosity_Elem[i] = 0.0;

		// then get the values from nodes
		for (j = 0; j < m_Elem->GetNodesNumber(false); j++)
		{
			idx_Node = m_Elem->GetNodeIndex(j); // get the connected nodes;
			m_Node = m_pcs->m_msh->nod_vector[idx_Node];
			number_of_nodes = m_Elem->GetNodesNumber(false);
			// m_porosity_Elem[i] += m_porosity[idx_Node] / number_of_nodes; // this is arithmetric mean
			// here we use harmonic mean, as porosity is used for permeability/diffusivity changes....flux in the
			// element is strongly influenced by the minimum values
			// this is for harmonic mean
			m_porosity_Elem[i] += 1.0 / m_porosity_Node[idx_Node] / number_of_nodes;
		}

		// cout << "error in ConvPorosity" << endl;
		if (m_porosity_Elem[i] > 1.e6)
			m_porosity_Elem[i] = 1.e6;
		m_porosity_Elem[i] = 1.0 / m_porosity_Elem[i];

		pormin = min(pormin, m_porosity_Elem[i]);
		pormax = max(pormax, m_porosity_Elem[i]);

		// push back porosities
		SetPorosityValue_MT(i, m_porosity_Elem[i], i_timestep);
	}
	cout << "min, max porosity: " << pormin << " " << pormax << endl;
}

#ifdef USE_MPI_BRNS
void REACT_BRNS::GetBRNSResult_MPI(void)
{
	// Retrieve the values from MPI buffer to the main memory
	MPI_Allreduce(rt_BRNS_buf, rt_BRNS, nNodes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(pre_ts_Conc_buf, pre_ts_Conc, num_Comp * nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void REACT_BRNS::CleanMPIBuffer(void)
{
	long in;
	for (in = 0; in < nNodes; in++)
		rt_BRNS_buf[in] = 0;
	for (in = 0; in < num_Comp * nNodes; in++)
		pre_ts_Conc_buf[in] = 0.0;
}
#endif
#endif // end of BRNS
