/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * rf_REACT_GEM.cpp
 * Haibing Shao 25.03.08
 * haibing.shao@ufz.de
 * GEM Reaction Package
 * based on the PSI node-GEM source code
 * using the node-GEM code from Paul Sherrer Institute (PSI)
 * current maintainer: Georg Kosakowski
 * georg.kosakowski@psi.ch
 * last major changes: 01 April 2010 start of doxygen documentation
 *	       03 April 2010 start of cleaning code and reducing memory
 *consumption/communication for ! 05.10.2010  extend doxygen documentation,
 *kintetics with solid solutions and richards flow coupling
 *             22. Feb. 2012 extended for multi-threading, needs boost library
 *and fixed GEMS3K kernel June 2012 adjusted Richards flow coupling and added
 *density dependence Code description
 *
 * The files rf_REACT_GEM.cpp and rf_REACT_GEM.h contain the main core modules
 *of the OpenGeosys - GEMIPM2K coupling. GEMIPM2K is the calculation kernel of
 *GEMS-PSI (http://gems.web.psi.ch). GEMS-PSI executables for various platforms
 *are freely
 * availabe for download.
 * The the kernel GEMIPM2K source code is available on request.
 * GEMIPM2K is currently coupled in a non-iterave sequential way to groundwater
 *flow and transport. The coupling to the Richards flow module is under
 *development.
 *
 */

// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI_GEMS)
#include "mpi.h"  //Parallel Computing Support
#include "par_ddc.h"
// HS 07.01.2008: Comment the following 2 lines on LiClus.
// int size;
// int myrank;
#endif

#include "msh_elem.h"
#include "msh_node.h"
#include "rf_REACT_GEM.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"
#include "rfmat_cp.h"
// GeoSys-FEMLib for Gauss points
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
// -----------------------
#include "files0.h"
// Headers for shuffling
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <cstring>
#include <iostream>
//
#ifdef _WIN32
#include "direct.h"  // on win32 and win64 platform
#else
#include "stdlib.h"
#include "unistd.h"  // on unix/linux platform
// use for time measurements
#include <sys/time.h>
#endif
using namespace std;

#ifdef GEM_REACT
/**
   REACT_GEM() is the main constructor. It initializes several variables with
   default values.
 */
REACT_GEM::REACT_GEM(void)
{
    REACT_GEM::dch_input_file_path = "calcite-dch.dat";
    REACT_GEM::ipm_input_file_path = "calcite-ipm.dat";
    REACT_GEM::dbr_input_file_path = "calcite-dbr-0-0000.dat";
    REACT_GEM::dbr_bc_input_file_path = "calcite-dbr-0-0001.dat";
    REACT_GEM::init_input_file_path = "calcite-init.dat";
    REACT_GEM::init_input_file_sig = "-init.dat";

    nIC = 0;
    nDC = 0;
    nPH = 0;
    nPS = 0;
    nNodes = 0;
    nElems = 0;
    idx_water = -1;
    idx_oxygen = -1;
    idx_hydrogen = -1;

    initialized_flag = 0;
    heatflag = 0;
    flowflag = 0;
    flag_porosity_change = 1;        // 0-not coupled;1=coupled;
    min_possible_porosity = 1.e-4;   // minimum porostiy in case of changing
                                     // porosity: avoid zero porosity
    max_possible_porosity = 0.9999;  // max porosity
    flag_coupling_hydrology = 1;     // 0-not coupled;1=coupled;
    m_gem_temperature = 298.15;  // default gem temperature is 25°C, in Kelvin!
    m_gem_pressure = 1.0e+5;     // default pressure 1 bar
    flag_iterative_scheme = 0;   // 0-not iteration;1=iteration;
    flag_disable_gems = 0;       // always calculate gems
    // flag for different iterative scheme
    // 0 - sequential non-iterative scheme
    // 1 - standard iterative scheme
    // 2 - symetric iterative scheme
    // 3 - strang splitting scheme
    flag_transport_b = 1;  // default is transport full speciation; 1: transport
                           // only dissolved components of b vector
    gem_mass_scale = 1.0e-0;  // GEMS default mass scaling parameter
    flag_calculate_boundary_nodes =
        0;  // set to zero if boundary nodes should not change the porosity
    m_max_failed_nodes = 5;  // default number of max allowed nodes to fail
    m_diff_gems = 0.0;

    Fem_Ele_Std = NULL;

    gem_nThread = 1;  // default number of threads
    string tinit_path = " ";

    boost::barrier* gem_barrier_start;
    boost::barrier* gem_barrier_finish;
    boost::thread* gemThread;
    boost::mutex rwmutex, getnode_mutex;
    // the next two definitions are set to default values such that they can be
    // used in serial version
    myrank = 0;
    mysize = 1;
}

REACT_GEM::~REACT_GEM(void)
{
    if (initialized_flag > 0)
    {
        delete[] m_xDC;
        delete[] m_gam;
        delete[] m_xPH;
        delete[] m_aPH;
        delete[] m_bSP;
        delete[] m_vPS;
        delete[] m_mPS;
        delete[] m_bPS;
        delete[] m_xPA;
        delete[] m_dul;
        delete[] m_dll;
        delete[] m_uIC;
        delete[] m_bIC;
        delete[] m_bIC_dummy;
        delete[] m_rMB;
        delete[] m_xDC_pts;
        delete[] m_soluteB_pts;
        delete[] m_bIC_pts;
        delete[] m_xDC_MT_delta;
        delete[] m_xDC_Chem_delta;
        delete[] m_NodeHandle;

        delete[] m_NodeStatusCH;
        delete[] m_IterDone;
        delete[] m_IterDoneIndex;
        delete[] m_IterDoneCumulative;
        delete[] m_T;
        delete[] m_P;
        delete[] m_Vs;
        delete[] m_Ms;
        delete[] m_Gs;
        delete[] m_Hs;
        delete[] m_IC;
        delete[] m_pH;
        delete[] m_pe;
        delete[] m_Eh;
        delete[] m_porosity;
        delete[] m_porosity_initial;
        delete[] m_excess_water;
        delete[] m_excess_gas;
        delete[] m_Node_Volume;
        delete[] m_gas_volume;
        delete[] m_fluid_volume;
        delete[] m_fluid_density;
        delete[] m_soluteB;
        delete[] m_soluteB_buff;
        // delete MPI buffer--------
        delete[] m_NodeHandle_buff;
        delete[] m_NodeStatusCH_buff;
        delete[] m_IterDone_buff;

        delete[] m_Vs_buff;
        delete[] m_Ms_buff;
        delete[] m_Gs_buff;
        delete[] m_Hs_buff;
        delete[] m_IC_buff;
        delete[] m_pH_buff;
        delete[] m_pe_buff;
        delete[] m_Eh_buff;

        delete[] m_xDC_buff;
        delete[] m_xPH_buff;
        delete[] m_xPA_buff;
        delete[] m_excess_water_buff;
        delete[] m_excess_gas_buff;
        delete[] m_porosity_buff;
        delete[] m_boundary;
        delete[] m_gas_volume_buff;
        delete[] m_fluid_volume_buff;
        delete[] m_fluid_density_buff;
        delete[] m_dul_buff;
        delete[] m_dll_buff;
        delete[] m_xDC_pts_buff;
        delete[] m_xDC_MT_delta_buff;
        delete[] m_xDC_Chem_delta_buff;
        delete[] m_aPH_buff;
        delete[] m_bIC_buff;
        delete[] m_bIC_dummy_buff;
        delete[] m_porosity_Elem_buff;
        delete[] m_porosity_Elem;
        // -------------------------

        delete[] mol_phase;
        delete[] omega_phase;
        delete[] omega_components;

        delete[] dmdt;
        delete[] omega_phase_buff;       // this we need for kinetics
        delete[] mol_phase_buff;         // this we need for kinetics
        delete[] omega_components_buff;  // this we need for kinetics

        delete[] dmdt_buff;

        m_flow_pcs = NULL;
        m_kin.clear();
    }
}

/**
   short REACT_GEM::Init_Nodes ( string Project_path )

   Initialization of the GEM TNode Class

   Here we read the files needed as input for initializing GEMIPM2K
   The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)

 */
short REACT_GEM::Init_Nodes(string Project_path)
{
    long ii = 0, i = 0, in = 0;
// make sure we have the correct values
#if defined(USE_PETSC)
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
    MPI_Comm_size(PETSC_COMM_WORLD, &mysize);

    glob_NodesNumber_Linear = GetGlobalNodeNumber_MT();
    loc_NodesNumber_Linear = GetLocalNodeNumber_MT();
    NodesNumber_Linear = GetNodeNumber_MT();
    gem_glob_buff =
        new PetscScalar[glob_NodesNumber_Linear];  // m_size is defined in the
                                                   // petsc part...this
                                                   // re-usable
    // buffer has size of vector in global equation system
    gem_glob_x0 = new PetscScalar[glob_NodesNumber_Linear];
    gem_glob_x1 = new PetscScalar[glob_NodesNumber_Linear];
#endif

    TNode* m_Node;
    // DATABR structure for exchange with GEMIPM
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;  // pointer to DATABR

    // Creating TNode structure accessible trough node pointer
    // Here we read the files needed as input for initializing GEMIPM2K
    // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)

    m_Node = new TNode();

    dCH = m_Node->pCSD();
    dBR = m_Node->pCNode();

    if (Load_Init_File(Project_path, m_Node))
    {
        // The init file is successfully loaded
        // Getting direct access to DataCH structure in GEMIPM2K memory

        dBR->NodeStatusCH = NEED_GEM_AIA;
        m_Node->GEM_run(false);
        // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

        // The init file is successfully loaded
        // Getting direct access to DataCH structure in GEMIPM2K memory
        // iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

        // Extracting data bridge array sizes
        nIC = dCH->nICb;  // Num of Independent components
        nDC = dCH->nDCb;  // Num of Chemical species in the reactive part
        nPH = dCH->nPHb;  // Num of Phases
        nPS = dCH->nPSb;  // Num of multicomponent phases; ASSERT(nPS < nPH)

        // get the index of water
        idx_water = FindWater_xDC(m_Node);
        // imediately check
        if (idx_water == -1)
            return 1;
        // get index for H and OH
        idx_oxygen = Findoxygen_bIC(m_Node);
        // imediately check
        if (idx_oxygen == -1)
            return 1;
        idx_hydrogen = Findhydrogen_bIC(m_Node);
        // imediately check
        if (idx_hydrogen == -1)
            return 1;

        heatflag = GetHeatFlag_MT();  // Get heatflag
        flowflag = GetFlowType_MT();  // Get flow flag

        // get m_flow_pcs already, then check the flag:
        if (flag_coupling_hydrology == 1)
            // need to couple to flow process;
            // mark the flag
            m_flow_pcs->flag_couple_GEMS = 1;
        else
            m_flow_pcs->flag_couple_GEMS = 0;  // default is set to not-coupled!

        gem_pressure_flag = 0;  // use not a user defined pressure
        // Get number of Nodes
        nNodes = GetNodeNumber_MT();
        // Get number of Elems
        nElems = GetElemNumber_MT();
        // kg44 03 april 2010 first we take all variables that we need for all
        // nodes!

        // Allocating work memory for FMT part (here only chemical variables)
        m_NodeHandle = new long[nNodes];
        m_NodeHandle_buff = new long[nNodes];

        m_NodeStatusCH = new long[nNodes];
        m_NodeStatusCH_buff = new long[nNodes];

        m_IterDone = new long[nNodes];
        m_IterDone_buff = new long[nNodes];

        m_IterDoneCumulative = new long[nNodes];
        m_IterDoneIndex = new long[nNodes];

        // MPI Buffer Variable---------------

        m_boundary = new int[nNodes];  // this marks boundary nodes with fixed
                                       // concentrations!?
        m_T = new double[nNodes];
        m_P = new double[nNodes];

        m_Vs = new double[nNodes];
        m_Vs_buff = new double[nNodes];

        m_Ms = new double[nNodes];
        m_Ms_buff = new double[nNodes];

        m_Gs = new double[nNodes];
        m_Gs_buff = new double[nNodes];

        m_Hs = new double[nNodes];
        m_Hs_buff = new double[nNodes];

        m_IC = new double[nNodes];
        m_IC_buff = new double[nNodes];

        m_pH = new double[nNodes];
        m_pH_buff = new double[nNodes];

        m_pe = new double[nNodes];
        m_pe_buff = new double[nNodes];

        m_Eh = new double[nNodes];
        m_Eh_buff = new double[nNodes];

        m_porosity = new double[nNodes];
        m_porosity_buff = new double[nNodes];
        m_porosity_initial = new double[nNodes];
        m_volumes_initial = new double[nNodes * nPH];

        m_excess_water = new double[nNodes];
        m_excess_water_buff = new double[nNodes];

        m_excess_gas = new double[nNodes];
        m_excess_gas_buff = new double[nNodes];

        m_Node_Volume = new double[nNodes];

        m_fluid_volume = new double[nNodes];
        m_fluid_volume_buff = new double[nNodes];

        m_fluid_density = new double[nNodes];
        m_fluid_density_buff = new double[nNodes];

        m_gas_volume = new double[nNodes];
        m_gas_volume_buff = new double[nNodes];

        m_porosity_Elem = new double[nElems];
        m_porosity_Elem_buff = new double[nElems];

        m_soluteB = new double[nNodes * nIC];
        m_soluteB_buff = new double[nNodes * nIC];

        m_bIC = new double[nNodes * nIC];
        m_bIC_buff = new double[nNodes * nIC];

        m_bIC_dummy = new double[nNodes * nIC];
        m_bIC_dummy_buff = new double[nNodes * nIC];

        m_dul = new double[nNodes * nDC];
        m_dll = new double[nNodes * nDC];
        m_dul_buff = new double[nNodes * nDC];
        m_dll_buff = new double[nNodes * nDC];

        m_xDC = new double[nNodes * nDC];
        m_xDC_buff = new double[nNodes * nDC];

        m_aPH = new double[nNodes * nPH];  // surface area for surface species
                                           // ..input to GEMS!
        m_aPH_buff = new double[nNodes * nPH];
        m_xPH = new double[nNodes * nPH];  // amount of carrier...used for smart
                                           // initial aproximation
        m_xPH_buff = new double[nNodes * nPH];
        m_xPA = new double[nNodes * nPS];
        m_xPA_buff = new double[nNodes * nPS];

        m_bSP =
            new double[nNodes * nIC];  // Bulk composition of all solids, moles
                                       // [nIC] ...not yet buffered via MPI,
        // because we do not use the data yet
        // ----------------------------------
        // this is for kinetics
        Kinetic_GEMS m_kin;  // new kinetic vector

        omega_phase = new double[nNodes * nPH];
        omega_phase_buff = new double[nNodes * nPH];
        mol_phase = new double[nNodes * nPH];
        mol_phase_buff = new double[nNodes * nPH];
        dmdt = new double[nNodes * nPH];
        dmdt_buff = new double[nNodes * nPH];
        omega_components = new double[nNodes * nDC];
        omega_components_buff = new double[nNodes * nDC];

        m_xDC_pts = new double[nNodes * nDC];
        m_soluteB_pts = new double[nNodes * nIC];
        m_bIC_pts = new double[nNodes * nIC];
        m_xDC_MT_delta = new double[nNodes * nDC];
        m_xDC_Chem_delta = new double[nNodes * nDC];

        m_xDC_pts_buff = new double[nNodes * nDC];
        m_xDC_MT_delta_buff = new double[nNodes * nDC];
        m_xDC_Chem_delta_buff = new double[nNodes * nDC];

        // kg44 03 april 2010 ...from here on, most data is only necessary once
        // (check if this is needed for all nodes!)
        m_rMB =
            new double[nNodes * nIC];  // buffer removed...we do not need the
                                       // mass balance residuals globally
        m_uIC = new double[nNodes * nIC];  // chemical potentials...in current
                                           // code not needed, but maybe later?

        m_gam = new double[nNodes * nDC];

        m_vPS = new double[nNodes * nPS];
        m_mPS = new double[nNodes * nPS];
        m_bPS = new double[nNodes * nIC * nPS];
        m_bPS_buff = new double[nNodes * nIC * nPS];

        m_ICNL = new char[nIC][MaxICN];  // List of IC names in the system,
                                         // [nIC]  of MaxICN length
        m_DCNL = new char[nDC][MaxDCN];  // List of DC names in the system,
                                         // [nDC] of MaxDCN length
        m_PHNL =
            new char[nDC]
                    [MaxPHN];  // List of Phase names  [nPH]  of MaxPHN length

        // now fill the arrays for vtk output

        m_ICNL = dCH->ICNL;

        m_DCNL = dCH->DCNL;

        m_PHNL = dCH->PHNL;

        // ------------------------
        for (in = 0; in < nNodes; in++)
        {
            m_boundary[in] = 0;  // cout << "boundary init ok"<<"\n";
            m_NodeHandle[in] = 0;
            m_NodeStatusCH[in] = 0;
            m_IterDone[in] = 0;
            m_IterDoneCumulative[in] = 0;
            m_IterDoneIndex[in] = in;  // initial values order of Nodes
            m_NodeHandle_buff[in] = 0;
            m_NodeStatusCH_buff[in] = 0;
            m_IterDone_buff[in] = 0;

            m_T[in] = 298.15;  // equivalent to 25°C
            m_P[in] = 1.0e+5;
            m_Vs[in] = 0.0;
            m_Ms[in] = 0.0;
            m_Gs[in] = 0.0;
            m_Hs[in] = 0.0;
            m_IC[in] = 0.0;
            m_pH[in] = 0.0;
            m_pe[in] = 0.0;
            m_Eh[in] = 0.0;
            m_porosity[in] = 0.0;
            m_porosity_initial[in] = 0.0;
            m_fluid_volume[in] = 0.0;
            m_fluid_density[in] = 1000.0;
            m_gas_volume[in] = 0.0;

            m_Vs_buff[in] = 0.0;
            m_Ms_buff[in] = 0.0;
            m_Gs_buff[in] = 0.0;
            m_Hs_buff[in] = 0.0;
            m_IC_buff[in] = 0.0;
            m_pH_buff[in] = 0.0;
            m_pe_buff[in] = 0.0;
            m_Eh_buff[in] = 0.0;
            m_porosity_buff[in] = 0.0;
            m_fluid_volume_buff[in] = 0.0;
            m_fluid_density_buff[in] = 0.0;
            m_gas_volume_buff[in] = 0.0;

            m_excess_water[in] = 0.0;
            m_excess_gas[in] = 0.0;

            m_Node_Volume[in] = REACT_GEM::GetNodeAdjacentVolume(in);
            m_excess_water_buff[in] = 0.0;
            m_excess_gas_buff[in] = 0.0;

            for (ii = 0; ii < nIC; ii++)
            {
                m_bIC[in * nIC + ii] = 0.0;
                m_bIC_dummy[in * nIC + ii] = 0.0;

                m_soluteB[in * nIC + ii] = 0.0;
                m_soluteB_buff[in * nIC + ii] = 0.0;
                m_soluteB_pts[in * nIC + ii] = 0.0;
                m_bIC_pts[in * nIC + ii] = 0.0;
                m_rMB[in * nIC + ii] = 0.0;
                m_uIC[in * nIC + ii] = 0.0;

                m_bIC_buff[in * nIC + ii] = 0.0;
                m_bIC_dummy_buff[in * nIC + ii] = 0.0;
            }

            for (ii = 0; ii < nDC; ii++)
            {
                m_xDC[in * nDC + ii] = 0.0;
                m_gam[in * nDC + ii] = 0.0;
                m_dul[in * nDC + ii] =
                    1.0e+10;  // this should be a large number, because after
                              // scaling to 1kg in Gems it should be 1.e+6
                m_dll[in * nDC + ii] = 0.0;  // zero is ok
                m_xDC_pts[in * nDC + ii] = 0.0;
                m_xDC_MT_delta[in * nDC + ii] = 0.0;
                m_xDC_Chem_delta[in * nDC + ii] = 0.0;

                m_xDC_buff[in * nDC + ii] = 0.0;
                m_dul_buff[in * nDC + ii] = 0.0;
                m_dll_buff[in * nDC + ii] = 0.0;
                m_xDC_pts_buff[in * nDC + ii] = 0.0;
                m_xDC_MT_delta_buff[in * nDC + ii] = 0.0;
                m_xDC_Chem_delta_buff[in * nDC + ii] = 0.0;
                omega_components[in * nDC + ii] = 0.0;
                omega_components_buff[in * nDC + ii] = 0.0;
            }

            for (ii = 0; ii < nPH; ii++)
            {
                m_aPH[in * nPH + ii] = 0.0;
                m_xPH[in * nPH + ii] = 0.0;

                m_xPH_buff[in * nPH + ii] = 0.0;
                m_aPH_buff[in * nPH + ii] = 0.0;

                omega_phase[in * nPH + ii] = 0.0;
                omega_phase_buff[in * nPH + ii] = 0.0;
                mol_phase[in * nPH + ii] = 0.0;
                mol_phase_buff[in * nPH + ii] = 0.0;

                dmdt[in * nPH + ii] = 0.0;
                dmdt_buff[in * nPH + ii] = 0.0;
                m_volumes_initial[in * nPH + ii] = 0.0;
            }

            for (ii = 0; ii < nPS; ii++)
            {
                m_vPS[in * nPS + ii] = 0.0;
                m_mPS[in * nPS + ii] = 0.0;
                m_xPA[in * nPS + ii] = 0.0;
                m_xPA_buff[in * nPS + ii] = 0.0;
            }

            for (ii = 0; ii < nIC; ii++)
                for (int jj = 0; jj < nPS; jj++)
                {
                    m_bPS[in * ii * nPS + jj] = 0.0;
                    m_bPS_buff[in * ii * nPS + jj] = 0.0;
                }
        }

        for (in = 0; in < nElems; in++)
        {
            m_porosity_Elem[in] = 0.0;
            m_porosity_Elem_buff[in] = 0.0;
        }
        nNodes = GetNodeNumber_MT();
        nElems = GetElemNumber_MT();
        if (nPH < 2 && flowflag == 3)
        {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            if (myrank == 0 /*should be set to root*/)
#endif
                cout << "Richards flow used and only the water pahse is "
                        "defined in GEMS...please add a gas phase.\n";
            exit(1);
        }

        string tinit_path = Project_path;
        //	  if ( myrank == 0 /*should be set to root*/ )
        {
            // here the boost barriers are defined....
            gem_barrier_finish = (new boost::barrier((gem_nThread + 1)));
            gem_barrier_start = (new boost::barrier((gem_nThread + 1)));
            gemThread =
                (new boost::thread[gem_nThread]);  // each mpi task creates
                                                   // gem_nThread workers

            for (i = 0; i < gem_nThread; ++i)  // here we create the threads!
            {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                rwmutex.lock();
                cout << "Creating GEMS-worker thread " << i << " for MPI task "
                     << myrank << "\n";
                rwmutex.unlock();
#else
                rwmutex.lock();
                cout << "Creating GEMS-worker thread " << i << "\n";
                rwmutex.unlock();
#endif
                gemThread[i] = boost::thread(boost::bind(
                    &REACT_GEM::gems_worker, this, (int)i, Project_path));
            }
        }
        return 0;  // successed
    }
    return 1;  // something went wrong
}

short REACT_GEM::Init_RUN(string Project_path)
{
    nNodes = GetNodeNumber_MT();
    nElems = GetElemNumber_MT();
    long ii = 0, in = 0, i = 0, j;
    //        CompProperties *m_cp = NULL;
    CRFProcess* this_pcs;
    double BCValue = 0.0;
    string cstr;

    rwmutex.lock();  // make sure we do not interfere with the threads

    TNode* m_Node;
    // DATABR structure for exchange with GEMIPM
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;  // pointer to DATABR

    // Creating TNode structure accessible trough node pointer
    m_Node = new TNode();
    // Here we read the files needed as input for initializing GEMIPM2K
    // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)
    if (Load_Init_File(Project_path, m_Node))
    {
        // The init file is successfully loaded
        // Getting direct access to DataCH structure in GEMIPM2K memory
        dCH = m_Node->pCSD();
        if (!dCH)
        {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            MPI_Finalize();  // make sure MPI exits
#endif
            exit(1);
        }

        // Getting direct access to work node DATABR structure which
        // exchanges data between GEMIPM and FMT parts
        dBR = m_Node->pCNode();
        if (!dBR)
        {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            MPI_Finalize();  // make sure MPI exits
#endif
            exit(1);
        }
        // run GEMS once
        dBR->NodeStatusCH = NEED_GEM_AIA;
        m_Node->GEM_run(false);
    }
    else
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
        return 5;
    }

    for (i = 0; i < nNodes;
         i++)  // after running GEMS once we initialize the OGS-GEMS data arrays

        // initialize the arrays
        REACT_GEM::GetReactInfoFromGEM(
            i, m_Node);  // get the data even if GEMS failed...this is necessary
                         // for
    // initializing the arrays for example m_aPH
    // unfortunately aPH is passed to GEMS and needs correct input for working
    // with sorption

    for (ii = 0; ii < (int)m_kin.size();
         ii++)  // this loop is for identifying kinetically controlled phases
    {
        m_kin[ii].phase_number = -1;
        m_kin[ii].dc_counter = 0;  // this is the starting dc
        // phase numbers are not yet set...do it!
        for (j = 0; j < nPH; j++)
        {
            cstr.assign(dCH->PHNL[j]);
            //   cout << m_kin[ii].phase_name  << " gems phase: " << cstr << "
            //   "<< " component start number " << m_kin[ii].dc_counter << "\n";
            if (m_kin[ii].phase_name == cstr)
            {
                m_kin[ii].phase_number = j;
                break;
            }
            // add the number of components to the counter..after test for
            // phasename...as this counter gives the starting position of the
            // dependent components
            m_kin[ii].dc_counter += dCH->nDCinPH[j];
        }

        if (m_kin[ii].phase_number < 0 || m_kin[ii].phase_number >= nPH)
        {
            cout << " GEMS: Error in Phase kinetics..check input for "
                 << m_kin[ii].phase_name << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            MPI_Finalize();  // make sure MPI exits
#endif
            exit(1);
        }
        else
            cout << "GEM Kinetics phase number:  " << m_kin[ii].phase_name
                 << " in phase number " << m_kin[ii].phase_number << "\n";
    }
    // Marking species with fixed concentrations (in boundary nodes)
    // first we check the boundary nodes with fixed concentrations
    // this is adopted from rf_REACT_BRNS...we only look for the first species,
    // as with GEMS we should define boundary conditions for ALL species
    this_pcs = NULL;
    // Get the pointer to the proper PCS.
    this_pcs = PCSGet("MASS_TRANSPORT");
    for (i = 0; i < nNodes; i++)
    {
        if (this_pcs)
        {
            // BC printing
            if (IsThisPointBCIfYesStoreValue(i, this_pcs, BCValue))
                // If this node is on the fixed boudnary for this component
                //                cout << "Node " << i <<", Comp " <<
                //                this_pcs->pcs_primary_function_name[0] << "
                //                ,Value " << BCValue << " is boundary node"
                //                <<"\n";
                m_boundary[i] = 1;
            else
                // If this node is NOT on the fixed boudnary for this component
                m_boundary[i] = 0;
        }
        else  // not getting the pointer to the proper PCS.
        {
            cout << this_pcs->pcs_primary_function_name[0]
                 << "!!! In InitGEMS, can not find corresponding PCS for "
                    "checking boundary conditions! "
                 << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            MPI_Finalize();  // make sure MPI exits
#endif
            exit(1);
        }
    }  // end loop over all nodes
    // here we have to make sure that the flow process has all initial data to
    // do the saturation! synchronize pressures should be enough!

    flowflag = GetFlowType_MT();  // Get flow flag

    // we have to make sure saturation is properly defined..otherwise we can not
    // calculate concentrations properly kg44 16.05.2013 something wrong here?
    // saturation seems not properly defined for Richards flow and restart
    if (flowflag == 3)  // this is only for Richards flow
    {
        // TODO for petsc probably need to synchronize the pressure field here
        // (as shadow nodes are not initializee via surfaces/volumes/plylines
        // etc)
        cout << "GEM-INIT: ";
        if (m_flow_pcs->saturation_switch == true)
        {
            cout << "CalcSaturationRichards "
                 << "\n";
            // is true here correct?
            m_flow_pcs->CalcSaturationRichards(1, true);
        }  // JOD
        else
        {
            // WW
            m_flow_pcs->CalcSecondaryVariablesUnsaturatedFlow(true);
            cout << " CalcSecondaryVariablesUnsaturatedFlow"
                 << "\n";
        }

        // TODO for petsc we may need to synchronize saturations ...have to
        // check!
    }

    // get the restart data specific for gems
    if (m_flow_pcs->GetRestartFlag() == FiniteElement::READ ||
        m_flow_pcs->GetRestartFlag() == FiniteElement::READ_WRITE)
    {
        if (!ReadReloadGem())
        {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            MPI_Finalize();  // make sure MPI exits
#endif

            exit(1);  // we save m_bic and m_soluteb as is, therefore these are
                      // concentrations -> call
                      // concentration_to:mass
        }
        /*       for ( i=0 ; i < nNodes ; i++ )  Better not here.....
           {
             REACT_GEM::ConcentrationToMass ( i,1); // now it should be
           identical to the normal start: total b vectors!
           }  */
    }

    //    else  test: we do this always! ...should be safe ;-)
    {
        GetInitialReactInfoFromMassTransport(
            1);  // get the initial values from MD ...IC are total B vectors -
                 // last time
        // step! this is not necessary for restart
        cout << "Attentione GEMS users: Initial kinetics calculated without "
                "restart! This probably kills kinetics, as "
                "phases in m_xDC are not yet properly initialized!"
             << "\n";
        cout << "No upper or lower constrains set during equilibration!...If "
                "your setup requires constrains, please "
                "contact georg.kosakowski@psi.ch"
             << "\n";
    }

    delete m_Node;     // get rid of this GEMS instance
    rwmutex.unlock();  // now it is save to release the lock
// from here the gems threads are responsible for GEMS

// distribute the data
#ifdef USE_MPI_GEMS
    // MPI initialization.
    MPI_Bcast(&nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
    //  if ( myrank == 0 /*should be set to root*/ )
    {
        //    cout << "main waits for start barrier " << "\n";

        gem_barrier_start->wait();  // give the workers the start signal
        //    cout << "main passed start barrier " << "\n";

        gem_barrier_finish
            ->wait();  // wait for the workers to finish the initial stuff
    }
#if defined(USE_MPI_GEMS)
    // For MPI scheme, gather the data here.
    REACT_GEM::GetGEMResult_MPI();
    REACT_GEM::CleanMPIBuffer();
#endif

    for (in = 0; in < nNodes; in++)  // set the correct boundary conditions
    {
        REACT_GEM::SetBValue_MT(in, 0, &(m_soluteB[in * nIC]));  // old timestep
        REACT_GEM::SetBValue_MT(in, 1, &(m_soluteB[in * nIC]));  // new timestep
    }

    // always calculate porosity at the beginning
    ConvPorosityNodeValue2Elem(
        1);  // new timestep: update element porosity and push back values
    ConvPorosityNodeValue2Elem(0);  // old timestep: fill with current values!
    CopyCurBPre();
    rwmutex.lock();
    cout << "Initial Running GEM  to get the correct porosities successful. "
         << "\n";
    rwmutex.unlock();

    return 0;
}

string REACT_GEM::Get_Init_File_Path(void)
{
    return init_input_file_path;
}

string REACT_GEM::Get_IPM_File_Path(void)
{
    return ipm_input_file_path;
}

string REACT_GEM::Get_DBR_File_Path(void)
{
    return dbr_input_file_path;
}

string REACT_GEM::Get_DCH_File_Path(void)
{
    return dch_input_file_path;
}

int REACT_GEM::Set_IPM_FILE_PATH(string m_path)
{
    REACT_GEM::ipm_input_file_path = m_path;
    return 0;
}

int REACT_GEM::Set_DBR_FILE_PATH(string m_path)
{
    REACT_GEM::dbr_input_file_path = m_path;
    return 0;
}

int REACT_GEM::Set_DCH_FILE_PATH(string m_path)
{
    REACT_GEM::dch_input_file_path = m_path;
    return 0;
}

int REACT_GEM::Set_Init_File_Path(string m_path)
{
    REACT_GEM::init_input_file_path = m_path;
    return 0;
}

bool REACT_GEM::Load_Init_File(string m_Project_path, TNode* m_Node)
{
    string init_path;
    char* buffer = NULL;

    init_path = m_Project_path.append(REACT_GEM::init_input_file_path);

#ifdef _WIN32
    // keep this on windows
    if (init_path.rfind("\\") == string::npos)
#else
    // keep this on linux
    if (init_path.rfind("/") == string::npos)
#endif
    {
#ifdef _WIN32
        if ((buffer = _getcwd(NULL, 0)) == NULL)
#else
        if ((buffer = getcwd(NULL, 0)) == NULL)
#endif
            perror("_getcwd error");
        else
        {
#ifdef _WIN32
            init_path.insert(0, "\\");  // keep this on window
#else
            init_path.insert(0, "/");  // keep this on linux
#endif
            init_path.insert(0, buffer);
        }
    }

    if (buffer)
        free(buffer);

    if (m_Node->GEM_init(init_path.c_str()))
        return 0;  // error occured during reading the files
    else
        return 1;  // read init file successed
}

short REACT_GEM::GetInitialReactInfoFromMassTransport(int timelevel)
{
    heatflag = GetHeatFlag_MT();
    flowflag = GetFlowType_MT();
    REACT_GEM::nNodes = GetNodeNumber_MT();

    for (long node_i = 0; node_i < nNodes; node_i++)
    {
        // get temperature from MT
        m_T[node_i] = REACT_GEM::GetTempValue_MT(node_i, timelevel);

        // get pressure from MT
        m_P[node_i] = REACT_GEM::GetPressureValue_MT(node_i, timelevel);
        // get Independent and dependent Component value from MT
        if ((flag_transport_b == 1) &&
            (m_flow_pcs->GetRestartFlag() == FiniteElement::WRITE ||
             m_flow_pcs->GetRestartFlag() == FiniteElement::NO_IO))
            REACT_GEM::GetBValue_MT(
                node_i, timelevel,
                m_bIC +
                    node_i *
                        nIC);  // do this not for restart...overwrites values!!!
    }
#if defined(USE_PETSC)
    // arrays are filled, now we should synchronize the values
    SynchronizeData(
        m_P);  // this also overwrites default values for shadow nodes!
    SynchronizeData(
        m_T);  // this also overwrites default values for shadow nodes!
    // for b vector we need to copy data into working vector first...with use
    // the b_buffer for this_pcs in is node...i is vector component!
    for (long i = 0; i < nIC; i++)
    {
        for (long in = 0; in < loc_NodesNumber_Linear; in++)
            m_bIC_buff[in] = m_bIC[in * nIC + i];  // copy only local nodes
        SynchronizeData(m_bIC_buff);
        for (long in = 0; in < nNodes; in++)
            m_bIC[in * nIC + i] = m_bIC_buff[in];  // copy also shadow nodes
    }
#endif
    return 0;
}

short REACT_GEM::GetReactInfoFromMassTransport(int timelevel)
{
    heatflag = GetHeatFlag_MT();
    flowflag = GetFlowType_MT();
    REACT_GEM::nNodes = GetNodeNumber_MT();

    for (long node_i = 0; node_i < nNodes; node_i++)
    {
        // get temperature from MT
        m_T[node_i] = REACT_GEM::GetTempValue_MT(node_i, timelevel);

        // get pressure from MT
        m_P[node_i] = REACT_GEM::GetPressureValue_MT(node_i, timelevel);
        // get Independent and dependent Component value from MT
        if (flag_transport_b == 1)
            REACT_GEM::GetBValue_MT(node_i, timelevel,
                                    m_soluteB + node_i * nIC);
        // Convert to mole values
        // if(conv_concentration == 1)	REACT_GEM::ConcentrationToMass( node_i
        // );

        // Setting Solid Phase Component // HS: Solid does not move.
        // REACT_GEM::GetSoComponentValue_MT(node_i, timelevel, m_xPH+node_i*nPH
        // );
    }

    return 0;
}

short REACT_GEM::SetReactInfoBackMassTransport(int timelevel)
{
    heatflag = GetHeatFlag_MT();
    flowflag = GetFlowType_MT();
    REACT_GEM::nNodes = GetNodeNumber_MT();

    for (long in = 0; in < nNodes; in++)
    {
        // Setting Temperature // disabled by HS. temperature is NOT the output
        // from chemistry. REACT_GEM::SetTempValue_MT(in,timelevel,m_T[in]);

        // Setting Pressure // disabled by HS. pressure is NOT the output from
        // chemistry. REACT_GEM::SetPressureValue_MT(in,timelevel,m_P[in]);

        // if (m_pcs->m_msh->nod_vector[in]->onBoundary() == false) {
        // Setting Independent Component
        if (m_NodeStatusCH[in] == OK_GEM_AIA ||
            m_NodeStatusCH[in] == OK_GEM_SIA ||
            m_NodeStatusCH[in] == BAD_GEM_AIA ||
            m_NodeStatusCH[in] == BAD_GEM_SIA)
            if (flag_transport_b == 1)
                REACT_GEM::SetBValue_MT(in, timelevel, &(m_soluteB[in * nIC]));

        // Set the extra water as source/sink term; not for boundary nodes

        if (flag_coupling_hydrology > 0 && !m_boundary[in])
            REACT_GEM::SetSourceSink_MT(in, dt /*in sec*/);
    }
#if defined(USE_MPI_GEMS)
    if (flag_coupling_hydrology > 0)
        m_flow_pcs->SetSTWaterGemSubDomain(
            myrank);  // necessary for domain decomposition

#endif
    if (flag_porosity_change > 0)
        ConvPorosityNodeValue2Elem(0);  // old timestep :copy current values to
                                        // old timestep before updating porosity
    if (flag_porosity_change > 0)
        ConvPorosityNodeValue2Elem(timelevel);  // new timestep :update element
                                                // porosity and push back values
    return 0;
}

void REACT_GEM::GetReactInfoFromGEM(long in, TNode* m_Node)
{
    m_Node->GEM_to_MT(m_NodeHandle[in],
                      m_NodeStatusCH[in],
                      m_IterDone[in],
                      m_Vs[in],
                      m_Ms[in],
                      m_Gs[in],
                      m_Hs[in],
                      m_IC[in],
                      m_pH[in],
                      m_pe[in],
                      m_Eh[in],
                      m_rMB + in * nIC,
                      m_uIC + in * nIC,
                      m_xDC + in * nDC,
                      m_gam + in * nDC,
                      m_xPH + in * nPH,
                      m_vPS + in * nPS,
                      m_mPS + in * nPS,
                      m_bPS + in * nPS * nIC,
                      m_xPA + in * nPS,
                      m_aPH + in * nPH,
                      m_bSP + in * nIC);
    // old kernel version m_bPS+in*nPS*nIC, m_xPA+in*nPS, m_aPH+in*nPH );
}

void REACT_GEM::SetReactInfoBackGEM(long in, TNode* m_Node)
{
    // Setting input data for GEMIPM

    //  for (i=0;i<nIC;i++) cout << m_bIC[in*nIC+i] << "\n";

    if (flag_transport_b == 1)  // here we insert the actual B vector

        m_Node->GEM_from_MT(m_NodeHandle[in],
                            m_NodeStatusCH[in],
                            m_T[in],
                            m_P[in],
                            m_Vs[in],
                            m_Ms[in],
                            m_bIC + in * nIC,
                            m_dul + in * nDC,
                            m_dll + in * nDC,
                            m_aPH + in * nPH);
    //	cout << m_xDC+in*nDC << "\n";
    // set charge to zero
    m_Node->pCNode()->bIC[nIC - 1] = 0.0;
}

short REACT_GEM::Run_MainLoop()
{
    if (flag_disable_gems)
        return 0;  // do nothing if GEMS calculations are disabled
    max_kinetic_timestep = 1.0e+99;  // restrict time step for kinetics

#ifdef USE_MPI_GEMS
    // MPI initialization.
    // So here is going to distribute the task.
    MPI_Bcast(&nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
    //  if ( myrank == 0 /*should be set to root*/ )
    {
        gem_barrier_start->wait();
        // cout << "main "<< " passed for start barrier "<<"\n";

        // gem_condition.notify_all();
        // check mutexes if all threads are finished....better use barrier?

        // wait for barrier
        //  cout << "main "<< " waiting for finish barrier "<<"\n";
        gem_barrier_finish->wait();
        // cout << "main "<< " passed for finish barrier "<<"\n";
    }
#ifdef USE_MPI_GEMS
    // For MPI scheme, gather the data here.
    REACT_GEM::GetGEMResult_MPI();
    REACT_GEM::CleanMPIBuffer();
#endif
    rwmutex.lock();  // avoid mutual exclusion in the MPI version
    //    cout << "DEBUG failed: max kinetic time step " << max_kinetic_timestep
    //    << "\n"; cout << " GEM  run successful. "  << "\n";
    rwmutex.unlock();

    return 0;
}

int REACT_GEM::GetHeatFlag_MT(void)
{
    CRFProcess* m_pcs = NULL;
    // heat transport
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //                if ( m_pcs->pcs_type_name.compare ( "HEAT_TRANSPORT" )
        //                == 0 ) {
        // TF
        if (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)
            return 1;
    }
    return 0;
}

int REACT_GEM::GetFlowType_MT(void)
{
    CRFProcess* m_pcs = NULL;
    // flow type
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //                if ( m_pcs->pcs_type_name.compare ( "GROUNDWATER_FLOW"
        //                ) ==0 ) {
        if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
        {
            m_flow_pcs = m_pcs;
            return 1;
            //                } else if ( m_pcs->pcs_type_name.compare (
            //                "LIQUID_FLOW" ) ==0 ) {
        }
        else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
        {
            m_flow_pcs = m_pcs;
            return 2;
            //                } else if ( m_pcs->pcs_type_name.compare (
            //                "RICHARDS_FLOW" ) ==0 ) {
        }
        else if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
        {
            m_flow_pcs = m_pcs;
            return 3;
            //                } else if ( m_pcs->pcs_type_name.compare (
            //                "MULTI_PHASE_FLOW" ) ==0 ) {
        }
        else if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        {
        }
        m_flow_pcs = m_pcs;
        return 4;
    }
    return 0;
}

long REACT_GEM::GetNodeNumber_MT(void)
{
    long number;
    CRFProcess* m_pcs = NULL;
    //------------read number of nodes--------------
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //		if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            number = (long)m_pcs->m_msh->GetNodesNumber(false);
            return number;
        }
    }
    //------------end of reading number of nodes----
    return 0;
}

#if defined(USE_PETSC)
long REACT_GEM::GetGlobalNodeNumber_MT(void)
{
    long number;
    CRFProcess* m_pcs = NULL;
    //------------read number of nodes--------------
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //		if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            number = (long)m_pcs->m_msh->getNumNodesGlobal();
            return number;
        }
    }
    //------------end of reading number of nodes----
    return 0;
}
long REACT_GEM::GetLocalNodeNumber_MT(void)
{
    long number;
    CRFProcess* m_pcs = NULL;
    //------------read number of nodes--------------
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //		if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            number = (long)m_pcs->m_msh->getNumNodesLocal();
            return number;
        }
    }
    //------------end of reading number of nodes----
    return 0;
}
#endif

long REACT_GEM::GetElemNumber_MT(void)
{
    long number;
    CRFProcess* m_pcs = NULL;
    //------------read number of elems--------------
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
        //                ==0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            number = (long)m_pcs->m_msh->ele_vector.size();
            return number;
        }
    }
    return 0;
}

double REACT_GEM::GetTempValue_MT(long node_Index, int timelevel)
{
    int indx;
    double temp;
    CRFProcess* m_pcs = NULL;

    if (heatflag == 1)
    {
        m_pcs = PCSGet("HEAT_TRANSPORT");

        indx = m_pcs->GetNodeValueIndex("TEMPERATURE1") + timelevel;
        temp = m_pcs->GetNodeValue(node_Index, indx);
    }
    else
        temp = m_gem_temperature;
    return temp;
}
short REACT_GEM::SetTempValue_MT(long node_Index, int timelevel, double temp)
{
    int indx;
    CRFProcess* m_pcs = NULL;
    if (heatflag == 1)
    {
        m_pcs = PCSGet("HEAT_TRANSPORT");

        indx = m_pcs->GetNodeValueIndex("TEMPERATURE1") + timelevel;
        m_pcs->SetNodeValue(node_Index, indx, temp);

        // sysT[i] = m_pcs->GetNodeValue(i, indx1);
        // if (sysT0[i] <273.15) sysT0[i] += 273.15;  //ToDo �C->K
        // if (sysT[i] <273.15) sysT[i] += 273.15;  //ToDo �C->K
        return 1;
    }
    else
        return 0;
}

double REACT_GEM::GetPressureValue_MT(long node_Index, int timelevel)
{
    // Get pressure value
    double pressure;
    int indx;
    pressure = 0.0;

    CFluidProperties* m_FluidProp;
    m_FluidProp = MFPGet("LIQUID");

    if (flowflag > 0)
    {
        switch (flowflag)
        {
            case 1:  // for "GROUNDWATER_FLOW";

                if (gem_pressure_flag < 1)
                {
                    // just set to 1.0 bar.
                    pressure = m_gem_pressure;
                    break;
                }
                else
                {
                    indx = m_flow_pcs->GetNodeValueIndex("HEAD");
                    // The unit of HEAD is in meters
                    pressure =
                        m_flow_pcs->GetNodeValue(node_Index, indx + timelevel);
                    // cout << pressure << " " << " timelevel "<<timelevel;
                    // change the pressure unit from hydraulic head to Pa.
                    pressure =
                        Pressure_M_2_Pa(pressure, m_FluidProp->Density());
                    // y cout << pressure << " density" <<
                    // m_FluidProp->Density()<<"\n"; add atmospheric pressure
                }

                if (pressure <= 0.0 /*valcumm suction in groundwater is not so
                                       realistic*/
                    || pressure > 1.0e+15 /*some very high pressure*/
                )
                {
                    // then set it to 1.0 bar = 1.0e5 Pa;
                    cout << " high pressure " << pressure << "\n";
                    pressure = 1.0e+05;
                }
                break;

            case 2:  // for "LIQUID_FLOW", not tested!!!

                indx = m_flow_pcs->GetNodeValueIndex("PRESSURE1") + timelevel;
                // The unit of HEAD is in meters
                pressure = m_flow_pcs->GetNodeValue(node_Index, indx);

                // change the pressure unit from meters of water to bar.
                // pressure = Pressure_M_2_Bar ( pressure ,
                // m_FluidProp->Density() ); add atmospheric pressure
                if (pressure < 0.0 /*valcumm suction in groundwater is not so
                                      realistic*/
                    || pressure > 1.0e+15 /*some very high pressure*/
                )
                    pressure = 1.0e+5;  // then set it to 1.0 bar;
                break;
            case 3:  // for "RICHARDS_FLOW", not tested!!!
                pressure = m_gem_pressure;
                //	indx = m_flow_pcs->GetNodeValueIndex ( "PRESSURE1" )
                //+timelevel; 	pressure = m_flow_pcs->GetNodeValue ( node_Index,
                //indx ); // The unit of HEAD is in meters

                // change the pressure unit from meters of water to bar.
                //	pressure = Pressure_M_2_Bar ( pressure ,
                //m_FluidProp->Density() );
                // add atmospheric pressure
                //	pressure +=1.0;
                if (pressure < 0.0 /*valcumm suction in groundwater is not so
                                      realistic*/
                    || pressure > 1.0e+15 /*some very high pressure*/
                )
                    pressure = 1.0e+5;  // then set it to 1.0 bar;
                break;
            case 4:  // MULTIPHASE ....not tested
                indx = m_flow_pcs->GetNodeValueIndex("PRESSURE1");
                // The unit of HEAD is in meters
                pressure =
                    m_flow_pcs->GetNodeValue(node_Index, indx + timelevel);

                // change the pressure unit from meters of water to bar.
                // pressure = Pressure_M_2_Bar ( pressure ,
                // m_FluidProp->Density() ); add atmospheric pressure

                if (pressure < 0.0 /*valcumm suction in groundwater is not so
                                      realistic*/
                    || pressure > 1.0e+15 /*some very high pressure*/
                )
                    pressure = 1.0e+5;  // then set it to 1.0 bar;
                break;
            default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                if (myrank == 0 /*should be set to root*/)
#endif
                    cout << "Error: Not implemented for the flow in GEM case!!!"
                         << "\n";
                pressure = 1.0e+05;
                break;
        }  // end of switch case;
    }      // end of if (flow_flag);
    else
        // if no valid flow pcs existing;
        cout << "Warning: No valid flow process!!"
             << "\n";
    return pressure;
}

short REACT_GEM::SetPressureValue_MT(long node_Index, int timelevel,
                                     double pressure)
{
    CRFProcess* m_pcs = NULL;
    CFluidProperties* m_FluidProp;
    m_FluidProp = MFPGet("LIQUID");

    // Set pressure value
    int indx;
    indx = 0;
    if (flowflag > 0)
    {
        switch (flowflag)
        {
            case 1:
                m_pcs = PCSGet("GROUNDWATER_FLOW");
                pressure = Pressure_Pa_2_M(pressure, m_FluidProp->Density());
                indx = m_pcs->GetNodeValueIndex("HEAD") + timelevel;
                m_pcs->SetNodeValue(node_Index, indx, pressure);
                break;
            case 2:
                m_pcs = PCSGet("LIQUID_FLOW");
                indx = m_pcs->GetNodeValueIndex("PRESSURE1") + timelevel;

                m_pcs->SetNodeValue(node_Index, indx, pressure);
                break;
            case 3:
                // do nothing for Richards flow
                break;
            //				m_pcs = PCSGet ( "RICHARDS_FLOW" );
            //				indx = m_pcs->GetNodeValueIndex ( "PRESSURE1" )
            //+timelevel; 				pressure = Pressure_Bar_2_Pa ( pressure );
            //				m_pcs->SetNodeValue ( node_Index, indx, pressure );
            //				break;
            case 4:
                m_pcs = PCSGet("MULTI_PHASE_FLOW");
                indx = m_pcs->GetNodeValueIndex("PRESSURE1") + timelevel;

                m_pcs->SetNodeValue(node_Index, indx, pressure);
            default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                if (myrank == 0 /*should be set to root*/)
#endif
                    cout << "Error: Not implemented for the flow in GEM case!!!"
                         << "\n";
                break;
        }
    }
    else
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        if (myrank == 0 /*should be set to root*/)
#endif
            cout << "Warning: No valid flow process!!"
                 << "\n";
        return 0;
    }
    return 1;
}

short REACT_GEM::GetDCValue_MT(long node_Index, int timelevel, double* m_DC,
                               double* m_DC_pts, double* m_DC_MT_delta)
{
    string str;
    double /*DC_MT_pre,*/ DC_MT_cur;
    CRFProcess* m_pcs = NULL;
    int i = -1;

    for (size_t j = 0; j < pcs_vector.size(); j++)
    {
        m_pcs = pcs_vector[j];  // dangerous!!
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            // if ( m_pcs->m_msh->nod_vector[node_Index]->onBoundary() == false
            // ) // do not update values for boundary node?
            i += +1;
            str = m_pcs->pcs_primary_function_name[0];

            // Get previous iteration mass transport concentration value
            // DC_MT_pre = m_pcs->GetNodeValue (
            // node_Index,m_pcs->GetNodeValueIndex ( str ) +0 ); Get current
            // iteration mass transport concentration value
            DC_MT_cur = m_pcs->GetNodeValue(
                node_Index, m_pcs->GetNodeValueIndex(str) + timelevel);
            *(m_DC + i) = DC_MT_cur;
        }
    }

    return 1;
}

short REACT_GEM::GetBValue_MT(long node_Index, int timelevel, double* m_soluteB)
{
    string str;
    double /*DC_MT_pre,*/ DC_B_cur;
    CRFProcess* m_pcs = NULL;
    long i = -1;  // set it to minus 1 as we add 1 before we use it
    for (size_t j = 0; j < pcs_vector.size(); j++)
    {
        m_pcs = pcs_vector[j];  // dangerous!!
        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
        //                == 0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            // if ( m_pcs->m_msh->nod_vector[node_Index]->onBoundary() == false
            // ) // do not update values for boundary node?

            str = m_pcs->pcs_primary_function_name[0];
            // here we could add some type of test in order to make sure we get
            // the correct process!
            i += 1;  // add one to counter
            // Get previous iteration mass transport concentration value
            // DC_MT_pre = m_pcs->GetNodeValue (
            // node_Index,m_pcs->GetNodeValueIndex ( str ) +0 ); Get current
            // iteration mass transport concentration value
            DC_B_cur = m_pcs->GetNodeValue(
                node_Index, m_pcs->GetNodeValueIndex(str) + timelevel);

            *(m_soluteB + i) = DC_B_cur;
        }
    }

    return 1;
}

// same functionality as GetComponentValue_MT
double REACT_GEM::GetDCValueSpecies_MT(long node_Index, int timelevel, int iDc)
{
    string str;
    double /*DC_MT_pre,*/ DC_MT_cur = 0.0;
    size_t i = 0;  // counter for processes
    CRFProcess* m_pcs = NULL;
    // first find out the number for the first Mass_transport process
    while (i < pcs_vector.size())
    {
        m_pcs = pcs_vector[i];  // dangerous!!
        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
        //                == 0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            m_pcs = pcs_vector[iDc + i];  // dangerous!! ... this is now the
                                          // process we are looking for
            //        if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==
            //        0 ) {
            if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
            {
                str = m_pcs->pcs_primary_function_name[0];
                if (str.compare("pH") != 0 && str.compare("pe") != 0 &&
                    str.compare("Eh") != 0 && str.compare("NodePorosity") != 0)
                {
                    // Get previous iteration mass transport concentration value
                    // DC_MT_pre = m_pcs->GetNodeValue (
                    // node_Index,m_pcs->GetNodeValueIndex ( str ) +0 ); Get
                    // current iteration mass transport concentration value
                    getnode_mutex.lock();
                    DC_MT_cur = m_pcs->GetNodeValue(
                        node_Index,
                        m_pcs->GetNodeValueIndex(str) +
                            timelevel);  // KG44 I am not quite
                    // sure, but it looks as
                    // this is not thread
                    // safe...lets look it
                    // with a special mutex
                    getnode_mutex.unlock();
                }
                else
                {
                    cout
                        << "Error in GetDCValueSpecies_MT ... return zero value"
                        << "\n";
                    DC_MT_cur = 0.0;
                }
            }
            else
            {
                cout << "Error in GetDCValueSpecies_MT ... return zero value"
                     << "\n";
                DC_MT_cur = 0.0;
            }

            return DC_MT_cur;
            break;
        }
        i += 1;
    }
    // something went wrong...we exit the program
    cout << "error in GetDCValueSpecies_MT "
         << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
    MPI_Finalize();  // make sure MPI exits
#endif
    exit(1);
}

short REACT_GEM::GetSoComponentValue_MT(long node_Index, int timelevel,
                                        double* m_Phase, TNode* m_Node)
{
    string str;
    int x_Component = 0;
    CRFProcess* m_pcs = NULL;
    for (size_t i = 0; i < pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
        //                == 0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            x_Component = -1;

            // get the name of compound from MT;
            str = m_pcs->pcs_primary_function_name[0];
            // get the index of certain compound, -1: no match
            x_Component = m_Node->Ph_name_to_xDB(str.c_str());
            if (x_Component > -1)
                *(m_Phase + x_Component) = m_pcs->GetNodeValue(
                    node_Index, m_pcs->GetNodeValueIndex(str) + timelevel);
            else
            {
                // DisplayErrorMsg("Error: Corresponding Component NOT FOUND in
                // GEM part!!"); return 0;
            }
        }
    }
    // DisplayErrorMsg("Error: MASS TRANSPORT NOT FOUND!!");
    return 1;
}
/*
   short REACT_GEM::SetDCValue_MT ( long node_Index, int timelevel, double* m_DC
   ) // This routine does not work properly!!!!!!!!!!!!!
   {
    CRFProcess* m_pcs = NULL;
    string str;
    for ( int i=0; i < nDC ; i++ )
    {

        m_pcs = pcs_vector[i+1]; // DANGEROUS: this does not work with more than
   1 process before

        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
   == 0 ) { if ( m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT )
        {
            str = m_pcs->pcs_primary_function_name[0];
            if ( str.compare ( "pH" ) != 0 && str.compare ( "pe" ) != 0 &&
   str.compare ( "Eh" ) != 0 && str.compare ( "NodePorosity" ) != 0 )
            {
                if ( flag_iterative_scheme > 0 )
                {
                    if ( CPGetMobil ( m_pcs->GetProcessComponentNumber() ) > 0 )
                    {
                        // m_pcs->eqs->b[node_Index] +=
   m_xDC_Chem_delta[node_Index*nDC+i] / dt ;

                        m_pcs->SetNodeValue ( node_Index ,
   m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_DC+i ) );

                    }
                    else
                    {

                        m_pcs->SetNodeValue ( node_Index ,
   m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_DC+i ) );

                    }
                }
                else
                {

                    m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex
   ( str ) +timelevel , * ( m_DC+i ) );

                }
            }
        }
    }

    return 1;
   }
 */
short REACT_GEM::SetBValue_MT(long node_Index, int timelevel, double* m_soluteB)
{
    CRFProcess* m_pcs = NULL;
    string str;
    long i = -1;
    for (size_t j = 0; j < pcs_vector.size(); j++)
    {
        m_pcs = pcs_vector[j];  // not good!

        //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" )
        //                == 0 ) {
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            str = m_pcs->pcs_primary_function_name[0];
            i += 1;
            if (flag_iterative_scheme > 0)
            {
                if (CPGetMobil(m_pcs->GetProcessComponentNumber()) > 0)
                    // m_pcs->eqs->b[node_Index] +=
                    // m_xDC_Chem_delta[node_Index*nDC+i] / dt ;

                    m_pcs->SetNodeValue(
                        node_Index, m_pcs->GetNodeValueIndex(str) + timelevel,
                        *(m_soluteB + i));

                else

                    m_pcs->SetNodeValue(
                        node_Index, m_pcs->GetNodeValueIndex(str) + timelevel,
                        *(m_soluteB + i));
            }
            else

                m_pcs->SetNodeValue(node_Index,
                                    m_pcs->GetNodeValueIndex(str) + timelevel,
                                    *(m_soluteB + i));
        }
    }

    return 1;
}

// i_timestep 0: old timestep 1: new timestep
int REACT_GEM::SetPorosityValue_MT(long ele_Index, double m_porosity_Elem,
                                   int i_timestep)
{
    int idx;
    idx = -1;

    CRFProcess* m_pcs = NULL;

    if (flowflag > 0)
    {
        switch (flowflag)
        {
            case 1:
                m_pcs = PCSGet("GROUNDWATER_FLOW");
                idx = m_pcs->GetElementValueIndex("POROSITY");
                // set new porosity;
                m_pcs->SetElementValue(ele_Index, idx + i_timestep,
                                       m_porosity_Elem);
                break;
            case 2:
                m_pcs = PCSGet("LIQUID_FLOW");
                idx = m_pcs->GetElementValueIndex("POROSITY");
                // always write into the new step
                m_pcs->SetElementValue(ele_Index, idx + i_timestep,
                                       m_porosity_Elem);
                break;
            case 3:
                m_pcs = PCSGet("RICHARDS_FLOW");
                idx = m_pcs->GetElementValueIndex("POROSITY");
                // always write into the new step
                m_pcs->SetElementValue(ele_Index, idx + i_timestep,
                                       m_porosity_Elem);
                break;
            case 4:  // kg44: do we have to update POROSITY_IL and POROSITY_SW?
                m_pcs = PCSGet("MULTI_PHASE_FLOW");
                idx = m_pcs->GetElementValueIndex("POROSITY1");
                // always write into the new step
                m_pcs->SetElementValue(ele_Index, idx + i_timestep,
                                       m_porosity_Elem);
                break;
            default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                if (myrank == 0 /*should be set to root*/)
#endif
                    cout << "Error: Not implemented for the flow in GEM case!!!"
                         << "\n";
                break;
        }
    }
    return 1;
}

int REACT_GEM::SetSourceSink_MT(long in, double time_step_size /*in sec*/)
{
    Water_ST_GEMS m_st;

    switch (flowflag)
    {
        case 1:  // groundwater flow
            m_st.index_node = in;
            m_st.water_st_value = m_excess_water[in] / time_step_size;
            // normalize with node volume
            m_st.water_st_value *= m_Node_Volume[in];
            m_flow_pcs->Water_ST_vec.push_back(m_st);
            return 1;
            break;
        case 2:  // liquid flow
            m_st.index_node = in;
            m_st.water_st_value = m_excess_water[in] / time_step_size;
            // normalize with node volume
            m_st.water_st_value *= m_Node_Volume[in];
            m_flow_pcs->Water_ST_vec.push_back(m_st);
            return 1;
            break;
        case 3:  // Richards flow
            m_st.index_node = in;
            m_st.water_st_value = m_excess_water[in] / time_step_size;
            // normalize with node volume
            m_st.water_st_value *= m_Node_Volume[in];
            m_flow_pcs->Water_ST_vec.push_back(m_st);
            return 1;
            break;
        case 4:  // multiphase flow...works with case 1 ...pressure saturation
                 // scheme
            m_st.index_node = in;
            m_st.water_st_value = m_excess_water[in] / time_step_size;
            // normalize with node volume
            m_st.water_st_value *= m_Node_Volume[in];
            m_flow_pcs->Water_ST_vec.push_back(m_st);
            return 1;
            break;
        default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            if (myrank == 0 /*should be set to root*/)
#endif
                cout << "Error: Not implemented for the flow in GEM case!!!"
                     << "\n";
            break;
    }
    return 0;
}

int REACT_GEM::FindWater_xDC(TNode* m_Node)
{
    // initialization
    int rt = -1;
    int i;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // loop over all the xDC names, and find the one that is water
    for (i = 0; i < nDC; i++)
        if (dCH->ccDC[i] == 'W')
        {
            rt = i;
            return rt;
        }
    return rt;
}

int REACT_GEM::Findhydrogen_bIC(TNode* m_Node)
{
    // initialization
    int rt = -1;
    int i;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // loop over all the bIC names, and find the one that is hydrogen
    for (i = 0; i < nIC; i++)
        if (dCH->ccIC[i] == 'h')
        {
            rt = i;
            return rt;
        }
    return rt;
}

int REACT_GEM::Findoxygen_bIC(TNode* m_Node)
{
    // initialization
    int rt = -1;
    int i;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // loop over all the bIC names, and find the one that is oxygen
    for (i = 0; i < nIC; i++)
        if (dCH->ccIC[i] == 'o')
        {
            rt = i;
            return rt;
        }
    return rt;
}

double REACT_GEM::Pressure_Pa_2_Bar(double Pre_in_Pa)
{
    return Pre_in_Pa / 1.0e+5;
}

double REACT_GEM::Pressure_Bar_2_Pa(double Pre_in_Bar)
{
    return Pre_in_Bar * 1.0e+5;
}

double REACT_GEM::Pressure_M_2_Bar(double Pre_in_M, double flu_density)
{
    return Pre_in_M * 9.81 * flu_density / 1.0e+5;
}

double REACT_GEM::Pressure_Bar_2_M(double Pre_in_Bar, double flu_density)
{
    return Pre_in_Bar * 1.0e5 / 9.81 / flu_density;
}

double REACT_GEM::Pressure_M_2_Pa(double Pre_in_M, double flu_density)
{
    return Pre_in_M * 9.81 * flu_density;
}

double REACT_GEM::Pressure_Pa_2_M(double Pre_in_Pa, double flu_density)
{
    return Pre_in_Pa / 9.81 / flu_density;
}

double REACT_GEM::GetNodeAdjacentVolume(long Idx_Node)
{
    double volume;
    long Idx_Ele;
    int number_of_nodes;
    volume = 0.0;
    number_of_nodes = 0;

    MeshLib::CNode* m_ogsNode;
    MeshLib::CElem* m_Elem;

    // get the pointer to current node;
    m_ogsNode = m_flow_pcs->m_msh->nod_vector[Idx_Node];

    // loop over all the elements that adjacent to this node;
    for (int i = 0; i < (long)m_ogsNode->getConnectedElementIDs().size(); i++)
    {
        // get the index of current element;
        Idx_Ele = m_ogsNode->getConnectedElementIDs()[i];

        // get the pointer of this element;
        m_Elem = m_flow_pcs->m_msh->ele_vector[Idx_Ele];

        // get the number of nodes in this element;
        // given argument "false" means giving node number instead of Gauss
        // points;
        number_of_nodes = m_Elem->GetNodesNumber(false);

        // taking part of volume from this element;
        volume += m_Elem->GetVolume() / number_of_nodes;
    }

    return volume;
}

// i_timestep: 0: old timestep 1: new timestep
void REACT_GEM::ConvPorosityNodeValue2Elem(int i_timestep)
{
    long i, idx_Node, group;
    int j, number_of_nodes;
    double pormin = 2.0, pormax = 0.0;
    double distance, weight, sum_weights;
    MeshLib::CElem* m_Elem;
    MeshLib::CNode* m_ogsNode;
    CMediumProperties* m_mmp = NULL;

    for (i = 0; i < nElems; i++)
    {
        distance = weight = sum_weights = 0.0;
        m_Elem = m_flow_pcs->m_msh->ele_vector[i];
        group = m_Elem->GetPatchIndex();
        m_mmp = mmp_vector[group];
        if (i_timestep ==
            0)  // only do something if we have model 15...gems changes porosity
        {
            // push back porosities
            if (!(m_mmp->porosity_model == 15))
                m_porosity_Elem[i] =
                    m_mmp->Porosity(i, 1);  // get element porosity
            SetPorosityValue_MT(i, m_porosity_Elem[i], i_timestep);
        }
        else if (m_mmp->porosity_model ==
                 15)  // only do something if porosity model 15
        {
            // first set the parameters to zero;
            m_porosity_Elem[i] = 0.0;
            number_of_nodes = (int)m_Elem->GetNodesNumber(false);
            // then get the values from nodes
            for (j = 0; j < number_of_nodes; j++)
            {
                // get the connected nodes;
                idx_Node = m_Elem->GetNodeIndex(j);
                m_ogsNode = m_flow_pcs->m_msh->nod_vector[idx_Node];

                // calculate distance between the node and the barycentre
                double const* gravity_centre(m_Elem->GetGravityCenter());
                double const* const pnt(m_ogsNode->getData());
                distance =
                    (gravity_centre[0] - pnt[0]) * (gravity_centre[0] - pnt[0]);
                distance +=
                    (gravity_centre[1] - pnt[1]) * (gravity_centre[1] - pnt[1]);
                distance +=
                    (gravity_centre[2] - pnt[2]) * (gravity_centre[2] - pnt[2]);
                distance = sqrt(distance);

                // Weight of each face depending on distance
                weight = (1.0 / distance);
                // Sum of weights
                sum_weights += weight;

                // this is arithmetric mean
                m_porosity_Elem[i] += m_porosity[idx_Node] / number_of_nodes;
                // here we use harmonic mean, as porosity is used for
                // permeability/diffusivity changes....flux in the element is
                // strongly influenced by the minimum values cout << " porosity
                // " << idx_Node <<" "<<m_porosity[idx_Node] << "\n";
                // m_porosity_Elem[i] += 1.0/m_porosity[idx_Node] ; // this is
                // for harmonic mean
            }
            //	cout << " Porosity: Element "<< i << " " <<m_porosity_Elem[i] <<
            //"\n";

            // m_porosity_Elem[i] = (double) number_of_nodes /
            // m_porosity_Elem[i];

            // upper limit of porosity
            if (m_porosity_Elem[i] >= max_possible_porosity)
                m_porosity_Elem[i] = max_possible_porosity;
            // lower limit of porosity..
            if (m_porosity_Elem[i] <= min_possible_porosity)
                m_porosity_Elem[i] = min_possible_porosity;

            pormin = min(pormin, m_porosity_Elem[i]);
            pormax = max(pormax, m_porosity_Elem[i]);

            // push back porosities
            SetPorosityValue_MT(i, m_porosity_Elem[i], i_timestep);
        }
        else
            m_porosity_Elem[i] = m_mmp->Porosity(i, 1);  // get element porosity
                                                         // debug
        // cout << " GEMS3K DEBUG elment material group porosity " << i << " "
        // << m_mmp->porosity_model << " " << m_porosity_Elem[i] << "\n";
    }
    if (i_timestep == 1)
        cout << "min, max porosity: " << pormin << " " << pormax
             << "\n";  // only output for current time step....old time step
                       // will give wrong values
}

double REACT_GEM::FluidDensity(long elem, int gaussnode)
{
    long idx_Node;
    int number_of_nodes, j, i;
    MeshLib::CElem* m_Elem;
    MeshLib::CNode* m_ogsNode;
    double density;
    double distance, weight, sum_weights;
    int size_m = 20;  // assigned to the value in CFiniteElementStd(CRFProcess
                      // *Pcs, const int C_Sys_Flad, const int order=1);
    double* NodalVal_BG;

    density = 1000.0;  // default backup value

    distance = weight = sum_weights = 0.0;

    //  cout << "GEMS DEBUG: FluidDensity " << density << "elem " << elem << "
    //  gaussnode " << gaussnode<< "\n";

    m_Elem = m_flow_pcs->m_msh->ele_vector[elem];
    number_of_nodes = (int)m_Elem->GetNodesNumber(false);

    if (gaussnode > -1)  // gauss point like in InterpolatPorpertytoGausPoints
    {
        if (!Fem_Ele_Std)
        {
            cout << "DEBUG REACTGEM: fluiddensity from gauss node failed: "
                    "could not get Fem_Ele_Std"
                 << "\n";
            return density;  // make sure call to interpolate does not fail
        }
        NodalVal_BG = new double[size_m];  // BG
        // Get gauss point data
        // GetGaussData(gp, gp_r, gp_s, gp_t);
        // fkt = GetGaussData(GPIndex, gp_r, gp_s, gp_t);
        // Compute the shape function for interpolation within element
        // ComputeShapefct(1);
        // read density from nodes
        for (i = 0; i < number_of_nodes; i++)
        {
            idx_Node = m_Elem->GetNodeIndex(i);
            NodalVal_BG[i] = m_fluid_density[idx_Node];
        }
        // Interpolate density from nodes to gauss point
        density = Fem_Ele_Std->interpolate(NodalVal_BG);
    }
    else  // arithmetric average on element
    {
        // get the values from nodes
        for (j = 0; j < number_of_nodes; j++)
        {
            // get the connected nodes;
            idx_Node = m_Elem->GetNodeIndex(j);
            m_ogsNode = m_flow_pcs->m_msh->nod_vector[idx_Node];
            // calculate distance between the node and the barycentre
            double const* gravity_centre(m_Elem->GetGravityCenter());
            double const* const pnt(m_ogsNode->getData());
            distance =
                (gravity_centre[0] - pnt[0]) * (gravity_centre[0] - pnt[0]);
            distance +=
                (gravity_centre[1] - pnt[1]) * (gravity_centre[1] - pnt[1]);
            distance +=
                (gravity_centre[2] - pnt[2]) * (gravity_centre[2] - pnt[2]);
            distance = sqrt(distance);

            // Weight of each face depending on distance
            weight = (1.0 / distance);
            // Sum of weights
            sum_weights += weight;
            // Density
            // this is arithmetric mean
            density += m_fluid_density[idx_Node] * weight;
        }
        density = density / sum_weights;
    }

    //  cout << "GEMS DEBUG: FluidDensity " << density << "\n";
    return density;
}

int REACT_GEM::CalcPorosity(long in, TNode* m_Node)
{
    int k;        // must be int for Ph_Volume !
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    m_porosity[in] = 0.0;
    for (k = 0; k < dCH->nPHb; k++)
        if ((dCH->ccPH[k] == 's') || (dCH->ccPH[k] == 'x') ||
            (dCH->ccPH[k] == 'd'))
            m_porosity[in] += m_Node->Ph_Volume(k);

    // normalized by m_Vs
    //	m_porosity[in] = 1.0 - m_porosity[in] / ( m_Vs[in] * 1.0e-6 /*convert to
    //cm3 here*/) ;
    // kg44 this is the correct way to do it
    m_porosity[in] = 1.0 - (m_porosity[in]);

    //	cout <<" porosity:" << m_porosity[in] << " node: "<< in <<"\n";
    // checking whether out of bounary. only needed if calc porosity is done

    if (m_porosity[in] >= max_possible_porosity)
        m_porosity[in] = max_possible_porosity;  // upper limit of porosity
    if (m_porosity[in] <= min_possible_porosity)
        m_porosity[in] =
            min_possible_porosity;  // lower limit of porosity..	     cout << "
                                    // skal factor " << skal_faktor << "
    // excess water volume " << m_excess_water[in] ;
    // cout <<" porosity:" << m_porosity[in] << " node: "<< in <<" Vs
    // "<<m_Vs[in]<<"\n";
    return 1;
}

int REACT_GEM::MassToConcentration(
    long in, int i_failed,
    TNode* m_Node)  // attention second argument is not timestep. I
// is a flag that indicates if we deal with
// failed nodes!!!!...do not get data from GEMS
// for this nodes
{
    // converting the value from moles to the value in mol/m^3 water.
    long i, j, k, ii;
    int idx;
    double gas_volume, fluid_volume;
    double skal_faktor = 0.0, skal_faktor_gas = 0.0, saturation = 1.0;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    // get the fluid volume
    if (i_failed)
    {
        fluid_volume = m_fluid_volume[in];
        gas_volume = m_gas_volume[in];
    }
    else
    {
        fluid_volume = 0.0;
        gas_volume = 0.0;
        for (k = 0; k < dCH->nPHb; k++)
        {
            if (dCH->ccPH[k] == 'a')
                fluid_volume += m_Node->Ph_Volume(k);
            if (dCH->ccPH[k] == 'g')
                gas_volume += m_Node->Ph_Volume(k);

            // * 1.0e-6;    // transform cm3 -> m3 !!
        }
        //     cout << " node " << in << " fluid_volume " << fluid_volume << "
        //     gas volume " << gas_volume << " total-gas: " <<
        //     m_Vs[in]-gas_volume << "\n";
    }

    if ((fluid_volume <= 0.0))
    {
        cout << "OGSGEM MassToConcentration: fluid volume negative or zero "
             << fluid_volume << " node " << in << " " << m_fluid_volume[in]
             << "\n";
        m_Node->GEM_write_dbr("dbr_for_crash_node_fluid_volume.txt");
        cout << " skal factor " << skal_faktor << " excess water volume "
             << m_excess_water[in];
        m_Node->GEM_print_ipm("ipm_for_crash_node_fluid_volume.txt");

#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif
        exit(1);
    }
    if ((gas_volume < 0.0))
    {
        cout << "OGSGEM MassToConcentration: gas volume negative" << gas_volume
             << " node " << in << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif

        exit(1);
    }

    // store volumes for reuse during concentration_to_mass
    m_fluid_volume[in] = fluid_volume;
    m_gas_volume[in] = gas_volume;

    // calculate excess water: volume of fluid phase bigger/smaller than
    // porosity; requires updated porosity (calc porosity before calculation of
    // execc water)
    CRFProcess* m_pcs = NULL;
    switch (flowflag)
    {
        case 1:  // groundwater flow	     cout << " skal factor " <<
                 // skal_faktor << " excess water volume " <<
            // m_excess_water[in] ;
            m_excess_water[in] = m_fluid_volume[in] - m_porosity[in];
            skal_faktor =
                (m_porosity[in]) /
                m_fluid_volume[in];  // scale water to available volume
            m_fluid_volume[in] = m_porosity[in];
            break;
        case 2:  // liquid flow
            m_excess_water[in] = m_fluid_volume[in] - m_porosity[in];
            skal_faktor =
                (m_porosity[in]) /
                m_fluid_volume[in];  // scale water to available volume
            m_fluid_volume[in] = m_porosity[in];
            break;
        case 3:  // Richards flow: Saturation is directly calculated from the
                 // pressure, therefore add a source-sink term
            // to induce pressure/saturation changes
            m_pcs = PCSGet("RICHARDS_FLOW");
            idx = m_pcs->GetNodeValueIndex("SATURATION1");

            if (m_fluid_volume[in] >
                m_porosity[in])  // fluid volume exceeds available pore space
            {
                m_excess_water[in] = m_fluid_volume[in] - m_porosity[in];
                skal_faktor =
                    m_porosity[in] /
                    m_fluid_volume[in];  // mol amount of water in first phase
                m_fluid_volume[in] = m_porosity[in];
                m_pcs->SetNodeValue(
                    in, idx, 1.0);  // instead change saturation accordingly;
                //    cout << " skal factor " << skal_faktor << " excess water
                //    volume " << m_excess_water[in] ;
            }
            else
            {
                m_excess_water[in] =
                    m_fluid_volume[in] -
                    m_porosity[in] * m_pcs->GetNodeValue(in, idx + 1);
                skal_faktor = m_porosity[in] *
                              m_pcs->GetNodeValue(in, idx + 1) /
                              m_fluid_volume[in];  // KG44 18.10.2012  is this
                                                   // really correct??
                m_fluid_volume[in] =
                    m_porosity[in] * m_pcs->GetNodeValue(in, idx + 1);

                //	     m_pcs->SetNodeValue ( in, idx, m_fluid_volume[in] /
                //m_porosity[in]); // instead change
                // saturation accordingly; this is done always.
            }
            //  m_excess_gas[in] = m_gas_volume[in]- m_porosity[in]* ( 1.0 -
            //  m_pcs->GetNodeValue ( in,idx+1 ) ); // in contrast to fluid this
            //  is not really strongly coupled ...we calculate the stuff but do
            //  not couple it to gas flow (yet) skal_faktor_gas= (
            //  m_porosity[in]* ( 1.0-m_pcs->GetNodeValue ( in,idx+1 ) ) )
            //  /m_gas_volume[in]; m_gas_volume[in] = m_porosity[in]* (
            //  1-m_pcs->GetNodeValue ( in,idx+1 ) );

            break;
        case 4:  // multiphase flow...works with case 1 ...pressure saturation
                 // scheme
            m_pcs = PCSGet("MULTI_PHASE_FLOW");
            idx = m_pcs->GetNodeValueIndex("SATURATION0");
            m_excess_water[in] =
                m_fluid_volume[in] -
                m_porosity[in] * m_pcs->GetNodeValue(in, idx + 1);
            m_excess_gas[in] =
                m_gas_volume[in] -
                m_porosity[in] * (1.0 - m_pcs->GetNodeValue(in, idx + 1));
            if (m_fluid_volume[in] > 0.0)
                m_xDC[in * nDC + idx_water] *=
                    m_pcs->GetNodeValue(in, idx + 1) * m_porosity[in] /
                    m_fluid_volume[in];
            else
                m_xDC[in * nDC + idx_water] = 0.0;

            break;
        default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            if (myrank == 0 /*should be set to root*/)
#endif
                cout << "Error: Not implemented for the flow in GEM case!!!"
                     << "\n";
            break;
    }
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
//	if ( fabs ( m_excess_water_buff[in] ) >= 0.01 ) cout << "node "<< in <<"
//m_excess_water" << m_excess_water_buff[in]
//<<"\n";
//	if ( fabs ( m_excess_gas_buff[in] ) >= 0.01 ) cout << "node "<< in <<"
//m_excess_gas" << m_excess_water_buff[in]
//<<"\n";
#else
//	if ( fabs ( m_excess_water[in] ) >= 0.01 ) cout << "node "<< in <<"
//m_excess_water " << m_excess_water[in] <<"\n"; 	if ( fabs ( m_excess_gas[in] )
//>= 0.01 ) cout << "node "<< in <<" m_excess_gas " << m_excess_water[in]
//<<"\n";
#endif
    // if positive, then source, otherwise sink.
    // change amount of fluid .... at the moment only water!!!!!!

    // if ( fabs ( m_excess_water[in] ) >= 1.e-10 ) cout << "node "<< in
    // <<"m_excess_water" << m_excess_water[in]
    // <<"\n";
    // do we need to scale this accoring to the node volume?..I think so: it is
    // done while transering this to GEMS
    if (flowflag < 3)
    {
        for (j = 0; j < nIC; j++)
        {
            i = in * nIC + j;
            ii = in * nPS * nIC + 0 * nIC +
                 j;  // corresponding index of first phase (fluid) for m_bPS
            m_bIC[i] -= m_bPS[ii];  // B vector without solute
            if (flag_coupling_hydrology)
                m_bPS[ii] *=
                    skal_faktor;  // completely newly scaled first phase ...only
                                  // for coupling with hydraulics
            else
            {  // this is necessary for not coupling with hydrology, but
               // porosity change or water is used by gems....
                if (idx_hydrogen == j)
                    m_bPS[ii] *= skal_faktor;
                if (idx_oxygen == j)
                    m_bPS[ii] *= skal_faktor;
            }
        }
        // for ( j=0 ; j <= idx_water; j++ )
        j = idx_water;
        {
            i = in * nDC + j;
            m_xDC[i] *= skal_faktor;  // newly scaled xDC including water
                                      // /excluding rest
        }

        // here we do not need to account for different flow processes ....
        //****************************************************************************************
        if (flag_transport_b == 1)
            for (j = 0; j < nIC; j++)
            {
                i = in * nIC + j;
                ii = in * nPS * nIC + 0 * nIC +
                     j;  // corresponding index of first phase (fluid) for m_bPS

                // correct for water and transform bPS into concentrations
                // carrier for zero(first phase)  is normally water!
                if (idx_hydrogen == j)
                    m_soluteB[i] =
                        m_bPS[ii] - (2.0 * m_xDC[in * nDC + idx_water]);
                else if (idx_oxygen == j)
                    m_soluteB[i] = m_bPS[ii] - m_xDC[in * nDC + idx_water];
                else
                    m_soluteB[i] = m_bPS[ii];

                m_soluteB[i] /=
                    m_fluid_volume[in];  // now these are the concentrations
            }

    }                        // flowflag < 3 is finished
    else if (flowflag == 3)  // richards flow ....no corrections in gas phase
                             // for Richards flow!
    {
        for (j = 0; j < nIC; j++)
        {
            i = in * nIC + j;
            ii = in * nPS * nIC + 0 * nIC +
                 j;  // corresponding index of first phase (fluid) for m_bPS
            m_bIC[i] -= (m_bPS[ii]);      // B vector without solute
            if (flag_coupling_hydrology)  // in any case this accounts only for
                                          // fully saturated conditions

                m_bPS[ii] *=
                    skal_faktor;  // completely newly scaled first phase ...only
                                  // for coupling with hydraulics
            else
            {  // this is necessary for not coupling with hydrology, but
               // porosity change or water is used by gems....
                if (idx_hydrogen == j)
                    m_bPS[ii] *= skal_faktor;
                if (idx_oxygen == j)
                    m_bPS[ii] *= skal_faktor;
            }
        }
        // for ( j=0 ; j <= idx_water; j++ )
        j = idx_water;
        {
            i = in * nDC + j;
            m_xDC[i] *= skal_faktor;  // newly scaled xDC for water
        }

        // here we do not need to account for different flow processes .... as
        // long as we work only with one wetting phase (water!!!) ...for liqid
        // CO2 we have to change something ....
        //****************************************************************************************
        if (flag_transport_b == 1)
            for (j = 0; j < nIC; j++)
            {
                i = in * nIC + j;
                ii = in * nPS * nIC + 0 * nIC +
                     j;  // corresponding index of first phase (fluid) for m_bPS
                // correct for water and transform bPS into concentrations
                // carrier for zero(first phase)  is normally water!
                if (idx_hydrogen == j)
                    m_soluteB[i] =
                        m_bPS[ii] - (2.0 * m_xDC[in * nDC + idx_water]);
                else if (idx_oxygen == j)
                    m_soluteB[i] = m_bPS[ii] - m_xDC[in * nDC + idx_water];
                else
                    m_soluteB[i] = m_bPS[ii];

                m_soluteB[i] /=
                    m_fluid_volume[in];  // now these are the concentrations
            }
    }
    return 1;
}

int REACT_GEM::ConcentrationToMass(long l /*idx of node*/, int i_timestep)
{
    // converting the value from mol/m^3 water to moles.
    long i, j, ii;
    double water_volume = 0.0, scale_factor = 1.0;
    int idx;
    string ErrorOut;

    //	cout <<"water volume " << water_volume << "\n";
    // I think here it has to be different than in MassToConcentration module,
    // as after hydraulics, the volume is limited to porosity As we work with an
    // unit volume of "one", it should be save to directly take porosity as
    // volume.... of course we have to differenciate between the different flow
    // models
    CRFProcess* m_pcs = NULL;
    switch (flowflag)
    {
        case 1:                                // groundwater flow
            water_volume = m_fluid_volume[l];  // kg44: reuse stored fluid
                                               // volume from last timestep!
            break;
        case 2:  // liquid flow
            water_volume = m_fluid_volume[l];
            break;
        case 3:  // Richards flow
            m_pcs = PCSGet("RICHARDS_FLOW");
            idx = m_pcs->GetNodeValueIndex("SATURATION1");
            scale_factor =
                m_porosity[l] *
                m_pcs->GetNodeValue(
                    l, idx + i_timestep);      // current volume of water phase
                                               // after hydraulic step
            water_volume = m_fluid_volume[l];  // old volume of water phase
                                               // after last GEMS step
            scale_factor /= water_volume;      // 1: same volume <1 decrease in
                                           // water volume >1 incresing volume
            m_fluid_volume[l] = water_volume;  // update fluid volume
            // cout << "conctomass " << m_fluid_volume[l] << " " << scale_factor
            // <<" " ;
            break;
        case 4:  // multiphase flow...works with case 1 ...pressure saturation
                 // scheme
            m_pcs = PCSGet("MULTI_PHASE_FLOW");
            idx = m_pcs->GetNodeValueIndex("SATURATION0");
            water_volume =
                m_porosity[l] * m_pcs->GetNodeValue(l, idx + i_timestep);
            break;
        default:
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
            if (myrank == 0 /*should be set to root*/)
#endif
                cout << "Error: Not implemented for the flow in GEM case!!!"
                     << "\n";
            break;
    }

    if ((water_volume < min_possible_porosity) || (water_volume > 1.0))
    {
        cout << "conctomass water volume " << water_volume << " at node: " << l
             << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif

        exit(1);
    }

    if (flowflag < 3)  // groundwater flow or fluid flow

        for (j = 0; j < nIC; j++)
        {
            i = l * nIC + j;

            // also check if xDC values is negative
            if (m_soluteB[i] < 0.0)
                // make sure it is greater/eq zero and take value from last
                // timestep ---> second argument is zero!
                m_soluteB[i] = fabs(GetDCValueSpecies_MT(l, 0, (int)j));
            m_soluteB[i] *= water_volume;
            // carrier for zero(first phase)  is normally water!
            if (idx_hydrogen == j)
                m_soluteB[i] += (2.0 * m_xDC[l * nDC + idx_water]);
            else if (idx_oxygen == j)
                m_soluteB[i] += m_xDC[l * nDC + idx_water];

            m_bIC[i] += m_soluteB[i];  // updated B vector for GEMS
            // here we check again if the vector if negative or smaller than a
            // minimum amount...if yes we add some stuff.... adding 10-6 Mol/m^3
            // should be save, as this corresponds aprox to 10-9 Mol/kg ..which
            // is well above the accuracy for most calculations
            if (m_bIC[i] <= 1.0e-8)
                m_bIC[i] = 1e-6;
            //                        cout <<  " i " << i << " " << m_bIC[i] <<
            //                        "\n";
        }
    else if (flowflag == 3)  // Richards flow
    {
        for (j = 0; j < nIC; j++)
        {
            i = l * nIC + j;
            ii = l * nPS * nIC + 0 * nIC +
                 j;  // corresponding index of first phase (fluid) for m_bPS
            // cout <<  " i " << i << " " << m_soluteB[i] << " " <<m_bPS[ii] <<
            // " scale factor " << scale_factor <<
            // "\n";

            // test for NaN!! ---seems necessary as nan run through the system
            //  this part is for debug reasons....
            if (!(m_soluteB[i] <= 1.0) && !(m_soluteB[i] > 1.0))
            {
                cout << " m_soluteB[i] " << m_soluteB[i] << " is NaN "
                     << " i " << i << " stop calculations "
                     << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                MPI_Finalize();  // make sure MPI exits
#endif
                exit(1);
            }
            // also check if xDC values is negative
            if (m_soluteB[i] < 0.0)
                // make sure it is greater/eq zero and take value from last
                // timestep ---> second argument is zero!
                m_soluteB[i] = fabs(GetDCValueSpecies_MT(l, 0, (int)j));
            m_soluteB[i] *=
                water_volume *
                scale_factor;  // match solutes with old fluid volume and
                               // account for hydraulic changes
            // carrier for zero(first phase)  is normally water!
            if (idx_hydrogen == j)
                m_soluteB[i] +=
                    (2.0 * m_xDC[l * nDC + idx_water]) * scale_factor;
            // match to volume after hydraulic step
            else if (idx_oxygen == j)
                m_soluteB[i] += m_xDC[l * nDC + idx_water] * scale_factor;

            //            m_bIC[i]+= ( m_soluteB[i]+m_bPS[ii] ); //updated B
            //            vector for GEMS ....no correction in gas phase for
            //            Richards flow!
            m_bIC[i] +=
                (m_soluteB[i]);  // updated B vector for GEMS ....no correction
                                 // in gas phase for Richards flow!

            // here we check again if the vector if negative or smaller than a
            // minimum amount...if yes we add some stuff.... adding 10-6 Mol/m^3
            // should be save, as this corresponds aprox to 10-9 Mol/kg ..which
            // is well above the accuracy for most calculations
            if (m_bIC[i] <= 1.0e-8)
                m_bIC[i] = 1e-6;
        }
    }
    return 1;
}

void REACT_GEM::CopyCurXDCPre(void)
{
    long i;
    for (i = 0; i < nNodes * nDC; i++)
        m_xDC_pts[i] = m_xDC[i];
}

void REACT_GEM::UpdateXDCChemDelta(void)
{
    long i;
    for (i = 0; i < nNodes * nDC; i++)
        m_xDC_Chem_delta[i] = m_xDC[i] - m_xDC_pts[i];
}

void REACT_GEM::CopyCurBPre(void)
{
    long i;
    for (i = 0; i < nNodes * nIC; i++)
    {
        m_soluteB_pts[i] = m_soluteB[i];
        m_bIC_pts[i] = m_bIC[i];
    }
}

double REACT_GEM::CalcSoluteBDelta(long in)
{
    long i;
    double dummy = 0;
    for (i = 0; i < nIC - 1; i++)
        dummy = max(dummy,
                    abs(m_soluteB[in * nIC + i] - m_soluteB_pts[in * nIC + i]) /
                        m_soluteB[in * nIC + i]);
    return dummy;
}

void REACT_GEM::RestoreOldSolution(long in)
{
    long i;
    for (i = 0; i < nIC - 1; i++)
    {
        m_soluteB[in * nIC + i] = m_soluteB_pts[in * nIC + i];
        m_bIC[in * nIC + i] = m_bIC_pts[in * nIC + i];
    }
    for (i = 0; i < nDC; i++)
        m_xDC[in * nDC + i] = m_xDC_pts[in * nDC + i];
}

bool GEMRead(string base_file_name, REACT_GEM* m_GEM_p)
{
    cout << "GEMRead"
         << "\n";
    char line[MAX_ZEILE];
    string sub_line;
    string line_string;
    ios::pos_type position;
    //========================================================================
    // file handling
    string gem_file_name;
    gem_file_name = base_file_name + GEM_FILE_EXTENSION;
    ifstream gem_file(gem_file_name.data(), ios::in);
    if (!gem_file.good())
    {
        cout << "! Error in GEMRead: No GEM data !"
             << "\n";
        return false;
    }
    gem_file.seekg(0L, ios::beg);
    //========================================================================
    // keyword loop
    while (!gem_file.eof())
    {
        gem_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != string::npos)
            return true;
        //----------------------------------------------------------------------
        // keyword found
        if (line_string.find("#GEM_PROPERTIES") != string::npos)
        {
            position = m_GEM_p->Read(&gem_file);
            gem_file.seekg(position, ios::beg);
        }  // keyword found
    }      // eof
    return true;
}

/** Read: This function reads the input file with ending ".gem"
 *
 *
 */
ios::pos_type REACT_GEM::Read(std::ifstream* gem_file)
{
    Kinetic_GEMS d_kin;  // dummy kinetic vector
    int j;
    // Initialization----------------
    string sub_line;
    string line_string;
    string delimiter(" ");
    bool new_keyword = false;
    string hash("#");
    ios::pos_type position;
    string sub_string;
    string dollar("$");
    string delimiter_type(":");
    std::stringstream in;
    // ------------------------------

    // Loop over all the key words----------------------
    while (!new_keyword)
    {
        position = gem_file->tellg();
        line_string = GetLineFromFile1(gem_file);
        if (line_string.size() < 1)
            break;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            break;
        }
        /// Key word "$GEM_INIT_FILE" : next line contains the name of the GEMS
        /// file with ending ".lst" that defines the setup, thermodynamic data
        /// and the numerical settings for the
        /// simulation........................
        if (line_string.find("$GEM_INIT_FILE") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> init_input_file_path;
            in.clear();
            continue;
        }

        /// Key word "$GEM_THREADS" : next line contains the number of GEMS3K
        /// worker threads to be started (for each MPI process). The minimum
        /// number of threads is 1!.......................
        if (line_string.find("$GEM_THREADS") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> gem_nThread;
            cout << "OGS-GEM gem read: number of threads: " << gem_nThread
                 << "\n";
            in.clear();
            continue;
        }

        /// Key word "$TRANSPORT_B": next line should be "1" to indicate the
        /// transport of total concentrations of independet species in fluid
        /// phase. The transport of all dependent species in the fluid phase is
        /// currently not supported anymore. .......................
        if (line_string.find("$TRANSPORT_B") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_transport_b;
            if (!flag_transport_b)
            {
                cout << "Transport of all species is currently not supported! "
                        "stop program"
                     << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                MPI_Finalize();  // make sure MPI exits
#endif

                exit(1);
            }
            in.clear();
            continue;
        }

        // ......................................................
        /// Key word "$FLAG_POROSITY_CHANGE": next line could be either "0" or
        /// "1" (recommended). "1": changes in the volume of solid phases cause
        /// the porosity to change   .......................
        if (line_string.find("$FLAG_POROSITY_CHANGE") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_porosity_change;
            in.clear();
            continue;
        }
        // ......................................................
        /// Key word "$MIN_POROSITY": next line should contain the minim
        /// porosity which has to be bigger than "0" e.g. "1e-6"
        /// .............................
        if (line_string.find("$MIN_POROSITY") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> min_possible_porosity;
            in.clear();
            continue;
        }
        // ......................................................
        /// Key word "$MAX_POROSITY": next line should contain the maximum
        /// porosity, recommended value: 1.0
        /// .............................
        if (line_string.find("$MAX_POROSITY") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> max_possible_porosity;
            in.clear();
            continue;
        }
        // ......................................................
        /// Key word "$FLAG_COUPLING_HYDROLOGY": If the system volume exceeds
        /// the reference volume (1 m^3) due to volume changes of solids and/or
        /// fluid phases the excess fluid is substracted from the volume. It is
        /// reinjected into the system during the next time step via source/sink
        /// time in the flow equation. ..................
        if (line_string.find("$FLAG_COUPLING_HYDROLOGY") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_coupling_hydrology;
            in.clear();
            continue;
        }
        // ......................................................
        /// Key word "$ITERATIVE_SCHEME": not yet implemented, set to zero or
        /// remove keyword from input file
        /// .........................
        if (line_string.find("$ITERATIVE_SCHEME") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_iterative_scheme;
            in.clear();
            continue;
        }
        // ......................................................
        /// Key word "$GEM_CALCULATE_BOUNDARY_NODES": By setting this to "1" it
        /// is possible to switch on GEMS calculations at boundary nodes (e.g.
        /// kinetically controlled porosity change at boundary nodes). This
        /// normally requires the definition of time dependent concentration
        /// values via time curves (see documentation of .bc files how to do
        /// this). ..................
        if (line_string.find("$CALCULATE_BOUNDARY_NODES") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_calculate_boundary_nodes;
            in.clear();
            continue;
        }
        // ......................................................
        // ......................................................
        /// Key word "$TEMPERATURE_GEM": with this keyword it is possible to
        /// change the overall temperature of the system (default 298.15
        /// K).........................
        if (line_string.find("$TEMPERATURE_GEM") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> m_gem_temperature;
            in.clear();
            continue;
        }
        /// Key word "$PRESSURE_GEM": with this keyword it is possible to
        /// disable the coupling between fluid pressures and GEMS calculations:
        /// set the desired value for GEMS calculations
        /// .........................
        if (line_string.find("$PRESSURE_GEM") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> m_gem_pressure;
            in.clear();
            gem_pressure_flag = 1;
            continue;
        }

        // ......................................................
        /// Key word "$MAX_FAILED_NODES" limits the number of failed nodes of
        /// GEMS calcuations for each MPI process. Run is exited if the number
        /// is exceeded..........................
        if (line_string.find("$MAX_FAILED_NODES") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> m_max_failed_nodes;
            in.clear();
            continue;
        }
        /// Key word "$MY_SMART_GEMS": if the summarized change in b-vector is
        /// smaller than this value, GEMS calculations are not done and previous
        /// values are taken.  .........................
        if (line_string.find("$MY_SMART_GEMS") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> m_diff_gems;
            in.clear();
            continue;
        }  // ......................................................
        /// Key word "$DISABLE_GEMS": gems calculations are not done / and no
        /// update of porosity etc...usefull for creating initial files if
        /// boundary conditions cause oszillations.........................
        if (line_string.find("$FLAG_DISABLE_GEM") != string::npos)
        {
            // subkeyword found
            in.str(GetLineFromFile1(gem_file));
            in >> flag_disable_gems;
            in.clear();
            continue;
        }  // ......................................................

        // kg44 26.11.2008 read in parameters for kinetics and GEM
        /** $KINETIC_GEM
         *       first line: general kinetic parameters
         * Name of phase (from GEMS setup) which is controlled
         * number of kinetic model
         * number of species activites used
         * double E_acid,E_neutral,E_base;  activation energies
         * double k_acid, k_neutral,k_base;  dissolution/precipitation rate
         * constants double q1,p1,q2,q3,p2,p3; exponents for omega name of first
         * species double n_1, n_2,n_3; exponents for acidic, neutral and base
         * cases for species one name of next species
         * .....three more exponents for each additional species
         * next line: reactive surface area
         * number of model for calculation of reactive surface area
         * model parameters for reactive surface area: e.g. specific reactive
         * surface area example: $KINETIC_GEM       ; Quartz phase no 45 (starts
         * at 0) wie Fernandez et al 2009 Quartz 1  1   0.0  0.0  0.0  0.0
         * -13.99 -16.29   0.0 1.0 0.0 1.0 1.0 1.0   H+  0.0  0.0 -0.5    ; this
         * are the kinetic parameters
         *   1 1.0e4    ; this are the parameters for reactive surface area
         */
        if (line_string.find("$KINETIC_GEM") != string::npos)
        {
            in.str(GetLineFromFile1(gem_file));
            in >> d_kin.phase_name >> d_kin.kinetic_model;
            if (d_kin.kinetic_model >= 1 && d_kin.kinetic_model <= 5)
            {
                cout << " found kinetics " << d_kin.kinetic_model << " "
                     << d_kin.phase_name << "\n";
                in >> d_kin.n_activities;
                if (d_kin.n_activities > 10)
                {
                    cout << "To many dependent species for GEM kinetic model "
                         << d_kin.n_activities << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                    MPI_Finalize();  // make sure MPI exits
#endif
                    exit(1);
                }
                //                cout <<" activities " << n_activities << "\n";
                // first general kinetic parameters
                //	0,1,2  double E_acid,E_neutral,E_base; // activation
                //energies
                in >> d_kin.kinetic_parameters[0] >>
                    d_kin.kinetic_parameters[1] >> d_kin.kinetic_parameters[2];
                //			cout << kinetic_parameters[0] << kinetic_parameters[1]
                //<< kinetic_parameters[1]<<"\n";

                //      3-5  double k_acid, k_neutral,k_base; //
                //      dissolution/precipitation rate constants
                in >> d_kin.kinetic_parameters[3] >>
                    d_kin.kinetic_parameters[4] >> d_kin.kinetic_parameters[5];
                //			cout << kinetic_parameters[3] << kinetic_parameters[4]
                //<< kinetic_parameters[5]<<"\n";

                //      6-11  double q1,p1,q2,q3,p2,p3; // exponents for omega
                in >> d_kin.kinetic_parameters[6] >>
                    d_kin.kinetic_parameters[7] >>
                    d_kin.kinetic_parameters[8] >>
                    d_kin.kinetic_parameters[9] >>
                    d_kin.kinetic_parameters[10] >>
                    d_kin.kinetic_parameters[11];
                for (j = 0; j < d_kin.n_activities; j++)
                {
                    in >> d_kin.active_species[j];
                    //				cout << active_species[j] ;
                    //      12,13,14  double n_1, n_2,n_3; // exponents for
                    //      acidic, neutral and base cases for species one

                    in >> d_kin.kinetic_parameters[j + 12] >>
                        d_kin.kinetic_parameters[j + 13] >>
                        d_kin.kinetic_parameters[j + 14];
                }
            }
            in.clear();
            // next line is surface area
            in.str(GetLineFromFile1(gem_file));
            in >> d_kin.surface_model;
            cout << d_kin.surface_model << "\n";
            if (d_kin.surface_model >= 1 ||
                d_kin.surface_model <
                    4)  // surface model 1, 2 and 3....only one parameter...
            {
                in >> d_kin.surface_area[0];  // surface: m*m / mol
                cout << "surface area " << d_kin.surface_area[0] << "\n";
            }
            else if (d_kin.kinetic_model ==
                     4)  // model 4 is kinetic a la crunch
            {
                in >> d_kin.surface_area[0];  // surface: m*m / mol
                cout << "mimic crunch: surface area " << d_kin.surface_area[0]
                     << "\n";
            }
            else if (d_kin.kinetic_model ==
                     5)  // model 5 is for solid solutions, needed parameters
                         // are the number of
            // endmembers and the surface area for each endmember
            {
                // next line are SS parameters
                in.str(GetLineFromFile1(gem_file));
                in >> d_kin.ss_endmembers;
                try
                {
                    d_kin.ss_scaling = new double[d_kin.ss_endmembers];
                }
                catch (bad_alloc)
                {
                    cout << "Reading Gems input: problem while allocating "
                            "memory for Solid Solution scaling parameters"
                         << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                    MPI_Finalize();  // make sure MPI exits
#endif
                    exit(1);
                }

                for (j = 0; j < d_kin.ss_endmembers; j++)
                    in >> d_kin.ss_scaling[j];
            }
            in.clear();

            // push back vector
            m_kin.push_back(d_kin);
        }  // subkeyword found
    }
    // End of looping over all the key words-----------
    return position;
}

/** calculate Reaction rates at a node in for all solids (for which kinetics is
 * defined)! kg44 30.07.2009 in: node temp: temperature
 */
int REACT_GEM::CalcReactionRate(long in, double temp, TNode* m_Node)
{
    int idx = 0, i, ii;
    long j, k;
    double rrn = 0.0, rrb = 0.0, rra = 0.0, sa = 0.0;
    double R = 8.31451070;  // molar gas konstant [J K-1 mol-1]
    double
        aa = 1.0,
        ab = 1.0,
        ac = 1.0;  // activity products ...species are input from material file
    double sactivity;  // dummy variable for extracting activities
    const char* species;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
        return 0;

    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
        return 0;
    /**
       int kinetic_model;  // only 1 = GEMS implemented right now
         int n_activities;  // number of species for activities
         string active_species[10];  // name for species ...maximum 10 names
       double kinetic_parameters[32];
       0,1,2  double E_acid,E_neutral,E_base; // activation energies
       3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate
       constants at standart konditions 6-11  double p1,q1,p2,q2,p3,q3; //
       exponents for omega 12,13,14  double n_1, n_2,n_3; // exponents for
       acidic neutral and base cases for species one append for each species
       another set of n_1, n_2  n_3 (up to 10 sets -> up to ten species)

     */

    // initialize in the first time step some variables...currently only
    // necessary for kinetic mocel 4 (implementation a la crunch)
    if (aktueller_zeitschritt < 1)
    {
        if (m_porosity_initial[in] < min_possible_porosity)
            m_porosity_initial[in] =
                m_porosity[in];  // make sure m_porosity_initial is not zero!
                                 // ...does not work
        // properly with RESTART!!!!!
        if (m_porosity_initial[in] < min_possible_porosity)
            m_porosity_initial[in] =
                m_porosity[in];  // secondary mineral ...see Cochepin et al 2008
        for (j = 0; j < nPH; j++)
        {
            m_volumes_initial[in * nPH + j] = m_Node->Ph_Volume(j);
            // secondary mineral ...see Cochepin et al 2008
            if (m_volumes_initial[in * nPH + j] <= 0.0)
                m_volumes_initial[in * nPH + j] = 0.01;
        }
    }

    // loop over all kinetic vectors and do for the defined phases and get rate
    // for each phase ... this algorithm assumes that we have correct ordering
    // of phases and components that belong into the phase

    for (ii = 0; ii < (int)m_kin.size(); ii++)
    {
        k = m_kin[ii].phase_number;

        if (m_kin[ii].kinetic_model >
            0)  // do it only if kinetic model is defined take model
        // kinetic_model==1 dissolution+precipitation kinetics
        // kinetic_model==2 only dissolution (no precipitation)aktueller
        // kinetic_mocel==3 only precipitation (no dissolution)
        {
            mol_phase[in * nPH + k] = 0.0;
            omega_phase[in * nPH + k] = 0.0;
            dmdt[in * nPH + k] = 0.0;

            if (m_kin[ii].kinetic_model ==
                5)  // special treatment of solid solution phases with e.g.
                    // vanselow convention
            {
                for (j = m_kin[ii].dc_counter;
                     j < m_kin[ii].dc_counter + dCH->nDCinPH[k]; j++)
                    // do not include surface complexation species!
                    if (!(dCH->ccDC[j] == '0') && !(dCH->ccDC[j] == 'X') &&
                        !(dCH->ccDC[j] == 'Y') && !(dCH->ccDC[j] == 'Z'))
                    {
                        // we need this later for solid solutions....
                        omega_components[in * nDC + j] = m_Node->DC_a(j);
                        // loop over all components of the phase
                        omega_phase[in * nPH + k] += m_Node->DC_a(j);

                        mol_phase[in * nPH + k] +=
                            (m_xDC[in * nDC + j] *
                             m_kin[ii].ss_scaling[j - m_kin[ii].dc_counter]);
                        // cout << "Debug kin omega phase "<<
                        // omega_phase[in*nPH+k] << " fraction component " <<
                        // omega_components[in*nDC+j] << "  mol phase " <<
                        // mol_phase[in*nPH+k] << "\n"; // debug
                    }
            }
            else  // normal behabviour for single component phases and SS which
                  // do have all the same endmember
            // characteristics
            {
                for (j = m_kin[ii].dc_counter;
                     j < m_kin[ii].dc_counter + dCH->nDCinPH[k]; j++)
                    // do not include surface complexation species!
                    if (!(dCH->ccDC[j] == '0') && !(dCH->ccDC[j] == 'X') &&
                        !(dCH->ccDC[j] == 'Y') && !(dCH->ccDC[j] == 'Z'))
                    {
                        //				omega_phase[k] += CalcSaturationIndex ( j,
                        //in,tempC,press ); // loop over all
                        // components of the phase
                        // we need this later for solid solutions....
                        omega_components[in * nDC + j] = m_Node->DC_a(j);
                        // loop over all components of the phase
                        omega_phase[in * nPH + k] += m_Node->DC_a(j);

                        mol_phase[in * nPH + k] += m_xDC[in * nDC + j];
                        // if (omega_phase[in*nPH+k] >1e6) cout << "Debug kin
                        // omega phase "<< omega_phase[in*nPH+k] << " fraction
                        // component " << omega_components[in*nDC+j] << "  mol
                        // phase " <<  mol_phase[in*nPH+k]
                        // << "\n"; // debug
                    }
            }

            //		        cout << omega_phase[k] << " " <<  mol_phase[k] <<
            //"\n"; // debug

            sa = REACT_GEM::SurfaceAreaPh(
                ii, in,
                m_Node);  // value for surface area in m^2...(specific surface
                          // area multiplied with volume of the phase)
            if (m_kin[ii].kinetic_model ==
                4)  // in the next part we try to mimic Crunchflow (at least
                    // partially)..could
            // be also done in CalcLimits if this makes the code easier to read
            {
                if (omega_phase[in * nPH + k] <
                    1.0)  // this is the dissolution case

                    sa *= pow(m_porosity[in] / m_porosity_initial[in] *
                                  m_Node->Ph_Volume(k) /
                                  m_volumes_initial[in * nPH + k],
                              0.66666666667);

                else if (omega_phase[in * nPH + k] >
                         1.0)  // this is the precipitation case

                    sa *= pow(m_porosity[in] / m_porosity_initial[in],
                              0.66666666667);
            }

            aa = 1.0;
            ab = 1.0;
            ac = 1.0;  // reset values for each phase!
            for (i = 0; i < m_kin[ii].n_activities; i++)
            {
                species = m_kin[ii].active_species[i].c_str();
                // loop over all the names in the list
                idx = m_Node->DC_name_to_xCH(species);
                if (idx < 0)
                {
                    if (m_porosity_initial[in] < min_possible_porosity)
                        m_porosity_initial[in] =
                            m_porosity[in];  // make sure m_porosity_initial is
                                             // not zero! ...does
                    // not work properly with RESTART!!!!!
                    if (m_porosity_initial[in] < min_possible_porosity)
                        m_porosity_initial[in] = min_possible_porosity;

                    cout << "failed CalcReactionRate: no DC-name "
                         << m_kin[ii].active_species[i] << " found"
                         << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                    MPI_Finalize();  // make sure MPI exits
#endif

                    exit(1);
                }
                else
                {
                    //      cout << "activities " <<aa << " " << ab << " " << ac
                    //      <<" " << temp << "\n";
                    sactivity = m_Node->DC_a(idx);  // extract activities (not
                                                    // activity coefficients!)
                    //                  cout << "Activity " << sactivity <<
                    //                  pow(10.0,sactivity)<< "\n";
                    aa *= pow(sactivity, m_kin[ii].kinetic_parameters[12 + i]);
                    ab *= pow(sactivity, m_kin[ii].kinetic_parameters[13 + i]);
                    ac *= pow(sactivity, m_kin[ii].kinetic_parameters[14 + i]);
                    //    *Y_m,     // Molalities of aqueous species and
                    //    sorbates [0:Ls-1] *Y_la,    // log activity of DC in
                    //    multi-component phases (mju-mji0) [0:Ls-1] *Y_w, //
                    //    Mass concentrations of DC in multi-component
                    //    phases,%(ppm)[Ls] *Gamma,   // DC activity
                    //    coefficients in molal or other phase-specific scale
                    //    [0:L-1]
                }
            }
            //      cout << "activities " <<aa << " " << ab << " " << ac <<" "
            //      << temp << "\n";
            // terms for each case
            rra = exp(-1.0 * m_kin[ii].kinetic_parameters[0] / R *
                      (1.0 / temp + 1.0 / 298.15)) *
                  aa *
                  pow((1.0 - pow(omega_phase[in * nPH + k],
                                 m_kin[ii].kinetic_parameters[6])),
                      m_kin[ii].kinetic_parameters[7]);

            rrn = exp(-1.0 * m_kin[ii].kinetic_parameters[1] / R *
                      (1.0 / temp + 1.0 / 298.15)) *
                  ab *
                  pow((1.0 - pow(omega_phase[in * nPH + k],
                                 m_kin[ii].kinetic_parameters[8])),
                      m_kin[ii].kinetic_parameters[9]);

            rrb = exp(-1.0 * m_kin[ii].kinetic_parameters[2] / R *
                      (1.0 / temp + 1.0 / 298.15)) *
                  ac *
                  pow((1.0 - pow(omega_phase[in * nPH + k],
                                 m_kin[ii].kinetic_parameters[10])),
                      m_kin[ii].kinetic_parameters[11]);

            // rate is scaled to the total amount available via the surface area
            dmdt[in * nPH + k] =
                -1.0 * sa *
                (pow(10.0, m_kin[ii].kinetic_parameters[3]) * rra +
                 pow(10.0, m_kin[ii].kinetic_parameters[4]) * rrn +
                 pow(10.0, m_kin[ii].kinetic_parameters[5]) * rrb);

            // test for NaN!! ---seems necessary as sometimes rra, rrn, rrb get
            // Inf! ---seems enough to test the upper limit---this test does not
            // resolve the real problem ;-)...probably pow(0.0,0.0) for
            // rra,rrn,rrb ?
            if (!(dmdt[in * nPH + k] <= 1.0) && !(dmdt[in * nPH + k] > 1.0))
            {
                cout << "failed " << m_kin[ii].phase_name << "at node " << in
                     << " dmdt " << dmdt << " is NaN "
                     << " sa " << sa << " rra " << rra << " rrn " << rrn
                     << " rrb " << rrb << "\n";
                cout << "mol_phase " << mol_phase[in * nPH + k] << " m_gam "
                     << m_gam[idx] << " dmdt " << dmdt[in * nPH + k]
                     << " omegaPhase " << omega_phase[in * nPH + k] << "\n";
                dmdt[in * nPH + k] = 0.0;  // no change!
            }
            // calculate max kinetic time step
            if (fabs(dmdt[in * nPH + k]) > 0.0 &&
                fabs(m_xDC[in * nDC + j]) > 0.0)
            {
                rwmutex.lock();  //
                max_kinetic_timestep =
                    min(max_kinetic_timestep,
                        fabs(m_xDC[in * nDC + j] / dmdt[in * nPH + k]));
                rwmutex.unlock();  // avoid mutual exclusion in the MPI version
            }
        }  // end if for kinetic model >1
    }

    return 1;
}

/**
 * REACT_GEM::CalcLimitsInitial ( long in )
 * This is part of the OGS-GEMS kinetic implementation. It should be called
 during initialization phase, when
 * no information from the previous timestep is available. All kinetically
 controlled phases are set to their
 * initial values by assigning dll (lower limit) and dul (upper limit) to the
 xDC values from a restart file or from IC files (in
 * case transport is done with full speciation).
 * In case the simulation starts and no xDC restart values are available, one
 should create such a restart file e.g. by conducting a
 * a equilibrium simulation with one time-step and - if necessary - adjucst the
 values for the kinetically controlled phases e.g.
 * with an editor or by scripts.

 */
int REACT_GEM::CalcLimitsInitial(long in, TNode* m_Node)
{
    long ii, k, j;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
        return 0;

    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
        return 0;
    for (j = 0; j < nDC; j++)
    {
        m_dll[in * nDC + j] = 0.0;      // set to zero
        m_dul[in * nDC + j] = 1.0e+10;  // very high number
    }
    // we test if restart flag is set....if not this will not work, as x_dc
    // might be not correct
    if (m_flow_pcs->GetRestartFlag() == FiniteElement::READ ||
        m_flow_pcs->GetRestartFlag() == FiniteElement::READ_WRITE)
    {
        for (ii = 0; ii < (int)m_kin.size(); ii++)
        {
            k = m_kin[ii].phase_number;

            if (m_kin[ii].kinetic_model >
                0)  // do it only if kinetic model is defined take model

                // kinetic_model==1 dissolution+precipitation kinetics
                // kinetic_model==2 only dissolution (no precipitation)
                // kinetic_mocel==3 only precipitation (no dissolution)

                for (j = m_kin[ii].dc_counter;
                     j < m_kin[ii].dc_counter + dCH->nDCinPH[k]; j++)
                {
                    if ((dCH->ccDC[j] == '0') || (dCH->ccDC[j] == 'X') ||
                        (dCH->ccDC[j] == 'Y') || (dCH->ccDC[j] == 'Z'))
                    {
                        m_dll[in * nDC + j] = 0.0;      // set to zero
                        m_dul[in * nDC + j] = 1.0e+10;  // very high number
                    }
                    else
                    {
                        m_dul[in * nDC + j] = m_xDC[in * nDC + j];
                        m_dll[in * nDC + j] = m_xDC[in * nDC + j];
                        if (m_dll[in * nDC + j] < 0.0)
                            m_dll[in * nDC + j] =
                                0.0;  // no negative masses allowed
                        if (m_dul[in * nDC + j] < 0.0)
                            m_dul[in * nDC + j] = 0.0;
                        if (m_dll[in * nDC + j] > m_dul[in * nDC + j])
                            m_dll[in * nDC + j] = m_dul[in * nDC + j];
                    }
                }
            // end kinetic model
        }          // end loop over phases
        return 1;  // end restart
    }
    return 1;
}

/** In this function we calculate the actual upper and lower metastability
 * constraints for the GEMS solution from the phase reaction rates (calculated
 * at the previous time step)
 */
int REACT_GEM::CalcLimits(long in, TNode* m_Node)
{
    double dummy;
    long ii, k, j;
    DATACH* dCH;  // pointer to DATACH
    DATABR* dBR;
    // Getting direct access to DataCH structure in GEMIPM2K memory
    dCH = m_Node->pCSD();
    if (!dCH)
        return 0;

    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    dBR = m_Node->pCNode();
    if (!dBR)
        return 0;
    for (j = 0; j < nDC; j++)
    {
        m_dll[in * nDC + j] = 0.0;      // set to zero
        m_dul[in * nDC + j] = 1.0e+10;  // very high number
    }

    for (ii = 0; ii < (int)m_kin.size(); ii++)
    {
        k = m_kin[ii].phase_number;
        // cout << " kinetics phase " << ii << " "  << m_kin[ii].kinetic_model
        // << "\n";
        if (m_kin[ii].kinetic_model > 0 &&
            m_kin[ii].kinetic_model <
                6)  // do it only if kinetic model is defined take model

            // kinetic_model==1 dissolution+precipitation kinetics
            // kinetic_model==2 only dissolution (no precipitation)
            // kinetic_mocel==3 only precipitation (no dissolution)
            for (j = m_kin[ii].dc_counter;
                 j < m_kin[ii].dc_counter + dCH->nDCinPH[k]; j++)
            {
                //     cout << "Kin debug " << in << " mol amount species " <<
                //     m_xDC[in*nDC+j] << " saturation phase "
                //     << omega_phase[in*nPH+k] << " saturation species" <<
                //     omega_components[in*nDC+j] << " mol fraction now" <<
                //     (m_xDC[in*nDC+j]/mol_phase[in*nPH+k] )<< "\n";
                // surface complexation species are not kinetically controlled
                // -- 0 is old way...X is new way in DCH files
                if ((dCH->ccDC[j] == '0') || (dCH->ccDC[j] == 'X') ||
                    (dCH->ccDC[j] == 'Y') || (dCH->ccDC[j] == 'Z'))
                {
                    m_dll[in * nDC + j] = 0.0;      // set to zero
                    m_dul[in * nDC + j] = 1.0e+10;  // very high number
                }
                else
                {
                    // cout << "Kin debug SS mol phase " << mol_phase[in*nPH+k]
                    // << " dmdt " << dmdt[in*nPH+k]*dt << " omega comp "
                    // <<omega_components[in*nDC+j] <<" omega phase " <<
                    // omega_phase[in*nPH+k]<< "\n";
                    dummy =
                        (mol_phase[in * nPH + k] + dmdt[in * nPH + k] * dt) *
                        omega_components[in * nDC + j] /
                        omega_phase[in * nPH + k];

                    if (m_kin[ii].kinetic_model ==
                        5)  // only for Solid solution models: rescale to change
                            // of endmember
                        // in case of Vanselow convenction or similar

                        dummy /= m_kin[ii].ss_scaling[j - m_kin[ii].dc_counter];
                    if (dummy > 1.0e10)
                        dummy = 1.0e10;

                    if (!(dummy <= 1.0) && !(dummy > 1.0))
                    {
                        // no change!
                        m_dul[in * nDC + j] = m_xDC[in * nDC + j];
                        m_dll[in * nDC + j] = m_xDC[in * nDC + j];
                    }
                    else if (m_xDC[in * nDC + j] >
                             dummy)  // This is the dissolution case
                    {
                        m_dul[in * nDC + j] = m_xDC[in * nDC + j];
                        m_dll[in * nDC + j] = dummy;
                    }
                    else  // This is the precipitation case
                    {
                        m_dll[in * nDC + j] = m_xDC[in * nDC + j];
                        m_dul[in * nDC + j] = dummy;
                    }
                    // do some corrections
                    // kinetic_model==2 only dissolution is kontrolled (free
                    // precipitation) kinetic_mocel==3 only precipitation is
                    // copntroleld (free dissolution)
                    if ((m_kin[ii].kinetic_model == 2) &&
                        (m_dul[in * nDC + j] > m_xDC[in * nDC + j]))
                        m_dul[in * nDC + j] =
                            1.0e+10;  // m_dul[in*nDC+j]= m_xDC[in*nDC+j];
                    if ((m_kin[ii].kinetic_model == 3) &&
                        (m_dll[in * nDC + j] < m_xDC[in * nDC + j]))
                        m_dll[in * nDC + j] =
                            0.0;  // m_dll[in*nDC+j]= m_xDC[in*nDC+j];
                    if ((m_xDC[in * nDC + j] < 1.0e-6) &&
                        (omega_phase[in * nPH + k] >= 1.0001) &&
                        (m_dul[in * nDC + j] < 1.0e-6))
                    {
                        m_dul[in * nDC + j] =
                            1.0e-6;  // allow some kind of precipitation...based
                                     // on saturation index
                        // for component value...here we set 10-6 mol per m^3
                        // ..which is maybe 10-10 per litre ...?
                        m_dll[in * nDC + j] = 0.0;
                    }
                    if (m_dll[in * nDC + j] > m_dul[in * nDC + j])
                        m_dll[in * nDC + j] =
                            m_dul[in * nDC +
                                  j];  // dll should be always lower than dul
                    // no negative masses allowed
                    if (m_dll[in * nDC + j] < 0.0)
                        m_dll[in * nDC + j] = 0.0;
                    // no negative masses allowed..give some freedom
                    if (m_dul[in * nDC + j] <= 0.0)
                    {
                        m_dul[in * nDC + j] = 1.0e-6;
                        m_dll[in * nDC + j] = 0.0;
                    }
                }
                // cout << "Kin debug for component no. " << j << " at node " <<
                // in << " m_xDC "  <<  m_xDC[in*nDC+j] << " m_dll, mdul " <<
                // m_dll[in*nDC+j] << " " << m_dul[in*nDC+j] << " diff " <<
                // m_dul[in*nDC+j]- m_dll[in*nDC+j] << "\n";
                //            if ((fabs((m_dul[in*nDC+j]-
                //            m_dll[in*nDC+j]))>0.0)) cout << "Kinetics for
                //            component no. "
                //            << j << " at node " << in << " m_xDC "  <<
                //            m_xDC[in*nDC+j] << " m_dll, mdul " <<
                //            m_dll[in*nDC+j] << " " << m_dul[in*nDC+j] << "
                //            diff " << m_dul[in*nDC+j]- m_dll[in*nDC+j]
                //            << "\n"; // give some debug output for kinetics
            }

    }  // end loop over phases

    return 1;
}

// simplest case....scaling with a specific surface area per volume mineral
// phasenr: index for phase compnr: index of component which belongs to the
// phase and for which a specific surface area is defined
double REACT_GEM::SurfaceAreaPh(long kin_phasenr, long in, TNode* m_Node)
{
    double surf_area = 0.0;

    // now it is volume of the phase ..with respect to the unit volume in m^3
    surf_area = m_Node->Ph_Volume(m_kin[kin_phasenr].phase_number);

    if (m_kin[kin_phasenr].surface_model == 1)
        // multiplication with specific surface area gives area
        surf_area *= m_kin[kin_phasenr].surface_area[0];
    else if (m_kin[kin_phasenr].surface_model == 2)
        // constant surface area
        surf_area = m_kin[kin_phasenr].surface_area[0];
    else if (m_kin[kin_phasenr].surface_model == 3)
        surf_area *= m_kin[kin_phasenr].surface_area[0] /
                     m_porosity[in];  // multiplication with specific surface
                                      // area and division by porosity
    else if (m_kin[kin_phasenr].surface_model == 4)
        surf_area = m_kin[kin_phasenr].surface_area[0] *
                    m_porosity[in];  // multiplication of specific surface area
                                     // and  porosity
    else
        surf_area = 0.0;  // no kinetics...safe solution

    // cout << "phase " << kin_phasenr << " area default " <<
    // m_kin[kin_phasenr].surface_area[0] << " model " <<
    // m_kin[kin_phasenr].surface_model <<  " surface area: " << surf_area <<
    // "\n";
    return surf_area;
}

#if defined(USE_MPI_GEMS)
void REACT_GEM::GetGEMResult_MPI(void)
{
    // Now gather the calculated
    // values------------------------------------------------------------------------------
    MPI_Allreduce(m_NodeHandle_buff, m_NodeHandle, nNodes, MPI_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_NodeStatusCH_buff, m_NodeStatusCH, nNodes, MPI_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_IterDone_buff, m_IterDone, nNodes, MPI_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(m_Vs_buff, m_Vs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_Ms_buff, m_Ms, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_Gs_buff, m_Gs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_Hs_buff, m_Hs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_IC_buff, m_IC, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_pH_buff, m_pH, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_pe_buff, m_pe, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_Eh_buff, m_Eh, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_porosity_buff, m_porosity, nNodes, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_fluid_volume_buff, m_fluid_volume, nNodes, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_fluid_density_buff, m_fluid_density, nNodes, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_gas_volume_buff, m_gas_volume, nNodes, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_excess_water_buff, m_excess_water, nNodes, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_excess_gas_buff, m_excess_gas, nNodes, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_bPS_buff, m_bPS, nNodes * nIC * nPS, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(m_bIC_buff, m_bIC, nNodes * nIC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_bIC_dummy_buff, m_bIC_dummy, nNodes * nIC, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_soluteB_buff, m_soluteB, nNodes * nIC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(m_xDC_buff, m_xDC, nNodes * nDC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(m_dll_buff, m_dll, nNodes * nDC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_dul_buff, m_dul, nNodes * nDC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_aPH_buff, m_aPH, nNodes * nPH, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_xPH_buff, m_xPH, nNodes * nPH, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_xPA_buff, m_xPA, nNodes * nPS, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(m_xDC_pts_buff, m_xDC_pts, nNodes * nDC, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(m_xDC_MT_delta_buff, m_xDC_MT_delta, nNodes * nDC, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(m_xDC_Chem_delta_buff, m_xDC_Chem_delta, nNodes * nDC,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(omega_phase_buff, omega_phase, nNodes * nPH, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(mol_phase_buff, mol_phase, nNodes * nPH, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(dmdt_buff, dmdt, nNodes * nPH, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(omega_components_buff, omega_components, nNodes * nDC,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // --------------------------------------------------------------------------------------------------------------
}

void REACT_GEM::CleanMPIBuffer(void)
{
    long in, ii, jj;
    for (in = 0; in < nNodes; in++)
    {
        m_IterDoneCumulative[in] += m_IterDone[in];
        m_NodeHandle_buff[in] = 0;
        m_NodeStatusCH_buff[in] = 0;
        m_IterDone_buff[in] = 0;

        m_Vs_buff[in] = 0.0;
        m_Ms_buff[in] = 0.0;
        m_Gs_buff[in] = 0.0;
        m_Hs_buff[in] = 0.0;
        m_IC_buff[in] = 0.0;
        m_pH_buff[in] = 0.0;
        m_pe_buff[in] = 0.0;
        m_Eh_buff[in] = 0.0;
        m_porosity_buff[in] = 0.0;
        m_fluid_volume_buff[in] = 0.0;
        m_fluid_density_buff[in] = 0.0;
        m_gas_volume_buff[in] = 0.0;

        m_excess_water_buff[in] = 0.0;
        m_excess_gas_buff[in] = 0.0;

        for (ii = 0; ii < nIC; ii++)
        {
            m_soluteB_buff[in * nIC + ii] = 0.0;
            m_bIC_buff[in * nIC + ii] = 0.0;
            m_bIC_dummy_buff[in * nIC + ii] = 0.0;
        }

        for (ii = 0; ii < nDC; ii++)
        {
            m_xDC_buff[in * nDC + ii] = 0.0;
            m_dul_buff[in * nDC + ii] = 0.0;
            m_dll_buff[in * nDC + ii] = 0.0;
            m_xDC_pts_buff[in * nDC + ii] = 0.0;
            m_xDC_MT_delta_buff[in * nDC + ii] = 0.0;
            m_xDC_Chem_delta_buff[in * nDC + ii] = 0.0;
            omega_components_buff[in * nDC + ii] = 0.0;
        }

        for (ii = 0; ii < nPH; ii++)
        {
            m_xPH_buff[in * nPH + ii] = 0.0;
            m_aPH_buff[in * nPH + ii] = 0.0;
            omega_phase_buff[in * nPH + ii] = 0.0;
            mol_phase_buff[in * nPH + ii] = 0.0;
            dmdt_buff[in * nPH + ii] = 0.0;
        }

        for (ii = 0; ii < nPS; ii++)
            m_xPA_buff[in * nPS + ii] = 0.0;

        for (ii = 0; ii < nIC; ii++)
            for (int jj = 0; jj < nPS; jj++)
                m_bPS_buff[in * ii * nPS + jj] = 0.0;
    }

    for (in = 0; in < nElems; in++)
        m_porosity_Elem_buff[in] = m_porosity_Elem[in];
}

void REACT_GEM::CopyToMPIBuffer(long in)
{
    long ii;
    m_NodeHandle_buff[in] = m_NodeHandle[in];
    m_NodeStatusCH_buff[in] = m_NodeStatusCH[in];
    m_IterDone_buff[in] = m_IterDone[in];

    m_Vs_buff[in] = m_Vs[in];
    m_Ms_buff[in] = m_Ms[in];
    m_Gs_buff[in] = m_Gs[in];
    m_Hs_buff[in] = m_Hs[in];
    m_IC_buff[in] = m_IC[in];
    m_pH_buff[in] = m_pH[in];
    m_pe_buff[in] = m_pe[in];
    m_Eh_buff[in] = m_Eh[in];
    m_porosity_buff[in] = m_porosity[in];
    m_fluid_volume_buff[in] = m_fluid_volume[in];
    m_fluid_density_buff[in] = m_fluid_density[in];

    m_gas_volume_buff[in] = m_gas_volume[in];

    m_excess_water_buff[in] = m_excess_water[in];
    m_excess_gas_buff[in] = m_excess_gas[in];

    for (ii = 0; ii < nIC; ii++)
    {
        m_bIC_buff[in * nIC + ii] = m_bIC[in * nIC + ii];
        m_soluteB_buff[in * nIC + ii] = m_soluteB[in * nIC + ii];

        m_bIC_dummy_buff[in * nIC + ii] = m_bIC_dummy[in * nIC + ii];
    }

    for (ii = 0; ii < nDC; ii++)
    {
        m_xDC_buff[in * nDC + ii] = m_xDC[in * nDC + ii];
        m_dul_buff[in * nDC + ii] = m_dul[in * nDC + ii];
        m_dll_buff[in * nDC + ii] = m_dll[in * nDC + ii];
        m_xDC_pts_buff[in * nDC + ii] = m_xDC_pts[in * nDC + ii];
        m_xDC_MT_delta_buff[in * nDC + ii] = m_xDC_MT_delta[in * nDC + ii];
        m_xDC_Chem_delta_buff[in * nDC + ii] = m_xDC_Chem_delta[in * nDC + ii];
        omega_components_buff[in * nDC + ii] = omega_components[in * nDC + ii];
    }

    for (ii = 0; ii < nPH; ii++)
    {
        m_xPH_buff[in * nPH + ii] = m_xPH[in * nPH + ii];
        m_aPH_buff[in * nPH + ii] = m_aPH[in * nPH + ii];
        omega_phase_buff[in * nPH + ii] = omega_phase[in * nPH + ii];
        mol_phase_buff[in * nPH + ii] = mol_phase[in * nPH + ii];
        dmdt_buff[in * nPH + ii] = dmdt[in * nPH + ii];
    }

    for (ii = 0; ii < nPS; ii++)
        m_xPA_buff[in * nPS + ii] = m_xPA[in * nPS + ii];
}

#endif  // end MPI

/** GetNodePorosityValue(long node_Index): function coming from GEMS
 * coupling...extract node based porosities..does only work with GEMS coupling
 * georg.kosakowski@psi.ch 02.11.2009
 * georg.kosakowski@psi.ch 22.02.2013
 */
double REACT_GEM::GetNodePorosityValue(long node_Index)
{
    double node_poros =
        -1.0;  // default value is negative...that should not happen

    node_poros = m_porosity[node_Index];

    return node_poros;
}

/** GetNodePorosityValueInitial(long node_Index): function coming from GEMS
 * coupling...extract initial node based porosities..does only work with GEMS
 * coupling georg.kosakowski@psi.ch 02.11.2009
 * georg.kosakowski@psi.ch 22.02.2013
 */
double REACT_GEM::GetNodePorosityValueInitial(long node_Index)
{
    double node_poros =
        -1.0;  // default value is negative...that should not happen

    node_poros = m_porosity_initial[node_Index];

    return node_poros;
}

// function coming from GEMS coupling...extract node based fluid phase
// densities..does only work with GEMS coupling a
// georg.kosakowski@psi.ch 15.06.2012

double REACT_GEM::GetNodeFluidDensityValue(long node_Index)
{
    double fluid_density =
        -1.0;  // default value is negative...that should not happen

    fluid_density = m_fluid_density[node_Index];

    return fluid_density;
}

// taken from rf_REACT_BRNS
int REACT_GEM::IsThisPointBCIfYesStoreValue(long index, CRFProcess* m_pcs,
                                            double& value)
{
    for (long p = 0; p < (int)m_pcs->bc_node_value.size(); ++p)
        if (index == m_pcs->bc_node_value[p]->msh_node_number)
        {
            value = m_pcs->bc_node_value[p]->node_value;
            return 1;  // Yes, found it.
        }

    return 0;
}

int REACT_GEM::WriteReloadGem()
{
    long i, j;

    string rank_str;
    rank_str = "0";
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
    int rank, msize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &msize);
    stringstream ss(stringstream::in | stringstream::out);
    ss.clear();
    ss.str("");
    ss << rank;
    rank_str = ss.str();
    ss.clear();
    // write out a separate time stamp into file
    string m_file_namet = FileName + "_" + "time" + "_gem_" + rank_str + ".asc";
#else  // if defined(USE_PETSC)
    string m_file_namet = FileName + "_" + "time" + "_gem.asc";
#endif

    ofstream ost(m_file_namet.c_str(), ios::trunc | ios::out);
    if (!ost.good())
    {
        cout << "Failure to open file: " << m_file_namet << "\n";
        abort();
    }
    ost.precision(15);  // 15 digits accuracy seems enough? more fields are
                        // filled up with random numbers!
    ost.setf(ios_base::scientific, ios_base::floatfield);
    ost << "time: " << aktuelle_zeit << "\n";
    ost.close();

    // node porosity ...necessary for restart
    // first test if m_poorosity exists (possible if no gem process exists)
    if (!m_porosity)
        return 2;
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
    string m_file_name =
        FileName + "_" + "m_porosity" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name = FileName + "_" + "m_porosity" + "_gem.asc";
#endif
    ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
    if (!os.good())
    {
        cout << "Failure to open file: " << m_file_name << "\n";
        abort();
    }
    os.precision(15);  // 15 digits accuracy seems enough? more fields are
                       // filled up with random numbers!
    os.setf(ios_base::scientific, ios_base::floatfield);

    for (i = 0; i < nNodes; i++)
    {
        os << m_porosity[i] << "  ";
        os << "\n";
    }
    os.close();
    // now write the fluid_volumes
    if (!m_fluid_volume)
        return 2;
#if defined(USE_PETSC)
    string m_file_name2 =
        FileName + "_" + "m_fluid_volume" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name2 = FileName + "_" + "m_fluid_volume" + "_gem.asc";
#endif
    ofstream os2(m_file_name2.c_str(), ios::trunc | ios::out);
    if (!os2.good())
    {
        cout << "Failure to open file: " << m_file_name2 << "\n";
        abort();
    }
    os2.precision(15);  // 15 digits accuracy seems enough? more fields are
                        // filled up with random numbers!
    os2.setf(ios_base::scientific, ios_base::floatfield);

    for (i = 0; i < nNodes; i++)
    {
        os2 << m_fluid_volume[i] << "  ";
        os2 << "\n";
    }
    os2.close();

    // write b vector
    if (!m_bIC)
        return 2;
#if defined(USE_PETSC)
    string m_file_name3 =
        FileName + "_" + "m_bIC" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name3 = FileName + "_" + "m_bIC" + "_gem.asc";
#endif
    ofstream os3(m_file_name3.c_str(), ios::trunc | ios::out);
    if (!os3.good())
    {
        cout << "Failure to open file: " << m_file_name3 << "\n";
        abort();
    }
    os3.precision(15);  // 15 digits accuracy seems enough? more fields are
                        // filled up with random numbers!
    os3.setf(ios_base::scientific, ios_base::floatfield);

    for (i = 0; i < nNodes; i++)
        for (j = 0; j < nIC; j++)
        {
            os3 << m_bIC[i * nIC + j] << "  ";
            os3 << m_soluteB[i * nIC + j] << "  ";
            os3 << "\n";
        }
    os3.close();
    // write xDC vector
    if (!m_xDC)
        return 2;
#if defined(USE_PETSC)
    string m_file_name4 =
        FileName + "_" + "m_xDC" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name4 = FileName + "_" + "m_xDC" + "_gem.asc";
#endif
    ofstream os4(m_file_name4.c_str(), ios::trunc | ios::out);
    if (!os4.good())
    {
        cout << "Failure to open file: " << m_file_name4 << "\n";
        abort();
    }
    os4.precision(15);  // 15 digits accuracy seems enough? more fields are
                        // filled up with random numbers!
    os4.setf(ios_base::scientific, ios_base::floatfield);

    for (i = 0; i < nNodes; i++)
        for (j = 0; j < nDC; j++)
        {
            os4 << m_xDC[i * nDC + j] << "  " << m_dll[i * nDC + j] << "  "
                << m_dul[i * nDC + j];
            os4 << "\n";
        }
    os4.close();

    return 1;
}

int REACT_GEM::ReadReloadGem()
{
    long i, j;

    string rank_str;
    rank_str = "0";
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
    int rank, msize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &msize);
    stringstream ss(stringstream::in | stringstream::out);
    ss.clear();
    ss.str("");
    ss << rank;
    rank_str = ss.str();
    ss.clear();
    string m_file_name =
        FileName + "_" + "m_porosity" + "_gem_" + rank_str + ".asc";
#else  // if defined(USE_PETSC)
    string m_file_name = FileName + "_" + "m_porosity" + "_gem.asc";
#endif
    ifstream is(m_file_name.c_str(), ios::in);
    if (!is.good())
    {
        cout << "Failure to open file: " << m_file_name << "\n";
        abort();
    }
    //
    // node porosity ...necessary for restart

    for (i = 0; i < nNodes; i++)
    {
        is >> m_porosity[i];
        is >> ws;
    }
    is.close();
#if defined(USE_PETSC)
    string m_file_name2 =
        FileName + "_" + "m_fluid_volume" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name2 = FileName + "_" + "m_fluid_volume" + "_gem.asc";
#endif
    ifstream is2(m_file_name2.c_str(), ios::in);
    if (!is2.good())
    {
        cout << "Failure to open file: " << m_file_name2 << "\n";
        abort();
    }

    for (i = 0; i < nNodes; i++)
    {
        is2 >> m_fluid_volume[i];
        is2 >> ws;
    }
    is2.close();
#if defined(USE_PETSC)
    string m_file_name3 =
        FileName + "_" + "m_bIC" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name3 = FileName + "_" + "m_bIC" + "_gem.asc";
#endif
    ifstream is3(m_file_name3.c_str(), ios::in);
    if (!is3.good())
    {
        cout << "Failure to open file: " << m_file_name3 << "\n";
        abort();
    }

    for (i = 0; i < nNodes; i++)
        for (j = 0; j < nIC; j++)
        {
            is3 >> m_bIC[i * nIC + j];
            is3 >> m_soluteB[i * nIC + j];
            is3 >> ws;
        }
    is3.close();
#if defined(USE_PETSC)
    string m_file_name4 =
        FileName + "_" + "m_xDC" + "_gem_" + rank_str + ".asc";
#else
    string m_file_name4 = FileName + "_" + "m_xDC" + "_gem.asc";
#endif
    ifstream is4(m_file_name4.c_str(), ios::in);
    if (!is4.good())
    {
        cout << "Failure to open file: " << m_file_name4 << "\n";
        abort();
    }

    for (i = 0; i < nNodes; i++)
        for (j = 0; j < nDC; j++)
        {
            is4 >> m_xDC[i * nDC + j];
            is4 >> m_dll[i * nDC + j];
            is4 >> m_dul[i * nDC + j];
            is4 >> ws;
        }
    is4.close();

    return 1;
}

#if defined(USE_PETSC)
/** WriteVTKGEMValues (fstream &vtk_file) ...for serial or normal mpi version it
 *appends GEMS node based data to vtk output for PETSC version it writes its own
 *vtk file
 *
 */
void REACT_GEM::WriteVTKGEMValuesPETSC(PetscViewer viewer)
{
    PetscScalar* xp;  // used for pointer
    PetscInt low, high;
    PetscInt count;
    double bdummy;
    // MeshLib::CFEMesh *mesh = fem_msh_vector[0];
    // const int nn = mesh->getNumNodesGlobal(); //global number of nodes
    // ..without shadow nodes cout << "DEBUG nNodes nn" << nNodes << " " << nn
    // <<" \n"; test vtk output
    PETSc_Vec x;  //
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, nNodes, PETSC_DECIDE);
    VecSetFromOptions(x);  //
    // get range of local variables
    VecGetOwnershipRange(x, &low, &high);
    VecGetLocalSize(x, &count);
    // get local part of vectors
    VecGetArray(x, &xp);
    // write node based data
    PetscViewerSetFormat(
        viewer, PETSC_VIEWER_ASCII_VTK);  // will cause a failure if cell data
                                          // is written before this

    //******************************** B Vector ********************
    for (int j = 0; j < nIC; j++)
    {
        const char* fieldname = m_ICNL[j];
        PetscObjectSetName((PetscObject)x, fieldname);
        //....................................................................
        for (int k = 0; k < count; k++)
        {
            bdummy = m_soluteB[k * nIC + j] *
                     m_fluid_volume[k];  // soluteB contains volume based
                                         // concentrations
            // now we have to add water
            if (idx_hydrogen == j)
                bdummy +=
                    (2.0 *
                     m_xDC[k * nDC + idx_water]);  // carrier for zero(first
                                                   // phase)  is normally water!
            else if (idx_oxygen == j)
                bdummy += m_xDC[k * nDC + idx_water];

            bdummy += m_bIC[k * nIC + j];  // add the solids
            xp[k] = bdummy;
        }
        VecView(x, viewer);
    }
    //***********************    speciation vector
    //*********************************
    for (int i = 0; i < nDC; i++)
    {
        const char* fieldname = m_DCNL[i];
        PetscObjectSetName((PetscObject)x, fieldname);
        //....................................................................
        for (int j = 0; j < nNodes; j++)
            xp[j] = m_xDC[j * nDC + i];
        VecView(x, viewer);
    }

    // ******************* eh, pe, pH, Nodeporosity ****************************
    PetscObjectSetName((PetscObject)x, "Ph");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_pH[j];
    VecView(x, viewer);

    PetscObjectSetName((PetscObject)x, "pe");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_pe[j];
    VecView(x, viewer);

    PetscObjectSetName((PetscObject)x, "Eh");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_Eh[j];
    VecView(x, viewer);
    // ***************************Nodeporosity, fluidvolume, excess volume, node
    // volume

    PetscObjectSetName((PetscObject)x, "NodePorosity");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_porosity[j];
    VecView(x, viewer);

    PetscObjectSetName((PetscObject)x, "FluidVolume");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_fluid_volume[j];
    VecView(x, viewer);

    PetscObjectSetName((PetscObject)x, "ExcessVolume");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_excess_water[j];
    VecView(x, viewer);

    PetscObjectSetName((PetscObject)x, "NodeVolume");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_Vs[j];
    VecView(x, viewer);

    // *********************Fluiddensity  ****************************

    PetscObjectSetName((PetscObject)x, "FluidDensity");
    //....................................................................
    for (int j = 0; j < nNodes; j++)
        xp[j] = m_fluid_density[j];
}

#endif

void REACT_GEM::WriteVTKGEMValues(fstream& vtk_file)
{
    long i, j, k;
    double bdummy = 0.0;
    // this is point data

    for (j = 0; j < nIC; j++)
    {
        vtk_file << "SCALARS " << m_ICNL[j] << " double 1"
                 << "\n";
        vtk_file << "LOOKUP_TABLE default"
                 << "\n";
        //....................................................................
        for (k = 0; k < nNodes; k++)
        {
            bdummy = m_soluteB[k * nIC + j] *
                     m_fluid_volume[k];  // soluteB contains volume based
                                         // concentrations
            // now we have to add water
            if (idx_hydrogen == j)
                bdummy +=
                    (2.0 *
                     m_xDC[k * nDC + idx_water]);  // carrier for zero(first
                                                   // phase)  is normally water!
            else if (idx_oxygen == j)
                bdummy += m_xDC[k * nDC + idx_water];

            bdummy += m_bIC[k * nIC + j];       // add the solids
            vtk_file << " " << bdummy << "\n";  // and output
        }
    }

    // loop over speciation vector!
    for (i = 0; i < nDC; i++)
    {
        vtk_file << "SCALARS " << m_DCNL[i] << " double 1"
                 << "\n";
        vtk_file << "LOOKUP_TABLE default"
                 << "\n";
        //....................................................................
        for (j = 0; j < nNodes; j++)
            vtk_file << " " << (float)m_xDC[j * nDC + i] << "\n";
    }
    // eh, pe, pH, Nodeporosity
    vtk_file << "SCALARS "
             << " pH "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_pH[j] << "\n";
    vtk_file << "SCALARS "
             << " pe "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_pe[j] << "\n";
    vtk_file << "SCALARS "
             << " Eh "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_Eh[j] << "\n";
    vtk_file << "SCALARS "
             << " NodePorosity "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_porosity[j] << "\n";
    vtk_file << "SCALARS "
             << " Fluidvolume "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_fluid_volume[j] << "\n";
    vtk_file << "SCALARS "
             << " ExcessVolume "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_excess_water[j] << "\n";
    vtk_file << "SCALARS "
             << " NodeVolume "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_Vs[j] << "\n";
    vtk_file << "SCALARS "
             << " FluidDensity "
             << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";
    //....................................................................
    for (j = 0; j < nNodes; j++)
        vtk_file << " " << m_fluid_density[j] << "\n";
}

/** gems_worker:
 * This function creates and initializes one GEMS3K kernel. It may be spawned
 * several times. Calculations are coordinated with boost barrier commands.
 */
void REACT_GEM::gems_worker(int tid, string m_Project_path)
{
    nNodes = GetNodeNumber_MT();
    nElems = GetElemNumber_MT();
    long j, in, node_fail = 0, repeated_fail = 0;
    double oldvolume;

    // collection of time data for rough performance analysis
    double time_gem_total, time_fraction, tdummy, tdummy1, tdummy2, twait,
        tstart, twaittotal;
    time_gem_total = 0.0;
    time_fraction = 0.0;
    tdummy = 0.0;
    twaittotal = 0.0;
    tstart = GetTimeOfDayDouble();

    // create data necessary for running gems in thread
    TNode* t_Node;

    // DATABR structure for exchange with GEMIPM this should be local data for
    // each thread
    DATACH* tdCH;  // pointer to DATACH
    DATABR* tdBR;  // pointer to DATABR

    string tinit_path;
    char* buffer = NULL;

    tinit_path = m_Project_path.append(REACT_GEM::init_input_file_path);
    // int tid=1;
    //    string tinit_path = "./BC-dat.lst";

#ifdef _WIN32
    // keep this on windows
    if (tinit_path.rfind("\\") == string::npos)
#else
    // keep this on linux
    if (tinit_path.rfind("/") == string::npos)
#endif
    {
#ifdef _WIN32
        if ((buffer = _getcwd(NULL, 0)) == NULL)
#else
        if ((buffer = getcwd(NULL, 0)) == NULL)
#endif
            perror("_getcwd error");
        else
        {
#ifdef _WIN32
            tinit_path.insert(0, "\\");  // keep this on window
#else
            tinit_path.insert(0, "/");  // keep this on linux
#endif
            tinit_path.insert(0, buffer);
        }
    }

    if (buffer)
        free(buffer);

    rwmutex.lock();  // first lock for reading gem init files

    t_Node = new TNode();
    if (t_Node->GEM_init(tinit_path.c_str()))
    {
        // error occured during reading the files
        cout << "Hello World! I failed!!!! It's me, thread " << tid << " "
             << tinit_path.c_str() << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif

        exit(1);
    }

    tdCH = t_Node->pCSD();
    if (!tdCH)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif

        exit(1);
    }
    // Getting direct access to work node DATABR structure which
    // exchanges data between GEMIPM and FMT parts
    tdBR = t_Node->pCNode();
    if (!tdBR)
    {
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
        MPI_Finalize();  // make sure MPI exits
#endif

        exit(1);
    }
    // run GEMS once
    tdBR->NodeStatusCH = NEED_GEM_AIA;
    t_Node->GEM_run(false);

    //    t_Node->GEM_write_dbr ( "dbr_for_crash_node_init_thread1.txt" );
    //    t_Node->GEM_print_ipm ( "ipm_for_crash_node_init_thread1.txt" );

    rwmutex.unlock();
    // now we set the loop variables
    int mycount, mystart;
#if defined(USE_MPI_GEMS)
    mystart = tid + (myrank * gem_nThread);
    mycount = gem_nThread * mysize;

//    rwmutex.lock();   //avoid mutual exclusion in the MPI version
//    cout << "GEMS3K MPI Processe / Thread: " << myrank << " " << tid << "
//    mystart,mycont: " << mystart << " " << mycount << "\n"; rwmutex.unlock();

// here "myrank" is the index of the CPU Processes, and "mysize" is the number
// of CPU Processes
#else
    mycount = gem_nThread;  // make sure the loop over the nodes counts only the
                            // threads!
    mystart = tid;
    int myrank = 0;
    int mysize = 1;
#endif

    //    cout << "thread " << tid << " waits for start barrier " << "\n";
    gem_barrier_start->wait();  // wait for init run start
    //    cout << "thread " << tid << " passed start barrier " << "\n";
    tdummy = GetTimeOfDayDouble();

    for (in = mystart; in < nNodes;
         in += mycount)  // myrank ist defined vi USE_MPI_GEMS
    {
        //       rwmutex.lock();
        //        cout << "GEMS3K MPI Processe / Thread: " << myrank << " " <<
        //        tid << " in " << in << "\n"; rwmutex.unlock();
        // everything is stored in concentrations for restart ...moved it to
        // here from init_gems
        if (m_flow_pcs->GetRestartFlag() == FiniteElement::READ ||
            m_flow_pcs->GetRestartFlag() == FiniteElement::READ_WRITE)
            // Convert from concentration
            REACT_GEM::ConcentrationToMass(
                in, 1);  // I believe this is save for MPI
        // this we have already

        // now we calculate kinetic constraints for GEMS!
        if (m_flow_pcs->GetRestartFlag() == FiniteElement::WRITE ||
            m_flow_pcs->GetRestartFlag() == FiniteElement::NO_IO)
            REACT_GEM::CalcLimitsInitial(
                in, t_Node);  // kg44 16.05.2013 new version, restart files
                              // contain upper and lower limits
        // Manipulate some kinetic contstraints for special initial conditions
        // get all nodes of mesh
        //	const std::vector<MeshLib::CNode*>& msh_nodes
        //(m_flow_pcs->m_msh->getNodeVector()); 	double const* const coords
        //(msh_nodes[in]->getData()); 	if (coords[0] <=0.0)
        //m_dll[in*nDC+nDC-3]=27300.0;  // here set monocorn... 	if ((coords[0]
        //>0.0) && (coords[0]<10.0)) m_dll[in*nDC+nDC-6]= 2.869412109375e+04; //
        //here set gravel
        //        if (coords[0] >=10.0) m_dll[in*nDC+nDC-3]=2.532753515625e+04;
        //        // here set monocorn again

        // Get data
        REACT_GEM::SetReactInfoBackGEM(
            in, t_Node);  // thcc1p.pcs.initialis is necessary, otherwise the
                          // correct data is not available

        // Order GEM to run
        tdBR->NodeStatusCH = NEED_GEM_AIA;
        m_NodeStatusCH[in] = t_Node->GEM_run(false);
        //                t_Node->GEM_write_dbr (
        //                "dbr_for_crash_node_init_thread_1.txt" );

        REACT_GEM::GetReactInfoFromGEM(
            in, t_Node);  // test case..get the data even if GEMS failed

        if (!(m_NodeStatusCH[in] == OK_GEM_AIA ||
              m_NodeStatusCH[in] == OK_GEM_SIA))
        {
            rwmutex.lock();
            //    cout << "Initial GEMs run first pass failed at node " << in ;
            //    cout  << " repeat calculations and change kinetic
            //    constraintsby 0.1%" << "\n";
            rwmutex.unlock();
            // change a bit the kinetic constraints -> make the system less
            // stiff
            for (j = 0; j < nDC; j++)
            {
                m_dll[in * nDC + j] =
                    0.999 * m_dll[in * nDC + j] - 1.0e-6;  // make smaller
                if (m_dll[in * nDC + j] < 0.0)
                    m_dll[in * nDC + j] = 0.0;
                m_dul[in * nDC + j] =
                    1.001 * m_dul[in * nDC + j] + 1.0e-6;  // make bigger
            }
            REACT_GEM::SetReactInfoBackGEM(in, t_Node);  // needs to be done to

            // run GEMS again
            tdBR->NodeStatusCH = NEED_GEM_AIA;
            m_NodeStatusCH[in] = t_Node->GEM_run(false);

            if ((m_NodeStatusCH[in] == ERR_GEM_AIA ||
                 m_NodeStatusCH[in] == ERR_GEM_SIA))
            {
                /*
                          rwmutex.lock();
                                cout << " Error: Init Loop failed when running
                   GEM on Node #" << in << "." << "\n"; cout << "Returned Error
                   Code: " << m_NodeStatusCH[in] << "\n"; t_Node->GEM_write_dbr
                   ( "dbr_for_crash_node_init_thread.txt" );
                                t_Node->GEM_print_ipm (
                   "ipm_for_crash_node_init_thread.txt" ); rwmutex.unlock(); #if
                   defined(USE_MPI_GEMS) || defined(USE_PETSC) MPI_Finalize();
                   //make sure MPI exits #endif

                                exit ( 1 );
                 */
                rwmutex.lock();
                cout << "error: Initial GEMs run after Read GEMS gives bad "
                        "result..proceed in any case. node: "
                     << in << "\n";
                rwmutex.unlock();
            }
            else if ((m_NodeStatusCH[in] == BAD_GEM_AIA ||
                      m_NodeStatusCH[in] == BAD_GEM_SIA))
            {
                rwmutex.lock();
                cout << "error: Initial GEMs run after Read GEMS gives bad "
                        "result..proceed in any case. node: "
                     << in << "\n";
                rwmutex.unlock();
            }
            else
            {
                //              rwmutex.lock();
                //              cout << " sucess with second try.... "<<  "\n";
                //              rwmutex.unlock();
            }
        }  // end loop if initial gems run fails

        // Get data also for volumes!
        REACT_GEM::GetReactInfoFromGEM(
            in, t_Node);  // this we need also for restart
        // calculate density of fluid phase, which is normally the first phase
        m_fluid_density[in] = m_mPS[in * nPS + 0] / m_vPS[in * nPS + 0];

        // we do not need the second pass for complete restart
        if (m_flow_pcs->GetRestartFlag() == FiniteElement::WRITE ||
            m_flow_pcs->GetRestartFlag() == FiniteElement::NO_IO)
            // scale data so that second pass gives the normalized volume of
            // 1m^3
            if (m_Vs[in] >= 0.0)
            {  // this should be not done for restart,  decoupled porosity runs
               // do not change volumes and the other runs
                // should be ok
                for (j = 0; j < nIC; j++)
                    m_bIC[in * nIC + j] /=
                        m_Vs[in];  // This is then for b vector
                for (j = 0; j < nDC; j++)
                {
                    m_dll[in * nDC + j] /= m_Vs[in];
                    m_dul[in * nDC + j] /= m_Vs[in];
                }
            }
        // end if for restart

        REACT_GEM::SetReactInfoBackGEM(
            in, t_Node);  // this is necessary, otherwise the correct data is
                          // not available

        // Order GEM to run
        tdBR->NodeStatusCH = NEED_GEM_AIA;

        m_NodeStatusCH[in] = t_Node->GEM_run(false);

        if ((m_NodeStatusCH[in] == ERR_GEM_AIA ||
             m_NodeStatusCH[in] == ERR_GEM_SIA))
        {
            rwmutex.lock();
            cout << "Initial GEMs run second pass failed at node " << in;
            cout << " repeat calculations and change kinetic constraintsby 0.1%"
                 << "\n";
            rwmutex.unlock();
            // change a bit the kinetic constraints -> make the system less
            // stiff
            for (j = 0; j < nDC; j++)
            {
                m_dll[in * nDC + j] =
                    0.999 * m_dll[in * nDC + j] - 1.0e-6;  // make smaller
                if (m_dll[in * nDC + j] < 0.0)
                    m_dll[in * nDC + j] = 0.0;
                m_dul[in * nDC + j] =
                    1.001 * m_dul[in * nDC + j] + 1.0e-6;  // make bigger
            }
            REACT_GEM::SetReactInfoBackGEM(in, t_Node);  // needs to be done to
            //            m_Node->GEM_write_dbr ( "dbr_for_crash_node_fail1.txt"
            //            );

            // run GEMS again
            tdBR->NodeStatusCH = NEED_GEM_AIA;
            m_NodeStatusCH[in] = t_Node->GEM_run(false);
            //            m_Node->GEM_write_dbr ( "dbr_for_crash_node_fail2.txt"
            //            );
            if ((m_NodeStatusCH[in] == ERR_GEM_AIA ||
                 m_NodeStatusCH[in] == ERR_GEM_SIA))
            {
                rwmutex.lock();
                cout << " Error: Init Loop second pass failed when running GEM "
                        "on Node #"
                     << in << "."
                     << "\n";
                cout << "Returned Error Code: " << m_NodeStatusCH[in] << "\n";
                t_Node->GEM_write_dbr("dbr_for_crash_node_init2.txt");
                rwmutex.unlock();
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                MPI_Finalize();  // make sure MPI exits
#endif
                exit(1);
            }
            else if (m_NodeStatusCH[in] == BAD_GEM_AIA ||
                     m_NodeStatusCH[in] == BAD_GEM_SIA)
            {
                rwmutex.lock();
                cout << "error: Initial GEMs run after Read GEMS gives bad "
                        "result..proceed in any case. node"
                     << in << "\n";
                rwmutex.unlock();
            }

            else
            {
                rwmutex.lock();
                cout << " sucess with second try.... "
                     << "\n";
                rwmutex.unlock();
            }
        }  // end loop if initial gems run fails

        REACT_GEM::GetReactInfoFromGEM(in, t_Node);

        //        } // end if for restart

        // calculate the chemical porosity
        if (m_flow_pcs->GetRestartFlag() == FiniteElement::WRITE ||
            m_flow_pcs->GetRestartFlag() == FiniteElement::NO_IO)
            REACT_GEM::CalcPorosity(in,
                                    t_Node);  // during init it should be always
                                              // done, except for restart !!!

        REACT_GEM::CalcReactionRate(
            in, m_T[in], t_Node);  // moved it after porosity calculation
                                   // because of chrunchflow kinetics (model 4)!
        // Convert to concentration: should be done always...

        REACT_GEM::MassToConcentration(in, 0,
                                       t_Node);  // concentrations are now based
                                                 // on the fluid-gas volumes...

#if defined(USE_MPI_GEMS)
        REACT_GEM::CopyToMPIBuffer(in);  // copy data to MPI buffer in any case
#endif
    }  // end for loop for all nodes

    gem_barrier_finish
        ->wait();  // init run finished ...now the init routine takes over
    time_gem_total += (GetTimeOfDayDouble() - tdummy);
    time_fraction = time_gem_total / GetTimeOfDayDouble();

    // main loop which is synchronized via run_main
    for (;;)
    {
        gem_barrier_start
            ->wait();  // central barrier for synchonizing the main runs
        tdummy = GetTimeOfDayDouble();

        repeated_fail = 0;  // set this to zero for each new run_main

        for (in = mystart; in < nNodes;
             in += mycount)  // myrank ist defined vi USE_MPI_GEMS
        {
            if ((!flag_calculate_boundary_nodes && m_boundary[in]) ||
                (CalcSoluteBDelta(in) <
                 m_diff_gems))  // do this only without calculation if  on a
                                // boundary or differences very small!!!
            {
#if defined(USE_MPI_GEMS)
                REACT_GEM::CopyToMPIBuffer(
                    in);  // copy old values to buffer..otherwise we loose them
#endif
            }
            else  // calculate if not on a boundary or transport causes
                  // differences in solutes
            {
                // Convert from concentration
                REACT_GEM::ConcentrationToMass(
                    in, 1);  // I believe this is save for MPI
                // now we calculate kinetic constraints for GEMS!
                REACT_GEM::CalcLimits(in, t_Node);
                // Get data

                REACT_GEM::SetReactInfoBackGEM(
                    in, t_Node);  // this should be also save for MPI
                // take values from old B volume for comparison
                oldvolume = m_Vs[in];

                // Order GEM to run
                tdBR->NodeStatusCH = NEED_GEM_AIA;  // first try without simplex
                                                    // using old solution
                m_NodeStatusCH[in] = t_Node->GEM_run(false);

                if (!(m_NodeStatusCH[in] == OK_GEM_AIA ||
                      m_NodeStatusCH[in] == OK_GEM_SIA ||
                      m_NodeStatusCH[in] == BAD_GEM_AIA ||
                      m_NodeStatusCH[in] == BAD_GEM_SIA) ||
                    (((abs(oldvolume - tdBR->Vs) / oldvolume) > 0.1) &&
                     (flowflag !=
                      3)))  // not for Richards flow  // ups...failed..try again
                            // with changed kinetics
                {
                    //#if !defined(USE_MPI_GEMS) && !defined(USE_PETSC)
                    rwmutex.lock();  // KG44 try to avoid mutual exclusion at
                                     // least in the parallel version, as this
                    // might slow down execution
                    cout << "Error: Main Loop failed when running GEM on Node #"
                         << in << "."
                         << " Returned Error Code: " << m_NodeStatusCH[in];
                    cout << " or GEM weird result at node " << in << " volume "
                         << tdBR->Vs << " old volume " << oldvolume;
                    cout << " repeat calculations and change kinetic "
                            "constraintsby 0.1%"
                         << "\n";
                    rwmutex.unlock();
                    //#endif
                    // change a bit the kinetic constraints -> make the system
                    // less stiff
                    for (j = 0; j < nDC; j++)
                    {
                        m_dll[in * nDC + j] = 0.999 * m_dll[in * nDC + j] -
                                              1.0e-6;  // make smaller
                        if (m_dll[in * nDC + j] < 0.0)
                            m_dll[in * nDC + j] = 0.0;
                        m_dul[in * nDC + j] = 1.001 * m_dul[in * nDC + j] +
                                              1.0e-6;  // make bigger
                    }
                    REACT_GEM::SetReactInfoBackGEM(
                        in, t_Node);  // needs to be done to for update dll dul

                    // run GEMS again

                    tdBR->NodeStatusCH = NEED_GEM_AIA;
                    m_NodeStatusCH[in] = t_Node->GEM_run(false);
                }

                // test for bad GEMS and for volume changes bigger than 10%
                // ...maximum 5 failed nodes per process.....

                if (!(m_NodeStatusCH[in] == OK_GEM_AIA ||
                      m_NodeStatusCH[in] == OK_GEM_SIA ||
                      m_NodeStatusCH[in] == BAD_GEM_AIA ||
                      m_NodeStatusCH[in] == BAD_GEM_SIA) ||
                    (((abs(oldvolume - tdBR->Vs) / oldvolume) > 0.1) &&
                     (flowflag != 3))  // not for Richards flow
                )
                {
                    //#if !defined(USE_MPI_GEMS) && !defined(USE_PETSC)
                    rwmutex.lock();
                    cout << "Error: Main Loop failed when running GEM on Node #"
                         << in << "."
                         << " Returned Error Code: " << m_NodeStatusCH[in];
                    cout << " or GEM weird result at node " << in << " volume "
                         << tdBR->Vs << " old volume " << oldvolume;
                    cout << " continue with last good solution for this node"
                         << "\n";
                    //                    t_Node->GEM_write_dbr (
                    //                    "dbr_for_crash_node_fail.txt" );
                    //                    t_Node->GEM_print_ipm (
                    //                    "ipm_for_crash_node_fail.txt" );
                    rwmutex.unlock();
                    //#endif
                    // exit ( 1 );
                    node_fail = 1;
                    repeated_fail += 1;
                    if (repeated_fail > m_max_failed_nodes)
                    {
                        rwmutex.lock();
                        cout << "GEFMS: " << repeated_fail
                             << "nodes failed this timestep, check chemical "
                                "system!"
                             << "\n";
                        rwmutex.unlock();
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
                        MPI_Finalize();  // make sure MPI exits
#endif
                        exit(1);
                    }  // we do not tolerate more than three failed nodes ->
                       // something wrong with chemistry/time step?
                }      // end loop if initial gems run fails
                else

                    // this is if gem run is ok
                    REACT_GEM::GetReactInfoFromGEM(
                        in,
                        t_Node);  // from here on we have buffer values if
                                  // GEMS_MPI is defined

                if (node_fail <
                    1)  // this needs to be done with buffer variables for MPI
                {
                    // CALC kintetic constrains
                    REACT_GEM::CalcReactionRate(
                        in, m_T[in],
                        t_Node);  // check for kinetics is done in the
                                  // subroutine for each species separately

                    // calculate the chemical porosity
                    if (flag_porosity_change == 1)
                        REACT_GEM::CalcPorosity(in, t_Node);

                    // Convert to concentration  ..
                    REACT_GEM::MassToConcentration(in, 0, t_Node);

                    // calculate density of fluid phase, which is normally the
                    // first phase
                    m_fluid_density[in] =
                        m_mPS[in * nPS + 0] / m_vPS[in * nPS + 0];
                }
                else
                {
                    RestoreOldSolution(in);
                    node_fail = 0;
                }
#if defined(USE_MPI_GEMS)
                REACT_GEM::CopyToMPIBuffer(
                    in);  // copy data to MPI buffer in any case
#endif
            }  // end if check for boundary node12
        }      // end for loop for all nodes

        twait = GetTimeOfDayDouble();  // this if for performance check...init
                                       // waiting time

        gem_barrier_finish->wait();  // barrier for synchonizing threads

        tdummy1 = GetTimeOfDayDouble();  // now get and calculate the times for
                                         // GEMS execution
        time_gem_total += (tdummy1 - tdummy);
        tdummy2 = (tdummy1 - tdummy);
        time_fraction = time_gem_total / (tdummy1 - tstart);

        // this is for performance check .....
        twait -= tdummy1;
        twait *= -1.0;
        twaittotal += twait;
//        rwmutex.lock();
//        cout << "GEMS3K MPI Processe / Thread: " << myrank << " " << tid << "
//        waiting time: " << twait << "\n"; rwmutex.unlock();
// write some time values for performance ...only for first thread..do only for
// all threads if
#if !defined(USE_MPI_GEMS) && !defined(USE_PETSC)
        if (tid == 0 &&
            myrank == 0)  // try to avoid mutual exclusion in MPI version
        {
            rwmutex.lock();
            cout << "GEMS3K MPI Processes / Threads: " << mysize << " "
                 << gem_nThread << " mid/tid " << myrank << " " << tid
                 << " total time: " << time_gem_total
                 << " s, this dt for GEMS: " << tdummy2
                 << " total fraction in GEMS: " << time_fraction
                 << " idle time: " << twaittotal << "\n";
            rwmutex.unlock();
        }
#endif
    }  // end of for loop....
}

double REACT_GEM::GetTimeOfDayDouble()
{
    double dtime;

#ifdef WIN32
    dtime = 1.0;  // plugin dummy value as long as sys/time.h for windoof is not
                  // available
#else
    timeval current;
    gettimeofday(&current, 0);
    dtime = 1.0 * current.tv_sec + 1.0e-6 * current.tv_usec;
#endif

    return dtime;
}

#if defined(USE_PETSC)
// The following routine is NOT thread safe, therefore do not call within
// threads!!!!!!!!!!!!!!!!!!!!!!!!!!!
void REACT_GEM::SynchronizeData(
    double* data)  // pass pointer to the vector we would like to synchronize...
{
    long i;
    CRFProcess* this_pcs;
    this_pcs =
        PCSGet("MASS_TRANSPORT");  // first mass transport process...for gems
                                   // coupling we should always have this
    // gem_glob_buff[m_size] is re-used to distribute values....
    // first fill data into global buffer...then get back data...
    // size of global buffer is m_size

    int receivecount;
    PetscInt low, high, otherlow;
    MPI_Status status;
    PetscInt count;
    int tag = 9999;

    // clean buffer, otherwise allreduce will not work
    for (i = 0; i < glob_NodesNumber_Linear; i++)
    {
        gem_glob_buff[i] = 0.0;  // put data into array
        gem_glob_x0[i] = 0.0;    // put data into array
    }
    // arrays are filled, now we should synchronize the values
    //    for ( long node_i=0; node_i < nNodes ; node_i++ )
    //    {
    //      cout << "before i, value, global_indes " << node_i << " " <<
    //      data[node_i] << " " << this_pcs->m_msh->Eqs2Global_NodeIndex[node_i]
    //      << "\n" ;
    //    }
    //    cout << "loc, nodes and glob nodes number " << loc_NodesNumber_Linear
    //    << " " << nNodes << " " << glob_NodesNumber_Linear << "\n";
    //  put data into buffer array ...only local data without shadow nodes
    for (i = 0; i < loc_NodesNumber_Linear; i++)
        gem_glob_buff[this_pcs->m_msh->Eqs2Global_NodeIndex[i]] =
            data[i];  // put data into array

    // now fill gem_glob_x0
    MPI_Allreduce(gem_glob_buff, gem_glob_x0, glob_NodesNumber_Linear,
                  MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    // put data back into data array, but now with shadow nodes
    for (i = 0; i < NodesNumber_Linear; i++)
        data[i] = gem_glob_x0[this_pcs->m_msh->Eqs2Global_NodeIndex[i]];

    //    for ( long node_i=0; node_i < nNodes ; node_i++ )
    //    {
    //      cout << "after i, value " << node_i << " " << data[node_i] << "\n";
    //    }
}

#endif

#endif
