/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   rf_kinreact.cpp

   KINETIC REACTIONS

   FEMLib-Object KinReact

   Programming:
   01/2004    Dirk Schaefer      Original Implementation
   02/2006    Sebastian Bauer    Adaption to C++ Class structure, new FEM
 concept 05/2007    Dirk Schaefer      Addition of NAPL-dissolution
 ***************************************************************************/

#include <cfloat>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "display.h"
#include "StringTools.h"
#include "files0.h"
#include "makros.h"
#include "msh_lib.h"
#include "rf_kinreact.h"
#include "mathlib.h"
#include "timer.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"
#include "tools.h"
#include "rf_react.h"
#include "rf_react_int.h"
#include "PhysicalConstant.h"

#ifdef OGS_FEM_CAP  // CAP_REACT  //CB merge CAP 0311
#include "rf_react_cap.h"
#endif
#if defined(USE_MPI) || defined(USE_MPI_KRC)
//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END
#include "mpi.h"  //Parallel Computing Support
#include "par_ddc.h"
#include "SplitMPI_Communicator.h"
#endif

//#include "msh_mesh.h"
using namespace std;
using SolidProp::CSolidProperties;
using namespace Math_Group;

const double residual = 1.e-20;

vector<CKinReact*> KinReact_vector;  // declare instance CKinReact_vector
vector<CKinReactData*>
    KinReactData_vector;           // declare instance CKinReact_vector
vector<CKinBlob*> KinBlob_vector;  // declare extern instance of class Blob
// CB _drmc_
vector<MicrobeData*>
    MicrobeData_vector;  // declare extern instance of class Blob

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b)                   \
    (dmaxarg1 = (a), dmaxarg2 = (b), \
     (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define DMIN(a, b)                   \
    (dmaxarg1 = (a), dmaxarg2 = (b), \
     (dmaxarg1) < (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

/* Constructor for MonodSubstruct */
MonodSubstruct::MonodSubstruct(void)
{
    species = "blubb";
    speciesnumber = -1;
    concentration = -1.0e9;
    order = -99.0;
    isotopecouplenumber = -1;  // CB isotope fractionation
    threshhold = false;
    threshConc = -1.0e9;
    threshOrder = -99;
}

/* Destructor for MonodSubstruct */
MonodSubstruct::~MonodSubstruct(void) {}

// CB _drmc_
MicrobeData::MicrobeData(void)
{
    MicrobeName = "Blupp";
    MonodReactionName = "Bluxx";
    steepness = -1.0;
    k_incr = 0.0;
    k_decr = 0.0;
    G0 = 0;
    dt = 0;
    decayrate = 0.0;
    _drmc_ = false;
    MonodReactionIdx = -1;

    Gibbs.clear();
    _drmc_level.clear();
}

MicrobeData::~MicrobeData(void) {}

double MicrobeData::GetGibbsEnergy(long node)
{
    double R = 8.314472;   // Gas constant
    double T = 298.25;     // K
    double G = 0, dG = 0;  // dummy
    int MonodSpecies = 0;
    double MonodConcentration = 1.0;
    int MonodOrder = 0;
    double Ctot, C = 0;
    // int Isotopespecies = 0;
    double Monodterm = 1.0;
    // double Inhibition = 1.0;
    // double Rateconstant = 1.0;
    // int NumberMonod = 0;
    double scoeff = 0;
    double unitfactor_l = 0.001;

    CKinReact* m_kr = NULL;
    m_kr = KinReact_vector[MonodReactionIdx];
    CKinReactData* m_krd = NULL;
    m_krd = KinReactData_vector[0];

    if (REACTINT_vec.size() > 0)  // Get the Temperature
        T = REACTINT_vec[0]->GetTemperature(node);

    // G = dG * -µmax*Prod[Ai/(Ai+KAi)]
    //   = [dG0 + R*T*ln(Q)] * -µmax*Prod[Ai/(Ai+KAi)]

    // calculate dG = dG0 + R*T*ln(Q)
    for (int i = 0; i < m_kr->number_reactionpartner; i++)
    {
        MonodSpecies = m_kr->ProdStochhelp[i]->speciesnumber;
        scoeff = m_kr->ProductionStoch[MonodSpecies];
        C = cp_vec[MonodSpecies]->getProcess()->GetNodeValue(
            node, m_krd->sp_varind[MonodSpecies]);
        if (C <= 0)
            C = 1e-10;
        dG += scoeff * log(1 * C * unitfactor_l);
    }
    dG *= R * T;
    dG += m_kr->dG0;

    // dS/dt/B term = -µmax*Prod[Ai/(Ai+KAi)]
    for (int i = 0; i < m_kr->number_monod; i++)
    {
        MonodSpecies = m_kr->monod[i]->speciesnumber;
        MonodConcentration = m_kr->monod[i]->concentration;
        MonodOrder = (int)m_kr->monod[i]->order;  // CB higher order Monod terms
        C = cp_vec[MonodSpecies]->getProcess()->GetNodeValue(
            node, m_krd->sp_varind[MonodSpecies]);
        Ctot = C;  // standard case, no iso fractionation
        // CB Isotope fractionation: here Ctot=C_light+C_heavy must be passed to
        // fct Monod()
        if ((m_kr->typeflag_iso_fract == 1) &&
            (m_kr->monod[i]->isotopecouplenumber >= 0))
        {
            MonodSpecies = m_kr->monod[i]->isotopecouplenumber;
            Ctot = C + cp_vec[MonodSpecies]->getProcess()->GetNodeValue(
                           node, m_krd->sp_varind[MonodSpecies]);
        }
        Monodterm *= m_kr->Monod(MonodConcentration, C, Ctot,
                                 MonodOrder);  // new formulation
        // now multiply by additional Threshhold Term, if present,
        // technically this is the same as an additional Monod term for the same
        // MonodSpecies usually of higher order and lower Monod-Conc =
        // threshConc
        if (m_kr->monod[i]->threshhold == true)
        {
            MonodConcentration =
                m_kr->monod[i]->threshConc;  // CB Threshhold concentration
            MonodOrder = (int)m_kr->monod[i]
                             ->threshOrder;  // CB higher order Monod terms
            Monodterm *=
                m_kr->Monod(MonodConcentration, Ctot, Ctot,
                            MonodOrder);  // C should in any case be total C
        }
    }
    for (int i = 0; i < m_kr->number_inhibit; i++)
    {
        MonodSpecies = m_kr->inhibit[i]->speciesnumber;
        MonodConcentration = m_kr->inhibit[i]->concentration;
        C = cp_vec[MonodSpecies]->getProcess()->GetNodeValue(
            node, m_krd->sp_varind[MonodSpecies]);
        Monodterm *= m_kr->Inhibition(MonodConcentration, C);
    }

    // put together  G = dG * ds/dt/X ;
    // division by X not necessary here, as we calculated the rate explicitly
    G = -m_kr->rateconstant * Monodterm * dG;
    // dG = 200000; // dummy

    return G;
}

/**************************************************************************
 FEMLib-Method:
 Task: OBJ configure function
 Programing:
 05/2012 CB Implementation
 **************************************************************************/
void MicrobeConfig(void)
{
    // create vector for _drmc_ level and switsh function and initialization
    // with value
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    if (m_msh == NULL)
    {
        cout << "No mesh in KRConfig"
             << "\n";
        exit(1);
    }

    for (size_t i = 0; i < MicrobeData_vector.size(); i++)
    {
        MicrobeData* m_md(MicrobeData_vector[i]);
        for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
        {
            m_md->Gibbs.push_back(1.0);        // Vorbelegung mit theta=-1
            m_md->_drmc_level.push_back(1.0);  // Vorbelegung mit S=-1
        }
    }
}

/**************************************************************************
 FEMLib-Method:
 Task: OBJ read function for CKinBlob-Structure
 Programing:
 02/2007 DS Implementation
 **************************************************************************/
bool MicrobeData::Read(ifstream* rfd_file)
{
    char line[MAX_ZEILE];
    string line_string, line_str1, s_geo_type, s_geo_name;
    string hash("#"), dollar("$");
    string species;

    bool new_keyword = false;  //, OK = true, new_subkeyword = false;
    long index;
    std::stringstream in;

    //========================================================================
    while (!new_keyword)
    {
        index = rfd_file->tellg();
        //    if(!rfd_file->getline(line,MAX_ZEILE)) break;
        if (!GetLineFromFile(line, rfd_file))
            break;
        line_string = line;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
                                     // n�chste keyword weg
            break;
        }
        /* Keywords nacheinander durchsuchen */
        //....................................................................
        // subkeyword found
        if (line_string.find("$MICROBENAME") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> MicrobeName;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$_drmc__PARAMETERS") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> steepness >> k_incr >> k_decr >> G0 >> decayrate;
            // check for read in
            if ((steepness != -1.0) && (k_incr != 0.0) && (k_decr != 0.0))
                _drmc_ = true;
            else
            {
                Display::DisplayMsgLn(" ERROR reading Microbe _drmc_ Terms  - skipping");
                return false;
            }
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$MONOD_REACTION_NAME") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> MonodReactionName;
            in.clear();
        }

    }  // end while keyword

    return true;
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinReact constructor
   Programing:
   02/2006 SB Implementation
   06/2007 DS NAPL-dissolution added
***************************************************************************/
CKinReact::CKinReact(void)
{
    name = "NULL";
    type = "NULL";
    number_reactionpartner = 0;
    reactionpartner.clear();
    stochmet.clear();
    rateconstant = 0.0;
    rateorder = 1.0;  // CB default changed from 0 to 1, 06/2012
    number_monod = 0;
    number_inhibit = 0;
    number_production = 0;
    number_isotope_couples = 0;  // CB isotope fractionation
    monod.clear();
    inhibit.clear();
    production.clear();
    bacteria_name = "NULL";
    bacteria_number = -1;
    ProductionStoch.clear();
    grow = -1;
    specif_cap = -1.0;
    isoenfac = 0;      // CB Isotope fractionation
    degType = "NULL";  // CB Isotope fractionation
    T_dependence = false;
    T_model = 0;
    T_params.clear();

    // CB _drmc_
    _drmc_ = false;
    MicrobeData_idx = -1;
    DormType = "NULL";
    DormTypeIdx = -1;
    dG0 = 0;
    Yieldcoefficient = 1;  // microbial yield coeff, is 1 as standard case

    // CB Not this particular reaction on specified GEO-Objects
    switched_off_node.clear();

    ProdStochhelp.clear();
    ex_param.clear();
    ex_species.clear();
    ex_species_names.clear();
    exSurfaceID = -1;
    exType = "NULL";
    // NAPL-dissolution
    blob_name = "NULL";
    blob_ID = -1;
    Csat_pure = 0.;
    current_Csat = 0.;
    Density_NAPL = 0.;
    Current_Csat.clear();
    //
    typeflag_monod = 0;
    typeflag_exchange = 0;
    typeflag_exchange_linear = 0;
    typeflag_exchange_langmuir = 0;
    typeflag_exchange_freundlich = 0;
    typeflag_napldissolution = 0;
    typeflag_iso_fract = 0;        // CB isotope fractionation
    typeflag_mineralkinetics = 0;  // CB
    typeflag_gasdissolution = 0;
    // Minkin
    mechvec.clear();
    number_Mech = 0;
    minSpeciesIdx.clear();
    Am.clear();  // initial specific reactive mineral surface area
    Km.clear();  // Equilibrium constant (at node) vector (if Km_uniform ==
                 // false)
    Cminini.clear();
    Km_uniform = true;
    Km_CHEMAPP = false;
    Km_HKF = false;
    Am_constant = true;
    Am_model = 0;
    Am_ini = 0;
    Eta = 1.0;
    Theta = 1.0;
    precip_baseterm_only = false;
    precipfactor = 1.0;
    precipexponent = 1.0;
    scale_rate = false;
    lagneau = false;
    OmegaThreshhold =
        1e-10;  // if |1-Omega| < OmegaThreshhold, set rate to zero
    MinKinRateCoeff = 0;
    ExplicitMinKinRateCoeff = false;
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinReact destructor
   Programing:
   02/2006 SB Implementation
***************************************************************************/
CKinReact::~CKinReact(void) {}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   02/2004 SB Implementation
**************************************************************************/
bool KRRead(const std::string& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
    // File handling
    std::string krc_file_name(file_base_name + KRC_FILE_EXTENSION);
    std::ifstream krc_file(krc_file_name.data(), ios::in);
    if (!krc_file.good())
        return false;

    CKinReact *m_kr = NULL, *m_kr1 = NULL;
    CKinReactData* m_krd = NULL;
    CKinBlob* m_bp = NULL;
    MicrobeData* m_md = NULL;
    char line[MAX_ZEILE];
    string sub_line;
    string line_string;
    string database_file_name;
    int found;
    size_t length;
    ios::pos_type position;
    string m_name, sp_name;

    KRCDelete();
    //========================================================================
    // Keyword loop
    cout << "KinReact Read"
         << "\n";
    while (!krc_file.eof())
    {
        krc_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != string::npos)
            break;
        //----------------------------------------------------------------------
        // keyword found  // Read kinetic Reactions
        if (line_string.find("#REACTION") != string::npos)
        {
            m_kr = new CKinReact();
            m_kr->Read(&krc_file, geo_obj, unique_name);
            m_kr->number = (int)KinReact_vector.size();
            KinReact_vector.push_back(m_kr);
        }  // keyword found
        // keyword found  // Read KinReactData
        if (line_string.find("#KINREACTIONDATA") != string::npos)
        {
            m_krd = new CKinReactData();
            m_krd->Read(&krc_file, geo_obj, unique_name);
            KinReactData_vector.push_back(m_krd);
        }  // keyword found
        // keyword found  // Read BlobProperties
        if (line_string.find("#BLOB_PROPERTIES") != string::npos)
        {
            m_bp = new CKinBlob();
            m_bp->Read(&krc_file, geo_obj, unique_name);
            KinBlob_vector.push_back(m_bp);
        }  // keyword found
        // CB _drmc_
        if (line_string.find("#MICROBE_PROPERTIES") != string::npos)
        {
            m_md = new MicrobeData();
            if (m_md->Read(&krc_file))
                MicrobeData_vector.push_back(m_md);
        }  // keyword found

    }  // eof
    // Close input file
    krc_file.close();

    if (!m_krd)
    {
        cout << " No keyword #KINREACTIONDATA specified - setting defaults"
             << "\n";
        m_krd = new CKinReactData();
        KinReactData_vector.push_back(m_krd);
    }

    //========================================================================

    /* Check, if data has to be completed from database file */
    //  cout << "checking database" << "\n";
    if (!KinReact_vector.empty())
    {
        // File handling, open reaction database file
        database_file_name = "reactions.dbf";
        ifstream dbf_file(database_file_name.data(), ios::in);
        //	if (!dbf_file.good()) cout << " No database file found " << "\n";

        // go through all reactions in reaction vector
        length = KinReact_vector.size();
        for (size_t i = 0; i < length; i++)
        {
            m_kr = KinReact_vector[i];
            // if no type is given in input file, only name, then read from
            // database file
            if (m_kr->type == "NULL")
            {
                found = 0;
                dbf_file.seekg(0L, ios::beg);  // rewind database file
                while (!dbf_file.eof() && found == 0)
                {
                    if (!GetLineFromFile(line, &dbf_file))
                        break;
                    line_string = line;
                    // keyword found
                    if (line_string.find("#REACTION") != string::npos)
                    {
                        m_kr1 = new CKinReact();
                        position = m_kr1->Read(&dbf_file, geo_obj, unique_name);
                        // check if this is the right one
                        if (m_kr->name == m_kr1->name)
                        {
                            // Insert in Reaction vector and remove old reaction
                            // (replacement)
                            KinReact_vector.insert(
                                KinReact_vector.begin() + i + 1, m_kr1);
                            KinReact_vector.erase(KinReact_vector.begin() + i);
                            found = 1;
                        }
                        else
                            // not the right equation:
                            delete m_kr1;
                    }  // end if(line_str..)  keyword found
                }      // eof dbf-file

                if (found == 0)
                {
                    // reaction not complete in input file and also not
                    // specified in dbf_file
                    cout << " ERROR! Reaction " << m_kr->name
                         << "  not found in database file"
                         << "\n";
                    cout << " Reaction removed from set of reaction "
                         << "\n";
                    // remove reaction from reaction vector
                    KinReact_vector.erase(KinReact_vector.begin() + i);
                }
            }  // end of if(m_kr=NULL
        }      // end for
        // close database file
        dbf_file.close();
    }  // end if kinreaction_vec.size() >
    /* end data base file read for reactions  */

    //========================================================================

    /* check reaction data consistency */
    std::cout << " Checking reaction data consistency for "
              << KinReact_vector.size() << " specified reactions "
              << "\n";
    length = KinReact_vector.size();
    for (size_t i = 0; i < length; i++)
    {
        m_kr = KinReact_vector[i];
        if (!m_kr->CheckReactionDataConsistency())
        {
            cout << " ERROR! Reaction " << m_kr->name
                 << "  removed from set of reactions due to data inconsistency"
                 << "\n";
            KinReact_vector.erase(KinReact_vector.begin() + i);
        }
    }  // end for(i=0;.. consistency check
    return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ configure function
   Programing:
   02/2004 SB Implementation
   05/2007 DS NAPL dissolution added
   07/2009 CB Isotope fractionation
   08/2009 CB Reaction deactivation
**************************************************************************/
void KRConfig(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
    CKinReact* m_kr = NULL;  //, *m_kr1=NULL;
    CKinReactData* m_krd = NULL;
    MinkinMech* m_mech = NULL;
    MicrobeData* m_md = NULL;
    int k;  // i, length, j, ;
    // int idummy; //idx,
    long ll, lm;
    string m_name, sp_name;
    CompProperties* m_cp = NULL;
    CRFProcess* m_pcs = NULL;
    vector<long> nodes_vector;
    CMediumProperties* m_mat_mp = NULL;
    // long group;
    double ww;  // foc, , w;

    string dummy;
    bool ok = true;
    std::vector<double> helpvec;

    // CB reaction deactivation
    // int annode_idx, lll, llll, lllll, nnodpneigh, duplicate,nn, nnodsneigh
    // int dnele_idx; // ,  pneighnod_idx, snele_idx,
    //  int dnnode, snnode, duplicate2;
    MeshLib::CElem* m_dnele = NULL;
    MeshLib::CElem* m_snele = NULL;
    MeshLib::CNode* m_dnnod = NULL;
    std::vector<int> ReactNeighborNodes;
    vec<long> secnnodesindices(8);
    vec<long> primnnodesindices(8);

    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    if (m_msh == NULL)
    {
        cout << "No mesh in KRConfig"
             << "\n";
        exit(0);
    }
    //========================================================================
    if (!KinReactData_vector.empty())
    {
        m_krd = KinReactData_vector[0];
        if (m_krd == NULL)
            cout << " Error - no m_krd data "
                 << "\n";
        // Set up additional data structures for calculations

        // Set number of reactions
        m_krd->NumberReactions = (int)KinReact_vector.size();
        size_t length(cp_vec.size());

        // Check if all reaction partners are specified as processes
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            // Check presence of bacteria species
            if (m_kr->bacteria_name.compare("NULL") != 0)
            {
                m_pcs = PCSGet("MASS_TRANSPORT", m_kr->bacteria_name);
                if (m_pcs == NULL)
                {
                    cout << " Warning: Component " << m_kr->bacteria_name
                         << " specified in KinReact as biomass but not given "
                            "as transport process "
                         << "\n";
                    ok = false;
                }
            }
            // Check if monod species are specified as components
            for (int i = 0; i < m_kr->number_monod; i++)
            {
                m_pcs = PCSGet("MASS_TRANSPORT", m_kr->monod[i]->species);
                if (m_pcs == NULL)
                {
                    cout << " Warning: Component " << m_kr->monod[i]->species
                         << " specified in KinReact as monod species but not "
                            "given as transport process "
                         << "\n";
                    ok = false;
                }
            }
            // Check if inhibition species are specified as components
            for (int i = 0; i < m_kr->number_inhibit; i++)
            {
                m_pcs = PCSGet("MASS_TRANSPORT", m_kr->inhibit[i]->species);
                if (m_pcs == NULL)
                {
                    cout << " Warning: Component " << m_kr->inhibit[i]->species
                         << " specified in KinReact as inhibition species but "
                            "not given as transport process "
                         << "\n";
                    ok = false;
                }
            }
            // Check productionstoch
            for (size_t i = 0; i < m_kr->ProdStochhelp.size(); i++)
            {
                m_pcs =
                    PCSGet("MASS_TRANSPORT", m_kr->ProdStochhelp[i]->species);
                if (m_pcs == NULL)
                {
                    cout << " Warning: Component "
                         << m_kr->ProdStochhelp[i]->species
                         << " specified in KinReact as produced species but "
                            "not given as transport process "
                         << "\n";
                    ok = false;
                }
            }
            if (m_kr->type.compare("exchange") == 0)
                for (int i = 0; i < m_kr->number_reactionpartner; i++)
                {
                    m_pcs = PCSGet("MASS_TRANSPORT", m_kr->reactionpartner[i]);
                    if (m_pcs == NULL)
                    {
                        cout << " Warning: Component "
                             << m_kr->reactionpartner[i]
                             << " specified in KinReact as reaction partner "
                                "but not given as transport process "
                             << "\n";
                        ok = false;
                    }
                }
            // check isotope couples
            // todo CB isotope fract
            for (int i = 0; i < m_kr->number_isotope_couples; i++)
            {
                m_pcs = PCSGet("MASS_TRANSPORT", m_kr->Isotope_light);
                if (m_pcs == NULL)
                {
                    cout << " Warning: Component " << m_kr->Isotope_light
                         << " specified in KinReact in isotope couple but not "
                            "given as transport process "
                         << "\n";
                    ok = false;
                }
            }
            // check mechanism terms
            for (int i = 0; i < (int)m_kr->number_Mech; i++)
            {  // todo CB isotope fract
                m_mech = m_kr->mechvec[i];
                for (k = 0; k < m_mech->no_mechSpec; k++)
                {  // todo CB isotope fract
                    m_pcs =
                        PCSGet("MASS_TRANSPORT", m_mech->mechSpeciesNames[k]);
                    if (m_pcs == NULL)
                    {
                        cout << " Warning: Component "
                             << m_mech->mechSpeciesNames[k]
                             << " specified in KinReact in Mechanism term but "
                                "not given as transport process "
                             << "\n";
                        ok = false;
                    }
                }
            }

        }  // loop over m_krd->NumberReactions
        //
        if (ok == false)
        {
            cout << " Components missing, Stopping"
                 << "\n";
            cout.flush();
            exit(1);
        }

        // Set vector is_a_bacterium
        //   cout << " Length of cp_vec: " << length << "\n";
        for (size_t j = 0; j < length; j++)
            m_krd->is_a_bacterium.push_back(0);  // initialize
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("monod") == 0)
                for (size_t i = 0; i < length; i++)
                {
                    m_cp = cp_vec[i];
                    if (m_kr->bacteria_name.compare(m_cp->compname) == 0)
                        m_krd->is_a_bacterium[i] = 1;
                }
        }
        // Set Vector ProductionStoch for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if ((m_kr->type.compare("monod") == 0) ||
                (m_kr->type.compare("Mineralkinetics") == 0))
            {
                for (size_t i = 0; i < length; i++)
                    // initialize
                    m_kr->ProductionStoch.push_back(0.0);

                // Get Stochiometry
                for (size_t k = 0; k < m_kr->ProdStochhelp.size(); k++)
                {
                    std::string const& m_name(m_kr->ProdStochhelp[k]->species);
                    for (size_t i = 0; i < length; i++)
                        if (m_name.compare(cp_vec[i]->compname) == 0)
                        {
                            m_kr->ProductionStoch[i] =
                                m_kr->ProdStochhelp[k]->concentration;
                            m_kr->ProdStochhelp[k]->speciesnumber = i;
                        }
                }
            }  // if type == monod or type == Mineralkinetics
        }      // vector ProductionStoch

        // Set numbers for monod species for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("monod") == 0)
                for (size_t i = 0; i < m_kr->monod.size(); i++)
                {
                    std::string const& m_name(m_kr->monod[i]->species);
                    for (size_t k = 0; k < length; k++)
                        if (m_name.compare(cp_vec[k]->compname) == 0)
                            m_kr->monod[i]->speciesnumber = k;
                }
        }  // monod substructure numbers

        // CB Isotope fractionation
        // Set number for isotope couple in monod reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->degType.compare("isotope_fractionation") == 0)
                for (size_t i = 0; i < m_kr->monod.size(); i++)
                {
                    std::string m_name(m_kr->monod[i]->species);
                    if (m_name.compare(m_kr->Isotope_light) == 0)
                    {
                        m_name = m_kr->Isotope_heavy;
                        for (size_t k = 0; k < length; k++)
                            if (m_name.compare(cp_vec[k]->compname) == 0)
                            {
                                m_kr->monod[i]->isotopecouplenumber = k;
                                break;
                            }
                    }
                    else if (m_name.compare(m_kr->Isotope_heavy) == 0)
                    {
                        m_name = m_kr->Isotope_light;
                        for (size_t k = 0; k < length; k++)
                            if (m_name.compare(cp_vec[k]->compname) == 0)
                            {
                                m_kr->monod[i]->isotopecouplenumber = k;
                                break;
                            }
                    }
                }
        }  // monod isotope substructure numbers

        // Set numbers for inhibition species for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("monod") == 0)
                for (size_t i = 0; i < m_kr->inhibit.size(); i++)
                {
                    std::string const& m_name(m_kr->inhibit[i]->species);
                    for (size_t k = 0; k < length; k++)
                        if (m_name.compare(cp_vec[k]->compname) == 0)
                            m_kr->inhibit[i]->speciesnumber = k;
                }
        }  // inhibition substructure numbers

        // Set bacteria number for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("monod") == 0)
            {
                std::string const& m_name(m_kr->bacteria_name);
                for (size_t k = 0; k < length; k++)
                {
                    m_cp = cp_vec[k];
                    if (m_name.compare(m_cp->compname) == 0)
                        m_kr->bacteria_number = k;
                }
            }
        }  // bacteria numbers

        // Set MicrobeData for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->_drmc_)
            {
                // Set MicrobeData group index
                m_name = m_kr->Microbename;
                for (size_t k = 0; k < MicrobeData_vector.size(); k++)
                {
                    if (m_name.compare(MicrobeData_vector[k]->MicrobeName) == 0)
                        m_kr->MicrobeData_idx = k;
                }
                if (m_kr->MicrobeData_idx == -1)
                {
                    cout << " Microbe name " << m_name
                         << " does not match MicrobeData! Stopping"
                         << "\n";
                    cout.flush();
                    exit(1);
                }
                // set _drmc_ type index
                if (m_kr->DormType.compare("GROWTH") == 0)
                    m_kr->DormTypeIdx = 0;
                else if (m_kr->DormType.compare("DEACTIVATION") == 0)
                    m_kr->DormTypeIdx = 1;
                else if (m_kr->DormType.compare("REACTIVATION") == 0)
                    m_kr->DormTypeIdx = 2;
                else if (m_kr->DormType.compare("DECAY") == 0)
                    m_kr->DormTypeIdx = 3;
                else
                {
                    cout << " Unknown _drmc_ Reaction Type " << m_kr->DormType
                         << "! Stopping"
                         << "\n";
                    cout.flush();
                    exit(1);
                }
                //
            }
        }  // MicrobeData

        // Set flags type_monod and type_exchange
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            //			dummy = m_kr->type;
            if (m_kr->type.compare("monod") == 0)  // CB Isotope fractionation
            {
                m_kr->typeflag_monod = 1;
                if (m_kr->degType.compare("isotope_fractionation") == 0)
                    m_kr->typeflag_iso_fract = 1;
            }
            if (m_kr->type.compare("exchange") == 0)
            {
                m_kr->typeflag_exchange = 1;
                if (m_kr->exType.compare("linear") == 0)
                    m_kr->typeflag_exchange_linear = 1;
                if (m_kr->exType.compare("langmuir") == 0)
                    m_kr->typeflag_exchange_langmuir = 1;
                if (m_kr->exType.compare("freundlich") == 0)
                    m_kr->typeflag_exchange_freundlich = 1;
            }
            if (m_kr->type.compare("NAPLdissolution") == 0)
                m_kr->typeflag_napldissolution = 1;
            if (m_kr->type.compare("Mineralkinetics") == 0)
                m_kr->typeflag_mineralkinetics = 1;
        }  // typeflags

        // set index number for MicrobeData
        for (int j = 0; j < (int)MicrobeData_vector.size(); j++)
        {
            m_md = MicrobeData_vector[j];
            for (int k = 0; k < m_krd->NumberReactions; k++)
            {
                m_kr = KinReact_vector[k];
                if (m_md->MonodReactionName.compare(m_kr->name) == 0)
                {
                    m_md->MonodReactionIdx = k;
                    break;
                }
            }
            if (m_md->MonodReactionIdx == -1)
            {  // failure
                cout << " Unknown Monod Reaction for _drmc_ model: "
                     << m_md->MonodReactionName << " ! Stopping"
                     << "\n";
                cout.flush();
                exit(1);
            }
        }

        // exchange reactions
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("exchange") == 0)
                // move names from equation to vector with names
                for (size_t i = 0; i < m_kr->reactionpartner.size(); i++)
                {
                    std::string const& m_name(m_kr->reactionpartner[i]);
                    m_kr->ex_species_names.push_back(m_name);
                    //	find species numbers for soecies names
                    for (size_t k = 0; k < length; k++)
                        if (m_name.compare(cp_vec[k]->compname) == 0)
                            m_kr->ex_species.push_back(k);
                }
        }  // exchange species numbers

        //#ds NAPLdissolution reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("NAPLdissolution") == 0)
            {
                // move names from equation to vector with names
                for (size_t i = 0; i < m_kr->reactionpartner.size(); i++)
                {
                    std::string const& m_name(m_kr->reactionpartner[i]);
                    m_kr->ex_species_names.push_back(m_name);
                    //	find species numbers for species names
                    for (size_t k = 0; k < length; k++)
                        if (m_name.compare(cp_vec[k]->compname) == 0)
                            m_kr->ex_species.push_back(k);
                }
                // Prepare vector Current_Csat for each node
                for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
                    m_kr->Current_Csat.push_back(0);
            }

            // set index numbers for blob_id
            CKinBlob* m_kb = NULL;
            for (size_t i = 0; i < KinBlob_vector.size(); i++)
            {
                m_kb = KinBlob_vector[i];
                if (m_kr->blob_name.compare(m_kb->name) == 0)
                {
                    m_kr->blob_ID = i;
                    if (m_kb->gas_dissolution_flag)
                        m_kr->typeflag_gasdissolution = 1;
                }
            }  // blob vector size
        }      // NAPLdissolution species numbers

        // CB Mineralkinetics reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("Mineralkinetics") == 0)
            {
                if (m_krd->activity_model > 0)
                {
                    // initialize help vector for activties, do only once
                    if (helpvec.size() == 0)
                        for (size_t i = 0; i < length;
                             i++)  // length = no of species
                            helpvec.push_back(0.0);
                    // Activity coefficients, do only once
                    if (m_krd->ActivityCoefficients.size() == 0)
                        for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
                            m_krd->ActivityCoefficients.push_back(helpvec);
                    // Ionic Strengths, do only once
                    if (m_krd->IonicStrengths.size() == 0)
                        for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
                            m_krd->IonicStrengths.push_back(0);
                }
                // Prepare Km vector, for each node 1 value, per reaction
                if (m_kr->Km_uniform == false)
                {
                    ww = m_kr->Km[0];
                    for (size_t l = 1; l < m_msh->nod_vector.size();
                         l++)  // start at idx 1, as first value has been set
                               // already in CKinReact::READ
                        m_kr->Km.push_back(ww);
                }
                // Prepare Am vector, for each node 1 value, per reaction
                if (m_kr->Am_constant == false)
                {
                    ww = m_kr->Am[0];
                    for (size_t l = 1; l < m_msh->nod_vector.size();
                         l++)  // start at idx 1, as first value has been set
                               // already in CKinReact::READ
                        m_kr->Am.push_back(ww);
                }
                // Set Omega treshholds from MinReactData
                m_kr->OmegaThreshhold = m_krd->OmegaThresh;
            }
        }  // Mineralkinetics

        // Set numbers for Mineral species for each reaction
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("Mineralkinetics") == 0)
            {
                // Mineral index
                m_name = m_kr->mineral_name;
                for (size_t k = 0; k < length; k++)
                {
                    if (m_name.compare(cp_vec[k]->compname) == 0)
                    {
                        m_kr->mineral_number = k;
                        if (cp_vec[k]->transport_phase != 1)
                            cout << "Warning: Mineral is not defined as solid "
                                    "phase species, phase # "
                                 << cp_vec[k]->transport_phase << "\n";
                        break;
                    }
                }
                // water index (this needs to be defined for reactions involving
                // H2O as a species
                for (size_t k = 0; k < length; k++)
                {
                    if ((cp_vec[k]->compname.compare("H2O") == 0) ||
                        (cp_vec[k]->compname.compare("H2O_liquid") == 0) ||
                        (cp_vec[k]->compname.compare("water_liquid") == 0))
                    {
                        m_kr->water_number = k;
                        break;
                    }
                }
                // Mineral species
                for (int i = 0; i < (int)m_kr->number_reactionpartner; i++)
                {
                    m_name = m_kr->reactionpartner[i];
                    for (size_t k = 0; k < length; k++)
                        if (m_name.compare(cp_vec[k]->compname) == 0)
                            m_kr->minSpeciesIdx.push_back(k);
                }
                // Mechanism term species
                for (size_t i = 0; i < m_kr->mechvec.size(); i++)
                {
                    m_mech = m_kr->mechvec[i];
                    for (int l = 0; l < m_kr->mechvec[i]->no_mechSpec; l++)
                    {
                        m_name = m_kr->mechvec[i]->mechSpeciesNames[l];
                        for (size_t k = 0; k < length; k++)
                            if (m_name.compare(cp_vec[k]->compname) == 0)
                                m_kr->mechvec[i]->mechSpeciesIdx.push_back(k);
                    }
                }
            }
        }  // Mineral Kinetics & substructure numbers

        // check for necessary data in Mineral kinetics
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->type.compare("Mineralkinetics") == 0)
            {
                // surface area data
                if (m_kr->Am_constant == false)
                {
                    if (m_kr->Am_model > 1)
                    {
                        if (cp_vec[m_kr->mineral_number]->molar_weight == 0)
                        {
                            cout
                                << " Warning for surface area model of Mineral "
                                << m_kr->mineral_name << ":"
                                << "\n";
                            cout << "   No molar_weight defined in component "
                                    "properties."
                                 << "\n";
                            ok = false;
                        }
                        if (cp_vec[m_kr->mineral_number]->mineral_density == 0)
                        {
                            cout
                                << " Warning for surface area model of Mineral "
                                << m_kr->mineral_name << ":"
                                << "\n";
                            cout << "   No mineral_density defined in "
                                    "component properties."
                                 << "\n";
                            ok = false;
                        }
                    }
                }
                // node porosities
                // water concentrations
                // etc
            }
        }
        // set special case of Lagneau-Benchmark
        if (m_krd->lagneau)
        {
            for (int j = 0; j < m_krd->NumberReactions; j++)
            {
                KinReact_vector[j]->lagneau = true;
            }
        }
        //....................................................................
        if (m_krd->scale_dcdt)
        {  // scale dcdt vector in derivs for stabilization of ODE solver
            for (int j = 0; j < m_krd->NumberReactions; j++)
                KinReact_vector[j]->scale_rate = true;
        }
        // exchange reactions numbers for m_krd
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if (m_kr->exType.compare("linear") == 0)
                m_krd->NumberLinear++;
            if (m_kr->exType.compare("langmuir") == 0)
                m_krd->NumberLangmuir++;
            if (m_kr->exType.compare("freundlich") == 0)
                m_krd->NumberFreundlich++;
            if (m_kr->type.compare("monod") == 0)
                m_krd->NumberMonod++;
            if (m_kr->type.compare("NAPLdissolution") == 0)
                m_krd->NumberNAPLdissolution++;
            if (m_kr->type.compare("Mineralkinetics") == 0)
                m_krd->NumberMineralkinetics++;
        }
        m_krd->NumberMicrobeData = (int)MicrobeData_vector.size();

        // set up vectors sp_pcs and sp_varind
        for (size_t j = 0; j < length; j++)
        {
            std::string const& sp_name(cp_vec[j]->compname);
            // HS, PCSGet not needed any more.
            // m_pcs = PCSGet("MASS_TRANSPORT", sp_name);
            // new style:
            m_pcs = cp_vec[cp_name_2_idx[sp_name]]->getProcess();
            if (m_pcs != NULL)
                // HS, not needed any more
                // idummy = PCSGetPCSIndex("MASS_TRANSPORT", sp_name);
                // m_krd->sp_pcsind.push_back(idummy);
                // new timelevel
                m_krd->sp_varind.push_back(m_pcs->GetNodeValueIndex(sp_name) +
                                           1);
            //	cout << " PCS: " << j << ", " << sp_name << ", " << idummy << ",
            //" << m_pcs->GetNodeValueIndex(sp_name)
            //+ 1 << "\n";
        }

        // CB isotope fractionation
        // modify rateconstant for Monod-type reactions and isotope
        // fractionation
        for (int j = 0; j < m_krd->NumberReactions; j++)
        {
            m_kr = KinReact_vector[j];
            if ((m_kr->typeflag_monod == 1) && (m_kr->typeflag_iso_fract == 1))
                m_kr->rateconstant *= (1 + m_kr->isoenfac / 1000);
        }
        /********************************************************/
        // check global requirements for kinetic reactions:
        // check if porosities for phases are set
        if (m_krd->NumberMonod > 0)
            for (size_t j = 0; j < mmp_vector.size(); j++)
                if (mmp_vector[j]->vol_bio < MKleinsteZahl)
                    std::cout << "Warning: Porosity of bio-phase is 0.0 ! "
                                 "Change Settings in *.mmp file "
                              << "\n";
        //#ds  k= m_krd->NumberLinear + m_krd->NumberFreundlich +
        // m_krd->NumberLangmuir;
        int k = m_krd->NumberLinear + m_krd->NumberFreundlich +
                m_krd->NumberLangmuir + m_krd->NumberNAPLdissolution;
        if (k > 0)
            for (size_t j = 0; j < mmp_vector.size(); j++)
                if (mmp_vector[j]->vol_mat < MKleinsteZahl)
                    std::cout << "Warning: Porosity of solid phase is 0.0 ! "
                                 "Change Settings in *.mmp file "
                              << "\n";

        /********************************************************/
        // Set up vector is_a_CCBC
        CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
        if (m_msh == NULL)
        {
            cout << "No mesh in KRConfig"
                 << "\n";
            exit(1);
        }

        const size_t mesh_node_vector_size(m_msh->nod_vector.size());
        // Initialize vector is_a_CCBC
        for (size_t l = 0; l < mesh_node_vector_size; l++)
            m_krd->is_a_CCBC.push_back(false);
        // Go through specified geometry elements
        std::string s_geo_name, s_geo_type;
        size_t s_geo_id;

        // 1) AllowReactions: This allows reactions only on specific geometries

        // First, set all nodes in vector is_a_CCBC inactive, if
        // $ALLOW_REACTIONS was set
        if (m_krd->AllowReactGeoName.size() > 0)
            for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
                m_krd->is_a_CCBC[l] = true;

        // then, activate the nodes on geometries allowed for reactions
        for (size_t j = 0; j < m_krd->AllowReactGeoName.size(); j++)
        {
            s_geo_name = m_krd->AllowReactGeoName[j];
            s_geo_type = m_krd->AllowReactGeoType[j];
            s_geo_id = m_krd->AllowReactGeoID[j];
            //------------------------------------------------------------------
            if (s_geo_type.compare("POINT") == 0)
            {
                const std::vector<GEOLIB::Point*>* pnt_vec(
                    geo_obj.getPointVec(unique_name));
                lm = m_msh->GetNODOnPNT((*pnt_vec)[m_krd->AllowReactGeoID[j]]);
                m_krd->is_a_CCBC[lm] = false;
            }
            //------------------------------------------------------------------
            if (s_geo_type.compare("POLYLINE") == 0)
            {
                // CGLPolyline *ply (polyline_vector[s_geo_id]);
                std::vector<GEOLIB::Polyline*> const* const ply_vec(
                    geo_obj.getPolylineVec(unique_name));
                GEOLIB::Polyline const* const ply((*ply_vec)[s_geo_id]);
                if (ply)
                {
                    // if (ply->getType() == 100)         //WW
                    //   m_msh->GetNodesOnArc(ply, nodes_vector);
                    // else
                    m_msh->GetNODOnPLY(ply, nodes_vector);
                    for (size_t i = 0; i < nodes_vector.size(); i++)
                    {
                        ll = nodes_vector[i];
                        lm = ll;  //+ShiftInNodeVector;
                        m_krd->is_a_CCBC[lm] = false;
                    }
                }
            }  // if(POLYLINE)
            //------------------------------------------------------------------
            if (s_geo_type.compare("SURFACE") == 0)
            {
                Surface* m_surface = NULL;
                m_surface = GEOGetSFCByName(s_geo_name);
                if (m_surface)
                {
                    m_msh->GetNODOnSFC(
                        m_surface, nodes_vector);  // kg44 if used together with
                                                   // PETSC one maybe needs to
                    // set the for_ic flag...see rf_id_new.cpp
                    for (size_t i = 0; i < nodes_vector.size(); i++)
                    {
                        ll = nodes_vector[i];
                        lm = ll;  //+ShiftInNodeVector;
                        m_krd->is_a_CCBC[lm] = false;
                    }
                }
            }
        }  // end of for(j=0;j<m_krd->NoReactGeoName.size()...
        // test output
        /*  cout << " Vector KinReactData::is_a_CCBC: " << "\n";
         for(l=0; l< (long)m_msh->nod_vector.size();l++)
         if(m_krd->is_a_CCBC[l] == true) cout << l <<", ";
         cout << "\n";
         */
        nodes_vector.clear();

        // 2) NoReactions, this suppresses reactions on specific geometries
        //    can be used in combination with 1)
        for (size_t j = 0; j < m_krd->NoReactGeoName.size(); j++)
        {
            s_geo_name = m_krd->NoReactGeoName[j];
            s_geo_type = m_krd->NoReactGeoType[j];
            s_geo_id = m_krd->NoReactGeoID[j];
            //------------------------------------------------------------------
            if (s_geo_type.compare("POINT") == 0)
            {
                const std::vector<GEOLIB::Point*>* pnt_vec(
                    geo_obj.getPointVec(unique_name));
                if (pnt_vec)
                    m_krd->is_a_CCBC[m_msh->GetNODOnPNT(
                        (*pnt_vec)[m_krd->NoReactGeoID[j]])] = true;
            }
            //------------------------------------------------------------------
            if (s_geo_type.compare("POLYLINE") == 0)
            {
                CGLPolyline* ply(polyline_vector[s_geo_id]);

                std::vector<GEOLIB::Polyline*> const* const ply_vec(
                    geo_obj.getPolylineVec(unique_name));
                GEOLIB::Polyline const* const polyline((*ply_vec)[s_geo_id]);
                double msh_min_edge_length(m_msh->getMinEdgeLength());
                m_msh->setMinEdgeLength(ply->epsilon);
                m_msh->GetNODOnPLY(polyline, nodes_vector);
                m_msh->setMinEdgeLength(msh_min_edge_length);
                for (size_t i = 0; i < nodes_vector.size(); i++)
                    //					m_krd->is_a_CCBC[nodes_vector[i]+ShiftInNodeVector]
                    //= true;
                    m_krd->is_a_CCBC[nodes_vector[i]] = true;
                nodes_vector.clear();
            }  // if(POLYLINE)
            //------------------------------------------------------------------
            if (s_geo_type.compare("SURFACE") == 0)
            {
                Surface* m_surface = NULL;
                m_surface = GEOGetSFCByName(s_geo_name);
                if (m_surface)
                {
                    m_msh->GetNODOnSFC(
                        m_surface, nodes_vector);  // kg44 if used together with
                                                   // PETSC one maybe needs to
                    // set the for_ic flag...see rf_id_new.cpp
                    for (size_t i = 0; i < nodes_vector.size(); i++)
                    {
                        long ll = nodes_vector[i];  //+ShiftInNodeVector;
                        m_krd->is_a_CCBC[ll] = true;
                    }
                }
            }
        }  // end of for(j=0;j<m_krd->NoReactGeoName.size()...
        // test output
        /*  cout << " Vector KinReactData::is_a_CCBC: " << "\n";
         for(l=0; l< (long)m_msh->nod_vector.size();l++)
         if(m_krd->is_a_CCBC[l] == true) cout << l <<", ";
         cout << "\n";
         */
        nodes_vector.clear();

        /********************************************************/
        // Set up vectors switched_off_node for individual reactions
        m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
        if (m_msh == NULL)
        {
            std::cout << "No mesh in KRConfig"
                      << "\n";
            exit(1);
        }
        // for all reactions
        for (int k = 0; k < m_krd->NumberReactions; k++)
        {
            m_kr = KinReact_vector[k];
            if (!m_kr->NotThisReactGeoName.empty())
            {
                // Initialize vector is_a_CCBC
                for (size_t l = 0; l < mesh_node_vector_size; l++)
                    m_kr->switched_off_node.push_back(false);
                // Go through specified geometry elements
                for (size_t j = 0; j < m_kr->NotThisReactGeoName.size(); j++)
                {
                    s_geo_name = m_kr->NotThisReactGeoName[j];
                    s_geo_type = m_kr->NotThisReactGeoType[j];
                    s_geo_id = m_kr->NotThisReactGeoID[j];
                    //------------------------------------------------------------------
                    if (s_geo_type.compare("POINT") == 0)
                    {
                        const std::vector<GEOLIB::Point*>* pnt_vec(
                            geo_obj.getPointVec(unique_name));
                        m_kr->switched_off_node[m_msh->GetNODOnPNT(
                            (*pnt_vec)[m_kr->NotThisReactGeoID[j]])] = true;
                    }
                    //------------------------------------------------------------------
                    if (s_geo_type.compare("POLYLINE") == 0)
                    {
                        std::vector<GEOLIB::Polyline*> const* const ply_vec(
                            geo_obj.getPolylineVec(unique_name));
                        GEOLIB::Polyline const* const polyline(
                            (*ply_vec)[s_geo_id]);
                        m_msh->GetNODOnPLY(polyline, nodes_vector);
                        for (size_t i = 0; i < nodes_vector.size(); i++)
                        {
                            long ll = nodes_vector[i];  //+ShiftInNodeVector;
                            m_krd->is_a_CCBC[ll] = true;
                        }
                    }

                    if (s_geo_type.compare("SURFACE") == 0)
                    {
                        Surface* m_surface = NULL;
                        m_surface = GEOGetSFCByName(s_geo_name);
                        if (m_surface)
                        {
                            m_msh->GetNODOnSFC(
                                m_surface,
                                nodes_vector);  // kg44 if used together with
                                                // PETSC one maybe
                            // needs to set the for_ic flag...see
                            // rf_id_new.cpp
                            for (size_t i = 0; i < nodes_vector.size(); i++)
                            {
                                long ll =
                                    nodes_vector[i];  //+ShiftInNodeVector;
                                m_kr->switched_off_node[ll] = true;
                            }
                        }
                    }
                    //------------------------------------------------------------------
                }  // loop over j NotThisReactGeoName elements
                nodes_vector.clear();
            }  // if NotThisReactGeoName.size > 0
        }      // loop over k nreactions

        /********************************************************/
        // Get foc average of connected elements
        MeshLib::CNode* m_nod = NULL;
        MeshLib::CElem* m_ele = NULL;
        for (size_t l = 0; l < mesh_node_vector_size; l++)
            m_krd->node_foc.push_back(1.0);
        for (size_t l = 0; l < mesh_node_vector_size; l++)
        {
            m_nod = m_msh->nod_vector[l];
            double ww(0.0);
            double w(0.0);
            for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
            {
                m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
                m_mat_mp = mmp_vector[m_ele->GetPatchIndex()];
                ww += m_mat_mp->foc * m_ele->GetVolume();
                w += m_ele->GetVolume();
            }
            double foc(ww / w);
            //// CB Brand project
            // if(foc>0.1)
            //  foc = 1;
            // else
            //  foc = 1e-100;
            m_krd->node_foc[l] = foc;
        }

        /********************************************************/
        // CB Set neighborhood relations for reaction deactivation switch
        if ((m_krd->ReactDeactEpsilon) < 0)
            m_krd->ReactDeactFlag = false;
        else
            m_krd->ReactDeactFlag = true;
        if (m_krd->ReactDeactFlag)
        {
            for (size_t l = 0; l < mesh_node_vector_size; l++)
            {
                // initialize all nodes as active
                m_krd->ReactDeact.push_back(false);
                // initialize all reaction terms as zero
                m_krd->React_dCdT.push_back(0);
                m_nod = m_msh->nod_vector[l];  // get current node object

                // loop over PRIMARY NEIGHBOR elements
                // this routine takes longer, but gets all the nodes
                for (size_t ll = 0; ll < m_nod->getConnectedElementIDs().size();
                     ll++)
                {
                    // get index of direct neighbour element
                    size_t dnele_idx = m_nod->getConnectedElementIDs()[ll];
                    // get direct neighbour element object
                    m_dnele = m_msh->ele_vector[dnele_idx];
                    // get the neighbor element node indices
                    m_dnele->GetNodeIndeces(primnnodesindices);
                    // get the neighbor element number of nodes
                    int nnodpneigh = m_dnele->GetNodesNumber(false);
                    // loop over primary neighbour element number of nodes
                    for (int lll = 0; lll < nnodpneigh; lll++)
                    {
                        // get the current node index
                        long pneighnod_idx = primnnodesindices[lll];
                        // get the node object
                        m_dnnod = m_msh->nod_vector[pneighnod_idx];
                        // loop over the connected elements, this now includes
                        // the SECONDARY NEIGHBORS
                        for (size_t llll = 0;
                             llll < m_dnnod->getConnectedElementIDs().size();
                             llll++)
                        {
                            // get the current connected element indices
                            size_t snele_idx =
                                m_dnnod->getConnectedElementIDs()[llll];
                            // get secondary neighbour element object
                            m_snele = m_msh->ele_vector[snele_idx];
                            // get the neighbor element node indices
                            m_snele->GetNodeIndeces(secnnodesindices);
                            // get the neighbor element number of nodes
                            int nnodsneigh = m_snele->GetNodesNumber(false);
                            // loop over secondary neighbour element number of
                            // nodes
                            for (int lllll = 0; lllll < nnodsneigh; lllll++)
                            {
                                bool duplicate(false);  // flag, initialize
                                int annode_idx = secnnodesindices[lllll];
                                // check vs. all previously identified nodes to
                                // prevent duplicate entries
                                for (size_t nn = 0;
                                     nn < ReactNeighborNodes.size();
                                     nn++)
                                    // node has already been found
                                    if (annode_idx == ReactNeighborNodes[nn])
                                    {
                                        duplicate =
                                            true;  // set flag to "not fresh"
                                        break;     // skip rest of comparisons
                                    }
                                if (!duplicate)  // fresh node
                                    // push back node index
                                    ReactNeighborNodes.push_back(annode_idx);
                                else
                                    duplicate = false;  // was no fresh node,
                                                        // reinitialize flag
                            }
                        }
                    }
                }
                // Add local node neighborhood indices vector to global vector
                // This is the same for both above neighborhood search routines
                m_krd->ReactNeighborhood.push_back(ReactNeighborNodes);
                ReactNeighborNodes.clear();
            }
            // debug
            // for (nn=0;nn<m_krd->ReactNeighborhood.size();nn++){
            //  cout << nn << " " << m_krd->ReactNeighborhood[nn].size() << " "
            //  ; for (ll=0;ll<m_krd->ReactNeighborhood[nn].size();ll++)
            //    cout << m_krd->ReactNeighborhood[nn][ll] << " " ;
            //  cout << "\n";
            //}
            m_krd->concentrationmatrix = (double**)malloc(
                ((long)mesh_node_vector_size) * sizeof(double*));
            for (size_t ix = 0; ix < mesh_node_vector_size; ix++)
                m_krd->concentrationmatrix[ix] =
                    (double*)malloc((length) * sizeof(double));
            for (size_t ix = 0; ix < mesh_node_vector_size; ix++)
                for (size_t ixx = 0; ixx < length; ixx++)
                    m_krd->concentrationmatrix[ix][ixx] = 0;
        }  // if (m_krd->ReactDeactFlag)

        /********************************************************/
        if (m_krd->debugoutflag)
        {
            m_krd->debugoutfilename = FileName + "_KR_Debug_out.txt";
            m_krd->c_dumpfilename = FileName + "_KR_C_dump.txt";
        }
        /********************************************************/
        if (m_krd->sortnodes)
        {
            m_krd->substeps = new long[(long)m_msh->nod_vector.size()];
            m_krd->noderanks = new long[(long)m_msh->nod_vector.size()];
            for (long ix = 0; ix < ((long)m_msh->nod_vector.size()); ix++)
            {
                m_krd->noderanks[ix] = ix;
                m_krd->substeps[ix] = 0;
            }
        }
        /********************************************************/
    }  // if KinReactData_vector.size() > 0
}

/**************************************************************************
   FEMLib-Method:
   Task: Clear KinReact_vector
   Programing:
   02/2006 SB Implementation
   last modified:
**************************************************************************/
void KRCDelete(void)
{
    long i;
    int no_krc = (int)KinReact_vector.size();
    for (i = 0; i < no_krc; i++)
        delete KinReact_vector[i];
    KinReact_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB New Class structure, IO, C++, FEM
   06/2010 TF formated, removed unused variables, changed signature for new
GEOLIB data structures
**************************************************************************/
bool CKinReact::Read(std::ifstream* rfd_file, const GEOLIB::GEOObjects& geo_obj,
                     const std::string& unique_name)
{
    char line[MAX_ZEILE];
    std::string line_string, line_str1;
    bool new_keyword = false, new_subkeyword = false;
    std::stringstream in;
    MonodSubstruct *m_monod = NULL, *m_inhibit = NULL, *m_production = NULL;
    // MicrobeData* m_microbe = NULL;
    size_t index, index1;
    double dval;
    // CB 10/09
    std::string s_geo_type, s_geo_name;
    string thresh_species;
    double thresh_conc, thresh_ord;
    bool found = false;
    int i;
    double dummy = 0;
    // clear vectors
    monod.clear();
    inhibit.clear();
    production.clear();
    reactionpartner.clear();
    stochmet.clear();
    T_params.clear();
    mechvec.clear();        // CB minkin
    minSpeciesIdx.clear();  // CB minkin
    MinkinMech* m_mech = NULL;
    double expo = -99;
    string speciesname = "Nix";
    string dummystring;

    while (!new_keyword)
    {
        index = rfd_file->tellg();
        if (!GetLineFromFile(line, rfd_file))
            break;
        line_string = line;
        if (line_string.find("#") != string::npos)
        {
            new_keyword = true;
            rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
                                     // n�chste keyword weg
            break;
        }
        /* Keywords nacheinander durchsuchen */
        //....................................................................
        // subkeyword found
        if (line_string.find("$NAME") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> name;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$TYPE") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> type;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$BACTERIANAME") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> bacteria_name;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$EQUATION") != string::npos)
        {
            line_str1 = GetLineFromFile1(rfd_file);
            ReadReactionEquation(line_str1);
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$RATECONSTANT") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> rateconstant >> rateorder;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$GROWTH") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> grow;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$MONODTERMS") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_keyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste subkeyword weg
                    break;
                }
                m_monod = new MonodSubstruct();
                number_monod++;
                in.str(line_str1);
                in >> m_monod->species >> m_monod->concentration >>
                    m_monod->order;
                if ((m_monod->order != -99.0) &&
                    (m_monod->concentration != -1.0e9))  // check for read in
                    monod.push_back(m_monod);
                else
                {
                    Display::DisplayMsgLn(" ERROR reading Monod Terms  - skipping");
                    number_monod--;
                    delete m_monod;
                }
                in.clear();
            }
        //....................................................................
        // subkeyword found
        if (line_string.find("$THRESHHOLDTERMS") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_keyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste subkeyword weg
                    break;
                }
                in.str(line_str1);
                in >> thresh_species >> thresh_conc >> thresh_ord;
                // check for read in
                if ((thresh_ord != -99.0) && (thresh_conc != -1.0e9))
                {
                    found = false;
                    for (i = 0; i < number_monod; i++)
                        // find the respective monod species
                        if (thresh_species.compare(monod[i]->species) == 0)
                        {
                            monod[i]->threshhold = true;  // and store the data
                            monod[i]->threshConc = thresh_conc;
                            monod[i]->threshOrder = thresh_ord;
                            found = true;
                            break;
                        }
                    if (found == false)
                    {
                        cout << " WARNING: no matching MONOD SPECIES found in "
                                "reaction ";
                        cout << name << " for THRESHHOLD TERM SPECIES "
                             << thresh_species << "\n";
                    }
                }
                else
                    Display::DisplayMsgLn(" ERROR reading Threshhold Terms  - skipping");
                in.clear();
            }
        //....................................................................
        // subkeyword found
        if (line_string.find("$INHIBITIONTERMS") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_keyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste subkeyword weg
                    break;
                }
                in.str(line_str1);
                m_inhibit = new MonodSubstruct();
                number_inhibit++;
                in >> m_inhibit->species;
                in >> m_inhibit->concentration;
                in >> m_inhibit->order;
                if ((m_inhibit->order != -99.0) &&
                    (m_inhibit->concentration != -1.0e9))  // check for read in
                    inhibit.push_back(m_inhibit);
                else
                {
                    delete m_inhibit;  // CB_merge_0513 ??
                    Display::DisplayMsgLn(" ERROR reading Inhibition Terms  - skipping");
                    number_inhibit--;
                }

                in.clear();
            }
        //....................................................................
        // subkeyword found
        if (line_string.find("$PRODUCTIONTERMS") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_keyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste subkeyword weg
                    break;
                }
                in.str(line_str1);
                m_production = new MonodSubstruct();
                number_production++;
                in >> m_production->species >> m_production->concentration >>
                    m_production->order;
                if ((m_production->order != -99.0)
                    // check for read in
                    && (m_production->concentration != -1.0e9))
                    production.push_back(m_production);
                else
                {
                    Display::DisplayMsgLn(" ERROR reading Production Terms  - skipping");
                    number_production--;
                }
                in.clear();
            }

        //....................................................................
        // subkeyword found
        if (line_string.find("$PRODUCTIONSTOCH") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_subkeyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste subkeyword weg
                    break;
                }
                in.str(line_str1);
                m_production = new MonodSubstruct();
                in >> m_production->species >> m_production->concentration;
                ProdStochhelp.push_back(m_production);
                in.clear();
            }
        //....................................................................
        if (line_string.find("$BACTERIAL_YIELD") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Yieldcoefficient;
            in.clear();
        }
        //....................................................................
        // CB Iso routine for reading in the isotope couples
        // subkeyword found
        if (line_string.find("$ISOTOPE_FRACTIONATION") != string::npos)
        {
            number_isotope_couples++;
            in.str(GetLineFromFile1(rfd_file));
            in >> Isotope_light >> Isotope_heavy >> isoenfac;
            in.clear();
            degType = "isotope_fractionation";  // CB besser in KRConfig??
        }
        //....................................................................
        if (line_string.find("$BACTERIA_SPECIFIC_CAPACITY") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> specif_cap;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$TEMPERATURE_DEPENDENCE") != string::npos)
        {  // subkeyword found
            T_dependence = true;
            in.str(GetLineFromFile1(rfd_file));
            in >> T_model;
            //
            if (T_model == 0)
            {
                T_dependence = false;
                cout << " Warning - Model 0: Reaction assumed independent of "
                        "Temperature. "
                     << "\n";
            }
            else if ((T_model == 1) || (T_model == 2))
            {
                if (type.compare("monod") == 0)
                {
                    if (T_model == 1)
                    {  // Biomass growth, Ratkowsky2 : Zwietering et al.
                        in >> dummy;  // Ratkowsky Tmin
                        T_params.push_back(dummy);
                        in >> dummy;  // Ratkowsky Tmax
                        T_params.push_back(dummy);
                        in >> dummy;  // Ratkowsky b
                        T_params.push_back(dummy);
                        in >> dummy;  // Ratkowsky c
                        T_params.push_back(dummy);
                    }
                    else if (T_model == 2)
                    {
                        cout << " Unknown model for Temperature Dependence of "
                                "Reaction "
                             << name << " in CKinReact::Read(). "
                             << "\n";
                        cout << " Aborting now!"
                             << "\n";
                        exit(0);
                    }
                }
                else if (type.compare("NAPLdissolution") == 0)
                {  // NAPL dissolution, Chen et al (1) or Knauss et al (2) model
                    in >> dummy;  // Knauss A
                    T_params.push_back(dummy);
                    in >> dummy;  // Knauss B
                    T_params.push_back(dummy);
                    in >> dummy;  // Knauss C
                    T_params.push_back(dummy);
                }
                else if (type.compare("Mineralkinetics") == 0)
                {
                    cout
                        << " Warning in CKinReact::Read() - Arrhenius-Equation "
                           "is used for Temperature correction of "
                           "Mineral Reaction "
                        << name << "\n";
                    cout << " Input values for $TEMPERATURE_DEPENDENCE keyword "
                            "will be ignored. "
                         << "\n";
                }
                else
                {
                    cout << " Warning in CKinReact::Read() - Temperature "
                            "correction of "
                         << type << " Reaction " << name
                         << " is not implemented. "
                         << "\n";
                    cout << " Input values for $TEMPERATURE_DEPENDENCE keyword "
                            "will be ignored. "
                         << "\n";
                }
            }
            else
            {
                cout << " Unknown model for Temperature Dependence of Reaction "
                     << name << " in CKinReact::Read(). "
                     << "\n";
                cout << " Aborting now!"
                     << "\n";
                exit(0);
            }

            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$_drmc_") != string::npos)
        {
            _drmc_ = true;
            in.str(GetLineFromFile1(rfd_file));
            in >> Microbename >> DormType;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$STANDARD_GIBBS_ENERGY") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> dG0;
            in.clear();
        }

        //....................................................................
        // subkeyword found
        if (line_string.find("$EXCHANGE_PARAMETERS") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> dval;  // velocity
            ex_param.push_back(dval);
            in >> dval;  // kd
            ex_param.push_back(dval);
            if (exType.compare("langmuir") == 0)
                in >> exSurfaceID;
            if (exType.compare("freundlich") == 0)
            {
                in >> dval;
                ex_param.push_back(dval);
                dval = 1;
                in >> dval;
                ex_param.push_back(dval);
            }
            if (exType.compare("linear") == 0)
            {
                dval = 1;
                in >> dval;
                ex_param.push_back(dval);
            }

            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$SORPTION_TYPE") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> exType;
            in.clear();
        }

        //....................................................................
        // subkeyword found
        if (line_string.find("$NAPL_PROPERTIES") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> blob_name;  // >> Csat_pure >> Density_NAPL;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$REACTION_ORDER") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> rateorder;  // CB rateorder added 06/2012
            in.clear();
        }
        //....................................................................
        if (line_string.find("$MINERALNAME") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> mineral_name;
            in.clear();
        }
        if (line_string.find("$CHEMAPPNAME") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> chemapp_name;
            in.clear();
        }
        if (line_string.find("$EQUILIBRIUM_CONSTANT") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> dummystring >> Km_default;
            if (dummystring.compare("UNIFORM") == 0)
            {
                Km_uniform = true;
                Km_CHEMAPP = false;
                Km_HKF = false;
            }
            else if (dummystring.compare("CHEMAPP") == 0)
            {
                Km_uniform = false;
                Km_CHEMAPP = true;
                Km_HKF = false;
            }
            else if (dummystring.compare("HKF") == 0)
            {
                Km_uniform = false;
                Km_CHEMAPP = false;
                Km_HKF = true;
            }
            else
                cout << " Warning in CKinReact::Read - No valid model for "
                        "Mineralkinetics equilibrium constant Km."
                     << "\n";
            // anyway, pushback the input value
            Km.push_back(pow(10, Km_default));
            in.clear();
        }
        //....................................................................
        if (line_string.find("$RATE_EXPONENTS") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Theta >> Eta;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$REACTIVE_SURFACE_AREA") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Am_model >> Am_ini;
            if (Am_model == 0)
            {
                Am_constant = true;
            }
            else if (Am_model <= 5)
            {
                Am_constant = false;
            }
            else
                cout << " Warning in CKinReact::Read - No valid model for "
                        "Mineralkinetics reactive surface area Am."
                     << "\n";
            // anyway, pushback the input value
            Am.push_back(Am_ini);
            in.clear();
        }
        //....................................................................
        if (line_string.find("$PRECIPITATION_BY_BASETERM_ONLY") != string::npos)
        {  // subkeyword found
            precip_baseterm_only = true;
        }
        //....................................................................
        if (line_string.find("$PRECIPITATION_FACTOR") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> this->precipfactor;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$PRECIPITATION_EXPONENT") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> this->precipexponent;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$BASETERM") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            m_mech = new MinkinMech();
            number_Mech++;
            in >> m_mech->k25 >> m_mech->Eact;
            in.clear();
            // convert log(k25) --> k25
            m_mech->k25 = pow(10, m_mech->k25);
            mechvec.push_back(m_mech);
        }
        //....................................................................
        if (line_string.find("$MECHANISMTERM") != string::npos)
        {  // subkeyword found
            // before reading in k mechanistic species, store ratteconstant and
            // Eact
            m_mech = new MinkinMech();
            number_Mech++;
            in.str(GetLineFromFile1(rfd_file));
            in >> m_mech->k25 >> m_mech->Eact;
            in.clear();
            // convert log(k25) --> k25
            m_mech->k25 = pow(10, m_mech->k25);
            if ((m_mech->k25 != 0) /*&&(m_mech->Eact != 0)*/)
            {  // check for read in
                // now read mechanistic species
                while (!new_subkeyword)
                {
                    index1 = rfd_file->tellg();
                    line_str1 = GetLineFromFile1(rfd_file);
                    // Check for end of data block
                    if ((line_str1.find("#") != string::npos) ||
                        (line_str1.find("$") != string::npos))
                    {
                        if (line_str1.find("#") != string::npos)
                            new_keyword = true;
                        rfd_file->seekg(
                            index1);  // Dateipointer zurücksetzen, sonst ist
                                      // das nächste subkeyword weg
                        break;
                    }
                    // Here, read individual species data
                    in.str(line_str1);
                    in >> speciesname >> expo;
                    if ((speciesname.compare("NIX") != 0) && (expo != -99))
                    {  // check for read in
                        m_mech->no_mechSpec++;
                        m_mech->mechSpeciesNames.push_back(speciesname);
                        m_mech->mechSpeciesExpo.push_back(expo);
                        expo = -99;
                        speciesname = "NIX";
                    }
                    else
                        Display::DisplayMsgLn(
                            " ERROR reading Mechanism Species  - skipping");
                    in.clear();
                }
                in.clear();
                mechvec.push_back(m_mech);
            }  // successfull read
            else
            {
                Display::DisplayMsgLn(" ERROR reading Mechanism Terms  - skipping");
                number_Mech--;
                delete m_mech;
            }
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$SWITCH_OFF_GEOMETRY") != string::npos)
            while (!new_subkeyword)
            {
                index1 = rfd_file->tellg();
                line_str1 = GetLineFromFile1(rfd_file);
                // Check for end of data block
                if ((line_str1.find("#") != string::npos) ||
                    (line_str1.find("$") != string::npos))
                {
                    if (line_str1.find("#") != string::npos)
                        new_keyword = true;
                    rfd_file->seekg(
                        index1);  // Dateipointer zurueksetzen, sonst ist das
                                  // naehste subkeyword weg
                    break;
                }
                in.str(line_str1);
                in >> s_geo_type >> s_geo_name;
                NotThisReactGeoType.push_back(s_geo_type);
                NotThisReactGeoName.push_back(s_geo_name);
                size_t geo_obj_idx(std::numeric_limits<size_t>::max());

                if (s_geo_type.find("POINT") != std::string::npos)
                    // TF 06/2010 - get the point vector and set the geo_obj_idx
                    if (!((geo_obj.getPointVecObj(unique_name))
                              ->getElementIDByName(s_geo_name, geo_obj_idx)))
                    {
                        std::cerr
                            << "error in CKinReact::Read: point name not found!"
                            << "\n";
                        exit(1);
                    }
                if (s_geo_type.find("POLYLINE") != std::string::npos)
                    // TF 07/2010 - get the polyline vector and set the
                    // geo_obj_idx
                    if (!((geo_obj.getPolylineVecObj(unique_name))
                              ->getElementIDByName(s_geo_name, geo_obj_idx)))
                    {
                        std::cerr << "error in CKinReact::Read: polyline name "
                                     "not found!"
                                  << "\n";
                        exit(1);
                    }
                NotThisReactGeoID.push_back(geo_obj_idx);

                in.clear();
            }
    }
    return true;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Check consistency of reactions read in from input or database file
   Programing:
   02/2006 SB Implementation
   last modified:
**************************************************************************/
int CKinReact::CheckReactionDataConsistency()
{
    int i, ok = 1, length;
    bool found = false;
    string name, name1;
    //    CRFProcess *m_pcs = NULL;
    //    int no_processes =(int)pcs_vector.size();

    // Check number of reaction partners
    if (number_reactionpartner < 2)
    {
        // SB todo		ok=0;
        // SB todo		cout << " Less than two reaction partners in equation
        // found"
        // << "\n";
    }
    // Ratekonstant zero?
    if ((fabs(rateconstant) < MKleinsteZahl) && (type.compare("monod") == 0))
    {
        //#ds erstetzt statt:	if((fabs(rateconstant) < MKleinsteZahl) &&
        //(type.compare("exchange")!=0) ){
        ok = 0;
        cout << " Rateconstant is zero"
             << "\n";
    }
    // Number of partners higher than in arrays?
    length = (int)reactionpartner.size();  // Muss so gemacht werden, weil sonst
                                           // signed/unsigned warning kommt
    if (number_reactionpartner != length)
    {
        ok = 0;
        cout << " Not enough reaction partners"
             << "\n";
        for (i = 0; i < length; i++)
            cout << reactionpartner[i] << "  ";
        cout << "\n";
    }
    length = (int)stochmet.size();
    if (stochmet.size() != reactionpartner.size())
    {
        ok = 0;
        cout << " Not enough stochieometric coefficients for reaction partners"
             << "\n";
        for (i = 0; i < length; i++)
            cout << stochmet[i] << "  ";
        cout << "\n";
    }
    // check type
    if (type.compare("monod") != 0)
    {
        if (type.compare("exchange") != 0)
            if (type.compare("NAPLdissolution") != 0)
                if (type.compare("Mineralkinetics") != 0)
                {
                    ok = 0;
                    cout << "Unknown reaction type"
                         << "\n";
                }
    }

    /* Check Monod-, Inhibition and Production terms */
    length = (int)monod.size();
    if (length < number_monod)
        number_monod = (int)monod.size();
    length = (int)inhibit.size();
    if (length < number_inhibit)
        number_inhibit = (int)inhibit.size();
    length = (int)production.size();
    if (length < number_production)
        number_production = (int)production.size();

    if ((number_monod + number_inhibit + number_production == 0) &&
        (type == "monod"))
    {
        //		cout << "Warning: no monod terms specified for monod equation
        //"<<
        //"\n";
    }

    // Test concentrations and constants for > 0.0
    length = (int)monod.size();
    for (i = 0; i < length; i++)
        if (monod[i]->concentration < MKleinsteZahl)
        {
            cout << " Monod Concentration of reaction " << name
                 << " smaller than 0.0  "
                 << "\n";
            ok = 0;
        }
    length = (int)inhibit.size();
    for (i = 0; i < length; i++)
        if (inhibit[i]->concentration < MKleinsteZahl)
        {
            cout << " Inhibition Concentration of reaction " << name
                 << " smaller than 0.0  "
                 << "\n";
            ok = 0;
        }
    //#ds check NAPL dissolution reactions
    if (type.compare("NAPLdissolution") == 0)
    {
        // does blob-name exist?
        CKinBlob* m_kb = NULL;
        for (i = 0; i < (int)KinBlob_vector.size(); i++)
        {
            m_kb = KinBlob_vector[i];
            if (blob_name.compare(m_kb->name) == 0)
                found = true;
        }
        if (found == false)
        {
            ok = 0;
            cout << " Blob_Name " << blob_name << " defined in Reaction "
                 << name << " does not exist "
                 << "\n";
        }
        // Csat_pure >= 0
        // if (Csat_pure < -MKleinsteZahl)
        if (cp_vec[cp_name_2_idx[reactionpartner[0]]]->max_solubility <
            -MKleinsteZahl)
        {
            ok = 0;
            cout << " Invalid maximum solubility Csat_pure: "
                 << cp_vec[cp_name_2_idx[reactionpartner[0]]]
                        ->max_solubility /*Csat_pure*/
                 << " in Reaction " << name << "\n";
        }
        // Density_NAPL >0
        // if (Density_NAPL < MKleinsteZahl)
        if (cp_vec[cp_name_2_idx[reactionpartner[0]]]->molar_density <
            MKleinsteZahl)
        {
            ok = 0;
            cout << " Invalid NAPL density Density_NAPL: "
                 << cp_vec[cp_name_2_idx[reactionpartner[0]]]
                        ->molar_density /*Density_NAPL*/
                 << " in Reaction " << name << "\n";
        }
    }  // end type=NAPLdissolution

    if (type.compare("Mineralkinetics") == 0)
    {
        // is mineral defined as solid phase species?
        // TODO CB_merge_0513 ??
    }
    return ok;
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::Write(ofstream* rfe_file)
{
    int i, flag = 0, length;

    // Write Keyword
    *rfe_file << "#REACTION"
              << "\n";
    // Name of reaction
    *rfe_file << "$NAME"
              << "\n"
              << name << "\n";
    // Type of reaction
    *rfe_file << "$TYPE"
              << "\n"
              << type << "\n";
    // bacteria name
    if (type == "monod")
        *rfe_file << "$BACTERIANAME"
                  << "\n"
                  << bacteria_name << "\n";
    if (type == "exchange")
        *rfe_file << "$SORPTION_TYPE"
                  << "\n"
                  << exType << "\n";
    // ReactionEquation
    *rfe_file << "$EQUATION"
              << "\n";
    for (i = 0; i < number_reactionpartner; i++)
    {
        if (stochmet[i] < 0.0)  // left side of equation
        {
            if (i == 0)
                *rfe_file << " " << fabs(stochmet[i]) << " "
                          << reactionpartner[i];
            else
                *rfe_file << " + " << fabs(stochmet[i]) << " "
                          << reactionpartner[i];
        }
        if (stochmet[i] > 0 && (flag > 0))  // remaining right hand side
            *rfe_file << " + " << fabs(stochmet[i]) << " "
                      << reactionpartner[i];
        if (stochmet[i] > 0 &&
            (flag == 0))  // " = " Sign and first term on right hand side
        {
            *rfe_file << " = " << fabs(stochmet[i]) << " "
                      << reactionpartner[i];
            flag = 1;
        }
    }
    *rfe_file << "\n";
    // Rateconstant and order
    if (type == "monod")
    {
        *rfe_file << "$RATECONSTANT"
                  << "\n"
                  << rateconstant << "   " << rateorder << "\n";
        *rfe_file << "$GROWTH"
                  << "\n"
                  << grow << "\n";
        // Monod terms
        *rfe_file << "$MONODTERMS"
                  << "\n";  //<< number_monod << "\n";
        for (i = 0; i < number_monod; i++)
            *rfe_file << monod[i]->species << "  " << monod[i]->concentration
                      << "  " << monod[i]->order << "\n";
        // Inhibition terms
        *rfe_file << "$INHIBITIONTERMS"
                  << "\n";  // << number_inhibit << "\n";
        for (i = 0; i < number_inhibit; i++)
            *rfe_file << inhibit[i]->species << "  "
                      << inhibit[i]->concentration << "  " << inhibit[i]->order
                      << "\n";
        // Production Terms
        //*rfe_file << "$PRODUCTIONTERMS" << "\n" << number_production << "\n";
        // for(i=0;i<number_production;i++)
        //	*rfe_file << production[i]->species << "  " <<
        // production[i]->concentration << "  " << production[i]->order
        //<< "\n";
        // Production Terms
        length = (int)ProdStochhelp.size();
        *rfe_file << "$PRODUCTIONSTOCH"
                  << "\n";  // << length << "\n";
        for (i = 0; i < length; i++)
            *rfe_file << ProdStochhelp[i]->species << "  "
                      << ProdStochhelp[i]->concentration << "  "
                      << "\n";
        if (degType == "isotope_fractionation")
        {
            *rfe_file << "$ISOTOPE_FRACTIONATION"
                      << "\n";
            *rfe_file << Isotope_light << "  " << Isotope_heavy << "  "
                      << isoenfac << "\n";
        }
        //#ds output f�r NAPL-dissolution
        //	*rfe_file << "$NAPL_PROPERTIES" << "\n";
        //	*rfe_file << "blob_name " << blob_name << " Csat_pure " << Csat_pure
        //<< " Density_NAPL " << Density_NAPL <<
        //"\n";
    }
    if (type == "exchange")
    {
        *rfe_file << "EXCHANGE_PARAMETERS"
                  << "\n";
        length = (int)ex_param.size();
        for (i = 0; i < length; i++)
            *rfe_file << ex_param[i] << "  ";
        *rfe_file << "\n";
    }
    *rfe_file << "\n";
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class read function
   Programing:
   05/2004 SB Implementation
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::ReadReactionEquation(string line_string_all)
{
    string line_string, name, helpstring, c, calt;
    string substring, subsubstring;
    int indexlow, indexhigh, help1, strl, linelength, ih1, ih2, ih3,
        indexgleich, i, partners, p;
    double zahl, vz = -1.0;

    stochmet.clear();
    reactionpartner.clear();

    // Anf�ngliche Leerzeichen, Kommentar am Ende raus
    ih1 = (int)line_string_all.find_first_not_of(" ", 0);
    ih2 = (int)line_string_all.find_first_of(";", ih1);
    ih3 = (int)line_string_all.length();
    line_string = line_string_all.substr(ih1, ih2 - ih1);

    // Auf mehrfache Leerzeichen und Tabs �berpr�fen - ersetzt durch je ein
    // Leerzeichen
    ih3 = (int)line_string.length();
    for (i = 0; i < ih3; i++)
    {
        c = line_string[i];
        if (c == "	")
            c = " ";  // replace tab by space
        if ((c == " ") && (calt == " "))
            ih1++;  // nothing
        else
        {
            helpstring.append(c);
            calt = c;
        }
    }

    // Leerzeichen am Ende raus
    ih1 = (int)helpstring.find_first_not_of(" ", 0);
    ih2 = (int)helpstring.find_last_of(" ");
    ih3 = (int)helpstring.length();
    if (ih2 == ih3 - 1)  // Leerzeichen am Ende
        line_string = helpstring.substr(ih1, ih2 - ih1);
    else
        line_string = helpstring.substr(ih1, ih3 - ih1);

    // Count reaction partners by looking for " + "
    linelength = (int)line_string.length();
    indexhigh = 0;
    partners = 0;
    ih1 = (int)line_string.find(" = ");
    if (ih1 < 0)
        Display::DisplayMsgLn(" Error in keyword REACTION");
    // Exception handling
    while (indexhigh < linelength)
    {
        ih1 = (int)line_string.find(" + ", indexhigh);
        if (ih1 > 0)
        {
            partners++;
            indexhigh = ih1 + 3;
        }
        else
            // no " + " found
            indexhigh = linelength;
    }
    // Number of partners is 2 for " = " and one additional for each " + "
    number_reactionpartner = partners + 2;

    /* Zerlegen der Gleichung in Bl�che mit "Zahl Namen"*/

    indexhigh = 0;
    p = 0;  // Zaehlt partner hoch
    indexlow = indexhigh;
    linelength = (int)line_string.length();

    indexgleich = (int)line_string.find(" = ");

    while (indexhigh < linelength)
    {
        /* einen Block holen  - find next occurence of +, -, = */
        ih1 = (int)line_string.find(" + ", indexlow);
        if (ih1 < 0)
            ih1 = linelength;
        ih2 = (int)line_string.find(" - ", indexlow);
        if (ih2 < 0)
            ih2 = linelength;
        ih3 = (int)line_string.find(" = ", indexlow);
        if (ih3 < 0)
            ih3 = linelength;
        indexhigh = min(ih1, min(ih2, ih3));
        if (indexhigh > indexgleich)
            vz = 1.0;  // rechte Seite der Gleichung: Vorzeichenwechsel
        substring = line_string.substr(indexlow, indexhigh - indexlow);

        /* Leerzeichen drin ? */
        help1 = (int)substring.find(" ", 0);
        strl = (int)substring.length();
        if (help1 > 0)  // Leerzeichen gefunden
        {
            subsubstring = substring.substr(0, help1);
            //	if(subsubstring.find(".",0)>0) //double value given
            zahl = atof(subsubstring.data()) * vz;
            stochmet.push_back(zahl);
            subsubstring = substring.substr(help1 + 1, strl - help1 - 1);
            // reactionpartner[p] = subsubstring.data();
            //	sprintf(reactionpartner[p],"%s", subsubstring.data());
            reactionpartner.push_back(subsubstring);
        }  // kein Leerzeichen gefunden, substring enth�lt nur den Namen => Zahl
           // wird 1 gesetzt
        else
        {
            zahl = 1.0 * vz;
            stochmet.push_back(zahl);
            //	reactionpartner[p] = substring.data();
            //	sprintf(reactionpartner[p],"%s", subsubstring.data());
            reactionpartner.push_back(substring);
        }
        p++;  // Reactionpartner hochzaehlen
        // cout << zahl <<" "<< name << " ";

        indexlow = indexhigh + 3;  // add length of " + "
    }                              // end while

    if (p != number_reactionpartner)
        cout << "Error: Parser for kinetic equations found varying number of "
                "reaction partners"
             << "\n";
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinBlob constructor
   Programing:
   02/2007 DS Implementation
***************************************************************************/
CKinBlob::CKinBlob(void)
{
    // default= nonsense values for input, will be checked in CKinCheck()
    d50 = -1.;
    Sh_factor = -1.;
    Re_expo = -1.;
    Sc_expo = -1.;
    Geometry_expo = -1.;
    Mass = 0.;
    Volume = 0.;
    Masstransfer_k = 0.;
    current_Interfacial_area = 0.;
    BlobGeoType.clear();
    BlobGeoName.clear();
    Area_Value.clear();
    Interfacial_area.clear();
    Masstransfer_K.clear();
    gas_dissolution_flag = false;
    // New-Sherwood-Number
    shidx = 0;  // SP: Index for different Sherwood calculations
    NContent_expo = 0;
    WContent_expo = 0;
    Pe_expo = 0;
    D50_expo = 0;
    Poro_expo = 0;
    Delta_expo = 0;
    UI_expo = 0;
    NContent_ini_expo = 0;
    NContent_res_expo = 0;
    GSR_expo = 0;
    Length_expo = 0;
    Tort_expo = 0;
    Pfannkuch_constant = 0;
    beta_4 = 0;
    one_three = 0;
    two_three = 0;
    four_nine = 0;
    fife_nine = 0;
    elev_nine = 0;
    fife_three = 0;
    grain_expo = 0;
    modSherwood = false;
    Sherwood_model = false;

    dM = -1;
    dS = -1;
    UI = -1;
    grain_var = -1;
    NCont_ini = -1;
    NCont_res = -1;
    NCont_Sh1 = -1;
    NCont_Sh2 = -1;
    // mod_Reynolds = -1 ;
    GSR = -1;
    Length = -1;
    Tort = -1;
    // CB 12/09 new data structures for parallel computation
    // Exchangeterms.clear();
    // Volume_old.clear();
}

/***************************************************************************
   FEMLib-Method:
   Task: CKinBlob destructor
   Programing:
   02/2007 DS Implementation
***************************************************************************/
CKinBlob::~CKinBlob(void) {}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function for CKinBlob-Structure
   Programing:
   02/2007 DS Implementation
**************************************************************************/
bool CKinBlob::Read(ifstream* rfd_file, const GEOLIB::GEOObjects& geo_obj,
                    const std::string& unique_name)
{
    char line[MAX_ZEILE];
    string line_string, line_str1, s_geo_type, s_geo_name;
    string hash("#"), dollar("$");
    bool new_keyword = false, OK = true;
    size_t index, index1;
    double d_inivalue;
    std::stringstream in;

    // STP New-Sherwood-Number
    //========================================================================
    while (!new_keyword)
    {
        index = rfd_file->tellg();
        //    if(!rfd_file->getline(line,MAX_ZEILE)) break;
        if (!GetLineFromFile(line, rfd_file))
            break;
        line_string = line;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
                                     // n�chste keyword weg
            break;
        }
        /* Keywords nacheinander durchsuchen */
        //....................................................................
        // subkeyword found
        if (line_string.find("$NAME") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> name;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$D50") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> d50;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        // if (line_string.find("$CALC_SHERWOOD") != string::npos)
        //{
        //	in.str(GetLineFromFile1(rfd_file));
        //	in >> Sh_factor >> Re_expo >> Sc_expo;
        //	in.clear();
        //}
        //....................................................................
        // subkeyword found
        if (line_string.find("$DM") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> dM;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$DS") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> dS;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$UI") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> UI;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$NAPL_CONTENT_INI") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> NCont_ini;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$NAPL_CONTENT_RES") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> NCont_res;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$GRAIN_SPHERE_RATIO") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> GSR;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$TORTUOSITY") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Tort;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$LENGTH") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Length;  // Länge der NAPL-Quelle
            in.clear();
        }
        //....................................................................
        if (line_string.find("$CALC_SHERWOOD") != string::npos &&
            line_string.find("$CALC_SHERWOOD_MODIFIED") == string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Sh_factor >> Re_expo >> Sc_expo;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$CALC_SHERWOOD_MODIFIED") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> Pe_expo >> NContent_expo >> WContent_expo >> D50_expo >>
                Poro_expo >> Delta_expo >> UI_expo >> NContent_ini_expo >>
                NContent_res_expo >> GSR_expo >> Tort_expo >> Length_expo >>
                Pfannkuch_constant;
            in.clear();
            modSherwood = true;
        }
        //....................................................................
        if (line_string.find("$SHERWOOD_MODEL") !=
            string::npos)  // subkeyword found, SP: keyword for calculation of
                           // diffenrent Sherwood-models
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> shidx;  // Sh_factor >> Re_expo >> Sc_expo >> Pe_expo >>
                          // NContent_expo >> WContent_expo >> D50_expo >>
            // Poro_expo >> Delta_expo >> UI_expo >> NContent_ini_expo >>
            // NContent_res_expo >> GSR_expo >> Tort_expo >> Length_expo >>
            // Pfannkuch_constant;
            if (shidx == 0)  // SP: Standard Sherwood calculation after Bradford
                             // and Abriola, 2001
            {
                in >> Sh_factor >> Re_expo >> Sc_expo;
            }
            if (shidx ==
                1)  // SP: Sherwood calculation after Powers et al., 1994
            {
                in >> Sh_factor >> Re_expo >> Delta_expo >> UI_expo;
                Sherwood_model = true;
            }
            if (shidx ==
                2)  // SP: Sherwood calculation after Powers et al., 1992
            {
                in >> Sh_factor >> Re_expo >> D50_expo >> UI_expo;
                Sherwood_model = true;
            }
            if (shidx == 3)  // SP: Sherwood calculation after Wilson and
                             // Geankoplis, 1966
            {
                in >> Sh_factor;
            }
            if (shidx == 4)  // SP: Sherwood calculation after Pfannkuch, 1984
            {
                in >> Sh_factor >> Pfannkuch_constant;
            }
            if (shidx == 5)  // SP: Sherwood calculation after Miller, 1990
            {
                in >> Sh_factor >> Re_expo >> WContent_expo >> Sc_expo;
                Sherwood_model = true;
            }
            if (shidx ==
                6)  // SP: Sherwood calculation after Saba&Illangasekare, 2000
            {
                in >> Sh_factor >> Re_expo >> Sc_expo;
                Sherwood_model = true;
            }
            if (shidx == 7)  // SP: Sherwood calculation after Geller&Hunt, 1993
            {
                in >> Sh_factor >> Re_expo >> NContent_expo >>
                    NContent_ini_expo >> Poro_expo >> grain_expo;
                Sherwood_model = true;
            }
            in.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$GEOMETRY") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> Geometry_expo;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$GAS_DISSOLUTION") != string::npos)
            gas_dissolution_flag = true;
        // subkeyword found
        if (line_string.find("$INTERFACIAL_AREA") != string::npos)
            while (OK)
            {
                index1 = rfd_file->tellg();
                if (!GetLineFromFile(line, rfd_file))
                    break;
                line_str1 = line;
                if ((line_str1.find(hash) != string::npos) ||
                    (line_str1.find(dollar) != string::npos))
                {
                    OK = false;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste keyword weg
                    break;
                }
                in.str(line_str1);
                //		in >> s_geo_type >> s_geo_name >> d_inivalue ;
                in >> s_geo_type;

                BlobGeoType.push_back(s_geo_type);

                size_t geo_obj_idx(std::numeric_limits<size_t>::max());
                if (s_geo_type.compare("DOMAIN") == 0)
                {
                    BlobGeoName.push_back("domain");
                    in >> d_inivalue;
                    Area_Value.push_back(d_inivalue);
                }
                else
                {
                    in >> s_geo_name >> d_inivalue;
                    BlobGeoName.push_back(s_geo_name);

                    if (s_geo_type.find("POINT") != std::string::npos)
                        // TF 06/2010 - get the point vector and set the
                        // geo_obj_idx
                        if (!((geo_obj.getPointVecObj(unique_name))
                                  ->getElementIDByName(s_geo_name,
                                                       geo_obj_idx)))
                        {
                            std::cerr << "error in CKinBlob::Read: point name "
                                         "not found!"
                                      << "\n";
                            exit(1);
                        }

                    if (s_geo_type.find("POLYLINE") != std::string::npos)
                        // TF 06/2010 - get the polyline vector and set the
                        // geo_obj_idx
                        if (!((geo_obj.getPolylineVecObj(unique_name))
                                  ->getElementIDByName(s_geo_name,
                                                       geo_obj_idx)))
                        {
                            std::cerr << "error in CKinBlob::Read: polyline "
                                         "name not found!"
                                      << "\n";
                            exit(1);
                        }

                    Area_Value.push_back(d_inivalue);
                }
                BlobGeoID.push_back(geo_obj_idx);
                in.clear();
            }
    }  // end while keyword

    return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function for CKinBlob-Structure
   Programing:
   02/2007 DS Implementation
   06/2010 TF changed signature
    (in order to use the Write method for general output streams, for example
std::cout)
**************************************************************************/
void CKinBlob::Write(std::ostream& rfe_file) const
{
    rfe_file << "\n";
    rfe_file << "#BLOB_PROPERTIES"
             << "\n";
    rfe_file << "$NAME"
             << "\n"
             << name << "\n";
    rfe_file << "$D50"
             << "\n"
             << d50 << "\n";
    rfe_file << "$CALC_SHERWOOD"
             << "\n";
    rfe_file << " Sh_factor : " << Sh_factor << "\n";
    rfe_file << " Re_expo   : " << Re_expo << "\n";
    rfe_file << " Sc_expo   : " << Sc_expo << "\n";
    rfe_file << "$GEOMETRY"
             << "\n"
             << Geometry_expo << "\n";
    rfe_file << "$INTERFACIAL_AREA"
             << "\n";
    for (size_t j = 0; j < BlobGeoName.size(); j++)
    {
        rfe_file << " Geotype : " << BlobGeoType[j];
        rfe_file << "   Geoname : " << BlobGeoName[j];
        rfe_file << "   Value   : " << Area_Value[j] << "\n";
    }
    if (modSherwood == true)
    {
        rfe_file << "$CALC_SHERWOOD_MODIFIED "
                 << "\n";
        rfe_file << " UI                   : " << this->UI << "\n";
        rfe_file << " dM                   : " << this->dM << "\n";
        rfe_file << " NAPL_CONTENT_INI     : " << this->NCont_ini << "\n";
        rfe_file << " NAPL_CONTENT_RESIDUAL: " << this->NCont_res << "\n";
        rfe_file << " GRAIN_SPHERE_RATIO   : " << this->GSR << "\n";
        rfe_file << " TORTUOSITY           : " << this->Tort << "\n";
        rfe_file << " LENGTH               : " << this->Length << "\n";

        rfe_file << " Pe_expo           : " << this->Pe_expo << "\n";
        rfe_file << " NContent_expo     : " << this->NContent_expo << "\n";
        rfe_file << " WContent_expo     : " << this->WContent_expo << "\n";
        rfe_file << " d50_expo          : " << this->D50_expo << "\n";
        rfe_file << " n_expo            : " << this->Poro_expo << "\n";
        rfe_file << " delta_expo        : " << this->Delta_expo << "\n";
        rfe_file << " UI_expo           : " << this->UI_expo << "\n";
        rfe_file << " NCont_ini_expo    : " << this->NContent_ini_expo << "\n";
        rfe_file << " NCont_res_expo    : " << this->NContent_res_expo << "\n";
        rfe_file << " GSR_expo          : " << this->GSR_expo << "\n";
        rfe_file << " Tort_expo         : " << this->Tort_expo << "\n";
        rfe_file << " Length_expo       : " << this->Length_expo << "\n";
        rfe_file << " Pfannkuch_constant: " << this->Pfannkuch_constant << "\n";
    }

    rfe_file << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ configure function
   Programing:
   05/2007 DS Implementation
**************************************************************************/
void KBlobConfig(const GEOLIB::GEOObjects& geo_obj,
                 const std::string& unique_name)
{
    std::string s_geo_name, s_geo_type;
    size_t s_geo_id;
    std::vector<long> nodes_vector;

    // create vector for Interfacial_area and initialization with input value
    // (Interfacial_area[0]=input value)
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    if (m_msh == NULL)
    {
        cout << "No mesh in KRConfig"
             << "\n";
        exit(1);
    }

    for (size_t i = 0; i < KinBlob_vector.size(); i++)
    {
        CKinBlob* m_kb(KinBlob_vector[i]);
        for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
        {
            m_kb->Interfacial_area.push_back(-1.);  // Vorbelegung mit Area=-1
            m_kb->Masstransfer_K.push_back(-1.);    // Vorbelegung mit K=-1
        }
        for (size_t j = 0; j < m_kb->BlobGeoName.size(); j++)
        {
            s_geo_name = m_kb->BlobGeoName[j];
            s_geo_type = m_kb->BlobGeoType[j];
            s_geo_id = (m_kb->getBlobGeoID())[j];
            if (s_geo_type.compare("POINT") == 0)
            {
                // 06/2010 TF - switch to new GEOLIB - REMOVE CANDIDATE
                //				CGLPoint* m_geo_point = NULL; // make new GEO
                // point 				m_geo_point =
                // GEOGetPointByName(s_geo_name);//Get GEO
                // point by name 				if (m_geo_point) l =
                // m_msh->GetNODOnPNT(m_geo_point); // + ShiftInNodeVector; //
                // find MSH point
                // number
                // stored
                // in l

                const std::vector<GEOLIB::Point*>* pnt_vec(
                    geo_obj.getPointVec(unique_name));
                size_t msh_node_id =
                    m_msh->GetNODOnPNT((*pnt_vec)[(m_kb->getBlobGeoID())[j]]);

                m_kb->Interfacial_area[msh_node_id] = m_kb->Area_Value[j];
            }  // end if POINT

            if (s_geo_type.compare("POLYLINE") == 0)
            {
                std::vector<GEOLIB::Polyline*> const* const ply_vec(
                    geo_obj.getPolylineVec(unique_name));
                GEOLIB::Polyline const* const polyline((*ply_vec)[s_geo_id]);
                m_msh->GetNODOnPLY(polyline, nodes_vector);
                for (size_t k = 0; k < nodes_vector.size(); k++)
                    m_kb->Interfacial_area[nodes_vector[k]] =
                        m_kb->Area_Value[j];
            }  // end if POLYLINE

            if (s_geo_type.compare("SURFACE") == 0)
            {
                Surface* m_surface = NULL;
                m_surface = GEOGetSFCByName(s_geo_name);
                if (m_surface)
                {
                    m_msh->GetNODOnSFC(
                        m_surface, nodes_vector);  // kg44 if used together with
                                                   // PETSC one maybe needs to
                    // set the for_ic flag...see rf_id_new.cpp
                    for (size_t k = 0; k < nodes_vector.size(); k++)
                    {
                        //+ShiftInNodeVector;
                        size_t msh_node_id = nodes_vector[k];
                        m_kb->Interfacial_area[msh_node_id] =
                            m_kb->Area_Value[j];
                    }
                }
            }  // end if SURFACE

            if (s_geo_type.compare("DOMAIN") == 0)
                for (size_t k = 0; k < m_msh->nod_vector.size(); k++)
                    m_kb->Interfacial_area[k] = m_kb->Area_Value[j];
            // end if SURFACE
        }
    }  // end loop over i = Number of blobs

    nodes_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task: Check OBJ configuration for input errors
   Programing:
   05/2007 DS Implementation
**************************************************************************/
void KBlobCheck(void)
{
    int i;
    long k, l;

    CKinBlob* m_kb = NULL;
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt

    cout << "\n"
         << "Checking defined blob classes"
         << "\n";
    for (i = 0; i < (int)KinBlob_vector.size(); i++)
    {
        m_kb = KinBlob_vector[i];
        cout << m_kb->name << "\n";

        if (m_kb->d50 <= 0.)
        {
            cout << "  Warning: D50 <= 0., D50 is set to 1.E-4"
                 << "\n";
            m_kb->d50 = 1.E-4;
        }
        if (m_kb->Sh_factor < 0.)
        {
            cout << "  Warning: Sh_factor < 0., Sh_factor is set to 1.15"
                 << "\n";
            m_kb->Sh_factor = 1.15;
        }
        if (m_kb->Re_expo < 0.)
        {
            cout << "  Warning: Re_expo < 0., Re_expo is set to 0.654"
                 << "\n";
            m_kb->Re_expo = 0.654;
        }
        if (m_kb->Sc_expo < 0.)
        {
            cout << "  Warning: Sc_expo < 0., Sc_expo is set to 0.486"
                 << "\n";
            m_kb->Sc_expo = 0.486;
        }
        if (m_kb->Geometry_expo < 0.)
        {
            cout
                << "  Warning: Geometry_expo < 0., Geometry_expo is set to 0.66"
                << "\n";
            m_kb->Geometry_expo = 0.66;
        }

        k = 0;
        for (l = 0; l < (long)m_msh->nod_vector.size(); l++)
        {
            if (m_kb->Interfacial_area[l] < 0.)
            {
                cout << "  Warning: no value for node " << l
                     << " , Interfacial area is set to zero."
                     << "\n";
                m_kb->Interfacial_area[l] = 0.;
            }
            else
                k += 1;
        }
        cout << "  Values for interfacial Surfaces have been defined for " << k
             << " nodes by user."
             << "\n";
    }  // end loop i over blob classes
}

/**************************************************************************
   Reaction-Method:
   Task: Class constructor
   Programing:
   02/2006 SB Implementation
**************************************************************************/
CKinReactData::CKinReactData(void)
{
    SolverType = 1;
    relErrorTolerance = 1.0e-10;
    minTimestep = 1.0e-10;
    initialTimestep = 1.0e-10;
    NumberReactions = 0;
    NumberFreundlich = 0;
    NumberLangmuir = 0;
    NumberLinear = 0;
    NumberMonod = 0;
    NumberNAPLdissolution = 0;
    NumberMineralkinetics = 0;
    NumberMicrobeData = 0;
    usedt = -1.0;
    maxBacteriaCapacity = -1.0;
    minBactConcentration = 1e-30;
    minConcentrationMode = 0;
    minConcentrationThreshhold = 1e-30;
    minConcentrationSet = 1e-30;
    is_a_bacterium.clear();
    testoutput = false;  // true;
    maxSurfaces = 3;     // Why is max = 3 ?? CB
    exSurface.clear();
    // sp_index.clear();
    // kr_active_species = -1;
    sp_varind.clear();
    // sp_pcsind.clear(); //HS
    // das Surface-Array k�nnte auch dynamisch allokiert werden
    for (int i = 0; i < maxSurfaces; i++)
        exSurface.push_back(-1.0);
    NoReactGeoName.clear();
    NoReactGeoType.clear();
    NoReactGeoID.clear();
    is_a_CCBC.clear();
    AllowReactGeoName.clear();
    AllowReactGeoType.clear();
    AllowReactGeoID.clear();
    node_foc.clear();
    cpu_time_krc = 0;  // to measure cpu time

    IonicStrengths.clear();
    ActivityCoefficients.clear();
    activity_model = 0;

    ReactDeactMode = -1;
    ReactDeactEpsilon = -1.0;  // CB ReactDeact
    ReactDeactCThresh = -1.0;
    ReactDeactRelative = false;
    ReactDeactFlag = false;
    ReactDeactPlotFlag = 0;
    ReactDeact.clear();  // flags for individual nodes
    React_dCdT.clear();  // Sum of reaction rates for individual nodes
    ReactNeighborhood
        .clear();  // node indices of local neighborhood around individual nodes
    copy_concentrations = false;  // SB 01.2011
    copy_nodes.clear();           // SB 01.2011
    radial = false;               // default case
    batch = false;
    TwoDinThreeD = false;

    debugoutflag = false;
    scale_dcdt = false;
    lagneau = false;
    sortnodes = false;
    OmegaThresh = 1e-10;

    noototnodes = 0;
    noocalcnodes = 0;
}

/**************************************************************************
   Reaction-Method:
   Task: Class destructor
   Programing:
   02/2006 SB Implementation
**************************************************************************/
CKinReactData::~CKinReactData(void)
{
    if (sortnodes)
    {
        delete[] substeps;
        delete[] noderanks;
    }
}

/**************************************************************************
   FEMLib-Method:
   Task: Reaction class data read function
   Programing:
   02/2006 SB New Class structure, IO, C++, FEM
**************************************************************************/
bool CKinReactData::Read(ifstream* rfd_file, const GEOLIB::GEOObjects& geo_obj,
                         const std::string& unique_name)
{
    char line[MAX_ZEILE];
    string sub_line;
    string line_string, line_str1, help("1");
    string delimiter(" ");
    bool new_keyword = false, OK = true;
    string hash("#"), dollar("$");
    std::stringstream in;
    size_t index, index1;
    int /* count_surf, */ surf_id;
    string s_geo_type, s_geo_name;
    string dummy;
    bool allow = false;
    bool allow_not = false;

    //========================================================================
    while (!new_keyword)
    {
        index = rfd_file->tellg();
        if (!GetLineFromFile(line, rfd_file))
            break;
        line_string = line;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            rfd_file->seekg(index);  // Dateipointer zur�cksetzen, sonst ist das
                                     // n�chste keyword weg
            break;
        }
        /* Keywords nacheinander durchsuchen */
        //....................................................................
        // subkeyword found
        if (line_string.find("$SOLVER_TYPE") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> SolverType;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$RELATIVE_ERROR") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> relErrorTolerance;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$MIN_TIMESTEP") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> minTimestep;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$INITIAL_TIMESTEP") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> initialTimestep;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$BACTERIACAPACITY") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> maxBacteriaCapacity;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$MIN_BACTERIACONC") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> minBactConcentration;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$MIN_CONCENTRATION_REPLACE") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> minConcentrationMode >> minConcentrationThreshhold >>
                minConcentrationSet;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$SURFACES") != string::npos)
            while (OK)
            {
                index1 = rfd_file->tellg();
                if (!GetLineFromFile(line, rfd_file))
                    break;
                line_str1 = line;
                if ((line_str1.find(hash) != string::npos) ||
                    (line_str1.find(dollar) != string::npos))
                {
                    OK = false;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste keyword weg
                    break;
                }
                in.str(line_str1);
                in >> surf_id;
                if (surf_id > maxSurfaces - 1)
                {
                    cout
                        << " Warning in CKinReactData::Read: Specified surf_id "
                        << surf_id << " is equal"
                        << "\n";
                    cout << " or larger than max no. of surfaces maxSurfaces = "
                         << maxSurfaces << "."
                         << "\n";
                    cout << " First surf_id specified in KRC input file should "
                            "be 0."
                         << "\n";
                }
                in >> exSurface[surf_id];
                in.clear();
            }
        //....................................................................
        // subkeyword found
        allow = allow_not = false;
        if (line_string.find("$ALLOW_REACTIONS") != string::npos)
            allow = true;
        if (line_string.find("$NO_REACTIONS") != string::npos)
            allow_not = true;
        if ((allow == true) || (allow_not == true))
        {
            while (OK)
            {
                index1 = rfd_file->tellg();
                if (!GetLineFromFile(line, rfd_file))
                    break;
                line_str1 = line;
                if ((line_str1.find(hash) != string::npos) ||
                    (line_str1.find(dollar) != string::npos))
                {
                    OK = false;
                    rfd_file->seekg(
                        index1);  // Dateipointer zur�cksetzen, sonst ist das
                                  // n�chste keyword weg
                    break;
                }
                in.str(line_str1);
                in >> s_geo_type >> s_geo_name;
                if (allow_not == true)
                {
                    NoReactGeoType.push_back(s_geo_type);
                    NoReactGeoName.push_back(s_geo_name);
                }
                else if (allow == true)
                {
                    AllowReactGeoType.push_back(s_geo_type);
                    AllowReactGeoName.push_back(s_geo_name);
                }

                // TF 06/2010 - for the change string-ID to size_t-ID
                size_t geo_obj_idx(std::numeric_limits<size_t>::max());

                if (s_geo_type.find("POINT") != std::string::npos)
                {
                    // get the point vector and set the geo_obj_idx
                    if (!((geo_obj.getPointVecObj(unique_name))
                              ->getElementIDByName(s_geo_name, geo_obj_idx)))
                    {
                        std::cerr << "error in CKinReactData::Read: (type="
                                  << s_geo_type << "): " << s_geo_name
                                  << " point name not found!"
                                  << "\n";
                        exit(1);
                    }
                }
                if (s_geo_type.find("POLYLINE") != std::string::npos)
                {
                    // get the point vector and set the geo_obj_idx
                    if (!((geo_obj.getPolylineVecObj(unique_name))
                              ->getElementIDByName(s_geo_name, geo_obj_idx)))
                    {
                        std::cerr
                            << "error in CKinReactData::Read: polyline name "
                            << s_geo_name << " not found!"
                            << "\n";
                        exit(1);
                    }
                }

                //
                if (allow_not == true)
                    NoReactGeoID.push_back(geo_obj_idx);
                else if (allow == true)
                    AllowReactGeoID.push_back(geo_obj_idx);

                in.clear();
            }
        }

        //....................................................................
        if (line_string.find("$COPY_CONCENTRATIONS") != string::npos)
        {  // SB 01.2011
            std::cout << "COPY_CONCENTRATIONS allows copying concentrations "
                         "after reactions to neighbouring"
                      << " nodes for symmetric models:"
                      << "\n";
            std::cout << " - BATCH: 2 node quasi 0D models."
                      << "\n";
            std::cout << " - RADIAL symmetric quasi 2D hex models with y = 0 "
                         "as symmetry axis."
                      << "\n";
            std::cout << " - 2DIN3D: vertical x-z cross sections in full 3D "
                         "hex models with arbitrary y as symmetry axis."
                      << "\n";

            this->copy_concentrations = true;
            in.str(GetLineFromFile1(rfd_file));
            in >> dummy;
            in.clear();
            if (dummy.compare("RADIAL") == 0)
                radial = true;  // default case
            else if (dummy.compare("BATCH") == 0)
                batch = true;
            else if (dummy.compare("2DIN3D") == 0)
                TwoDinThreeD = true;
            else
            {
                cout << "Error in CKinReactData::Read: $COPY_CONCENTRATIONS: "
                        "No copy mode defined."
                     << "\n"
                     << " Define either RADIAL or BATCH or 2DIN3D. Terminating "
                        "now!"
                     << "\n";
                exit(1);
            }
        }
        //....................................................................
        if (line_string.find("$LAGNEAU_BENCHMARK") != string::npos)
        {  // CB 04.2011
            this->lagneau = true;
        }
        //....................................................................
        if (line_string.find("$SCALE_DCDT") != string::npos)
        {                             // CB 04.2011
            this->scale_dcdt = true;  // scale dcdt vector in derivs for
                                      // stabilization of ODE solver
        }
        //....................................................................
        if (line_string.find("$SORT_NODES") != string::npos)
        {                            // CB 04.2011
            this->sortnodes = true;  // scale dcdt vector in derivs for
                                     // stabilization of ODE solver
        }
        //.....................................................................
        if (line_string.find("$OMEGA_THRESHOLD") != string::npos)
        {  // CB 11.2011
            in.str(GetLineFromFile1(rfd_file));
            in >> OmegaThresh;  // set Minrate to 0 if |1-omega|< OmegaThresh
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$REACTION_DEACTIVATION") != string::npos)
        {
            in.str(GetLineFromFile1(rfd_file));
            in >> ReactDeactMode >> ReactDeactEpsilon >> ReactDeactPlotFlag >>
                dummy >> ReactDeactCThresh;
            in.clear();
            if (dummy.compare("RELATIVE") == 0)
                ReactDeactRelative = true;
            else
                ReactDeactRelative = false;
        }
        // subkeyword found
        if (line_string.find("$DEBUG_OUTPUT") != string::npos)
            debugoutflag = true;
        //....................................................................
        if (line_string.find("$ACTIVITY_MODEL") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(rfd_file));
            in >> activity_model;
            in.clear();
        }
        //....................................................................
    }
    usedt = DMAX(initialTimestep, minTimestep);
    return true;
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction data class write function
   Programing:
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReactData::Write(ofstream* rfe_file)
{
    int i;
    // Write Keyword
    *rfe_file << "\n";
    *rfe_file << "#KINREACTIONDATA"
              << "\n";
    *rfe_file << "$SOLVER_TYPE"
              << "\n"
              << SolverType << "\n";
    *rfe_file << "$RELATIVE_ERROR"
              << "\n"
              << relErrorTolerance << "\n";
    *rfe_file << "$MIN_TIMESTEP"
              << "\n"
              << minTimestep << "\n";
    *rfe_file << "$INITIAL_TIMESTEP"
              << "\n"
              << initialTimestep << "\n";
    *rfe_file << "$BACTERIACAPACITY"
              << "\n"
              << maxBacteriaCapacity << "\n";
    //*rfe_file << " max Surfaces : " << maxSurfaces << "\n";
    *rfe_file << "$SURFACES"
              << "\n";  //<< (int)exSurface.size() << "\n";
    for (i = 0; i < (int)exSurface.size(); i++)
        *rfe_file << i + 1 << "  " << exSurface[i] << "\n";
    *rfe_file << "NO_REACTIONS"
              << "\n";  // << (int)NoReactGeoName.size() << "\n";
    for (i = 0; i < (int)NoReactGeoName.size(); i++)
        *rfe_file << NoReactGeoType[i] << "  " << NoReactGeoName[i] << "\n";
    *rfe_file << "$REACTION_DEACTIVATION	"
              << "\n"
              << ReactDeactMode << " " << ReactDeactEpsilon << " "
              << ReactDeactPlotFlag << "\n";
    *rfe_file << "$DEBUG_OUTPUT	"
              << "\n"
              << debugoutflag << "\n";
    //*rfe_file << " Number of reactions: " << NumberReactions << "\n";
    //*rfe_file << " Number of linear exchange reactions: " << NumberLinear <<
    //"\n"; *rfe_file << " Number of freundlich exchange reactions: " <<
    // NumberFreundlich << "\n"; *rfe_file << " Number of langmuir exchange
    // reactions: " << NumberLangmuir << "\n"; *rfe_file << " is_a_bacterium: "
    //<< "\n";
    // for(i=0;i<is_a_bacterium.size();i++) *rfe_file << is_a_bacterium[i] << "
    // ";
    //*rfe_file << "\n";
    *rfe_file << "\n";
    /*
       *rfe_file << " usedt "<< usedt << "\n";
       *rfe_file << " Number Reactions "<< NumberReactions << "\n";
       *rfe_file << " is_a_bacterium " << "\n";
       for(int i=0; i < (int)is_a_bacterium.size();i++)
       *rfe_file <<  is_a_bacterium[i];
       *rfe_file << "\n";
     */
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction data class write function
   Programing:
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReactData::TestWrite(void)
{
    int i;  // CB
    // Write Keyword
    cout << "#KINREACTIONDATA"
         << "\n";
    cout << "$SOLVER_TYPE"
         << "\n"
         << SolverType << "\n";
    cout << "$RELATIVE_ERROR"
         << "\n"
         << relErrorTolerance << "\n";
    cout << "$MIN_TIMESTEP"
         << "\n"
         << minTimestep << "\n";
    cout << "$INITIAL_TIMESTEP"
         << "\n"
         << initialTimestep << "\n";
    cout << "$BACTERIACAPACITY"
         << "\n"
         << maxBacteriaCapacity << "\n";
    cout << "$REACTION_DEACTIVATION	"
         << "\n"
         << ReactDeactMode << " " << ReactDeactEpsilon << " "
         << ReactDeactPlotFlag << "\n";
    cout << "$DEBUG_OUTPUT	"
         << "\n"
         << debugoutflag << "\n";

    cout << "\n";
    cout << " usedt " << usedt << "\n";
    cout << " Number Reactions " << NumberReactions << "\n";
    cout << " is_a_bacterium " << (int)is_a_bacterium.size() << "\n";
    for (i = 0; i < (int)is_a_bacterium.size(); i++)
        cout << is_a_bacterium[i] << " ";
    cout << "\n";
    cout << " Exchange reactions : "
         << "\n";
    // ACHTUNG: die folgende Zeile wird ausgegeben BEVOR NumberXXX berechnet
    // wurde !
    cout << " Linear exchange : " << NumberLinear << "\n"
         << " Freundlich exchange : " << NumberFreundlich << "\n"
         << " Langmuir exchange : " << NumberLangmuir << "\n"
         << " NAPL dissolution : " << NumberNAPLdissolution << "\n";
    for (i = 0; i < (int)exSurface.size(); i++)
        cout << " " << exSurface[i];
    cout << "\n";
}

/**************************************************************************
   Reaction-Method:
   Task: KinReaction write function - echo of input values to rfe - file
   Programing:
   05/2004 SB Implementation
   02/2006 SB C++, IO
**************************************************************************/
bool KRWrite(std::string const& prot_name)
{
    CKinReact* m_kr = NULL;
    CKinReactData* m_krd = NULL;
    CKinBlob* m_kp = NULL;
    string rfe_file_name;

    if (KinReact_vector.empty() && KinBlob_vector.empty() &&
        KinReactData_vector.empty())
        return true;  // NW skip writing a file if no KinReaction data are
                      // loaded

    //========================================================================
    // File handling
    rfe_file_name = prot_name + "_echo";
    ofstream rfe_file(rfe_file_name.data(), ios::app);
    if (!rfe_file.good())
        return false;
    rfe_file.seekp(0L, ios::end);  // go to end
    //========================================================================
    rfe_file
        << "\n"
        << "; Reactions "
           "----------------------------------------------------------------- "
        << "\n";
    // Output all Reactions
    size_t length(KinReact_vector.size());
    for (size_t i = 0; i < length; i++)
    {
        m_kr = KinReact_vector[i];
        m_kr->Write(&rfe_file);
        if (KinReactData_vector[0]->testoutput)
            m_kr->TestWrite();
    }
    // Output all BlobProperties
    length = KinBlob_vector.size();
    for (size_t i = 0; i < length; i++)
    {
        m_kp = KinBlob_vector[i];
        m_kp->Write(rfe_file);
        m_kp->Write(std::cout);
    }
    // Output all kinetic reaction data
    if (!KinReactData_vector.empty())
        m_krd = KinReactData_vector[0];
    if (m_krd != NULL)
    {
        m_krd->Write(&rfe_file);
        if (KinReactData_vector[0]->testoutput)
            m_krd->TestWrite();
    }
    rfe_file.close();

    return true;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ExecuteKinReact()                                 */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Haupt-Subroutine zur Berechnung der mikrobiellen kinetischen Reaktionen*/
/* Aufruf in void ExecuteReactions(void)  (rf_react.cpp)                  */
/*                                                                        */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ExecuteKinReact(void)
{
    const size_t nnodes(
        fem_msh_vector[0]->nod_vector.size());  // SB: ToDo hart gesetzt

    // CB Reaction deactivation for this time step
    if (ReactDeactFlag)
    {
        if (aktueller_zeitschritt > 2)
            ReactionDeactivation(
                nnodes);  // Check if nodes should be deactivated or activated
                          // for this time step
        ClockTimeVec[0]->StopTime("ReactDeact");  // CB time
        ClockTimeVec[0]->StartTime();             // CB time
    }
    if (debugoutflag)
    {
        debugoutstr.setf(ios::scientific, ios::floatfield);
        debugoutstr.precision(12);
        debugoutstr.open(debugoutfilename.c_str());
    }
    if (time_vector.size() > 0)
        dt = time_vector[0]->CalcTimeStep();
    // This function  prepares
    // - flow velocities of PS_GLOBAL
    // at all nodes
    if (this->NumberNAPLdissolution > 0)
        NAPLDissolutionPreprocessing();
    // This function  prepares
    //  - Ionic streghts
    //  - activity coefficients
    // at all nodes and for all species
    if (this->NumberMineralkinetics > 0)
        PreprocessMinKin();
    if (this->NumberMicrobeData > 0)
        this->PreprocessMicrobe_drmc_(dt);

    if ((dt > 1.E-20) && (aktueller_zeitschritt > 0))
    {
        size_t save_node(0);
        int save_nok = 0, save_nbad = 0;
        double usedttmp = 1.E+30;
        /* Einstellungen Gleichungsl�ser f�r alle Knoten gleich */
        /* relative Genauigkeit des Gleichungsl�sers (eps< fehler/c0) */
        double eps = relErrorTolerance;
        /* min zulaessiger Zeitschritt*/
        double hmin = minTimestep;

        double usedtneu = 0.0;
        long node;
        size_t it, node_idx = 0;
        int nok = 0, nbad = 0, totsteps = 0;
        double* Concentration;  // HS: concentration of all result
        int nComponents;
        bool dry = false;
        // CFEMesh* m_msh = fem_msh_vector[0];            //SB: ToDo hart
        // gesetzt CTimeDiscretization *m_tim = NULL;
        int tt = 0, tsteps = 1;
        // cout << " ExecuteKineticReactions" <<  "\n";

        nComponents = (int)cp_vec.size();
        nComponents++;

        Concentration = new double[nComponents * nnodes];
        for (it = 0; it < (nComponents)*nnodes; it++)
            Concentration[it] = 0.0;
        if (this->sortnodes)
        {
            for (it = 0; it < nnodes; it++)
                substeps[it] = 0;
        }
#if defined(USE_MPI) || defined(USE_MPI_KRC)
        if (this->NumberReactions > 0)
            cout << " ExecuteKineticReactions MPI"
                 << "\n";
        double* Concentration_buff;
        Concentration_buff = new double[nComponents * nnodes];
        for (it = 0; it < (long)(nComponents)*nnodes; it++)
            Concentration_buff[it] = 0.0;
#else
        if (this->NumberReactions > 0)
            cout << " ExecuteKineticReactions"
                 << "\n";
#endif

        // Get the reaction interface data for checking for dried out nodes from
        // eclipse coupling --> rateflag = 0
        REACTINT* m_rei = NULL;
        if (REACTINT_vec.size() > 0)
            m_rei = REACTINT_vec[0];

        if (NumberReactions > 0)
        {
            unsigned count = 0;
#if defined(USE_MPI) || defined(USE_MPI_KRC)
            // MPI initialization.
            // So here is going to distribute the task.
            long nNodes = (long)nnodes;
            MPI_Bcast(&nNodes, 1, MPI_LONG, 0, comm_DDC);
            // here "myrank" is the index of the CPU Processes, and "mysize" is
            // the number of CPU Processes
            for (it = myrank; it < nnodes; it += mysize)
            {
#else
            for (it = 0; it < nnodes; it++)
            {
#endif
                //      substepping
                for (tt = 0; tt < tsteps; tt++)
                {
                    if (sortnodes)  // execute reactions at nodes in sequence
                                    // according to no of substeps in last time
                        // step
                        node = noderanks[it];
                    else
                        node = it;
                    // cout << node << "\n";

                    // check for dryout
                    dry = false;  // by default, all nodes are wet
                    if (m_rei)
                    {
                        if (m_rei->s_water_limit)
                            if (m_rei->dried_out_nodes[node])
                                dry = true;
                    }
                    // no reactions at Concentration BCs
                    if (is_a_CCBC[node] == true)
                    {
                    }
                    // CB no reactions at deactivated nodes
                    else if ((ReactDeactFlag) && (ReactDeact[node] == true))
                    {
                        noototnodes++;
                    }
                    // dryout node
                    else if (dry)
                    {
                    }
                    else
                    {
#if defined(USE_MPI) || defined(USE_MPI_KRC)
                        Biodegradation(Concentration_buff + node * nComponents,
                                       node, eps, hmin, &usedtneu, &nok, &nbad,
                                       tt, tsteps);
#else
#ifdef OGS_FEM_CAP
                        if (tt > 0)
                        {  // update equilibrium system by call to chemapp
                            m_rei->ReactionPreProcessing();
                            REACT_CAP_vec[0]->ExecuteReactionsChemApp(
                                -1, it);  // liquid
                        }
#endif
                        Biodegradation(Concentration + node * nComponents, node,
                                       eps, hmin, &usedtneu, &nok, &nbad, tt,
                                       tsteps);

#ifdef OGS_FEM_CAP
                        if (tt < tsteps - 1)
                        {  // update equilibrium system by call to chemapp
                            REACT_CAP_vec[0]->ExecuteReactionsChemApp(
                                1, it);  // solid+liquid
                        }
#endif

#endif
                        if (usedtneu < usedttmp)
                        {
                            usedttmp = usedtneu;
                            save_nok = nok;
                            save_nbad = nbad;
                            save_node = node;
                        }
                        count++;
                        totsteps = totsteps + nbad + nok;
                        noocalcnodes++;
                        noototnodes++;

                        if (sortnodes)
                            substeps[node] =
                                nbad + nok;  // no of substeps as indicator for
                                             // computation time

                        if (tsteps > 1)
                            usedt =
                                usedtneu;  // use last successfull timestep to
                                           // continue integration over substeps
                    }
                }  // tt
                   //

            }  // end for(node...

// cout << " Total no. of RK-Steps: " << totsteps << "\n";

// Get MPI results
#if defined(USE_MPI) || defined(USE_MPI_KRC)
            MPI_Allreduce(Concentration_buff, Concentration,
                          nnodes * nComponents, MPI_DOUBLE, MPI_SUM, comm_DDC);
#endif
            // post-processing
            // PostprocessNAPLDissolution(node);
            /**/
            // update results
            for (node_idx = 0; node_idx < nnodes; node_idx++)
            {
                // only update those nodes, where Biodegradation() was executed

                dry = false;  // by default, all nodes are wet
                if (m_rei)
                {
                    if (m_rei->s_water_limit)
                        if (m_rei->dried_out_nodes[node_idx])
                            dry = true;
                }

                if (is_a_CCBC[node_idx] == true)
                {
                }
                else if ((ReactDeactFlag) && (ReactDeact[node_idx] == true))
                {
                }
                else if (dry)
                {
                }
                else
                {
                    for (int sp = 0; sp < (nComponents - 1); sp++)
                    {
                        // Notlösung gegen das vollständige Absterben der
                        // Bakterien
                        if ((is_a_bacterium[sp]) &&
                            (Concentration[node_idx * nComponents + sp + 1] <
                             minBactConcentration))
                            Concentration[node_idx * nComponents + sp + 1] =
                                minBactConcentration;
                        // Konzentrationen aller Substanzen in Datenstruktur
                        // zurückschreiben
                        // pcs_vector[sp_pcsind[sp]]->SetNodeValue(node_idx,sp_varind[sp],Concentration[node_idx*nComponents+sp+1]);
                        // // CB HS update
                        cp_vec[sp]->getProcess()->SetNodeValue(
                            node_idx, sp_varind[sp],
                            Concentration[node_idx * nComponents + sp + 1]);
                        // save exchange term SB todo
                    }
                }
            }  // end of concentration update
            cout << "    Kinetic reactions executed at " << count << " of "
                 << nnodes << " nodes."
                 << "\n";
            cout << "    Total no of nodes with reactions so far:            "
                 << noototnodes << "\n";
            cout << "    Total no of nodes with reactions calculated so far: "
                 << noocalcnodes << "\n";

            // sort noderanks according to computation time
            if (sortnodes)
                SortIterations(substeps, noderanks, nnodes);
            // for(node_idx = 0; node_idx < nnodes; node_idx++)
            //  cout << substeps[node_idx] << " " << noderanks[node_idx] <<
            //  "\n";

            // CB Reaction deactivation for next time step
            if (ReactDeactFlag)
            {
                // cout << "    Kinetic reactions executed at " << count << " of
                // " << nnodes << " nodes." << "\n";
                if (ReactDeactMode ==
                    2)  // For mode 2 the C_new must be updated by C_old (C
                        // after eraction of last time step)
                    ReactDeactSetOldReactionTerms(nnodes);
                // Aromaticum(nnodes);
            }
        }  // end if(NumberRactions>0)

        if (usedttmp < usedt)
        {
            cout << "\n"
                 << "Kinetics in node " << save_node
                 << " limit integration step - nok: ";
            cout << save_nok << " nbad: " << save_nbad << "\n";
        }

        // update des zul�ssigen Integrationsschritts, verwendet beim Aufruf von
        // odeint kleinster Wert, der in einem der Knoten zu einer zuverl�ssigen
        // Integration gef�hrt hat konservative, aber stabile Annahme
        usedttmp = DMAX(usedttmp, hmin);
        usedt = DMIN(usedttmp, dt);
        // cout << "\n" << " Next suggested integration step " << usedt << "\n";

#if defined(USE_MPI) || defined(USE_MPI_KRC)
        // free_dvector(Concentration_buff,(long)1,(long)nComponents*nnodes);
        delete[] Concentration_buff;
#endif
        // free_dvector(Concentration,(long)1,(long)nComponents*nnodes);
        delete[] Concentration;

    }  // end if((dt>1.E-20)&&(aktueller_zeitschritt>0)){

    if (debugoutflag)
        debugoutstr.close();
    if (ReactDeactFlag)
        if (ReactDeactPlotFlag == 1)
            ReactDeactPlotFlagsToTec();
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Biodegradation(node, )                            */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet Bioabbau im Knoten node in dt                                */
/* Tr�gt berechnete Senke in ?? ein                                       */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: node (Knotennummer), tstart (Zeitpunkt Rechenbeginn)                */
/*    tend (Zeitpunkt Rechenende)                                         */
/*                                                                        */
/* Ergebnis:                                                              */
/* - void -                                                               */
/* Ergebnisse werden in Datenstruktur zur�ckgeschrieben                   */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/* 05/2007     DS         NAPL-dissolution added                          */
/*                                                                        */
/**************************************************************************/

void CKinReactData::Biodegradation(double* m_Conc, long node, double eps,
                                   double hmin, double* usedtneu, int* nok,
                                   int* nbad, int tt, int tsteps)
{
    // double *Concentration;
    double* newVolume = NULL;
    double* oldVolume = NULL;
    double* oldMass = NULL;
    double nexth = 0.;
    long sp;  //, timelevel;
    //  int nok=0, nbad=0, Number_of_Components;
    int Number_of_Components, nreactions, r, Sp1, blob, Number_of_blobs;
    double Csat_max = 0, DensityNAPL = 0, DensityAQ = 0, DiffusionAQ = 0,
           ViscosityAQ = 0, PoreVelocity = 0.0, d50 = 0, Reynolds = 0,
           Schmidt = 0, Sherwood = 0, NCont_Sh2 = 0;  // NCont_Sh1=0,
    double DarcyVelocity = 0.0;
    double mod_Reynolds;
    double beta_4, grain_var;  //, grain_expo;
    double one_three, two_three;
    double Peclet, Delta;
    double NAPLcontent = 1;
    double WATERcontent = 1;
    double Sat, Poro;
    double tstart, tend, dt = 0.0;
    double baditerations;
    CRFProcess* m_pcs = NULL;
    int idxS = 0;
    int shidx;
    CTimeDiscretization* m_tim = NULL;
    float garage = 0;
    bool success = false;
    // double *m_Conc_save;  // CB: for failure of ODE integration
    // this is a temp array used to feed C as float to odeint in case
    // integration failed m_Conc_save = new double[Number_of_Components];
    vector<double> tempstore;
    double variables[3];
    variables[0] = variables[1] = variables[2] = 0;
    bool Tflag = false;
    // double A = 0; //Daq(T) PCE
    // double B = 0; //Daq(T) PCE
    bool gassdiss = false;
    double Cmin = 0;

    CKinReact* m_kr = NULL;
    CKinBlob* m_kb = NULL;
    CKinReactData* m_krd = NULL;
    m_krd = KinReactData_vector[0];
    nreactions = m_krd->NumberReactions;
    // timelevel = 1;                                 // concentrations are in
    // new timelevel

    if (time_vector.size() > 0)
    {
        m_tim = time_vector[0];
        dt = m_tim->CalcTimeStep();
    }
    // Number_of_Components = kr_active_species; //
    Number_of_Components = (int)cp_vec.size();
    // Get storage for vector of concentrations
    // Concentration = dvector(1,Number_of_Components);

    /* Konzentrationen aller Substanzen aus Datenstruktur auslesen und in neuem
     * Array speichern */
    for (sp = 0; sp < Number_of_Components; sp++)
    {
        m_Conc[sp + 1] =
            cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
        // CB reset concentrations below a threshhold to a minimum value
        // 1e-19 and 0.0 by default
        if (minConcentrationMode > 0)
        {
            if (minConcentrationMode == 1)
                Cmin = fabs(m_Conc[sp + 1]);
            else if (minConcentrationMode == 2)
                Cmin = m_Conc[sp + 1];
            if (Cmin < minConcentrationThreshhold)
                m_Conc[sp + 1] = minConcentrationSet;
        }
        // SB todo - ist abewr gerade eh noch ein dummy
        // ExchangeTerm[sp]=TBCGetExchange(node,sp)/(dt);
    }

    // CB CLEAN speed up fix
    // if(m_Conc[12] < 5 ) return; // CO2 concentration
    if (debugoutflag)
    {
        if (node == 1)
        {
            debugoutstr << "node timestep " << flush;
            // debugoutstr << " Concentrations before odeint: " << "\n" << " "
            // << flush;
            for (sp = 0; sp < Number_of_Components; sp++)
                debugoutstr << cp_vec[sp]->getProcess()->nod_val_name_vector[0]
                            << " " << flush;
            debugoutstr << " SATURATION1 " << flush;

            // debugoutstr << pcs_vector[sp_pcsind[sp]]->nod_val_name_vector[0]
            // << " " << flush;  // CB HS update
            debugoutstr << "usedt nok nbad"
                        << "\n";
        }
        debugoutstr << node << " " << aktueller_zeitschritt << " " << flush;
        for (sp = 0; sp < Number_of_Components; sp++)
            debugoutstr << m_Conc[sp + 1] << " " << flush;
        debugoutstr << KinReact_vector[0]->GetSaturation(node, 1, 1) << " "
                    << flush;
        debugoutstr << usedt << " - -"
                    << "\n"
                    << flush;
    }
    //#ds
    /* PREPARE PARAMETERS FOR NAPL-DISSOLUION*/
    /* calculate Mass and Volume of blobs for this node */
    // first get T
    for (r = 0; r < (int)KinReact_vector.size(); r++)
    {
        if (KinReact_vector[r]->T_dependence == true &&
            KinReact_vector[r]->typeflag_napldissolution == 1)
        {
            Tflag = true;
            if (REACTINT_vec.size() > 0)  // Get the Temperature
                variables[1] = REACTINT_vec[0]->GetTemperature(node);
            break;
        }
    }

    if (Tflag ==
        false)  // CB, SP, 19.04.2013, get Temp. for DensityAQ calculations
        if ((mfp_vector[0]->getCompressibilityTModel() != 0) ||
            (mfp_vector[0]->diffusion_model == 2))
            if (REACTINT_vec.size() > 0)  // Get the Temperature
                variables[1] = REACTINT_vec[0]->GetTemperature(node);

    // set all mass(blob)=0, volume(blob)=0
    Number_of_blobs = (int)KinBlob_vector.size();
    if (Number_of_blobs > 0)
    {
        for (r = 0; r < Number_of_blobs; r++)
        {
            if (KinBlob_vector[r]->gas_dissolution_flag)
            {
                gassdiss = true;
                break;
            }
        }

        if (gassdiss == false)
            DensityAQ = mfp_vector[0]->Density(
                variables);  // should be returned as temperature dependent
                             // value, if required
        if (Tflag)
        {
            if (REACTINT_vec.size() > 0)  // Get the Temperature
                if (REACTINT_vec[0]->constanttemperature)
                {
                    variables[1] = REACTINT_vec[0]->GetTemperature(node);
                    DensityAQ =
                        REACTINT_vec[0]->LiquidDensity_Busch(variables[1]);
                }
        }

        // get storage for local data structures
        oldVolume = dvector(0, Number_of_blobs);
        oldMass = dvector(0, Number_of_blobs);
        newVolume = dvector(0, Number_of_blobs);
        // initialize
        for (r = 0; r < Number_of_blobs; r++)
            oldVolume[r] = oldMass[r] = newVolume[r] = 0.;
        // 1) Here, calculate current Mass and old Volume of NAPL for each blob
        //    - Mass is required for computing current Csat in 2)
        //    - old Volume is required for updating Interfacial Area in
        //    postprocessing
        if (KinBlob_vector[0]->gas_dissolution_flag == false)
        {
            for (r = 0; r < nreactions; r++)
            {
                m_kr = KinReact_vector[r];
                // CB new reaction switch for individual reactions
                if (m_kr->switched_off_node.size() > 0)
                    if (m_kr->switched_off_node[node] == true)
                        continue;
                if (m_kr->typeflag_napldissolution)  // dissolution reaction
                                                     // identified
                {
                    Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
                    blob = m_kr->blob_ID;
                    DensityNAPL = cp_vec[Sp1 - 1]->molar_density;
                    // DensityNAPL = m_kr->Density_NAPL; // CB: this should be
                    // obtained from comp properties
                    m_kb = KinBlob_vector[blob];  // pointer to blob-properties
                                                  // set in the reaction r
                    if (m_Conc[Sp1] > 0.)
                    {
                        // new local data structures
                        oldMass[blob] += DMAX(m_Conc[Sp1], 0.);
                        oldVolume[blob] += DMAX(m_Conc[Sp1], 0.) / DensityNAPL;
                        // Sb todo Achtung - das wird ja gar nicht
                        // zurückgespeichert...
                    }
                }
            }  // end for nreactions
        }

        // 2) Now, calculate current Csat for each species (i.e. reaction)
        // depending on Raoults law for this node.
        //    This requires previously calculated Mass.
        for (r = 0; r < nreactions; r++)
        {
            m_kr = KinReact_vector[r];
            // CB new reaction switch for individual reactions
            if (m_kr->switched_off_node.size() > 0)
                if (m_kr->switched_off_node[node] == true)
                    continue;
            if (m_kr->typeflag_napldissolution)  // dissolution reaction
                                                 // identifie
            {
                if (m_kr->typeflag_gasdissolution == 1)
                {
                    Csat_max = m_kr->GetMaxSolubility(node, 0);
                    m_kr->Current_Csat[node] =
                        Csat_max;  // keine NAPL-Masse vorhanden, NAPL-Bildung
                                   // möglich wenn c(singleSubstance) > Csat
                }
                else
                {
                    Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
                    Csat_max = m_kr->GetMaxSolubility(node, DensityAQ);
                    // Csat_max = cp_vec[Sp1-1]->max_solubility;
                    blob = m_kr->blob_ID;
                    if (oldMass[blob] > 0.)
                        m_kr->Current_Csat[node] =
                            Csat_max * DMAX(m_Conc[Sp1], 0.) / oldMass[blob];
                    else
                        // m_kr->current_Csat = Csat_max;              // keine
                        // NAPL-Masse vorhanden, NAPL-Bildung möglich wenn
                        // c(singleSubstance) > Csat
                        m_kr->Current_Csat[node] =
                            Csat_max;  // keine NAPL-Masse vorhanden,
                                       // NAPL-Bildung möglich wenn
                                       // c(singleSubstance) > Csat
                }
            }
        }  // end for nreactions

        // 3) Next, calculate current Masstransfer-coefficient k for this node
        //    First, get the pore velocity for mass transfer coefficient
        for (r = 0; r < nreactions; r++)
        {
            m_kr = KinReact_vector[r];
            // CB new reaction switch for individual reactions
            if (m_kr->switched_off_node.size() > 0)
                if (m_kr->switched_off_node[node] == true)
                    continue;
            if (m_kr->typeflag_napldissolution)
            {
                PoreVelocity = m_kr->GetNodePoreVelocity(node);
                DarcyVelocity = m_kr->GetNodeDarcyVelocity(
                    node);  // SP: added for Sherwood calculation after Powers
                            // et al. 1992
                break;      // do this only once
            }
        }
        // calculate current Masstransfer-coefficient
        if (Tflag == false)
            DiffusionAQ = mfp_vector[0]->PhaseDiffusion(
                variables);  // SP: Implemented to be temperature dependent
        if (gassdiss == false)
            ViscosityAQ = mfp_vector[0]->Viscosity(
                variables);  // in case of heat transport this T dependent
        // mfp_vector[0]->mode = tmod;
        if (Tflag)
        {
            if (REACTINT_vec.size() > 0)  // Get the Temperature
                // if(REACTINT_vec[0]->constanttemperature){
                variables[1] = REACTINT_vec[0]->GetTemperature(node);
            ViscosityAQ =
                REACTINT_vec[0]->LiquidViscosity_Yaws_1976(variables[1]);
            //}
        }

        // calculate current Masstransfer-coefficient
        m_pcs = PCSGetFlow();
        Poro = KinReact_vector[0]->GetPhaseVolumeAtNode(node, 1, 0);
        idxS = m_pcs->GetNodeValueIndex("SATURATION2") + 1;  // new timelevel
        NAPLcontent = MRange(0, m_pcs->GetNodeValue(node, idxS) * Poro, Poro);
        WATERcontent = Poro - NAPLcontent;

        for (r = 0; r < Number_of_blobs; r++)
        {
            m_kb = KinBlob_vector[r];
            d50 = m_kb->d50;
            shidx = m_kb->shidx;
            if (m_kb->gas_dissolution_flag == false)
            {
                DiffusionAQ = mfp_vector[0]->PhaseDiffusion(
                    variables);  // CB Todo: this should be a component property
                                 // => Sherwood is component dependent
                DensityAQ = mfp_vector[0]->Density(variables);
                // ViscosityAQ  = mfp_vector[0]->Viscosity();
                mod_Reynolds =
                    DensityAQ * PoreVelocity * d50 /
                    ViscosityAQ;  // modified Re' includes transport velo
                Reynolds = DensityAQ * DarcyVelocity * d50 /
                           ViscosityAQ;  // Re includes darcy velo
                Schmidt = ViscosityAQ / DiffusionAQ / DensityAQ;
                Peclet = PoreVelocity * WATERcontent * d50 / DiffusionAQ;
                Delta = d50 / m_kb->dM;
                grain_var = d50 / m_kb->dS;
                // NCont_Sh1    = m_kb->NCont_ini / m_kb->NCont_ini ;
                NCont_Sh2 = NAPLcontent / m_kb->NCont_ini;
                beta_4 = 0.518 + 0.114 * Delta + 0.1 + m_kb->UI;
                one_three = 1.00 / 3.00;
                two_three = 2.00 / 3.00;

                switch (shidx)
                {
                    case 0:  // Standard Sherwood calculation after Bradford and
                             // Abriola (2001)
                        Sherwood = m_kb->Sh_factor *
                                   pow(mod_Reynolds, m_kb->Re_expo) *
                                   pow(Schmidt, m_kb->Sc_expo);
                        break;
                    case 1:  // Sherwood calculation after Powers et al. (1994)
                        // Sherwood = m_kb->Sh_factor * pow(Reynolds,
                        // m_kb->Re_expo) * pow(Delta, m_kb->Delta_expo) *
                        // pow(m_kb->UI, m_kb->UI_expo) * pow(m_kb->NCont_ini,
                        // m_kb->Geometry_expo) ;
                        Sherwood = m_kb->Sh_factor *
                                   pow(mod_Reynolds, m_kb->Re_expo) *
                                   pow(Delta, m_kb->Delta_expo) *
                                   pow(m_kb->UI, m_kb->UI_expo) *
                                   pow(NCont_Sh2, beta_4);
                        break;
                    case 2:  // Sherwood calculation after Powers et al. (1992)
                        Sherwood = m_kb->Sh_factor *
                                   pow(Reynolds, m_kb->Re_expo) *
                                   pow(d50, m_kb->D50_expo) *
                                   pow(m_kb->UI, m_kb->UI_expo);
                        break;
                    case 3:  // Sherwood calculation after Wilson and Geankoplis
                             // (1966)
                        Sherwood = m_kb->Sh_factor * pow(Poro, -1) *
                                   pow(Peclet, one_three);
                        break;
                    case 4:  // Sherwood calculation after Pfannkuch (1984)
                        Sherwood = m_kb->Sh_factor + m_kb->Pfannkuch_constant *
                                                         pow(Peclet, two_three);
                        break;
                    case 5:  // Sherwood calculation after Miller et al. (1990)
                        Sherwood = m_kb->Sh_factor *
                                   pow(mod_Reynolds, m_kb->Re_expo) *
                                   pow(WATERcontent, m_kb->WContent_expo) *
                                   pow(Schmidt, m_kb->Sc_expo);
                        break;
                    case 6:  // Sherwood calculation after Saba&Illangasekare
                             // (2000)
                        Sherwood =
                            m_kb->Sh_factor * pow(mod_Reynolds, m_kb->Re_expo) *
                            pow(Schmidt, m_kb->Sc_expo) *
                            ((d50 * NAPLcontent) / (m_kb->Tort * m_kb->Length));
                        break;
                    case 7:  // Sherwood calculation after Geller&Hunt (1993)
                        Sherwood =
                            m_kb->Sh_factor * pow(mod_Reynolds, m_kb->Re_expo) *
                            pow(NAPLcontent, m_kb->NContent_expo) *
                            pow(m_kb->NCont_ini, m_kb->NContent_ini_expo) *
                            pow(Poro, m_kb->Poro_expo) *
                            pow(grain_var, m_kb->grain_expo);
                        break;
                    default:
                        break;
                }
            }
            else if (m_kb->gas_dissolution_flag == true)
                Sherwood =
                    d50 *
                    (4 / d50 +
                     pow((2 * PoreVelocity / PI / d50 / DiffusionAQ), 0.5));

            if (m_kb->Sherwood_model == false)
                m_kb->Masstransfer_K[node] =
                    Sherwood * DiffusionAQ / d50;  // k in m/s
            else if (m_kb->Sherwood_model ==
                     true)  // modified Sh'= k * (d50)^2 / Daq --> k =
                            // Sh*Daq/(d50)^2
                m_kb->Masstransfer_K[node] =
                    Sherwood * DiffusionAQ / (d50 * d50);  // k in m/s
        }

        // This step is necessary only for gas dissolution, as saturation
        // changes over time 4) Finally, save current Interfacial areas for this
        // node.
        for (r = 0; r < Number_of_blobs; r++)
        {
            m_kb = KinBlob_vector[r];
            if (m_kb->gas_dissolution_flag)
            {
                Poro = KinReact_vector[0]->GetPhaseVolumeAtNode(node, 1, 0);
                Sat = KinReact_vector[0]->GetSaturation(node, 1, 1);
                m_kb->Interfacial_area[node] =
                    12 * Sat * Poro / d50 * m_kb->Area_Value[0];
            }
        }
    }  //  if No of blobs > 0

    /**********************************************/
    /* ODE Integration */
    /**********************************************/
    tstart = DMAX(aktuelle_zeit - dt, 0.);
    tend = aktuelle_zeit;
    // calculate substep length
    tstart += dt * double(tt) / double(tsteps);
    tend -= dt * double(tsteps - (tt + 1)) / double(tsteps);
    /* Aufruf Gleichungsl�ser */
    /* eigentliche Rechenroutinen sind "derivs" und "jacobn" (namen fest
       vorgegeben), die vom Gleichungsl�ser aufgerufen werden */
    if (SolverType == 0)
    {
        Calc_linearized_rates(m_Conc, Number_of_Components, tend - tstart,
                              node);
        success = true;  // works always
        *nok = 1;
        *nbad = 0;
    }
    else
    {
        success =
            odeint(m_Conc, Number_of_Components, tstart, tend, eps, usedt, hmin,
                   &nexth, nok, nbad, derivs, stifbs, rkqs, node, SolverType);

        if (success == 0)
        {  // if odeint failed for some reason
            cout << "Warning in CKinReactData::Biodegradation: Odeint did not "
                    "converge at node "
                 << node << "."
                 << "\n";
            cout << "   -Trying increased relative error " << eps * 10 << "\n";
            success = odeint(m_Conc, Number_of_Components, tstart, tend,
                             eps * 10, usedt, hmin, &nexth, nok, nbad, derivs,
                             stifbs, rkqs, node, SolverType);
            if (success == 0)
            {
                cout << "  -Trying perturbation of initial time step usedt"
                     << "\n";
                success = odeint(m_Conc, Number_of_Components, tstart, tend,
                                 eps, usedt * 0.99, hmin, &nexth, nok, nbad,
                                 derivs, stifbs, rkqs, node, SolverType);
                if (success == 0)
                {
                    cout << "   -Trying conversion of C vector from double to "
                            "float."
                         << "\n";
                    for (sp = 0; sp < Number_of_Components; sp++)
                    {
                        tempstore.push_back(m_Conc[sp + 1]);
                        garage = float(m_Conc[sp + 1]);
                        m_Conc[sp + 1] = double(garage);
                    }
                    success =
                        odeint(m_Conc /*_save*/, Number_of_Components, tstart,
                               tend, eps, usedt, hmin, &nexth, nok, nbad,
                               derivs, stifbs, rkqs, node, SolverType);
                    // if(success) for(sp=0;sp<Number_of_Components;sp++)
                    // m_Conc[sp+1] = m_Conc_save[sp+1] ;
                    /*else*/ if (success == 0)
                    {
                        cout << "   -Trying perturbation of concentration "
                                "vector by 1 per mil"
                             << "\n";
                        for (sp = 0; sp < Number_of_Components; sp++)
                            m_Conc[sp + 1] =
                                tempstore[sp];  // restore original double value
                                                // precision
                        for (sp = 0; sp < Number_of_Components; sp++)
                            m_Conc[sp + 1] *= 1.01;  // pertube
                        success =
                            odeint(m_Conc, Number_of_Components, tstart, tend,
                                   eps, usedt, hmin, &nexth, nok, nbad, derivs,
                                   stifbs, rkqs, node, SolverType);
                        for (sp = 0; sp < Number_of_Components; sp++)
                            m_Conc[sp + 1] /= 1.01;  // restore Masses
                        if (success == 0)
                        {
                            cout << "   -Trying increased relative error "
                                 << eps * 100 << "\n";
                            success = odeint(m_Conc, Number_of_Components,
                                             tstart, tend, eps * 100, usedt,
                                             hmin, &nexth, nok, nbad, derivs,
                                             stifbs, rkqs, node, SolverType);
                            if (success == 0)
                                cout << "  ODE integration finally failed. "
                                        "Skipping this node."
                                     << "\n";
                        }
                    }
                }
            }
        }
    }
    if (success)
        baditerations = double(*nbad) / double(*nok + *nbad);
    else
        baditerations = 1;
    if (baditerations < 0.001)
    {
        // fehlerfreie Integration, zeitschritt kann vergr��ert werden
        if (nexth > usedt)
            *usedtneu = DMAX(nexth, usedt * 2.);
        else
            *usedtneu = usedt * 1.5;
    }
    else
    {
        // Integrationsfehler, zeitschritt beibehalten oder verkleinern
        if (*nbad == 1)
            *usedtneu = DMAX(nexth, usedt * 1.10);
        else if (*nok > *nbad * 2)
            *usedtneu = DMAX(nexth, usedt * 1.01);
        else
            *usedtneu = DMAX(nexth, usedt / 5.);
    }

    // update results, now outside in ExecuteKinReact

    if (debugoutflag)
    {
        // debugoutstr << " Concentrations after odeint: " << "\n" << " " <<
        // flush;
        debugoutstr << node << " " << aktueller_zeitschritt << " " << flush;
        for (sp = 0; sp < Number_of_Components; sp++)
            debugoutstr << m_Conc[sp + 1] << " " << flush;
        debugoutstr << KinReact_vector[0]->GetSaturation(node, 1, 1) << " "
                    << flush;
        debugoutstr << nexth << " " << int(*nok) << " " << int(*nbad) << "\n"
                    << flush;
    }

    /**********************************************/
    /* Postprocessing NAPL Dissolution */
    /* calculate Interfacial areas for this node after dissolution for next time
     * step */
    /**********************************************/
    if (Number_of_blobs > 0)
    {
        if (KinBlob_vector[0]->gas_dissolution_flag == false)
        {
            for (r = 0; r < nreactions; r++)
            {
                m_kr = KinReact_vector[r];
                // CB new reaction switch for individual reactions
                if (m_kr->switched_off_node.size() > 0)
                    if (m_kr->switched_off_node[node] == true)
                        continue;
                if (m_kr->typeflag_napldissolution)  // dissolution reaction
                                                     // identified
                {
                    Sp1 = m_kr->ex_species[0] + 1;  // Sp1 = NAPL-species
                    blob = m_kr->blob_ID;
                    DensityNAPL = cp_vec[Sp1 - 1]->molar_density;
                    newVolume[blob] += DMAX(m_Conc[Sp1], 0.) / DensityNAPL;
                }
            }  // end for nreactions
            for (r = 0; r < Number_of_blobs; r++)
            {
                m_kb = KinBlob_vector[r];
                if ((newVolume[r] > 0.) && (oldVolume[r] > 0.))
                    m_kb->Interfacial_area[node] *=
                        pow((newVolume[r] / oldVolume[r]), m_kb->Geometry_expo);
                else
                    m_kb->Interfacial_area[node] =
                        1.E-20;  // residual interfacial area to allow
                                 // re-building of phase
            }
        }
        // Free storage
        free_dvector(newVolume, 0, Number_of_blobs);
        free_dvector(oldVolume, 0, Number_of_blobs);
        free_dvector(oldMass, 0, Number_of_blobs);
    }  //  if No of blobs > 0
    // free_dvector(Concentration,1,Number_of_Components);
    // free(m_Conc_save);
    tempstore.clear();
}

void CKinReactData::Calc_linearized_rates(double* m_Conc,
                                          long Number_of_Components,
                                          double deltaT, long node)
{
    double* dydx;

    dydx = dvector(1, Number_of_Components);
    derivs(deltaT, m_Conc, dydx, Number_of_Components, node, deltaT);
    for (long i = 0; i < Number_of_Components; i++)
        m_Conc[i + 1] += dydx[i + 1] * deltaT;
}

/*************************************************************************************/
/* Routine f�r Bulirsch-Stoer Gleichungsl�ser (steife ODEs) */
/* DS-TBC */
/*                                                                                   */
/* Berechnung der ersten Ableitungen �ber die Zeit dC/dt=... */
/*                                                                                   */
/* Input: */
/* t = aktuelle Zeit */
/* c[1..Number_of_Components] = aktuelle Konzentration der Substanzen in Zelle
 */
/* n =! Number_of_Components */
/*                                                                                   */
/* Output: */
/* dcdt[1..Number_of_Components] = aktuelle Ableitung der Konzentrationen nach
 * der   */
/*            Zeit, = Wachstum der Bakterien und Verbrauch der Substanzen */
/*                                                                                   */
/* Programmaenderungen: */
/* 05/2003     DS         Erste Version */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/* 05/2007     DS         NAPL-dissolution added */
/*                                                                                   */
/*************************************************************************************/

void derivs(double t, double c[], double dcdt[], int n, long node,
            double /*steplength*/)
{
    t = t;  // OK411
    int i, j, r, nreactions, BacteriaNumber;
    int Sp1, Sp2, surfaceID = -1, blob;  // phase,
    double BacteriaMass, BacGrowth, Yield, sumX = 0., maxkap;
    double porosity1, porosity2, exchange = 0.0, exch, kd, density1,
                                 saturation2, kadsorb, kdesorb, totalSurface,
                                 exponent, parameter,
                                 chochexp;  // OK411
    double ratefact = 1;
    double dt;
    double foc;
    //#ds
    //	int blob;
    double Csat;
    //	double occupiedSurface[m_krd->maxSurfaces+1];
    vector<double> occupiedSurface;
    // MinKin
    double MinGrowth, omega;
    int MineralNumber;
    double scale = 1;
    double scalemin = 1;
    // scale rates correctly
    std::vector<double> ratevector;
    std::vector<std::vector<double> > ratematrix;

    // phase = 0;
    CKinReact* m_kr = NULL;
    // OK411 CKinReact *m_kr1 = NULL;
    CKinBlob* m_kb = NULL;
    CKinReactData* m_krd = NULL;
    m_krd = KinReactData_vector[0];
    /* Anzahl der mikrobiellen Reaktionen aus Datenstruktur auslesen */
    nreactions = m_krd->NumberReactions;  // BioDegradation.NumberReactions;
    // if(m_krd->debugoutflag)
    //  m_krd->debugoutstr << " derivs" << "\n" << flush;

    CTimeDiscretization* m_tim = NULL;
    m_tim = time_vector[0];
    dt = m_tim->CalcTimeStep();

    if (m_krd->scale_dcdt)
    {
        for (i = 0; i < nreactions; i++)
            ratevector.push_back(0);
        for (i = 0; i < n; i++)
            ratematrix.push_back(ratevector);
        ratevector.clear();
    }
    /* reset array with derivatives */
    /* ACHTUNG, unterschiedliche Indizierung der Arrays, c[1..n]
     * BioDegradation[0..n-1] */
    for (i = 0; i < n; i++)
        // SBtodo		dcdt[i+1]=ExchangeTerm[i];
        dcdt[i + 1] = 0.0;

    /* calculate present bacteria capacity */
    sumX = 0.0;  // SB added
    maxkap = m_krd->maxBacteriaCapacity;
    if (maxkap > 1.E-30)
        for (i = 0; i < n; i++)
            if (m_krd->is_a_bacterium[i])
                sumX = sumX + c[i + 1];

    /**********************************************************************************************/
    /* Anzahl der mikrobiellen Reaktionen aus Datenstruktur auslesen */
    nreactions = m_krd->NumberReactions;  // BioDegradation.NumberReactions;
    /* loop over reactions dX/dt= nymax * X * monodterms * inhibitionterms */
    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;
        if (m_kr->typeflag_monod)
        {
            BacteriaNumber = m_kr->bacteria_number + 1;
            BacteriaMass = c[BacteriaNumber];
            porosity1 = m_kr->GetReferenceVolume(BacteriaNumber - 1, node);
            if (BacteriaMass > 1.E-40)
            {
                m_kr->currentnode = node;  // CB 19/10/09 This is eclusively for
                                           // Brand model to allow porosity in
                // Inhibition constant calculation

                // This is where growth rate is computed
                BacGrowth = m_kr->BacteriaGrowth(r, c, sumX, -1, node);
                if (m_kr->grow)
                    dcdt[BacteriaNumber] +=
                        BacGrowth * m_kr->EffectiveYield(
                                        node);  // CB 08/12 added a microbial
                                                // yield coefficient
                /* microbial consumption of substances */
                for (i = 0; i < n; i++)
                {
                    Yield =
                        m_kr->ProductionStoch[i];  // CB 08/12 this is the
                                                   // stoichiometry coefficient
                                                   // of the reaction
                    if (fabs(Yield) > 1.E-30)
                    {
                        porosity2 = m_kr->GetReferenceVolume(i, node);
                        dcdt[i + 1] +=
                            BacGrowth * Yield * porosity1 / porosity2;
                    }
                }
            }
        }  // type == monod
    }      // nreactions

    /**********************************************************************************************/
    /* Berechnung der Mineralkinetiken */
    if (m_krd->NumberMineralkinetics > 0)
    {
        for (r = 0; r < nreactions; r++)
        {
            m_kr = KinReact_vector[r];
            // CB new reaction switch for individual reactions
            if (m_kr->switched_off_node.size() > 0)
                if (m_kr->switched_off_node[node] == true)
                    continue;
            if (m_kr->typeflag_mineralkinetics)
            {
                // Equilibrium: when rate is zero, division by 0 (1-omega^theta)
                // in fct jacbn is a problem, therefore, do nothing here when in
                // in equilibrium
                omega = m_kr->Omega(c, node);
                if (fabs(1 - pow(omega, m_kr->Theta)) > 0)
                {
                    MineralNumber = m_kr->mineral_number + 1;
                    // This is where growth rate is computed
                    MinGrowth = m_kr->MinRate(c, node, dt, omega);
                    // cout << MinGrowth << " ";
                    // unit: --> mol/s/m³_solid
                    if (MinGrowth != 0)
                    {
                        dcdt[MineralNumber] += MinGrowth;
                        if (m_krd->scale_dcdt)
                            ratematrix[MineralNumber - 1][r] = MinGrowth;
                        // mineral species precipitation or dissolution (not for
                        // the actual mineral,
                        // ProductionStoch[MineralNumber-1]=0)
                        porosity1 =
                            m_kr->GetReferenceVolume(MineralNumber - 1, node);
                        for (i = 0; i < n; i++)
                        {
                            Yield = m_kr->ProductionStoch[i];
                            porosity2 = m_kr->GetReferenceVolume(i, node);
                            if (fabs(Yield) > 1.E-30)
                            {
                                // unit conversion: --> mol/s/m³_solid *
                                // m³sol/m³Aq * (m³w/m³Aq)^-1 = mol/s/m³w
                                dcdt[i + 1] +=
                                    MinGrowth * Yield * porosity1 / porosity2;
                                if (m_krd->scale_dcdt)
                                    ratematrix[i][r] += MinGrowth * Yield *
                                                        porosity1 / porosity2;
                            }
                        }
                    }
                }
            }  // type == mineralkinetics
        }      // nreactions
    }

    if (m_krd->scale_dcdt)
    {
        // cout << "derivs: dt=" << dt << "s, node=" << node << ": " << "\n";
        // for (i=0; i<n; i++)
        //    cout << c[i+1] << " " << dcdt[i+1]*dt << " " << dcdt[i+1] << " "
        //    << c[i+1]+dcdt[i+1]*dt << "\n";
        // cout << "\n";

        // limit total dc/dt to the present amount of species to avoid negative
        // concentrations
        scale = scalemin = 1;
        for (i = 0; i < n; i++)
        {
            // find the species, which go negative from sum of all rates
            if ((c[i + 1] + dcdt[i + 1] * dt < 0) && (dcdt[i + 1] * dt != 0))
            {
                scale = c[i + 1] / fabs(dcdt[i + 1] * dt);
                // then rescale rates of all reactions involving this species
                for (r = 0; r < nreactions; r++)  // go through all reactions
                    if (ratematrix[i][r] !=
                        0)  // is species i is involved, rescale the reaction?
                        for (j = 0; j < n; j++)
                        {
                            if (m_krd->SolverType ==
                                0)  // for linearized rate calculation
                                ratematrix[j][r] *= scale;
                            else
                                ratematrix[j][r] *= scale * 0.9;
                        }
                // then recompute dcdt
                for (j = 0; j < n; j++)
                {  // species by species
                    dcdt[j + 1] = 0;
                    for (r = 0; r < nreactions; r++)
                        dcdt[j + 1] +=
                            ratematrix[j][r];  // sum up contributions of all
                                               // reactions
                }
            }
            // scalemin= DMIN(scalemin, scale); // scalemin is obsolete with
            // introduction of ratematrix
        }
        // if(scalemin<1)
        //  for (i=0; i<n; i++) dcdt[i+1]*= (scalemin*0.9);

        // cout << "derivs: dt=" << dt << "s, node=" << node << ": " << "\n";
        // for (i=0; i<n; i++)
        //    cout << c[i+1] << " " << dcdt[i+1]*dt << " " << dcdt[i+1] << " "
        //    << c[i+1]+dcdt[i+1]*dt << "\n";
        // cout << "\n";
    }
    /**********************************************************************************************/
    /* Berechnung der Austauschprozesse */

    // calculate already occupied surface
    if (m_krd->NumberLangmuir > 0)
    {
        // Initialise Surfaces for langmuir isotherms
        for (i = 0; i < m_krd->maxSurfaces; i++)
            occupiedSurface.push_back(0.0);

        for (r = 0; r < nreactions; r++)
        {
            m_kr = KinReact_vector[r];
            // CB new reaction switch for individual reactions
            if (m_kr->switched_off_node.size() > 0)
                if (m_kr->switched_off_node[node] == true)
                    continue;
            if ((m_kr->typeflag_exchange) && (m_kr->typeflag_exchange_langmuir))
            {
                Sp1 = m_kr->ex_species[0] + 1;
                surfaceID = m_kr->exSurfaceID;
                occupiedSurface[surfaceID] += c[Sp1];
            }
        }
    }  // if NumberLangmuir > 0

    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;
        if (m_kr->typeflag_exchange)
        {
            /* linearer Austausch ggf. mit kd */
            if (m_kr->typeflag_exchange_linear)
            {
                // Matrix
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                // geloest
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
                ratefact = m_kr->ex_param[2];

                exch = m_kr->ex_param[0];
                kd = m_kr->ex_param[1];

                if (fabs(kd) < MKleinsteZahl)
                {
                    // no kd, exchange between two species in solution

                    // slow down desorption rate by constant factor
                    if ((c[Sp2] - c[Sp1]) < 0)
                        exch *= ratefact;

                    exchange = exch * (c[Sp2] - c[Sp1]);
                    dcdt[Sp1] += exchange / porosity1;
                    dcdt[Sp2] += -exchange / porosity2;
                }
                else
                {
                    // with kd, exchange between matrix (mol/kg) and solution
                    // (mol/l)
                    foc = m_krd->node_foc[node];

                    if (foc > MKleinsteZahl)
                        kd = kd * foc;
                    // else
                    //  kd = 0;

                    // slow down desorption rate by constant factor
                    if ((kd * c[Sp2] - c[Sp1]) < 0)
                        exch *= ratefact;

                    exchange = exch * (c[Sp2] * kd - c[Sp1]);
                    /* Die Abfrage verringert die Desorptionsgeschwindigkeit,
                     * wenn absehbar ist, dass Csorbiert im Negativen landet */
                    if (-exchange * dt > c[Sp1])
                        exchange = -c[Sp1] / dt;

                    dcdt[Sp1] += exchange;
                    dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
                }
            }  // ende if exType == linear

            /* Freundlich Kinetik */
            if (m_kr->typeflag_exchange_freundlich)
            {
                // Matrix
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                // geloest
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

                // if both species are in the same phase, density multiplication
                // below is not required
                if (cp_vec[Sp1 - 1]->transport_phase ==
                    cp_vec[Sp2 - 1]->transport_phase)
                    density1 = 1;

                exponent = m_kr->ex_param[2];
                parameter = m_kr->ex_param[1];
                exch = m_kr->ex_param[0];
                ratefact = m_kr->ex_param[3];

                if (c[Sp2] > residual)
                    // no linearisation required
                    chochexp = pow(c[Sp2], exponent);
                else
                    // linearisation required due to instability of c^x if
                    // c<residual
                    chochexp =
                        (1. - exponent) * pow(residual, exponent) +
                        exponent * pow(residual, (exponent - 1)) * c[Sp2];

                exchange = exch * (parameter * chochexp - c[Sp1]);
                // slow down desorption rate by constant factor
                if (exchange < 0)
                    exchange *= ratefact;
                /* Die Abfrage verringert die Desorptionsgeschwindigkeit, wenn
                 * absehbar ist, dass Csorbiert im Negativen landet */
                if (-exchange * dt > c[Sp1])
                    exchange = -c[Sp1] / dt;

                dcdt[Sp1] += exchange;
                dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
            }  // if freundlich

            /* Langmuir Kinetik */
            if (m_kr->typeflag_exchange_langmuir)
            {
                // Surfaces were initialized above
                //	for (i=0; i<nexchange; i++)	{
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);

                kadsorb = m_kr->ex_param[0];
                kdesorb = m_kr->ex_param[1];
                // SB_langmuir		surfaceID	= m_kr->exSurfaceID;
                // //Exchange.Langmuir[i].SurfaceID;
                totalSurface = m_krd->exSurface[m_kr->exSurfaceID];

                //#ds ACHTUNG hier muss sicher gestellt sein, dass Sp1 die
                // adsorbierte und Sp2 die gel�ste Species ist !
                exchange = kadsorb *
                               (totalSurface - occupiedSurface[surfaceID]) *
                               c[Sp2] -
                           kdesorb * c[Sp1];

                /* Die Abfrage verringert die Desorptionsgeschwindigkeit, wenn
                 * absehbar ist, dass Csorbiert im Negativen landet */
                if (-exchange * dt > c[Sp1])
                    exchange = -c[Sp1] / dt;

                dcdt[Sp1] += exchange;
                dcdt[Sp2] += -exchange * porosity1 / porosity2 * density1;
                //	}
            }  // ende if exType == langmuir
        }      // if type == exchange
    }          // for r

    //#ds
    /**********************************************************************************************/
    /* Berechnung NAPL-L�sung */
    /**********************************************************************************************/

    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;

        if (m_kr->typeflag_napldissolution)
        {
            Sp1 = m_kr->ex_species[0] +
                  1;  // Exchange.Linear[i].Species1;    Sp1 muss NAPL sein
            //		porosity1	= m_kr->GetPorosity(Sp1-1,node);
            //		density1	= m_kr->GetDensity(Sp1-1,node);
            ////GetDensity(phase);
            blob = m_kr->blob_ID;
            m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
                                          // the reaction r

            Sp2 = m_kr->ex_species[1] +
                  1;  // Exchange.Linear[i].Species2;    Sp2 = mobile Phase
            // CB this includes the saturation
            porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
            //#ds TODO	    saturation2 = ??
            saturation2 = 1.;

            // CB now access data for each node directly from node vectors
            Csat = m_kr->Current_Csat[node];  // Csat externally calculated in
                                              // Function Biodegradation
            exch = m_kb->Masstransfer_K[node] *
                   m_kb->Interfacial_area[node];  // k * A externally calculated
                                                  // in Function Biodegradation
            exponent = m_kr->rateorder;

            exchange = 0;
            if (exch > 0)
                exchange =
                    exch * pow((Csat - c[Sp2]),
                               exponent);  // CB introduced rateorder 06/2012

            /* Die Abfrage verringert die L�sungsgeschwindigkeit, wenn absehbar
               ist, dass CNAPL im Negativen landet Verhindert negative
               CNAPL-Konzentrationen */
            if (exchange * dt > c[Sp1])
            {
                // k1 =  c[Sp1] / dt;
                // k2 =  (Csat - c[Sp2]) * porosity2 * saturation2 / dt;
                exchange = c[Sp1] / dt;
                // exchange = DMIN(k1, k2);
                // cout << "Warning in derivs: exchange * dt > c[Sp1]" << "\n";
            }
            //         if ((exchange > 0) && (exchange * dt + c[Sp2] > Csat))
            //            exchange = (Csat - c[Sp2]) / dt /5;
            if (c[Sp1] > MKleinsteZahl)  // no dissolution or NAPL mass present
                if ((exchange < 0.) || (c[Sp1] > MKleinsteZahl))
                {
                    dcdt[Sp1] += -exchange;  // NAPL
                    // concentration in solution, refers to mass / volume of
                    // water
                    dcdt[Sp2] += exchange / porosity2 / saturation2;
                }
        }  // ifNAPL dissolution
    }      // loop ofer reactions r
    //#ds

    occupiedSurface.clear();
    ratematrix.clear();
}

/**************************************************************************/
/* Phase bestimmen, in der sich Substanz i befindet                       */
/* DS-TBC                                                                 */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/**************************************************************************/
int CKinReact::GetPhase(int species)
{
    int phase = -1;
    CompProperties* cp_m = NULL;

    cp_m = cp_vec[species];
    if (cp_m != NULL)
        phase = cp_m->transport_phase;
    else
        cout << " Error: component does not exist !"
             << "\n";
    return phase;
}

/*****************************************************************************************/
/* Calculate the reference volume of a phase at a node */
/* DS-TBC */
/* 02/2006     SB         Introduced new C++ concept, Data structures */
/* 08/2008     DS         Consider saturation of water phase in case of
 * multiphase flow  */
/* 09/2009     CB         Heterogeneous Porosities update */
/*****************************************************************************************/
double CKinReact::GetReferenceVolume(int comp, long index)
{
    // Get process
    CRFProcess const* const pcs(cp_vec[comp]->getProcess());
    double theta(pcs->m_num->ls_theta), refvol = 0;
    long phase = cp_vec[comp]->transport_phase;

    if (phase == 0)
    {
        int idx;
        refvol = GetPhaseVolumeAtNode(index, theta, phase);
        // water phase, reference volume might be less than pore space in case
        // of multiphase or richards flow
        // --> Get node saturation of mobile (water) phase, required for all
        // exchange processes
        double saturation = 1.0;  // default
        CRFProcess* pcs_flow = PCSGetFlow();
        if (pcs_flow->getProcessType() == FiniteElement::PS_GLOBAL)
        {
            idx = pcs_flow->GetNodeValueIndex(
                "SATURATION1");  // Sat of water phase
            saturation = pcs_flow->GetNodeValue(index, idx);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
        {
            if (pcs_flow->pcs_type_number == 0)
                // this is the saturation equation
                pcs_flow = pcs_vector[pcs_flow->pcs_number + 1];
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = pcs_flow->GetNodeValue(index, idx);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::RICHARDS_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = pcs_flow->GetNodeValue(index, idx);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = pcs_flow->GetNodeValue(index, idx);
        }
        refvol *= saturation;
    }
    else if (phase == 1)  // solid
    {
        refvol = GetPhaseVolumeAtNode(index, theta, phase);
    }
    else if (phase == 2)  // bio phase
        refvol = GetPhaseVolumeAtNode(index, theta, phase);
    else if (phase == 3)  // NAPL phase (refers to REV)
        refvol = 1.0;
    else
        cout << " Warning: Unknown phase index " << phase
             << " in CKinReact::GetReferenceVolume()."
             << "\n";

    return refvol;
}

/*****************************************************************************************/
/* Calculate the reference volume of a phase at a node */
/* DS-TBC */
/* 02/2006     SB         Introduced new C++ concept, Data structures */
/* 08/2008     DS         Consider saturation of water phase in case of
 * multiphase flow  */
/* 09/2009     CB         Heterogeneous Porosities update */
/*****************************************************************************************/
double CKinReact::GetSaturation(long node, int theta, int phase)
{
    double saturation = 1;  // refvol = 0.0,
    CRFProcess* pcs_flow = PCSGetFlow();
    int idx;

    //
    if (phase == 0)
    {  // water
        if (pcs_flow->getProcessType() == FiniteElement::PS_GLOBAL)
        {
            idx = pcs_flow->GetNodeValueIndex(
                "SATURATION2");  // Sat of NAPL phase
            saturation = 1.0 - pcs_flow->GetNodeValue(node, idx + theta);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::RICHARDS_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = pcs_flow->GetNodeValue(node, idx + theta);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = pcs_flow->GetNodeValue(node, idx + theta);
        }
    }
    //
    else if (phase == 1)  // NAPL, gas
    {
        if (pcs_flow->getProcessType() == FiniteElement::PS_GLOBAL)
        {
            idx = pcs_flow->GetNodeValueIndex(
                "SATURATION2");  // Sat of NAPL phase
            saturation = pcs_flow->GetNodeValue(node, idx + theta);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::RICHARDS_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = 1.0 - pcs_flow->GetNodeValue(node, idx + theta);
        }
        else if (pcs_flow->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        {
            // Sat of water phase
            idx = pcs_flow->GetNodeValueIndex("SATURATION1");
            saturation = 1.0 - pcs_flow->GetNodeValue(node, idx + theta);
        }
    }

    return saturation;
}

/**************************************************************************/
/* Return the volume fraction of a particular phase at a node             */
/* 0 pore space, 1 solid phase, 2 bio phase                               */
/* DS-TBC                                                                 */
/* 09/2009     CB         Introduced new C++ concept, Data structures     */
/* 09/2011	TF hanged access to coordinates of mesh node,
    - substituted access to mesh_element from pointer to direct access into the
   vector
    - made the mesh node a const pointer
    - made the pointer to the mesh const, made the mesh itself const
    - substituted pow(x,2) by x*x
    - reduced scope of loop variable i
    - moved declaration of gravity center to appropriate place			  */
/**************************************************************************/
double CKinReact::GetPhaseVolumeAtNode(long node_number, double theta,
                                       int phase)
{
    CFEMesh const* const mesh(fem_msh_vector[0]);  // SB: ToDo hart gesetzt

    size_t idx = 0;
    long elem;  // OK411
    double distance, weight, sum_w = 0;
    double vol = 0, poro = 0;

    // if node poro vector is used
    if (phase == 0 || phase == 1)
    {
        REACTINT* m_rei = NULL;
        if (REACTINT_vec.size() > 0)
        {
            m_rei = REACTINT_vec[0];
            vol = m_rei->node_porosity[node_number];
            if (phase == 1)  // solid phase
                vol = 1.0 - vol;
            return vol;
        }
    }
    // get Indices for phase 1 or 2, only if heterogeneous porosity model = 11,
    // i.e. vol_mat_model = vol_bio_model = 2
    long group = 0;  // SB todo group = m_ele->GetPatchIndex(); Todo CB
    CMediumProperties* m_mat_mp(mmp_vector[group]);
    if (m_mat_mp->vol_bio_model == 2 && m_mat_mp->vol_mat_model == 2)
    {
        switch (phase)
        {
            case 1:  // solid phase
                // Get VOL_MAT index
                for (idx = 0; idx < mesh->mat_names_vector.size(); idx++)
                    if (mesh->mat_names_vector[idx].compare("VOL_MAT") == 0)
                        break;
                break;
            case 2:  // bio phase
                // Get VOL_BIO index
                for (idx = 0; idx < mesh->mat_names_vector.size(); idx++)
                    if (mesh->mat_names_vector[idx].compare("VOL_BIO") == 0)
                        break;
                break;
            default:
                break;
        }
    }

    // Get node coordinates
    MeshLib::CNode const* const node(mesh->nod_vector[node_number]);
    double const* const coord = node->getData();  // Coordinates(coord);

    for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
    {
        // initialize for each connected element
        distance = weight = poro = 0;
        // Get the connected element
        elem = node->getConnectedElementIDs()[el];  // element index
        MeshLib::CElem const* const m_ele(mesh->ele_vector[elem]);
        // get the phase volume of current element elem
        group = 0;  // group = m_ele->GetPatchIndex(); Todo CB
        m_mat_mp = mmp_vector[group];
        switch (phase)
        {
            case 0:  // pore space
                // CB Now provides also heterogeneous porosity, model 11
                poro = m_mat_mp->Porosity(elem, theta);
                break;
            case 1:                                // solid phase
                if (m_mat_mp->vol_mat_model == 1)  // homogeneous
                    poro = m_mat_mp->vol_mat;
                else if (m_mat_mp->vol_mat_model == 2)  // CB heterogeneous
                    poro = m_ele->mat_vector(idx);
                else
                    cout << "Warning! No valid VOL_MAT model in "
                            "CKinReact::GetPhaseVolumeAtNode, vol_mat_model ="
                         << m_mat_mp->vol_mat_model << "\n";
                break;
            case 2:                                // bio phase
                if (m_mat_mp->vol_bio_model == 1)  // homogeneous
                    poro = m_mat_mp->vol_bio;
                else if (m_mat_mp->vol_bio_model == 2)  // CB heterogeneous
                    poro = m_ele->mat_vector(idx);
                else
                    cout << "Warning! No valid VOL_BIO model in "
                            "CKinReact::GetPhaseVolumeAtNode, vol_bio_model ="
                         << m_mat_mp->vol_bio_model << "\n";
                break;
            case 3:  // NAPL phase (refers to REV)
                poro = 1.0;
                break;
            default:
                cout << "Error in CKinReact::GetPhaseVolumeAtNode: no valid "
                        "phase"
                     << "\n";
                break;
        }
        // calculate distance node <-> element center of gravity
        double const* grav_c(m_ele->GetGravityCenter());
        for (size_t i = 0; i < 3; i++)
            distance +=
                (coord[i] - grav_c[i]) *
                (coord[i] - grav_c[i]);  // pow((coord[i] - grav_c[i]), 2);
        // linear inverse distance weight = 1/(distance)
        distance =
            sqrt(distance);  // for quadratic interpolation uncomment this line
        weight = (1 / distance);
        sum_w += weight;
        // add the weighted phase volume
        vol += poro * weight;
    }  // loop over connected elements

    // normalize weighted sum by sum_of_weights sum_w
    vol *= 1 / sum_w;

    return vol;
}

/**************************************************************************/
/* Max Solubility of a NAPL component                                     */
/* 05/2012     CB                                                         */
/**************************************************************************/
double CKinReact::GetMaxSolubility(long node, double density)
{
    double Csatmax = 0;
    int Sp1 = ex_species[0] + 1;  // Sp1 = NAPL-species
    double Temperature = 0;
    double logK = 0;
    double R = 8.314472;  // gas constant = R 8.314472(15) J/K/mol

    if (typeflag_gasdissolution == 1)
    {
        if (REACTINT_vec.size() > 0)
            Csatmax = REACTINT_vec[0]->GetCO2SolubilityDuan(node);
        else
        {
            cout << " Warning in CKinReact::GetMaxSolubility(): "
                 << "\n";
            cout << " No Reaction Interface REI defined, but required for CO2 "
                    "solubility calculation."
                 << "\n";
            cout << " Using constant value from cp_data."
                 << "\n";
            Csatmax = cp_vec[Sp1 - 1]->max_solubility;  // temporary substitute
        }
        return Csatmax;
    }

    if (T_dependence == false)
        Csatmax = cp_vec[Sp1 - 1]->max_solubility;
    else if (T_dependence == true)
    {
        if (REACTINT_vec.size() > 0)  // Get the Temperature
            Temperature = REACTINT_vec[0]->GetTemperature(node);
        if (T_model == 1)
        {  // Chen Model : R*ln(K) = A + B/T + C*ln(T); based on molar fraction
           // mol/molw
            logK = (T_params[0] + T_params[1] / Temperature +
                    T_params[2] * log(Temperature)) /
                   R;             // nat log
            Csatmax = exp(logK);  // mol/mol wasser
            Csatmax *=
                (1 / 0.01801528) * density;  // mol/molw * molw/kgw * kgw/m³
        }
        else if (T_model == 2)
        {  // Knauss Model : R*ln(K) = A + B/T + C*ln(T); based on molality
           // mol/kgw
            logK = (T_params[0] + T_params[1] / Temperature +
                    T_params[2] * log(Temperature)) /
                   R;             // nat log
            Csatmax = exp(logK);  // mol/kgw  wasser
            Csatmax *= density;   // mol/kgw * kgw/m³
        }
        else
        {
            cout << " Unknown model for Temperature Dependence in "
                    "CKinReact::GetMaxSolubility(). "
                 << "\n";
            cout << " Aborting now!"
                 << "\n";
            exit(0);
        }
    }
    return Csatmax;
}

/**************************************************************************/
/* Yield coefficient of microbe growth, including _drmc_ effects        */
/* 08/2012     CB                                                         */
/**************************************************************************/
double CKinReact::EffectiveYield(long node)
{
    double Yieldeff = Yieldcoefficient;

    if (DormTypeIdx == 0)
    {
        MicrobeData* m_md = NULL;
        m_md = MicrobeData_vector[MicrobeData_idx];
        Yieldeff *= (1 - m_md->G0 / m_md->Gibbs[node]);
        return MRange(0.0, Yieldeff, Yieldcoefficient);
    }

    return Yieldeff;
}

/**************************************************************************/
/* Dichte einer Phase bestimmen                                        */
/* DS-TBC                                                                 */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/**************************************************************************/

double CKinReact::GetDensity(int comp, long index)
{
    double dens = 0.;
    long group, phase;
    // WW CRFProcess *m_pcs = NULL;

    group = index;  // avoid warning

    cout << comp << " " << cp_vec.size() << " ";
    cout << cp_vec[comp]->transport_phase << endl;

    phase = cp_vec[comp]->transport_phase;

    // Get material properties of element
    group = 0;  // m_pcs->m_msh->ele_vector[index]->GetPatchIndex(); //SB todo

    if (phase == 0)  // component in water phase, return liquid density

        dens = mfp_vector[0]->Density();
    // component in solid phase, return solid phase density
    else if (phase == 1)
        dens = msp_vector[group]->Density();
    else if (phase == 2)  // component in bio phase, return 1 as density
        dens = 1;

    return dens;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: BacteriaGrowth(reaction)                          */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet mikrobielles Wachstum aufgrund Reaktion reaction             */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: reaction - Nummer der Wachstumsreaktion                             */
/*    Concentration - Array aktuelle Konzentrationen aus Gleichungsl�ser  */
/*    sumX - Bakteriendichte f�r Ber�cksichtigung MaxKapazit�t            */
/*                                                                        */
/* Ergebnis:                                                              */
/* R: Growth - Wachstum der zur Reaktion gehoerigen Bakteriengruppe       */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/*                                                                        */
/**************************************************************************/

double CKinReact::BacteriaGrowth(int r, double* c, double sumX, int exclude,
                                 long node)
{
    r = r;  // OK411
    int i, BacteriaNumber, MonodSpecies, InhibitionSpecies, Isotopespecies,
        NumberMonod, NumberInhibition;
    double Growth, BacteriaMass, maxVelocity, maxkap;
    double MonodConcentration, InhibitionConcentration, C;
    double Ctot = 0;  // CB Isotope fractionation
    double MonodOrder;
    double Temperature = 0;
    double theta = 0;
    CKinReactData* m_krd = NULL;
    MicrobeData* m_md = NULL;

    m_krd = KinReactData_vector[0];

    BacteriaNumber = bacteria_number + 1;
    maxVelocity = rateconstant;
    BacteriaMass = c[BacteriaNumber];
    if (T_dependence)
    {
        if (this->T_model == 0)
        {
            // do nothing
        }
        else
        {
            if (REACTINT_vec.size() > 0)  // Get the Temperature
                Temperature = REACTINT_vec[0]->GetTemperature(node);
            if (T_model == 1)
            {  // Ratkowski2
                if (Temperature < T_params[0] || Temperature > T_params[1])
                    maxVelocity = 0;
                else
                {
                    maxVelocity =
                        1.0 / 3600 *
                        pow((T_params[2] * (Temperature - T_params[0])),
                            2.0);  // [b*(T-Tmin)]^2
                    maxVelocity *=
                        (1.0 - exp(T_params[3] *
                                   (Temperature -
                                    T_params[1])));  // *(1-exp(C(T-Tmax)))
                }
            }
            else
            {
                cout << " Unknown model for growth dependent on Temperature in "
                        "CKinReact::BacteriaGrowth(). "
                     << "\n";
                cout << " Aborting now!"
                     << "\n";
                exit(0);
            }
        }
    }

    if (_drmc_)
    {
        m_md = MicrobeData_vector[MicrobeData_idx];
        theta = 1 / (exp((m_md->G0 - m_md->Gibbs[node]) /
                         (m_md->steepness * m_md->G0)) +
                     1);
        if (theta > 1.0)
            theta = 1.0;
        switch (DormTypeIdx)
        {
            case 0:  // GROWTH
                maxVelocity *= theta;
                break;
            case 1:  // DEACTIVATION
                maxVelocity *= (1 - theta);
                break;
            case 2:  // REACTIVATION
                maxVelocity *= theta * m_md->_drmc_level[node];
                break;
            case 3:  // DECAY
                maxVelocity += m_md->_drmc_level[node] * m_md->decayrate;
                break;
            default:
                cout << " Warning in CKinReact::BacteriaGrowth(): Unknown "
                        "_drmc_ reaction type!"
                     << "\n";
                break;
        }
    }

    Growth = maxVelocity * BacteriaMass;

    /* Hemmung durch Bakteriendichte sumX nur bei Wachstum */
    if ((sumX > 1.E-30) && (maxVelocity > 0.) && (specif_cap < 1e-30))
    {  // specif_cap is used to limit growth only to the specific bacteria
        /* Max. Bakterienkapazit�t aus Datenstruktur auslesen */
        maxkap = m_krd->maxBacteriaCapacity;
        Growth *= maxkap / (sumX + maxkap);
    }
    // Limit growth by specific bacteria density and not by total density
    else if (specif_cap > 1e-30)
    {
        Growth *= specif_cap / (BacteriaMass + specif_cap);
    }

    /*  FOR-Schleife �ber vorhandene Monodterme */

    NumberMonod = number_monod;
    for (i = 0; i < NumberMonod; i++)
    {
        /* Möglichkeit zum Weglassen Monodterm für partielle Ableitungen
         * erforderlich */
        if (i != exclude)
        {
            MonodSpecies = monod[i]->speciesnumber;
            MonodConcentration = monod[i]->concentration;
            MonodOrder = monod[i]->order;  // CB higher order Monod terms
            C = c[MonodSpecies + 1];
            Ctot = C;  // standard case, no iso fractionation
            // CB Isotope fractionation: here Ctot=C_light+C_heavy must be
            // passed to fct Monod()
            if ((typeflag_iso_fract == 1) &&
                (monod[i]->isotopecouplenumber >= 0))
            {
                Isotopespecies = monod[i]->isotopecouplenumber;
                Ctot = C + c[Isotopespecies + 1];
            }
            // Growth = Growth * Monod(MonodConcentration,C);  // old
            // formulation without isotope fractionation
            Growth *= Monod(MonodConcentration, C, Ctot,
                            MonodOrder);  // new formulation
            // now multiply by additional Threshhold Term, if present,
            // technically this is the same as an additional Monod term for the
            // same MonodSpecies usually of higher order and lower Monod-Conc =
            // threshConc
            if (monod[i]->threshhold == true)
            {
                MonodConcentration =
                    monod[i]->threshConc;  // CB Threshhold concentration
                MonodOrder =
                    monod[i]->threshOrder;  // CB higher order Monod terms
                Growth *= Monod(MonodConcentration, Ctot, Ctot,
                                MonodOrder);  // C should in any case be total C
            }
        }
    }  //  for NumberMonod

    /*	FOR-Schleife �ber vorhandene Inhibitionsterme */
    NumberInhibition = number_inhibit;
    for (i = 0; i < NumberInhibition; i++)
    {
        InhibitionSpecies = inhibit[i]->speciesnumber;
        InhibitionConcentration = inhibit[i]->concentration;
        // ATTENTION!!!
        // CB 16/10/09 fix of Fe3 inhibition concentration for Brand model,
        // this parameter depends on porosity as in Min3P Fe3 is in solid phase
        // and Inhibition concentrations are expressed in terms of Volume
        // fraction Vol_Fe3/Vol_BulkAquifer [m?m�] while in Geosys, Fe3 is
        // considered an immobile aqueous species and concentrations were
        // converted to [mol/L_water] by C_gs = C_min3p * rho_Fe3 / molweight /
        // porosity for inhibition concentrations, division by porosity is still
        // required
        //    if(inhibit[i]->species.compare("Fe3")==0) {
        //      InhibitionConcentration *= 1/GetPhaseVolumeAtNode(currentnode,
        //      1, 0);
        //    } // CB  further changes in derivs (1), jacbn (2), class
        //    CKinReact{}
        C = c[InhibitionSpecies + 1];
        Growth *= Inhibition(InhibitionConcentration, C);
    }

    return Growth;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Monod(MC, C)                                      */
/* DS-TBC                                                                 */
/*                                                                        */
/* Aufgabe:                                                               */
/* Berechnet modifizierten Monod-Term                                     */
/*                                                                        */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)                 */
/* E: MC - MonodConcentration                                             */
/*    C  - Concentration of substance                                     */
/*                                                                        */
/* Ergebnis:                                                              */
/* A: Monod - Wert des Monodterms                                         */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 05/2003     DS         Erste Version                                   */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                        */
/**************************************************************************/
double CKinReact::Monod(double MC, double C, double Ctot, double order)
// double CKinReact::Monod ( double MC, double C )
{
    double Monod;
    // CKinReactData *m_krd = KinReactData_vector[0];
    if (C > 0.)
        // Monod=C/(MC+C);        // normaler Monodterm
        Monod = pow((C / (MC + Ctot)),
                    order);  // CB higher order Monod terms --> factor order for
                             // partial derivatives
    else
        /* linearisierter Term fuer c<=0, with very small slope due to high
         * concentrations */
        Monod = pow(((C / 1000.) / MC), order);  // CB higher order Monod terms

    return Monod;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Inhibition(IC, C)
   DS-TBC                                                                 */
/* Aufgabe:
   Berechnet modifizierten Inhibitions-Term
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: IC - InhibitionConcentration
   C  - Concentration of substance
 */
/* Ergebnis:
   A: Inhibition - Wert des Inhibitionsterms
 */
/* Programmaenderungen:
   05/2003     DS         Erste Version									  */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                          */
/**************************************************************************/

double CKinReact::Inhibition(double IC, double C)
{
    double Inhibition;
    //  CKinReactData *m_krd = KinReactData_vector[0];

    if (C > 0.)
        /* normaler Inhibitionsterm */
        Inhibition = IC / (IC + C);
    else
        /* linearisierter Term fuer C<=0 */
        Inhibition = 1.0;  // Inhibition=1.-C/IC;   CB changed due to stimulance
                           // of growth for neg conc.

    return Inhibition;
}

/*************************************************************************************/
/* Routine f�r Bulirsch-Stoer Gleichungsl�ser (steife ODEs) */
/* DS-TBC */
/*                                                                                   */
/* Berechnung der Jacobi-Matrix, partielle Ableitungen von dC/dt aus Funktion */
/* derive nach */
/* a) der Zeit d2C/d2t = dfdt[] = 0. f�r mikrobiellen Abbau */
/* b) den Konzentrationen d2C/dt*dC = dfdc[][] */
/*                                                                                   */
/* Input: */
/* t = aktuelle Zeit */
/* c[1..Number_of_Components] = aktuelle Konzentration der Substanzen in Zelle
 */
/* n =! Number_of_Components */
/*                                                                                   */
/* Output: */
/* dfdt[1..Number_of_Components] = aktuelle Ableitung der Funktion nach der Zeit
 */
/* dfdc[1..Number_of_Components][1..Number_of_Components] = aktuelle Ableitung
 */
/*             der Funktion nach den Konzentrationen */
/*                                                                                   */
/* Programmaenderungen: */
/* 05/2003     DS         Erste Version */
/* 02/2006     SB         Introduced new C++ concept, Data structures     */
/*                                                                                   */
/*************************************************************************************/

void jacobn(double t, double c[], double dfdt[], double** dfdc, int n,
            long node)
{
    t = t;  // OK411
    int i, j, r, nreactions, BacteriaNumber, NumberMonod, MonodSpecies,
        NumberInhibition, InhibitionSpecies;
    int Sp1, Sp2, SpX, surfaceID = -1, surfaceID2, blob;
    double maxkap, speccap, BacteriaMass, sumX = 0., BacGrowth, maxVelocity,
                                          *d2X_dtdS;
    double CMonodSpecies, MonodConcentration, CInhibitionSpecies,
        InhibitionConcentration, Yield;
    double porosity1, porosity2, exch, kd, density1, saturation2, kadsorb;
    double kdesorb, totalSurface, adsorb, exponent, parameter;
    // SBtodo	double occupiedSurface[maxSurfaces+1];
    vector<double> occupiedSurface;
    int IsotopeSpecies;
    double Ciso;
    double MonodOrder;
    double ThreshOrder;
    double ThreshConc;
    double Yieldcff = 0;

    double ratefact = 1;

    int MineralNumber, MineralSpecies, MechSpecies, no_mech, no_mechspec;
    double MineralGrowth, CMineralSpecies, StoechMinSpecies, CMechSpecies,
        ExpoMechSpecies, omega, dXdt_Mech;
    double theta, eta, eta2, area, Betrag, BetragArg;
    CKinReact *m_kr = NULL, *m_kr1 = NULL;  // OK411 *m_kr2=NULL;
    CKinBlob* m_kb = NULL;
    CKinReactData* m_krd = NULL;
    MinkinMech* m_mech = NULL;
    //	CMediumProperties *m_mat_mp = NULL;
    double foc;
    double dt = 0;

    m_krd = KinReactData_vector[0];
    // what's the current time step size?
    if (m_krd->NumberMineralkinetics > 0)
    {
        CTimeDiscretization* m_tim = NULL;
        m_tim = time_vector[0];
        dt = m_tim->CalcTimeStep();
    }

    // if(m_krd->debugoutflag)
    //  m_krd->debugoutstr << " jacobn" << "\n" << flush;

    /* Hilfsvektor f�r partielle Ableitung des Bakterienwachstums nach Species S
     */
    d2X_dtdS = dvector(1, n);

    /* weitere Ableitungen nach t dfdt[] alle null */
    /* Ableitungen nach c dfdc[][] werden inkrementiv berechnet, also erst alles
     * null setzen */
    /* ACHTUNG, unterschiedliche Indizierung der Arrays, c[1..n]
     * BioDegradation[0..n-1] */
    for (i = 0; i < n; i++)
    {
        dfdt[i + 1] = 0.;
        for (j = 0; j < n; j++)
            dfdc[i + 1][j + 1] = 0.;
    }

    /* calculate present bacteria capacity */
    maxkap = m_krd->maxBacteriaCapacity;
    // F�r Berechnung der Ableitungen f�r den Fall dass eine maximale Kapazit�t
    // ber�cksichtigt werden muss Muss sein, weil Ableitungen h�here Ordnung
    // haben (Bakterienmasse steckt auch in Kapazit�tsgleichung)
    sumX = 0.;  // added CB
    if (maxkap > 1.E-30)
    {
        for (i = 0; i < n; i++)
            if (m_krd->is_a_bacterium[i])
            {
                BacteriaMass = c[i + 1];
                sumX += BacteriaMass;
            }
    }

    /* Anzahl der mikrobiellen Reaktionen aus Datenstruktur auslesen */
    nreactions = m_krd->NumberReactions;

    /* loop over reactions dX/dt= nymax * X * monodterms * inhibitionterms */
    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;

        if (m_kr->typeflag_monod)
        {
            BacteriaNumber = m_kr->bacteria_number + 1;
            BacteriaMass = c[BacteriaNumber];
            speccap = m_kr->specif_cap;

            if (BacteriaMass > 1.E-40)
            {
                /* Ableitungen werden aus dX/dt = BacGrowth berechnet */
                // sumX is different for case with (>0) or without (==0) maxkap

                m_kr->currentnode =
                    node;  // CB This is eclusively for Brand model to allow
                           // porosity in Inhibition
                // constant calculation
                Yieldcff = m_kr->EffectiveYield(node);
                BacGrowth = m_kr->BacteriaGrowth(r, c, sumX, -1, node);
                for (i = 0; i < n; i++)
                    d2X_dtdS[i + 1] = 0.;

                // Berechnung der Bakterien Ableitungen für den Fall dass eine
                // maximale Kapazität berücksichtigt werden muss Muss sein, weil
                // Ableitungen höhere Ordnung haben (Bakterienmasse steckt auch
                // in Kapazitätsgleichung) This is for the case that growth of
                // any bacteria is limited by total bacteria density sumX &
                // maxkap
                if ((maxkap > 1.E-30) &&
                    (speccap <
                     1e-30))  // no bacteria specific capacity term defined
                {
                    maxVelocity = m_kr->rateconstant;
                    /* Wachstumsterm, ber�cksichtige Kapazit�tsterm */
                    if (maxVelocity > 1.E-30)
                    {
                        // Erst Berechnen der Ableitungen des Bakterienwachstums
                        // nach allen anderen Substanzen
                        //   Ableitung nach der Bakterienmasse (mit
                        //   Ber�cksichtigung Kapazit�tsterm) d2Xi / dt*dXi =
                        //   BacGrowth *(sumx+maxkap-Xi) / (Xi*(sumx+maxkap))
                        d2X_dtdS[BacteriaNumber] =
                            BacGrowth * (sumX + maxkap - BacteriaMass) /
                            (BacteriaMass * (sumX + maxkap)) *
                            Yieldcff;  // CB 08/12 added yield coeff
                        for (i = 0; i < n; i++)
                        {
                            if (m_krd->is_a_bacterium[i] &&
                                (i + 1 != BacteriaNumber))
                                // Ableitung nach den anderen Bakterienmassen im
                                // Kapazit�tsterm
                                //   d2Xi / dt*dXj = BacGrowth / -(sumx+maxkap)
                                d2X_dtdS[i + 1] =
                                    BacGrowth / -(sumX + maxkap) *
                                    Yieldcff;  // CB 08/12 added yield coeff
                        }
                    } /* Sterbeterm, grunds�tzlich keine Ber�cksichtigung des
                         Kapazit�tsterms */
                    else
                    { /* Sterbeterm, grundsätzlich keine Berücksichtigung des
                         Kapazitätsterms */
                        /* d2Xi / dt*dXi = BacGrowth / Xi */
                        d2X_dtdS[BacteriaNumber] = BacGrowth / BacteriaMass;
                        /* d2Xi / dt*dXj = 0 for decay */
                    }
                }
                // This is for the case that growth of a specific bacteria is
                // limited by the SPECIFIC bacteria density BacteriaMass &
                // specif_cap
                else if (speccap > 1e-30)
                {  // bacteria specific capacity term defined
                    maxVelocity = m_kr->rateconstant;
                    /* Wachstumsterm, berücksichtige Kapazitätsterm */
                    if (maxVelocity > 1.E-30)
                    {
                        // Erst Berechnen der Ableitungen des Bakterienwachstums
                        // nach allen anderen Substanzen
                        //   Ableitung nach der Bakterienmasse (mit
                        //   Berücksichtigung d. spezif. Kapazitätsterms) d2Xi /
                        //   dt*dXi = BacGrowth *(specif_cap) /
                        //   (Xi*(Xi+specif_cap))
                        d2X_dtdS[BacteriaNumber] =
                            BacGrowth * (speccap) /
                            (BacteriaMass * (BacteriaMass + speccap)) *
                            Yieldcff;  // CB 08/12 added yield coeff
                        // Ableitung nach den anderen Bakterienmassen im
                        // Kapazitätsterm non existant
                    }
                    else
                    { /* Sterbeterm, grundsätzlich keine Berücksichtigung des
                         Kapazitätsterms */
                        /* d2Xi / dt*dXi = BacGrowth / Xi */
                        d2X_dtdS[BacteriaNumber] = BacGrowth / BacteriaMass;
                        /* d2Xi / dt*dXj = 0 for decay */
                    }
                }
                // Berechnung der Bakterien Ableitungen für den Fall dass KEINE
                // maximale Kapazität berücksichtigt werden muss without
                // capacity term, this is identical for growth and decay term;
                // but not if yield coeff is considered (CB)
                else
                {  // maxkap = 0
                    /* d2Xi / dt*dXi = BacGrowth / Xi */
                    d2X_dtdS[BacteriaNumber] =
                        BacGrowth / BacteriaMass *
                        Yieldcff;  // CB 08/12 added yield coeff = 1 for decay;
                }

                /* Schleife f�r Ableitungen nach Substanzen in Monodtermen;
                 * unabh�ngig von maxkap */
                NumberMonod = m_kr->number_monod;
                for (i = 0; i < NumberMonod; i++)
                {
                    // d2X / dt*dS_j =      S_j = monod-species
                    //   S_j may be ZERO or below !
                    MonodSpecies = m_kr->monod[i]->speciesnumber + 1;
                    MonodConcentration = m_kr->monod[i]->concentration;
                    MonodOrder = m_kr->monod[i]->order;
                    CMonodSpecies = c[MonodSpecies];
                    if (CMonodSpecies > 1.E-20)
                    {
                        // S_j > 0, normal Monod Term used
                        //   divide BacGrowth through Monod-Term of spec. j
                        //   and multiplicate with partial derivate of
                        //   Monod-Term

                        // In case of isotope fractionation of substrate Si,Sj
                        // (i,j=l,h ; light,heavy)
                        // - the partial derivative d2X/dtdSi is different:
                        //   d2X_dtdS[Si] = BacGrowth *
                        //   (MonodConcentration+CisotopePartner)
                        //           / CMonodSpecies /
                        //           (MonodConcentration+CMonodSpecies+CisotopePartner)
                        // - an additional partial derivative d2X/dtdSj with
                        // respect to the isotope partner Sj appears
                        //   and must be accounted for and added in the vector
                        //   d2X_dtdS[...] for the current reaction:
                        //   d2X_dtdS[Sj] = BacGrowth * (-1) /
                        //   (MonodConcentration+CMonodSpecies+CisotopePartner)...

                        if ((m_kr->monod[i]->species == m_kr->Isotope_heavy) ||
                            (m_kr->monod[i]->species == m_kr->Isotope_light))
                        {
                            IsotopeSpecies =
                                m_kr->monod[i]->isotopecouplenumber + 1;
                            Ciso = c[IsotopeSpecies];
                            // this is the term for the Monod-Species
                            d2X_dtdS[MonodSpecies] =
                                BacGrowth * MonodOrder *
                                (MonodConcentration + Ciso) / CMonodSpecies /
                                (MonodConcentration + CMonodSpecies +
                                 Ciso);  // no isofrac
                            // now get the partial derivative d2x/dtdSj with
                            // respect to isotope partner
                            d2X_dtdS[IsotopeSpecies] =
                                BacGrowth * MonodOrder * (-1) /
                                (MonodConcentration + CMonodSpecies + Ciso);
                            // If a threshhold term exists for the Monod
                            // species, the partial derivative is modified for
                            // this species
                            //   - the ODE with Monod and threshold term for
                            //   Monod species C and Isotope fractionation is:
                            //     dX/dt = my*R * [Cl/(Cl+Ch+K)]^n *
                            //     [(Cl+Ch)/(Cl+Ch+T)]^m
                            //      - with n and K the Order and
                            //      Monod-concentration of the Monod term
                            //      - with m and T the Order and
                            //      Threshhold-concentration of the threshhold
                            //      term
                            //   - The partial derivative with respect to Cl
                            //   (taking into account Division by the Monod and
                            //   Threshhold term) is
                            //     d2X/dtdCl = my*R * [n*(Ch+K)/Cl/(Cl+Ch+K) +
                            //     p*T/(Cl+Ch)/(Cl+Ch+T)]
                            //   - The partial derivative with respect to Ch
                            //   (taking into account Division by the Monod and
                            //   Threshhold term) is
                            //     d2X/dtdCh = my*R * [-n/(Cl+Ch+K) +
                            //     p*T/(Cl+Ch)/(Cl+Ch+T)]
                            //   - The latter term is identical for both derivs
                            //   hence must be added to the previously
                            //   calculated first terms
                            if (m_kr->monod[i]->threshhold == true)
                            {
                                ThreshConc = m_kr->monod[i]->threshConc;
                                ThreshOrder = m_kr->monod[i]->threshOrder;
                                d2X_dtdS[MonodSpecies] +=
                                    BacGrowth * ThreshOrder * ThreshConc /
                                    (CMonodSpecies + Ciso) /
                                    (ThreshConc + CMonodSpecies +
                                     Ciso);  // isofrac
                                d2X_dtdS[IsotopeSpecies] +=
                                    BacGrowth * ThreshOrder * ThreshConc /
                                    (CMonodSpecies + Ciso) /
                                    (ThreshConc + CMonodSpecies +
                                     Ciso);  // isofrac
                            }
                        }
                        else  // no isofrac standard case:
                              // d2X_dtdS[MonodSpecies] = BacGrowth *
                              // MonodConcentration /
                        // CMonodSpecies / (MonodConcentration+CMonodSpecies);
                        // // no isofrac
                        {
                            d2X_dtdS[MonodSpecies] =
                                BacGrowth * MonodOrder * MonodConcentration /
                                CMonodSpecies /
                                (MonodConcentration +
                                 CMonodSpecies);  // no isofrac
                            // If a threshhold term exists for the Monod
                            // species, the partial derivative is modified for
                            // this species
                            //   - the ODE with Monod and threshold term for
                            //   Monod species C is:
                            //     dX/dt = my*R * [C/(C+K)]^n * [C/(C+T)]^m
                            //      - with n and K the Order and
                            //      Monod-concentration of the Monod term
                            //      - with m and T the Order and
                            //      Threshhold-concentration of the threshhold
                            //      term
                            //   - The partial derivative with respect to C
                            //   (taking into account Division by the Monod and
                            //   Threshhold term) is
                            //     d2X/dtdC = my*R * [n*K/C/(C+K) + p*T/C/(C+T)]
                            //   - The latter term hence must be added to the
                            //   previously calculated first term
                            if (m_kr->monod[i]->threshhold == true)
                            {
                                ThreshConc = m_kr->monod[i]->threshConc;
                                ThreshOrder = m_kr->monod[i]->threshOrder;
                                d2X_dtdS[MonodSpecies] +=
                                    BacGrowth * ThreshOrder * ThreshConc /
                                    CMonodSpecies /
                                    (ThreshConc + CMonodSpecies);  // no isofrac
                            }
                        }
                    }
                    else if (CMonodSpecies < -1.E-20)
                    {
                        /* S_j << 0, linear Monod Term used */
                        // d2X_dtdS[MonodSpecies] = BacGrowth / CMonodSpecies ;
                        d2X_dtdS[MonodSpecies] =
                            BacGrowth * MonodOrder / CMonodSpecies /
                            1000;  // Changed monod term with smaller slope CB
                        if (m_kr->monod[i]->threshhold == true)
                        {
                            ThreshConc = m_kr->monod[i]->threshConc;
                            ThreshOrder = m_kr->monod[i]->threshOrder;
                            d2X_dtdS[MonodSpecies] +=
                                BacGrowth * ThreshOrder * ThreshConc /
                                CMonodSpecies /
                                (ThreshConc + CMonodSpecies);  // no isofrac
                        }
                    }  // Todo CB isofrac special case necessary? Threshhold
                       // terms??
                    else
                    {
                        // S_j near 0 numerically instable
                        //   recompute BacGrowth without S_j
                        //   (hope, that will only sometimes occur)

                        m_kr->currentnode =
                            node;  // CB 19/10/09 This is eclusively for Brand
                                   // model to allow porosity in
                        // Inhibition constant calculation

                        d2X_dtdS[MonodSpecies] =
                            m_kr->BacteriaGrowth(r, c, sumX, MonodSpecies,
                                                 node) /
                            MonodConcentration;
                    }
                }  // for NumberMonod

                /* Schleife f�r Ableitungen nach Substanzen in
                 * Inhibitionstermen, unabh�ngig von maxkap */
                NumberInhibition = m_kr->number_inhibit;
                for (i = 0; i < NumberInhibition; i++)
                {
                    // d2X / dt*dS_j =      S_j = inhibition-species
                    //   S_j may be Zero without any problem
                    InhibitionSpecies = m_kr->inhibit[i]->speciesnumber + 1;
                    InhibitionConcentration = m_kr->inhibit[i]->concentration;

                    // ATTENTION!!!
                    // CB 16/10/09 fix of Fe3 inhibition concentration for Brand
                    // model, this parameter depends on porosity as in Min3P Fe3
                    // is in solid phase and Inhibition concentrations are
                    // expressed in terms of Volume fraction
                    // Vol_Fe3/Vol_BulkAquifer [m?m�]
                    // while in Geosys, Fe3 is considered an immobile aqueous
                    // species and concentrations were converted to
                    // [mol/L_water] by C_gs = C_min3p * rho_Fe3 / molweight /
                    // porosity for inhibition concentrations, division by
                    // porosity is still required
                    //        if(m_kr->inhibit[i]->species.compare("Fe3")==0) {
                    //          InhibitionConcentration *=
                    //          1/m_kr->GetPhaseVolumeAtNode(node, 1, 0);
                    //        } // CB  further changes in derivs (1), jacbn (2),
                    //        class CKinReact{}

                    CInhibitionSpecies = c[InhibitionSpecies];
                    if (CInhibitionSpecies > 0.)
                    {
                        // S_j > 0, normal Inhibition Term used
                        //   divide BacGrowth through Inhibition-Term of spec. j
                        //   and multiplicate with partial derivate of Inhi-Term
                        d2X_dtdS[InhibitionSpecies] =
                            -BacGrowth /
                            (InhibitionConcentration + CInhibitionSpecies);
                    }
                    else
                    {
                        /* S_j <= 0, linear Inhibition Term used */
                        // d2X_dtdS[InhibitionSpecies] = - BacGrowth /
                        // (InhibitionConcentration-CInhibitionSpecies);// CB
                        // changed as in next line
                        d2X_dtdS[InhibitionSpecies] =
                            -BacGrowth /
                            (InhibitionConcentration);  // CB changed due to
                                                        // stimulance of growth
                                                        // for neg conc.
                    }
                }

                /* transfer partial derivatives to dfdc-Array of equation solver
                 */
                if (m_kr->grow)
                {
                    for (i = 0; i < n; i++)
                        /* transfer der berechneten Ableitungen f�r die
                         * Bakteriengruppe */
                        dfdc[BacteriaNumber][i + 1] += d2X_dtdS[i + 1];
                }

                /* Berechnung der Ableitungen f�r die vom Bakteriellen Wachstum
                 * abh�ngigen Substanzen, unabh�ngig von maxkap */
                /* d2S_j / dt*dS_k = yield(j) * d2X/dt*dS_k */
                porosity1 = m_kr->GetReferenceVolume(BacteriaNumber - 1, node);
                for (i = 0; i < n; i++)
                {
                    Yield = m_kr->ProductionStoch[i];
                    if (fabs(Yield) > 1.E-30)
                    {
                        porosity2 = m_kr->GetReferenceVolume(i, node);
                        for (j = 0; j < n; j++)
                        {
                            if (fabs(d2X_dtdS[j + 1]) > 1.E-30)
                                dfdc[i + 1][j + 1] += d2X_dtdS[j + 1] * Yield *
                                                      porosity1 / porosity2;
                        }
                    }
                }

            }  // Ende if BacteriaMass > 1e-40*/
        }      // Ende if type monod
    }          // Ende Schleife �ber nreactions*/

    // Berechnung der Ableitungen der Mineralkinetiken
    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;

        // Mineral kinetics
        if (m_kr->typeflag_mineralkinetics)
        {
            // first, compute quantities that are used several times below
            omega = m_kr->Omega(c, node);
            theta = m_kr->Theta;
            BetragArg = 1 - pow(omega, theta);
            Betrag = fabs(BetragArg);
            // Equilibrium: when rate is zero, division by 0 (1-omega^theta) is
            // a problem, therefore, do only when not in equilibrium
            if (Betrag > 0)
            {
                // get quantities
                eta = m_kr->Eta;
                eta2 = m_kr->precipexponent;
                MineralNumber = m_kr->mineral_number + 1;
                // This is the mineral growth rate, i.e. the basic derivative
                // dM/dt
                MineralGrowth = m_kr->MinRate(c, node, dt, omega);
                // initialize vector before calculating partial derivs
                for (i = 0; i < n; i++)
                    d2X_dtdS[i + 1] = 0.0;
                //(1)
                //
                // d2M_dtdCmin derivative wrt mineral species in IAP of rate law
                // d2M_dtdCmin = - K * Am * eta * theta * Stoech_i / C_i *
                // |1-Om^theta|^(eta-1) * Om^theta
                //
                // This is complete, if the mineral species is NOT a mech
                // species at the same time; otherwise, additional terms will
                // apear in 2nd deriv, these are automaticaly computed below in
                // loop over mechterms and mechspecies (2) Anyway, we calculate
                // these partial derivatives based on dM/dt = MineralGrowth.
                // StoechCoeff sign in factor2 must be changed because it is
                // defined negative for species which in the IAP are forming the
                // mineral (there, they are positive, See fct. Omega)

                // todo: check root function derivative

                for (i = 0; i < int(m_kr->minSpeciesIdx.size()); i++)
                {
                    MineralSpecies = m_kr->minSpeciesIdx[i] + 1;
                    CMineralSpecies = c[MineralSpecies];
                    StoechMinSpecies =
                        m_kr->ProductionStoch[MineralSpecies - 1];
                    d2X_dtdS[MineralSpecies] =
                        -fabs(MineralGrowth);  // get dM/dt and make it negative
                    d2X_dtdS[MineralSpecies] /=
                        Betrag;  // reduce power ^eta -> ^(eta-1)
                    d2X_dtdS[MineralSpecies] *=
                        eta * theta *
                        (-StoechMinSpecies);  // factor1 from derivative
                    if (BetragArg <= -1)      // strong supersaturation
                        d2X_dtdS[MineralSpecies] *=
                            1 / eta2;  //   --> slow down precipitation rate by
                                       //   root fct. (default = 1)
                    if (CMineralSpecies > 1e-25 ||
                        CMineralSpecies < -1e-25)  // factor2 from derivative
                        d2X_dtdS[MineralSpecies] *= 1 / CMineralSpecies;
                    else  // for almost zero concentration
                        d2X_dtdS[MineralSpecies] *= 1 / 1e-25;
                    d2X_dtdS[MineralSpecies] *=
                        pow(omega, theta);  // factor3 from derivative
                }

                //(2)
                //
                // d2M_dtdCmech derivative wrt mechanism species in rate
                // constant of rate law d2M_dtdCmech = 0 + 0 + ... + K_mech*A*
                // |1-Om^theta|^eta * MechExpo_i/Cmech_i
                //
                // calculate these partial derivatives here without using
                // previously calculated dM/dt as only the a term containing
                // Cmech remains from sum after taking derivative. This partial
                // derviatives adds up to a mineral species p.d. from (1), when
                // it is a mech species at the same time; however, no special
                // treatment is required first, get some quantities

                // todo: check root function derivative

                no_mech = int(m_kr->mechvec.size());

                // get Surface Area
                if (m_kr->Am_constant)
                    area = m_kr->Am[0];
                else
                    area = m_kr->Am[node];
                porosity1 = m_kr->GetReferenceVolume(MineralNumber - 1, node);
                for (i = 0; i < no_mech; i++)
                {
                    m_mech = m_kr->mechvec[i];
                    no_mechspec = m_mech->no_mechSpec;
                    // the next is identical for all MechSpecies partial derivs
                    dXdt_Mech = pow(Betrag, eta);  // = |1-Om^theta|^eta
                    if (BetragArg <= -1)           // strong supersaturation
                        dXdt_Mech = pow(
                            dXdt_Mech,
                            1 / eta2);  //   --> slow down precipitation rate by
                                        //   root function (default = 1)
                    dXdt_Mech *= m_mech->Mech(c, node) * area /
                                 porosity1;  //   *K_mech*A0/(1-n)  ;  unit: -->
                                             //   mol/s/m³_solid
                    if (BetragArg <= 0)      // supersaturation
                        dXdt_Mech *=
                            m_kr->precipfactor;  //   --> slow down
                                                 //   precipitation rate by
                                                 //   factor (default = 1)
                    for (j = 0; j < no_mechspec; j++)
                    {
                        MechSpecies = m_mech->mechSpeciesIdx[j] + 1;
                        CMechSpecies = c[MechSpecies];
                        // pH!!!
                        // if(m_mech->mechSpeciesNames[j].compare("pH")==0) //
                        // pH is -log(gamma*H+)
                        //  CMechSpecies  = pow(10,-CMechSpecies);
                        ExpoMechSpecies = m_mech->mechSpeciesExpo[j];
                        // the next is specific to MechSpecies_j p. deriv
                        // but severeal mechanisms_i may contain the same mech
                        // species, and it may ALSO be a mineral species from
                        // (1) at the same time therefore += results
                        if (CMechSpecies > 1e-25 || CMechSpecies < -1e-25)
                            d2X_dtdS[MechSpecies] +=
                                dXdt_Mech * ExpoMechSpecies /
                                CMechSpecies;  // * MechExpo_i/Cmech_i
                        else
                            d2X_dtdS[MechSpecies] +=
                                dXdt_Mech * ExpoMechSpecies /
                                1e-25;  // * MechExpo_i/Cmech_i
                    }
                }

                //(3)
                //
                // d2M_dtdM derivative wrt mineral mass is zero always,
                // d2M_dtdM = 0
                //
                // ... as activity of mineral in IAP of rate law is 1 always per
                // definition
                d2X_dtdS[MineralNumber] = 0.0;

                // now transfer partial derivatives to dfdc-Array of equation
                // solver rows i of J contain a single species' derivatives wrt
                // all other species d2Ci/dtdCj
                for (j = 0; j < n; j++)  // for the Mineral, first row of J
                    dfdc[MineralNumber][j + 1] +=
                        d2X_dtdS[j + 1];  // unit: --> mol/s/m³_solid

                //(4)
                //
                // derivatives of produced or consumed species, the Mineral
                // species d2S_j / dt*dS_k = yield(j) * d2X/dt*dS_k
                //
                // These are simply coupled to p.d. of Mineral wtr. to all other
                // species d2X_dtdS via the stoechiometric coefficients;
                // Attention: this MUST be zero for the mineral itself, as
                // otherwise, the p.d. is added here again
                porosity1 = m_kr->GetReferenceVolume(MineralNumber - 1, node);
                for (i = 0; i < n; i++)
                {  // these are the rows i
                    Yield =
                        m_kr->ProductionStoch[i];  // = 0 for Mineral (must), <0
                                                   // for min species consumed,
                                                   // >0 for
                    // species produced by precipitation
                    if (fabs(Yield) > 1.E-30)
                    {
                        porosity2 = m_kr->GetReferenceVolume(i, node);
                        for (j = 0; j < n; j++)
                        {  // these are the j column entries of a row i, Yi is
                           // the same in 1 row
                            if (fabs(d2X_dtdS[j + 1]) > 1.E-30)
                                // unit conversion: --> mol/s/m³_solid *
                                // m³sol/m³Aq * (m³w/m³Aq)^-1 = mol/s/m³w
                                dfdc[i + 1][j + 1] += d2X_dtdS[j + 1] * Yield *
                                                      porosity1 / porosity2;
                        }
                    }
                }

            }  // if (fabs(1-pow(omega, theta)) > 0) --> no equiliubrium
        }      // if(m_kr->typeflag_mineralkinetics)
    }
    /**********************************************************************************************/
    /* Berechnung der Ableitungen der Austauschprozesse */
    // calculate already occupied surfaces first
    if (m_krd->NumberLangmuir > 0)
    {
        // Initialise Surfaces for langmuir isotherms
        for (i = 0; i < m_krd->maxSurfaces; i++)
            occupiedSurface.push_back(0.0);

        for (r = 0; r < nreactions; r++)
        {
            m_kr = KinReact_vector[r];
            // CB new reaction switch for individual reactions
            if (m_kr->switched_off_node.size() > 0)
                if (m_kr->switched_off_node[node] == true)
                    continue;
            if ((m_kr->getType().compare("exchange") == 0) &&
                (m_kr->typeflag_exchange_langmuir))
            {
                Sp1 = m_kr->ex_species[0] + 1;
                surfaceID = m_kr->exSurfaceID;
                occupiedSurface[surfaceID] += c[Sp1];
            }
        }
    }  // if NumberLangmuir > 0

    /* Berechnung der Ableitungen der Austauschprozesse */
    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;

        if (m_kr->typeflag_exchange)
        {
            /* linearer Austausch mit kd */
            if (m_kr->typeflag_exchange_linear)
            {
                // Matrix
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                // geloest
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
                ratefact = m_kr->ex_param[2];

                exch = m_kr->ex_param[0];
                kd = m_kr->ex_param[1];

                if (fabs(kd) < MKleinsteZahl)
                {
                    // slow down desorption rate by constant factor
                    if ((c[Sp2] - c[Sp1]) < 0)
                        exch *= ratefact;

                    // no kd, exchange between two species in solution
                    dfdc[Sp1][Sp1] += -exch / porosity1;
                    dfdc[Sp1][Sp2] += exch / porosity1;
                    dfdc[Sp2][Sp1] += exch / porosity2;
                    dfdc[Sp2][Sp2] += -exch / porosity2;
                }
                else
                {
                    // with kd, exchange between matrix (mol/kg) and solution
                    // (mol/l)
                    foc = m_krd->node_foc[node];
                    if (foc > MKleinsteZahl)
                        kd = kd * foc;
                    // else
                    //  kd = 0;
                    // slow down desorption rate by constant factor
                    if ((kd * c[Sp2] - c[Sp1]) < 0)
                        exch *= ratefact;
                    dfdc[Sp1][Sp1] += -exch;
                    dfdc[Sp1][Sp2] += exch * kd;
                    dfdc[Sp2][Sp1] += exch * porosity1 / porosity2 * density1;
                    dfdc[Sp2][Sp2] +=
                        -exch * kd * porosity1 / porosity2 * density1;
                }
            }  // linear
            /* Langmuir Kinetik */
            if (m_kr->typeflag_exchange_langmuir)
            {
                // calculated already occupied surfaces above
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
                kadsorb = m_kr->ex_param[0];
                kdesorb = m_kr->ex_param[1];
                // SB_langmuir		surfaceID	= m_kr->exSurfaceID;
                // //Exchange.Langmuir[i].SurfaceID;
                totalSurface = m_krd->exSurface[m_kr->exSurfaceID];
                //      occupiedSurface is calculated just above
                adsorb = kadsorb * (totalSurface - occupiedSurface[surfaceID]);
                dfdc[Sp1][Sp1] += -kdesorb;
                dfdc[Sp2][Sp1] += kdesorb * porosity1 / porosity2 * density1;
                dfdc[Sp1][Sp2] += adsorb;
                dfdc[Sp2][Sp2] += -adsorb * porosity1 / porosity2 * density1;
                // additional derivatives due to occupied surface
                for (j = 0; j < nreactions; j++)
                {
                    m_kr1 = KinReact_vector[j];
                    if (m_kr1->getType().compare("exchange") == 0)
                        if (m_kr1->typeflag_exchange_langmuir)
                        {
                            SpX = m_kr1->ex_species[0] + 1;
                            surfaceID2 = m_kr1->exSurfaceID;
                            if (surfaceID == surfaceID2)
                            {
                                dfdc[Sp1][SpX] += -kadsorb * c[Sp2];
                                dfdc[Sp2][SpX] += kadsorb * c[Sp2] * porosity1 /
                                                  porosity2 * density1;
                            }
                        }
                }
            }  // end if langmuir

            /* Freundlich Kinetik */
            if (m_kr->typeflag_exchange_freundlich)
            {
                Sp1 = m_kr->ex_species[0] + 1;
                porosity1 = m_kr->GetReferenceVolume(Sp1 - 1, node);
                density1 = m_kr->GetDensity(Sp1 - 1, node);
                Sp2 = m_kr->ex_species[1] + 1;
                porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
                exponent = m_kr->ex_param[2];
                parameter = m_kr->ex_param[1];
                exch = m_kr->ex_param[0];
                ratefact = m_kr->ex_param[3];

                // if both species are in the same phase, density multiplication
                // below is not required
                if (cp_vec[Sp1 - 1]->transport_phase ==
                    cp_vec[Sp2 - 1]->transport_phase)
                    density1 = 1;

                // slow down DEsorption rate by constant factor
                if ((parameter * pow(c[Sp2], exponent) - c[Sp1]) < 0)
                    exch *= ratefact;

                if (c[Sp2] > residual)
                    // no linearisation required
                    adsorb = exch * parameter * exponent *
                             pow(c[Sp2], (exponent - 1.0));
                else
                    // linearisation required due to instability of c^x if
                    // c<residual
                    adsorb = exch * parameter * exponent *
                             pow(residual, (exponent - 1.0));

                dfdc[Sp1][Sp1] += -exch;
                dfdc[Sp2][Sp1] += exch * porosity1 / porosity2 * density1;
                dfdc[Sp1][Sp2] += adsorb;
                dfdc[Sp2][Sp2] += -adsorb * porosity1 / porosity2 * density1;
            }  // end freundlich
        }      // end if exchange
    }          // end loop over reactions

    //#ds
    /**********************************************************************************************/
    /* NAPL-dissolution */
    /**********************************************************************************************/
    for (r = 0; r < nreactions; r++)
    {
        m_kr = KinReact_vector[r];
        // CB new reaction switch for individual reactions
        if (m_kr->switched_off_node.size() > 0)
            if (m_kr->switched_off_node[node] == true)
                continue;

        if (m_kr->typeflag_napldissolution)
        {
            /* NAPL-L�sung */
            Sp1 = m_kr->ex_species[0] +
                  1;  // Exchange.Linear[i].Species1; should be NAPL
            //		porosity1	= m_kr->GetReferenceVolume(Sp1-1,node);
            blob = m_kr->blob_ID;
            m_kb = KinBlob_vector[blob];  // pointer to blob-properties set in
                                          // the reaction r
            Sp2 = m_kr->ex_species[1] +
                  1;  // Exchange.Linear[i].Species2; should be dissolved
            // CB this includes the saturation
            porosity2 = m_kr->GetReferenceVolume(Sp2 - 1, node);
            //#ds TODO	    saturation2 = ??
            saturation2 = 1.;
            exponent = m_kr->rateorder;
            exch = m_kb->Masstransfer_K[node] * m_kb->Interfacial_area[node];
            exch *= exponent *
                    pow((c[Sp2]),
                        (exponent - 1.0));  // CB introduced rateorder 06/2012

            // Remark: In function derivs the dissolution velocity can be
            // reduced, if NAPL concentration falls below zero.
            //         The partial derivatives should change, too, but this is
            //         not considered here. However, it works fine

            // Case 1 : reduced rate in derivs
            if ((m_kb->Masstransfer_K[node] * m_kb->Interfacial_area[node] *
                 pow((m_kr->Current_Csat[node] - c[Sp2]), exponent) * dt) >
                c[Sp1])
            {
                //		dfdc[Sp1][Sp2] = 0     derivatives for
                // water-concentration is zero in case of reduced rate in
                // derivs
                //		dfdc[Sp2][Sp2] = 0

                if ((m_kr->Current_Csat[node] < c[Sp2]) ||
                    (c[Sp1] > MKleinsteZahl))
                {  // no dissolution or NAPL mass present
                    // dCNAPL/dt = CNAPL/dt
                    // d2CNAPL / dt dCNAPL = 1/dt
                    // d2Cmob  / dt dCmob = -k*A/n/Sw
                    dfdc[Sp1][Sp2] += -1 / dt;
                    dfdc[Sp2][Sp2] += 1 / dt / porosity2 / saturation2;
                }
            }
            // Case 2 : original rate in derivs
            else
            {
                //		dfdc[Sp1][Sp1] = 0     derivatives for
                // NAPL-concentration always zero 		dfdc[Sp2][Sp1] = 0

                // if ( (m_kr->current_Csat < c[Sp2]) || (c[Sp1] >
                // MKleinsteZahl) ) {         // no dissolution or NAPL mass
                // present
                if ((m_kr->Current_Csat[node] < c[Sp2]) ||
                    (c[Sp1] >
                     MKleinsteZahl))  // no dissolution or NAPL mass present
                {
                    // d2CNAPL / dt dCmob = k*A
                    // d2Cmob  / dt dCmob = -k*A/n/Sw
                    dfdc[Sp1][Sp2] += exch;
                    dfdc[Sp2][Sp2] += -exch / porosity2 / saturation2;
                }
            }
        }  // NAPL-dissolution
    }      // loop over reactions r
    //#ds

    free_dvector(d2X_dtdS, 1, n);
}

/**************************************************************************
   Reaction-Method:
   Task: Reaction class test output function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
   02/2006 SB Adapted to new FEM structure
**************************************************************************/
void CKinReact::TestWrite(void)
{
    int i, length, flag = 0;

    // Write Keyword
    cout << "888888888888888888888888888888888888888888888888888888888888888888"
            "8888888888888888888888888888888888888"
         << "\n";
    cout << " Test Output "
         << "\n";
    cout << "888888888888888888888888888888888888888888888888888888888888888888"
            "8888888888888888888888888888888888888"
         << "\n";
    cout << "#REACTION"
         << "\n";
    // Name of reaction
    cout << "$NAME"
         << "\n"
         << name << "\n";
    // Type of reaction
    cout << "$TYPE"
         << "\n"
         << type << "\n";
    // bacteria name
    cout << "$BACTERIANAME"
         << "\n"
         << bacteria_name << "\n";
    // ReactionEquation
    cout << "$EQUATION"
         << "\n";
    for (i = 0; i < number_reactionpartner; i++)
    {
        if (stochmet[i] < 0.0)  // left side of equation
        {
            if (i == 0)
                cout << " " << fabs(stochmet[i]) << " " << reactionpartner[i];
            else
                cout << " + " << fabs(stochmet[i]) << " " << reactionpartner[i];
        }
        if (stochmet[i] > 0 && (flag > 0))  // remaining right hand side
            cout << " + " << fabs(stochmet[i]) << " " << reactionpartner[i];
        if (stochmet[i] > 0 &&
            (flag == 0))  // " = " Sign and first term on right hand side
        {
            cout << " = " << fabs(stochmet[i]) << " " << reactionpartner[i];
            flag = 1;
        }
    }
    cout << "\n";
    // Rateconstant and order
    cout << "$RATEKONSTANT"
         << "\n"
         << rateconstant << "   " << rateorder << "\n";
    cout << "$GROWTH"
         << "\n"
         << grow << "\n";
    // Monod terms
    cout << "$MONODTERMS"
         << "\n"
         << number_monod << "\n";
    for (i = 0; i < number_monod; i++)
        cout << monod[i]->species << "  " << monod[i]->concentration << "  "
             << monod[i]->order << "\n";
    // Inhibition terms
    cout << "$INHIBITIONTERMS"
         << "\n"
         << number_inhibit << "\n";
    for (i = 0; i < number_inhibit; i++)
        cout << inhibit[i]->species << "  " << inhibit[i]->concentration << "  "
             << inhibit[i]->order << "\n";
    // Production Terms
    cout << "$PRODUCTIONTERMS"
         << "\n"
         << number_production << "\n";
    for (i = 0; i < number_production; i++)
        cout << production[i]->species << "  " << production[i]->concentration
             << "  " << production[i]->order << "\n";
    // ProductionStochhelp Terms
    cout << "$PRODUCTIONSTOCH"
         << "\n"
         << (int)ProdStochhelp.size() << "\n";
    for (i = 0; i < (int)ProdStochhelp.size(); i++)
        cout << ProdStochhelp[i]->species << "  "
             << ProdStochhelp[i]->concentration << "\n";
    // exchange
    cout << "$EXCHANGE_PARAMETERS"
         << "\n"
         << (int)ex_param.size() << "\n";
    for (i = 0; i < (int)ex_param.size(); i++)
        cout << ex_param[i] << "  ";
    cout << "\n";

    cout << "\n";

    cout << "number_reactionpartner " << (int)number_reactionpartner << "\n";
    cout << "bacteria_number " << (int)bacteria_number << "\n";
    cout << "grow " << grow << "\n";
    length = (int)ProductionStoch.size();
    cout << "length ProductionStoch: " << length << "\n";
    for (i = 0; i < length; i++)
        cout << (int)ProductionStoch[i] << " ";
    cout << "\n";

    length = (int)ex_species.size();
    cout << "length exSpecies: " << length << "\n";
    for (i = 0; i < length; i++)
        cout << ex_species_names[i] << " " << ex_species[i] << "\n";
    cout << "\n";
    cout << " sorption type : " << exType << "\n";

    // Test output
}

/**************************************************************************
Reaction-Method:
Task: returns true if NAPL dissolution is modeled
Programing:
08/2008 CB Implementation
**************************************************************************/
bool KNaplDissCheck(void)
{
    bool NAPLdiss = false;

    CKinReactData* m_krd = NULL;
    if (KinReactData_vector.size() == 0)
        return NAPLdiss;

    else
    {
        m_krd = KinReactData_vector[0];
        if (m_krd->NumberNAPLdissolution > 0)
            NAPLdiss = true;

        return NAPLdiss;
    }
}
/**************************************************************************
Reaction-Method:
Task: returns true if mineral kinetics is modeled

Programing:
12/2010 CB Implementation
**************************************************************************/
bool KMinKinCheck(void)
{
    bool MinKin = false;

    CKinReactData* m_krd = NULL;

    if (KinReactData_vector.size() == 0)
        return MinKin;

    else
    {
        m_krd = KinReactData_vector[0];
        if (m_krd->NumberMineralkinetics > 0)
            MinKin = true;
        return MinKin;
    }
}

/**************************************************************************/
// Task: Calculates and sets the NAPL density for all nodes
//      in case of NAPL dissolution
// Programing:
//   01/2008   CB   Implementation
/**************************************************************************/
void KNaplCalcDensity(void)
{
    CRFProcess* m_pcs = NULL;
    CFluidProperties* m_mfp = NULL;
    // CFEMesh* m_msh = NULL;

    double density = 0;
    // double viscosity = 0;
    long i;
    int ndx_density_phase = -1;
    long nnodes = 0;

    m_pcs = PCSGetFlow();
    ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY2");

    nnodes = int(m_pcs->m_msh->nod_vector.size());
    m_mfp = mfp_vector[m_pcs->pcs_type_number + 1];
    m_mfp->mode = 1;

    for (i = 0; i < nnodes; i++)
    {
        density = CalcNAPLDens(i);

        // CB: Density2 reintroduced as secondary variable of PS_GLOBAL,
        m_pcs->SetNodeValue(i, ndx_density_phase, density);
    }
    m_mfp->mode = 0;
}

/**************************************************************************
Task: Postprocessing function calculates the NAPL saturation after
      Flow, Transport and kinetic NAPL dissolution

Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void CalcNewNAPLSat()
{
    long i, j, k, l;
    int idx0 = 0;
    int idx1, idx2, idxC;
    long nnodes, nNAPLcomps;
    double conc, conc2, rho_N_new, rho_N_old, rho_N_fluid;
    double mass, volume;
    double satu_N_new, satu_N_old;

    vector<int> pcs_napl_comps_vector;
    vector<double> molar_weights_vector;
    vector<double> molar_densities_vector;

    string var_name;
    int no_processes = (int)pcs_vector.size();
    CRFProcess* m_pcs = NULL;
    CRFProcess* n_pcs = NULL;
    CFluidProperties* m_mfp = NULL;

    m_pcs = PCSGetFlow();

    bool psg = false;
    if (m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
        psg = true;

    nnodes = (long)fem_msh_vector[0]->nod_vector.size();

    // CB PS_GLOBAL this needs update
    m_mfp = mfp_vector[m_pcs->pcs_type_number];
    // m_mfp->mode = 1; // CB ??

    // Get indices of node value variables for phase 2
    if (psg)
        idx0 = m_pcs->GetNodeValueIndex("SATURATION2");  // old timelevel
    idx1 =
        m_pcs->GetNodeValueIndex("DENSITY2");  // CB PS_GLOBAL this needs update
    idx2 = m_pcs->GetNodeValueIndex("SATURATION1");  // old timelevel

    i = j = k = l = 0;
    no_processes = (int)pcs_vector.size();

    // get the parameters
    for (i = 0; i < no_processes; i++)
    {
        n_pcs = pcs_vector[i];
        // if(n_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0)
        if (n_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            j = n_pcs->GetProcessComponentNumber();
            if (cp_vec[j]->transport_phase == 3)
            {  // is in napl
                pcs_napl_comps_vector.push_back(
                    i);  // store processes that are NAPL dissolution components
                if (KinBlob_vector[0]->gas_dissolution_flag == false)
                {
                    molar_weights_vector.push_back(
                        cp_vec[j]->molar_weight);  // get the corresponding
                                                   // molar weight
                    molar_densities_vector.push_back(
                        cp_vec[j]->molar_density);  // get the corresponding
                                                    // densities
                }
                l++;
            }
        }
    }
    nNAPLcomps = l;

    // do the update
    for (i = 0; i < nnodes; i++)
    {
        conc = rho_N_new = rho_N_old = rho_N_fluid = satu_N_new = satu_N_old =
            mass = volume = conc2 = 0;

        if (KinBlob_vector[0]->gas_dissolution_flag == false)
        {
            // determine the old NAPL DENSITY after flow / transport step
            for (j = 0; j < nNAPLcomps; j++)
            {
                l = pcs_napl_comps_vector[j];
                idxC = pcs_vector[l]->GetNodeValueIndex(
                    pcs_vector[l]->pcs_primary_function_name[0]);
                conc = pcs_vector[l]->GetNodeValue(i, idxC);  // old timelevel
                mass += conc * molar_weights_vector[j];
            }
            if (mass > 0)
                rho_N_old = mass;  // [kg/mÂ³REV]
            else
                rho_N_old = 0;
            // determine the new NAPL density rho_N_neu at current node
            mass = 0;
            for (j = 0; j < nNAPLcomps; j++)
            {
                l = pcs_napl_comps_vector[j];
                idxC = pcs_vector[l]->GetNodeValueIndex(
                    pcs_vector[l]->pcs_primary_function_name[0]);
                conc = pcs_vector[l]->GetNodeValue(
                    i, idxC + 1);  // +1 new timelevel
                if (fabs(conc) < 1e-19)
                    conc = 0.0;
                mass += conc * molar_weights_vector[j];
                volume +=
                    conc /
                    molar_densities_vector[j];  // this is required calculating
                                                // the napl fluid density
            }
            if (mass > 0)
                rho_N_new = mass;  // [kg/mÂ³REV]
            else
                rho_N_new = 0;

            // get the old SATURATION2 of NAPL after flow / transport step
            if (psg)  // PS_GLOBAL
                satu_N_old = m_pcs->GetNodeValue(i, idx0 + 1);
            else  // MULTI_PHASE_FLOW
                satu_N_old = 1.0 - m_pcs->GetNodeValue(i, idx2 + 1);
            // calculate new NAPL Saturation: dSatu = Satu(t+dt)-Satu(t) =
            // Satu(t)*(1-rho(t)/rho(t+dt))
            if (satu_N_old * rho_N_new * rho_N_old > 0)
                satu_N_new = satu_N_old +
                             satu_N_old * (1 - rho_N_old / rho_N_new);  // fast
            // satu_N_new = satu_N_old + satu_N_old * (rho_N_new / rho_N_old -
            // 1);   //slow
            else
                satu_N_new = satu_N_old;
            if (satu_N_new < 0)
            {
                cout << " Warning in fct CalcNewNAPLSat: NAPL-Sat: "
                     << satu_N_new << "\n";
                satu_N_new = MRange(0.0, satu_N_new, 1.0);
            }

            // set new SATURATION2
            // if(satu_N_new > 0.00001)  satu_N_new = 1-0.95; // 0.985
            // CB 14.01.09 Hansen+Kueper BM

            // m_pcs->SetNodeValue(i, idx0, satu_N_new);  // idx0 for timelevel
            // 0 ??
            if (psg)
                m_pcs->SetNodeValue(i, idx0 + 1,
                                    satu_N_new);  // idx0+1 for timelevel 1 ??
            // set new SATURATION1
            // m_pcs->SetNodeValue(i, idx2, 1-satu_N_new);     // idx2 for
            // timelevel 0 ??
            m_pcs->SetNodeValue(i, idx2 + 1,
                                1 - satu_N_new);  // idx2+1 for timelevel 1 ??
            // finally determine the new napl fluid density
            if (mass * volume > 0)
                rho_N_fluid =
                    mass / volume;  // [kg/mÂ³N] = [kg/mÂ³REV] / [mÂ³N/mÂ³REV]
            else
                rho_N_fluid = m_mfp->Density();  // use NAPL phase fluid density
                                                 // as defined in .mfp
            // set new DENSITY2
            m_pcs->SetNodeValue(i, idx1, rho_N_fluid);
            // CB PS_GLOBAL this needs update
        }
        else
        {  // gas_dissolution
            // get the old SATURATION2 of NAPL after flow / transport step
            if (psg)  // PS_GLOBAL
                satu_N_old = m_pcs->GetNodeValue(i, idx0);
            else  // MULTI_PHASE_FLOW only has water satu
                satu_N_old =
                    1.0 - m_pcs->GetNodeValue(
                              i, idx2 + 1);  // idx2+1 as in case of Eclipse
                                             // coupling, Sats are transient.
            // hence the sat before KRC is available only at "new" TL
            // get the old and new CO2 conc in gas phase after flow / transport
            // step
            for (j = 0; j < nNAPLcomps; j++)
            {
                l = pcs_napl_comps_vector[j];
                idxC = pcs_vector[l]->GetNodeValueIndex(
                    pcs_vector[l]->pcs_primary_function_name[0]);
                conc += pcs_vector[l]->GetNodeValue(i, idxC);  // old timelevel
                conc2 +=
                    pcs_vector[l]->GetNodeValue(i, idxC + 1);  // new timelevel
            }
            satu_N_new = satu_N_old * (conc2 / conc);  // simple scaling
            if (psg)  // Set new NAPL SATURATION2
                m_pcs->SetNodeValue(i, idx0 + 1, satu_N_new);
            // Set new water SATURATION1 in any case
            m_pcs->SetNodeValue(i, idx2 + 1, 1 - satu_N_new);
        }
    }
    pcs_napl_comps_vector.clear();
    molar_weights_vector.clear();
    molar_densities_vector.clear();

    // m_pcs->WriteAllVariables();
}

/**************************************************************************
Task: Postprocessing function calculates the NAPL saturation after
      Flow, Transport and kinetic NAPL dissolution

Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void CalcNewPhasePressure()
{
    long i;  //, j=0;
    double Cgo, Cgn, Clo, Cln, TT, Po, Pn, H2O, rho;
    double Cgn2, Cln2;
    double Vg = 0, Vl = 1;   //, delta = 0;
    double NumberReactions;  //, scale=0;

    // int idx0=0;
    int idxCg, idxCl, idxP, idxS, sp;
    long nnodes;

    double poro = 0, Satu = 1;

    CRFProcess* m_pcs = NULL;
    m_pcs = PCSGetFlow();
    if (m_pcs->getProcessType() != FiniteElement::MULTI_PHASE_FLOW)
    {
        if (m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
            // this is for NAPL only, for CO2 use the routine below
            CalcNewNAPLSat();
        return;
    }

    // this is only for CO2 dissolution

    CRFProcess* g_pcs = NULL;
    CRFProcess* l_pcs = NULL;

    idxP = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;    // Pnw new TL
    idxS = m_pcs->GetNodeValueIndex("SATURATION1") + 1;  // Sw new TL

    NumberReactions = KinReact_vector.size();

    // Get dissolution Reaction
    CKinReact* m_kr = NULL;
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_napldissolution)
            break;
    }
    if (m_kr->typeflag_napldissolution != 1)
        return;

    // CKinReactData *m_krd = NULL;
    // m_krd = KinReactData_vector[0];

    nnodes = (long)fem_msh_vector[0]->nod_vector.size();

    // do the Pressure update
    for (i = 0; i < nnodes; i++)
    {
        Cgo = Cgn = Clo = Cln = TT = Po = Pn = H2O = rho = Cgn2 = Cln2 = 0;
        // Get Sw
        Satu = m_pcs->GetNodeValue(i, idxS);

        // Get P, T, H2O
        if (REACTINT_vec.size() > 0)
        {
            Po = REACTINT_vec[0]->GetPressure(i);
            TT = REACTINT_vec[0]->GetTemperature(i);
            poro = REACTINT_vec[0]->node_porosity[i];
            H2O = REACTINT_vec[0]->water_conc[i] * Satu * poro;
        }

        // Get concentrations
        // C in NAPL / gas phase
        sp = m_kr->ex_species[0];
        g_pcs = cp_vec[sp]->getProcess();  // true: new->old
        idxCg = g_pcs->GetNodeValueIndex(m_kr->ex_species_names[0]);
        Cgo = g_pcs->GetNodeValue(i, idxCg);
        Cgn = g_pcs->GetNodeValue(i, idxCg + 1);

        sp = m_kr->ex_species[1];
        l_pcs = cp_vec[sp]->getProcess();
        idxCl = l_pcs->GetNodeValueIndex(m_kr->ex_species_names[1]);
        Clo = l_pcs->GetNodeValue(i, idxCl);
        Cln = l_pcs->GetNodeValue(i, idxCl + 1);
        Clo *= poro * Satu;
        Cln *= poro * Satu;

        // delta = fabs(fabs(Clo - Cln) - fabs(Cgo - Cgn));
        // if( delta > 1e-6)
        //  cout << "!!! " << Cgo << " " << Cgn << " " << Clo << " " << Cln << "
        //  " << ((Cgo - Cgn)+(Clo - Cln)) << " "
        //  << poro << " " << Satu << "\n";
        //  //cout << "Warning: Wrong Mass balance at node " << i << " in
        //  CalcNewPhasePressure. dC = " << delta << "\n";

        // P & S update
        Pn = Po;

        VLE_CalcNewPressure(TT, Pn, Vg, Vl, Cgo, Cgn, Clo, Cln, H2O, H2O, rho);
        Pn *= 1.0e+5;  // bar --> Pa
        m_pcs->SetNodeValue(i, idxP, Pn);

        Vl *= 1e-6 / poro;
        Vg *= 1e-6 / poro;
        Satu = 1 - Vg / (Vg + Vl);
        m_pcs->SetNodeValue(i, idxS, Satu);

        // as pressure and vol have changed, density and thus gas conc. is
        // different now
        Cgn2 = rho * 1000 * poro * (1 - Satu) / (44.009 / 1000);
        g_pcs->SetNodeValue(i, idxCg + 1, Cgn2);
        // put the residual to water phase
        Cln2 = (Cln /*+ (Cgn-Cgn2)*/) / (poro * Satu);
        l_pcs->SetNodeValue(i, idxCl + 1, Cln2);
    }
}

/**************************************************************************/
// Task: Calculates the NAPL mass flux across Model boundary
//      in case of NAPL infiltration
// Programing:
//   01/2008   CB   Implementation
/**************************************************************************/
void CalcNAPLCompMasses()
{
    long i, j, k, l;
    int idxC, idx, idxS;
    long nnodes, nelenodes, nNAPLcomps, nDisscomps;
    double conc = 0;
    double sat = 0;
    double maxconc = 0;
    double mass, volume;
    ofstream _dump;
    double time;
    int no_processes;
    const double* coord;

    double xmax = 0;
    double ymax = 0;
    double ymin = 1e9;
    double area;

    long idxVx;  //, idxVy, idxVz;
    // double vel_nod[3];
    double vel_nod;

    vector<int> pcs_napl_comps_vector;
    vector<double> conc_napl_comps_vector;
    vector<int> pcs_comps_vector;
    vector<double> conc_comps_vector;
    vector<double> mf_vector;

    MeshLib::CNode* m_nod = NULL;
    MeshLib::CElem* m_ele = NULL;
    CFEMesh* m_msh = NULL;
    m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt

    string var_name;
    CRFProcess* m_pcs = NULL;
    CRFProcess* n_pcs = NULL;
    CTimeDiscretization* m_tim = NULL;
    m_tim = time_vector[0];
    time = m_tim->time_current;

    i = j = k = l = 0;
    no_processes = (int)pcs_vector.size();
    nnodes = (long)fem_msh_vector[0]->nod_vector.size();
    volume = fem_msh_vector[0]->ele_vector[0]->GetVolume();

    m_pcs = PCSGetFlow();
    // get the indices of velocity of flow process
    idxVx = m_pcs->GetNodeValueIndex("VELOCITY_X1");
    // idxVy = m_pcs->GetNodeValueIndex("VELOCITY_Y1");
    // idxVz = m_pcs->GetNodeValueIndex("VELOCITY_Z1");
    idxS = m_pcs->GetNodeValueIndex("SATURATION1") + 1;  // new timelevel

    bool mpf = false;
    if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        mpf = true;

    // prepare indices and data size structures
    for (i = 0; i < no_processes; i++)
    {
        n_pcs = pcs_vector[i];
        // if(n_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0)
        if (n_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            j = n_pcs->GetProcessComponentNumber();
            if (cp_vec[j]->transport_phase == 3)
            {  // is in napl
                pcs_napl_comps_vector.push_back(
                    i);  // store processes that are NAPL dissolution components
                conc_napl_comps_vector.push_back(
                    0);  // get the corresponding molar weight
                l++;
            }
            else if (cp_vec[j]->transport_phase == 0)
            {
                pcs_comps_vector.push_back(i);
                conc_comps_vector.push_back(0);
                mf_vector.push_back(0);
                k++;
            }
        }
    }
    nNAPLcomps = l;
    nDisscomps = k;

    // integrate NAPL masses
    for (j = 0; j < nNAPLcomps; j++)
    {
        conc = mass = 0;
        l = pcs_napl_comps_vector[j];
        idxC = pcs_vector[l]->GetNodeValueIndex(
                   pcs_vector[l]->pcs_primary_function_name[0]) +
               1;  // new timelevel
        for (i = 0; i < nnodes; i++)
        {
            conc = pcs_vector[l]->GetNodeValue(i, idxC);
            sat = m_pcs->GetNodeValue(i, idxS);  // should be napl sat
            if (mpf)  // in case of MULTI_PHASE_FLOW, it is water sat
                sat = 1 - sat;
            if (sat > 1e-3)
                mass += conc * volume;
        }
        conc_napl_comps_vector[j] = mass;
    }

    // get downgrad  model boundary x coordinate
    for (i = 0; i < nnodes; i++)
    {
        m_nod = m_msh->nod_vector[i];
        coord = m_nod->getData();
        if (coord[0] > xmax)
            xmax = coord[0];
    }

    // get ele flux area
    m_ele = m_msh->ele_vector[0];
    area = m_ele->GetFluxArea();
    nelenodes = m_ele->GetNodesNumber(false);
    for (j = 0; j < nelenodes; j++)
    {
        idx = m_ele->GetNodeIndex(j);
        m_nod = m_msh->nod_vector[idx];
        coord = m_nod->getData();
        if (coord[2] > ymax)
            ymax = coord[2];
        if (coord[3] < ymin)
            ymin = coord[2];
    }
    area *= (ymax - ymin);

    // maximum C
    for (j = 0; j < nDisscomps; j++)
    {
        conc = maxconc = 0;
        l = pcs_comps_vector[j];
        idxC = 1 + pcs_vector[l]->GetNodeValueIndex(
                       pcs_vector[l]->pcs_primary_function_name[0]);
        for (i = 0; i < nnodes; i++)
        {
            conc = pcs_vector[l]->GetNodeValue(i, idxC);  // new timelevel
            if (conc > maxconc)
                maxconc = conc;
        }
        conc_comps_vector[j] = maxconc;
    }

    // Mass flux across model boundary
    for (j = 0; j < nDisscomps; j++)
    {
        mass = conc = 0;
        for (i = 0; i < nnodes; i++)
        {
            m_nod = m_msh->nod_vector[i];
            coord = m_nod->getData();
            if (coord[0] == xmax)
            {
                // Get the velocity components
                vel_nod = m_pcs->GetNodeValue(i, idxVx);
                // vel_nod[1] = m_pcs->GetNodeValue(i, idxVy);
                // vel_nod[2] = m_pcs->GetNodeValue(i, idxVz);
                // Get the concentration
                l = pcs_comps_vector[j];
                idxC = 1 + pcs_vector[l]->GetNodeValueIndex(
                               pcs_vector[l]->pcs_primary_function_name[0]);
                conc = pcs_vector[l]->GetNodeValue(i, idxC);  // new timelevel
                // mass flux
                mass += conc * vel_nod;
            }
        }
        mf_vector[j] = mass * area;
    }

    // output to Tecfile
    ifstream in;
    string file;
    file = FileName + "_Masses_.tec";
    _dump.setf(ios::scientific, ios::floatfield);
    _dump.precision(12);
    if (m_tim->step_current == 1)
    {
        _dump.open(file.c_str(), ios::out);
        // header
        _dump << "TITLE = \"Masses and concentrations over time\" "
              << "\n";
        _dump << "VARIABLES  = \"TIME [a]\",\"TIME [s]\"";
        for (j = 0; j < nNAPLcomps; j++)
        {
            l = pcs_napl_comps_vector[j];
            _dump << "\"" << pcs_vector[l]->pcs_primary_function_name[0]
                  << "\",";
        }
        for (j = 0; j < nDisscomps; j++)
        {
            l = pcs_comps_vector[j];
            _dump << "\"" << pcs_vector[l]->pcs_primary_function_name[0]
                  << "\",";
        }
        for (j = 0; j < nDisscomps; j++)
        {
            l = pcs_comps_vector[j];
            _dump << "\"" << pcs_vector[l]->pcs_primary_function_name[0]
                  << "_mf\",";
        }
        _dump << "\n" << flush;
        _dump << "ZONE T=\"Source Zone\""
              << "\n";
    }
    else
        _dump.open(file.c_str(), ios::app);
    // data
    _dump << time / 86400 / 365 << " " << time << " ";
    for (j = 0; j < nNAPLcomps; j++)
        _dump << conc_napl_comps_vector[j] << " ";
    for (j = 0; j < nDisscomps; j++)
        _dump << conc_comps_vector[j] << " ";
    for (j = 0; j < nDisscomps; j++)
        _dump << mf_vector[j] << " ";
    _dump << "\n";
    _dump.close();

    // clear vectors
    pcs_napl_comps_vector.clear();
    conc_napl_comps_vector.clear();
    pcs_comps_vector.clear();
    conc_comps_vector.clear();
    mf_vector.clear();
}

/**************************************************************************/
// Task: Calculates the NAPL density for a node
//      in case of NAPL dissolution
// Programing:
//   01/2008   CB   Implementation
/**************************************************************************/
double CalcNAPLDens(int node)
{
    long i, j, k, l;
    int idxC;  // idx1
    int nNAPLcomps;
    double conc, rho_N_new, mass, volume;
    string var_name;

    vector<int> pcs_napl_comps_vector;
    vector<double> molar_weights_vector;
    vector<double> molar_densities_vector;

    int no_processes = (int)pcs_vector.size();
    CRFProcess* m_pcs = NULL;
    CRFProcess* n_pcs = NULL;
    CFluidProperties* m_mfp = NULL;

    m_pcs = PCSGetFlow();

    m_mfp =
        mfp_vector[m_pcs->pcs_type_number + 1];  // CB ToDo: check, if this also
                                                 // applies for MULTI_PHASE_FLOW
    // m_mfp = mfp_vector[m_pcs->pcs_type_number];
    // m_mfp->mode = 1; // CB ??

    i = j = k = l = 0;

    // collect parameters and concentrations for NAPL-components
    for (i = 0; i < no_processes; i++)
    {
        n_pcs = pcs_vector[i];
        // if(n_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0)
        if (n_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            j = n_pcs->GetProcessComponentNumber();
            if (cp_vec[j]->transport_phase == 3)
            {  // is in NAPL component
                pcs_napl_comps_vector.push_back(
                    i);  // store processes that are NAPL dissolution components
                molar_weights_vector.push_back(
                    cp_vec[j]
                        ->molar_weight);  // get the corresponding molar weight
                molar_densities_vector.push_back(
                    cp_vec[j]
                        ->molar_density);  // get the corresponding densities
                l++;
            }
        }
    }
    nNAPLcomps = l;

    // determine the NAPL density rho_N_neu at current node
    conc = rho_N_new = mass = volume = 0;
    for (j = 0; j < nNAPLcomps; j++)
    {
        l = pcs_napl_comps_vector[j];
        idxC = pcs_vector[l]->GetNodeValueIndex(
            pcs_vector[l]->pcs_primary_function_name[0]);
        conc = pcs_vector[l]->GetNodeValue(node, idxC);
        if (fabs(conc) < 1e-19)
            conc = 0.0;
        mass +=
            conc *
            molar_weights_vector[j];  // [kg/m³REV] = [molN/m³REV] * [kg/molN]
        volume +=
            conc / molar_densities_vector[j];  // [m³N/m³REV] = [molN/m³REV] /
                                               // [molN/m³REV]
    }
    if (mass * volume > 0)
        rho_N_new = mass / volume;  // [kg/m³N] = [kg/m³REV] / [m³N/m³REV]
    else
        rho_N_new =
            m_mfp
                ->Density();  // use NAPL phase fluid density as defined in .mfp
    // This may be not possible if new spatially distributed density model will
    // be implemented.. What to do in case of no NAPL? Set a default density?

    pcs_napl_comps_vector.clear();
    molar_weights_vector.clear();
    molar_densities_vector.clear();

    return rho_N_new;
}

/**************************************************************************
   Reaction-Method:
   Task: returns the Porevelocity of the mobile (water) phase in case of a NAPL
   dissolution model and TwoPhaseFlow; v at the node as a inverse distance
   weighted mean of the connecting elements velocities
   Programing:
   08/2008 CB Implementation
   10/2010 TF changed access to process type
   09/2011 TF changed access to coordinates of mesh node,
    - substituted access to mesh_element from pointer to direct access into the
vector
    - made the mesh node a const pointer
    - made the pointer to the mesh const, made the mesh itself const
    - substituted pow(x,2) by x*x
    - reduced scope of some loop variables
**************************************************************************/
double CKinReact::GetNodePoreVelocity(long node_number)
{
    CRFProcess* pcs(PCSGetFlow());
    // CFEMesh const* const msh(fem_msh_vector[0]); //SB: ToDo hart gesetzt

    // long group;
    // long el, elem;
    long idxVx, idxVy, idxVz, idxs1;
    double vel_nod[3];  //, coord[3], vel_ele[3];
    // double distance, sum_w, weight;
    // double* grav_c;
    double PoreVel(0), poro(0), satu = 1.0;  // default
    double theta = pcs->m_num->ls_theta;

    // Get node saturation of mobile (water) phase
    if (pcs->getProcessType() == FiniteElement::PS_GLOBAL)
    {
        idxs1 = pcs->GetNodeValueIndex("SATURATION1");  // Sat of water phase
        satu = pcs->GetNodeValue(node_number, idxs1);
    }
    else if (pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
    {
        idxs1 = pcs->GetNodeValueIndex("SATURATION1");  // Sat of water phase
        satu = pcs->GetNodeValue(node_number, idxs1);
    }
    else if (pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
    {
        if (pcs->pcs_type_number == 0)
            // this is the saturation equation
            pcs = pcs_vector[pcs->pcs_number + 1];
        // Sat of water phase
        idxs1 = pcs->GetNodeValueIndex("SATURATION1");
        satu = pcs->GetNodeValue(node_number, idxs1);
    }
    else if (pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
    {
        // Sat of water phase
        idxs1 = pcs->GetNodeValueIndex("SATURATION1");
        satu = pcs->GetNodeValue(node_number, idxs1);
    }

    // Get the porosity as average from neighbouring elements from
    poro = GetPhaseVolumeAtNode(node_number, theta, 0);

    // initialize data structures
    for (int i = 0; i < 3; i++)
        vel_nod[i] = 0;
    PoreVel = 0;

    // get the indices of velocity of flow process
    idxVx = pcs->GetNodeValueIndex("VELOCITY_X1");
    idxVy = pcs->GetNodeValueIndex("VELOCITY_Y1");
    idxVz = pcs->GetNodeValueIndex("VELOCITY_Z1");
    // Get the velocity components
    vel_nod[0] = pcs->GetNodeValue(node_number, idxVx);
    vel_nod[1] = pcs->GetNodeValue(node_number, idxVy);
    vel_nod[2] = pcs->GetNodeValue(node_number, idxVz);
    // absolute value of velocity vector
    for (int i = 0; i < 3; i++)
        PoreVel += pow(vel_nod[i], 2);
    PoreVel = sqrt(PoreVel);
    // divide by porosity and saturation to obtain transport velocity
    PoreVel /= (poro * satu);

    return PoreVel;
}

/**************************************************************************
 Reaction-Method:
 Task: returns the Darcyvelocity of the mobile (water) phase in case of a NAPL
 dissolution model and TwoPhaseFlow; v at the node as a inverse distance
 weighted mean of the connecting elements velocities
 Programing:
 08/2008 CB Implementation
 10/2010 TF changed access to process type
 05/2013 SP added calculation of Darcyvelocity for Sherwoodmodel: Powers et al.
 1992
 **************************************************************************/
double CKinReact::GetNodeDarcyVelocity(long node)
{
    CRFProcess* m_pcs = NULL;

    long i;
    long idxVx, idxVy, idxVz;
    double vel_nod[3];
    double DarcyVel;

    m_pcs = PCSGetFlow();

    // initialize data structures
    for (i = 0; i < 3; i++)
        vel_nod[i] = 0;
    DarcyVel = 0;

    // get the indices of velocity of flow process
    idxVx = m_pcs->GetNodeValueIndex("VELOCITY_X1");
    idxVy = m_pcs->GetNodeValueIndex("VELOCITY_Y1");
    idxVz = m_pcs->GetNodeValueIndex("VELOCITY_Z1");
    // Get the velocity components
    vel_nod[0] = m_pcs->GetNodeValue(node, idxVx);
    vel_nod[1] = m_pcs->GetNodeValue(node, idxVy);
    vel_nod[2] = m_pcs->GetNodeValue(node, idxVz);
    // absolute value of velocity vector
    for (i = 0; i < 3; i++)
        DarcyVel += pow(vel_nod[i], 2);
    DarcyVel = sqrt(DarcyVel);

    return DarcyVel;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactionDeactivation()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Deactivates kinetic reaction calculation at individual nodes based on  */
/* evaluation of reaction rates of the previous time step in local        */
/* neighborhoods around nodes, or based on comparison to previous         */
/*   concentrations at the nod                                            */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 08/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactionDeactivation(long nonodes)
{
    long node;
    int react_t =
        10;  // reactions are calculated at all nodes every react_t timesteps
    size_t sp;
    const size_t Number_of_Components(cp_vec.size());
    if (ReactDeactRelative)
        cout << " Reaction deactivation mode " << ReactDeactMode
             << " with relative concentration change! "
             << "\n";
    else
        cout << " Reaction deactivation  mode " << ReactDeactMode
             << " with absolute concentration change! "
             << "\n";

    // reactivate all nodes every n time steps
    if (((aktueller_zeitschritt) % react_t) ==
        0)  // when evaluated before loop over nodes, i.e. prepare for this time
            // step

        // if(((aktueller_zeitschritt+1) % react_t ) == 0) {  // when evaluated
        // after  loop over nodes, i.e. prepare for next time step
        for (node = 0; node < nonodes; node++)
        {
            ReactDeact[node] = false;
            if (ReactDeactMode != 3)  // for all timesteps, anyway prepare the
                                      // concentrationmatrix for the next time
                // step, i.e. save the current concentrations after Transport
                for (sp = 0; sp < Number_of_Components; sp++)
                    concentrationmatrix[node][sp] =
                        cp_vec[sp]->getProcess()->GetNodeValue(node,
                                                               sp_varind[sp]);
        }
    // CB Now check if node can be deactivated for the next time step
    else
    {
        // only for mode 1: first calculate for each node the sum of rates and
        // store in vector, which then is accessed in next loop over nodes, to
        // evaluate the sum of sum of rates over the neighbours
        if (ReactDeactMode == 1)
            for (node = 0; node < nonodes; node++)
            {
                React_dCdT[node] = 0;
                for (sp = 0; sp < Number_of_Components; sp++)
                {
                    // This is C of last time step after reactions, i.e. the old
                    // time level for this time step
                    double Concentration =
                        cp_vec[sp]->getProcess()->GetNodeValue(
                            node, (sp_varind[sp] - 1));
                    // This is C of previous time step after transport only,
                    // which was stored in matrix
                    double Concentration_old = concentrationmatrix[node][sp];
                    // relative concentration change
                    double maxi = 1;
                    if (ReactDeactRelative)
                        maxi = DMAX(Concentration_old, Concentration);
                    if (maxi < ReactDeactCThresh)
                        React_dCdT[node] += 0.0;
                    else
                    {
                        React_dCdT[node] +=
                            fabs((Concentration_old - Concentration) / maxi) /
                            dt;  // normalized by current local concentration
                        // and now prepare concentrationmatrix for next time
                        // step, i.e. save current concentrations after
                        // Transport
                        concentrationmatrix[node][sp] =
                            cp_vec[sp]->getProcess()->GetNodeValue(
                                node, sp_varind[sp]);
                    }
                }
                // cout << nod << " " << React_dCdT[nod] << "\n";
            }

        // this is the check, if a node may be deactivated for this time step
        for (node = 0; node < nonodes; node++)
        {
            double sumReact_dCdT = 0;

            switch (ReactDeactMode)
            {
                case 1:  // loop over no of connected nodes and their respective
                         // neighbours
                    for (size_t nn = 0; nn < ReactNeighborhood[node].size();
                         nn++)
                    {
                        int node_idx = ReactNeighborhood[node][nn];
                        sumReact_dCdT += React_dCdT[node_idx];
                    }
                    if (ReactDeactRelative)
                        sumReact_dCdT /= double(ReactNeighborhood[node].size());
                    break;
                case 2:  // compare with C after transport of last timestep,
                         // loop over all components
                    for (sp = 0; sp < Number_of_Components; sp++)
                    {
                        // This is C of current time step after transport
                        double Concentration =
                            cp_vec[sp]->getProcess()->GetNodeValue(
                                node, sp_varind[sp]);
                        // This is C of previous time step after transport only
                        double Concentration_old =
                            concentrationmatrix[node][sp];
                        double maxi = 1;
                        if (ReactDeactRelative)
                            maxi = DMAX(Concentration_old, Concentration);
                        if (maxi > ReactDeactCThresh)
                            sumReact_dCdT += fabs(
                                (Concentration - Concentration_old) / maxi);
                        // and now prepare the concentrationmatrix for the next
                        // time step, i.e. save the current concentrations after
                        // Transport
                        concentrationmatrix[node][sp] = Concentration;
                    }
                    break;
                case 3:  // compare with C after reaction of last timestep, loop
                         // over all components
                    for (sp = 0; sp < Number_of_Components; sp++)
                    {
                        // This is C of current time step after transport
                        double Concentration =
                            cp_vec[sp]->getProcess()->GetNodeValue(
                                node, sp_varind[sp]);
                        // This is C of previous time step after transport &
                        // reaction
                        double Concentration_old =
                            cp_vec[sp]->getProcess()->GetNodeValue(
                                node, (sp_varind[sp] - 1));
                        double maxi = 1;
                        if (ReactDeactRelative)
                            maxi = DMAX(Concentration_old, Concentration);
                        if (maxi > ReactDeactCThresh)
                            sumReact_dCdT += fabs(
                                (Concentration - Concentration_old) / maxi);
                    }
                    break;
                default:
                    break;
            }

            // check if deactivation criterion is met
            if (ReactDeactRelative)
                sumReact_dCdT /= double(Number_of_Components);
            if (sumReact_dCdT < ReactDeactEpsilon)
                ReactDeact[node] = true;  // negligible change, deactivate the
                                          // node for the next time step
            else
                ReactDeact[node] = false;  // sufficient change, reactivate the
                                           // node for the next time step
        }
    }  // else

    // Resets the reaction rates vector, only in case of model 1
    if (ReactDeactMode == 1)
        for (node = 0; node < nonodes; node++)
            React_dCdT[node] = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactReset_dCdT()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Sets the C after reaction of last time step as new C after reaction    */
/* for deactivated nodes */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactDeactSetOldReactionTerms(long nonodes)
{
    long node;
    int sp;
    int Number_of_Components = (int)cp_vec.size();
    double Concentration;

    for (node = 0; node < nonodes; node++)
    {
        if (ReactDeact[node] == true)
            for (sp = 0; sp < Number_of_Components; sp++)
            {
                // Get the C after reactions of last time step (old time level,
                // index = 0)
                Concentration = cp_vec[sp]->getProcess()->GetNodeValue(
                    node, sp_varind[sp] - 1);
                // Set this C as the new concentration after reaction
                cp_vec[sp]->getProcess()->SetNodeValue(node, sp_varind[sp],
                                                       Concentration);
            }
    }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactPlotFlagsToTec()                        */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Prints flags for reaction deactivation in tecplot format               */
/* for all nodes                                                          */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void CKinReactData::ReactDeactPlotFlagsToTec()
{
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    std::string eleType;

    if (m_msh->getNumberOfLines() > 0)
        eleType = "QUADRILATERAL";
    if (m_msh->getNumberOfQuads() > 0)
        eleType = "QUADRILATERAL";
    if (m_msh->getNumberOfHexs() > 0)
        eleType = "BRICK";
    if (m_msh->getNumberOfTris() > 0)
        eleType = "QUADRILATERAL";
    if (m_msh->getNumberOfTets() > 0)
        eleType = "TETRAHEDRON";
    if (m_msh->getNumberOfPrisms() > 0)
        eleType = "BRICK";

    if (NumberReactions > 0)
    {
        std::string filename(FileName + "_Deactivated_nodes.tec");
        std::ofstream aus;
        if (aktueller_zeitschritt == 1)
            aus.open(filename.c_str());
        else
            aus.open(filename.c_str(), ios::app);

        const size_t nnodes(m_msh->nod_vector.size());
        const size_t nele(m_msh->ele_vector.size());

        aus << "VARIABLES = "
            << "\"x\""
            << " "
            << "\"y\""
            << " "
            << "\"z\""
            << "\"active\""
            << "\n";
        aus << "ZONE T="
            << "\"aktueller_zeitschritt=" << aktueller_zeitschritt << "\"";
        aus << ", N=" << nnodes << ", E=" << nele
            << " F=FEPOINT, ET=" << eleType << "\n";

        for (size_t i = 0; i < nnodes; i++)
        {
            double const* coord =
                (m_msh->nod_vector[i])->getData();  // Coordinates(coord);
            aus << coord[0] << " " << coord[1] << " " << coord[2] << " ";
            if (is_a_CCBC[i] == true)
                aus << 0 << "\n";
            else if (ReactDeact[i] == true)
                aus << 0 << "\n";
            else
                aus << 1 << "\n";
        }
        for (size_t i = 0; i < nele; i++)
            m_msh->ele_vector[i]->WriteIndex_TEC(aus);

        aus.close();
    }
}

/**************************************************************************
 Reaction-Method:
 Task: This function calculates decay of Aromaticum bacteria species
       independent from other kinetic reactions,
       reqiured for some model tests only
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void CKinReactData::Aromaticum(long nonodes)
{
    long node;
    double conc, lambda = 0.0;  // OK411
    // int pcsindex = 0;
    CRFProcess* m_pcs = NULL;
    int varindex = 0;
    int nospec = (int)sp_varind.size();

    for (int sp = 0; sp < nospec; sp++)
        if (cp_vec[sp]->getProcess()->nod_val_name_vector[0].compare(
                "Aromaticum") == 0)
        {
            // pcsindex = sp_pcsind[sp];
            m_pcs = cp_vec[sp]->getProcess();
            varindex = sp_varind[sp];
            break;
        }
    for (int sp = 0; sp < nospec; sp++)
        if (cp_vec[sp]->compname.compare("Aromaticum") == 0)
        {
            lambda = cp_vec[sp]->decay_model_values[0];
            break;
        }

    for (node = 0; node < nonodes; node++)
    {
        conc = m_pcs->GetNodeValue(node, varindex);
        conc = conc * exp(-lambda * dt);
        m_pcs->SetNodeValue(node, varindex, conc);
    }
}

/**************************************************************************
 Reaction-Method:
 Task: This function prepares
   - node velocities
   - vector of nodes for which NAPL dissolution shall be switched
     off once saturation drops below hard coded thereshhold
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void CKinReactData::NAPLDissolutionPreprocessing()
{
    CKinReact* m_kr = NULL;
    CRFProcess* m_pcs = NULL;
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt

    int i, j, sp;  // idx0,
    bool mpf = false;
    bool psg = false;

    // To access nodal pore velocities , these may be calculated for all nodes
    // beforehands
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_napldissolution)
        {
            m_pcs = PCSGetFlow();
            if (m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
                psg = true;
            else if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
                mpf = true;
            // if(m_pcs->pcs_type_name.compare("PS_GLOBAL")!=0)
            if (!psg && !mpf)
                cout << "WARNING in CKinReact::NAPLDissolutionPreprocessing: "
                        "Flow process cannot be used"
                     << "\n";
            // CB copied from Problem::PostCouplingLoop()
            else if (m_pcs->cal_integration_point_value)  // WW
                m_pcs->Extropolation_GaussValue();
            break;  // do only once
        }
    }

    // In case of Eclipse coupling and MULTI_PHASE_FLOW set data to old time
    // level
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_napldissolution)
        {
            // C in NAPL / gas phase
            sp = m_kr->ex_species[0];
            cp_vec[sp]->getProcess()->CopyTimestepNODValues(
                true);  // true: new->old
            // saturation
            m_pcs->CopyTimestepNODValues(true);
        }
    }

    // Switch off nodes with too low NAPL saturation
    // now set up vectors switched_off_node for individual reactions
    // if (psg)
    //  idx0 = m_pcs->GetNodeValueIndex("SATURATION2"); // old timelevel
    // check all reactions
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];

        if (m_kr->typeflag_napldissolution)
        {
            // Initialize vector if not yet done in KRConfig()
            if (m_kr->NotThisReactGeoType.size() == 0)
                for (j = 0; j < (long)m_msh->nod_vector.size(); j++)
                    m_kr->switched_off_node.push_back(false);
            // Go through all nodes
            // for(j=0;j<(long)m_msh->nod_vector.size();j++){
            // if(m_pcs->GetNodeValue(j, idx0) < 1e-13)
            //   m_kr->switched_off_node[j] = true;
            //}
        }  // reaction is NAPL diss

    }  // loop over k nreactions
}

/**************************************************************************
 Reaction-Method:
 Task: This function prepares
   - node velocities
   - vector of nodes for which NAPL dissolution shall be switched
     off once saturation drops below hard coded thereshhold
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void CKinReactData::PreprocessMicrobe_drmc_(double steplength)
{
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    MicrobeData* m_md = NULL;

    double theta, G, S;
    int i, j;

    for (i = 0; i < NumberMicrobeData; i++)
    {
        m_md = MicrobeData_vector[i];
        for (j = 0; j < (long)m_msh->nod_vector.size(); j++)
        {
            // here, get the gibbs energy for the reaction at beginning of time
            // step
            G = m_md->GetGibbsEnergy(j);
            m_md->Gibbs[j] = G;
            theta =
                1 / (exp((m_md->G0 - G) / (m_md->steepness * m_md->G0)) + 1);
            if (theta > 1.0)
                theta = 1.0;
            // variable dt is updated only after the first time step and is
            // initialized with zero before the first time step, update of S in
            // the first time step hence is Snew = Sold * exp(0) this is an
            // analytical approximtion of dS/dt =
            // [kincr*theta*(1-theta_s)-kdecr*(1-theta)] * S
            S = m_md->_drmc_level[j] * exp((m_md->k_incr * theta * (1 - theta) -
                                            m_md->k_decr * (1 - theta)) *
                                           m_md->dt);
            if (S > 1.0)
                S = 1.0;
            m_md->_drmc_level[j] = S;
        }
        // store time step length for next time step update
        m_md->dt = steplength;
    }
}

/**************************************************************************
 Reaction-Method:
 Task: This function prepares
   - activity cooeffcient matrix (coupling to Chemapp reqiuired for model 3)
   - Equilibrium constant matrix (coupling to Chemapp reqiuired in case of
     non-uniform models CAP or HKF)
   - reactive surface areas
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void CKinReactData::PreprocessMinKin()
{
    int i, j, k = 0, nnode = 0, ncomp, widx = 0;  // node,
    double val, conc, A0, poro;
    // double fact;
    double A = 0, B = 0, a0, I, Imin, Imax, dens;
    double T = 298.15;
    double TC;
    double phi, MV, epsi;
    bool warning = false;
    double unitfactor_l = 1;

    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    CKinReact* m_kr = NULL;
    nnode = (long)m_msh->nod_vector.size();
    ncomp = (int)cp_vec.size();
// CB merge CAP 0311
#ifdef OGS_FEM_CAP  // CAP_REACT
    REACT_CAP* m_cap;
#endif

    REACTINT* m_rei = NULL;
    if (REACTINT_vec.size() > 0)
        m_rei = REACTINT_vec[0];

    // CB merge CAP 0311
    // get ChemApp data structure
    if (activity_model == 3)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        if (REACT_CAP_vec.size() > 0)
            m_cap = REACT_CAP_vec[0];
#else
        cout << " Warning: in PreprocessMinKin: Activity model = 3, but CAP is "
                "not defined. "
             << "\n";
#endif
    }
    else
    {
        for (i = 0; i < NumberReactions; i++)
        {
            m_kr = KinReact_vector[i];
            if (m_kr->typeflag_mineralkinetics && m_kr->Km_uniform == false)
            {
#ifdef OGS_FEM_CAP  // CAP_REACT
                m_cap = REACT_CAP_vec[0];
                break;
#else
                cout << " Warning: in PreprocessMinKin: Km is not uniform, but "
                        "CAP is not defined. "
                     << "\n";
#endif
            }
        }
    }

    // get water index
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_mineralkinetics)
        {
            widx = m_kr->water_number;
            break;
        }
    }

    if (activity_model == 1 || activity_model == 2)
    {
        // first calculate ionic strength at each node
        // I = 0.5SUM(z_i^2 * Moles_i / MassWater) = 0.5SUM(z_i^2 * C_i)
        Imin = 1000;
        Imax = 0;
        for (j = 0; j < nnode; j++)
        {
            IonicStrengths[j] = 0;
            // get the unit conversion factor: mol/m³l --> molality
            if (m_rei)
            {
                if (m_rei->unitconversion)
                    unitfactor_l = MOLH2OPERKG / m_rei->water_conc[j];
            }

            for (i = 0; i < ncomp; i++)
            {
                val = cp_vec[i]->valence;  //  Valenz;
                // Attention: next conc should be molality [mol/kg_h2o], so
                // miltiply with conversion factor
                conc = cp_vec[i]->getProcess()->GetNodeValue(j, sp_varind[i]) *
                       unitfactor_l;
                // conc = pcs_vector[sp_pcsind[i]]->GetNodeValue(j,sp_varind[i])
                // * unitfactor_l; // CB HS update
                if (cp_vec[i]->compname.compare("pH") == 0)
                    conc = pow(
                        10,
                        -conc *
                            unitfactor_l);  // this is inacurate as it should be
                                            // molality of H+, not the activity.
                IonicStrengths[j] += 0.5 * val * val * conc;
            }
            // store min and max I in model domain
            if (Imin > IonicStrengths[j])
                Imin = IonicStrengths[j];
            if (Imax < IonicStrengths[j])
                Imax = IonicStrengths[j];
        }  // nnodes
        // check activity model vs. min / max ionic strength
        switch (activity_model)
        {
            case 0: /* constant activity = 1*/
                break;
            case 1: /* DH*/
                if (Imax > 0.1)
                    Display::DisplayMsgLn(
                        "Warning in CKinReactData::PreprocessMinkin(): I > 0.1 "
                        "--> use of DH activity model is not "
                        "recomended!");
                break;
            case 2: /* Davies*/
                if (Imin < 0.1)
                    Display::DisplayMsgLn(
                        "Warning in CKinReactData::PreprocessMinkin(): I < 0.1 "
                        "--> use of Davies activity model is not "
                        "recomended!");
                break;
            case 3: /* CHEMAPP*/
                break;
            default:
                Display::DisplayMsgLn(
                    "Unknown activity model in "
                    "CKinReactData::PreprocessMinkin()!");
                break;
        }

        // then calculate T dependent parameters A & B, from Groundwater
        // Geochemistry, Broder Merkel Todo: get T
        TC = T - PhysicalConstant::CelsiusZeroInKelvin;
        dens = 1 -
               (TC - 3.9863) * (TC - 3.9863) * (TC + 288.9414) / 508929.2 /
                   (TC + 68.12963) +
               0.011445 * exp(-374.3 / TC);
        epsi = 2727.586 + 0.6224107 * T - 466.9151 * log(T) - 52000.87 / T;
        A = 1.82483e6 * pow(dens, 0.5) * pow(epsi * T, -1.5);
        B = 50.2916 * pow(dens, 0.5) * pow(epsi * T, -0.5);
    }

    // now calculate activity coefficients for all species at all nodes
    switch (activity_model)
    {
        case 0: /* constant activity = 1*/
            // do nothing
            break;
        case 1: /* DH*/
            for (j = 0; j < ncomp; j++)
            {
                if (cp_vec[j]->transport_phase == 0)
                {
                    val = double(cp_vec[j]->valence);  //  Valenz;
                    a0 = cp_vec[j]->a_zero;            //  a0;
                    for (i = 0; i < nnode; i++)
                    {
                        I = IonicStrengths[i];
                        if (I > 0)
                        {
                            if (fabs(val) > 0)  // charged species
                                ActivityCoefficients[i][j] =
                                    pow(10, -A * val * val * sqrt(I) /
                                                (1 + B * a0 * sqrt(I)));
                            else
                            {  // uncharged species
                                if (j == widx)
                                    ActivityCoefficients[i][j] =
                                        1;  // set water activity to 1
                                else
                                    ActivityCoefficients[i][j] = pow(
                                        10,
                                        0.1 *
                                            I);  // as in PHREEQC, manual, p 157
                            }
                        }
                        else
                        {
                            ActivityCoefficients[i][j] = 1;
                            if (i == 0 && j == 0)
                                cout << "Warning in KinReactData: zero Ionic "
                                        "strength in PreprocessMinKin; using "
                                        "gamma "
                                        "= 1.0"
                                     << "\n";
                        }
                    }
                }
            }
            break;
        case 2: /* DAVIS*/
            for (j = 0; j < ncomp; j++)
            {
                if (cp_vec[j]->transport_phase == 0)
                {
                    val = double(cp_vec[j]->valence);  //  Valenz;
                    for (i = 0; i < nnode; i++)
                    {
                        I = IonicStrengths[i];
                        if (I > 0)
                        {
                            if (fabs(val) > 0)  // charged species
                                ActivityCoefficients[i][j] = pow(
                                    10,
                                    -A * val * val *
                                        (sqrt(I) / (1 + sqrt(I)) - 0.3 * I));
                            else
                            {  // uncharged species
                                if (j == widx)
                                    ActivityCoefficients[i][j] =
                                        1;  // set water activity to 1
                                else
                                    ActivityCoefficients[i][j] = pow(
                                        10,
                                        0.1 *
                                            I);  // as in PHREEQC, manual, p 157
                            }
                        }
                        else
                        {
                            ActivityCoefficients[i][j] = 1;
                            if (i == 0 && j == 0)
                                cout << "Warning in KinReactData: zero Ionic "
                                        "strength in PreprocessMinKin; using "
                                        "gamma "
                                        "= 1.0"
                                     << "\n";
                        }
                    }
                }
            }
            break;
        case 3:     /* CHEMAPP*/
#ifdef OGS_FEM_CAP  // CAP_REACT
            for (j = 0; j < ncomp; j++)
            {  // CB merge CAP 0311
                // first, find the matching CAP species idx k
                if (REACT_CAP_vec.size() > 0)
                {
                    for (k = 0; k < m_cap->mass_num; k++)
                    {  // Attention, m_cap->mass_num can be less than ncomp
                        if (cp_vec[j]->compname.compare(
                                m_cap->species_name[k]) == 0)
                            break;
                    }  // CB merge CAP 0311
                }      // CB merge CAP 0311
                for (i = 0; i < nnode; i++)
                {
                    // Chemapp now returns activity coefficients gamma instead
                    // of activities.

                    // then, get the concentration of species j at node i
                    // Attention: next conc should be molality [mol/kg_h2o]
                    // conc=pcs_vector[sp_pcsind[j]]->GetNodeValue(i,sp_varind[j]);
                    conc =
                        cp_vec[j]->getProcess()->GetNodeValue(i, sp_varind[j]);
                    if (cp_vec[j]->compname.compare("pH") == 0)
                        conc = pow(10,
                                   -conc);  // this is inacurate as it should be
                                            // molality of H+, not the activity.
                    // finally, get activity from vector node_ac and divide by
                    // conc to get gamma
                    if (REACT_CAP_vec.size() > 0)
                        if (m_cap->node_ac[i].size() > 0)
                        {  // CB merge CAP 0311
                            if (/*conc>0 && */ k < m_cap->mass_num)
                                ActivityCoefficients[i][j] =
                                    m_cap->node_ac[i][k] /*/ conc*/;
                            else
                                ActivityCoefficients[i][j] = 1.0;
                        }
                        else
                        {  // CB merge CAP 0311
                            ActivityCoefficients[i][j] = 1.0;
                            if (j == 0 && is_a_CCBC[i] == false)
                            {
                                cout << " Warning in "
                                        "CKinReactData::PreprocessMinKin():"
                                     << "\n";
                                cout << "    No Chemapp-activity coefficient "
                                        "data available for node "
                                     << i << "\n";
                            }
                        }  // CB merge CAP 0311
                }
            }
#else
            cout << " Warning: in PreprocessMinKin: Km is not uniform, but CAP "
                    "is not defined. "
                 << "\n";
#endif
            break;
        default:
            Display::DisplayMsgLn(
                "Unknown activity model in CKinReactData::PreprocessMinkin()!");
            break;
    } /* switch */

    // Set K_m as function of P and T in vector m_kr->Km[node]; get data from
    // ChemApp
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_mineralkinetics && m_kr->Km_uniform == false)
        {
#ifdef OGS_FEM_CAP  // CAP_REACT   //CB merge CAP 0311
            if (m_kr->Km_CHEMAPP == true)
            {  // find the matching CAP Kinetic reaction
                for (j = 0; j < int(m_cap->Kin_Reactions.size()); j++)
                {
                    if (m_kr->chemapp_name.compare(
                            m_cap->Kin_Reactions[j].type) == 0)
                    {  // found
                        for (k = 0; k < nnode; k++)
                        {
                            if (m_cap->node_logK[k].size() > 0)
                                m_kr->Km[k] = pow(10, -m_cap->node_logK[k][j]);
                            else if (is_a_CCBC[k] == false)
                            {  // what, if no data has been calculated by
                               // chemapp??
                                cout << " Warning in "
                                        "CKinReactData::PreprocessMinkin:"
                                     << "\n";
                                cout << "    No Chemapp-log(K) data available "
                                        "for reaction "
                                     << m_kr->chemapp_name;
                                cout << " at node " << k << "\n";
                                if (k > 0)
                                    m_kr->Km[k] =
                                        m_kr->Km[k - 1];  // use neighbour node
                                                          // value
                                else
                                    m_kr->Km[k] =
                                        m_kr->Km_default;  // use default input
                                                           // value
                            }
                        }
                        break;
                    }
                }
            }
            else if (m_kr->Km_HKF == true)
            {  // find the matching HKF Kinetic reaction
                for (j = 0; j < int(m_cap->Kin_HKF_Reactions.size()); j++)
                {
                    if (m_kr->chemapp_name.compare(
                            m_cap->Kin_HKF_Reactions[j].type) == 0)
                    {  // found
                        for (k = 0; k < nnode; k++)
                        {
                            if (m_cap->node_HKF_logK[k].size() > 0)
                                m_kr->Km[k] =
                                    pow(10, m_cap->node_HKF_logK[k][j]);
                            else if (is_a_CCBC[k] == false)
                            {  // what, if no data has been calculated by
                               // chemapp??
                                cout << " Warning in "
                                        "CKinReactData::PreprocessMinkin:"
                                     << "\n";
                                cout << "    No HKF-log(K) data available for "
                                        "reaction "
                                     << m_kr->chemapp_name;
                                cout << " at node " << k << "\n";
                                if (k > 0)
                                    m_kr->Km[k] =
                                        m_kr->Km[k - 1];  // use neighbour node
                                                          // value
                                else
                                    m_kr->Km[k] =
                                        m_kr->Km_default;  // use default input
                                                           // value
                            }
                        }
                        break;
                    }
                }
            }
#else
            cout << " Warning: in PreprocessMinKin: Km is not uniform, but CAP "
                    "is not defined. "
                 << "\n";
#endif
        }
    }  // CB merge CAP 0311

    // set mineral reactive surface areas
    // they are updated each time step by the remaing mass fraction C(t)/C0
    // before first zime step, prepare Am/C0, this is done just once at the
    // beginning of simulation
    if (aktueller_zeitschritt == 1)
    {
        for (i = 0; i < NumberReactions; i++)
        {
            m_kr = KinReact_vector[i];

            if (m_kr->typeflag_mineralkinetics && m_kr->Am_constant == false)
            {
                for (k = 0; k < nnode; k++)
                {
                    // conc =
                    // pcs_vector[sp_pcsind[m_kr->mineral_number]]->GetNodeValue(k,sp_varind[m_kr->mineral_number]);
                    // // CB HS update
                    conc =
                        cp_vec[m_kr->mineral_number]
                            ->getProcess()
                            ->GetNodeValue(k, sp_varind[m_kr->mineral_number]);
                    switch (m_kr->Am_model)
                    {
                        case 1:
                            // A0/C0  to multiply with C(t) before each time
                            // step
                            A0 = m_kr->Am[0];  // save area
                            if (conc > 0)
                                m_kr->Am[k] = A0 / conc;
                            break;
                        case 2:
                            // store thi initial mineral concentration for
                            // surface updating in postprocessing
                            m_kr->Cminini.push_back(conc);
                            break;
                        case 3:
                            // do nothing
                            break;
                        case 4:
                            // do nothing
                            break;
                        case 5:
                            // do nothing
                            break;
                        default:
                            // do nothing
                            break;
                    }
                }
            }
        }
    }

    // here, update each time step by the remaing mass fraction C(t)/C0
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_mineralkinetics && m_kr->Am_constant == false)
        {
            MV = cp_vec[m_kr->mineral_number]->molar_weight /
                 cp_vec[m_kr->mineral_number]->mineral_density;

            for (k = 0; k < nnode; k++)
            {
                // A0/C0 here is multiplied with C(t) before current time step
                // conc =
                // pcs_vector[sp_pcsind[m_kr->mineral_number]]->GetNodeValue(k,sp_varind[m_kr->mineral_number]);
                // // CB HS update
                conc = cp_vec[m_kr->mineral_number]->getProcess()->GetNodeValue(
                    k, sp_varind[m_kr->mineral_number]);

                phi = 0;

                switch (m_kr->Am_model)
                {
                    case 1:
                        if (conc > 0)
                            m_kr->Am[k] *= conc;
                        if (m_kr->Am[k] < 1e-6)
                        {
                            m_kr->Am[k] =
                                1e-6;  // minimum allowed surface area to allow
                                       // regrowth of a mineral
                            warning = true;
                            // node = k;
                        }
                        break;
                    case 2:
                        // do nothing here, update is done in preprocessing
                        // function
                        break;
                    case 3:
                        // do nothing here, update is done in preprocessing
                        // function
                        break;
                    case 4:
                        // do nothing here, update is done in preprocessing
                        // function
                        break;
                    case 5:
                        // SA is proportional to mineral volume fraction
                        if (conc > 0)
                        {
                            if (REACTINT_vec.size() > 0)
                                poro = REACTINT_vec[0]->node_porosity[k];
                            else
                                poro = m_kr->GetPhaseVolumeAtNode(
                                    k, 1, 0);  // no poro update when RACINT
                                               // does not exist
                            phi = conc * MV * (1 - poro);
                        }
                        phi = DMAX(1e-7,
                                   phi);  // guarantee a minimum surface area
                        // else  phi = 1e-6;
                        m_kr->Am[k] =
                            phi * m_kr->Am_ini *
                            cp_vec[m_kr->mineral_number]->molar_weight / MV;
                        break;
                    default:
                        // do nothing here, update is done in preprocessing
                        // function
                        break;
                }
            }
        }
    }
    if (warning)
        cout << " Warning: A mineral concentration is 0.0 at node " << k << "\n"
             << " --> setting surface area to minimum of 1e-8."
             << "\n";

    CTimeDiscretization* m_tim = NULL;
    m_tim = time_vector[0];
    int steps = int(m_tim->time_step_vector.size());
    bool plot = false;
    if (aktueller_zeitschritt == 1 || aktueller_zeitschritt % 10 == 0 ||
        steps == int(aktueller_zeitschritt))
        plot = true;

    // debug output
    if (debugoutflag && activity_model > 0 && plot == true)
    {
        c_dumpfilename = FileName + "_Activities.dump";
        c_dump.setf(ios::scientific, ios::floatfield);
        c_dump.precision(12);
        c_dump.open(c_dumpfilename.c_str());
        // header
        c_dump << "Node IonicStrength ";
        for (i = 0; i < ncomp; i++)
            c_dump << cp_vec[i]->compname << " ";
        c_dump << "\n";
        // data
        for (int j = 0; j < nnode; j++)
        {
            c_dump << j << " " << IonicStrengths[j] << " ";
            for (i = 0; i < ncomp; i++)
                c_dump << ActivityCoefficients[j][i] << " ";
            c_dump << "\n" << flush;
        }
        c_dump.close();
    }
    if (debugoutflag && plot == true)
    {
        c_dumpfilename = FileName + "_EquilibriumConstants.dump";
        c_dump.setf(ios::scientific, ios::floatfield);
        c_dump.precision(12);
        c_dump.open(c_dumpfilename.c_str());
        // header
        c_dump << "Node ";
        for (i = 0; i < this->NumberReactions; i++)
        {
            m_kr = KinReact_vector[i];
            if (m_kr->typeflag_mineralkinetics != 0)
                c_dump << m_kr->getName() << " ";
        }
        c_dump << "\n";
        // data
        for (int j = 0; j < nnode; j++)
        {
            c_dump << j << " ";
            for (i = 0; i < this->NumberReactions; i++)
            {
                m_kr = KinReact_vector[i];
                if (m_kr->typeflag_mineralkinetics != 0)
                {
                    if (m_kr->Km_uniform == false)
                        c_dump << m_kr->Km[j] << " ";
                    else
                        c_dump << m_kr->Km[0] << " ";
                }
            }
            c_dump << "\n" << flush;
        }
        c_dump.close();
    }

    return;
}

/**************************************************************************
 Reaction-Method:
 Task: This function updates mineral reactive surface areas after reactions
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void CKinReactData::PostprocessMinKin()
{
    int i, k, nnode;
    double conc0 = 0;
    double phi0, phi1, n0, n1, conc1, mv;
    bool outflag = false;

    phi0 = phi1 = n0 = n1 = conc0 = mv = 0;

    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    nnode = (long)m_msh->nod_vector.size();
    CKinReact* m_kr = NULL;

    // update or reset mineral reactive surface areas for next time step
    for (i = 0; i < NumberReactions; i++)
    {
        m_kr = KinReact_vector[i];
        if (m_kr->typeflag_mineralkinetics && m_kr->Am_constant == false)
        {
            outflag = true;
            for (k = 0; k < nnode; k++)
            {
                // this function is executed before new time level Cs were
                // copied to old time level so access the old time level
                // concentrazion here by index-1 conc0 =
                // pcs_vector[sp_pcsind[m_kr->mineral_number]]->GetNodeValue(k,sp_varind[m_kr->mineral_number]-1);
                // conc1 =
                // pcs_vector[sp_pcsind[m_kr->mineral_number]]->GetNodeValue(k,sp_varind[m_kr->mineral_number]);
                conc0 =
                    cp_vec[m_kr->mineral_number]->getProcess()->GetNodeValue(
                        k, sp_varind[m_kr->mineral_number] - 1);
                conc1 =
                    cp_vec[m_kr->mineral_number]->getProcess()->GetNodeValue(
                        k, sp_varind[m_kr->mineral_number]);
                // CB HS update

                switch (m_kr->Am_model)
                {
                    case 1:  // restore A0/C0 from A0*C(t)/C0 by division with
                             // C(t) after transport
                        if (conc0 > 0)
                            m_kr->Am[k] /= conc0;
                        else
                            m_kr->Am[k] =
                                1e-6;  // minimum allowed surface area to allow
                                       // regrowth of a mineral
                        break;
                    case 2:  // for primary minerals, volume fraction weighted
                        if (REACTINT_vec.size() > 0)
                        {
                            n0 = REACTINT_vec[0]
                                     ->node_ini_porosity[k];  // here plug in
                                                              // n_ini at node
                            n1 = REACTINT_vec[0]->node_porosity[k];
                        }
                        else
                            n0 = n1 = m_kr->GetPhaseVolumeAtNode(
                                k, 1, 0);  // no poro update when RACINT does
                                           // not exist
                        if (conc0 > conc1)
                        {  // dissolution
                            mv = cp_vec[m_kr->mineral_number]->molar_weight /
                                 cp_vec[m_kr->mineral_number]->mineral_density;
                            phi0 = m_kr->Cminini[k] * (1 - n0) * mv;
                            phi1 = conc1 * (1 - n1) * mv;
                        }
                        else  // precipitation
                            phi0 = phi1 = 1;
                        m_kr->Am[k] = m_kr->Am_ini *
                                      pow((phi1 / phi0 * n1 / n0), 0.66667);
                        break;
                    case 3:  // for secondary minerals, only poro fraction
                             // weighted
                        if (REACTINT_vec.size() > 0)
                        {
                            n0 = REACTINT_vec[0]
                                     ->node_porosity[k];  // here plug in n_ini
                                                          // at node
                            n1 = REACTINT_vec[0]->node_porosity[k];
                        }
                        else
                            n0 = n1 = m_kr->GetPhaseVolumeAtNode(
                                k, 1, 0);  // no poro update when RACINT does
                                           // not exist
                        if (conc0 > conc1)
                        {  // dissolution
                            mv = cp_vec[m_kr->mineral_number]->molar_weight /
                                 cp_vec[m_kr->mineral_number]->mineral_density;
                            phi1 = conc1 * (1 - n1) * mv;
                        }
                        else  // precipitation
                            phi1 = 1;
                        m_kr->Am[k] =
                            m_kr->Am_ini * pow(phi1 * n1 / n0, 0.66667);
                        break;
                    case 4:  // for secondary minerals, specific Area
                             // dissolution model
                        if (REACTINT_vec.size() > 0)
                        {
                            n0 = REACTINT_vec[0]
                                     ->node_ini_porosity[k];  // here plug in
                                                              // n_ini at node
                            n1 = REACTINT_vec[0]->node_porosity[k];
                        }
                        else
                            n0 = n1 = m_kr->GetPhaseVolumeAtNode(
                                k, 1, 0);  // no poro update when RACINT does
                                           // not exist
                        if (conc0 > conc1)
                        {  // dissolution
                            mv = cp_vec[m_kr->mineral_number]->molar_weight /
                                 cp_vec[m_kr->mineral_number]->mineral_density;
                            phi1 = conc1 * (1 - n1) * mv;
                            m_kr->Am[k] =
                                m_kr->Am_ini * phi1 *
                                cp_vec[m_kr->mineral_number]->molar_weight / mv;
                        }
                        else
                        {  // precipitation as model 3
                            phi1 = 1;
                            m_kr->Am[k] =
                                m_kr->Am_ini * pow(phi1 * n1 / n0, 0.66667);
                        }
                        break;
                    case 5:
                        // do nothing
                        break;
                    default:
                        // do nothing
                        break;
                }
            }
        }
    }

    CTimeDiscretization* m_tim = NULL;
    m_tim = time_vector[0];
    int steps = int(m_tim->time_step_vector.size());
    bool plot = false;
    if (aktueller_zeitschritt == 1 || aktueller_zeitschritt % 10 == 0 ||
        steps == int(aktueller_zeitschritt))
        plot = true;

    // debug output
    if (debugoutflag && outflag && plot)
    {
        ofstream aus;
        string file = FileName + "_reactivesurfaceareas.dump";
        aus.setf(ios::scientific, ios::floatfield);
        aus.precision(12);
        aus.open(file.c_str());

        for (i = 0; i < nnode; i++)
        {
            if (i == 0)
            {
                aus << "node ";
                for (int j = 0; j < NumberReactions; j++)
                {
                    m_kr = KinReact_vector[j];
                    if (m_kr->typeflag_mineralkinetics &&
                        m_kr->Am_constant == false)
                        aus << m_kr->mineral_name << " ";
                }
                aus << "\n";
            }
            aus << i << " ";
            // data
            for (int j = 0; j < NumberReactions; j++)
            {
                m_kr = KinReact_vector[j];
                if (m_kr->typeflag_mineralkinetics &&
                    m_kr->Am_constant == false)
                    aus << m_kr->Am[i] << " ";
            }
            aus << "\n";
        }
        aus.close();
    }

    return;
}

/**************************************************************************
 Reaction-Method:
 Task: This function calculates the IAP/Kequi term of the Lasaga rate law
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double CKinReact::Omega(double* c, int node)
{
    // if(node == 68)
    //  cout << node ;
    int i, idx;
    double omega = 1;
    double K_pt = 0;
    double conc = 0;
    double unitfactor_l = 1;

    CKinReactData* m_krd = NULL;
    m_krd = KinReactData_vector[0];

    if (REACTINT_vec.size() > 0)
    {
        REACTINT* m_rei = NULL;
        m_rei = REACTINT_vec[0];
        if (m_rei->unitconversion)
            unitfactor_l = MOLH2OPERKG / m_rei->water_conc[node];
    }
    // get K equi at given P and T
    if (Km_uniform == true)
        K_pt = Km[0];  // KR variable
    else if (Km_uniform == false)
        K_pt = Km[node];  // KR variable

    // ion activity product
    // -ProductionStoch[idx] is exponent in IAP; sign must be changed here
    // negative, because the kinetic reaction is formulated for precipitation
    // --> removal of dissolved species from solution --> negative StoechCoeff,
    // but they are positive in IAP for educts of a mineral
    for (i = 0; i < number_reactionpartner; i++)
    {
        idx = minSpeciesIdx[i];
        conc = c[idx + 1];
        if (c[idx + 1] < 0)
            conc = 1e-20;
        if (idx == mineral_number)
            omega /= 1;  // Mineral activity = 1 by definition
        else if (idx == water_number)
        {
            if (m_krd->activity_model >
                0)  // in case of water, this contains the activity directly
                omega *= pow(m_krd->ActivityCoefficients[node][idx],
                             -ProductionStoch[idx]);
            else
                omega /= 1;
        }
        else
        {
            // if(reactionpartner[i].compare("pH")!=0){
            if (m_krd->activity_model > 0)
                omega *= pow(m_krd->ActivityCoefficients[node][idx] * conc *
                                 unitfactor_l,
                             -ProductionStoch[idx]);
            else
                omega *= pow(1 * conc * unitfactor_l, -ProductionStoch[idx]);
            //}
            // else // pH is -log(gamma*H+)
            //  omega *= pow( exp(-conc), -ProductionStoch[idx]);
        }
    }
    // divide by equilibrium constant
    if (number_reactionpartner == 0)
        omega = 0;
    else
        omega /= K_pt;
    return omega;
}

/**************************************************************************
 Reaction-Method:
 Task: This function calculates the mineral dissolution / precipitation
    rate as a function of the rate constant Ktot, the surface area
    and the IAP
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double CKinReact::MinRate(double* c, int node, double dt, double Om)
{
    double rate = 0, maxrate = 0;
    double ktot = 0, conc = 0, poro1 = 1, poro2 = 1;
    double /*Om,*/ H1, H2, Betrag, vol, A = 0;
    int i = 0;

    // get Surface Area
    if (Am_constant)
        A = Am[0];
    else
        A = Am[node];

    // IAP
    // Om = Omega(c, node);
    H1 = pow(Om, Theta);
    // DMIN(H1, 11.0);
    // total rate
    if (H1 < 1 &&
        c[mineral_number + 1] <=
            1.1e-7)  // undersaturated, but no mineral present, no dissolution
        rate = 0;
    else
    {  // supersat.-->precip, or undersat. & mineral present
        H2 = fabs(1 - H1);
        if (H2 < OmegaThreshhold)  // very close to equilibrium, do nothing
            rate = 0;
        else
        {
            Betrag = pow(H2, Eta);
            vol = GetReferenceVolume(mineral_number, node);
            if (ExplicitMinKinRateCoeff)
                ktot = this->MinKinRateCoeff;
            else
                ktot = -MinRateConstant(
                    c, node,
                    (1 - H1));  // default negative rate: dissolution of mineral
            // sign function:, positive rate for supersat, and slow it down by
            // precipfactor (defaulfactor = 1)
            if ((1 - H1) <= 0)
                rate = -1.0 * precipfactor * ktot * A /
                       vol;  // precipfactor used to slow down precipitation by
                             // constant factor
            else
                rate = ktot * A / vol;
            // some modification for Lagneau Benchmark probably required here

            // Lagneau Benchmark
            // bool lagneau = false;
            if (lagneau)
                rate = -1.0 * A * vol * ktot / vol;
            // if(node==1) cout << " Attention: Lagneau-Benchmark with modified
            // mineral rate";

            // in case of strong supersaturation, reduce precipitation rate by
            // root function
            if ((1 - H1) <= -1.0)
                rate *= pow(Betrag, 1 / precipexponent);
            else
                rate *= Betrag;
        }
    }

    // limit rates to the present amount of mineral and dissolved species to
    // avoid overshoot
    if (scale_rate && rate != 0)
    {
        maxrate = rate;
        poro1 = GetReferenceVolume(mineral_number,
                                   node);  // mineral (solid phase volume)
        if (rate < 0)
        {  // dissolution
            // only limited by left hand side reaction equation species, i.e.
            // Stcoef>0 or mineral itself
            for (i = 0; i < int(ProductionStoch.size()); i++)
            {
                if (ProductionStoch[i] > 0 || i == mineral_number)
                {
                    conc = DMAX(0, c[i + 1]);  // do not evaluate negative conc
                    if (mineral_number == i)   // this is the mineral
                        maxrate =
                            -conc / dt;  // max amount that can be dissolved
                    else if (ProductionStoch[i] > 0)
                    {  // this is other species consumed by dissolution, Stcoef
                       // > 0
                        poro2 = GetReferenceVolume(
                            i, node);  // dissolved phase volume
                        maxrate =
                            -conc / dt * poro2 / poro1 / ProductionStoch[i];
                    }
                    // update maxrate,the larger, the more negative, so find
                    // maximum here
                    rate = DMAX(rate, maxrate);
                }
            }
        }
        else if (rate > 0)
        {  // precipitation, not limited by mineral concentration
            for (i = 0; i < int(ProductionStoch.size()); i++)
            {
                if (ProductionStoch[i] < 0)
                {  // this is species consumed by precipitation, Stcoef < 0
                    conc = DMAX(0, c[i + 1]);  // do not evaluate negative conc
                    poro2 =
                        GetReferenceVolume(i, node);  // dissolved phase volume
                    maxrate =
                        -conc / dt * poro2 / poro1 /
                        ProductionStoch[i];  // - before conc is necessary here
                }
                // update maxrate, the larger, the more positive, so find
                // minimum here
                rate = DMIN(rate, maxrate);
            }
        }
    }

    // if(scale_rate){
    //  maxrate = rate;
    //  poro1 = GetReferenceVolume(mineral_number,node); // mineral (solid phase
    //  volume) if(rate<0){ // dissolution
    //    //only limited by left hand side reaction equation species, i.e.
    //    Stcoef>0 or mineral itself for(i=0;i<int(minSpeciesIdx.size());i++){
    //      conc = DMAX(0,c[minSpeciesIdx[i]+1]);         // do not evaluate
    //      negative conc if(mineral_number==minSpeciesIdx[i])          // this
    //      is the mineral
    //        maxrate = -conc / dt;                      // max amount that can
    //        be dissolved
    //      else if(ProductionStoch[minSpeciesIdx[i]]>0){ // this is other
    //      species consumed by dissolution, Stcoef > 0
    //        poro2	= GetReferenceVolume(minSpeciesIdx[i],node);         //
    //        dissolved phase volume maxrate = -conc / dt * poro2 / poro1 /
    //        ProductionStoch[minSpeciesIdx[i]];
    //      }
    //      // update maxrate,the larger, the more negative, so find maximum
    //      here rate = DMAX(rate,maxrate);
    //    }
    //  }
    //  else if (rate>0){ // precipitation, not limited by mineral concentration
    //    for(i=0;i<int(minSpeciesIdx.size());i++){
    //      conc = DMAX(0,c[minSpeciesIdx[i]+1]);         // do not evaluate
    //      negative conc if(ProductionStoch[minSpeciesIdx[i]]<0){      // this
    //      is species consumed by precipitation, Stcoef < 0
    //        poro2	= GetReferenceVolume(minSpeciesIdx[i],node);         //
    //        dissolved phase volume maxrate = -conc / dt * poro2 / poro1 /
    //        ProductionStoch[minSpeciesIdx[i]]; // - before conc is necessary
    //        here
    //      }
    //      // update maxrate, the larger, the more positive, so find minimum
    //      here rate = DMIN(rate,maxrate);
    //    }
    //  }
    //}
    if (rate > 0 && rate < 1e-30)
        rate = 0;
    if (rate < 0 && rate > -1e-30)
        rate = 0;
    return rate;
    // return rate/double(KinReactData_vector[0]->NumberMineralkinetics);
}

/**************************************************************************
 Reaction-Method:
 Task: This function resturns the mineral disoolution / precipitation
    rate constant as a sum over all reaction mechanisms
    Ktot = Sum(Mech)
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double CKinReact::MinRateConstant(double* c, int node, double disso)
{
    int i;
    double ktot = 0;

    // calc the rate constant
    for (i = 0; i < number_Mech; i++)
    {
        if (disso < 0 && precip_baseterm_only == true)
        {  // for precipitation, consider only base mechanism, if set in
           // inputfile
            if (mechvec[i]->no_mechSpec == 0)
                ktot += mechvec[i]->Mech(c, node);
        }
        else
        {  // consider all mechanisms
            ktot += mechvec[i]->Mech(c, node);
        }
    }

    return ktot;
}

/**************************************************************************/
// Functions for class MinKinMech
/**************************************************************************/

/**************************************************************************
 Reaction-Method:
 Task: Constructor of MinkinMech class
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
MinkinMech::MinkinMech(void)
{
    no_mechSpec = 0;
    mechSpeciesIdx.clear();
    mechSpeciesNames.clear();
    mechSpeciesExpo.clear();
    Eact = -1;
    k25 = -1;
}
MinkinMech::~MinkinMech(void) {}

/**************************************************************************
 Reaction-Method:
 Task: This function returns the contribution of a single mineral
    precipitation / dissolution mechanism to the total rate constant
    Ktot = Sum(Mech)
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double MinkinMech::Mech(double* c, int node)
{
    // the base term is treated as a mechanism as well
    // Todo: add some function to get T

    int i, idx, actmodel;
    double kmech = 0, conc = 0;
    double R = 8.314472;  // gas constant = R 8.314472(15) J/K/mol
    CKinReactData* m_krd = NULL;
    m_krd = KinReactData_vector[0];
    double unitfactor_l = 1;
    double T = 298.15;

    if (REACTINT_vec.size() > 0)
    {
        REACTINT* m_rei = NULL;
        m_rei = REACTINT_vec[0];
        if (m_rei->unitconversion)
            unitfactor_l = MOLH2OPERKG / m_rei->water_conc[node];
        T = m_rei->GetTemperature(node);
    }
    actmodel = m_krd->activity_model;

    // Todo: add some function to get T from heat transport pcs

    // Arrhenius
    kmech = k25 * exp(-Eact / R * (1 / T - 1 / 298.15));
    // multiply by mechanism activities term
    for (i = 0; i < no_mechSpec; i++)
    {
        idx = mechSpeciesIdx[i];
        conc = c[idx + 1];
        if (conc < 0)
            conc = 1.0e-20;
        if (mechSpeciesNames[i].compare("pH") != 0)
        {
            if (actmodel > 0)
                kmech *= pow((m_krd->ActivityCoefficients[node][idx] * conc *
                              unitfactor_l),
                             mechSpeciesExpo[i]);
            else
                kmech *= pow(
                    (1 * conc * unitfactor_l),
                    mechSpeciesExpo[i]);  // all activity coefficients are 1.0
        }
        else
        {  // pH is -log(gamma*H+)
            if (actmodel > 0)
                kmech *= pow(pow(10, m_krd->ActivityCoefficients[node][idx] *
                                         (-conc * unitfactor_l)),
                             mechSpeciesExpo[i]);
            else
                kmech *= pow(pow(10, 1 * (-conc * unitfactor_l)),
                             mechSpeciesExpo[i]);
        }
    }

    return kmech;
}

/**************************************************************************
 Reaction-Method:
 Task: Copy concentrations after reaction step to symmetric nodes
        Implemented only for radial models with one hex element in theta
 direction Then it is sufficient to simulate flow and conservative transport at
 all nodes, but reactions only at one side of the element. Then copy to other
 side of element. Programing:
 //SB 01.2011 SB First implementation
 **************************************************************************/
void CKinReactData::CopyConcentrations(void)
{
    long i, l, k, nnodes, nele;
    size_t j;
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    CRFProcess* m_pcs = NULL;
    // CNode* m_nod0 = NULL, *m_nod1 = NULL, *m_nod2 = NULL, *m_nod3 = NULL;
    MeshLib::CElem* m_ele = NULL;
    vec<MeshLib::CNode*> ele_nodes(8);
    Math_Group::vec<long> vec_nod_index(8);
    const double *coord0, *coord1, *coord2, *coord3;
    bool err_flag =
        true;  // if true, all is OK, if false, do not use this option
    int copy_direction = -1;
    long i_copy_from, i_copy_to;
    int nidx0;
    long nidx = 0;

    if (copy_concentrations == false)
        return;

    cout << " CopyConcentrations "
         << "\n";
    // get number of nodes
    nnodes = (long)m_msh->nod_vector.size();
    nele = (long)m_msh->ele_vector.size();

    // first time called, Check if mesh is appropriate
    if (copy_nodes.size() < 1)
    {
        // BATCH
        if (batch)
        {
            if (nnodes == 2 && nele == 1 &&
                m_msh->ele_vector[0]->GetElementType() == 1)
            {
                // for a single line element
                if (is_a_CCBC[0] == true && is_a_CCBC[1] == false)
                    batch = true;
                else if (is_a_CCBC[0] == false && is_a_CCBC[1] == true)
                    batch = true;
                else
                {
                    cout << "\n"
                         << " Copy Concentrations only for a single (batch) "
                            "line element with one deactivated node!"
                         << " - Option switched off "
                         << "\n";
                    batch = copy_concentrations = false;
                }
            }
        }
        // RADIAL
        else if ((radial) || (TwoDinThreeD))
        {
            // Only for hex-elements
            for (i = 0; i < nele; i++)
            {
                if (m_msh->ele_vector[i]->GetElementType() != 3)
                {
                    cout << "\n"
                         << " Copy Cons only for hex elements ! - Option "
                            "switched off "
                         << "\n";
                    radial = TwoDinThreeD = copy_concentrations = false;
                    break;
                }
            }

            // Now test for radial flow model along x-Axis
            // Test is if model is symmetric to x-Axis; not tested is, if it is
            // one column
            m_ele = m_msh->ele_vector[0];  // get first element
            m_ele->GetNodes(ele_nodes);    // get element nodes
            coord0 = ele_nodes[0]
                         ->getData();  // Get coordinates of first four element
                                       // nodes ("upper" part of element)
            coord1 = ele_nodes[1]->getData();  // CB 03/11
            coord2 = ele_nodes[2]->getData();
            coord3 = ele_nodes[3]->getData();

            /*       *3-----*2
                     |\     |\
                     | *0---+-*1
                     *-|----* |
                      \|     \|
                       *------*        */

            // check orientation of element: all should have same z-value
            if ((coord0[2] != coord1[2]) || (coord0[2] != coord2[2]) ||
                (coord0[2] != coord3[2]))
            {
                cout << " Error: Not identical z coordinates of element nodes "
                        "of element 0 "
                     << "\n";
                err_flag = false;
            }
            // node0 and node3 should have same x value, as well as node1 and
            // node2
            if ((coord0[0] != coord3[0]) || (coord1[0] != coord2[0]))
            {
                cout << " Error: Not identical x coordinates of element nodes "
                        "(0/3) or (1/2) of element 0 "
                     << "\n";
                err_flag = false;
            }
            // this is different for radial and TwoDinThreeD
            if (radial)
            {
                // node0 and node3 should have same y=-y value, as well as node1
                // and node2
                if ((coord0[1] != -1.0 * coord3[1]) ||
                    (coord1[1] != -1.0 * coord2[1]))
                {  // this is just for radial flow models
                    cout << " Error: Not y_i = -y_j for element nodes (0/3) or "
                            "(1/2) of element 0 "
                         << "\n";
                    err_flag = false;
                }
            }
            else if (TwoDinThreeD)
            {
                // node0 and node1 should have same y value, as well as node2
                // and node3
                if ((coord0[1] != coord1[1]) || (coord2[1] != coord3[1]))
                {  // this is for regular hex grid flow models
                    cout << " Error: Not y_i = -y_j for element nodes (0/3) or "
                            "(1/2) of element 0 "
                         << "\n";
                    err_flag = false;
                }
            }
            // everything fine so far?
            if (err_flag == false)
            {
                cout << "\n"
                     << " Copy Cons - Grid not appliccable ! - Option switched "
                        "off "
                     << "\n";
                radial = TwoDinThreeD = copy_concentrations = false;
            }

            // Now check, if the correct nodes are switched off and determine
            // copy direction, only for radial
            err_flag = true;  // reset for reuse
            if (radial)
            {
                if (is_a_CCBC[ele_nodes[0]->GetIndex()] == true)
                {  // element node 0 is not calculated
                    if (is_a_CCBC[ele_nodes[1]->GetIndex()] != true)
                        err_flag = false;
                    if (is_a_CCBC[ele_nodes[2]->GetIndex()] != false)
                        err_flag = false;
                    if (is_a_CCBC[ele_nodes[3]->GetIndex()] != false)
                        err_flag = false;
                    copy_direction = 1;
                }
                if (is_a_CCBC[ele_nodes[3]->GetIndex()] == true)
                {  // element node 0 is not calculated
                    if (is_a_CCBC[ele_nodes[2]->GetIndex()] != true)
                        err_flag = false;
                    if (is_a_CCBC[ele_nodes[1]->GetIndex()] != false)
                        err_flag = false;
                    if (is_a_CCBC[ele_nodes[0]->GetIndex()] != false)
                        err_flag = false;
                    copy_direction = 0;
                }
                // everything fine so far?
                if (err_flag == false)
                {
                    cout << "\n"
                         << " Copy Cons - Please switch off one model side for "
                            "reaction calculation ! - Option "
                            "switched off "
                         << "\n";
                    radial = TwoDinThreeD = copy_concentrations = false;
                }
            }
        }  // ((radial) || (TwoDinThreeD)){

    }  // end of mesh and geometry checks

    // now quit, if somethings wrong
    if (copy_concentrations == false)
        return;

    // first time called, build data structures and set up index vector
    if (copy_nodes.size() < 1)
    {
        // Make vector for all nodes for mapping
        for (i = 0; i < nnodes; i++)
            copy_nodes.push_back(i);  // map each node to itself by default

        // Now, prepare copy node relations
        if (batch)
        {
            for (i = 0; i < nnodes; i++)
                if (is_a_CCBC[i] == true && i > 0)
                    copy_nodes[i] = i - 1;
                else if (is_a_CCBC[i] == true)
                    copy_nodes[i] = i + 1;
        }  // batch
        else if (radial)
        {
            // Go through all elements and set copy values respectively
            for (i = 0; i < nele; i++)
            {
                m_ele = m_msh->ele_vector[i];  // get element
                m_ele->GetNodeIndeces(vec_nod_index);
                if (copy_direction == 1)
                {  // element nodes 0, 1, 4, 5 not calculated
                    // copy concentrations 3 -> 0
                    i_copy_from = vec_nod_index[3];
                    i_copy_to = vec_nod_index[0];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 2 -> 1
                    i_copy_from = vec_nod_index[2];
                    i_copy_to = vec_nod_index[1];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 7 -> 4
                    i_copy_from = vec_nod_index[7];
                    i_copy_to = vec_nod_index[4];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 6 -> 5
                    i_copy_from = vec_nod_index[6];
                    i_copy_to = vec_nod_index[5];
                    copy_nodes[i_copy_to] = i_copy_from;
                }
                if (copy_direction == 0)
                {  // element nodes 2, 3, 6, 7 not calculated
                    // copy concentrations 0 -> 3
                    i_copy_from = vec_nod_index[0];
                    i_copy_to = vec_nod_index[3];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 1 -> 2
                    i_copy_from = vec_nod_index[1];
                    i_copy_to = vec_nod_index[2];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 4 -> 7
                    i_copy_from = vec_nod_index[4];
                    i_copy_to = vec_nod_index[7];
                    copy_nodes[i_copy_to] = i_copy_from;
                    // copy concentrations 5 -> 6
                    i_copy_from = vec_nod_index[5];
                    i_copy_to = vec_nod_index[6];
                    copy_nodes[i_copy_to] = i_copy_from;
                }
                if (i == 0)
                {
                    cout << "\n"
                         << "copy_nodes[i] = "
                         << "\n";
                    for (int ii = 0; ii < 8; ii++)
                        cout << "          " << ii << " = " << copy_nodes[ii]
                             << "\n";
                }
            }
        }  // radial
        else if (TwoDinThreeD)
        {
            for (i = 0; i < nnodes; i++)
            {
                if (is_a_CCBC[i] == false)
                {  // reactions are quantified at this node
                    // get node coordinates
                    coord0 = m_msh->nod_vector[i]->getData();
                    // check all neighbour nodes
                    for (j = 0;
                         j < m_msh->nod_vector[i]->getNumConnectedNodes();
                         j++)
                    {
                        nidx = m_msh->nod_vector[i]->getConnectedNodes()[j];
                        coord1 = m_msh->nod_vector[nidx]->getData();
                        // check orientation of nodes: should have same x + z
                        // and y of neighbour should be larger
                        if ((coord0[0] == coord1[0]) &&
                            (coord0[2] == coord1[2]) && (coord0[1] < coord1[1]))
                        {
                            // Found!
                            // Check if reaction is really switched off at this
                            // neighbour
                            if (is_a_CCBC[nidx] ==
                                true)  // set the neighbour index
                                // copy_nodes[i]=nidx;
                                copy_nodes[nidx] = i;
                            break;
                        }
                    }  // done
                    // now check for success:
                    if (copy_nodes[nidx] == nidx)
                        cout << "Warning in CKinReactData::CopyConcentrations:!"
                             << " No Neighbour found for copy nodes " << i
                             << "\n";
                }  // is_a_CCBC[i] == false -> node with reactions
            }      // node loop
        }          // TwoDinThreeD

    }  // if first time

    // Now, Copy concentrations in all cases
    // Go through all transport processes and copy values for all nodes
    // Where to put this? -> ToDo CB in new concept
    // this is the same for batch or radial or 2Din3D
    for (i = 0; i < (int)pcs_vector.size(); i++)
    {
        m_pcs = pcs_vector[i];
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {  // if it is a MASS_TRANSPORT process
            l = 0;
            // get index of concentration
            for (j = 0; j < (size_t)m_pcs->pcs_number_of_primary_nvals; j++)
            {
                nidx0 = m_pcs->GetNodeValueIndex(
                            m_pcs->pcs_primary_function_name[j]) +
                        1;  // new time level
                for (k = 0; k < nnodes; k++)
                {
                    i_copy_to = k;
                    i_copy_from = copy_nodes[i_copy_to];
                    if (i_copy_to != i_copy_from)
                    {  // only do for differing node indicess
                        m_pcs->SetNodeValue(
                            i_copy_to, nidx0,
                            m_pcs->GetNodeValue(i_copy_from, nidx0));
                        l++;
                    }
                }
            }
        }
    }
    cout << "     at " << l << " of " << nnodes << " nodes."
         << "\n";
    // end function
}

// this is a "slow" bubblesort, but it gives back the indexes of our sort
void CKinReactData::SortIterations(long* iterations, long* indexes, long len)
{
    long i, j, iTmp;
    // Set the index to match the current sort
    // order of the names
    for (i = 0; i < len; i++)
        indexes[i] = i;

    // Bubblesort. Sort the indexes
    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
        {
            if (iterations[indexes[i]] < iterations[indexes[j]])
            {
                // swap the indexes
                iTmp = indexes[i];
                indexes[i] = indexes[j];
                indexes[j] = iTmp;
            }
        }
    }
}
