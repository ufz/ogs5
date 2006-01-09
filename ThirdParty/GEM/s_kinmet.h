//-------------------------------------------------------------------
// $Id: s_kinmet.h 725 2012-10-02 15:43:37Z kulik $
//
// Stub Declaration of TKinMet class for further development of 
//  kinetics/metastability models
//
// Copyright (C) 2012  D.Kulik, B.Thien, G.Kosakowski
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------
//
#ifndef S_KINMET_H
#define S_KINMET_H

//#include "s_fgl.h"

const int   MAXDCNAME =      16, MAXPHNAME = 16;

struct KinMetData {

    char  KinModCod_;   // Type code of the kinetic/metastability model, see KINMETTYPECODES
    char  UptModCod_;   // Type code of Tr uptake model (solution/sorption phases only), see UPTMTYPECODES
    char  PhasNam_[MAXPHNAME];      // Phase name (for specific built-in models)

    long int NSpec_; // Number of components in the phase
    long int nlPh_;  // number of linked phases (cf. lPh), default 0.
    long int nFace_; // number of kinetically different faces (on solid phase surface), default 1, up to 6.
    long int nReg_;  // number of kinetic regions (and catalyzing aqueous species)
    long int nrpC_;  // number of kinetic rate constants and coefficients
    long int numpC_; // number of uptake model parameter coefficients (per end member)
    long int kins_;  // kinetic state: -1 dissolution +1 precipitation ...
//    long int Nres_;  // reserved

    double T_k_;         // Temperature, K (initial)
    double P_bar_;       // Pressure, bar (initial)
    double kTau_;        // current time, s (initial)
    double kdT_;         // current time step (initial)
    double IS_;          // Effective molal ionic strength of aqueous electrolyte
    double pH_;          // pH of aqueous solution
    double pe_;          // pe of aqueous solution
    double Eh_;          // Eh of aqueous solution, V

    double sSA_;  // Specific surface area of the phase, m2/g, default: 0.
    double sgw_;  // Standard mean surface energy of solid-aqueous interface, J/m2
    double sgg_;  // Standard mean surface energy of gas-aqueous interface, J/m2
    double rX0_;  // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
    double hX0_;  // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
    //
    double sVp_;  // Specific pore volume of phase, m3/g (default: 0)
    double frSA_; // reactive fraction of surface area (def. 1)
    double nPh_;  // current amount of this phase, mol (read-only)
    double mPh_;  // current mass of this phase, g (read-only)
    double vPh_;  // current volume of this phase, cm3 (read-only)
    double sAPh_;  // current surface of this phase, m2
    double OmPh_;  // phase stability index (lg scale), input

    double *nPul_; // upper restriction to this phase amount, mol (calculated here)
    double *nPll_; // lower restriction to this phase amount, mol (calculated here)
    double *sGP_;  // surface free energy of the phase, J (YOF*PhM)

    double *fSAk;     // Pointer to input fractions of surface area of solid occupied by different faces [nFaces]
    double *arKrpc_;  // pointer to input array of kinetic rate constants [nReg*nFaces][nrpC]
    double *arUmpc_;  // pointer to input array of uptake model coefficients [nComp][numpC]

    char  (*SM_)[MAXDCNAME];  // pointer to the list of DC names in the phase [NComp] read-only
    char  arDCC_;     // pointer to the classifier of DCs involved in sorption phase [NComp] read-only

    long int *xFaces;    // indexes of faces per set of kinetic region parameters [nReg*nFaces], default 0.
    long int *arjCrDC_;  // input array of DC indexes used in rate regions [nReg*nFaces], default 0.

    double *arym_;    // Pointer to molalities of all species in MULTI (provided), read-only
    double *arla_;    // Pointer to lg activities of all species in MULTI (provided), read-only

    double *arnx_;    // Pointer to mole amounts of phase components (provided) [NComp] read-only

    double *arnxul_;     // Vector of upper kinetic restrictions to nx, moles [L]  (DUL) direct access output
    double *arnxll_;     // Vector of lower kinetic restrictions to nx, moles [L]  (DLL) direct access output

    double *arWx_;       // Species (end member) mole fractions ->NSpecies
    double *arVol_;      // molar volumes of end-members (species) cm3/mol ->NSpecies

};

// Class describing a mineral face precipitation/dissolution data (Lasaga's general form)
class TFaceKinetics
{
protected:

     long int nReg; // number of kinetic regions applying to this face
     long int nrC;  // number of kinetic rate constants and coefficients (max. 32)
     long int jrDC; // index of catalytic (rate-controlling) species
     long int iRes; // reserved

     double
        nDC[], // amounts of catalytic species [nReg]
        cDC[], // concentrations (molality) of catalytic species [nReg]
        aDC[],   // activity (fugacity) of catalytic species [nReg]

        krc[][], // kinetic rate constants (used differently in subclasses) [nReg][nrC]

        fSAf, // fraction of surface area under this face (default 1)

        aft[],  // affinity term (real scale) [nReg]
        arf[],  // Arrhenius factor [nReg]
        cat[],  // catalytic term  [nReg]

        kc[],   // rate constant (involving all corrections) in mol/m2/s [nReg]
        rr[],   // rate for this region (output) [nReg]
        rmol,   // rate for the whole face (output) in mol/s
        velo,   // velocity of face growth (positive) or dissolution (negative) nm/s
        ;

        double RegRateCon( double T_p, long int iReg );  // calculates the region rate constant
        double RegRate( double T_p, long int iReg );  // calculates the region mole rate

public:

    // Generic constructor
    TFaceKinetics( long int nReg_p, long int nrC_p, long int jrDC_p, double* nDC_p, double* cDC_p, double* aDC_p,
              double *krc_p,  double fsaf_p );

    // Destructor
    virtual ~TFaceKinetics();

    double FaceRate( );  // calculates the face mole rate
    double FaceVelo( );  // calculates the face growth/dissolution velocity

};


class TKinMet  // Base class for kinetics and metastability models
{
    protected:
        char  KinMT;   // Type code of the kinetic/metastability model, see KINMETTYPECODES
        char  SpliMd;  // Splitting model code for dissolution
        char  SpliMp;  // Splitting model code for precipitation
        char  UptMT;   // Type code of Tr uptake model (solution/sorption phases only), see UPTMTYPECODES
        char PhaseName[MAXPHASENAME+1];    // Phase name (for specific built-in models)

        long int NComp,       // Number of components in the phase
        nlPh,  // number of linked phases (cf. lPh), default 0.
        nFaces,// number of kinetically different faces (on solid phase surface), default 1, up to 6.
        nReg,  // number of kinetic regions (and catalyzing aqueous species) default 1
        nrpC,  // number of kinetic rate constants and coefficients default 2
        numpC, // number of uptake model parameter coefficients (per end member) default 0
        kinst, // status: 0: equilibrium; -1 dissolution; +1 precipitation; ...
        iRes,  // reserved
        ;
        double
              kTau,        // current time, s (initial)
              kdT,         // current time step (initial)
              Tk,    	    // Temperature, K
              Pbar,  	    // Pressure, bar
              R_CONST,    // Gas constant, J/K/mol
              IS,         // Effective molal ionic strength of aqueous electrolyte
              pH,         // pH of aqueous solution
              pe,         // pe of aqueous solution
              Eh,         // Eh of aqueous solution, V
              ;
        double sSA,  // Specific surface area of the phase, m2/g, default: 0.
               sgw, // Standard mean surface energy of solid-aqueous interface, J/m2
               sgg, // Standard mean surface energy of gas-aqueous interface, J/m2
               rX0,    // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
               hX0,    // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
        //
               sVp,  // Specific pore volume of phase, m3/g (default: 0)
               frSA,   // reactive fraction of surface area (def. 1)

               nPh,  // current amount of this phase, mol (read-only)
               mPh,  // current mass of this phase, g (read-only)
               vPh,  // current volume of this phase, cm3 (read-only)
               sAPh,  // current surface of this phase, m2
               OmPh,  // phase stability index (lg scale), input
            //
               nPul, // upper restriction to this phase amount, mol (calculated here)
               nPll, // lower restriction to this phase amount, mol (calculated here)
               sGP,  // surface free energy of the phase, J (YOF*PhM)
        ;

        double
        fSAf[],    // Pointer to input fractions of surface area of solid occupied by different faces [nFaces]
        arKrpc[][],  // entry pointer to array of kinetic rate constants [nReg*nFaces][nrpC]
        arUmpc[][],  // entry pointer to array of uptake model coefficients [nComp][numpC]
        ;
        char  (*arSM)[MAXDCNAME];  // pointer to the list of DC names in the phase [NComp] read-only
        char  *arDCC;   // pointer to the classifier of DCs involved in this phase [NComp] read-only

        long int xFaces[];   // indexes of faces per set of kinetic region parameters [nReg*nFaces], default 0.
        long int arjCrDC[];  // pointer to input array of DC indexes used in rate regions [nFaces*nReg]

        long int arPhLin[][];  // indexes of linked phases and link type codes [nlPh*2] read-only

        double *arym;     // molalities of all species in MULTI, read-only
        double *arla;     // lg activities of all species in MULTI, read-only (for catalytic terms)
        double *arnx;     // Pointer to mole amounts of phase components (provided) [NComp] read-only
        double *arWx;     // Species (end member) mole fractions (provided) ->NSpecies
        double *arVol;    // molar volumes of end-members (species) (provided) cm3/mol ->NSpecies

        TFaceKinetics arFk[]; // work dyn array of kinetic faces data (model-specific) [nFaces]

        double spcfu[];    // work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]
        double spcfl[];    // work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]

        double *arnxul;    // Vector of upper kinetic restrictions to nx, moles [L]  (DUL) direct access output
        double *arnxll;    // Vector of lower kinetic restrictions to nx, moles [L]  (DLL) output access output

        public:
                // Generic constructor
                TKinMet( KinMetData kmd );

                // Destructor
                virtual ~TKinMet();

                virtual long int SetTime()
                {
                        return 0;
                };

                virtual long int PTparam()
                {
                        return 0;
                };

                virtual long int RateMod()
                {
                        return 0;
                };

                virtual long int SplitMod()
                {
                        return 0;
                };

                virtual long int UptakeMod()
                {
                        return 0;
                };

                // modified faces area fractions fSAf and specific surface area of the phase
                virtual double ModifiedFSA ( double *fSAf_p );

                // set new system state
                long int UpdatePT ( double T_k, double P_bar );

                // set new time and time step
                bool ResetTime( double Tau, double dTau );

                // Checking dimensions and re-allocation, if necessary
                bool testSizes( long int NComps, long int NlPhs, long int NFaces, long int NRegs,
                                long int NrpCs, long int numpCs, char kModCode, char uModCode );

                // getting phase name
                void GetPhaseName( const char *PhName );

};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - derived classes - - - - - - - - - - - - - -

class TPalandri: public TKinMet  // Generic dissolution/precipitation models following Lasaga 1998
{
    protected:

 // specific stuff for Palandri form

        // internal functions
        void alloc_internal();
        void free_internal();

public:
        // Constructor
        TPalandri( KinMetData kmd, /* specific params */ );

        // Destructor
        ~TPalandri();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        // Calculates splitting to end members
        long int SplitMod();

        // Calculates uptake rates
        long int UptakeMod();

};


class TWolthers: public TKinMet  // Uptake model
{
    protected:

 // specific stuff for Palandri form

        // internal functions
        void alloc_internal();
        void free_internal();

public:
        // Constructor
        TWolthers( KinMetData kmd, /* specific params */ );

        // Destructor
        ~TWolthers();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        // Calculates uptake rates
        long int UptakeMod();


};

class TLaidler: public TKinMet  // Laidler d/p and uptake model
{
    protected:

 // specific stuff for Laidler form

        // internal functions
        void alloc_internal();
        void free_internal();

public:

        // Constructor
        TLaidler( KinMetData kmd, /* specific params */ );

        // Destructor
        ~TLaidler();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        // Calculates uptake rates
        long int UptakeMod();


};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif // S_KINMET_H
