//-------------------------------------------------------------------
// $Id: ms_multi_format.cpp 829 2013-05-01 16:14:07Z kulik $
//
/// \file ms_multi_format.cpp
/// Implementation of writing/reading IPM, DCH and DBR text I/O files
//
// Copyright (c) 2006-2012 S.Dmytriyeva,D.Kulik
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
//-------------------------------------------------------------------

#include "io_arrays.h"
#include "m_param.h"
#include <iomanip>

#ifdef IPMGEMPLUGIN
  istream& f_getline(istream& is, gstring& str, char delim);
#endif

//bool _comment = true;
const char *_GEMIPM_version_stamp = " GEMS3K v.3.2 r.830 (rc) ";

//===================================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here. The fourth field contains
// the text of the comment for this data object, optionally written into the
// text-format output IPM file.
//
outField MULTI_static_fields[8] =  {
  { "pa_PE" , 0 , 0, 0, "# PE: Flag for using electroneutrality condition in GEM IPM calculations (1 or 0)" },
  { "PV" ,    0 , 0, 0, "\n# PV: Flag for the volume balance constraint (on Vol IC) for indifferent equilibria at P_Sat (0 or 1)" },
  { "PSOL" ,  0 , 0, 0, "\n# PSOL: Total number of DCs in liquid hydrocarbon phases (0; reserved)" },
  { "PAalp" , 0 , 0, 0, "# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases " },
  { "PSigm" , 0 , 0, 0, "\n# PSigm: Flag for using (+) or ignoring (-) specific surface free energies  " },
  { "Lads" ,  0 , 0, 0, "# Lads: Total number of Dependent Components in sorption phases included into this system" },
  { "FIa" ,   0 , 0, 0, "# FIa: Number of sorption phases included in this system (0 if no sorption phases )" },
  { "FIat" ,  0 , 0, 0, "# FIat: Maximum number of surface types per adsorption phase (if FIa > 0, set FIat = 6)" }
};

outField MULTI_dynamic_fields[70] =  {
// write/read dynamic (array) data to/from the text-format IPM file
   {  "sMod",  1 , 0, 0, "# sMod: Codes for TSolMod built-in  models of mixing in multicomponent phases [nPS*6]" },
   {  "LsMod", 1 , 0, 0, "\n# LsMod: Dimensions of TSolMod <IPxPH> and <PMc> data arrays [nPS*3]. In each row (for phase):"
                         "\n# [0] number of interaction parameters (rows in <IPx>); [1] max. parameter order (columns in <IPx>);"
                         "\n# [2] number of coefficients per interaction parameter in <PMc> array" },
   {  "LsMdc", 1 , 0, 0, "\n# LsMdc: Dimensions of TSolMod <DMc> and <MoiSN> arrays [nPS*3]: In each row (for phase):"
                         "\n# [0] number of parameters per component; [1] 0; [2] 0. For multi-site (sublattice) models: "
                         "\n#   [1] number of sublattices nS; [2] total number of moieties nM acting in sublattice sites" },
   {  "B",     1 , 0, 0, "# B: Full total bulk composition (vector b), moles [nIC] (will be partially re-written from DBR files)" },
   {  "DCCW",  0 , 0, 0, "# DCCW: internal DC class codes [nDC], will be reset automatically from DCH file content" },
               // DCCW is placeholder - something else can be used here, if needed
   {  "Pparc", 0 , 0, 0, "# Pparc: Partial pressures or fugacities of pure Dependent Components [nDC] (reserved)" },
   {  "fDQF",  0 , 0, 0, "\n# fDQF: DQF parameters of end members or pure gas fugacities, (J/mol/(RT) [nDC]" },
   {  "lnGmf", 0 , 0, 0, "\n# lnGmf: Natural logarithms of DC activity coefficients used at Simplex LP approximation only [nDC]" },
   {  "RLC",   0 , 0, 0, "# RLC: Code of metastability constraints for DCs {L U B (default)} [nDC]" },
   {  "RSC",   0 , 0, 0, "\n# RSC: Units of metastability/kinetic constraints for DCs {M} moles [nDC]" },
   {  "DLL",   0 , 0, 0, "\n# DLL: Lower metastability constraints on DC amounts <xDC>, moles [nDC] (default: 0)" },
   {  "DUL",   0 , 0, 0, "\n# DUL: Upper metastability constraints on DC amounts <xDC>, moles [nDC] (default: 1e6)" },
   {  "Aalp",  0 , 0, 0, "# Aalp: Specific surface areas of phases, m2/g [nPH]" },
   {  "Sigw",  0 , 0, 0, "\n# Sigw: Specific surface free energy for phase-water interface, J/m2 [nPH] (reserved)" },
   {  "Sigg",  0 , 0, 0, "\n# Sigg: Specific surface free energy for phase-gas interface, J/m2 (not yet used) [nPH]" },
   {  "YOF",   0 , 0, 0, "\n# YOF: Surface free energy parameter for phases in J/g (to accomodate for variable phase composition)  [nPH]" },
   {  "Nfsp",  1 , 0, 0, "\n# Nfsp: Fractions of the sorbent specific surface area allocated to surface types [nPS*6]" },
   {  "MASDT", 1 , 0, 0, "\n# MASDT: Total maximum site  density per surface type, mkmol/g [nPS*6]" },
   {  "C1",    1 , 0, 0, "\n# C1: Inner capacitance density parameter C1, F/m2 [nPS*6]" },
   {  "C2",    1 , 0, 0, "\n# C2: Outer capacitance density parameter C2, F/m2 [nPS*6]" },
   {  "C3",    1 , 0, 0, "\n# C3: Third capacitance density parameter C3, F/m2 [nPS*6]" },
   {  "pCh",   1 , 0, 0, "\n# pCh: Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [nPS*6]" },
   {  "SATX",  1 , 0, 0, "\n# SATX: Setup of surface sites and species (will be applied separately within each sorption phase) [Lads*4]"
                        "\n# [0] surface type; [1] sorbent emd member; [2] surface site in surf. type; [3] surface EDL plane" },
   {  "MASDJ", 1 , 0, 0, "\n# MASDJ: Parameters of surface species in surface complexation models [Lads*6]"
                         "\n# [0] max site density mmol/g; [1] species charge allocated to 0 plane;"
                         "\n# [2] species charge allocated to beta -or third plane; [3] Frumkin isotherm interaction parameter;"
                         "\n# [4] denticity or coordination number CN; [5] reserved isoterm parameter" },
   {  "SCM",   1 , 0, 0, "\n# SCM: Classifier of built-in electrostatic models applied to surface types in sorption phases [nPS*6]" },
   {  "SACT",  1 , 0, 0, "\n# SACT: Classifier of applied SACT equations (isotherm corrections) [Lads]" },
   {  "DCads", 1 , 0, 0, "\n# DCads: Classifier of DCs involved in sorption phases [Lads]" },
               // static: GEM IPM v3 numerical tolerances
   { "pa_DB" , 0 , 0, 0,  "# DB: Minimum amount of IC in the bulk composition, moles (except charge Zz) { 1e-17 }"},
   { "pa_DHB", 0 , 0, 0,  "\n# DHB: Maximum allowed relative mass balance residual for ICs { 1e-13 } " },
   { "pa_EPS", 0 , 0, 0,  "\n# EPS: Tolerance of the SolveSimplex() balance residual for ICs { 1e-10 } " },
   { "pa_DK",  0 , 0, 0,  "\n# DK: Tolerance for the Dikin's criterion of IPM convergence { 1e-6 } " },
   { "pa_DF" , 0 , 0, 0,  "\n# DF: Tolerance DF of the stability criterion for a lost phase to be inserted to mass balance { 0.01 } " },
   { "pa_DP",  0 , 0, 0,  "\n# DP: Maximal number of iterations in MassBalanceRefinement MBR() procedure { 130 }"  },
   { "pa_IIM", 0 , 0, 0,  "\n# IIM: Maximum allowed number of iterations in one main GEM IPM descent run { 7000 }" },
   { "pa_PD" , 0 , 0, 0,  "\n# PD: Mode of calculation of DC activity coefficients { 2 } " },
   { "pa_PRD" , 0 , 0, 0, "\n# PRD: Disable (0) or activate (-4 or less- max.dec.exp.for DC amount correction) SpeciationCleanup() { -5 }" },
   { "pa_AG" ,  0 , 0, 0, "\n# AG: Smoothing parameter 1 for non-ideal primal chemical potential increments (-1 to +1) { 1.0 }" },
   { "pa_DGC" , 0 , 0, 0, "\n# DGC: Smoothing parameter 2- exponent in smoothing function (-1 to +1) { 0 or 0.001 for adsorption }" },
   { "pa_PSM" , 0 , 0, 0, "\n# PSM: Level of diagnostic messages { 0- disabled (no ipmlog file); 1- default; 2-including warnings }" },
   { "pa_GAR" , 0 , 0, 0, "# GAR: Activity coefficient for major (M) species in solution phases at Simplex LP AIA { 1 }"  },
   { "pa_GAH" , 0 , 0, 0, "# GAH: Activity coefficient for minor (J) species in solution phases at Simplex LP AIA { 1000 }" },
   { "pa_DS",   0 , 0, 0, "\n# DS: Cutoff minimum amount of stable phase in GEM IPM primal solution, moles { 1e-20 }" },
   { "pa_XwMin" , 0 , 0, 0, "# XwMin: Cutoff mole amount of water-solvent for aqueous phase elimination { 1e-13 }" },
   { "pa_ScMin" , 0 , 0, 0, "# ScMin: Cutoff mole amount of solid sorbent for sorption phase elimination { 1e-13 }" },
   { "pa_DcMin" , 0 , 0, 0, "# DcMin: Cutoff mole amount for elimination of DC (species) in multi-component phase { 1e-33 }" },
   { "pa_PhMin" , 0 , 0, 0, "# PhMin: Cutoff mole amount for elimination of solution phases other than aqueous { 1e-20 }" },
   { "pa_ICmin" , 0 , 0, 0, "# ICmin: Cutoff effective molal ionic strength for calculation of aqueous activity coefficients { 1e-5 }" },
   { "pa_PC" ,   0 , 0, 0,  "\n# PC: Mode of Phase Selection: 1 old (Select-2), 2 new (PSSC), default { 2 }" },
   { "pa_DFM" ,  0 , 0, 0,  "# DFM: Tolerance for stability criterion for a phase to be eliminated from mass balance { 0.01 } " },
   { "pa_DFYw" , 0 , 0, 0,  "# DFYw: Insertion mole amount for water-solvent at Simplex()->MBR() bridge { 1e-5 }"},
   { "pa_DFYaq", 0 , 0, 0,  "# DFYaq: Insertion mole amount for aqueous species at Simplex()->MBR() bridge { 1e-5 }"  },
   { "pa_DFYid", 0 , 0, 0,  "\n# DFYid: Insertion mole amount for DCs of ideal solution phases at Simplex()->MBR() bridge { 1e-5 }" },
   { "pa_DFYr" , 0 , 0, 0,  "# DFYr: Insertion mole amount for major DCs in solution phases at Simplex()->MBR()bridge { 1e-5 }" },
   { "pa_DFYh" , 0 , 0, 0,  "# DFYh: Insertion mole amount for junior DCs in solution phases Simplex()->MBR() bridge{ 1e-5 }" },
   { "pa_DFYc" , 0 , 0, 0,  "# DFYc: Insertion mole amount for single-component phase at Simplex()->MBR() bridge { 1e-5 }" },
   { "pa_DFYs",  0 , 0, 0,  "# DFYs: Insertion mole amount for single-component phase in PSSC() algorithm { 1e-6 }" },
   { "pa_DW",    0 , 0, 0,  "# DW: Activate (1) or disable (0) error condition on maximum number of MBR() iterations DP { 1 }" },
   { "pa_DT",    0 , 0, 0,  "# DT: use DHB as relative maximum mass balance cutoff for all ICs (0, default); or for major ICs:"
                            "\n# decimal exponent (<-6) applied to DHB cutoff; (1) use DHB also as an absolute cutoff { 1 }" },
   { "pa_GAS",   0 , 0, 0,  "\n# GAS: Threshold for primal-dual chemical potential difference used in SpeciationCleanup() { 0.001 }" },
   { "pa_DG",    0 , 0, 0,  "# Total number of moles used in internal re-scaling of the system (disabled if < 1e-4) { 1000 }" },
   { "pa_DNS" ,  0 , 0, 0,  "# DNS: Standard surface number density, nm-2 for calculating activity of surface species { 12.05 }" },
   { "pa_IEPS" , 0 , 0, 0,  "# IEPS: Tolerance for calculation of surface activity coefficient terms for surface species { 0.001 }" },
   { "pKin" ,    0 , 0, 0,  "\n# pKin: Flag for using metastability constraints on DC amounts in primal GEM solution { 1 } " },
   { "pa_DKIN" , 0 , 0, 0,  "# DKIN: Tolerance for non-trivial metastability constraints on DC amounts, moles { 1e-10 } " },
   { "mui" ,     0 , 0, 0,  "\n# mui: IC indices in parent RMULTS IC list (not used in standalone GEMS3K)" },
   { "muk" ,     0 , 0, 0,  "\n# muk: Phase indices in parent RMULTS Phase list (not used in standalone GEMS3K)" },
   { "muj" ,     0 , 0, 0,  "\n# muj: DC indices in parent RMULTS DC list (not used in standalone GEMS3K)" },
   { "pa_PLLG" , 0 , 0, 0,  "# pa_PLLG: Tolerance for checking divergence in IPM dual solution, 1 to 32001 { 30000 }, 0 disables" },
   { "tMin" ,    0 , 0, 0,  "# tMin: Type of thermodynamic potential to minimize (reserved)" },
   { "dcMod",    0 , 0, 0,  "# dcMod: Codes for PT corrections of DC thermodynamic data [nDC] (reserved)" }
};



//===================================================================

/// Writing structure MULTI (GEM IPM work structure)
void TMulti::to_text_file_gemipm( const char *path, bool addMui,
		bool with_comments, bool brief_mode )
{
  SPP_SETTING *pa = paTProfil;
   bool _comment = with_comments;
   char PAalp;
   char PSigm;

#ifndef IPMGEMPLUGIN
   PAalp = TSyst::sm->GetSY()->PAalp;
   PSigm = TSyst::sm->GetSY()->PSigm;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
#endif
  fstream ff( path, ios::out );
  ErrorIf( !ff.good() , path, "Fileopen error");
  TPrintArrays  prar1( 8, MULTI_static_fields, ff);
  TPrintArrays  prar( 70, MULTI_dynamic_fields, ff);

// set up array flags for permanent fields
   if( !( pm.FIs > 0 && pm.Ls > 0 ) )
   {
     prar.setNoAlws( (long int)(f_sMod ) );
     prar.setNoAlws( f_LsMod );
     prar.setNoAlws( f_LsMdc );
   }
   if( PSigm == S_OFF )
   {
     prar.setNoAlws( f_Sigw);
     prar.setNoAlws( f_Sigg);
   }
   if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
   { /* ADSORPTION AND ION EXCHANGE */
     prar.setNoAlws( f_Nfsp);
     prar.setNoAlws( f_MASDT);
     prar.setNoAlws( f_C1 );
     prar.setNoAlws( f_C2 );
     prar.setNoAlws( f_C3 );
     prar.setNoAlws( f_pCh );
     prar.setNoAlws( f_SATX );
     prar.setNoAlws( f_MASDJ );
     prar.setNoAlws( f_SCM );
     prar.setNoAlws( f_SACT );
     prar.setNoAlws( f_DCads );
   }

if( _comment )
{  ff << "# " << _GEMIPM_version_stamp << endl << "# File: " << path << endl;
   ff << "# Comments can be marked with # $ ; as the first character in the line" << endl;
   ff << "# IPM text input file for the internal GEM IPM-3 kernel data" << endl;
   ff << "# (should be read after the DCH file and before DBR files)" << endl << endl;
   ff << "# ID key of the initial chemical system definition" << endl;
}
  ff << "\"" << pm.stkey << "\"" << endl;

 if( _comment )
     ff << "\n## (1) Flags that affect memory allocation";

 if(!brief_mode || pa->p.PE != pa_.p.PE )
   prar1.writeField(f_pa_PE, pa->p.PE, _comment, false  );

 //   ff << "# Do not know if this stuff is really necessary" << endl;
 //   ff << "# 'GWAT'         55.50837344" << endl;
 //   ff << left << setw(12) << "<GWAT> " <<  right << setw(8) << pm.GWAT << endl;

 prar1.writeField(f_PV, pm.PV, _comment, brief_mode  );
 prar1.writeField(f_PSOL, pm.PSOL, _comment, brief_mode  );
 if( _comment )
   ff << "\n\n# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases ";
 ff << endl << left << setw(12) << "<PAalp> " <<  right << setw(6) <<
    "\'" << PAalp << "\'" << endl;
 if( _comment )
  ff << "\n# PSigm: Flag for using (+) or ignoring (-) specific surface free energies  " << endl;
 ff << left << setw(12) << "<PSigm> " <<  right << setw(6) <<
    "\'" << PSigm << "\'" << endl;

  if( !brief_mode || pm.FIat > 0 || pm.Lads > 0 )
  { if( _comment )
      ff << "\n## (2) Dimensionalities that affect memory allocation";
    prar1.writeField(f_Lads, pm.Lads, _comment, false  );
    prar1.writeField(f_FIa, pm.FIa, _comment, false  );
    prar1.writeField(f_FIat,  pm.FIat, _comment, false  );

//   ff << left << setw(12) << "<sitNc> " <<  right << setw(8) << pm.sitNcat << endl;
//   ff << left << setw(12) << "<sitNa> " <<  right << setw(8) << pm.sitNan << endl;
    } // brief_mode

ff << "\n\n<END_DIM>\n";

// static data not affected by dimensionalities
   if( _comment )
   {  ff << "\n## (3) Numerical controls and tolerances of GEM IPM-3 kernel" << endl;
      ff << "#      - Need to be changed only in special cases (see gems3k_ipm.html)";
   }
   if( !brief_mode ||pa->p.DB != pa_.p.DB )
      prar.writeField(f_pa_DB, pa->p.DB, _comment, false  );
   if( !brief_mode ||pa->p.DHB != pa_.p.DHB )
      prar.writeField(f_pa_DHB, pa->p.DHB, _comment, false  );
   if( !brief_mode ||pa->p.EPS != pa_.p.EPS )
       prar.writeField(f_pa_EPS, pa->p.EPS, _comment, false  );
   if( !brief_mode ||pa->p.DK != pa_.p.DK )
       prar.writeField(f_pa_DK, pa->p.DK, _comment, false  );
   if( !brief_mode ||pa->p.DS != pa_.p.DS )
       prar.writeField(f_pa_DS,  pa->p.DS, _comment, false  );
   if( !brief_mode ||pa->p.DF != pa_.p.DF )
       prar.writeField(f_pa_DF, pa->p.DF, _comment, false  );
   if( !brief_mode ||pa->p.DFM != pa_.p.DFM )
       prar.writeField(f_pa_DFM,  pa->p.DFM, _comment, false  );
   if(!brief_mode || pa->p.DP != pa_.p.DP )
       prar.writeField(f_pa_DP,  pa->p.DP, _comment, false  );
   if(!brief_mode || pa->p.IIM != pa_.p.IIM )
       prar.writeField(f_pa_IIM,  pa->p.IIM, _comment, false  );
   if(!brief_mode || pa->p.PD != pa_.p.PD )
       prar.writeField(f_pa_PD,  pa->p.PD, _comment, false  );
   if(!brief_mode || pa->p.PRD != pa_.p.PRD )
       prar.writeField(f_pa_PRD,  pa->p.PRD, _comment, false  );
   if(!brief_mode || pa->p.AG != pa_.p.AG )
       prar.writeField(f_pa_AG,  pa->p.AG, _comment, false  );
   if(!brief_mode || pa->p.DGC != pa_.p.DGC )
       prar.writeField(f_pa_DGC,  pa->p.DGC, _comment, false  );
   if(!brief_mode || pa->p.PSM != pa_.p.PSM )
       prar.writeField(f_pa_PSM,  pa->p.PSM, _comment, false  );
   if(!brief_mode || pa->p.GAR != pa_.p.GAR )
       prar.writeField(f_pa_GAR,  pa->p.GAR, _comment, false  );
   if(!brief_mode || pa->p.GAH != pa_.p.GAH )
       prar.writeField(f_pa_GAH,  pa->p.GAH, _comment, false  );

   if(!brief_mode)
    if( _comment )
     {  ff << "\n\n# X*Min: Cutoff amounts for elimination of unstable species ans phases from mass balance";
     }

   if(!brief_mode || pa->p.XwMin != pa_.p.XwMin )
       prar.writeField(f_pa_XwMin,  pa->p.XwMin, _comment, false  );
   if(!brief_mode || pa->p.ScMin != pa_.p.ScMin )
       prar.writeField(f_pa_ScMin,  pa->p.ScMin, _comment, false  );
   if(!brief_mode || pa->p.DcMin != pa_.p.DcMin )
       prar.writeField(f_pa_DcMin,  pa->p.DcMin, _comment, false  );
   if(!brief_mode || pa->p.PhMin != pa_.p.PhMin )
       prar.writeField(f_pa_PhMin,  pa->p.PhMin, _comment, false  );
   if(!brief_mode || pa->p.ICmin != pa_.p.ICmin )
       prar.writeField(f_pa_ICmin,  pa->p.ICmin, _comment, false  );
   if(!brief_mode || pa->p.PC != pa_.p.PC )
       prar.writeField(f_pa_PC,  pa->p.PC, _comment, false  );

   if( _comment )
      ff << "\n# DFY: Insertion mole amounts used after the LPP AIA and in PhaseSelection() algorithm" << endl;

   if(!brief_mode || pa->p.DFYw != pa_.p.DFYw )
       prar.writeField(f_pa_DFYw,  pa->p.DFYw, _comment, false  );
   if(!brief_mode || pa->p.DFYaq != pa_.p.DFYaq )
       prar.writeField(f_pa_DFYaq,  pa->p.DFYaq, _comment, false  );
   if(!brief_mode || pa->p.DFYid != pa_.p.DFYid )
       prar.writeField(f_pa_DFYid,  pa->p.DFYid, _comment, false  );
   if(!brief_mode || pa->p.DFYr != pa_.p.DFYr )
       prar.writeField(f_pa_DFYr,  pa->p.DFYr, _comment, false  );
   if(!brief_mode || pa->p.DFYh != pa_.p.DFYh )
       prar.writeField(f_pa_DFYh,  pa->p.DFYh, _comment, false  );
   if(!brief_mode || pa->p.DFYc != pa_.p.DFYc )
       prar.writeField(f_pa_DFYc,  pa->p.DFYc, _comment, false  );
   if(!brief_mode || pa->p.DFYs != pa_.p.DFYs )
       prar.writeField(f_pa_DFYs,  pa->p.DFYs, _comment, false  );

   if( _comment )
     ff << "\n# Tolerances and controls of the high-precision IPM-3 algorithm ";

   if(!brief_mode || pa->p.DW != pa_.p.DW )
       prar.writeField(f_pa_DW,  pa->p.DW, _comment, false  );
   if(!brief_mode || pa->p.DT != pa_.p.DT )
       prar.writeField(f_pa_DT,  pa->p.DT, _comment, false  );
   if(!brief_mode || pa->p.GAS != pa_.p.GAS )
       prar.writeField(f_pa_GAS,  pa->p.GAS, _comment, false  );
   if(!brief_mode || pa->p.DG != pa_.p.DG )
       prar.writeField(f_pa_DG,  pa->p.DG, _comment, false  );
   if(!brief_mode || pa->p.DNS != pa_.p.DNS )
       prar.writeField(f_pa_DNS, pa->p.DNS, _comment, false  );
   if(!brief_mode || pa->p.IEPS != pa_.p.IEPS )
       prar.writeField(f_pa_IEPS, pa->p.IEPS, _comment, false  );
  prar.writeField(f_pKin, pm.PLIM, _comment, brief_mode  );
  if(!brief_mode || pa->p.DKIN != pa_.p.DKIN )
       prar.writeField(f_pa_DKIN, pa->p.DKIN, _comment, false  );
  if(!brief_mode || pa->p.PLLG != pa_.p.PLLG )
       prar.writeField(f_pa_PLLG, pa->p.PLLG, _comment, false  );
  if(!brief_mode || pm.tMin != G_TP_ )
       prar.writeField(f_tMin, pm.tMin, _comment, false  );

//dynamic arrays
if( pm.FIs > 0 && pm.Ls > 0 )
{
  if( _comment )
     ff << "\n\n## (4) Initial data for multicomponent phases (see DCH file for dimension nPHs)";
  prar.writeArrayF(  f_sMod, pm.sMod[0], pm.FIs, 6L, _comment, brief_mode );

long int LsModSum;
long int LsIPxSum;
long int LsMdcSum;
long int LsMsnSum;
long int LsSitSum;
getLsModsum( LsModSum, LsIPxSum );
getLsMdcsum( LsMdcSum, LsMsnSum, LsSitSum );

   prar.writeArray(  f_LsMod, pm.LsMod, pm.FIs*3, 3L, _comment, brief_mode);

  if(LsIPxSum )
  {
   if( _comment )
      ff << "\n\n# IPxPH: Index lists (in TSolMod convention) for interaction parameters of non-ideal solutions";
   prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
  }
  if(LsModSum )
   {
     if( _comment )
        ff << "\n\n# PMc: Tables (in TSolMod convention) of interaction parameter coefficients  for non-ideal solutions";
    prar.writeArray(  "PMc", pm.PMc,  LsModSum);
   }
   prar.writeArray(  f_LsMdc, pm.LsMdc, pm.FIs*3, 3L, _comment, brief_mode);
   if(LsMdcSum )
   {   if( _comment )
          ff << "\n\n# DMc: Tables (in TSolMod convention) of  parameter coefficients for dependent components";
    prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);
   }
   if(LsMsnSum )
   {   if( _comment )
          ff << "\n\n# MoiSN:  end member moiety / site multiplicity number tables (in TSolMod convention) ";
    prar.writeArray(  "MoiSN", pm.MoiSN,  LsMsnSum);
   }
} // sMod
  if( _comment )
    ff << "\n\n## (5) Data arrays which are provided neither in DCH nor in DBR files";
  prar.writeArray(  f_B, pm.B,  pm.N, -1L, _comment, brief_mode);

  if( _comment )
     ff << "\n\n# Initial data for DCs - see DATACH file for dimensions nDC, nDCs";
  prar.writeArray(  f_Pparc, pm.Pparc,  pm.L, -1L, _comment, brief_mode);
  //  ff << "\n\n# This is not necessary - can be calculated from G0 ???????????";
  // prar.writeArray(  "G0", pm.G0,  pm.L);
  prar.writeArray(  f_fDQF, pm.fDQF,  pm.L, -1L, _comment, brief_mode);
  prar.writeArray(  f_lnGmf, pm.lnGmf,  pm.L, -1L, _comment, brief_mode);

  if( _comment )
     ff << "\n\n# (6) Metastability constraints on DC amounts from above (DUL) and below (DLL)";
   prar.writeArrayF(  f_RLC, pm.RLC, pm.L, 1L, _comment, brief_mode );
   prar.writeArrayF(  f_RSC, pm.RSC, pm.L, 1L, _comment, brief_mode );
   prar.writeArray(  f_DLL, pm.DLL, pm.L, -1L, _comment, brief_mode);
   prar.writeArray(  f_DUL, pm.DUL,  pm.L, -1L, _comment, brief_mode);

   if( _comment )
     ff << "\n\n# (7) Initial data for Phases" << endl;
   prar.writeArray(  f_Aalp, pm.Aalp,  pm.FI, -1L, _comment, brief_mode);
   if( PSigm != S_OFF )
   {
     prar.writeArray(  f_Sigw, pm.Sigw,  pm.FI, -1L, _comment, brief_mode);
     prar.writeArray(  f_Sigg, pm.Sigg,  pm.FI, -1L, _comment, brief_mode);
   }
   prar.writeArray(  f_YOF, pm.YOF,  pm.FI, -1L, _comment, brief_mode);

   if( pm.FIat > 0 &&  pm.FIs > 0 )
    { // ADSORPTION AND ION EXCHANGE
      if( _comment )
        ff << "\n\n# (8) Initial data for sorption phases";

      prar.writeArray(  f_Nfsp, &pm.Nfsp[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_MASDT, &pm.MASDT[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_C1, &pm.XcapA[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_C2, &pm.XcapB[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_C3, &pm.XcapF[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_pCh, &pm.Xetaf[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
      prar.writeArray(  f_SATX, &pm.SATX[0][0], pm.Lads*4, 4L, _comment, brief_mode);
      prar.writeArray(  f_MASDJ, &pm.MASDJ[0][0], pm.Lads*DFCN, (long int)DFCN, _comment, brief_mode);
      prar.writeArrayF(  f_SCM, pm.SCM[0], pm.FIs, pm.FIat, _comment, brief_mode );
      prar.writeArrayF(  f_SACT, pm.SATT, pm.Lads, 1L, _comment, brief_mode );
      prar.writeArrayF(  f_DCads, pm.DCC3, pm.Lads, 1L, _comment, brief_mode );
    }

   //if(!brief_mode || prar.getAlws("dcMod" ))
   prar.writeArrayF(  f_dcMod, pm.dcMod[0], pm.L, 6L, _comment, brief_mode );

 /*
   outArray( ff, "Vol", pm.Vol,  pm.L);
   outArray( ff, "G0", pm.G0,  pm.L);
   outArray( ff, "PUL", pm.PUL,  pm.L);
   outArray( ff, "PLL", pm.PLL,  pm.L);
   outArray( ff, "lnGam", pm.lnGam,  pm.L);
   outArray( ff, "F0", pm.F0,  pm.L);
*/

 if( addMui && !brief_mode )
 {
   prar.writeArray(  f_mui, pm.mui,  pm.N, -1L, _comment, brief_mode);
   prar.writeArray(  f_muk, pm.muk,  pm.FI, -1L, _comment, brief_mode);
   prar.writeArray(  f_muj, pm.muj,  pm.L, -1L, _comment, brief_mode);
 }

 ff << endl;
 if( _comment )
   ff << "\n# End of file" << endl;

}

/// Reading structure MULTI (GEM IPM work structure)
void TMulti::from_text_file_gemipm( const char *path,  DATACH  *dCH )
{
  SPP_SETTING *pa = paTProfil;
  long int ii, nfild, len;

   //static values
   char PAalp;
   char PSigm;

#ifdef IPMGEMPLUGIN
   set_def();
#endif
  //mem_set( &pm.N, 0, 38*sizeof(long int));
  //mem_set( &pm.TC, 0, 55*sizeof(double));
  // get sizes from DATACH
  pm.TC = pm.TCc = 25.;
  pm.T = pm.Tc =298.15;
  pm.P = pm.Pc = 1.;
  pm.N = pm.NR = dCH->nIC;
  pm.L = dCH->nDC;
  pm.FI = dCH->nPH;
  pm.FIs = dCH->nPS;
  //
  pm.Ls = 0; //dCH->nDCs;
  for( ii=0; ii<dCH->nPS; ii++)
  {
    pm.Ls += dCH->nDCinPH[ii];
    if( dCH->ccPH[ii] == 'a' )
     pm.LO = pm.Ls-1;
    if( dCH->ccPH[ii] == 'g' || dCH->ccPH[ii] == 'p' || dCH->ccPH[ii] == 'f')
      pm.PG = dCH->nDCinPH[ii];
  }

  // setup default constants
  pa->p.PE =  pm.E = 1;
  pm.PV = 0;
  pm.PSOL = 0;
  PAalp = '+';
  PSigm = '+';
  pm.Lads = 0;
  pm.FIa = 0;
  pm.FIat = 0; //6
  pm.PLIM  = 1;

  // reads sizes and constants from txt file
  fstream ff( path, ios::in );
  ErrorIf( !ff.good() , path, "Fileopen error");

// static data
   TReadArrays  rdar( 8, MULTI_static_fields, ff);
   gstring str;
   rdar.skipSpace();
   f_getline( ff, str, '\n');
   copyValues( pm.stkey, (char * )str.c_str(), EQ_RKLEN );

   nfild = rdar.findNext();
   while( nfild >=0 )
   {
     switch( nfild )
     {
       case f_pa_PE: rdar.readArray("pa_PE" , &pa->p.PE, 1);
                 pm.E = pa->p.PE;
              break;
       case f_PV: rdar.readArray("PV" , &pm.PV, 1);
              break;
       case f_PSOL: rdar.readArray("PSOL" , &pm.PSOL, 1);
              break;
       case f_PAalp: rdar.readArray("PAalp" , &PAalp, 1, 1);
              break;
       case f_PSigm: rdar.readArray("PSigm" , &PSigm, 1, 1);
              break;
       case f_Lads: rdar.readArray("Lads" , &pm.Lads, 1);
              break;
       case f_FIa: rdar.readArray("FIa" , &pm.FIa, 1);
              break;
       case f_FIat: rdar.readArray("FIat" , &pm.FIat, 1);
              break;
    }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from the MULTI structure";
    Error( "Error", ret);
  }

#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#else
   dyn_new();
#endif

// get dynamic data from DATACH file
  for( ii=0; ii<dCH->nPH; ii++)
    pm.L1[ii] = dCH->nDCinPH[ii];

  for( ii=0; ii<dCH->nIC*dCH->nDC; ii++)
    pm.A[ii] = dCH->A[ii];

  if( pm.EZ )
  { long int iZ=-1;
    for(  ii=0; ii<dCH->nDC; ii++ )
     if( dCH->ccIC[ii] == IC_CHARGE )
         break;
    if( ii< dCH->nDC )
    { iZ = ii;
      for( ii=0; ii<dCH->nDC; ii++)
          pm.EZ[ii] = pm.A[pm.N*ii+iZ];
    }
  }

  for( ii=0; ii< dCH->nIC; ii++ )
  { pm.Awt[ii]  = dCH->ICmm[ii]*1e3;
    fillValue(pm.SB[ii], ' ', MaxICN );
    len = strlen(dCH->ICNL[ii]);
    //len = min(  len,MaxICN);
    copyValues( pm.SB[ii], dCH->ICNL[ii], min(len,(long int)MAXICNAME));
    pm.SB[ii][MaxICN] = dCH->ccIC[ii];
    pm.ICC[ii] =  dCH->ccIC[ii];
  }

if( fabs(dCH->DCmm[0]) < 1e-32 )  // Restore DCmm if skipped from the DCH file
  for( long int jj=0; jj< dCH->nDC; jj++ )  // Added by DK on 03.03.2007
  {
    dCH->DCmm[jj] = 0.0;
    for( ii=0; ii< dCH->nIC; ii++ )
       dCH->DCmm[jj] += dCH->ICmm[ii]*dCH->A[jj*dCH->nIC+ii];
  }

  for( ii=0; ii< dCH->nDC; ii++ )
  {
    pm.MM[ii] = dCH->DCmm[ii]*1e3;
    pm.DCC[ii] = dCH->ccDC[ii];
    len =strlen(dCH->DCNL[ii]);
    //len = min(  len,MaxDCN);
    copyValues( pm.SM[ii], dCH->DCNL[ii], min(len,(long int)MAXDCNAME) );
  }

  for( ii=0; ii< dCH->nPH; ii++ )
  {
     len =strlen(dCH->PHNL[ii]);
     //len = min(  len,MaxPHN);
          fillValue( pm.SF[ii], ' ', MAXPHNAME+MAXSYMB );
          copyValues( pm.SF[ii]+MAXSYMB, dCH->PHNL[ii], min(len,(long int)MAXPHNAME) );
     pm.SF[ii][0] = dCH->ccPH[ii];
     pm.PHC[ii] = dCH->ccPH[ii];
  }

// !!!!  copyValues( pm.DCCW, dCH->ccDCW, dCH->nDC);
  // set up DCCW
  ConvertDCC();

//reads dynamic values from txt file
   TReadArrays  rddar( 70, MULTI_dynamic_fields, ff);

// set up array flags for permanent fields

 if( !( pm.FIs > 0 && pm.Ls > 0 ) )
 {
   rddar.setNoAlws( (long int)(f_sMod ));
   rddar.setNoAlws( f_LsMod );
   rddar.setNoAlws( f_LsMdc );
 }
 if( PSigm == S_OFF )
 {
   rddar.setNoAlws( f_Sigw );
   rddar.setNoAlws( f_Sigg );
 }
 if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
 { /* ADSORPTION AND ION EXCHANGE */
   rddar.setNoAlws( f_Nfsp );
   rddar.setNoAlws( f_MASDT );
   rddar.setNoAlws( f_C1 );
   rddar.setNoAlws( f_C2 );
   rddar.setNoAlws( f_C3 );
   rddar.setNoAlws( f_pCh );
   rddar.setNoAlws( f_SATX );
   rddar.setNoAlws( f_MASDJ );
   rddar.setNoAlws( f_SCM );
   rddar.setNoAlws( f_SACT );
   rddar.setNoAlws( f_DCads );
 }

  // Read dynamic arrays
  nfild = rddar.findNext();
  while( nfild >=0 )
  {
    switch( nfild )
    { case f_sMod: if( !pm.sMod )
                Error( "Error", "Array sMod is not used in this problem");
              rddar.readArray( "sMod" , pm.sMod[0], pm.FIs, 6 );
              break;
      case f_LsMod:{ if( !pm.LsMod )
                Error( "Error", "Array LsMod is not used in this problem");
              rddar.readArray( "LsMod" , pm.LsMod, pm.FIs*3) ;
              long int LsModSum;
              long int LsIPxSum;
              getLsModsum( LsModSum, LsIPxSum );
              if(LsIPxSum )
              { rddar.readNext( "IPxPH");
#ifdef IPMGEMPLUGIN
              if(!pm.IPx )
                  pm.IPx = new long int[LsIPxSum];
#else
                 pm.IPx = (long int *)aObj[ o_wi_ipxpm ].Alloc(LsIPxSum, 1, L_);
#endif
                rddar.readArray( "IPxPH", pm.IPx,  LsIPxSum);
              }
              if(LsModSum )
              { rddar.readNext( "PMc");
#ifdef IPMGEMPLUGIN
              if(!pm.PMc )
                  pm.PMc = new double[LsModSum];
#else
               pm.PMc = (double *)aObj[ o_wi_pmc].Alloc( LsModSum, 1, D_);
#endif
                rddar.readArray( "PMc", pm.PMc,  LsModSum);
              }
              break;
             }
      case f_LsMdc: { if( !pm.LsMdc )
                   Error( "Error", "Array LsMdc not used in this problem");
                rddar.readArray( "LsMdc" , pm.LsMdc, pm.FIs*3 );
                long int LsMdcSum;
                long int LsMsnSum;
                long int LsSitSum;
                getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
                if(LsMdcSum )
                { rddar.readNext( "DMc");
#ifdef IPMGEMPLUGIN
                if(!pm.DMc )
                     pm.DMc = new double[LsMdcSum];
#else
                pm.DMc = (double *)aObj[ o_wi_dmc].Alloc( LsMdcSum, 1, D_ );
#endif
                rddar.readArray( "DMc", pm.DMc,  LsMdcSum);
                }
              if(LsMsnSum )
              { rddar.readNext( "MoiSN");
#ifdef IPMGEMPLUGIN
              if(!pm.MoiSN )
                   pm.MoiSN = new double[LsMsnSum];
              pm.SitFr = new double[LsSitSum];
#else
              pm.MoiSN = (double *)aObj[ o_wi_moisn].Alloc( LsMsnSum, 1, D_ );
              pm.SitFr  = (double *)aObj[ o_wo_sitfr ].Alloc( LsSitSum, 1, D_ );
#endif
              fillValue( pm.SitFr, 0., LsSitSum );
              rddar.readArray( "MoiSN", pm.MoiSN,  LsMsnSum);
              }
              break;
            }
      case f_B: rddar.readArray( "B", pm.B,  pm.N);
              break;
      case f_DCCW: rddar.readArray( "DCCW", pm.DCCW,  pm.L, 1);
              break;
      case f_Pparc: rddar.readArray( "Pparc", pm.Pparc,  pm.L);
              break;
      case f_fDQF: rddar.readArray( "fDQF", pm.fDQF,  pm.L);
              break;
      case f_lnGmf: rddar.readArray( "lnGmf", pm.lnGmf,  pm.L);
              break;
      case f_RLC: rddar.readArray( "RLC", pm.RLC, pm.L, 1 );
              break;
      case f_RSC: rddar.readArray( "RSC", pm.RSC, pm.L, 1 );
              break;
      case f_DLL: rddar.readArray( "DLL", pm.DLL,  pm.L);
              break;
      case f_DUL: rddar.readArray( "DUL", pm.DUL,  pm.L);
              break;
      case f_Aalp: rddar.readArray( "Aalp", pm.Aalp,  pm.FI);
              break;
      case f_Sigw: if( !pm.Sigw )
                Error( "Error", "Array Sigw not used in this problem");
              rddar.readArray( "Sigw", pm.Sigw,  pm.FI);
              break;
      case f_Sigg: if( !pm.Sigg )
                Error( "Error", "Array Sigg not used in this problem");
              rddar.readArray( "Sigg", pm.Sigg,  pm.FI);
              break;
      case f_YOF: rddar.readArray( "YOF", pm.YOF,  pm.FI);
              break;
      case f_Nfsp: if( !pm.Nfsp )
                Error( "Error", "Array Nfsp not used in this problem");
              rddar.readArray( "Nfsp", &pm.Nfsp[0][0], pm.FIs*pm.FIat);
              break;
      case f_MASDT: if( !pm.MASDT )
                Error( "Error", "Array MASDT not used in this problem");
              rddar.readArray( "MASDT", &pm.MASDT[0][0], pm.FIs*pm.FIat);
              break;
      case f_C1: if( !pm.XcapA )
                Error( "Error", "Array XcapA not used in this problem");
              rddar.readArray( "C1", &pm.XcapA[0][0], pm.FIs*pm.FIat);
              break;
      case f_C2: if( !pm.XcapB )
                Error( "Error", "Array XcapB not used in this problem");
              rddar.readArray( "C2", &pm.XcapB[0][0], pm.FIs*pm.FIat);
              break;
      case f_C3: if( !pm.XcapF )
                Error( "Error", "Array XcapF not used in this problem");
              rddar.readArray( "C3", &pm.XcapF[0][0], pm.FIs*pm.FIat);
              break;
      case f_pCh: if( !pm.Xetaf )
                Error( "Error", "Array Xetaf not used in this problem");
              rddar.readArray( "pCh", &pm.Xetaf[0][0], pm.FIs*pm.FIat);
              break;
      case f_SATX: if( !pm.SATX )
                Error( "Error", "Array SATX not used in this problem");
              rddar.readArray( "SATX", &pm.SATX[0][0], pm.Lads*4);
              break;
      case f_MASDJ: if( !pm.MASDJ )
                Error( "Error", "Array MASDJ not used in this problem");
              rddar.readArray( "MASDJ", &pm.MASDJ[0][0], pm.Lads*DFCN);
              break;
      case f_SCM: if( !pm.SCM )
                Error( "Error", "Array SCM not used in this problem");
              rddar.readArray( "SCM", pm.SCM[0], pm.FIs, pm.FIat );
              break;
      case f_SACT: if( !pm.SATT )
                Error( "Error", "Array SATT not used in this problem");
              rddar.readArray( "SACT", pm.SATT, pm.Lads, 1 );
              break;
      case f_DCads: if( !pm.DCC3 )
                Error( "Error", "Array DCC3 not used in this problem");
               rddar.readArray( "DCads", pm.DCC3, pm.Lads, 1 );
               break;
      case f_pa_DB: rddar.readArray( "pa_DB" , &pa->p.DB, 1);
               break;
      case f_pa_DHB: rddar.readArray("pa_DHB", &pa->p.DHB, 1);
               break;
      case f_pa_EPS: rddar.readArray("pa_EPS" , &pa->p.EPS, 1);
               break;
      case f_pa_DK: rddar.readArray("pa_DK" , &pa->p.DK, 1);
               break;
      case f_pa_DF: rddar.readArray("pa_DF" , &pa->p.DF, 1);
               break;
      case f_pa_DP: rddar.readArray("pa_DP", &pa->p.DP, 1);
               break;
      case f_pa_IIM: rddar.readArray("pa_IIM", &pa->p.IIM, 1);
               break;
      case f_pa_PD: rddar.readArray("pa_PD" , &pa->p.PD, 1);
               break;
      case f_pa_PRD: rddar.readArray("pa_PRD" , &pa->p.PRD, 1);
               break;
      case f_pa_AG: rddar.readArray("pa_AG" , &pa->p.AG, 1);
               break;
      case f_pa_DGC: rddar.readArray("pa_DGC" , &pa->p.DGC, 1);
               break;
      case f_pa_PSM: rddar.readArray("pa_PSM" , &pa->p.PSM, 1);
               break;
      case f_pa_GAR: rddar.readArray("pa_GAR" , &pa->p.GAR, 1);
               break;
      case f_pa_GAH: rddar.readArray("pa_GAH" , &pa->p.GAH, 1);
               break;
      case f_pa_DS: rddar.readArray("pa_DS", &pa->p.DS, 1);
               break;
      case f_pa_XwMin: rddar.readArray("pa_XwMin" , &pa->p.XwMin, 1);
               break;
      case f_pa_ScMin: rddar.readArray("pa_ScMin" , &pa->p.ScMin, 1);
               break;
      case f_pa_DcMin: rddar.readArray("pa_DcMin" , &pa->p.DcMin, 1);
               break;
      case f_pa_PhMin: rddar.readArray("pa_PhMin" , &pa->p.PhMin, 1);
               break;
      case f_pa_ICmin: rddar.readArray("pa_ICmin" , &pa->p.ICmin, 1);
               break;
      case f_pa_PC: rddar.readArray("pa_PC" , &pa->p.PC, 1);
               break;
      case f_pa_DFM: rddar.readArray("pa_DFM" , &pa->p.DFM, 1);
               break;
      case f_pa_DFYw: rddar.readArray("pa_DFYw" , &pa->p.DFYw, 1);
               break;
      case f_pa_DFYaq: rddar.readArray("pa_DFYaq" , &pa->p.DFYaq, 1);
               break;
      case f_pa_DFYid: rddar.readArray("pa_DFYid" , &pa->p.DFYid, 1);
               break;
      case f_pa_DFYr: rddar.readArray("pa_DFYr" , &pa->p.DFYr, 1);
               break;
      case f_pa_DFYh: rddar.readArray("pa_DFYh" , &pa->p.DFYh, 1);
               break;
      case f_pa_DFYc: rddar.readArray("pa_DFYc" , &pa->p.DFYc, 1);
               break;
      case f_pa_DFYs: rddar.readArray("pa_DFYs", &pa->p.DFYs, 1);
               break;
      case f_pa_DW: rddar.readArray("pa_DW", &pa->p.DW , 1);
               break;
      case f_pa_DT: rddar.readArray("pa_DT", &pa->p.DT , 1);
               break;
      case f_pa_GAS: rddar.readArray("pa_GAS", &pa->p.GAS, 1);
               break;
      case f_pa_DG: rddar.readArray("pa_DG" , &pa->p.DG, 1);
               break;
      case f_pa_DNS: rddar.readArray("pa_DNS" , &pa->p.DNS, 1);
               break;
      case f_pa_IEPS: rddar.readArray("pa_IEPS" , &pa->p.IEPS, 1);
               break;
      case f_pKin: rddar.readArray("pKin" , &pm.PLIM, 1);
               break;
      case f_pa_DKIN: rddar.readArray("pa_DKIN" , &pa->p.DKIN, 1);
               break;
      case f_mui: rddar.readArray("mui" , pm.mui, pm.N);
               break;
      case f_muk: rddar.readArray("muk" , pm.muk, pm.FI);
               break;
      case f_muj: rddar.readArray("muj" , pm.muj, pm.L);
               break;
      case f_pa_PLLG: rddar.readArray("pa_PLLG" , &pa->p.PLLG, 1);
               break;
      case f_tMin: rddar.readArray("tMin" , &pm.tMin, 1);
             break;
       case f_dcMod:   rddar.readArray( "dcMod" , pm.dcMod[0], pm.L, 6 );
                  break;
    }
    nfild = rddar.findNext();
  }
 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from the MULTY structure";
    Error( "Error", ret);
  }
}

//=============================================================================
// ms_multi_format.cpp

