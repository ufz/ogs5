/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
// using namespace std;

#ifndef OGS_FEM_CAP  // CAP_REACT
#ifdef UNIX
#define LI long
#define LIP long*
#define LNT long
#define DB double
#define DBP double*
#define CHP char*
#define CMT extern int
#define VDP void*
#define ftnlen long /* FORTRAN string length type */
#else
#define LI long            /* unsigned int		*/
#define LIP long*          /* unsigned int*	*/
#define LNT long           /* unsigned int		*/
#define DB double          /* double		*/
#define DBP double*        /* double*		*/
#define CHP char*          /* char*		*/
#define CMT void __stdcall /* void __stdcall	*/
#define VDP void*          /* void*		*/
#endif
#endif

extern int CAP_MODE, CAP_Time, CAP_Node, CAP_icount, CAP_count;
extern std::vector<std::vector<std::string> > CHEM_STATE, CHEM_STATE_AP;
extern std::vector<std::vector<int> > PHASE_CONSTI_MAP;

extern bool CAP_check_file(void);
extern int OGS_keyword_check(std::string in_file_new, std::string in_file_old,
                             std::string ext_name);

extern int CAP_tqini(LIP NOERR);
extern int CAP_tqopen(CHP FILE, LI LUN, LIP NOERR);
extern int CAP_tqclos(LI LUN, LIP NOERR);
extern int CAP_tqgio(CHP OPTION, LIP IVAL, LIP NOERR);
extern int CAP_tqcio(CHP OPTION, LI IVAL, LIP NOERR);
extern int CAP_tqrfil(LIP NOERR);
extern int CAP_tqgsu(CHP OPTION, CHP UNIT, LIP NOERR);
extern int CAP_tqcsu(CHP OPTION, CHP UNIT, LIP NOERR);
extern int CAP_tqinsc(CHP NAME, LIP INDEXS, LIP NOERR);
extern int CAP_tqgnsc(LI INDEXS, CHP NAME, LIP NOERR);
extern int CAP_tqnosc(LIP NSCOM, LIP NOERR);
extern int CAP_tqstsc(LI INDEXS, DBP STOI, DBP WMASS, LIP NOERR);
extern int CAP_tqcsc(CHP NAME, LIP NOERR);
extern int CAP_tqinp(CHP NAME, LIP INDEXP, LIP NOERR);
extern int CAP_tqgnp(LI INDEXP, CHP NAME, LIP NOERR);
extern int CAP_tqnop(LIP NPHASE, LIP NOERR);
extern int CAP_tqinpc(CHP NAME, LI INDEXP, LIP INDEXC, LIP NOERR);
extern int CAP_tqgnpc(LI INDEXP, LI INDEXC, CHP NAME, LIP NOERR);
extern int CAP_tqnopc(LI INDEXP, LIP NPCON, LIP NOERR);
extern int CAP_tqstpc(LI INDEXP, LI INDEXC, DBP STOI, DBP WMASS, LIP NOERR);
extern int CAP_tqgsp(LI INDEXP, CHP OPTION, LIP NOERR);
extern int CAP_tqcsp(LI INDEXP, CHP OPTION, LIP NOERR);
extern int CAP_tqgspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR);
extern int CAP_tqcspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR);
extern int CAP_tqsetc(CHP OPTION, LI INDEXP, LI INDEX, DB VAL, LIP NUMCON,
                      LIP NOERR);
extern int CAP_tqremc(LI NUMCON, LIP NOERR);
extern int CAP_tqsttp(CHP IDENTS, DBP VALS, LIP NOERR);
extern int CAP_tqstca(CHP IDENTS, LI INDEXP, LI INDEXC, DB VAL, LIP NOERR);
extern int CAP_tqstec(CHP OPTION, LI INDEXP, DB VAL, LIP NOERR);
extern int CAP_tqstrm(CHP IDENTS, LIP NOERR);
extern int CAP_tqce(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
extern int CAP_tqcel(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
extern int CAP_tqclim(CHP OPTION, DB VAL, LIP NOERR);
extern int CAP_tqgetr(CHP OPTION, LI INDEXP, LI INDEX, DBP VAL, LIP NOERR);
extern int CAP_tqgdpc(CHP OPTION, LI INDEXP, LI INDEXC, DBP VAL, LIP NOERR);
extern int CAP_tqshow(LIP NOERR);
extern int CAP_tqerr(CHP MESS, LIP NOERR);

extern int CAP_tqcprt(LIP NOERR);
extern int CAP_tqvers(LIP NVERS, LIP NOERR);
extern int CAP_tqsize(LIP NA, LIP NB, LIP NC, LIP ND, LIP NE, LIP NF, LIP NG,
                      LIP NH, LIP NI, LIP NJ, LIP NK, LIP NOERR);
extern int CAP_tqmodl(LI INDEXP, CHP NAME, LIP NOERR);
extern int CAP_tqstxp(CHP IDENTS, CHP OPTION, DBP VAL, LIP NOERR);
extern int CAP_tqlite(LIP LITE, LIP NOERR);
extern int CAP_tqrbin(LIP NOERR);
extern int CAP_tqmap(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT,
                     LIP NOERR);
extern int CAP_tqmapl(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT,
                      LIP NOERR);
extern int CAP_tqpcis(LI INDEXP, LI INDEXC, LIP ISPERM, LIP NOERR);
extern int CAP_tqopna(CHP FILE, LI LUN, LIP NOERR);
extern int CAP_tqopnb(CHP FILE, LI LUN, LIP NOERR);
extern int CAP_tqnosl(LI INDEXP, LIP NSUBL, LIP NOERR);
extern int CAP_tqnolc(LI INDEXP, LI INDEXL, LIP NSLCON, LIP NOERR);
extern int CAP_tqinlc(CHP NAME, LI INDEXP, LI INDEXL, LIP INDEXC, LIP NOERR);
extern int CAP_tqgnlc(LI INDEXP, LI INDEXL, LI INDEXC, CHP NAME, LIP NOERR);
extern int CAP_tqgtlc(LI INDEXP, LI INDEXL, LI INDEXC, DBP VAL, LIP NOERR);

/*
extern int CAP_tqgopn (CHP FILE,LI LUN,CHP FFORM,CHP FSTAT,CHP FACC,LI RECL,
        LIP IOSTAT,LIP NOERR);
*/
extern int CAP_tqbond(LI INDEXP, LI INDEXA, LI INDEXB, LI INDEXC, LI INDEXD,
                      DBP VAL, LIP NOERR);
extern int CAP_tqused(LIP NA, LIP NB, LIP NC, LIP ND, LIP NE, LIP NF, LIP NG,
                      LIP NH, LIP NI, LIP NJ, LIP NK, LIP NOERR);
extern int CAP_tqgtrh(LIP TFHVER,
                      CHP TFHNWP,
                      LIP TFHVNW,
                      CHP TFHNRP,
                      LIP TFHVNR,
                      LIP TFHDTC,
                      LIP TFHDTE,
                      CHP TFHID,
                      CHP TFHUSR,
                      CHP TFHREM,
                      LIP NOERR);
extern int CAP_tqopnt(CHP FILE, LI LUN, LIP NOERR);
extern int CAP_tqrcst(LIP NOERR);
extern int CAP_tqgtid(CHP ID, LIP NOERR);
extern int CAP_tqgtnm(CHP NAME, LIP NOERR);
extern int CAP_tqgtpi(CHP PID, LIP NOERR);
extern int CAP_tqwstr(CHP OPTION, CHP CTXT, LIP NOERR);
extern int CAP_tqgted(LIP EDMON, LIP EDYEAR, LIP NOERR);
extern int CAP_tqgthi(CHP HASPT, LIP HASPID, LIP NOERR);
extern int CAP_tqcen(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
extern int CAP_tqcenl(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
extern int CAP_tqwasc(CHP FILE, LIP NOERR);
extern int CAP_tqcdat(LI I1, LI I2, LI I3, LI I4, LI I5, DB VAL, LIP NOERR);
extern int CAP_tqchar(LI INDEXP, LI INDEXC, DBP VAL, LIP NOERR);
extern int CAP_tqcnsc(LI INDEXS, CHP NAME, LIP NOERR);
extern int CAP_tqlpar(LI INDEXP, CHP OPTION, LIP NOPAR, CHP CHRPAR, LIP LGTPAR,
                      LIP NOERR);
extern int CAP_tqgpar(LI INDEXP, CHP OPTION, LI INDEXX, LIP NOEXPR, LIP NVALA,
                      DBP VALA, LIP NOERR);
