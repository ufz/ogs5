/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: tools.h
 */
/* Aufgabe:
   verschiedene Funktionen, die von verschiedenen Modulen gebraucht
   werden und keine Adaptivitaet voraussetzen (sie aber teilweise
   unterstuetzen)
 */
/**************************************************************************/

#ifndef tools_INC
#define tools_INC

#include "rf_mmp_new.h"                           //MB
#include "rf_pcs.h"
#include <sstream>

typedef struct                                    /* fuer Kurven (Stuetzstellen) */
{
	double punkt;
	double wert;
} StuetzStellen;

typedef struct                                    /* fuer Kurven (Kurven) */
{
	long anz_stuetzstellen;
	/* Anzahl der Stuetzstellen */
	StuetzStellen* stuetzstellen;
	/* Feld mit den eingelesenen Stuetzstellen */
} Kurven;

//NB
//extern double GetMatrixValue (double, double, std::string, int*);
double GetMatrixValue(double var1, double var2, std::string caption, int *gueltig);
extern double GetCurveValue ( int kurve, int methode, double punkt, int* gueltig);
extern double GetCurveValueInverse ( int kurve, int methode, double wert, int* gueltig);
extern double GetCurveDerivative(int kurve, int methode, double punkt, int* gueltig);
extern double GetCurveInverseDerivative ( int kurve, int methode, double wert, int* gueltig);
extern Kurven* kurven;                            /* Feld mit Kurven */
extern int anz_kurven;                            /* Anzahl der Kurven */

/******************************************************/
/* C1.11 Miscellaneous                                 */
/******************************************************/

/* Deklarationen fuer Schluesselwort #FRACTURE_APERTURE_DISTRIBUTION */
extern double* fracture_aperture_array;
/* Feld mit Knoten, die ausgegeben werden sollen */
extern long fracture_aperture_anz;                /* Feldgroesse */

/* Baut Element-zu-Knoten-Verzeichnis auf (auch adaptiv) */
extern void ConstructElemsToNodesList ( void );

/* belegt Loesungsvektor vor --> NULLE_ERGEBNIS */
extern void PresetErgebnis ( double* ergebnis, int nidx );

/* berechnet Anfangszeitschritt --> automat. Zeitschrittsteuerung (aTM) */
extern double StartTimeStep ( int dtidx );
/* Ermittelt den zulaessigen Courant-Zeitschritt pro Element */
extern double CalcCourantTimeStep ( long index, long ndx, double dc);

/* Ermittelt den zulaessigen Courant-Zeitschritt im System */
extern double CalcSystemCourantTimeStep ( long ndx, double dc);

/* Sucht den zulaessigen Courant/Neumann-Zeitschritt im System */
extern double GetSystemCourantNeumannTimeStep ( long ndx, int dtidx, double acknowledge );

/* Testet, ob sich Knotenwerte im Element geaendert haben */
extern int TestElementDirtyness ( long index, long ndx1, long ndx2, double acknowledge);
/* Ermittelt das Vorzeichen */
extern int Signum(double);
/* Bildet den arithmetischen Mittel einer Elementgroesse durch Interpolation
   der zugehoerigen Knotenwerte */
double InterpolateElementNodesValues ( long index, long idx );
extern int FctCurves ( char* data, int found, FILE* f );
//SB
//extern int FctReadHeterogeneousPermeabilityField(char* name_file);
//SB
long GetNearestElement(double x,double y,double z, int* help);
//SB
long GetNearestHetVal(long EleIndex, CFEMesh*, long no_values, double** invals);
//MB
double GetAverageHetVal(long EleIndex, CFEMesh*, long no_values, double** invals);
extern double GetHetValue(int,char*);             //SB

typedef struct
{
	int nof;                              //number of field values
	char* filename;                       // Name of input file for heterogeneous values
	char** names;                         //Names of field variables
	double* convertfact;                  //conversion factors
}hetfields;                                       //SB

extern hetfields* hf;                             //SB
                                                  //SB
extern hetfields* Createhetfields(int n, char* name_file);
//SB
extern void set_hetfields_name(hetfields* hf,int i, char* name);
//SB
extern char* get_hetfields_name(hetfields* hf,int i);
extern int get_hetfields_number(hetfields* hf);   //SB
                                                  //MB
extern int FctReadHeterogeneousFields(char* name_file, CMediumProperties*);
extern long DampOscillations(int ndx1,
                             int oscil_damp_method,
                             double* oscil_damp_parameter,
                             double (* NodeCalcLumpedMass)(long));
extern int GetLineFromFile(char*, std::ifstream*);

typedef struct
{
	long row;
	long col;
	double** m;
} DMATRIX;

extern DMATRIX* CreateDoubleMatrix(long row, long col);
extern void DestroyDoubleMatrix(DMATRIX* dm);
extern double FindMin (std::vector<double>Vec);   //NB 4.9.05
extern double FindMax (std::vector<double>Vec);   //NB 4.9.05

/**
 * Author: NB 4.9.05
 * Finds all real roots of a third grade polynomial in the form:
 * P(x) = x^3 + px^2 + qx + r
 * employing the <a href="http://de.wikipedia.org/wiki/Kubische_Gleichung#Analytische_Bestimmung_reeller_L.C3.B6sungen">Cardano method</a>
 * @param p coefficient of cubic equation
 * @param q coefficient of cubic equation
 * @param r coefficient of cubic equation
 * @param t vector contains the roots
 */
extern void NsPol3 (double p, double q, double r, std::vector<double>* t);

typedef struct
{
	long row;
	long col;
	long** m;
} LMATRIX;

extern LMATRIX* CreateLongMatrix(long row, long col);
extern void DestroyLongMatrix(LMATRIX* lm);

#ifdef MFC                                        //WW
extern void CURRead(string);                      //WW
extern string ext_file_name;
#endif

/// Release the memory of arrary allocated by using 'new'. WW. 31.08.2010
template<class num> void  DeleteArray(num* an_array)
{
	if(an_array)
		delete [] an_array;
	an_array = NULL;
}
#endif
