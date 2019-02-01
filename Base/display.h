/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: display.h
 */
/* Aufgabe:
   Enthaelt alle Funktionen fuer Standard Ein- und Ausgabe (Bildschirm,
   Tastatur)
 */
/**************************************************************************/

#ifndef display_INC

#define display_INC
/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */
#include <cstdio>
#include <cstring>
//#include <ctype.h>

namespace Display
{
/*JT: Send output message*/
void ScreenMessage(const char *format , ... );
void ScreenMessageNoMPIRank(const char *format , ... );

/* Gibt Programm-Abspann aus */
void DisplayMsg(const char* s);
/* Schreibt Zeichenkette ohne Zeilenvorschub auf Standardausgabe */
void DisplayMsgLn(const char* s);
/* Schreibt Zeichenkette mit Zeilenruecklauf auf Standardausgabe */
void DisplayDouble(double x, int i, int j);
/* Schreibt Double-Wert ohne Zeilenvorschub auf Standardausgabe */
void DisplayLong(long x);
/* Schreibt Vektor auf Standardausgabe */
// OK411 extern void DisplayDoubleMatrix ( double *mat, long m, long n, char
// *text );
/* Schreibt Matrix auf Standardausgabe */
void DisplayErrorMsg(const char* s);
}
#endif
