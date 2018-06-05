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

/*JT: Send output message*/
extern void ScreenMessage(const char* message);
/* Gibt Programm-Abspann aus */
extern void DisplayMsg(const char* s);
/* Schreibt Zeichenkette ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayMsgLn(const char* s);
/* Schreibt Zeichenkette mit Zeilenvorschub auf Standardausgabe */
extern void DisplayMsgCR(const char* s);
/* Schreibt Zeichenkette mit Zeilenruecklauf auf Standardausgabe */
extern void DisplayDouble(double x, int i, int j);
/* Schreibt Double-Wert ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayLong(long x);
/* Schreibt Long-Wert ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayDoubleVector(double* vec, long grad, char* text);
/* Schreibt Vektor auf Standardausgabe */
// OK411 extern void DisplayDoubleMatrix ( double *mat, long m, long n, char *text );
/* Schreibt Matrix auf Standardausgabe */
extern void DisplayErrorMsg(const char* s);
/* Schreibt Fehlermeldung auf Standardausgabe */
extern void DisplayTimeMsg(const char* s, double d);
/* Schreibt Laufzeitmeldung auf Standardausgabe */
#endif
