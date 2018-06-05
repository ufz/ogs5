/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: break.c
 */
/* Aufgabe:
   Ctrl-C - Ueberwachung
 */
/* Programmaenderungen:
   06/1994     MSR        Erste Version
 */
/**************************************************************************/

#include "break.h"
#include "makros.h"

void BreakFunc(int sig);
int abbruch = 0;

/**************************************************************************/
/* ROCKFLOW - Funktion: NoBreak
 */
/* Aufgabe:
   schaltet Ctrl-C aus
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
 */
/**************************************************************************/
void NoBreak(void)
{
	signal(SIGINT, SIG_IGN);
	abbruch = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StandardBreak
 */
/* Aufgabe:
   schaltet Ctrl-C - interpretation auf Standard
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
 */
/**************************************************************************/
void StandardBreak(void)
{
	signal(SIGINT, SIG_DFL);
	abbruch = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: SaveBreak
 */
/* Aufgabe:
   speichert eventuelles Ctrl-C in abbruch
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
 */
/**************************************************************************/
void SaveBreak(void)
{
	signal(SIGINT, BreakFunc);
	abbruch = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ClearBreak
 */
/* Aufgabe:
   setzt abbruch auf 0
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
 */
/**************************************************************************/
void ClearBreak(void)
{
	abbruch = 0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: BreakFunc
 */
/* Aufgabe:
   setzt abbruch auf 1 bei Ctrl-C
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E in sig: Signal; muss sein, ist aber anscheinend redundant
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
 */
/**************************************************************************/
void BreakFunc(int sig)
{
	if (sig == SIGINT)
		abbruch = 1;
}
