/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: display.c
 */
/* Aufgabe:
   Enthaelt alle Funktionen fuer Standard Ein- und Ausgabe (Bildschirm,
   Tastatur)
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   10/1999     AH         Warnung entfernt
   01/2002     MK         Umleitung der DisplayX-Funktionen in MSG-Datei
                          Ausnahmen: DisplayStartMsg/DisplayEndMsg                                                  */
/**************************************************************************/
#include "display.h"
#include "makros.h"

extern FILE* OpenMsgFile(void);
extern void CloseMsgFile(FILE*);

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayMsg
 */
/* Aufgabe:
   Schreibt Zeichenkette ohne Zeilenvorschub auf Standardausgabe
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayMsg ( const char* s )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"%s",s);
	CloseMsgFile(f);
}

/**************************************************************************
 Task: Output message to screen. Helps to remove so many IFDEFS
 Programming:
  03/2012 JT
**************************************************************************/
void ScreenMessage(const char* message)
{
#ifdef USE_MPI
	if(myrank > 0)
		return;
#endif
	printf("%s",message);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayMsgLn
 */
/* Aufgabe:
   Schreibt Zeichenkette mit Zeilenvorschub auf Standardausgabe,
   beginnt immer erst nach 12 Zeichen Einrueckung.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayMsgLn ( const char* s )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"%s\n            ",s);
	fflush(stdout);
	CloseMsgFile(f);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayMsgCR
 */
/* Aufgabe:
   Schreibt Zeichenkette mit Zeilenruecklauf auf Standardausgabe,
   beginnt immer erst nach 12 Zeichen Einrueckung.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   08/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayMsgCR ( const char* s )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"%s\r            ",s);
	CloseMsgFile(f);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayDouble
 */
/* Aufgabe:
   Schreibt double-Wert ohne Zeilenvorschub formatiert auf
   Standardausgabe. Wird fuer beide Formatangaben 0 angegeben,
   wird im Standardformat geschrieben.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double x: double-Wert
   E int i: Gesamtstellenzahl
   E int j: Nachkommastellen
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   12/1995     cb         E-Format
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayDouble ( double x, int i, int j )
{
	FILE* f;
	f = OpenMsgFile();
	if ((i == 0) && (j == 0))
		/* printf("%f",x); */
		fprintf(f,"%g",x);
	else
		fprintf(f,"% *.*g",i,j,x);
	CloseMsgFile(f);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayLong
 */
/* Aufgabe:
   Schreibt long-Wert ohne Zeilenvorschub im Standardformat auf
   Standardausgabe.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long x: long-Wert
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayLong ( long x )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"%ld",x);
	CloseMsgFile(f);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayDoubleVector
 */
/* Aufgabe:
   Schreibt Vektor auf Standardausgabe
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *vec: Zeiger auf Vektor
   E long grad: Laenge des Vektors
   E char *text: Ueberschrift ueber Ausgabe
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   08/1994     Hans Herrmann    Erste Version
   12/1994     MSR              Von mathlib nach display portiert
   01/2002     MK               Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayDoubleVector ( double* vec, long grad, char* text )
{
	FILE* f;
	long i;
	DisplayMsgLn("");
	DisplayMsgLn(text);
	f = OpenMsgFile();
	for (i = 0; i < grad; i++)
	{
		fprintf(f,"| %+e |",vec[i]);
		fprintf(f,"%s\n            ","");
		fflush(stdout);
	}
	CloseMsgFile(f);
	DisplayMsgLn("");
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayErrorMsg
 */
/* Aufgabe:
   Schreibt Fehlermeldung auf Standardausgabe.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette (Fehlermeldung)
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayErrorMsg ( const char* s )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"\n!!!!!!!!  %s\n\n            ",s);
	CloseMsgFile(f);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayTimeMsg
 */
/* Aufgabe:
   Schreibt Laufzeitmeldung auf Standardausgabe.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette (Fehlermeldung)
   E double d: Zeitwert
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   07/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayTimeMsg ( const char* s, double d )
{
	FILE* f;
	f = OpenMsgFile();
	fprintf(f,"\n            %s%20ld s\n\n            ",s,((long) d));
	CloseMsgFile(f);
}
