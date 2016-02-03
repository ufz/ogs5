/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ogs_display.h"

#include "BuildInfo.h"

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayStartMsg
 */
/* Aufgabe:
   Gibt Eroeffnungsbildschirm aus
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2010     NW         Automatic centering of the version information
 */
/**************************************************************************/
void DisplayStartMsg ( void )
{
	int i, pad_len;
	char buf[128];

	printf("\n");
	printf("          ###################################################\n");
	printf("          ##                                               ##\n");
	printf("          ##               OpenGeoSys-Project              ##\n");
	printf("          ##                                               ##\n");
	printf("          ##  Helmholtz Center for Environmental Research  ##\n");
	printf("          ##    UFZ Leipzig - Environmental Informatics    ##\n");
	printf("          ##                  TU Dresden                   ##\n");
	printf("          ##              University of Kiel               ##\n");
	printf("          ##            University of Edinburgh            ##\n");
	printf("          ##         University of Tuebingen (ZAG)         ##\n");
	printf("          ##       Federal Institute for Geosciences       ##\n");
	printf("          ##          and Natural Resources (BGR)          ##\n");
	printf("          ##  German Research Centre for Geosciences (GFZ) ##\n");
	printf("          ##                                               ##\n");

	//align the version information to center of the line
	printf("          ## ");
	sprintf(buf, "Version %s  Date %s", BuildInfo::OGS_VERSION.c_str(), BuildInfo::OGS_DATE.c_str());
	pad_len = 45 - (int)strlen(buf);
	for (i = 0; i < pad_len / 2; i++)
		printf(" ");
	printf("%s", buf);
	for (i = 0; i < pad_len - pad_len / 2; i++)
		printf(" ");
	printf(" ##\n");

	printf("          ##                                               ##\n");
	printf("          ###################################################\n");
	printf("\n          File name (without extension): ");
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayEndMsg
 */
/* Aufgabe:
   Gibt Programm-Abspann aus.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
void DisplayEndMsg ( void )
{
	printf("\n          Programm beendet!\n\n\n");
}
