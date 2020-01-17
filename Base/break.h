/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: break.h
 */
/* Aufgabe:
   Ctrl-C - Ueberwachung
 */
/**************************************************************************/

#ifndef break_INC

#define break_INC
/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */
#include <signal.h>

/* Deklarationen */
extern void NoBreak(void);
/* schaltet Ctrl-C aus */
extern void StandardBreak(void);
/* schaltet Ctrl-C - interpretation auf Standard */
extern void SaveBreak(void);
/* speichert eventuelles Ctrl-C in abbruch */
extern void ClearBreak(void);
/* setzt abbruch auf 0 */

/* Weitere externe Objekte */

extern int abbruch;
/* 1, wenn Ctrl-C gedrueckt wurde */
#endif
