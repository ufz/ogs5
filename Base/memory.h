/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: memory.h
 */
/* Aufgabe:
   Ersetzt malloc, free und realloc
   Zugehoerige Schalter:
     MEMORY_TEST_IN_TIME    : Pruefen der Speicherallokierung/-freigaben
                              waehrend der Laufzeit.
     MEMORY_STR     : Zusaetzliche Uebergabe einer Aufrufkennung bei jedem
                      memory-Aufruf, um Aufrufstellen zu lokalisieren.
     MEMORY_REALLOC : Ersetzt Realloc durch Malloc und Free und speichert
                      Daten um
 */
/**************************************************************************/

#ifndef memory_INC

#define memory_INC
/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */

#ifndef MEMORY_STR
extern void* Malloc(long bytes);
/* Ersetzt malloc */
extern void* Free(void* block);
/* Ersetzt free, muss in der Form >> x = Free(x); << benutzt werden */
extern void* Realloc(void* block, long bytes);
/* Ersetzt realloc, muss in der Form >> x = Realloc(x,i); << benutzt werden */
#else
extern void* MAlloc(long bytes, char* datei, int zeile);
extern void* FRee(void* block, char* datei, int zeile);
extern void* REalloc(void* block, long bytes, char* datei, int zeile);
/* Funktionen wie oben, nur mit Aufrufstellenkennzeichnung */
#endif
#endif
