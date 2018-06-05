/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: memory.c

   Aufgabe:
   Ersetzt malloc, free und realloc
   Zugehoerige Schalter:
   MEMORY_TEST_IN_TIME  : Erstellt Bilanz der Anforderungen waehrend der
   Laufzeit (zeitintensiv!).
   MEMORY_SHOW_USAGE    : Zeigt aktuellen Bedarf bei jeder Anforderung an
   MEMORY_FLUSH         : Schreibt Protokoll ungepuffert (zeitintensiv!)
   MEMORY_REALLOC       : Ersetzt Realloc durch Malloc und Free und speichert
   Daten um
   MEMORY_STR           : Zusaetzliche Uebergabe einer Aufrufkennung bei jedem
   memory-Aufruf, um Aufrufstellen zu lokalisieren.
   MEMORY_ALLOCATION_TEST_SUCCESS
   : Prueft, ob eine Speicheranforderung erfolgreich
   absolviert wurde

   Programmaenderungen:
   12/1994     MSR        Erste Version
   10/1995     MSR        Erweiterung im Hinblick auf Auswertung der
   Ausgabedatei und Aufrufstellenlokalisierung
   5/1997     C.Thorenz  Erste Version des "in-time" Speichertests
   9/1997     C.Thorenz  Sichereres Handling fuer NULL-Zeiger und
   Zero-Groesse
   10/1998    C.Thorenz  Ueberpruefe den Erfolg der Aktionen
   7/1999    C.Thorenz  Bugfix

**************************************************************************/

#include "makros.h"
#include "memory.h"
#include <cstdio>
#include <cstdlib>

/* Interne (statische) Deklarationen */

#ifdef MEMORY_TEST_IN_TIME
#ifndef MEMORY_STR
typedef struct
{
	long address;
	long size;
} Memory_Table;

#else
typedef struct
{
	long address;
	long size;
	char* file;
	int line;
} Memory_Table;
#endif

static Memory_Table* memtab = NULL;

static long memory_list_size = 0;
static long memory_alloced = 0;
static long memory_max_alloced = 0;
#endif

/**************************************************************************
   ROCKFLOW - Funktion: Malloc (MAlloc)

   Aufgabe:
   Belegt Speicher auf dem Heap.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long bytes: Groesse des zu belegenden Speichers in Byte
   E char *datei: OPTIONAL, Name der aufrufenden Datei (--> MEMORY_STR)
   E int zeile: OPTIONAL, Zeile, aus der aufgerufen wird (--> MEMORY_STR)

   Ergebnis:
   Zeiger auf den reservierten Speicherbereich (NULL bei Fehler)

   Programmaenderungen:
   12/1994     MSR        Erste Version
   5/1997     C.Thorenz  Erste Version des "in-time" Speichertests
   9/1997     C.Thorenz  Sichereres Handling fuer NULL-Zeiger und
   Zero-Groesse

**************************************************************************/
#ifndef MEMORY_STR
void* Malloc(long bytes)
#else
void* MAlloc(long bytes, char* datei, int zeile)
#endif
{
	void* address = NULL;
#ifdef MEMORY_TEST_IN_TIME
	long i;
#endif

#ifdef MEMORY_MANAGEMENT_NOT_ANSI_COMPLIANT
	/* Laut ANSI nicht noetig, zur Sicherheit eingefuegt. CT */
	if (bytes == 0)
		return NULL;
#endif

	address = (void*)malloc((size_t)bytes);
#ifdef MEMORY_SHOW_USAGE
	printf("address %8x %ld  , ", address, (long)address);
#endif

#ifdef MEMORY_TEST_IN_TIME
	if (memtab)
	{
		/* Nach freiem Platz in der Speichertabelle suchen */
		i = 0;
		while ((i < memory_list_size) && (memtab[i].address > 0))
			i = i + 1;
		if (i == memory_list_size)
		{
			/* Kein freier Platz wurde gefunden, der Ueberhang wird benutzt.
			   Fuer den naechsten Durchlauf muss ein neuer Platz erzeugt werden. */
			memory_list_size = memory_list_size + 1;
			memtab = (Memory_Table*)realloc(memtab, (memory_list_size + 1) * sizeof(Memory_Table));
			memtab[memory_list_size].address = 0;
			memtab[memory_list_size].size = 0;
		}
		memory_alloced = memory_alloced + bytes;

		memtab[i].address = (long)address;
		memtab[i].size = (long)bytes;

#ifdef MEMORY_STR
		memtab[i].file = datei;
		memtab[i].line = zeile;
#endif
	}
#ifdef MEMORY_SHOW_USAGE
#ifdef MEMORY_STR
	printf("Malloc aufgerufen von %s Zeile %d.  ", datei, zeile);
#endif
	printf("Bytes im Speicher: %ld \n", memory_alloced);
#endif
	if (memory_alloced > memory_max_alloced)
		memory_max_alloced = memory_alloced;
#endif

#ifdef MEMORY_ALLOCATION_TEST_SUCCESS
	if ((bytes) && (address == 0)) /* Angeforderter Speicher wurde nicht geliefert */
	{
#ifdef MEMORY_STR
		printf("Malloc aufgerufen von %s Zeile %d. \n", datei, zeile);
#endif
		printf("Malloc ist fehlgeschlagen!!!\n");
	}
#endif

	return address;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Free (FRee)
 */
/* Aufgabe:
   Gibt belegten Speicher wieder frei.
   Die Zuweisung des Ergebnisses ist sinnvoll.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E void *block: Zeiger auf den freizugebenden Speicherbereich
   E char *datei: OPTIONAL, Name der aufrufenden Datei (--> MEMORY_STR)
   E int zeile: OPTIONAL, Zeile, aus der aufgerufen wird (--> MEMORY_STR)
 */
/* Ergebnis:
   NULL
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   5/1997     C.Thorenz  Erste Version des "in-time" Speichertests
   9/1997     C.Thorenz  Sichereres Handling fuer NULL-Zeiger und
   Zero-Groesse
 */
/**************************************************************************/
#ifndef MEMORY_STR
void* Free(void* block)
#else
void* FRee(void* block, char* datei, int zeile)
#endif
{
#ifdef MEMORY_TEST_IN_TIME
	long i;
#endif

#ifdef MEMORY_MANAGEMENT_NOT_ANSI_COMPLIANT
	/* Laut ANSI nicht noetig, zur Sicherheit eingefuegt. CT */
	if (block == NULL)
		return NULL;
#endif

#ifdef MEMORY_TEST_IN_TIME
	if (memtab)
	{
		/* Nach "Block" in der Speichertabelle suchen */
		i = memory_list_size - 1;
		/* Von hinten anfangen muesste schneller gehen, da sich die "alten"
		   Addressen (hoffentlich) am Anfang sammeln. */

		while (i >= 0)
		{
			if (memtab[i].address == (long)block)
			{
				/* Richtiger Block wurde gefunden! Kann freigegeben werden. */

				memory_alloced = memory_alloced - memtab[i].size;

				memtab[i].address = 0;
				memtab[i].size = 0;
				i = 0;
			}
			i = i - 1;
		}

#ifdef MEMORY_SHOW_USAGE
#ifdef MEMORY_STR
		printf("Free aufgerufen von %s Zeile %d.  ", datei, zeile);
#endif
		printf("Bytes im Speicher: %ld \n", memory_alloced);
#endif
	}
#endif

	free(block);

	return NULL;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Realloc (REalloc)
 */
/* Aufgabe:
   Aendert die Groesse eines Speicherbereiches auf dem Heap und
   uebernimmt dort zuvor gespeicherte Daten. Die Zuweisung des
   Ergebnisses ist zwingend.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E void *block: Zeiger auf den Speicherbereich
   E long bytes: Neue Groesse des Speicherbereichs in Byte
   E char *datei: OPTIONAL, Name der aufrufenden Datei (--> MEMORY_STR)
   E int zeile: OPTIONAL, Zeile, aus der aufgerufen wird (--> MEMORY_STR)
 */
/* Ergebnis:
   Zeiger auf den neuen Speicherbereich
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   5/1997     C.Thorenz  Erste Version des "in-time" Speichertests
   9/1997     C.Thorenz  Sichereres Handling fuer NULL-Zeiger und
   Zero-Groesse
 */
/**************************************************************************/
#ifndef MEMORY_STR
void* Realloc(void* block, long bytes)
#else
void* REalloc(void* block, long bytes, char* datei, int zeile)
#endif
{
	static void* address;
#ifdef MEMORY_REALLOC
	static char* a;
	static char* b;
#endif

#ifdef MEMORY_TEST_IN_TIME
	/* Es wird in den folgenden verschiedenen "#ifdef"-Abschnitten
	   sehr wuest mit den Eingangsgroessen hantiert, darum lege ich sie
	   erstmal ab ... */
	long old_block;
	long new_bytes;
	long i;

	old_block = (long)block;
	new_bytes = (long)bytes;
#endif

#ifdef MEMORY_MANAGEMENT_NOT_ANSI_COMPLIANT
	/* Laut ANSI nicht noetig, zur Sicherheit eingefuegt. CT */
	if (block == NULL)
#ifndef MEMORY_STR
		return Malloc(bytes);
#else
		return MAlloc(bytes, datei, zeile);
#endif
	if (bytes == 0)
#ifndef MEMORY_STR
		return Free(block);
#else
		return FRee(block, datei, zeile);
#endif
#endif

#ifdef MEMORY_REALLOC
#ifndef MEMORY_STR
	address = Malloc(bytes);
#else
	address = MAlloc(bytes, datei, zeile);
#endif
	if (!(block == NULL))
	{
		a = (char*)address;
		b = (char*)block;
		{
			/* beim Vergroessern wird ein laengerer Bereich gelesen,
			   als eigentlich belegt war !
			   --> nur mit eigener Speichertabelle zu aendern */
			register long i;
			i = 0l;
			while (i < bytes)
			{
				a[i] = b[i];
				i++;
			}
		}
	}
#ifndef MEMORY_STR
	block = Free(block);
#else
	block = FRee(block, datei, zeile); /* vorher Free msr 0796 */
#endif
#else // ifdef MEMORY_REALLOC
	address = (void*)realloc(block, (size_t)bytes);
#ifdef MEMORY_SHOW_USAGE
	printf("address %8x %ld  , ", address, (long)address);
#endif
#endif

#ifdef MEMORY_ALLOCATION_TEST_SUCCESS
	if ((bytes) && (address == 0)) /* Angeforderter Speicher wurde nicht geliefert */
	{
#ifdef MEMORY_STR
		printf("Malloc aufgerufen von %s Zeile %d. \n", datei, zeile);
#endif
		printf("Malloc ist fehlgeschlagen!!!\n");
	}
#endif

#ifdef MEMORY_TEST_IN_TIME
	if (memtab)
	{
		/* Nach "Block" in der Speichertabelle suchen */
		i = memory_list_size - 1;
		/* Von hinten anfangen muesste schneller gehen, da sich die "alten"
		   Addressen (hoffentlich) am Anfang sammeln. */

		while (i >= 0)
		{
			if (memtab[i].address == old_block)
			{
				/* Richtiger Block wurde gefunden! Kann freigegeben werden. */
				memory_alloced = memory_alloced - memtab[i].size + new_bytes;
				memtab[i].address = (long)address;
				memtab[i].size = (long)new_bytes;
#ifdef MEMORY_STR
				memtab[i].file = datei;
				memtab[i].line = zeile;
#endif
				i = 0;
			}
			i = i - 1;
		}
#ifdef MEMORY_SHOW_USAGE
#ifdef MEMORY_STR
		printf("Realloc aufgerufen von %s Zeile %d.  ", datei, zeile);
#endif
		printf("Bytes im Speicher: %ld \n", memory_alloced);
#endif
		if (memory_alloced > memory_max_alloced)
			memory_max_alloced = memory_alloced;
	}
#endif

	return address;
}
