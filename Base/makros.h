/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: makros.h
 */
/* Aufgabe:
   Diese Datei wird von allen ROCKFLOW-Modulen (*.c - Dateien !!!)
   importiert. Sie enthaelt globale Preprozessor-Definitionen und
   importiert weitere, (fast) ueberall benoetigte Header.

 */
/**************************************************************************/

#ifndef makros_INC
#define makros_INC

/* Global benoetigte Header */
//#include <stdlib.h>
/* Speicherverwaltung */
//#include <string.h>
#include <string>

/* Zeichenketten */
//#include <float.h>
/* Floating-Point */

/* ROCKFLOW-Version */
// LB: renamed ROCKFLOW_VERSION to OGS_VERSION and moved the #define to
// Base/Configure.h.in. Please set the version in the top-level CMakeLists.txt!!
// (see sources/CMakeLists.txt)

/* Definitionen von Makros zur Steuerung der bedingten Compilierung */
#define SWITCHES
/* Ausgabe der Schalterstellungen zu Beginn des Programms */
#ifdef MSVCPP6
#pragma warning(disable : 4786)
#endif

/* Laufzeitausgaben */
#define TESTTIME
#ifdef TESTTIME
#define TIMER_ROCKFLOW 0
#endif

/**********************************************************************/
/* Speicher */

#ifndef NO_ERROR_CONTROL /* Wird ggf. im Makefile gesetzt */
#define ERROR_CONTROL
/* Fehlertests (Feldgrenzen, Existenz o.ae.), die bei sauberen Netzen und
   einwandfrei funktionierendem Programm nichts bringen und nur Laufzeit
   kosten */
#endif

#define MEMORY_MANAGEMENT_NOT_ANSI_COMPLIANT
/*  Bei einigen Compilern werden
     malloc(0)
     realloc(NULL,xxx)
     realloc(xxx,0)
     free(NULL)
    nicht ANSI-gerecht gehandhabt. Mit diesem Schalter wird
    ANSI-Verhalten gewaehrleistet. */

#define noMEMORY_ALLOCATION_TEST_SUCCESS
/*  Prueft, ob eine Speicheranforderung erfolgreich absolviert wurde */

#define noMEMORY_TEST_IN_TIME /* fuer Versions-Speichertest */
/* Erstellt waehrend der Laufzeit eine Bilanz des allockierten
   und wieder freigegebenen Speichers. Sehr Zeitintensiv !!! */

#define noMEMORY_STR /* fuer Versions-Speichertest */
/* Gibt Informationen bei Memory-Funktionen zur
   Aufrufstellenlokalisation. Funktioniert nur zusammen mit
   MEMORY_TEST_IN_TIME. Sehr Speicherintensiv!!! */

#define noMEMORY_SHOW_USAGE
/* Gibt bei MEMORY_TEST_IN_TIME Informationen ueber jeden
   Malloc/Realloc/Free-Vorgang aus. */

#define noMEMORY_REALLOC
/* Ersetzt Realloc durch Malloc und Free und speichert um */

#ifndef MEMORY_TEST_IN_TIME
#ifdef MEMORY_SHOW_USAGE
#undef MEMORY_SHOW_USAGE
#endif
#ifndef MEMORY_ALLOCATION_TEST_SUCCESS
#ifdef MEMORY_STR
#undef MEMORY_STR
#endif
#endif
#endif

#ifdef MEMORY_STR
#define Malloc(a) MAlloc(a, __FILE__, __LINE__)
#define Free(a) FRee(a, __FILE__, __LINE__)
#define Realloc(a, b) REalloc(a, b, __FILE__, __LINE__)
/* Ersetzt Malloc, Free und Realloc in allen *.c-Dateien durch den
   erweiterten Aufruf mit Dateiname und Zeilennummer */
#endif

/**********************************************************************/
/* Daten */

#define noENABLE_ADT
/* Listen, Baeumen, etc sind dann moeglich */

/* Definitionen der Feldgroessen */
#define ELEM_START_SIZE 200l
/* Minimale Groesse des Elementverzeichnisses */
#define ELEM_INC_SIZE 1000l
/* Bei Erreichen von ELEM_START_SIZE wird das Elementverzeichnis
   automatisch um ELEM_INC_SIZE vergroessert */
#define NODE_START_SIZE 200l
/* Minimale Groesse des Knotenverzeichnisses */
#define NODE_INC_SIZE 1000l
/* Bei Erreichen von NODE_START_SIZE wird das Knotenverzeichnis
   automatisch um NODE_INC_SIZE vergroessert */
#define EDGE_START_SIZE 0l
/* Minimale Groesse des Kantenverzeichnisses */
#define EDGE_INC_SIZE 1000l
/* Bei Erreichen von EDGE_START_SIZE wird das Kantenverzeichnis
   automatisch um EDGE_INC_SIZE vergroessert */
#define PLAIN_START_SIZE 0l
/* Minimale Groesse des Flaechenverzeichnisses */
#define PLAIN_INC_SIZE 1000l
/* Bei Erreichen von PLAIN_START_SIZE wird das Flaechenverzeichnis
   automatisch um PLAIN_INC_SIZE vergroessert */

/**********************************************************************/
/* Protokolle */

/* Definitionen der Dateinamen-Erweiterungen */
#define TEXT_EXTENSION ".rfd"
/* Dateinamen-Erweiterung fuer Text-Eingabedatei */
#define PROTOCOL_EXTENSION ".rfe"
/* Dateinamen-Erweiterung fuer Text-Protokolldatei */
#define RF_INPUT_EXTENSION ".rfi"
/* Dateinamen-Erweiterung fuer RF-Input-Dateien */
#define RF_OUTPUT_EXTENSION ".rfo"
/* Dateinamen-Erweiterung fuer RF-Output-Dateien */
#define RF_MESSAGE_EXTENSION ".msg"
/* Dateinamen-Erweiterung fuer RF-Output-Dateien */
#define RF_SAVE_EXTENSION1 ".sv1"
/* Dateinamen-Erweiterung fuer RF-Sicherheitskopien */
#define RF_SAVE_EXTENSION2 ".sv2"
/* Dateinamen-Erweiterung fuer RF-Sicherheitskopien */
#define MESH_GENERATOR_EXTENSION ".rfm"
/* Dateinamen-Erweiterung fuer Text-Eingabedatei (Netzgenerator) */
#define MESH_GENERATOR_PROTOCOL_EXTENSION ".rfg"
/* Dateinamen-Erweiterung fuer Text-Protokolldatei (Netzgenerator) */
#define INVERSE_EXTENSION ".rfv" /* ah inv */
/* Dateinamen-Erweiterung fuer Text-Eingabedatei (Inverses Modellieren) */
#define INVERSE_PROTOCOL_EXTENSION ".rfp"
/* Dateinamen-Erweiterung fuer Text-Eingabedatei (Inverses Modellieren) */
#define CHEM_REACTION_EXTENSION ".pqc"
/* Dateinamen-Erweiterung fuer Text-Eingabedatei (Chemical reaction) */
#define CHEMAPP_REACTION_EXTENSION ".chm"
#define REACTION_EXTENSION_CHEMAPP ".cap" // DL/SB 11.2008
#define TEC_FILE_EXTENSION ".tec"
#define VTK_FILE_EXTENSION ".vtk" // GK
#define CSV_FILE_EXTENSION ".csv"

#define noTESTFILES
/* RFD-File Datenbank testen */

#define noEXT_RFD
/* Eingabeprotokoll ausfuehrlich kommentieren */
#define EXT_RFD_MIN
/* Eingabeprotokoll kommentieren, nur gefundene Schluesselworte */
#ifdef EXT_RFD
#undef EXT_RFD_MIN
#endif

/* Format der Double-Ausgabe ueber FilePrintDouble --> txtinout */
#define FORMAT_DOUBLE
#define FPD_GESAMT 4
#define FPD_NACHKOMMA 14

/**********************************************************************/
/* C1.2 Numerik */
#define noTESTTAYLOR
/* nur zu Testzwecken: Taylor-Galerkin-Verfahren nach Donea (1D) */

/**********************************************************************/
/* C1.4 Loeser */
#define noNULLE_ERGEBNIS
/* Nullt den Ergebnisvektor vor Aufruf des Loesers; ansonsten wird er
   mit den Ergebniswerten des letzten Zeitschritts vorbelegt. */
#define noSOLVER_SHOW_ERROR
/* Anzeigen des Iterationsfehlers */
#define noSOLVER_SHOW_RESULTS
/* Anzeigen der Iterationswerte */
#define RELATIVE_EPS
/* Abbruchkriterium bei CG-Loesern wird nicht als absolute Schranke
   benutzt, sondern mit der Norm der rechten Seite multipliziert */
/* Benutzte Normen bei Abbruchkriterien der CG-Loeser mit Speichertechnik:
     MVekNorm1 : Spaltensummennorm
     MVekNorm2 : euklidische Norm
     MVekNormMax : Maximumnorm
   Diese Normen muessen hier fuer die verschiedenen Loeser eingetragen
   werden !!!
 */
#define VEKNORM_BICG MVekNorm2
/* Norm fuer SpBICG-Loeser */
#define VEKNORM_BICGSTAB MVekNorm2
/* Norm fuer SpBICGSTAB-Loeser */
#define VEKNORM_QMRCGSTAB MVekNorm2
/* Norm fuer SpQMRCGSTAB-Loeser */
/*ahb*/
#define VEKNORM_CG MVekNorm2
/* Norm fuer SpCG-Loeser */
#define NORM 2
/* Norm fuer alle Objekte */
#if NORM == 0
#define VEKNORM MVekNormMax
#elif NORM == 1
#define VEKNORM MVekNorm1
#else
#define VEKNORM MVekNorm2
#endif
/*ahe*/

/**********************************************************************/
/* C1.9 Adaption */
#define noTEST_ADAPTIV
/* nur zu Testzwecken */
#define noREF_STATIC
/* erlaubt an ungefaehrlichen Stellen statische Variablen in Rekursionen.
   --> refine.c */

/**********************************************************************/
/* C1.10 Grafik */
#define no__RFGRAF
/* Bindet zur Zeit die Grafik-Funktionen unter X11 */

/**********************************************************************/
/* PCS / C++ */
#define PCS_OBJECTS
#define PCS_NUMBER_MAX 30
#define DOF_NUMBER_MAX 6 // JT: max # dof's per process
#define MAX_FLUID_PHASES 2 // JT: max # fluid phases
#define noPCS_NOD
//#define GLI // KR
#define noWINDOWS

/**********************************************************************/
/* Parallelization */
#define noPARALLEL
#define noCHEMAPP // MX
#define noREACTION_ELEMENT
#define noSX
#define noMPI
#define noOPEN_MP
/* Definitionen von Konstanten, die bei manchen Compilern benoetigt werden */
#ifndef NULL
#define NULL ((void*)0)
#endif
#ifndef TRUE
#define TRUE (0 == 0)
#endif
#ifndef FALSE
#define FALSE (1 == 0)
#endif
#ifndef PI
#define PI 3.14159265358979323846
#endif

/* Feste Zahlen fuer Genauigkeitspruefungen etc. */
#define Mdrittel (1.0 / 3.0)
#define MKleinsteZahl DBL_EPSILON
#define MFastNull DBL_MIN
#define MSqrt2Over3 sqrt(2.0 / 3.0)

/* CBLAS oder MKL_CBLAS verwenden? Wenn ja, wo? */
#define noCBLAS
#define noMKL_CBLAS

/* Includes fuer CBLAS oder MKL */
#ifdef MKL_CBLAS
#include <mkl_cblas.h>
#define CBLAS
#else
#ifdef CBLAS
#include <cblas.h>
#endif
#endif

/* Wo verwenden? */
#ifdef CBLAS
#define CBLAS_M2MatVek
#define CBLAS_MSkalarprodukt
#define CBLAS_MMultMatMat
#endif

/* UMF-Pack-Loeser verwenden? */
#ifdef UMFPACK31
#define UMFPACK
#endif
#ifdef UMFPACK40
#define UMFPACK
#endif

// min und max
#define Max(A, B) ((A) > (B) ? (A) : (B))

/*
   #ifndef min
   #define min(A,B) ((A) < (B) ? (A) : (B))
   #endif
   #ifndef max
   #define max(A,B) ((A) > (B) ? (A) : (B))
   #endif
 */
#define MAX_ZEILE 2048
/* max. Laenge einer UCD-Zeile; bei Leseproblemen vergroessern */

// enum DIS_TYPES {CONSTANT,LINEAR};

extern std::string FileName;
extern std::string FilePath; // WW

#define RESET_4410 // H2_ELE test

//---- MPI Parallel --------------
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_BRNS) || defined(USE_MPI_KRC) || defined(USE_PETSC)
extern int mysize; // WW
extern int myrank;
#endif
//---- MPI Parallel --------------

#endif
