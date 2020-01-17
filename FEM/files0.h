/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef strings_INC
#define strings_INC

#define FORMAT_DOUBLE
#define FPD_GESAMT 4
#define FPD_NACHKOMMA 14

/* Andere oeffentlich benutzte Module */
#include <cstdio>
#include <string>
//#include <fstream>

typedef int (*FctTestInt)(int*, FILE*);
typedef int (*FctTestLong)(long*, FILE*);
typedef int (*FctTestFloat)(float*, FILE*);
typedef int (*FctTestDouble)(double*, FILE*);
typedef int (*FctTestString)(char*, FILE*);
/* Funktionsprototypen zum Testen von eingelesenen Werten auf
   zulaessige Wertebereiche, mit Protokolldatei; Ergebnis: 0 bei Fehler,
   der zum Abbruch fuehren soll, sonst 1. Wenn korrigiert wurde: Protokoll
   schreiben und Ergebnis 1 */

extern int LineFeed(FILE* f);
/* Schreibt Zeilenvorschub in Textdatei */
extern int FilePrintString(FILE* f, const char* s);
/* Schreibt Zeichenkette in Textdatei */
extern int FilePrintInt(FILE* f, int x);
/* Schreibt Integer-Wert in Textdatei */
extern int FilePrintLong(FILE* f, long x);
/* Schreibt Long-Wert in Textdatei */
extern int FilePrintDouble(FILE* f, double x);
/* Schreibt Double-Wert in Textdatei */
extern int FilePrintIntWB(FILE* f, int x);
/* Schreibt Integer-Wert ohne fuehrende Leerzeichen in Textdatei */
extern int FilePrintLongWB(FILE* f, long x);
/* Schreibt Long-Wert ohne fuehrende Leerzeichen in Textdatei */
extern int FilePrintDoubleWB(FILE* f, double x);
/* Schreibt Double-Wert ohne fuehrende Leerzeichen in Textdatei */
extern int FileCommentDoubleMatrix(FILE* f, double* d, int adjust, long m,
                                   long n);
/* Schreibt Double-m x n-Matrix als Kommentar mit Einrueckung adjust
   in Textdatei */

extern int StrReadSubSubKeyword(char* x, char* s, int beginn, int* found,
                                int* ende);
/* Liest Sub-Sub-Keyword aus String */
extern int StrReadSubKeyword(char* x, char* s, int beginn, int* found,
                             int* ende);
/* Liest Sub-Keyword aus String */

extern int StrReadInt(int* x, char* s, FILE* f, FctTestInt func, int* pos);
/* Liest Integer-Wert aus String und schreibt Protokoll in Datei */
extern int StrReadLong(long* x, char* s, FILE* f, FctTestLong func, int* pos);
/* Liest Long-Wert aus String und schreibt Protokoll in Datei */
extern int StrReadFloat(float* x, char* s, FILE* f, FctTestFloat func,
                        int* pos);
/* Liest Float-Wert aus String und schreibt Protokoll in Datei */
extern int StrReadDouble(double* x, char* s, FILE* f, int* pos);
/* Liest Double-Wert aus String und schreibt Protokoll in Datei */
extern int StrReadString(char** x, char* s, FILE* f,
                         /*FctTestString func,*/ int* pos);
/* Liest Zeichenkette aus String und schreibt Protokoll in Datei */
extern int StrReadStr(char* x, char* s, FILE* f,
                      /*FctTestString func,*/ int* pos);
/* Liest Zeichenkette aus String und schreibt Protokoll in Datei */
extern int StrTestInt(char* s);
/* Testet, ob in s noch ein Integer kommt; 0:=nein */
extern int StrTestLong(char* s);
/* Testet, ob in s noch ein LongInt kommt; 0:=nein */
extern int StrTestFloat(char* s);
/* Testet, ob in s noch ein Float kommt; 0:=nein */
extern int StrTestDouble(char* s);
/* Testet, ob in s noch ein Double kommt; 0:=nein */
extern int StrTestString(char* s);
/* Testet, ob in s noch ein String kommt; 0:=nein */
extern int StrTestHash(char* s, int* pos);
/* Testet, ob in s ein # folgt; 0:=nein; 1:=ja, pos liefert Position
   nach dem # */
extern int StrTestDollar(char* s, int* pos);
/* Testet, ob in s ein $ folgt; 0:=nein; 1:=ja, pos liefert Position
   nach dem $ */
extern int StrTestInv(char* s, int* pos);
/* Testet, ob in s ein ? folgt; 0:=nein; 1:=ja, pos liefert Position
   nach dem ? */
extern int StrReadLongNoComment(long* x, char* s, FILE* f, FctTestLong func,
                                int* pos);
extern int StrReadStrNoComment(char* x, char* s, FILE* f, FctTestString func,
                               int* pos);

/* zusaetzliche Lesefunktionen fuer RF-SHELL */
extern int StringReadFloat(float* x, char* s, int* pos);
extern int StringReadDouble(double* x, char* s, int* pos);
extern int StringReadInt(int* x, char* s, int* pos);
extern int StringReadLong(long* x, char* s, int* pos);
extern int StringReadStr(char** x, char* s, int* pos);

/* allgemeine Dummy-Testfunktionen */
// extern int TFInt ( int *x, FILE *f );
// extern int TFLong ( long *x, FILE *f );
// extern int TFFloat ( float *x, FILE *f );
// extern int TFDouble ( double *x, FILE *f );
// extern int TFDoubleNew (char *s, FILE *f );
extern int TFString(char* x, FILE* f);

/* Liest Zeichenkette von Standardeingabe */
extern char* StrUp(const char* s);
/* wandelt s in Grossbuchstaben um */
extern char* StrDown(char* s);
/* wandelt s in Kleinbuchstaben um */

extern void GetRFINodesData();
/* Weitere externe Objekte */
#define MAX_NAME 80

/*MX*/
extern int StrOnlyReadStr(char* x, char* s, FILE* f,
                          /*FctTestString func,*/ int* pos);

extern std::string get_sub_string(const std::string&, const std::string&, int,
                                  int*);
extern void remove_white_space(std::string*);
// extern std::string get_sub_string2(std::string buffer,std::string
// delimiter,std::string cut_string);
extern std::string get_sub_string2(const std::string&, const std::string&,
                                   std::string*);
extern bool SubKeyword(const std::string&);
extern bool Keyword(const std::string&);

// CC move here
extern std::string GetLineFromFile1(std::ifstream*);
// SB
extern std::string GetUncommentedLine(std::string);
// extern std::string NumberToString(long);
extern void is_line_empty(std::string*);  // OK
#endif
