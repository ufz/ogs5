/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: files0.c
 */
/* Aufgabe:
   Enthaelt die uebergeordneten Datei- Ein- und Ausgaberoutinen, sowie
   das Speichern der Durchbruchskurven.
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   02/1999     CT         Bugfix: Anpassung fuer mit 0 beginnende Elemente.
                          Kennung fuer Version im Kopf. Versionsabhaengiges
                          lesen.
   02/1999     CT         Korrekturen fuer herausgeloeste Materialgruppen
   03/1999     AH         Korrekturen nicht noetig, da die Materialzuweisung
                          nicht mehr auf dieser Ebene erfolgt. Die Abfrage ist
                          jetzt auskommentiert. Sie war hier nur um
   Kompatibilitaeten zum Rockflow aufzubewahren. 03/1999     CT anz_matxx,
   start_mat entfernt 02/2000     CT         Restart wieder hergestellt 07/2000
   AH         Vorbereitungen zum HGM. 9/2000     CT         Neu:
   RefreshNodeOutputData, Warnungen beseitigt 10/2001     AH         Inverse
   Modellierung Trennung und Anpassung (CreateFileData) In DestroyFileData
   Datenfeld auskommentiert. Neues Konzept fuer Datenbank-Verwaltung ah inv
   01/2002     MK         DisplayMsgX-Umleitung in *.msg-Datei: OpenMsgFile
   08/2002     MK         GetPathRFDFile
   ConfigFileData aus CreateFileData herausgeloest
   03/2003     RK         Quellcode bereinigt, Globalvariablen entfernt
 */
/**************************************************************************/

#include "display.h"
#include "memory.h"
//#include "makros.h"
//#ifndef NEW_EQS //WW. 07.11.2008
//#include "solver.h"
//#endif
//#include "rf_pcs.h"
//#include "rf_mmp_new.h"
#include "FileTools.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"
#include "rfmat_cp.h"
#include "tools.h"
//#include "rf_pcs.h"
#include "rf_out_new.h"
//#include "rf_tim_new.h"
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
//#include "rf_num_new.h"
#include "rf_fct.h"             //OK
#include "rf_fluid_momentum.h"  // PCH
#include "rf_kinreact.h"
#include "rf_random_walk.h"  // PCH
#include "rf_react.h"
#include "rf_react_int.h"
// CB2406 #ifdef OGS_FEM_CAP // CAP_REACT
// CB_merge_0513
#include "rf_react_cap.h"

#ifdef CHEMAPP
#include "eqlink.h"  //MX
#endif
#include "fct_mpi.h"
/* Tools */
//#include "mathlib.h"
//#include "femlib.h"
/* GeoLib */
//#include "geo_lib.h"
#include "files0.h"
// MSHLib
//#include "msh_lib.h"
//#include "gs_project.h"
/* Dateinamen */
char* crdat = NULL;     /*MX*/
char* file_name = NULL; /* dateiname */
static char* msgdat = NULL;

#define RFD_FILE_EXTENSION ".rfd"  // OK
#ifndef MFC                        // WW
void CURRead(std::string);         // OK
#endif
std::ios::pos_type CURReadCurve(std::ifstream*);  // OK
void CURWrite();                                  // OK

#define KEYWORD '#'
#define SUBKEYWORD '$'

// GEOLIB
//#include "GEOObjects.h"

// FileIO
#include "OGSIOVer4.h"
#include "readNonBlankLineFromInputStream.h"
//#include "FEMIO.h"

using namespace std;

static bool isValidTextFileFormat(const std::string& basename,
                                  const std::string& fext)
{
    const std::string fname(basename + fext);
    if (!IsFileExisting(fname))
        return true;
#ifdef _WIN32
    const bool is_win32 = true;
#else
    const bool is_win32 = false;
#endif
    if (is_win32 == HasCRInLineEnding(fname))
    {
        return true;
    }
    else
    {
        if (is_win32)
            std::cout << "*** ERROR: Detect UNIX file format " << fname.data()
                      << std::endl;
        else
            std::cout << "*** ERROR: Detect Windows file format "
                      << fname.data() << std::endl;
        return false;
    }
}

static bool checkFormatOfInputFiles(const std::string& basename)
{
    bool valid = true;
    valid &= isValidTextFileFormat(basename, ".gli");
    valid &= isValidTextFileFormat(basename, ".msh");
    valid &= isValidTextFileFormat(basename, ".pcs");
    valid &= isValidTextFileFormat(basename, ".ic");
    valid &= isValidTextFileFormat(basename, ".bc");
    valid &= isValidTextFileFormat(basename, ".st");
    valid &= isValidTextFileFormat(basename, ".mfp");
    valid &= isValidTextFileFormat(basename, ".msp");
    valid &= isValidTextFileFormat(basename, ".mmp");
    valid &= isValidTextFileFormat(basename, ".mcp");
    valid &= isValidTextFileFormat(basename, ".out");
    valid &= isValidTextFileFormat(basename, ".tim");

    return valid;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReadData
 */
/* Aufgabe:
   Liest Daten aus den Eingabedateien ein
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *dateiname: Dateiname ohne Extension
 */
/* Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1
 */
/* Programmaenderungen:
   07/1996    MSR     Erste Version
   07/2000    AH      Eingabe fuer Netzgenerator
   10/2001    AH      Trennung und Anpassung (CreateFileData)
                      Inverse Modellierung
                      Neues Konzept fuer Datenbank-Verwaltung
   10/2002   OK   DATCreateFileNames(dateiname)
                  DATDestroyFileNames()

   last modified: OK 16.10.2002
 */
/**************************************************************************/
int ReadData(const char* dateiname,
             GEOLIB::GEOObjects& geo_obj,
             std::string& unique_name)
{
#if defined(USE_MPI)  // WW
    if (myrank == 0)
    {
#endif
        std::cout << "\n";
        std::cout << "---------------------------------------------"
                  << "\n";
        std::cout << "Data input:"
                  << "\n";
#if defined(USE_MPI)  // WW
    }
#endif
    /* Dateinamen generieren */
    // OK  DATCreateFileNames(dateiname);
    static int datlen;
    datlen = (int)strlen(dateiname) + 5;
    crdat = (char*)Malloc(datlen); /*MX*/
    /*MX*/
    crdat = strcat(strcpy(crdat, dateiname), CHEM_REACTION_EXTENSION);
    msgdat = (char*)Malloc(datlen);
    msgdat = strcat(strcpy(msgdat, dateiname), RF_MESSAGE_EXTENSION);
    FILE* f = NULL;
    if ((f = fopen(msgdat, "r")) == NULL) /* MSG-Datei existiert nicht */
        msgdat = (char*)Free(msgdat);
    else
    {
        fclose(f);
        if ((f = fopen(msgdat, "a")) ==
            NULL) /* MSG-Schreibzugriff nicht moeglich */
            msgdat = (char*)Free(msgdat);
        else
            fclose(f);
    }
    //----------------------------------------------------------------------
    // Check line ending of input files
    if (!checkFormatOfInputFiles(dateiname))
    {
        Display::ScreenMessage("terminate this program");
        exit(0);
    }
    //----------------------------------------------------------------------
    // Read GEO data
    GEOLIB_Read_GeoLib(dateiname);

    std::string geo_file_name(dateiname);
    geo_file_name += ".gli";
    std::vector<std::string> file_read_errors;
    FileIO::readGLIFileV4(
        geo_file_name, &geo_obj, unique_name, file_read_errors);

    //----------------------------------------------------------------------
    // Read object data
    PCSRead(dateiname);
    MFPRead(dateiname);
    // HS PCS immediately followed by the MCP read
    CPRead(dateiname);  // SB:GS4
    BCRead(dateiname, geo_obj, unique_name);
    STRead(dateiname, geo_obj, unique_name);
    ICRead(dateiname, geo_obj, unique_name);
    OUTRead(dateiname, geo_obj, unique_name);
    TIMRead(dateiname);

    MSPRead(dateiname);
    MMPRead(dateiname);
    REACINTRead(dateiname);  // CB new reaction interface
    RCRead(dateiname);
    REACT_CAP_Read(dateiname,
                   geo_obj,
                   unique_name);  // DL/SB 11/2008 ChemASpp inteface new

    KRRead(dateiname, geo_obj, unique_name);
    KRWrite(dateiname);
#ifdef CHEMAPP
    CHMRead(dateiname);  // MX for CHEMAPP
#endif
    NUMRead(dateiname);

    FEMDeleteAll();  // KR moved from FEMRead()
    std::vector<CFEMesh*> mesh_vec;
    FEMRead(dateiname, mesh_vec, &geo_obj, &unique_name);
    if (!mesh_vec.empty())  // KR
    {
        fem_msh_vector.insert(fem_msh_vector.end(),
                              mesh_vec.begin(),
                              mesh_vec.end());  // re-inserted by KR
        CompleteMesh();                         // WW
    }

    // SBOK4209 MSHWrite(dateiname);
    // PCTRead is bounded by msh
    PCTRead(dateiname);  // PCH
    FMRead(dateiname);   // PCH
    FCTRead(dateiname);  // OK
    CURRead(dateiname);  // OK
// CURWrite(); //OK
#ifdef USE_PETSC
    FCT_MPI::FCTCommRead(dateiname);
#endif
    //----------------------------------------------------------------------
    // Read Excel/CVS data
    // PNTPropertiesRead(dateiname);

    msgdat = (char*)Free(msgdat);

    if (!mesh_vec.empty())
        return 100;  // magic number?

    return 1;
}

/**************************************************************************
   GeoSys-Method: FEMOpen->RFDOpen
   Task:
   Programing:
   11/2003 OK Implementation
   08/2004 OK PCS2
   01/2005 OK Boolean type
**************************************************************************/
bool RFDOpen(std::string file_name_base)
{
    (void)file_name_base;
    return false;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: *OpenMsgFile *CloseMsgFile
 */
/* Aufgabe:
   Oeffnet MSG-Datei fuer Display-Umleitung
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
 */
/* Ergebnis:
   - FILE -
 */
/* Programmaenderungen:
   12/2001     MK        Erste Version
 */
/**************************************************************************/
FILE* OpenMsgFile()
{
    FILE* f = NULL;
    if (msgdat)
    {
        if ((f = fopen(msgdat, "a")) == NULL)
        {
            f = stdout;
            fprintf(f,
                    "\n!!!!!!!!  %s\n\n            ",
                    "Fehler: Schreibzugriff auf Message-Protokolldatei nicht "
                    "moeglich!!");
        }
    }
    else
        f = stdout; /* Dateiname existiert nicht */
    return f;
}

void CloseMsgFile(FILE* f)
{
    if (f != stdout)
        if (fclose(f))
            Display::DisplayErrorMsg(
                "Fehler: Message-Protokolldatei konnte nicht geschlossen "
                "werden !!");
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void PRJRead(std::string base_file_name)
{
    char line[MAX_ZEILE];
    string sub_line;
    string line_string;
    ios::pos_type position;
    //========================================================================
    // file handling
    string rfd_file_name;
    rfd_file_name = base_file_name + FCT_FILE_EXTENSION;
    std::ifstream rfd_file(rfd_file_name.data(), std::ios::in);
    if (!rfd_file.good())
        return;
    rfd_file.seekg(0L, std::ios::beg);
    //========================================================================
    // keyword loop
    std::cout << "RFDRead"
              << "\n";
    while (!rfd_file.eof())
    {
        rfd_file.getline(line, MAX_ZEILE);
        project_title = line;
    }  // eof
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void CURRead(std::string base_file_name)
{
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //----------------------------------------------------------------------
    StuetzStellen* stuetz = NULL;
    anz_kurven = 1;
    stuetz = (StuetzStellen*)Malloc(sizeof(StuetzStellen));
    stuetz[0].punkt = 1.0;
    stuetz[0].wert = 1.0;
    kurven = (Kurven*)Malloc(sizeof(Kurven));
    kurven[anz_kurven - 1].anz_stuetzstellen = 1;
    kurven[anz_kurven - 1].stuetzstellen = stuetz;
    //----------------------------------------------------------------------
    // file handling
    std::string cur_file_name;
    cur_file_name = base_file_name + RFD_FILE_EXTENSION;
    std::ifstream cur_file(cur_file_name.data(), std::ios::in);
    if (!cur_file.good())
        return;
    cur_file.seekg(0L, std::ios::beg);
    //========================================================================
    // keyword loop

    Display::ScreenMessage("CURRead\n");

    while (!cur_file.eof())
    {
        cur_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != std::string::npos)
            return;
        //----------------------------------------------------------------------
        // keyword found
        if (line_string.find("#CURVE") != std::string::npos)
        {
            position = CURReadCurve(&cur_file);
            cur_file.seekg(position, std::ios::beg);
        }  // keyword found
    }      // eof
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
std::ios::pos_type CURReadCurve(std::ifstream* cur_file)
{
    bool new_keyword = false;
    std::string hash("#");
    std::string line_string;
    ios::pos_type position;
    std::stringstream line_stream;
    int anz = 0;
    double d1, d2;
    StuetzStellen* stuetz = NULL;
    //----------------------------------------------------------------------
    while (!new_keyword)
    {
        position = cur_file->tellg();
        // OK    cur_file->getline(buffer,MAX_ZEILE);
        // OK    line_string = buffer;
        line_string = GetLineFromFile1(cur_file);
        if (line_string.size() < 1)
            continue;
        //....................................................................
        // Test next keyword
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            continue;
        }
        //--------------------------------------------------------------------
        if (line_string.find(";") != string::npos)
            continue;
        //--------------------------------------------------------------------
        // DATA
        // OK    cur_file->seekg(position,ios::beg);
        // OK    *cur_file >> d1 >> d2;
        line_stream.str(line_string);
        line_stream >> d1 >> d2;
        anz++;
        stuetz = (StuetzStellen*)Realloc(stuetz, (anz * sizeof(StuetzStellen)));
        stuetz[anz - 1].punkt = d1;
        stuetz[anz - 1].wert = d2;
        line_stream.clear();
        //--------------------------------------------------------------------
    }
    //----------------------------------------------------------------------
    if (anz >= 1l)
    {
        anz_kurven++;
        kurven = (Kurven*)Realloc(kurven, (anz_kurven * sizeof(Kurven)));
        kurven[anz_kurven - 1].anz_stuetzstellen = anz;
        kurven[anz_kurven - 1].stuetzstellen = stuetz;
    }
    return position;
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void CURWrite()
{
    //========================================================================
    // File handling
    std::string fct_file_name = "test.cur";
    std::fstream fct_file(fct_file_name.c_str(), ios::trunc | ios::out);
    fct_file.setf(ios::scientific, ios::floatfield);
    fct_file.precision(12);
    if (!fct_file.good())
        return;
    fct_file << "GeoSys-CUR: Functions "
                "------------------------------------------------"
             << "\n";
    //========================================================================
    int j;
    StuetzStellen stuetz;
    for (int i = 0; i < anz_kurven; i++)
    {
        fct_file << "#CURVES"
                 << "\n";
        for (j = 0; j < kurven[i].anz_stuetzstellen; j++)
        {
            stuetz = kurven[i].stuetzstellen[j];
            fct_file << stuetz.punkt << " " << stuetz.wert << "\n";
        }
    }
    fct_file << "#STOP";
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetLineFromFile1
 */
/* Aufgabe:
   Liest aus dem Eingabefile *ein die n�chste Zeile
   F�ngt die Zeile mit ";" an oder ist sie leer, wird sie ausgelassen
   R�ckgabe ist ist ein string mit dem Zeileninhalt ab dem ersten
   Nicht-Leerzeichen bis zum ersten Auftreten des Kommentartzeichens ";"
 */
/* Programmaenderungen:
    05/2004     SB  First Version
 */
/*  09/2005     CC move from fem to geo
 **************************************************************************/
string GetLineFromFile1(ifstream* ein)
{
    //	return readNonBlankLineFromInputStream(*ein);
    string line, zeile = "";
    int fertig = 0, i = 0, j = 0;
    char zeile1[MAX_ZEILEN];
    line = "";  // WW
    //----------------------------------------------------------------------
    while (fertig < 1)
    {
        if (ein->getline(zeile1, MAX_ZEILEN))  // Zeile lesen
        {
            line = zeile1;  // character in string umwandeln
            i = (int)line.find_first_not_of(
                " ", 0);  // Anf�ngliche Leerzeichen �berlesen, i=Position des
                          // ersten Nichtleerzeichens im string
            j = (int)line.find(
                ";",
                i);  // Nach Kommentarzeichen ; suchen. j = Position des
                     // Kommentarzeichens, j=-1 wenn es keines gibt.
            if (j != i)
                fertig =
                    1;  // Wenn das erste nicht-leerzeichen ein Kommentarzeichen
                        // ist, zeile �berlesen. Sonst ist das
            // eine Datenzeile
            if ((i != -1))
                zeile = line.substr(
                    i, j - i);  // Ab erstem nicht-Leerzeichen bis
                                // Kommentarzeichen rauskopieren in neuen
            // substring, falls Zeile nicht leer ist
            i = (int)zeile.find_last_not_of(
                " ");  // Suche nach dem letzten Zeichen, dass kein Leerzeichen
                       // ist
            if (i >= 0)
            {
                //		  line.clear(); // = "";
                line = zeile.substr(
                    0, i + 1);  // Leerzeichen am Ende rausschneiden
                //		  zeile.clear(); // = "";
                zeile = line;
            }
        }
        else  // end of file found

            fertig = 1;
    }  // end while(...)
    //----------------------------------------------------------------------
    return zeile;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintString
 */
/* Aufgabe:
   Schreibt Zeichenkette ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E char *s: Zeichenkette
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintString(FILE* f, const char* s)
{
    if ((int)fprintf(f, "%s", s) != (int)strlen(s))
        return 0;
    return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintInt
 */
/* Aufgabe:
   Schreibt Integer-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E int x: Integer-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintInt(FILE* f, int x)
{
    if (fprintf(f, " %i ", x) < 0)
        return 0;
    return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintLong
 */
/* Aufgabe:
   Schreibt Long-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E long x: Long-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
 */
/**************************************************************************/
int FilePrintLong(FILE* f, long x)
{
    if (fprintf(f, " %ld ", x) < 0)
        return 0;
    return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: FilePrintDouble
 */
/* Aufgabe:
   Schreibt Double-Wert ohne Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
   E double x: Double-Wert
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   05/1994     MSR        Erste Version
   12/1995     cb         E-Format
 */
/**************************************************************************/
int FilePrintDouble(FILE* f, double x)
{
#ifdef FORMAT_DOUBLE
    if (fprintf(f, " % #*.*g ", FPD_GESAMT, FPD_NACHKOMMA, x) < 0)
        return 0;
#else
    if (fprintf(f, " % #g ", x) < 0)
        return 0;
#endif
    return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadDouble
 */
/* Aufgabe:
   Liest Double-Wert aus String und schreibt Protokoll in Datei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R double *x: gelesener Double-Wert
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestDouble func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrReadDouble(double* x, char* s, FILE* f, int* pos)
{
    *x = 0.0;
    if (sscanf(s, " %lf%n", x, pos) <= 0)
    {
        *pos = 0; /* nichts sinnvolles gelesen */
        fprintf(f,
                "\n %f      *** Fehler: Kein Wert eingelesen (double) !!!\n",
                *x);
        return 0;
    }
    else
    {
        /* CT: Protokolformat geaendert */
        if ((fabs(*x) < 100000.) && (fabs(*x) >= 0.1))
            fprintf(f, " %f ", *x);
        else
            fprintf(f, " %e ", *x);

        return 1;
    }
}

/**************************************************************************
   STRLib-Method: SubKeyword
   Task:
   Programing:
   09/2004 OK Implementation
   last modification:
**************************************************************************/
bool SubKeyword(const std::string& line)
{
    if (line.find(SUBKEYWORD) != std::string::npos)
        return true;
    else
        return false;
}

/**************************************************************************
   STRLib-Method: SubKeyword
   Task:
   Programing:
   09/2004 OK Implementation
   last modification:
**************************************************************************/
bool Keyword(const std::string& line)
{
    if (line.find(KEYWORD) != std::string::npos)
        return true;
    else
        return false;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrUp
 */
/* Aufgabe:
   wandelt Zeichenkette in Grossbuchstaben um
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   X char *s: umzuwandelnde Zeichenkette
 */
/* Ergebnis:
   umgewandelte Zeichenkette
 */
/* Programmaenderungen:
   03/1994   MSR   Erste Version
 */
/**************************************************************************/
char* StrUp(const char* s)
{
    int i;
    int l = (int)strlen(s);
    char* tmp = new char[l];
    strcpy(tmp, s);
    for (i = 0; i < l; i++)
        if (islower((int)s[i]))
            tmp[i] = (char)toupper((int)s[i]);
    return tmp;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StringReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Zeiger auf Adresse der gelesenen Zeichenkette; eine vorher
               vorhandene wird geloescht, es sollten nur mit malloc erzeugte
               Zeichenketten verwendet werden.
   E char *s: Zeichenkette, aus der gelesen werden soll
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   06/1999   OK   aus StrReadString
 */
/**************************************************************************/
int StringReadStr(char** x, char* s, int* pos)
{
    *x = NULL;
    //  *x = (char *) Malloc(256);
    *x = new char[256];  // CC
    *x[0] = '\0';
    if (sscanf(s, " %s%n", *x, pos) <= 0)
    {
        int a = (int)strlen(*x);  // CC
        // delete[] *x;//CC
        *x = new char[a + 1];  // CC
        //*x = (char *) Realloc(*x,((int)strlen(*x)+1));
        *pos = 0; /* nichts sinnvolles gelesen */
        return 0;
    }
    else
        return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: LineFeed
 */
/* Aufgabe:
   Schreibt Zeilenvorschub in Textdatei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E FILE *f: Dateizeiger
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int LineFeed(FILE* f)
{
    if (fprintf(f, "\n") < 0)
        return 0;
    return 1;
}

// int TFDouble ( double *x, FILE *f )
//{
//   return 1;
//}

// int TFString ( char* x, FILE* f )
//{
//	return 1;
//}

/**************************************************************************
   STRLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   last modification:
**************************************************************************/
void remove_white_space(std::string* buffer)
{
    int pos = 0;
    while (pos >= 0)
    {
        pos = (int)buffer->find_first_of(" ");
        if (pos < 0)
            break;
        buffer->erase(pos, 1);
    }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String und schreibt Protokoll in Datei
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Adresse der gelesenen Zeichenkett, es muss vorher Speicher
               allokiert werden
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestString func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   08/2000     CT        Erste Version
 */
/**************************************************************************/
int StrReadStr(char* x, char* s, FILE* f, /*FctTestString func,*/ int* pos)
{
    //   int test;
    x[0] = '\0';
    if (sscanf(s, " %s%n", x, pos) <= 0)
    {
        *pos = 0; /* nichts sinnvolles gelesen */
        fprintf(
            f, "\n %s      *** Fehler: Kein Wert eingelesen (string) !!!\n", x);
        return 0;
    }
    else
    {
        //      test = func(x,f);
        //      fprintf(f,"%s ",x);
        //      return test;
        fprintf(f, "%s ", x);
        return 1;
    }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrTestDouble
 */
/* Aufgabe:
   Testet, ob in s noch ein Double kommt;
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   0: nein; 1: ja
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrTestDouble(char* s)
{
    double i;
    if (sscanf(s, " %lf", &i) <= 0)
        return 0;
    else
        return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrTestHash
 */
/* Aufgabe:
   Testet, ob in s ein # folgt
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
   R int *pos: gelesene Zeichen bis nach dem # (wenn gefunden)
 */
/* Ergebnis:
   0: nein; 1: ja
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
int StrTestHash(char* s, int* pos)
{
    int p;
    char h[256];
    if (sscanf(s, " %s%n", h, &p) <= 0)
        return 0;
    else
    {
        if (strcmp(h, "#") == 0)
        {
            *pos = p;
            return 1;
        }
        else
            return 0;
    }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrOnlyReadStr
 */
/* Aufgabe:
   Liest Zeichenkette aus String aber schreibt Protokoll in Datei nicht
   (for Phreeqc read function)
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char **x: Adresse der gelesenen Zeichenkett, es muss vorher Speicher
               allokiert werden
   E char *s: Zeichenkette, aus der gelesen werden soll
   E FILE *f: Dateizeiger fuer Protokolldatei
   E FctTestString func: Funktionszeiger auf die Funktion, die den
                         eingelesenen Wert auf Gueltigkeit testet
   R int *pos: Anzahl der bisher gelesenen Zeichen
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   06/2003     MX        Erste Version
 */
/**************************************************************************/
/*MX*/
int StrOnlyReadStr(char* x,
                   char* s,
                   FILE* /*f*/,
                   /*FctTestString func,*/ int* pos)
{
    //   int test;

    x[0] = '\0';
    if (sscanf(s, " %s%n", x, pos) <= 0)
    {
        *pos = 0; /* nichts sinnvolles gelesen */
        return 0;
    }
    else
        //      test = func(x,f);
        //      return test;
        return 1;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: StrReadSubKeyword
 */
/* Aufgabe:
   Liest ein Sub-Keyword (eingeleitet mit "$") aus Keyword-String. Nur bis zum
   naechsten Hash (naechstes Keyword) oder Stringende
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R char *sub :  Zeichenkette, in die das Sub-Keywords geschrieben wird
   E char *s   : Zeichenkette, aus der gelesen werden soll
   E int begin : Erstes Zeichen, ab dem gesucht werden soll
   R int *found: Erstes Zeichen des Sub-Keywords im String s
   R int *ende : Letztes Zeichen des Sub-Keywords im String s, oder
                 Beginn des naechsten Keywords, oder Ende des Strings

 */
/* Ergebnis:
   1 bei gefundenem Sub-Keyword, sonst 0
 */
/* Programmaenderungen:
   08/2000 C.Thorenz  Erste Version
 */
/**************************************************************************/
int StrReadSubKeyword(char* sub, char* s, int beginn, int* found, int* ende)
{
    int i, xi = 0;

    *found = -1;
    *ende = (int)strlen(s);

    for (i = beginn; i < (int)strlen(s); i++)
    {
        if (s[i] == '$')
        {
            if (*found < 1)
                /* Anfang des Sub-Keywords merken */
                *found = i;
            else
            {
                /* Ende des Sub-Keywords merken (neues Sub-Keyword folgt) */
                *ende = i;
                break;
            }
        }

        if (s[i] == '#')
        {
            /* Ende des Sub-Keywords merken (neues Keyword folgt) */
            *ende = i;
            break;
        }

        if (*found >= 0)
        {
            sub[xi] = s[i];
            xi++;
        }
    }

    if (*found >= 0)
        sub[xi] = '\0';

    return *found >= 0;
}

/**************************************************************************
   STRLib-Method: get_sub_string
   Task: sub_string between pos1 and delimiter
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
std::string get_sub_string(const std::string& buffer,
                           const std::string& delimiter,
                           int pos1,
                           int* pos2)
{
    int pos = 0;
    std::string empty_string("");
    // string sub_string_this;
    *pos2 = (int)buffer.find(delimiter, pos1);
    if (*pos2 < 0)
        return empty_string;
    while (*pos2 <= pos1)
    {
        pos1++;
        *pos2 = (int)buffer.find(delimiter, pos1);
        if (*pos2 < 0)
        {
            *pos2 = (int)buffer.size();
            break;
        }
        if (pos1 >= (int)buffer.size())
            break;
    }
    string sub_string_this = buffer.substr(pos1, *pos2);
    while (pos >= 0)
    {
        pos = (int)sub_string_this.find_first_of(" ");
        if (pos < 0)
            break;
        sub_string_this.erase(pos, 1);
    }
    return sub_string_this;
}

/**************************************************************************
   STRLib-Method: get_sub_string2
   Task:
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
std::string get_sub_string2(const std::string& buffer,
                            const std::string& delimiter,
                            std::string* tmp)
{
    int pos2 = (int)buffer.find_first_of(delimiter);
    std::string sub_string = buffer.substr(0, pos2);
    *tmp = buffer.substr(pos2 + delimiter.size());
    return sub_string;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetUncommentedLine
 */
/* Aufgabe:
   R�ckgabe ist ist ein string mit dem Zeileninhalt ab dem ersten
   Nicht-Leerzeichen bis zum ersten Auftreten des Kommentartzeichens ";"
   Abgeleitet aus GetLineFromFile1() */
/* Programmaenderungen:
    06/2009     SB  First Version
 **************************************************************************/
std::string GetUncommentedLine(std::string line)
{
    std::string zeile = "";
    int i = 0, j = 0;
    //----------------------------------------------------------------------
    i = (int)line.find_first_not_of(
        " ", 0);  // Anf�ngliche Leerzeichen �berlesen, i=Position des ersten
                  // Nichtleerzeichens im string
    j = (int)line.find(";",
                       i);  // Nach Kommentarzeichen ; suchen. j = Position des
                            // Kommentarzeichens, j=-1 wenn es keines gibt.
    if ((i != -1))
        zeile =
            line.substr(i, j - i);  // Ab erstem nicht-Leerzeichen bis
                                    // Kommentarzeichen rauskopieren in neuen
    // substring, falls Zeile nicht leer ist
    i = (int)zeile.find_last_not_of(
        " ");  // Suche nach dem letzten Zeichen, dass kein Leerzeichen ist
    if (i >= 0)
    {
        line = zeile.substr(0, i + 1);  // Leerzeichen am Ende rausschneiden
        zeile = line;
    }

    return zeile;
}
