/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#if defined(WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include <makros.h>

#ifdef OGS_FEM_CAP  // CAP_REACT
#include "cacint.h"
#endif

#include "rf_react.h"
#include "rf_react_int.h"
#include "CAP_IO.h"
#include "IO.h"
using namespace std;

int CAP_MODE, CAP_Time, CAP_Node, CAP_icount, CAP_count;
vector<vector<std::string> > CHEM_STATE, CHEM_STATE_AP;
vector<vector<int> > PHASE_CONSTI_MAP;

char getPathSepatator()  // WW
{
#if defined(WIN32)
    return '\\';
#else
    return '/';
#endif
}

void write_file(void)
{
#ifdef OGS_FEM_CAP  // CAP_REACT
    stringstream ss;
    ofstream out;
    string file_name, path_name, copy_file;
    char pcname[TQSTRLEN], name_element[TQSTRLEN];
    long int np, npc, nsc, noerror;
    long int ip, ipc, isc;
    DB value;
    ip = ipc = isc = np = npc = nsc = noerror = 0;

    ss << FilePath << "_CHEM_STATE";
    ss >> path_name;
    ss.clear();

    if (CAP_count == 1)
    {
#if defined(WIN32)
        mkdir(path_name.c_str());

        copy_file = "copy " + FileName + ".* " + path_name;
        system(copy_file.c_str());

#else
        mkdir(path_name.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
    }

#if defined(WIN32)
    ss << FilePath << "_CHEM_STATE\\T" << CAP_Time << "N" << CAP_Node << "C"
       << CAP_icount;
#else
    ss << FilePath << "_CHEM_STATE\/T" << CAP_Time << "N" << CAP_Node << "C"
       << CAP_icount;
#endif
    ss >> file_name;
    // cout << file_name << "\n";
    out.open(file_name.c_str(), ios::out);
    out.clear();
    out << "phase constituent name a ia ac G";
    tqnop(&np, &noerror);
    for (ip = 1; ip <= np; ip++)
    {
        tqnopc(ip, &npc, &noerror);
        for (ipc = 1; ipc <= npc; ipc++)
        {
            out << endl;
            tqgnpc(ip, ipc, pcname, &noerror);
            out << ip << " " << ipc << " " << pcname;
            tqgetr("a", ip, ipc, &value, &noerror);
            out << " " << value;
            tqgetr("ia", ip, ipc, &value, &noerror);
            out << " " << value;
            tqgetr("ac", ip, ipc, &value, &noerror);
            out << " " << value;
            tqgdpc("G", ip, ipc, &value, &noerror);
            out << " " << value;
        }
    }
    out.close();
    file_name += "_AP";
    out.open(file_name.c_str(), ios::out);
    out.clear();
    out << "idx element amount";
    tqnosc(&nsc, &noerror);
    for (isc = 1; isc <= nsc; isc++)
    {
        out << endl;
        tqgnsc(isc, name_element, &noerror);
        out << isc << " " << name_element;
        tqgetr("ap", 2, isc, &value, &noerror);
        out << " " << value;
    }
    out.close();
#endif
}

bool CAP_check_file(void)
{
    stringstream ss;
    ifstream in;
    string file_name;
    ss.clear();
    bool success = false;
    int f;

    ss << FilePath << "_CHEM_STATE" << getPathSepatator() << "T" << CAP_Time
       << "N" << CAP_Node << "C" << CAP_icount;

    ss >> file_name;

    in.clear();
    in.open(file_name.c_str(), ios::in);
    if (in.good())
        success = true;

    in.close();
    in.clear();

    string proj_name, in_file_new, in_file_old;

    std::string locfilepath, locfilename;

    locfilepath = FilePath;
    locfilename = FileName;

    if (locfilepath == "")
    {
#if defined(WIN32)
        locfilepath = "";
        locfilename = locfilename;
#else
        locfilepath = "./";
        locfilename = "./" + locfilename;
#endif
    }
    proj_name = locfilename.substr(locfilepath.length());

    cout << proj_name << endl;
    cout << locfilename << endl;
    cout << locfilepath << endl;

    int i;
    vector<string> ext_name0, ext_name;

    ext_name0.clear();
    ext_name0.push_back("bc");
    ext_name0.push_back("cap");
    ext_name0.push_back("gli");
    ext_name0.push_back("ic");
    ext_name0.push_back("mcp");
    ext_name0.push_back("mfp");
    ext_name0.push_back("mmp");
    ext_name0.push_back("msh");
    ext_name0.push_back("msp");
    ext_name0.push_back("num");
    ext_name0.push_back("out");
    ext_name0.push_back("pcs");
    ext_name0.push_back("rei");
    ext_name0.push_back("rfd");
    ext_name0.push_back("rfe");
    ext_name0.push_back("st");
    ext_name0.push_back("tim");

    ext_name.clear();
    for (i = 0; i < (int)ext_name0.size(); i++)
    {
        in_file_new = locfilename + "." + ext_name0[i];
        in_file_old = locfilepath + getPathSepatator() + "_CHEM_STATE" +
                      getPathSepatator() + proj_name + "." + ext_name0[i];
        if (!IO::file_compare(in_file_new, in_file_old))
            ext_name.push_back(ext_name0[i]);
    }

    for (i = 0; i < (int)ext_name.size(); i++)
    {
        in_file_new = locfilename + "." + ext_name[i];
        in_file_old = locfilepath + getPathSepatator() + "_CHEM_STATE" +
                      getPathSepatator() + proj_name + "." + ext_name[i];
        f = OGS_keyword_check(in_file_new, in_file_old, ext_name[i]);
        if (f == 1)
            cout << in_file_new << " is removed " << endl;
        if (f == -1)
            cout << in_file_new << " is added " << endl;
    }

    return success;
}

int OGS_keyword_check(std::string in_file_new, std::string in_file_old,
                      std::string ext_name)
{
    ifstream in_new, in_old;
    stringstream ss;
    string instr_new, instr_old, instr, line_no, subkey;
    int i, j, res, icount;
    unsigned pos0, pos1;

    vector<string> key_new, key_old, subkey_new, subkey_old, line_new, line_old,
        pies;

    in_new.clear();
    in_new.open(in_file_new.c_str(), ios::in);
    in_old.clear();
    in_old.open(in_file_old.c_str(), ios::in);

    key_new.clear();
    key_old.clear();
    subkey_new.clear();
    subkey_old.clear();
    line_new.clear();
    line_old.clear();

    if (in_new.good() && in_old.good())
    {
        res = 0;
        getline(in_new, instr_new);
        getline(in_old, instr_old);
        icount = 1;

        while (1)
        {
            if (in_new.eof())
                break;
            icount++;
            ss.clear();
            ss << icount;
            ss >> line_no;

            getline(in_new, instr_new);
            pies = IO::string2vector(instr_new);
            if (pies.size() > 0 && pies[0].substr(0, 1) != ";")
            {
                if (pies[0].substr(0, 1) == "#")
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " +
                                pies[i];  // remove more blank in one line
                    instr = line_no + " " + instr;
                    key_new.push_back(instr);
                }
                else if (pies[0].substr(0, 1) == "$")
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " + pies[i];
                    subkey = instr;
                    instr = line_no + " " + instr;
                    subkey_new.push_back(instr);
                }
                else
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " + pies[i];
                    instr = line_no + " " + subkey + " " + instr;
                    line_new.push_back(instr);
                }
            }
        }

        icount = 1;
        while (1)
        {
            if (in_old.eof())
                break;
            icount++;
            ss.clear();
            ss << icount;
            ss >> line_no;

            getline(in_old, instr_old);
            pies = IO::string2vector(instr_old);
            if (pies.size() > 0 && pies[0].substr(0, 1) != ";")
            {
                if (pies[0].substr(0, 1) == "#")
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " +
                                pies[i];  // remove more blank in one line
                    instr = line_no + " " + instr;
                    key_old.push_back(instr);
                }
                else if (pies[0].substr(0, 1) == "$")
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " + pies[i];
                    subkey = instr;
                    instr = line_no + " " + instr;
                    subkey_old.push_back(instr);
                }
                else
                {
                    instr = "";
                    for (i = 0; i < (int)pies.size(); i++)
                        instr = instr + " " + pies[i];
                    instr = line_no + " " + subkey + " " + instr;
                    line_old.push_back(instr);
                }
            }
        }

        if (key_new.size() > key_old.size())
            cout << " add new # key word in " << ext_name << " file!" << endl;
        if (key_new.size() > key_old.size())
            cout << " remove # key word from " << ext_name << " file!" << endl;
        for (i = 0; i < (int)key_new.size(); i++)
        {
            for (j = 0; j < (int)key_old.size(); j++)
            {
                pos1 = key_new[i].find(" ");
                pos0 = key_old[j].find(" ");
                if (key_new[i].substr(pos1) == key_old[j].substr(pos0))
                {
                    key_new[i] = " ";
                    key_old[j] = " ";
                }
            }
        }

        for (i = 0; i < (int)subkey_new.size(); i++)
        {
            for (j = 0; j < (int)subkey_old.size(); j++)
            {
                pos1 = subkey_new[i].find(" ");
                pos0 = subkey_old[j].find(" ");
                if (subkey_new[i].substr(pos1) == subkey_old[j].substr(pos0))
                {
                    subkey_new[i] = " ";
                    subkey_old[j] = " ";
                }
            }
        }

        for (i = 0; i < (int)line_new.size(); i++)
        {
            for (j = 0; j < (int)line_old.size(); j++)
            {
                pos1 = line_new[i].find(" ");
                pos0 = line_old[j].find(" ");
                if (line_new[i].substr(pos1) == line_old[j].substr(pos0))
                {
                    line_new[i] = " ";
                    line_old[j] = " ";
                }
            }
        }

        key_new = IO::vector_reduce(key_new);
        key_old = IO::vector_reduce(key_old);
        subkey_new = IO::vector_reduce(subkey_new);
        subkey_old = IO::vector_reduce(subkey_old);
        line_new = IO::vector_reduce(line_new);
        line_old = IO::vector_reduce(line_old);

        for (i = 0; i < (int)key_new.size(); i++)
            cout << ext_name << " file: Line. " << key_new[i]
                 << "   -> not consistent !" << endl;
        for (i = 0; i < (int)key_old.size(); i++)
            cout << ext_name << " file: Line. " << key_old[i]
                 << "   -> changed !" << endl;
        for (i = 0; i < (int)subkey_new.size(); i++)
            cout << ext_name << " file: Line. " << subkey_new[i]
                 << "   -> not consistent !" << endl;
        for (i = 0; i < (int)subkey_old.size(); i++)
            cout << ext_name << " file: Line. " << subkey_old[i]
                 << "   -> changed !" << endl;
        for (i = 0; i < (int)line_new.size(); i++)
            cout << ext_name << " file: Line. " << line_new[i]
                 << "   -> not consistent !" << endl;
        for (i = 0; i < (int)line_old.size(); i++)
            cout << ext_name << " file: Line. " << line_old[i]
                 << "   -> changed !" << endl;
    }
    else
    {
        if (!in_new.good() && in_old.good())
            res = 1;
        else if (!in_old.good() && in_new.good())
            res = -1;
        else if (!in_new.good() && !in_old.good())
            res = 2;
    }
    return res;
}

void read_file(void)
{
    stringstream ss;
    ifstream in;
    string file_name, instr;
    vector<int> idx;
    long int ip, i;

    CHEM_STATE.clear();
    CHEM_STATE_AP.clear();

    ss.clear();
    ss << FilePath << "_CHEM_STATE" << getPathSepatator() << "T" << CAP_Time
       << "N" << CAP_Node << "C" << CAP_icount;
    ss >> file_name;

    // cout << file_name << "\n";
    // read file to vector
    in.clear();
    in.open(file_name.c_str(), ios::in);
    if (in.good())
    {
        while (1)
        {
            if (in.eof())
                break;
            getline(in, instr);
            // CHEM_STATE.push_back(REACT_PRQ::string2vector(instr));
            CHEM_STATE.push_back(REACTINT::string2vector(instr));
            std::cout << "CHEM_STATE pushback: size: " << CHEM_STATE.size()
                      << " CAP_Time: " << CAP_Time << " CaP_Node: " << CAP_Node
                      << " CAP_icount: " << CAP_icount << " f:" << file_name
                      << "\n";
        }
    }
    in.close();
    file_name += "_AP";
    in.clear();
    in.open(file_name.c_str(), ios::in);
    if (in.good())
    {
        while (1)
        {
            if (in.eof())
                break;
            getline(in, instr);
            // CHEM_STATE_AP.push_back(REACT_PRQ::string2vector(instr));
            CHEM_STATE_AP.push_back(REACTINT::string2vector(instr));
        }
    }
    in.close();

    if (CAP_count == 1)
    {
        PHASE_CONSTI_MAP.clear();
        idx.clear();
        ip = 0;
        for (i = 0; i < (int)CHEM_STATE.size(); i++)
        {
            if (atoi(CHEM_STATE[i][0].c_str()) > ip)
            {
                PHASE_CONSTI_MAP.push_back(idx);
                idx.clear();
                idx.push_back(i);
                ip = atoi(CHEM_STATE[i][0].c_str());
            }
            else
                idx.push_back(i);
        }
        PHASE_CONSTI_MAP.push_back(idx);
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqce(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
#else
int CAP_tqce(CHP, LI, LI, DBP, LIP)
#endif
{
    CAP_count++;
    CAP_icount++;
    if (CAP_MODE == 0)
    {               // normal chemapp
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqce(OPTION, INDEXP, INDEXC, VALS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 1)  // using chemapp to create chemical state files
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        int res = tqce(OPTION, INDEXP, INDEXC, VALS, NOERR);
        write_file();
        return res;
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        read_file();
        return 0;
    }
    else
        return -1;
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcel(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
#else
int CAP_tqcel(CHP, LI, LI, DBP, LIP)
#endif
{
    CAP_count++;
    CAP_icount++;
    if (CAP_MODE == 0)
    {               // normal chemapp
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcel(OPTION, INDEXP, INDEXC, VALS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 1)  // using chemapp to create chemical state files
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        int res = tqcel(OPTION, INDEXP, INDEXC, VALS, NOERR);
        write_file();
        return res;
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        read_file();
        return 0;
    }
    else
        return -1;  // CB should not occur
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcen(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
#else
int CAP_tqcen(CHP, LI, LI, DBP, LIP)
#endif
{
    CAP_count++;
    CAP_icount++;
    if (CAP_MODE == 0)
    {               // normal chemapp
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcen(OPTION, INDEXP, INDEXC, VALS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 1)  // using chemapp to create chemical state files
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        int res = tqcen(OPTION, INDEXP, INDEXC, VALS, NOERR);
        write_file();
        return res;
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        read_file();
        return 0;
    }
    else
        return -1;  // CB should not occur
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcenl(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
#else
int CAP_tqcenl(CHP, LI, LI, DBP, LIP)
#endif
{
    CAP_count++;
    CAP_icount++;
    if (CAP_MODE == 0)
    {               // normal chemapp
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcenl(OPTION, INDEXP, INDEXC, VALS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }                        // normal chemapp
    else if (CAP_MODE == 1)  // using chemapp to create chemical state files
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        int res = tqcenl(OPTION, INDEXP, INDEXC, VALS, NOERR);
        write_file();
        return res;
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        read_file();
        return 0;
    }
    else
        return -1;  // CB should not occur
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqgetr(CHP OPTION, LI INDEXP, LI INDEX, DBP VAL, LIP NOERR)
#else
int CAP_tqgetr(CHP OPTION, LI INDEXP, LI INDEX, DBP VAL, LIP)
#endif
{
    int i;
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {  // normal chemapp   or //using chemapp to create chemical state files
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqgetr(OPTION, INDEXP, INDEX, VAL, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        if (strcmp(OPTION, "ap") == 0)
            VAL[0] = atof(CHEM_STATE_AP[INDEX][2].c_str());
        else
        {
            for (i = 0; i < (int)CHEM_STATE[0].size(); i++)
                if (OPTION == CHEM_STATE[0][i])
                    break;
            if (INDEX > 0)
                VAL[0] = atof(
                    CHEM_STATE[PHASE_CONSTI_MAP[INDEXP][INDEX - 1]][i].c_str());
            else if (INDEX == 0)
                VAL[0] = atof(
                    CHEM_STATE[PHASE_CONSTI_MAP[INDEXP][INDEX]][i].c_str());
        }
        return 0;
    }
    else
        return -1;  // CB should not occur
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqgdpc(CHP OPTION, LI INDEXP, LI INDEXC, DBP VAL, LIP NOERR)
#else
int CAP_tqgdpc(CHP OPTION, LI INDEXP, LI INDEXC, DBP VAL, LIP)
#endif
{
    int i;
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {  // normal chemapp   or //using chemapp to create chemical state files
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqgdpc(OPTION, INDEXP, INDEXC, VAL, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else if (CAP_MODE == 2)  // file mode, without real chemapp
    {
        for (i = 0; i < (int)CHEM_STATE[0].size(); i++)
            if (OPTION == CHEM_STATE[0][i])
                break;
        if (INDEXC > 0)
            VAL[0] = atof(
                CHEM_STATE[PHASE_CONSTI_MAP[INDEXP][INDEXC - 1]][i].c_str());
        else if (INDEXC == 0)
            VAL[0] =
                atof(CHEM_STATE[PHASE_CONSTI_MAP[INDEXP][INDEXC]][i].c_str());
        return 0;
    }
    else
        return -1;  // CB should not occur
}

int CAP_tqini(LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqini(NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        read_file();
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqvers(LIP NVERS, LIP NOERR)
{
#else
int CAP_tqvers(LIP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqvers(NVERS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqopna(CHP FILE, LI LUN, LIP NOERR)
{
#else
int CAP_tqopna(CHP, LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqopna(FILE, LUN, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqopnb(CHP FILE, LI LUN, LIP NOERR)
{
#else
int CAP_tqopnb(CHP, LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqopnb(FILE, LUN, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqopnt(CHP FILE, LI LUN, LIP NOERR)
{
#else
int CAP_tqopnt(CHP, LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqopnt(FILE, LUN, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqclos(LI LUN, LIP NOERR)
{
#else
int CAP_tqclos(LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqclos(LUN, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqrfil(LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqrfil(NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}
int CAP_tqrbin(LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqrbin(NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqrcst(LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqrcst(NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqnosc(LIP NSCOM, LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqnosc(NSCOM, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        NSCOM[0] = (int)CHEM_STATE_AP.size() - 1;
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqgnsc(LI INDEXS, CHP NAME, LIP NOERR)
{
#else
int CAP_tqgnsc(LI, CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqgnsc(INDEXS, NAME, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqnop(LIP NPHASE, LIP NOERR)
{
#else
int CAP_tqnop(LIP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqnop(NPHASE, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqgnp(LI INDEXP, CHP NAME, LIP NOERR)
{
#else
int CAP_tqgnp(LI, CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqgnp(INDEXP, NAME, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqnopc(LI INDEXP, LIP NPCON, LIP NOERR)
{
#else
int CAP_tqnopc(LI, LIP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqnopc(INDEXP, NPCON, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqgnpc(LI INDEXP, LI INDEXC, CHP NAME, LIP NOERR)
{
#else
int CAP_tqgnpc(LI, LI, CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqgnpc(INDEXP, INDEXC, NAME, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqinp(CHP NAME, LIP INDEXP, LIP NOERR)
{
    size_t i;
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqinp(NAME, INDEXP, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        for (i = 0; i < CHEM_STATE.size(); i++)
        {
            if (strcmp(NAME, CHEM_STATE[i][2].c_str()) == 0)
                INDEXP[0] = atoi(CHEM_STATE[i][0].c_str());
        }
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqinpc(CHP NAME, LI INDEXP, LIP INDEXC, LIP NOERR)
{
    size_t i;
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqinpc(NAME, INDEXP, INDEXC, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        for (i = 0; i < CHEM_STATE.size(); i++)
        {
            if (strcmp(NAME, CHEM_STATE[i][2].c_str()) == 0 &&
                INDEXP == atoi(CHEM_STATE[i][0].c_str()))
                INDEXC[0] = atoi(CHEM_STATE[i][1].c_str());
        }
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcsp(LI INDEXP, CHP OPTION, LIP NOERR)
{
#else
int CAP_tqcsp(LI, CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcsp(INDEXP, OPTION, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR)
{
#else
int CAP_tqcspc(LI, LI, CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcspc(INDEXP, INDEXC, OPTION, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqstpc(LI INDEXP, LI INDEXC, DBP STOI, DBP WMASS, LIP NOERR)
{
#else
int CAP_tqstpc(LI, LI, DBP, DBP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqstpc(INDEXP, INDEXC, STOI, WMASS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcsc(CHP NAME, LIP NOERR)
{
#else
int CAP_tqcsc(CHP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcsc(NAME, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqinsc(CHP NAME, LIP INDEXS, LIP NOERR)
{
#else
int CAP_tqinsc(CHP, LIP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqinsc(NAME, INDEXS, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqcio(CHP OPTION, LI IVAL, LIP NOERR)
{
#else
int CAP_tqcio(CHP, LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqcio(OPTION, IVAL, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqremc(LI NUMCON, LIP NOERR)
{
#else
int CAP_tqremc(LI, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqremc(NUMCON, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

#ifdef OGS_FEM_CAP  // CAP_REACT
int CAP_tqsetc(CHP OPTION, LI INDEXP, LI INDEX, DB VAL, LIP NUMCON, LIP NOERR)
{
#else
int CAP_tqsetc(CHP, LI, LI, DB, LIP, LIP NOERR)
{
#endif
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqsetc(OPTION, INDEXP, INDEX, VAL, NUMCON, NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}

int CAP_tqshow(LIP NOERR)
{
    if (CAP_MODE == 0 || CAP_MODE == 1)
    {
#ifdef OGS_FEM_CAP  // CAP_REACT
        return tqshow(NOERR);
#else
        return -1;  // CB  should not occur
#endif
    }
    else
    {
        *NOERR = 0;
        return 0;
    }
}
