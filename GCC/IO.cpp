/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include <cstdlib>
#include <cmath>
#include "IO.h"
#include "species.h"

using namespace std;

IO::IO(void) {}
IO::~IO(void) {}

/*  ====database convert subroutin===  */
int IO::file2code(string datafile_in, string codefile_out)
{
    int i;
    ifstream in;
    ofstream out;
    string instr;
    vector<string> pies;

    cout << endl << " Reading " << datafile_in << " file... ... ";
    in.open(datafile_in.c_str(), ios::in);
    out.open(codefile_out.c_str(), ios::out);
    out << "#include <string>" << endl;
    out << "using namespace std;" << endl;
    out << "const string " << codefile_out << "[]= {" << endl;
    if (in.good())
    {
        // cout << " OK " << endl;
        while (1)
        {
            if (in.eof())
                break;
            getline(in, instr);
            pies = IO::string2vector(instr);
            if ((int)pies.size() > 0 && pies[0].substr(0, 1) != "*")
            {
                out << "\"";
                for (i = 0; i < (int)pies.size(); i++)
                    out << pies[i] << " ";
                out << "\"," << endl;
            }
        }
        out << "\"END\"};";
    }
    else
        cout << endl << " Warning ! Could not open data file." << endl;
    in.close();
    out.close();
    return 1;
}
/*=========================================*/

int IO::file2vector(string datafile_in, vector<string>& vector_out)
{
    int i;
    ifstream in;
    string instr, pies0;
    vector<string> pies;

    // cout << endl << " Reading " << datafile_in << " file... ... " ;
    in.open(datafile_in.c_str(), ios::in);
    if (in.good())
    {
        // cout << " OK. " << endl;
        while (1)
        {
            if (in.eof())
                break;
            getline(in, instr);
            pies = IO::string2vector(instr);
            if ((int)pies.size() > 0 && pies[0].substr(0, 1) != "*")
            {
                pies0 = "";
                for (i = 0; i < (int)pies.size(); i++)
                    pies0 += pies[i] + " ";
                vector_out.push_back(pies0);
            }
        }
        vector_out.push_back("END");
    }
    else
        cout << endl << " Warning ! Could not open data file." << endl;
    in.close();
    return 1;
}

// split string line to pieces, and store in a vector
vector<string> IO::string2vector(string line)
{
    stringstream in;
    string sp;
    vector<string> pies;
    pies.clear();
    in.str(line);
    while (1)
    {
        if (in.eof())
            break;
        sp = "";
        in >> sp;
        if (sp != "")
            pies.push_back(sp);
    }
    return pies;
}

// 0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
// convert chemical formula to index
vector<int> IO::formula2index(std::string formula)
{
    const string Chemical_Element[10] = {"Li", "Na", "K", "Mg", "Ca",
                                         "Cl", "S",  "C", "H",  "O"};
    int i, j, n, symb_asc, is_bracket;
    vector<int> id_iz, type_x;
    vector<int> bracket_ia, bracket_ib, bracket_iz, bracket_ia0,
        bracket_ib0;  // CB clear
    vector<int> element_ia, element_ib, element_iz;
    vector<int> number_ia, number_ib, number_iz;
    vector<string> element_name;

    id_iz.clear();
    type_x.clear();
    n = int(formula.size());

    element_ia.clear();
    element_ib.clear();
    element_iz.clear();
    element_name.clear();

    number_ia.clear();
    number_ib.clear();
    number_iz.clear();

    bracket_ia.clear();  // to store begin position of bracket
    bracket_ib.clear();  // to store end position of bracket
    bracket_iz.clear();
    bracket_ia0.clear();
    bracket_ib0.clear();

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    for (i = 0; i < n; i++)
    {
        symb_asc = int(formula[i]);

        if (symb_asc >= 65 && symb_asc <= 90)  // A-Z
            type_x.push_back(0);
        else if (symb_asc >= 97 && symb_asc <= 122)  // a-z
            type_x.push_back(1);
        else if (symb_asc >= 48 && symb_asc <= 57)  // 0-9
            type_x.push_back(2);
        else if (symb_asc == 40 || symb_asc == 91)  // [(
            type_x.push_back(3);
        else if (symb_asc == 41 || symb_asc == 93)  // ])
            type_x.push_back(4);
        else
            type_x.push_back(5);
    }
    type_x.push_back(-1);  // set string end mark

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    // search for chemical elements
    for (i = 0; i < n; i++)
    {
        if (type_x[i] == 0)
            if (type_x[i + 1] == 1)
            {
                element_name.push_back(formula.substr(i, 2));
                element_ia.push_back(i);
                element_ib.push_back(i + 1);
            }
            else
            {
                element_name.push_back(formula.substr(i, 1));
                element_ia.push_back(i);
                element_ib.push_back(i);
            }
        else if (type_x[i] == 2)
        {
            if (i > 0 && type_x[i - 1] != 2)
                number_ia.push_back(i);
            if (i > 0 && type_x[i] == 2 && type_x[i + 1] != 2)
                number_ib.push_back(i);
        }
        else if (type_x[i] == 3)
            bracket_ia0.push_back(i);
        else if (type_x[i] == 4)
            bracket_ib0.push_back(i);
    }

    // to check the bracket symb is correct or not, e.g.  (( ... ))
    for (i = 0; i < (int)bracket_ia0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ib0.size(); j++)
        {
            if (i == (int)bracket_ia0.size() - 1)
            {
                if (bracket_ib0[j] > bracket_ia0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ib0[j] > bracket_ia0[i] &&
                    bracket_ib0[j] < bracket_ia0[i + 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ia.push_back(bracket_ia0[i]);
    }
    for (i = 0; i < (int)bracket_ib0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ia0.size(); j++)
        {
            if (i == 0)
            {
                if (bracket_ia0[j] < bracket_ib0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ia0[j] < bracket_ib0[i] &&
                    bracket_ia0[j] > bracket_ib0[i - 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ib.push_back(bracket_ib0[i]);
    }

    for (i = 0; i < (int)number_ia.size(); i++)
        number_iz.push_back(
            atoi(formula.substr(number_ia[i], number_ib[i] - number_ia[i] + 1)
                     .c_str()));

    for (i = 0; i < (int)element_name.size(); i++)
        if (type_x[element_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (element_ib[i] + 1 == number_ia[j])
                    element_iz.push_back(number_iz[j]);
            }
        else
            element_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        if (type_x[bracket_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (bracket_ib[i] + 1 == number_ia[j])
                    bracket_iz.push_back(number_iz[j]);
            }
        else
            bracket_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        for (j = 0; j < (int)element_ia.size(); j++)
            if (element_ia[j] > bracket_ia[i] && element_ib[j] < bracket_ib[i])
                element_iz[j] *= bracket_iz[i];

    for (i = 0; i < 10; i++)
    {
        id_iz.push_back(0);
        for (j = 0; j < int(element_ia.size()); j++)
            if (Chemical_Element[i] == element_name[j])
                id_iz[i] += element_iz[j];
    }

    // for(i=0;i<(int)id_iz.size();i++)
    //	cout << " " << i << " " << id_iz[i] << endl;
    return id_iz;
}

// convert chemical formula to index, from total list
vector<int> IO::formula2index_total(std::string formula)
{
    int i, j, n, symb_asc, is_bracket;
    vector<int> id_iz, type_x;
    vector<int> bracket_ia, bracket_ib, bracket_iz, bracket_ia0,
        bracket_ib0;  // CB clear
    vector<int> element_ia, element_ib, element_iz;
    vector<int> number_ia, number_ib, number_iz;
    vector<string> element_name;

    // get Chemical Element list from species.h
    vector<string> Chemical_Element, pies;
    Chemical_Element.clear();
    for (i = 0; (size_t)i < (sizeof(ELEMENTS) / sizeof(ELEMENTS[0])); i++)
    {
        pies = IO::string2vector(ELEMENTS[i]);
        Chemical_Element.push_back(pies[2]);
    }

    id_iz.clear();
    type_x.clear();
    n = int(formula.size());

    element_ia.clear();
    element_ib.clear();
    element_iz.clear();
    element_name.clear();

    number_ia.clear();
    number_ib.clear();
    number_iz.clear();

    bracket_ia.clear();  // to store begin position of bracket
    bracket_ib.clear();  // to store end position of bracket
    bracket_iz.clear();
    bracket_ia0.clear();
    bracket_ib0.clear();

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    for (i = 0; i < n; i++)
    {
        symb_asc = int(formula[i]);

        if (symb_asc >= 65 && symb_asc <= 90)  // A-Z
            type_x.push_back(0);
        else if (symb_asc >= 97 && symb_asc <= 122)  // a-z
            type_x.push_back(1);
        else if (symb_asc >= 48 && symb_asc <= 57)  // 0-9
            type_x.push_back(2);
        else if (symb_asc == 40 || symb_asc == 91)  // [(
            type_x.push_back(3);
        else if (symb_asc == 41 || symb_asc == 93)  // ])
            type_x.push_back(4);
        else
            type_x.push_back(5);
    }
    type_x.push_back(-1);  // set string end mark

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    // search for chemical elements
    for (i = 0; i < n; i++)
    {
        if (type_x[i] == 0)
            if (type_x[i + 1] == 1)
            {
                element_name.push_back(formula.substr(i, 2));
                element_ia.push_back(i);
                element_ib.push_back(i + 1);
            }
            else
            {
                element_name.push_back(formula.substr(i, 1));
                element_ia.push_back(i);
                element_ib.push_back(i);
            }
        else if (type_x[i] == 2)
        {
            if (i > 0 && type_x[i - 1] != 2)
                number_ia.push_back(i);
            if (i > 0 && type_x[i] == 2 && type_x[i + 1] != 2)
                number_ib.push_back(i);
        }
        else if (type_x[i] == 3)
            bracket_ia0.push_back(i);
        else if (type_x[i] == 4)
            bracket_ib0.push_back(i);
    }

    // to check the bracket symb is correct or not, e.g.  (( ... ))
    for (i = 0; i < (int)bracket_ia0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ib0.size(); j++)
        {
            if (i == (int)bracket_ia0.size() - 1)
            {
                if (bracket_ib0[j] > bracket_ia0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ib0[j] > bracket_ia0[i] &&
                    bracket_ib0[j] < bracket_ia0[i + 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ia.push_back(bracket_ia0[i]);
    }
    for (i = 0; i < (int)bracket_ib0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ia0.size(); j++)
        {
            if (i == 0)
            {
                if (bracket_ia0[j] < bracket_ib0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ia0[j] < bracket_ib0[i] &&
                    bracket_ia0[j] > bracket_ib0[i - 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ib.push_back(bracket_ib0[i]);
    }

    for (i = 0; i < (int)number_ia.size(); i++)
        number_iz.push_back(
            atoi(formula.substr(number_ia[i], number_ib[i] - number_ia[i] + 1)
                     .c_str()));

    for (i = 0; i < (int)element_name.size(); i++)
        if (type_x[element_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (element_ib[i] + 1 == number_ia[j])
                    element_iz.push_back(number_iz[j]);
            }
        else
            element_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        if (type_x[bracket_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (bracket_ib[i] + 1 == number_ia[j])
                    bracket_iz.push_back(number_iz[j]);
            }
        else
            bracket_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        for (j = 0; j < (int)element_ia.size(); j++)
            if (element_ia[j] > bracket_ia[i] && element_ib[j] < bracket_ib[i])
                element_iz[j] *= bracket_iz[i];

    for (i = 0; i < (int)Chemical_Element.size(); i++)
    {
        id_iz.push_back(0);
        for (j = 0; j < int(element_ia.size()); j++)
            if (Chemical_Element[i] == element_name[j])
                id_iz[i] += element_iz[j];
    }

    // for(i=0;i<(int)id_iz.size();i++)
    //	cout << " " << i << " " << id_iz[i] << endl;
    return id_iz;
}

vector<int> IO::formula2index_define(
    std::string formula, vector<string> Chemical_Element)
{
    // const string Chemical_Element[10]={"Li", "Na", "K", "Mg", "Ca", "Cl",
    // "S", "C", "H", "O"};
    int i, j, n, symb_asc, is_bracket;
    vector<int> id_iz, type_x;
    vector<int> bracket_ia, bracket_ib, bracket_iz, bracket_ia0,
        bracket_ib0;  // CB clear
    vector<int> element_ia, element_ib, element_iz;
    vector<int> number_ia, number_ib, number_iz;
    vector<string> element_name;

    id_iz.clear();
    type_x.clear();
    n = int(formula.size());

    element_ia.clear();
    element_ib.clear();
    element_iz.clear();
    element_name.clear();

    number_ia.clear();
    number_ib.clear();
    number_iz.clear();

    bracket_ia.clear();  // to store begin position of bracket
    bracket_ib.clear();  // to store end position of bracket
    bracket_iz.clear();
    bracket_ia0.clear();
    bracket_ib0.clear();

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    for (i = 0; i < n; i++)
    {
        symb_asc = int(formula[i]);

        if (symb_asc >= 65 && symb_asc <= 90)  // A-Z
            type_x.push_back(0);
        else if (symb_asc >= 97 && symb_asc <= 122)  // a-z
            type_x.push_back(1);
        else if (symb_asc >= 48 && symb_asc <= 57)  // 0-9
            type_x.push_back(2);
        else if (symb_asc == 40 || symb_asc == 91)  // [(
            type_x.push_back(3);
        else if (symb_asc == 41 || symb_asc == 93)  // ])
            type_x.push_back(4);
        else
            type_x.push_back(5);
    }
    type_x.push_back(-1);  // set string end mark

    // 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
    // search for chemical elements
    for (i = 0; i < n; i++)
    {
        if (type_x[i] == 0)
            if (type_x[i + 1] == 1)
            {
                element_name.push_back(formula.substr(i, 2));
                element_ia.push_back(i);
                element_ib.push_back(i + 1);
            }
            else
            {
                element_name.push_back(formula.substr(i, 1));
                element_ia.push_back(i);
                element_ib.push_back(i);
            }
        else if (type_x[i] == 2)
        {
            if (i > 0 && type_x[i - 1] != 2)
                number_ia.push_back(i);
            if (i > 0 && type_x[i] == 2 && type_x[i + 1] != 2)
                number_ib.push_back(i);
        }
        else if (type_x[i] == 3)
            bracket_ia0.push_back(i);
        else if (type_x[i] == 4)
            bracket_ib0.push_back(i);
    }

    // to check the bracket symb is correct or not, e.g.  (( ... ))
    for (i = 0; i < (int)bracket_ia0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ib0.size(); j++)
        {
            if (i == (int)bracket_ia0.size() - 1)
            {
                if (bracket_ib0[j] > bracket_ia0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ib0[j] > bracket_ia0[i] &&
                    bracket_ib0[j] < bracket_ia0[i + 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ia.push_back(bracket_ia0[i]);
    }
    for (i = 0; i < (int)bracket_ib0.size(); i++)
    {
        is_bracket = 0;
        for (j = 0; j < (int)bracket_ia0.size(); j++)
        {
            if (i == 0)
            {
                if (bracket_ia0[j] < bracket_ib0[i])
                    is_bracket = 1;
            }
            else
            {
                if (bracket_ia0[j] < bracket_ib0[i] &&
                    bracket_ia0[j] > bracket_ib0[i - 1])
                    is_bracket = 1;
            }
        }
        if (is_bracket == 1)
            bracket_ib.push_back(bracket_ib0[i]);
    }

    for (i = 0; i < (int)number_ia.size(); i++)
        number_iz.push_back(
            atoi(formula.substr(number_ia[i], number_ib[i] - number_ia[i] + 1)
                     .c_str()));

    for (i = 0; i < (int)element_name.size(); i++)
        if (type_x[element_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (element_ib[i] + 1 == number_ia[j])
                    element_iz.push_back(number_iz[j]);
            }
        else
            element_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        if (type_x[bracket_ib[i] + 1] == 2)
            for (j = 0; j < (int)number_iz.size(); j++)
            {
                if (bracket_ib[i] + 1 == number_ia[j])
                    bracket_iz.push_back(number_iz[j]);
            }
        else
            bracket_iz.push_back(1);

    for (i = 0; i < (int)bracket_ia.size(); i++)
        for (j = 0; j < (int)element_ia.size(); j++)
            if (element_ia[j] > bracket_ia[i] && element_ib[j] < bracket_ib[i])
                element_iz[j] *= bracket_iz[i];

    for (i = 0; i < (int)Chemical_Element.size(); i++)
    {
        id_iz.push_back(0);
        for (j = 0; j < int(element_ia.size()); j++)
            if (Chemical_Element[i] == element_name[j])
                id_iz[i] += element_iz[j];
    }

    // for(i=0;i<(int)id_iz.size();i++)
    //	cout << " " << i << " " << id_iz[i] << endl;
    return id_iz;
}

bool IO::file_compare(std::string in_file_new, std::string in_file_old)
{
    ifstream in_new, in_old;
    string instr_new, instr_old;
    bool res = true;

    in_new.clear();
    in_new.open(in_file_new.c_str(), ios::in);
    in_old.clear();
    in_old.open(in_file_old.c_str(), ios::in);

    if (in_new.good() && in_old.good())
    {
        while (1)
        {
            if (in_new.eof())
                break;
            getline(in_new, instr_new);
            getline(in_old, instr_old);
            if (instr_new != instr_old)
            {
                res = false;
                break;
            }
        }
    }
    else
        res = false;

    return res;
}

vector<string> IO::vector_reduce(std::vector<string> in_vec)
{
    int i;
    vector<string> out_vec;
    out_vec.clear();
    for (i = 0; i < (int)in_vec.size(); i++)
    {
        if (IO::string2vector(in_vec[i]).size() > 0)
            out_vec.push_back(in_vec[i]);
    }
    return out_vec;
}
