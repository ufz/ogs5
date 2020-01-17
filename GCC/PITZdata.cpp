/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdlib>

#include "PITZdata.h"
#include "species.h"
#include "IO.h"
#include "VLE.h"

using namespace std;

PITZdata::PITZdata(void) {}
PITZdata::~PITZdata(void) {}

double PITZdata::charge(string N)
{
    int i, flag = 0;
    double res;
    vector<string> pies;
    for (i = 0; (size_t)i < (sizeof(SPECIES) / sizeof(SPECIES[0])); i++)
    {
        pies = IO::string2vector(SPECIES[i]);
        if (pies[0] == N)
        {
            res = atof(pies[1].c_str());
            flag = 1;
        }
    }
    if (flag == 0)
    {
        cout << " No Enough Species In The List " << endl;
        exit(0);
    }
    return res;
}

double PITZdata::pitzer_parameters(double T, double P, string param_switch)
{
    double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0, c6 = 0.0, c7 = 0.0,
           c8 = 0.0, c9 = 0.0, c10 = 0.0;
    double c11 = 0.0 /*, c12=0.0, c13=0.0, c14=0.0, c15=0.0, c16=0.0, c17=0.0,
                        c18=0.0, c19=0.0, c20=0.0*/
        ;
    // double c21=0.0, c22=0.0, c23=0.0, c24=0.0, c25=0.0, c26=0.0, c27=0.0,
    // c28=0.0, c29=0.0, c30=0.0; double DW;
    double res = 0.0;

    if (param_switch == "LAMN_CO2_CO2")
    {
        c1 = -8.603471564E-01;
        c2 = 3.297141654E-03;
        c3 = 6.309267405E+01;
        c4 = -4.098960500E-06;
        c5 = 1.529493614E+01;
        c6 = 6.506644253E-03;
        c7 = -9.637977140E-04;
        c8 = -3.238222665E-01;
        c9 = 1.599113719E-02;
        c10 = 0.0;
        c11 = -1.886733300E-05;
        return c1 + c2 * T + c3 / T + c4 * T * T + c5 / (630.0e0 - T) + c6 * P +
               c7 * P * log(T) + c8 * P / T + c9 * P / (630.0e0 - T) +
               c10 * P * P / pow((630.0e0 - T), 2) + c11 * T * log(P);
    }
    else if (param_switch == "LAM_Na_CO2")
    {
        c1 = -2.739092216E-01;
        c2 = 7.399855859E-04;
        c3 = 5.552132850E+01;
        c4 = 0.0;
        c5 = 0.0;
        c6 = 0.0;
        c7 = 0.0;
        c8 = 5.683638727E-03;
        c9 = -8.009093476E-04;
        c10 = 0.0;
        c11 = -1.745620270E-05;
        return c1 + c2 * T + c3 / T + c4 * T * T + c5 / (630.0e0 - T) + c6 * P +
               c7 * P * log(T) + c8 * P / T + c9 * P / (630.0e0 - T) +
               c10 * P * P / pow((630.0e0 - T), 2) + c11 * T * log(P);
    }
    else

        //+++++++++=====test again=======+++++++++++
        if (param_switch == "ZETA_NaCl_CO2")
    {  // not "ZETA_Na_Cl_CO2", only for CO2 activity coeff
        c1 = -1.665719188E-02;
        c2 = 1.391618600E-06;
        c3 = 0.0;
        c4 = 0.0;
        c5 = 0.0;
        c6 = 0.0;
        c7 = 0.0;
        c8 = -1.873812115E-03;
        c9 = -1.577400757E-03;
        c10 = 0.0;
        c11 = 0.0;
        return c1 + c2 * T + c3 / T + c4 * T * T + c5 / (630.0e0 - T) + c6 * P +
               c7 * P * log(T) + c8 * P / T + c9 * P / (630.0e0 - T) +
               c10 * P * P / pow((630.0e0 - T), 2) + c11 * T * log(P);
    }

    else
        return res;
}
