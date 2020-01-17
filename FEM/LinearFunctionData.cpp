/*
 * LinearFunctionData.cpp
 *
 *  Created on: Sep 1, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cstring>
#include <sstream>
#include <vector>

#include "LinearFunctionData.h"

// TF do this only in cpp files !!!
using namespace std;

LinearFunctionData::LinearFunctionData(ifstream& ins, int num_var)
    : _ndata(0), _subdom_index(NULL), _a0(NULL), _b0(NULL), _c0(NULL), _d0(NULL)
{
    bool is_sub_domain = false;
    if (num_var > 0)
    {
        is_sub_domain = true;
        _ndata = num_var;
    }
    else
        _ndata = 1;

    char* pch;
    char seps[] = "+\n";
    char seps1[] = "*";
    double f_buff;
    _a0 = new double[_ndata];
    _b0 = new double[_ndata];
    _c0 = new double[_ndata];
    _d0 = new double[_ndata];

    if (is_sub_domain)
        _subdom_index = new size_t[_ndata];

    for (size_t i = 0; i < _ndata; i++)
    {
        _a0[i] = _b0[i] = _c0[i] = _d0[i] = 0.0;

        vector<string> tokens;
        string str_buff;
        stringstream buff;
        size_t ibuf(0);

        if (is_sub_domain)
        {
            ins >> ibuf >> str_buff >> ws;
            _subdom_index[i] = ibuf;
        }
        else
            ins >> str_buff >> ws;

        pch = strtok(const_cast<char*>(str_buff.c_str()), seps);
        buff << pch;
        buff >> _a0[i];
        buff.clear();
        while (pch != NULL)
        {
            pch = strtok(NULL, seps);
            if (pch == NULL)
                break;
            string token = pch;
            tokens.push_back(token);
        }
        for (size_t k = 0; k < tokens.size(); k++)
        {
            pch = strtok(const_cast<char*>(tokens[k].c_str()), seps1);
            buff << pch;
            buff >> f_buff;
            buff.clear();
            pch = strtok(NULL, seps1);
            switch (pch[0])
            {
                case 'x':
                    _b0[i] = f_buff;
                    break;
                case 'y':
                    _c0[i] = f_buff;
                    break;
                case 'z':
                    _d0[i] = f_buff;
                    break;
            }
        }
        //		tokens.clear();
    }
}
/*!
   \brief Destrcutor of class LinearFunctionData

   WW 24.08.2011
 */
LinearFunctionData::~LinearFunctionData()
{
    delete[] _a0;
    delete[] _b0;
    delete[] _c0;
    delete[] _d0;
    if (_subdom_index)
        delete[] _subdom_index;

    _subdom_index = NULL;
    _a0 = _b0 = _c0 = _d0 = NULL;
}

/*!
   \brief Get Value of class LinearFunctionData

   WW 24.08.2011
 */
double LinearFunctionData::getValue(size_t dom_i, double x, double y,
                                    double z) const
{
    for (size_t i = 0; i < _ndata; i++)
        if (dom_i == _subdom_index[i])
            return _a0[i] + _b0[i] * x + _c0[i] * y + _d0[i] * z;

    return 0.;
}
/*!
   \brief Get Value of class LinearFunctionData

   WW 24.08.2011
 */
double LinearFunctionData::getValue(double x, double y, double z) const
{
    return _a0[0] + _b0[0] * x + _c0[0] * y + _d0[0] * z;
}
