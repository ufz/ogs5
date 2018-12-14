/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>

// using namespace std;

class PITZdata
{
private:
public:
    PITZdata(void);
    ~PITZdata(void);
    // method
    static double charge(std::string N);
    static double pitzer_parameters(double T, double P,
                                    std::string param_switch);
};
