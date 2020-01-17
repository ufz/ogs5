/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

class NR
{
private:
public:
    NR(void);
    ~NR(void);

    static double dfridr(double func(double), const double /*x*/);
    static double dfridrX(double func(double, double), const double /*x*/,
                          const double /*y*/);
    static double dfridrY(double func(double, double), const double /*x*/,
                          const double /*y*/);
};
