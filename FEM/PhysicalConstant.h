/**
 * \file PhysicalConstant.h
 *  Define physical constants
 *
 * \author Christoph Lehmann
 * \author Wenqing Wang
 *
 * \date 04-2015, 07-2015
 *
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef OGS_PHYSICAL_CONSTANT_H
#define OGS_PHYSICAL_CONSTANT_H

/**
 * Namespace containing physical constants
 * All members of this namespace should be given in SI standard units,
 * i.e. in terms of kg, m, s, K, mol, A, cd.
 */
namespace PhysicalConstant
{
/// Celsius zero in Kelvin
const double CelsiusZeroInKelvin = 273.15;

/**
  Ideal gas constant in SI standard units (J mol^-1 K^-1)

  Source:               http://physics.nist.gov/cgi-bin/cuu/Value?r
  Date visited:         2015-04-17
  Standard uncertainty: 0.000 0075 J mol^-1 K^-1
*/
const double IdealGasConstant = 8.3144621;

/**
 * Molar masses of certain elements and chemical compounds
 */
namespace MolarMass
{
/**
 * Water
 *
 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of O is in the range [15.999 03, 15.999 77] g/mol
 * and the molar mass of H is in the range [1.007 84, 1.008 11] g/mol
 */
const double Water = 0.018016; ///< kg mol^-1

/**
 * N_2
 *
 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of N is in the range [14.006 43, 14.007 28] g/mol
 */
const double N2 = 0.028013; ///< kg mol^-1

/**
 * O_2
 *
 * Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of O is in the range [15.999 03, 15.999 77] g/mol
 */
const double O2 = 0.032; ///< kg mol^-1

/**
 * Air
 *
 * Source:       http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
 *
 */
const double Air = 0.02897; ///< kg mol^-1
}

/**
  Specific gas constant (J/(kg K)) defined by
  \f[
  R_v = \frac{R}{M}
  \f]
  with \$M\$, the momlar mass
*/
namespace SpecificGasConstant
{
const double WaterVapour = IdealGasConstant / MolarMass::Water; //=461.504;
}
}
#endif // OGS_CONSTANTS_H
