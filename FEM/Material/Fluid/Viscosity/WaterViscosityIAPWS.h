/**
 *  \brief Viscosity model according to
 *         Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
 *         Water Substance
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterViscosityIAPWS.h
 *
 * Created on December 1, 2016, 1:41 PM
 */

#pragma once

#include <string>

namespace MaterialLib
{
namespace Fluid
{
/**
 * \brief A class for viscosity model that is defined by
 *        The International Association for the Properties of Water and Steam
 *        <a href="http://www.iapws.org/relguide/visc.pdf">IAPWS</a>
 *
 *        With the definition, the viscosity is a function of temperature and
 *        water density
 *
 *  \attention The critical enhancement, \f$\bar{\mu}_2\f$, which is significant
 *             only within the boundaries specified by
 *                 \f[ T (\mbox{in K}) \in (645.91, 650.77) \f]
 *             and
 *                 \f[ \rho (\mbox{in kg m}^{-3}) \in (245.8, 405.3)\f],
 *             is not considered.
 */
class WaterViscosityIAPWS
{
public:
	WaterViscosityIAPWS() {}

    /// Get model name.
    std::string getName()
    {
        return "IAPWS temperature-density dependent viscosity model";
    }

    /**  Get density value.
         \param rho   Density.
         \param T     Temperature.
    */
    static double getValue(const double T, const double rho);

    /**  Get density value.
         \param rho   Density.
         \param T     Temperature.
    */
    static double getdValuedT(const double T, const double rho);

	/**  Get density value.
         \param rho   Density.
         \param T     Temperature.
    */
    static double getdValuedRho(const double T, const double rho); 

    // Coefficients Hi and Hij are given in two static arrays in the cpp file.
};

}  // end namespace
}  // end namespace
