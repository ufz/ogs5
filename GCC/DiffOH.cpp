/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/* References
(1).Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A., and
Helgeson, H. C., 1992. Calculation of the Thermodynamic Properties of Aqueous
Species at High-Pressures and Temperatures- Effective Electrostatic Radii,
Dissociation-Constants and Standard Partial Molal Properties to 1000-Degrees-C
and 5-Kbar. Journal of the Chemical Society-Faraday Transactions 88, 803-826.

(2).Oelkers, E. H. and Helgeson, H. C., 1988. Calculation of the Thermodynamic
and Transport- Properties of Aqueous Species at High-Pressures and Temperatures
- Aqueous Tracer Diffusion- Coefficients of Ions to 1000-Degrees-C and 5-Kb.
    Geochimica Et Cosmochimica Acta 52, 63-85.

(4).Helgeson, H. C., 1992. Effects of complex formation in flowing fluids on the
hydrothermal solubilities of minerals as a function of fluid pressure and
temperature in the critical and supercritical regions of the system H2O.
    Geochimica Et Cosmochimica Acta 56, 3191-3207.
*/

#include <limits>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <vector>
#include "DiffOH.h"
#include "IAPWS-IF97.h"
#include "VLE.h"
#include "NR.h"
#include "IO.h"

using namespace std;

Diffusion::Diffusion(void)
{
    index1 = 0;
}

Diffusion::~Diffusion(void) {}
