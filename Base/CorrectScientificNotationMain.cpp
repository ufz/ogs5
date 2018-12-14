/**
 * 28/5/2010 LB Initial implementation
 *
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "StringTools.h"

int main(int argc, const char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: correctScientificNotation filename" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    correctScientificNotation(filename);
}
