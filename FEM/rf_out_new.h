/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: OUT
   Task: class implementation
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_out_new_INC
#define rf_out_new_INC

class COutput;
namespace GEOLIB
{
class GEOObjects;
}

extern std::vector<COutput*> out_vector;

extern std::string defaultOutputPath;

/**
 * read file that stores information about output
 * @param file_base_name base file name (without extension)
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if file reading was successful, else false
 */
bool OUTRead(const std::string& file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

extern void OUTWrite(std::string);
#define OUT_FILE_EXTENSION ".out"
extern void OUTData(double, const int step, bool force_output);
extern void OUTDelete();
extern COutput* OUTGet(const std::string&);
extern void OUTCheck(void); // new SB
extern COutput* OUTGetRWPT(const std::string&); // JT
#endif
