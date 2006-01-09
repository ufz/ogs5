/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib-Object: FCT
   Task: Functions
   Programing:
   01/2005 OK Implementation
**************************************************************************/
#ifndef rf_fct_INC
#define rf_fct_INC
// C++ STL
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

class CFunction
{
private:
public:
	CFunction(void);
	~CFunction(void);
	void Set(std::string,std::string,double);
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::fstream*);
	// add int method. WW
	double GetValue(double point, bool* valid, int method = 0);
public:
	std::string type_name;
	std::string geo_name;
	std::string geo_type_name;
	std::string dis_type_name;            //CMCD
	std::vector<std::string>variable_names_vector;
	std::vector<double*>variable_data_vector;
	//	int matrix_dimension_x; //NB4703
	//	int matrix_dimension_y; //NB4703
	std::vector<int>matrix_dimension;     //NB 4.8.01
	bool selected;
};

extern std::vector<CFunction*>fct_vector;
extern void FCTWrite(std::string);
extern void FCTRead(std::string);
extern void FCTReadTIMData(std::string);
extern CFunction* FCTGet(std::string);
extern CFunction* FCTGet(long);                   //YD

#define FCT_FILE_EXTENSION ".fct"
#endif
