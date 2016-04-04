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
   09/2008 NB Function with multiple arguments
**************************************************************************/
#include "files0.h"
#include "makros.h"
#include "rf_fct.h"

#include <cstdlib>

using namespace std;

extern void remove_white_space(std::string* buffer);
extern string GetLineFromFile1(std::ifstream*);

std::vector<CFunction*> fct_vector;

/**************************************************************************
   FEMLib-Method: CFunction
   Task: constructor
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
CFunction::CFunction(void)
{
	selected = true;
}

/**************************************************************************
   FEMLib-Method: CFunction
   Task: destructor
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
CFunction::~CFunction(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: get instance by name
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
CFunction* FCTGet(std::string fct_name)
{
	int i;
	int fct_vector_size = (int)fct_vector.size();
	CFunction* m_fct = NULL;
	for (i = 0; i < fct_vector_size; i++)
	{
		m_fct = fct_vector[i];
		if (fct_name.compare(m_fct->type_name) == 0) // OK
			return m_fct;
	}
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: read functionn
   Programing:
   03/2005 OK Implementation
   09/2008 NB, matrix functions
   last modification:
**************************************************************************/
std::ios::pos_type CFunction::Read(std::ifstream* fct_file)
{
	bool ok_true = true;
	bool new_keyword = false;
	// WW bool dim = false;
	int unsigned no_variable_names;
	// OK411 int unsigned i;
	int i;
	double* variable_data = NULL;
	double* j = 0;
	double* k = 0;
	double* l = 0;
	char buffer[MAX_ZEILE];
	string hash("#");
	string line_string;
	string test_string;
	string dollar("$");
	ios::pos_type position;
	std::stringstream line_stream;
	std::vector<double*> dim_x; // OK memory release check
	std::vector<double*> dim_y;
	std::vector<double*> dim_z;
	//----------------------------------------------------------------------
	while (!new_keyword)
	{
		position = fct_file->tellg();
		fct_file->getline(buffer, MAX_ZEILE);
		line_string = buffer;
		//....................................................................
		// Test next keyword
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// TYPE
		if (line_string.find("$TYPE") != string::npos) // subkeyword found
		{
			line_stream.str(GetLineFromFile1(fct_file));
			line_stream >> type_name;
			line_stream.clear();
		}
		//....................................................................
		// GEO_TYPE
		// subkeyword found
		if (line_string.find("$GEO_TYPE") != string::npos)
		{
			line_stream.str(GetLineFromFile1(fct_file));
			line_stream >> geo_type_name;
			if (geo_type_name.find("POINT") != string::npos)
			{
				line_stream >> geo_name;
				line_stream.clear();
			}
		}
		//....................................................................
		// DIS_TYPE //CMCD
		// subkeyword found
		if (line_string.find("$DIS_TYPE") != string::npos)
		{
			line_stream.str(GetLineFromFile1(fct_file));
			line_stream >> dis_type_name; // TRANSIENT
			line_stream.clear();
		}
		//....................................................................
		// VARIABLES
		// subkeyword found
		if (line_string.find("$VARIABLES") != string::npos)
		{
			line_stream.str(GetLineFromFile1(fct_file));
			while (ok_true)
			{
				line_stream >> line_string;
				if (line_string.compare(test_string) == 0)
				{
					line_stream.clear();
					break;
				}
				variable_names_vector.push_back(line_string);
				test_string = line_string;
			}
		}

		//....................................................................
		// MATRIX_DIMENSION     NB
		//... reads the number of colums and rows of the given matrix ........
		// subkeyword found
		if (line_string.find("$DIMENSION") != string::npos)
		{
			line_stream.str(GetLineFromFile1(fct_file));

			line_stream >> i;
			matrix_dimension.push_back(i);
			line_stream >> i;
			matrix_dimension.push_back(i);

			line_stream.clear();
			line_string = "";
		}
		//....................................................................
		// MATRIX       NB
		//...reads both arguments and depending function values in the given..
		//...matrix and saves them consecutively in variable_data_vector......
		// subkeyword found
		if (line_string.find("$MATRIX") != string::npos)
		{
			for (i = 0; i < matrix_dimension[0] * matrix_dimension[1]; i++)
			{
				line_stream.str(GetLineFromFile1(fct_file));
				line_string = line_stream.str();

				j = new double; // OK
				k = new double;
				l = new double;

				line_stream >> *j >> *k >> *l;
				if ((int)dim_x.size() < matrix_dimension[0])
					dim_x.push_back(j);
				if ((dim_y.size() > 0) && ((*k != *dim_y[dim_y.size() - 1])))
					dim_y.push_back(k);
				else if (dim_y.size() == 0)
					dim_y.push_back(k);
				dim_z.push_back(l);
				line_stream.clear();
			}

			// OK411
			for (i = 0; i < (int)dim_x.size(); i++)
				variable_data_vector.push_back(dim_x[i]);
			for (i = 0; i < (int)dim_y.size(); i++)
				variable_data_vector.push_back(dim_y[i]);
			for (i = 0; i < (int)dim_z.size(); i++)
				variable_data_vector.push_back(dim_z[i]);

			// OK411
			if ((int)variable_data_vector.size()
			    != matrix_dimension[0] * matrix_dimension[1] + matrix_dimension[0] + matrix_dimension[1])
			{
				std::cout << "FCT function: Error! The number of data does not correspond to the specified matrix size"
				          << "\n";
				break;
			}
			line_stream.clear();
			line_string = "";
		}

		//--------------------------------------------------------------------
		// DATA
		if (line_string.find("$DATA") != string::npos) // subkeyword found
		{
			no_variable_names = (int)variable_names_vector.size();
			while (ok_true)
			{
				line_string = GetLineFromFile1(fct_file);
				line_stream.str(line_string);
				if (line_string.find("$") != string::npos)
				{
					line_stream.clear();
					break;
				}
				if (line_string.find("#") != string::npos)
				{
					line_stream.clear();
					return position;
				}
				variable_data = new double[no_variable_names];
				for (i = 0; i < (int)no_variable_names; i++) // OK411
				{
					line_stream >> variable_data[i];
					variable_data_vector.push_back(variable_data);
				}
				line_stream.clear();
			}
		}
		//--------------------------------------------------------------------
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
void CFunction::Write(std::fstream* fct_file)
{
	int i, j;
	int variable_names_vector_size = (int)variable_names_vector.size();
	int variable_data_vector_size = (int)variable_data_vector.size();
	//--------------------------------------------------------------------
	// CHECK
	//--------------------------------------------------------------------
	// KEYWORD
	*fct_file << "#FUNCTION"
	          << "\n";
	//--------------------------------------------------------------------
	// TYPE
	*fct_file << " $TYPE"
	          << "\n";
	*fct_file << "  ";
	*fct_file << type_name << "\n";
	//--------------------------------------------------------------------
	// GEO_TYPE
	*fct_file << " $GEO_TYPE"
	          << "\n";
	*fct_file << "  ";
	*fct_file << geo_type_name << " " << geo_name << "\n";
	//--------------------------------------------------------------------
	// DATA
	*fct_file << " $VARIABLES"
	          << "\n";
	*fct_file << " ";
	for (i = 0; i < variable_names_vector_size; i++)
		*fct_file << " " << variable_names_vector[i];
	*fct_file << "\n";
	//--------------------------------------------------------------------
	// DATA
	*fct_file << " $DATA"
	          << "\n";
	for (i = 0; i < variable_data_vector_size; i++)
	{
		*fct_file << " ";
		for (j = 0; j < variable_names_vector_size; j++)
			*fct_file << " " << variable_data_vector[i][j];
		*fct_file << "\n";
	}
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task: master read functionn
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
void FCTRead(std::string base_file_name)
{
	CFunction* m_fct = NULL;
	char line[MAX_ZEILE];
	string sub_line;
	string line_string;
	ios::pos_type position;
	//========================================================================
	// file handling
	std::string fct_file_name;
	fct_file_name = base_file_name + FCT_FILE_EXTENSION;
	std::ifstream fct_file(fct_file_name.data(), std::ios::in);
	if (!fct_file.good())
		return;
	fct_file.seekg(0L, std::ios::beg);
	//========================================================================
	// keyword loop
	std::cout << "FCTRead"
	          << "\n";
	while (!fct_file.eof())
	{
		fct_file.getline(line, MAX_ZEILE);
		line_string = line;

		if (line_string.find("#STOP") != std::string::npos)
		{
			return;
		}

		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#FUNCTION") != std::string::npos)
		{
			m_fct = new CFunction();
			position = m_fct->Read(&fct_file);
			fct_vector.push_back(m_fct);
			fct_file.seekg(position, std::ios::beg);
		} // keyword found
	} // eof
}

/**************************************************************************
   FEMLib-Method:
   Task: master write functionn
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
void FCTWrite(std::string file_base_name)
// void MATWriteMediumProperties(fstream *mp_file)
{
	CFunction* m_fct = NULL;
	std::string sub_line;
	std::string line_string;
	//========================================================================
	// File handling
	std::string fct_file_name = file_base_name + FCT_FILE_EXTENSION;
	std::fstream fct_file(fct_file_name.c_str(), ios::trunc | ios::out);
	fct_file.setf(std::ios::scientific, std::ios::floatfield);
	fct_file.precision(12);
	if (!fct_file.good())
		return;
	fct_file << "GeoSys-FCT: Functions ------------------------------------------------"
	         << "\n";
	//========================================================================
	// FCT vect list
	int no_fct = (int)fct_vector.size();
	int i;
	for (i = 0; i < no_fct; i++)
	{
		m_fct = fct_vector[i];
		m_fct->Write(&fct_file);
	}
	fct_file << "#STOP";
	// mp_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   04/2005 OK MODE
**************************************************************************/
void FCTReadTIMData(std::string file_name_base)
{
	int counter;
	int pos1, pos2;
	double scale_factor = 1.;
	double* variable_data = NULL;
	char line[MAX_ZEILE];
	std::string sub_string;
	std::string line_string;
	std::string delimiter_type(";");
	std::string fct_type_name;
	std::string fct_x_name;
	std::string fct_x_unit;
	std::string fct_y_name;
	std::string fct_y_unit;
	CFunction* m_fct = NULL;
	std::string header_rows;
	std::string header_cols;
	int mode = -1;
	int i;
	double time;
	double value;
	//----------------------------------------------------------------------
	std::cout << "Read FCT properties from " << file_name_base << "\n";
	//----------------------------------------------------------------------
	// File handling
	// OK4105
	std::string csv_file_name = file_name_base + CSV_FILE_EXTENSION;
	std::ifstream csv_file(csv_file_name.data(), std::ios::in);
	if (!csv_file.good())
		return;
	csv_file.seekg(0L, std::ios::beg);
	//----------------------------------------------------------------------
	// Evaluate header
	csv_file.getline(line, MAX_ZEILE);
	line_string = line;
	//......................................................................
	// MODE check
	pos1 = 0;
	sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
	sub_string = line_string.substr(pos1, pos2 - pos1);
	header_rows = sub_string;
	pos1 = pos2 + 1;
	sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
	sub_string = line_string.substr(pos1, pos2 - pos1);
	header_cols = sub_string;
	if (header_cols.find("TIME") != std::string::npos)
		mode = 0;
	if (header_rows.find("TIME") != std::string::npos)
		mode = 1;
	//======================================================================
	// MODE 0 - ROW: PNT, COL: TIM
	if (mode == 0)
	{
		pos1 = 0;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		fct_type_name = sub_string;
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, "=", pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		fct_x_name = sub_string;
		pos1 = pos2 + 1;
		sub_string = line_string.substr(pos1, pos2 + fct_x_name.size() + 2 - pos1);
		fct_x_unit = sub_string;
		if (fct_x_unit.find("MONTH") != string::npos)
			scale_factor = 86400. * 30.;
		pos1 = pos2 + (int)fct_x_name.size() + 2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		fct_y_name = sub_string;
		//----------------------------------------------------------------------
		// Read header
		csv_file.getline(line, MAX_ZEILE);
		//----------------------------------------------------------------------
		// Read POINT, TIME ...
		while (!csv_file.eof())
		{
			csv_file.getline(line, MAX_ZEILE);
			line_string = line;
			if (line_string.empty())
				return;
			//....................................................................
			m_fct = new CFunction();
			m_fct->type_name = fct_type_name;
			m_fct->variable_names_vector.push_back(fct_x_name);
			m_fct->variable_names_vector.push_back(fct_y_name);
			//....................................................................
			// POINT
			pos1 = 0;
			sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
			sub_string = line_string.substr(pos1, pos2 - pos1);
			m_fct->geo_type_name = "POINT";
			m_fct->geo_name = sub_string;
			fct_vector.push_back(m_fct);
			//....................................................................
			// T = 0
			variable_data = new double[2];
			variable_data[0] = 0.0;
			variable_data[1] = 0.0;
			m_fct->variable_data_vector.push_back(variable_data);
			//....................................................................
			// TIME, VALUE
			counter = 0;
			while ((pos1 > -1) && (pos2 > -1))
			{
				counter++;
				pos1 = pos2 + 1;
				sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
				if (pos2 < 0)
					break;
				variable_data = new double[2];
				// TIME
				variable_data[0] = counter * scale_factor;
				// VALUE
				sub_string = line_string.substr(pos1, pos2 - pos1);
				variable_data[1] = strtod(sub_string.data(), NULL);
				m_fct->variable_data_vector.push_back(variable_data);
			}
			//....................................................................
		} // eof
		//----------------------------------------------------------------------
	} // MODE==0
	//======================================================================
	// MODE 1 - ROW: TIM, COL: PNT
	if (mode == 1)
	{
		csv_file.seekg(0L, std::ios::beg);
		csv_file.getline(line, MAX_ZEILE);
		pos1 = 0;
		sub_string = get_sub_string(header_rows, "=", pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		fct_x_name = sub_string;
		pos1 = pos2 + 1;
		sub_string = line_string.substr(pos1, pos2 + fct_x_name.size() + 1 - pos1);
		fct_x_unit = sub_string;
		if (fct_x_unit.find("MONTH") != string::npos)
			scale_factor = 86400. * 30.;
		if (fct_x_unit.find("YEAR") != string::npos)
			scale_factor = 86400. * 364.25;
		pos1 = pos2 + (int)fct_x_name.size() + 1 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		fct_y_name = sub_string;
		//--------------------------------------------------------------------
		// Read header
		csv_file.getline(line, MAX_ZEILE);
		line_string = line;
		pos1 = 0;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		// sub_string = line_string.substr(pos1,pos2-pos1);
		int no_fct = 0;
		while (pos2 > 0)
		{
			m_fct = new CFunction();
			m_fct->geo_type_name = "POINT";
			m_fct->type_name = file_name_base;
			m_fct->variable_names_vector.push_back(fct_x_name);
			m_fct->variable_names_vector.push_back(fct_y_name);
			pos1 = pos2 + 1;
			sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
			sub_string = line_string.substr(pos1, pos2 - pos1);
			m_fct->geo_name = sub_string;
			fct_vector.push_back(m_fct);
			no_fct++;
		}
		//--------------------------------------------------------------------
		// Read data: TIME, POINT ...
		while (!csv_file.eof())
		{
			csv_file.getline(line, MAX_ZEILE);
			line_string = line;
			if (line_string.empty())
				return;
			//..................................................................
			// TIM
			pos1 = 0;
			sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
			sub_string = line_string.substr(pos1, pos2 - pos1);
			time = strtod(sub_string.data(), NULL);
			for (i = 0; i < no_fct; i++)
			{
				m_fct = fct_vector[i];
				variable_data = new double[2];
				pos1 = pos2 + 1;
				sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
				sub_string = line_string.substr(pos1, pos2 - pos1);
				variable_data[0] = time * scale_factor;
				value = strtod(sub_string.data(), NULL);
				variable_data[1] = value;
				m_fct->variable_data_vector.push_back(variable_data);
			}
			//..................................................................
		} // eof
		//--------------------------------------------------------------------
	}
	//======================================================================
}

/**************************************************************************
   FEMLib-Method:
   Task: read functionn
   Programing:

   01/2005 OK Implementation
   04/2006 CMCD Transient function
   last modification:
**************************************************************************/
double CFunction::GetValue(double point, bool* valid, int method)
{
	long anz;
	anz = (long)variable_data_vector.size();
	if (anz == 0)
		return 1.0;
	*valid = true;
	register long i;
	long j;
	double value = 1.0;
	//----------------------------------------------------------------------
	//
	value = variable_data_vector[0][1];
	if (point < variable_data_vector[0][0])
	{
		*valid = false;
		return value;
	}
	if (point > variable_data_vector[anz - 1][0])
	{
		*valid = false;
		return variable_data_vector[anz - 1][1];
	}
	// Check what type of function, continuous or transient
	if (dis_type_name != "TRANSIENT")
	{
		//----------------------------------------------------------------------
		// Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend geordnet
		i = 1l;
		while (point > variable_data_vector[i][0])
		{
			i++;
			switch (method)
			{
				case 0: // Lineare Interpolation
					value = variable_data_vector[i - 1][1]
					        + (variable_data_vector[i][1] - variable_data_vector[i - 1][1])
					              / (variable_data_vector[i][0] - variable_data_vector[i - 1][0])
					              * (point - variable_data_vector[i - 1][0]);
					break;
				case 1: // Piece wise constant. FS//WW
					value = variable_data_vector[i - 1][1];
					break;
			}
		}
	}
	if (dis_type_name == "TRANSIENT")
	{
		long index = 0;
		long pp = 0;
		long np = 0;
		for (j = 0; j < anz; j += 2)
		{
			double temp1 = variable_data_vector[j][0];
			// OK      double temp2 = variable_data_vector[j][1];
			if (point > temp1)
				index++;
		}
		if (index % 2 == 1)
		{
			pp = (index - 1) * 2;
			np = index * 2;
			value = variable_data_vector[pp][1]
			        + (variable_data_vector[np][1] - variable_data_vector[pp][1])
			              / (variable_data_vector[np][0] - variable_data_vector[pp][0])
			              * (point - variable_data_vector[pp][0]);
			*valid = true;
		}
		if (index % 2 == 0)
		{
			*valid = false;
			return -1.0;
		}
	}
	return value;
}

/**************************************************************************
   FEMLib-Method:
   Task: get instance by name
   Programing:
   08/2006 YD Implementation
   last modification:
**************************************************************************/
CFunction* FCTGet(long fct_number)
{
	CFunction* m_fct = NULL;
	m_fct = fct_vector[fct_number];
	return m_fct;
}
