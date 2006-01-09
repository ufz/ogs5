/*
 * BoundaryConditionIO.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// Base
#include "readNonBlankLineFromInputStream.h"

// FileIO
#include "BoundaryConditionIO.h"
#include "GeoIO.h"
#include "ProcessIO.h"

namespace FileIO
{
//CBoundaryCondition* BoundaryConditionIO::read(std::istream& in_str,
//		GEOLIB::GEOObjects const& geo_obj, std::string const& unique_fname)
//{
//	ProcessType pcs_type (INVALID_PROCESS);
//	GeoInfo *geo_info (new GeoInfo);
//	std::string geo_name;
//
//	bool new_keyword = false;
//	// loop over components of boundary conditions
//	while (!new_keyword) {
//		std::string line_string (readNonBlankLineFromInputStream (in_str));
//		if (line_string.empty()) break;
//
//		std::string sub_string, strbuff;
//		int ibuff; //pos,
//		double dbuff; //WW
//		std::stringstream in;
//
//		if (line_string.find("#") != std::string::npos) {
//			new_keyword = true;
//			break;
//		}
//
//		if (line_string.find("$PCS_TYPE") != std::string::npos) {
//			ProcessIO::readProcessInfo (in_str, pcs_type);
//		}
//
////		if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			std::string tmp;
////			in >> tmp; // _pcs_pv_name;
////			if (this->_pcs_type == MASS_TRANSPORT) {
////				// HS set the pointer to MCP based on component name.
////				// a check whether this name is existing and unique.
////				if (cp_name_2_idx.count(tmp) == 1) {
////					setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
////					setProcessPrimaryVariable(CONCENTRATION);
////				} else {
////					DisplayErrorMsg(
////							"Error: In reading BC file, the input component names are not found in MCP file!!!");
////					exit(1);
////				}
////			} else {
////				setProcess(PCSGet(this->getProcessType()));
////				setProcessPrimaryVariable(convertPrimaryVariable(tmp));
////			}
////			in.clear();
////		}
//
////		// HS, this is new. later on we should stick to COMP_NAME, PRIMARY_VARIABLE support will be removed.
////		if (line_string.find("$COMP_NAME") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			std::string tmp;
////			in >> tmp; // _pcs_pv_name;
////			if (this->_pcs_type == MASS_TRANSPORT) {
////				// HS set the pointer to MCP based on component name.
////				// check whether this name is existing and unique.
////				if (cp_name_2_idx.count(tmp) == 1) {
////					setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
////					setProcessPrimaryVariable(CONCENTRATION);
////				} else {
////					DisplayErrorMsg(
////							"Error: In reading BC file, the input component names are not found in MCP file!!!");
////					exit(1);
////				}
////			}
////			in.clear();
////		}
//
//		//subkeyword found
//		if (line_string.find("$GEO_TYPE") != std::string::npos) {
//			GeoIO::readGeoInfo (geo_info, in_str, geo_name, geo_obj, unique_fname);
//		}
//
////		//PCH
////		if (line_string.find("$DIS_TYPE") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			in >> line_string; //sub_line
////			_periodic = false; // JOD
////
////			// Soure terms are assign to element nodes directly. 23.02.2009. WW
////			if (line_string.find("DIRECT") != std::string::npos) {
////				this->setProcessDistributionType(FiniteElement::DIRECT);
////				in >> fname;
////				fname = FilePath + fname;
////				in.clear();
////			}
////			if (line_string.find("CONSTANT") != std::string::npos) {
////				this->setProcessDistributionType(FiniteElement::CONSTANT);
////				in >> geo_node_value; //sub_line
////				in.clear();
////			}
////			if (line_string.find("LINEAR") != std::string::npos) {
////				this->setProcessDistributionType(FiniteElement::LINEAR);
////				// Distribued. WW
////				size_t nLBC;
////				in >> nLBC; //sub_line
////				in.clear();
////
////				//        sub_string = strtok(buffer,seps);
////				//        sub_string = strtok( NULL, seps );
////				//        int nLBC = atoi(sub_string.c_str());
////				for (size_t i = 0; i < nLBC; i++) {
////					in.str(GetLineFromFile1(in_str));
////					in >> ibuff >> dbuff >> strbuff;
////					in.clear();
////
////					//           *in_str>>ibuff>>dbuff;
////					_PointsHaveDistribedBC.push_back(ibuff);
////					_DistribedBC.push_back(dbuff);
////					if (strbuff.size() > 0) {
////						_PointsFCTNames.push_back(strbuff);
////						time_dep_interpol = true;
////					}
////				}
////				//        in_str->ignore(MAX_ZEILE,'\n');
////			}
////		}
//
//		// Time dependent function
//		//..Time dependent curve ............................................
//		// subkeyword found
////		if (line_string.find("$TIM_TYPE") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			in >> line_string;
////
////			if (line_string.find("CURVE") != std::string::npos) {
////				this->setProcessDistributionType(FiniteElement::CONSTANT);
////				in >> _curve_index;
////				in.clear();
////			}
////			continue;
////		}
//
//		// subkeyword found
////		if (line_string.find("$FCT_TYPE") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			in >> fct_name; //sub_line
////			in.clear();
////		}
//
//		//subkeyword found
////		if (line_string.find("$MSH_TYPE") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			in >> sub_string; //sub_line
////			_msh_type_name = "NODE";
////			if (sub_string.find("NODE") != std::string::npos) {
////				in >> _msh_node_number;
////				in.clear();
////			}
////		}
//
//		// subkeyword found
////		if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str)); // CONSTANT -21500.0
////			in >> line_string;
////			if (line_string.find("CONSTANT") != std::string::npos) {
////				this->setProcessDistributionType(FiniteElement::CONSTANT);
////				in >> geo_node_value;
////				in.clear();
////			}
////			in.str(GetLineFromFile1(in_str)); // 0.0 IF HEAD > 0.04
////			std::string pcs_pv_name_cond; // 07/2010 TF temp string
////			in >> node_value_cond >> line_string >> pcs_pv_name_cond
////					>> line_string >> condition;
////			in.clear();
////			in.str(GetLineFromFile1(in_str)); // PCS OVERLAND_FLOW
////			std::string pcs_type_name_cond;
////			in >> line_string >> pcs_type_name_cond;
////			in.clear();
////			conditional = true;
////		}
////		// NW
////		if (line_string.find("$EPSILON") != std::string::npos) {
////			in.str(GetLineFromFile1(in_str));
////			in >> epsilon;
////			in.clear();
////		}
//	}
//
////	this->setProcessType(convertProcessType(tmp));
//
//	return NULL;
//}

void BoundaryConditionIO::write(std::ostream& out,
                                CBoundaryCondition const& bc)
{
	// keyword
	out << "#BOUNDARY_CONDITION" << "\n";

	// process and primary variable
	out << "\t$PCS_TYPE" << "\n";
	out << "\t\t" << convertProcessTypeToString(bc.getProcessType())
	    << "\n";

	out << "\t$PRIMARY_VARIABLE" << "\n";
	out << "\t\t" << convertPrimaryVariableToString(
	        bc.getProcessPrimaryVariable()) << "\n";

	// geometry
	out << "\t$GEO_TYPE" << "\n";
	out << "\t" << bc.getGeoTypeAsString() << " " << bc.geo_name << "\n";

	// distribution type
	out << "\t$DIS_TYPE" << "\n";
	out << "\t\t" << convertDisTypeToString(bc.getProcessDistributionType());
	if (bc.getProcessDistributionType() == FiniteElement::CONSTANT)
		out << "\t\t" << bc.geo_node_value << "\n";
	else if (bc.getProcessDistributionType() == FiniteElement::LINEAR)
	{
		out << "\t\t" << bc._PointsHaveDistribedBC.size() << "\n";
		for (size_t i = 0; i < bc._PointsHaveDistribedBC.size(); i++)
			out << "\t\t" << bc._PointsHaveDistribedBC[i] << "  "
			    << bc._DistribedBC[i] << "\n";
	}
	else
		out << "\n";

	// function name
	if (!bc.fct_name.empty())
	{
		out << "\t$FCT_TYPE" << "\n";
		out << "\t\t" << bc.fct_name << "\n";
	}
}
} // end namespace FileIO
