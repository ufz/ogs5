/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Class: BC BoundaryConditions
   Task:
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#include "makros.h"

// C++ STL
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <iostream>

#include "display.h"
#include "memory.h"

// FileIO
#include "BoundaryConditionIO.h"
#include "GeoIO.h"
#include "ProcessIO.h"
#include "readNonBlankLineFromInputStream.h"

// GEOLib
//#include "geo_lib.h"
//#include "geo_sfc.h"

// GEOLIB
#include "GEOObjects.h"

// MSHLib
//#include "mshlib.h"
// FEMLib
extern void remove_white_space(std::string*);
//#include "problem.h"
#include "tools.h"
//#include "rf_node.h"
#include "rf_bc_new.h"
//#include "rf_pcs.h"
//#include "rf_fct.h"
#include "rfmat_cp.h"
//#include "geo_ply.h"
// MathLib
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "mathlib.h"

#include "BoundaryCondition.h"

#ifndef _WIN32
#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

double cputime(double x)
{
    struct rusage rsrc;
    double usr, sys;

    if (getrusage(RUSAGE_SELF, &rsrc) == -1)
    {
        perror("times");
        exit(1);
    }

    usr = rsrc.ru_utime.tv_sec + 1.0e-6 * rsrc.ru_utime.tv_usec;
    sys = rsrc.ru_stime.tv_sec + 1.0e-6 * rsrc.ru_stime.tv_usec;

    return usr + sys - x;
}
#endif

using namespace Display;

CBoundaryConditionNode::CBoundaryConditionNode() : node_value_offset(0.0)
{
    conditional = false;
    for (std::size_t i = 0; i < 3; i++)
        _normal_vector[i] = 0;
}

void CBoundaryConditionNode::SetNormalVector(double const* const normal_vector)
{
    _normal_vector[0] = normal_vector[0];
    _normal_vector[1] = normal_vector[1];
    _normal_vector[2] = normal_vector[2];
}

double const* CBoundaryConditionNode::GetNormalVector() const
{
    return this->_normal_vector;
}

/**************************************************************************
   FEMLib-Method:
   Task: destructor
   Programing:
   08/2011 WW Implementation
**************************************************************************/
void CBoundaryConditionNode::Read(std::istream& is)
{
    is >> geo_node_number;
    is >> msh_node_number;
    is >> CurveIndex;
    is >> node_value;
    is >> std::ws;
}
/**************************************************************************
   FEMLib-Method:
   Task: destructor
   Programing:
   08/2011 WW Implementation
**************************************************************************/
void CBoundaryConditionNode::Write(std::ostream& os) const
{
    std::string deli = "  ";
    os << geo_node_number << deli;
    os << msh_node_number << deli;
    os << CurveIndex << deli;
    os << node_value << deli;
    os << "\n";
}

//==========================================================================
std::list<CBoundaryCondition*> bc_list;
std::vector<std::string> bc_db_head;
std::list<CBoundaryConditionsGroup*> bc_group_list;
std::vector<CBoundaryCondition*> bc_db_vector;

/**************************************************************************
   FEMLib-Method:
   Task: BC constructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::CBoundaryCondition()
    : GeoInfo(),
      geo_name(""),
      _curve_index(-1),
      dis_linear_f(NULL),
      _time_period(NULL)
{
    this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
    // FCT
    conditional = false;
    time_dep_interpol = false;
    epsilon = -1;                     // NW
    time_contr_curve = -1;            // WX
    bcExcav = -1;                     // WX
    MatGr = -1;                       // WX
    NoDispIncre = -1;                 // WX:12.2012
    gradient_ref_depth = 0;           // CB
    gradient_ref_depth_value = 0;     // CB
    gradient_ref_depth_gradient = 0;  // CB
    _pressure_as_head_model = -1;
    _pressure_as_head_density = 0;
    _isConstrainedBC = false;
    _isSeepageBC = false;
    _isSwitchBC = false;
}

// KR: Conversion from GUI-BC-object to CBoundaryCondition
CBoundaryCondition::CBoundaryCondition(const BoundaryCondition* bc)
    : ProcessInfo(bc->getProcessType(), bc->getProcessPrimaryVariable(), NULL),
      GeoInfo(bc->getGeoType(), bc->getGeoObj()),
      DistributionInfo(bc->getProcessDistributionType()),
      _time_period(NULL)
{
    setProcess(PCSGet(this->getProcessType()));
    this->geo_name = bc->getGeoName();
    const std::vector<size_t> dis_nodes = bc->getDisNodes();
    const std::vector<double> dis_values = bc->getDisValues();

    if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
    {
        this->geo_node_value = dis_values[0];
    }
    else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
    {
        for (size_t i = 0; i < dis_values.size(); i++)
        {
            this->_PointsHaveDistribedBC.push_back(
                static_cast<int>(dis_nodes[i]));
            this->_DistribedBC.push_back(dis_values[i]);
        }
    }
    else
        std::cout << "Error in CBoundaryCondition() - DistributionType \""
                  << FiniteElement::convertDisTypeToString(
                         this->getProcessDistributionType())
                  << "\" currently not supported."
                  << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task: BC deconstructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::~CBoundaryCondition()
{
    // DIS
    node_number_vector.clear();
    geo_node_number = -1;
    geo_node_value = 0.0;

    // WW
    if (dis_linear_f)
        delete dis_linear_f;
    dis_linear_f = NULL;

    if (_time_period)
        delete _time_period;
}

const std::string& CBoundaryCondition::getGeoName() const
{
    return geo_name;
}

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   09/2004 OK POINTS method
   11/2004 MX stream string
**************************************************************************/
std::ios::pos_type CBoundaryCondition::Read(std::ifstream* bc_file,
                                            const GEOLIB::GEOObjects& geo_obj,
                                            const std::string& unique_fname,
                                            bool& valid)
{
    std::string line_string;
    bool new_keyword = false;
    std::ios::pos_type position;

    std::string sub_string, strbuff;
    int ibuff;     // pos,
    double dbuff;  // WW
    std::stringstream in;

    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        position = bc_file->tellg();
        line_string = readNonBlankLineFromInputStream(*bc_file);
        if (line_string.size() < 1)
            break;
        if (line_string.find("#") != std::string::npos)
        {
            new_keyword = true;
            break;
        }

        if (line_string.find("$PCS_TYPE") != std::string::npos)
            if (!FileIO::ProcessIO::readProcessInfo(*bc_file, _pcs_type))
                valid = false;

        if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            std::string tmp;
            in >> tmp;  // _pcs_pv_name;
            if (this->_pcs_type == FiniteElement::MASS_TRANSPORT)
            {
                // HS set the pointer to MCP based on component name.
                // a check whether this name is existing and unique.
                if (cp_name_2_idx.count(tmp) == 1)
                {
                    setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
                    setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
                }
                else
                {
                    DisplayErrorMsg(
                        "Error: In reading BC file, the input component names "
                        "are not found in MCP file!!!");
                    exit(1);
                }
            }
            else
            {
                setProcess(PCSGet(this->getProcessType()));
                setProcessPrimaryVariable(
                    FiniteElement::convertPrimaryVariable(tmp));
            }
            in.clear();
        }

        // HS, this is new. later on we should stick to COMP_NAME,
        // PRIMARY_VARIABLE support will be removed.
        if (line_string.find("$COMP_NAME") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            std::string tmp;
            in >> tmp;  // _pcs_pv_name;
            if (this->_pcs_type == FiniteElement::MASS_TRANSPORT)
            {
                // HS set the pointer to MCP based on component name.
                // check whether this name is existing and unique.
                if (cp_name_2_idx.count(tmp) == 1)
                {
                    setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
                    setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
                }
                else
                {
                    DisplayErrorMsg(
                        "Error: In reading BC file, the input component names "
                        "are not found in MCP file!!!");
                    exit(1);
                }
            }
            in.clear();
        }

        if (line_string.find("$GEO_TYPE") != std::string::npos)
            if (!FileIO::GeoIO::readGeoInfo(
                    this, *bc_file, geo_name, geo_obj, unique_fname))
                valid = false;

        // PCH
        if (line_string.find("$DIS_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> line_string;  // sub_line
            _periodic = false;  // JOD

            // Source terms are assign to element nodes directly. 23.02.2009. WW
            if (line_string.find("DIRECT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::DIRECT);
                in >> fname;
                fname = FilePath + fname;
                in.clear();
            }

            if (line_string.find("CONSTANT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> geo_node_value;  // sub_line
                in.clear();
            }

            if (line_string.find("SWITCH") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::SWITCH);
                this->_isSwitchBC = true;
                std::string tempst;
                // read process
                in >> tempst;
                _switchBC.switchProcessType =
                    FiniteElement::convertProcessType(tempst);
                if (_switchBC.switchProcessType ==
                    FiniteElement::INVALID_PROCESS)
                {
                    std::cerr
                        << "Invalid Process type in SWITCH BC! Exiting now."
                        << std::endl;
                    std::exit(0);
                }
                // ... and associated primary variable
                in >> tempst;
                _switchBC.switchPrimVar =
                    FiniteElement::convertPrimaryVariable(tempst);

                // ... which will give the simulated switch value
                if (!(in >> _switchBC.switchValue))
                {
                    std::cerr
                        << "No switch value given in SWITCH BC! Exiting now."
                        << std::endl;
                    std::exit(0);
                }
                // ... to switch from the on state (ie simulated switch value is
                // higher than this one:)
                if (!(in >> _switchBC.switchOnValue))
                {
                    std::cerr
                        << "No switch on value given in SWITCH BC! Exiting now."
                        << std::endl;
                    std::exit(0);
                }
                // ... or to switch to the off state (ie simulated switch value
                // is lower than this one:)
                if (!(in >> _switchBC.switchOffValue))
                {
                    std::cerr << "No switch off value given in SWITCH BC! "
                                 "Exiting now."
                              << std::endl;
                    std::exit(0);
                }
                in.clear();
            }
            // If a linear function is given. 25.08.2011. WW
            if (line_string.find("FUNCTION") != std::string::npos)
            {
                setProcessDistributionType(FiniteElement::FUNCTION);
                in.clear();
                dis_linear_f = new LinearFunctionData(*bc_file);
            }
            if (line_string.find("LINEAR") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::LINEAR);
                // Distribuded. WW
                size_t nLBC;
                in >> nLBC;  // sub_line
                in.clear();

                for (size_t i = 0; i < nLBC; i++)
                {
                    in.str(readNonBlankLineFromInputStream(*bc_file));
                    in >> ibuff >> dbuff >> strbuff;
                    in.clear();

                    //           *bc_file>>ibuff>>dbuff;
                    _PointsHaveDistribedBC.push_back(ibuff);
                    _DistribedBC.push_back(dbuff);
                    if (strbuff.size() > 0)
                    {
                        _PointsFCTNames.push_back(strbuff);
                        time_dep_interpol = true;
                    }
                }
                //        bc_file->ignore(MAX_ZEILE,'\n');
            }
        }

        if (line_string.find("GRADIENT") != std::string::npos)  // 6/2012  JOD
        {
            this->setProcessDistributionType(FiniteElement::GRADIENT);
            in >> gradient_ref_depth;
            in >> gradient_ref_depth_value;
            in >> gradient_ref_depth_gradient;
            in.clear();
        }

        // Time dependent function
        //..Time dependent curve ............................................
        if (line_string.find("$TIM_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> line_string;

            if (line_string.find("CURVE") != std::string::npos)
            {
                //				tim_type_name = "CURVE";
                // WW this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> _curve_index;
                in.clear();

                //        pos1=pos2+1;
                //        sub_string = get_sub_string(buffer,"  ",pos1,&pos2);
                //		_curve_index = atoi(sub_string.c_str());
            }
            continue;
        }

        if (line_string.find("$TIME_PERIOD") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            double t1, t2;
            in >> t1 >> t2;
            _time_period = new TimePeriod(t1, t2);
            continue;
        }

        if (line_string.find("$FCT_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> fct_name;  // sub_line
            in.clear();
        }

        if (line_string.find("$MSH_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> sub_string;  // sub_line
            _msh_type_name = "NODE";
            if (sub_string.find("NODE") != std::string::npos)
            {
                in >> _msh_node_number;
                in.clear();
            }
        }

        if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(
                *bc_file));  // CONSTANT -21500.0
            in >> line_string;
            if (line_string.find("CONSTANT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> geo_node_value;
                in.clear();
            }
            in.str(readNonBlankLineFromInputStream(
                *bc_file));                // 0.0 IF HEAD > 0.04
            std::string pcs_pv_name_cond;  // 07/2010 TF temp string
            in >> node_value_cond >> line_string >> pcs_pv_name_cond >>
                line_string >> condition;
            in.clear();
            in.str(readNonBlankLineFromInputStream(
                *bc_file));  // PCS OVERLAND_FLOW
            std::string pcs_type_name_cond;
            in >> line_string >> pcs_type_name_cond;
            in.clear();
            conditional = true;
        }

        if (line_string.find("$EPSILON") != std::string::npos)  // NW
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> epsilon;
            in.clear();
        }
        //....................................................................
        // aktive state of the bc is time controlled  WX
        if (line_string.find("$TIME_CONTROLLED_ACTIVE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> time_contr_curve;
            in.clear();
        }
        //....................................................................
        // bc for excated boundaries WX
        if (line_string.find("$EXCAVATION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> bcExcav >> MatGr;
            in.clear();
        }
        //....................................................................
        // NO DISPLACEMENT INCREMENT WX:12.2012
        if (line_string.find("$NO_DISP_INCREMENT") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> NoDispIncre;
            in.clear();
        }
        //....................................................................
        // bc for copying of primary variables; give geometry type and name SB
        // 09.2012
        if (line_string.find("$COPY_VALUE") != std::string::npos)
        {
            // CB_merge_0513 ??
            in.str(readNonBlankLineFromInputStream(*bc_file));
            // in.str(GetLineFromFile1(bc_file));
            in >> copy_geom >> copy_geom_name;
            in.clear();
        }
        //....................................................................
        if (line_string.find("$PRESSURE_AS_HEAD") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> _pressure_as_head_model;  // 0 -> calc pressure from density;
                                            // 1 -> calc pressure from *given*
                                            // density
            if (_pressure_as_head_model == 1)
                in >> _pressure_as_head_density;
            else if (_pressure_as_head_model < 0 || _pressure_as_head_model > 1)
            {
                std::cout << "Unsupported PRESSURE_AS_HEAD model "
                          << _pressure_as_head_model << std::endl;
                _pressure_as_head_model = -1;
            }
            in.clear();
        }
        //....................................................................
        if (line_string.find("$CONSTRAINED") != std::string::npos)
        {
            Constrained temp;

            _isConstrainedBC = true;
            in.str(readNonBlankLineFromInputStream(*bc_file));
            std::string tempst, tempst2;

            in >> tempst >>
                tempst2;  // VELOCITY and DIRECTION (positive/negative scalar
                          // product between velocity vector
            // and surface normal); or PROCESS_TYPE and associated
            // PRIMARY_VARIABLE
            if (tempst == "VELOCITY")
            {
                temp.constrainedVariable = convertConstrainedVariable(tempst);
                temp.constrainedDirection = convertConstrainedType(tempst2);
                temp.constrainedValue = std::numeric_limits<size_t>::max();
                temp.constrainedPrimVar = FiniteElement::INVALID_PV;
                temp.constrainedProcessType = FiniteElement::INVALID_PROCESS;
                if (!(temp.constrainedDirection == ConstrainedType::POSITIVE ||
                      temp.constrainedDirection == ConstrainedType::NEGATIVE))
                {
                    std::cout << "No valid constrainedDirection for "
                              << convertConstrainedVariableToString(
                                     temp.constrainedVariable)
                              << "(" << tempst2 << ")" << std::endl;
                    _isConstrainedBC = false;
                }

                in >> tempst;
                if (tempst == "STABLE")
                    temp._isConstrainedVelStable = true;

                if (getGeoType() != GEOLIB::SURFACE)
                    std::cout << "\n Warning! Make sure, that a velocity "
                                 "constrained BC is a SURFACE!"
                              << std::endl;
            }
            else
            {
                temp.constrainedProcessType =
                    FiniteElement::convertProcessType(tempst);
                if (!(temp.constrainedProcessType ==
                          FiniteElement::MASS_TRANSPORT ||
                      temp.constrainedProcessType ==
                          FiniteElement::HEAT_TRANSPORT ||
                      temp.constrainedProcessType ==
                          FiniteElement::LIQUID_FLOW ||
                      temp.constrainedProcessType ==
                          FiniteElement::RICHARDS_FLOW))
                {
                    _isConstrainedBC = false;
                    break;
                }

                temp.constrainedPrimVar =
                    FiniteElement::convertPrimaryVariable(tempst2);

                in >> temp.constrainedValue >>
                    tempst;  // Constrained Value; and constrain direction
                             // (greater/smaller than value)
                temp.constrainedDirection = convertConstrainedType(tempst);
                temp.constrainedVariable =
                    ConstrainedVariable::INVALID_CONSTRAINED_VARIABLE;
                if (!(temp.constrainedDirection == ConstrainedType::SMALLER ||
                      temp.constrainedDirection == ConstrainedType::GREATER))
                {
                    std::cout << "No valid constrainedDirection for "
                              << FiniteElement::convertProcessTypeToString(
                                     temp.constrainedProcessType)
                              << " (" << tempst << ")" << std::endl;
                    _isConstrainedBC = false;
                }

                in >> tempst;  // Seepage face option (set BC to constrained
                               // value, if calculated value > constrained
                // value)
                if (tempst == "SEEPAGE")
                {
                    if (temp.constrainedDirection == ConstrainedType::SMALLER)
                        temp._isSeepageBC = true;
                    else
                        std::cout
                            << "Seepage not used as constrained direction is "
                               "not set to SMALLER.\n Please check "
                               ".bc file."
                            << std::endl;
                }
            }
            if (_isConstrainedBC)
                this->_constrainedBC.push_back(temp);
            in.clear();
        }
        //....................................................................
    }
    return position;
}

///**************************************************************************
// FEMLib-Method: CBoundaryCondition::Write
// 02/2004 OK Implementation
// 07/2007 OK LINEAR
// 10/2008 OK NOD
// 06/2009 OK MSH_TYPE off
// **************************************************************************/
// void CBoundaryCondition::Write(std::fstream* rfd_file) const
//{
//   //KEYWORD
//   *rfd_file << "#BOUNDARY_CONDITION" << "\n";
//   //--------------------------------------------------------------------
//   //NAME+NUMBER
//   *rfd_file << " $PCS_TYPE" << "\n";
//   *rfd_file << "  " << convertProcessTypeToString(getProcessType()) << "\n";
//   *rfd_file << " $PRIMARY_VARIABLE" << "\n";
//   *rfd_file << "  " <<
//   convertPrimaryVariableToString(this->getProcessPrimaryVariable()) << "\n";
//   //--------------------------------------------------------------------
//   //GEO_TYPE
//   *rfd_file << " $GEO_TYPE" << "\n";
//   *rfd_file << "  ";
//   *rfd_file << getGeoTypeAsString() << " " << geo_name << "\n";
//
//   //--------------------------------------------------------------------
//   /*OK4910
//    //MSH_TYPE
//    if(msh_node_number>0){
//    *rfd_file << " $MSH_TYPE" << "\n";
//    *rfd_file << "  ";
//    *rfd_file << "NODE" << " " << msh_node_number << "\n";
//    }
//    */
//   //--------------------------------------------------------------------
//   //DIS_TYPE
//   *rfd_file << " $DIS_TYPE" << "\n";
//   *rfd_file << "  ";
//   *rfd_file << convertDisTypeToString(this->getProcessDistributionType());
//   //switch (dis_type_name[0]) {
//   //case 'C': // Constant
//   if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
//   {
//      *rfd_file << " " << geo_node_value;
//      *rfd_file << "\n";
//      //break;
//   }
//   //case 'L': // Linear
//   else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
//   {
//      *rfd_file << " " << _PointsHaveDistribedBC.size() << "\n";
//      for (size_t i = 0; i < _PointsHaveDistribedBC.size(); i++)
//      {
//         *rfd_file << "  " << _PointsHaveDistribedBC[i] << " ";
//         *rfd_file << "  " << _DistribedBC[i] << "\n";
//      }
//      //break;
//   }
//
//   //FCT
//   if (fct_name.length() > 0)                     //OK4108
//   {
//      *rfd_file << " $FCT_TYPE" << "\n";
//      *rfd_file << "  ";
//      *rfd_file << fct_name << "\n";
//   }
//}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Write
   Task: write function
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
void CBoundaryCondition::WriteTecplot(std::fstream* tec_file) const
{
    long i;
    CGLPolyline* m_polyline1 = NULL;
    CGLPolyline* m_polyline2 = NULL;
    // list<CGLPolyline*>::const_iterator p;
    std::vector<CGLPolyline*>::iterator p;
    Surface* m_surface = NULL;
    long no_points = 0;
    std::vector<CTriangle*> triangle_vector;

    *tec_file << "VARIABLES = X,Y,Z,V1"
              << "\n";

    if (getGeoType() == GEOLIB::SURFACE)
    {
        m_surface = GEOGetSFCByName(geo_name);  // CC
        if (m_surface)
            switch (m_surface->type)
            {
                case 2:
                    p = m_surface->polyline_of_surface_vector.begin();
                    while (p != m_surface->polyline_of_surface_vector.end())
                    {
                        m_polyline1 = *p;
                        ++p;
                        m_polyline2 = *p;
                        break;
                    }
                    no_points = (long)m_polyline1->point_vector.size();
                    /*
                       for(i=0;i<no_points-1;i++) {
                       m_triangle = new CTriangle;
                       m_triangle->x[0] = m_polyline1->point_vector[i]->x;
                       m_triangle->y[0] = m_polyline1->point_vector[i]->y;
                       m_triangle->z[0] = m_polyline1->point_vector[i]->z;
                       m_triangle->x[1] = m_polyline1->point_vector[i+1]->x;
                       m_triangle->y[1] = m_polyline1->point_vector[i+1]->y;
                       m_triangle->z[1] = m_polyline1->point_vector[i+1]->z;
                       m_triangle->x[2] = m_polyline2->point_vector[i+1]->x;
                       m_triangle->y[2] = m_polyline2->point_vector[i+1]->y;
                       m_triangle->z[2] = m_polyline2->point_vector[i+1]->z;
                       triangle_vector.push_back(m_triangle);
                       m_triangle = new CTriangle;
                       m_triangle->x[0] = m_polyline2->point_vector[i]->x;
                       m_triangle->y[0] = m_polyline2->point_vector[i]->y;
                       m_triangle->z[0] = m_polyline2->point_vector[i]->z;
                       m_triangle->x[1] = m_polyline2->point_vector[i+1]->x;
                       m_triangle->y[1] = m_polyline2->point_vector[i+1]->y;
                       m_triangle->z[1] = m_polyline2->point_vector[i+1]->z;
                       m_triangle->x[2] = m_polyline1->point_vector[i+1]->x;
                       m_triangle->y[2] = m_polyline1->point_vector[i+1]->y;
                       m_triangle->z[2] = m_polyline1->point_vector[i+1]->z;
                       triangle_vector.push_back(m_triangle);
                       }
                     */
                    break;
            }
    }

    long no_nodes = 2 * no_points;
    // long no_elements = triangle_vector.size();
    long no_elements = 2 * (no_points - 1);
    // Write
    *tec_file << "ZONE T = " << geo_name << ", "
              << "N = " << no_nodes << ", "
              << "E = " << no_elements << ", "
              << "F = FEPOINT"
              << ", "
              << "ET = TRIANGLE"
              << "\n";
    if (m_polyline1)
        for (i = 0; i < no_points; i++)
            *tec_file << m_polyline1->point_vector[i]->x << " "
                      << m_polyline1->point_vector[i]->y << " "
                      << m_polyline1->point_vector[i]->z << " "
                      << geo_node_value << "\n";

    if (m_polyline2)
        for (i = 0; i < no_points; i++)
            *tec_file << m_polyline2->point_vector[i]->x << " "
                      << m_polyline2->point_vector[i]->y << " "
                      << m_polyline2->point_vector[i]->z << " "
                      << geo_node_value << "\n";

    for (i = 0; i < no_points - 1; i++)
        *tec_file << i + 1 << " " << i + 1 + 1 << " " << no_points + i + 1
                  << "\n";
    for (i = 0; i < no_points - 1; i++)
        *tec_file << no_points + i + 1 << " " << no_points + i + 1 + 1 << " "
                  << i + 1 + 1 << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   05/2010 TF changes due to new GEOLIB integration, some improvements
**************************************************************************/
bool BCRead(std::string const& file_base_name,
            const GEOLIB::GEOObjects& geo_obj,
            const std::string& unique_name)
{
    char line[MAX_ZEILE];
    std::string line_string, bc_file_name;

    // File handling
    bc_file_name = file_base_name + BC_FILE_EXTENSION;

    std::ifstream bc_file(bc_file_name.data(), std::ios::in);
    if (!bc_file.good())
    {
        ScreenMessage("! Error in BCRead: No boundary conditions !\n");
        return false;
    }

    // Keyword loop
    ScreenMessage("BCRead ... ");
    while (!bc_file.eof())
    {
        bc_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != std::string::npos)
        {
            ScreenMessage("done, read %d boundary conditions.\n",
                          bc_list.size());
            return true;
        }
        if (line_string.find("#BOUNDARY_CONDITION") != std::string::npos)
        {
            CBoundaryCondition* bc(new CBoundaryCondition());
            bool valid(true);
            std::ios::pos_type position =
                bc->Read(&bc_file, geo_obj, unique_name, valid);
            if (valid)
                bc_list.push_back(bc);
            else
                delete bc;
            bc_file.seekg(position, std::ios::beg);
        }  // keyword found
    }      // eof
    return true;
}

/**************************************************************************
   FEMLib-Method: BCWrite
   Task: master write function
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/

void BCWrite(std::string const& base_file_name)
{
    std::string sub_line;
    std::string line_string;

    // File handling
    std::string bc_file_name(base_file_name + BC_FILE_EXTENSION);
    std::fstream bc_file(bc_file_name.data(), std::ios::trunc | std::ios::out);
    bc_file.setf(std::ios::scientific, std::ios::floatfield);
    bc_file.precision(12);
    // OK string tec_file_name = base_file_name + ".tec";
    // OK fstream tec_file (tec_file_name.data(),ios::trunc|ios::out);
    // OK tec_file.setf(ios::scientific,ios::floatfield);
    // OK tec_file.precision(12);
    if (!bc_file.good())
        return;
    bc_file.seekg(0L, std::ios::beg);  // rewind?
    bc_file << "GeoSys-BC: Boundary Conditions "
               "------------------------------------------------\n";
    //========================================================================
    // BC list
    std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
    while (p_bc != bc_list.end())
    {
        FileIO::BoundaryConditionIO::write(bc_file, *(*p_bc));
        ++p_bc;
    }
    bc_file << "#STOP";
    bc_file.close();
    // OK tec_file.close();
}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   07/2007 OK V2, global function
**************************************************************************/
// CBoundaryCondition* BCGet(const std::string &pcs_name, const std::string
// &geo_type_name,
//		const std::string &geo_name)
//{
//	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
//	while (p_bc != bc_list.end()) {
//		if (((*p_bc)->pcs_type_name.compare(pcs_name) == 0)
//				&& ((*p_bc)->geo_type_name.compare(geo_type_name) == 0)
//				&& ((*p_bc)->getGeoName().compare(geo_name) == 0))
//			return *p_bc;
//		++p_bc;
//	}
//	return NULL;
//}

/**************************************************************************
   GeoSys source term function:
   02/2009 WW Implementation
**************************************************************************/
inline void CBoundaryCondition::DirectAssign(CRFProcess* m_pcs,
                                             long ShiftInNodeVector)
{
    std::string line_string;
    std::stringstream in;
    long n_index;
    double n_val;
    CBoundaryConditionNode* m_node_value = NULL;

    //========================================================================
    // File handling
    std::ifstream d_file(fname.c_str(), std::ios::in);
    // if (!st_file.good()) return;

    if (!d_file.good())
    {
        std::cout
            << "! Error in direct node source terms: Could not find file:!\n"
            << fname << "\n";
        abort();
    }
    // Rewind the file
    d_file.clear();
    d_file.seekg(0L, std::ios::beg);
    //========================================================================
    while (!d_file.eof())
    {
        line_string = readNonBlankLineFromInputStream(d_file);
        if (line_string.find("#STOP") != std::string::npos)
            break;

        in.str(line_string);
        in >> n_index >> n_val;
        in.clear();
        //
        m_node_value = new CBoundaryConditionNode;
        m_node_value->conditional = false;
        m_node_value->msh_node_number = n_index + ShiftInNodeVector;
        m_node_value->geo_node_number = n_index;
        m_node_value->node_value = n_val;
        m_node_value->CurveIndex = _curve_index;
        m_pcs->bc_node.push_back(this);
        m_pcs->bc_node_value.push_back(m_node_value);
    }  // eof
}

/**************************************************************************
   GeoSys BC function:
   03/2009 WW Implementation
**************************************************************************/
inline void CBoundaryCondition::PatchAssign(long ShiftInNodeVector)
{
    std::string line_string;
    std::stringstream in;
    long n_index;
    std::vector<long> sfc_nodes;
    CBoundaryConditionNode* m_node_value = NULL;

    CRFProcess* pcs(PCSGet(convertProcessTypeToString(this->getProcessType())));
    Surface* surface(GEOGetSFCByName(geo_name));

    // File handling
    std::ifstream d_file(fname.c_str(), std::ios::in);

    if (!d_file.good())
    {
        std::cout
            << "! Error in direct node source terms: Could not find file:!\n"
            << fname << "\n";
        abort();
    }
    // Rewind the file
    d_file.clear();
    d_file.seekg(0L, std::ios::beg);

    while (!d_file.eof())
    {
        line_string = readNonBlankLineFromInputStream(d_file);
        if (line_string.find("#STOP") != std::string::npos)
            break;

        in.str(line_string);
        in >> n_index;
        in.clear();
        sfc_nodes.push_back(n_index);
    }

    if (surface)
        pcs->m_msh->GetNODOnSFC_PLY_XY(surface, sfc_nodes, true);

    for (size_t i = 0; i < sfc_nodes.size(); i++)
    {
        m_node_value = new CBoundaryConditionNode;
        m_node_value->conditional = false;
        n_index = sfc_nodes[i];
        m_node_value->msh_node_number = n_index + ShiftInNodeVector;
        m_node_value->geo_node_number = n_index;
        m_node_value->node_value = geo_node_value;
        m_node_value->CurveIndex = _curve_index;
        pcs->bc_node.push_back(this);
        pcs->bc_node_value.push_back(m_node_value);
    }  // eof
}

CBoundaryConditionsGroup::CBoundaryConditionsGroup(void)
{
    msh_node_number_subst = -1;  //
    time_dep_bc = -1;
}

CBoundaryConditionsGroup::~CBoundaryConditionsGroup(void)
{
    /*
       int group_vector_length = group_vector.size();
       int i;
       for(i=0;i<group_vector_length;i++)
       group_vector.pop_back();
     */
    //  group_vector.clear();
}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Set
   Task: set boundary conditions
   Programing:
   02/2004 OK Implementation
   09/2004 WW Interpolation of piecewise linear BC
   02/2005 OK MSH types
   03/2005 OK MultiMSH, PNT
   08/2005 WW Changes due to the new geometry finite element.
   12/2005 OK FCT
   04/2006 WW New storage
   09/2006 WW Move linear interpolation to new MSH structure
   12/2007 WW Linear distributed BC in a surface
   10/2008 WW/CB SetTransientBCtoNodes
   last modification:
**************************************************************************/
void CBoundaryConditionsGroup::Set(CRFProcess* pcs,
                                   int ShiftInNodeVector,
                                   const double value_offset,
                                   const std::string& this_pv_name)
{
    //	long number_of_nodes = 0;
    long i, j;  // WX
    long* nodes = NULL;
    std::vector<long> nodes_vector;
    std::vector<double> node_value;
    CGLPolyline* m_polyline = NULL;
    //	CBoundaryCondition *m_bc = NULL;
    CBoundaryConditionNode* m_node_value = NULL;
    group_name = _pcs_type_name;
    bool quadratic = false;
    bool cont = false;  // WW

    if (!this_pv_name.empty())
        _pcs_pv_name = this_pv_name;
    CFEMesh* m_msh = pcs->m_msh;
    // Tests //OK

    const std::size_t previous_size = pcs->bc_node_value.size();

    if (!m_msh)
        std::cout << "Warning in CBoundaryConditionsGroup::Set - no MSH data"
                  << "\n";
    // return;
    if (m_msh)  // WW
    {
        /// In case of P_U coupling monolithic scheme
        if (pcs->type == 41)  // WW Mono
        {
            // Deform
            if (_pcs_pv_name.find("DISPLACEMENT") != std::string::npos ||
                _pcs_pv_name.find("VELOCITY_DM") != std::string::npos)
                quadratic = true;
            else
                quadratic = false;
        }
        else if (pcs->type == 4)
            quadratic = true;
        else
            quadratic = false;
        pcs->m_msh->SwitchOnQuadraticNodes(quadratic);
    }

    FiniteElement::PrimaryVariable primary_variable(
        FiniteElement::convertPrimaryVariable(_pcs_pv_name));
    std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();

    clock_t start_time(clock());

    while (p_bc != bc_list.end())
    {
        CBoundaryCondition* bc(*p_bc);
        if (bc->time_dep_interpol)  // WW/CB
        {
            ++p_bc;
            continue;
        }

        if ((bc->getProcess() == pcs) &&
            (bc->getProcessPrimaryVariable() == primary_variable))
        {
            cont = false;

            //------------------------------------------------------------------
            if (bc->getExcav() > 0 || bc->getGeoType() == GEOLIB::GEODOMAIN)
            // WX: 01.2011. boundary conditions for excavation. 03.2011.
            // Material domain BC
            {
                // GEOGetNodesInMaterialDomain(m_msh,
                // m_bc->getExcavMatGr(),nodes_vector, quadratic);
                size_t ii;
                long Size;
                int nn = 0, Domain_MG = -1;
                bool exist;
                if (bc->getGeoType() == GEOLIB::GEODOMAIN)
                    Domain_MG = atoi(bc->geo_name.c_str());
                MeshLib::CElem* elem = NULL;
                nodes_vector.resize(0);
                for (ii = 0; ii < m_msh->ele_vector.size(); ii++)
                {
                    elem = m_msh->ele_vector[ii];
                    nn = elem->GetNodesNumber(quadratic);
                    if (elem->GetPatchIndex() ==
                            static_cast<size_t>(bc->getExcavMatGr()) ||
                        static_cast<int>(elem->GetPatchIndex()) == Domain_MG)
                    {
                        Size = (int)nodes_vector.size();
                        for (i = 0; i < nn; i++)
                        {
                            exist = false;
                            for (j = 0; j < Size; j++)
                                if (elem->GetNodeIndex(i) == nodes_vector[j])
                                {
                                    exist = true;
                                    break;
                                }
                            if (!exist)
                                nodes_vector.push_back(elem->GetNodeIndex(i));
                        }
                    }
                }

                for (i = 0; i < (long)nodes_vector.size();
                     i++)  // possible nodes
                // for(int j=0;
                // j<m_msh->nod_vector[nodes_vector[i]]->connected_elements.size();
                // j++){
                // if(m_msh->ele_vector[m_msh->nod_vector[nodes_vector[i]]->connected_elements[j]]->GetPatchIndex()==m_bc->getExcavMatGr())
                {
                    m_node_value = new CBoundaryConditionNode();
                    m_node_value->msh_node_number = -1;
                    m_node_value->msh_node_number =
                        nodes_vector[i] + ShiftInNodeVector;  // nodes[i];
                    m_node_value->geo_node_number =
                        nodes_vector[i];  // nodes[i];
                    m_node_value->node_value = bc->geo_node_value;
                    m_node_value->pcs_pv_name = _pcs_pv_name;  // YD/WW
                    m_node_value->CurveIndex = bc->getCurveIndex();
                    pcs->bc_node.push_back(bc);                  // WW
                    pcs->bc_node_value.push_back(m_node_value);  // WW
                    // j=m_msh->nod_vector[nodes_vector[i]]->connected_elements.size();
                    // //only on material group boundary is selected
                }
                ++p_bc;
                continue;
            }
            //------------------------------------------------------------------
            //-- 23.02.3009. WW
            if (bc->getProcessDistributionType() == FiniteElement::DIRECT)
            {
                bc->DirectAssign(pcs, ShiftInNodeVector);
                ++p_bc;
                continue;
            }
            //................................................................
            if (bc->getGeoType() == GEOLIB::POINT)
            {
                // 05.2012. WW
                long node_idx = m_msh->GetNODOnPNT(
                    static_cast<const GEOLIB::Point*>(bc->getGeoObj()));
                if (node_idx < 0)
                {
                    ++p_bc;
                    continue;
                }
                //-------
                m_node_value = new CBoundaryConditionNode;
                // Get MSH node number
                if (bc->getProcessDistributionType() ==
                        FiniteElement::CONSTANT ||
                    bc->getProcessDistributionType() == FiniteElement::SWITCH)
                    m_node_value->geo_node_number = node_idx;

                m_node_value->conditional = cont;
                m_node_value->CurveIndex = bc->getCurveIndex();

                // Get value from a linear function. 25.08.2011. WW
                if (bc->getProcessDistributionType() == FiniteElement::FUNCTION)
                {
                    //                  a_node =
                    //                  m_msh->nod_vector[m_node_value->geo_node_number];
                    //					m_node_value->node_value =
                    // bc->dis_linear_f->getValue(a_node->X(),a_node->Y(),a_node->Z());
                    double const* const coords(
                        m_msh->nod_vector[m_node_value->geo_node_number]
                            ->getData());
                    m_node_value->node_value = bc->dis_linear_f->getValue(
                        coords[0], coords[1], coords[2]);
                }
                else
                    m_node_value->node_value = bc->geo_node_value;

                m_node_value->msh_node_number =
                    m_node_value->geo_node_number + ShiftInNodeVector;  // WW
                // YD/WW
                m_node_value->pcs_pv_name = _pcs_pv_name;
                m_node_value->msh_node_number_subst = msh_node_number_subst;
                pcs->bc_node.push_back(bc);  // WW
                // WW
                pcs->bc_node_value.push_back(m_node_value);
            }
            //--------------------------------------------------------------------
            if (bc->getGeoType() == GEOLIB::POLYLINE)
            {
                // CC
                m_polyline = GEOGetPLYByName(bc->geo_name);
                // 08/2010 TF get the polyline data structure
                GEOLIB::Polyline const* ply(
                    static_cast<const GEOLIB::Polyline*>(bc->getGeoObj()));

                if (m_polyline)
                {
                    if (bc->getProcessDistributionType() ==
                            FiniteElement::CONSTANT ||
                        bc->getProcessDistributionType() ==
                            FiniteElement::SWITCH)
                    {
                        // 08/2010 TF
                        double msh_min_edge_length = m_msh->getMinEdgeLength();
                        m_msh->setMinEdgeLength(m_polyline->epsilon);
                        std::vector<size_t> my_nodes_vector;
                        m_msh->GetNODOnPLY(ply, my_nodes_vector);
                        m_msh->setMinEdgeLength(msh_min_edge_length);
                        nodes_vector.clear();
                        for (size_t k(0); k < my_nodes_vector.size(); k++)
                            nodes_vector.push_back(my_nodes_vector[k]);

                        // for some benchmarks we need the vector entries sorted
                        // by index
                        std::sort(nodes_vector.begin(), nodes_vector.end());

                        const size_t nodes_vector_size(nodes_vector.size());
                        for (size_t i(0); i < nodes_vector_size; i++)
                        {
                            m_node_value = new CBoundaryConditionNode();
                            m_node_value->msh_node_number =
                                nodes_vector[i] + ShiftInNodeVector;
                            m_node_value->geo_node_number = nodes_vector[i];
                            // dis_prop[0];
                            if (bc->getProcessDistributionType() ==
                                FiniteElement::SWITCH)
                                m_node_value->node_value =
                                    bc->getSwitchBC().switchOnValue;
                            else
                                m_node_value->node_value = bc->geo_node_value;
                            m_node_value->CurveIndex = bc->getCurveIndex();
                            // YD/WW
                            m_node_value->pcs_pv_name = _pcs_pv_name;
                            // SB copy values 09.2012
                            m_node_value->bc_node_copy_geom = bc->copy_geom;
                            m_node_value->bc_node_copy_geom_name =
                                bc->copy_geom_name;

                            // WW
                            pcs->bc_node.push_back(bc);
                            // WW
                            pcs->bc_node_value.push_back(m_node_value);
                            // WW group_vector.push_back(m_node_value);
                            // WW bc_group_vector.push_back(m_bc); //OK
                        }
                    }

                    // Get value from a linear function. 25.06.2011. WW
                    if (bc->getProcessDistributionType() ==
                        FiniteElement::FUNCTION)
                    {
                        // 08/2010 TF
                        double msh_min_edge_length = m_msh->getMinEdgeLength();
                        m_msh->setMinEdgeLength(m_polyline->epsilon);
                        std::vector<size_t> my_nodes_vector;
                        m_msh->GetNODOnPLY(ply, my_nodes_vector);
                        m_msh->setMinEdgeLength(msh_min_edge_length);

                        nodes_vector.clear();
                        for (size_t k(0); k < my_nodes_vector.size(); k++)
                            nodes_vector.push_back(my_nodes_vector[k]);

                        // for some benchmarks we need the vector entries sorted
                        // by index
                        std::sort(nodes_vector.begin(), nodes_vector.end());

                        const size_t nodes_vector_size(nodes_vector.size());
                        for (size_t i(0); i < nodes_vector_size; i++)
                        {
                            m_node_value = new CBoundaryConditionNode();
                            m_node_value->msh_node_number =
                                nodes_vector[i] + ShiftInNodeVector;
                            m_node_value->geo_node_number = nodes_vector[i];
                            //                            a_node =
                            //                            m_msh->nod_vector[m_node_value->geo_node_number];
                            //                            m_node_value->node_value
                            //                            =
                            //                            bc->dis_linear_f->getValue(a_node->X(),a_node->Y(),a_node->Z());
                            double const* const coords(
                                m_msh->nod_vector[m_node_value->geo_node_number]
                                    ->getData());
                            m_node_value->node_value =
                                bc->dis_linear_f->getValue(
                                    coords[0], coords[1], coords[2]);

                            m_node_value->CurveIndex = bc->getCurveIndex();
                            m_node_value->pcs_pv_name = _pcs_pv_name;
                            pcs->bc_node.push_back(bc);
                            pcs->bc_node_value.push_back(m_node_value);
                        }
                    }

                    if (bc->getProcessDistributionType() ==
                        FiniteElement::GRADIENT)  // 6/2012 JOD
                    {
                        m_msh->GetNODOnPLY(ply, nodes_vector);

                        for (size_t k(0); k < nodes_vector.size(); k++)
                        {
                            m_node_value = new CBoundaryConditionNode();
                            m_node_value->msh_node_number = -1;
                            m_node_value->msh_node_number =
                                nodes_vector[k] + ShiftInNodeVector;
                            m_node_value->geo_node_number = nodes_vector[k];
                            m_node_value->node_value =
                                bc->gradient_ref_depth_gradient *
                                    (bc->gradient_ref_depth -
                                     m_msh->nod_vector[nodes_vector[k]]
                                         ->getData()[2]) +
                                bc->gradient_ref_depth_value;
                            m_node_value->CurveIndex = bc->getCurveIndex();
                            m_node_value->pcs_pv_name = _pcs_pv_name;
                            m_node_value->msh_node_number_subst =
                                msh_node_number_subst;
                            pcs->bc_node.push_back(bc);
                            pcs->bc_node_value.push_back(m_node_value);
                        }
                    }

                    // WW / TF
                    if (bc->getProcessDistributionType() ==
                        FiniteElement::LINEAR)
                    {
                        double msh_min_edge_length = m_msh->getMinEdgeLength();
                        m_msh->setMinEdgeLength(m_polyline->epsilon);
                        std::vector<size_t> my_nodes_vector;
                        m_msh->GetNODOnPLY(ply, my_nodes_vector);
                        std::vector<double> nodes_as_interpol_points;
                        m_msh->getPointsForInterpolationAlongPolyline(
                            ply, nodes_as_interpol_points);
                        m_msh->setMinEdgeLength(msh_min_edge_length);

                        nodes_vector.clear();
                        for (size_t k(0); k < my_nodes_vector.size(); k++)
                            nodes_vector.push_back(my_nodes_vector[k]);

                        std::vector<double> interpolation_points;
                        std::vector<double> interpolation_values;
                        for (size_t i(0); i < bc->getDistribedBC().size(); i++)
                        {
                            for (size_t j = 0; j < ply->getNumberOfPoints();
                                 j++)
                                if (bc->getPointsWithDistribedBC()[i] ==
                                    (int)ply->getPointID(j))
                                {
                                    if (fabs(bc->getDistribedBC()[i]) <
                                        MKleinsteZahl)
                                        bc->getDistribedBC()[i] = 1.0e-20;
                                    interpolation_points.push_back(
                                        ply->getLength(j));
                                    interpolation_values.push_back(
                                        bc->getDistribedBC()[i]);
                                    break;
                                }
                        }
                        MathLib::PiecewiseLinearInterpolation(
                            interpolation_points,
                            interpolation_values,
                            nodes_as_interpol_points,
                            node_value);

                        for (size_t i = 0; i < nodes_vector.size(); i++)
                        {
                            m_node_value = new CBoundaryConditionNode();
                            m_node_value->msh_node_number = -1;
                            m_node_value->msh_node_number =
                                nodes_vector[i] + ShiftInNodeVector;
                            m_node_value->geo_node_number = nodes_vector[i];
                            m_node_value->node_value = node_value[i];
                            // YD/WW
                            m_node_value->pcs_pv_name = _pcs_pv_name;
                            m_node_value->CurveIndex = bc->getCurveIndex();
                            // WW
                            pcs->bc_node.push_back(bc);
                            // WW
                            pcs->bc_node_value.push_back(m_node_value);
                            // WW group_vector.push_back(m_node_value);
                            // WW bc_group_vector.push_back(bc); //OK
                        }
                        node_value.clear();
                    }
                    Free(nodes);
                }  // if(m_ply)
            }
            //------------------------------------------------------------------
            if (bc->getGeoType() == GEOLIB::SURFACE)
            {
                // CC10/05
                // 04/2011 TF get the GEOLIB::Surface data structure
                GEOLIB::Surface const* sfc(
                    static_cast<const GEOLIB::Surface*>(bc->getGeoObj()));

                Surface* m_surface = GEOGetSFCByName(bc->geo_name);
                if (m_surface)
                {
                    nodes_vector.clear();

                    //					m_msh->GetNODOnSFC(m_surface,
                    // nodes_vector); #ifndef NDEBUG
                    // GEOLIB::GEOObjects const& geo_obj(*
                    // m_msh->getGEOObjects()); std::string const&
                    // geo_project_name (* m_msh->getProjectName());
                    // std::string sfc_name;
                    //					geo_obj.getSurfaceVecObj(geo_project_name)->getNameOfElement(sfc,
                    // sfc_name); 					std::string
                    // debug_fname("MeshNodesOld-BC-" + sfc_name + ".gli");
                    // std::ofstream debug_out (debug_fname.c_str());
                    // debug_out << "#POINTS" << "\n"; 					for
                    //(size_t k(0); k<nodes_vector.size(); k++) {
                    // debug_out << k
                    //<< " " <<
                    //							GEOLIB::Point((m_msh->getNodeVector())[nodes_vector[k]]->getData())
                    //<< 							" $NAME " << nodes_vector[k]
                    //<<
                    //"\n";
                    //					}
                    //					debug_out << "#STOP" << "\n";
                    //					debug_out.close();
                    //#endif
                    std::vector<size_t> msh_nod_vec;
                    double computed_search_length = m_msh->getSearchLength();
                    if (bc->epsilon != -1)
                        m_msh->setSearchLength(bc->epsilon);
                    m_msh->GetNODOnSFC(sfc, msh_nod_vec);
                    m_msh->setSearchLength(computed_search_length);

#ifndef NDEBUG
#ifdef DEBUGMESHNODESEARCH
                    {
                        std::string const debug_fname(bc->geo_name +
                                                      "-FoundNodes.gli");
                        std::ofstream debug_out(debug_fname.c_str());
                        debug_out << "#POINTS\n";
                        for (size_t k(0); k < msh_nod_vec.size(); k++)
                        {
                            debug_out
                                << k << " "
                                << GEOLIB::Point(
                                       (m_msh->getNodeVector())[msh_nod_vec[k]]
                                           ->getData())
                                << " $NAME " << msh_nod_vec[k] << "\n";
                        }
                        debug_out << "#STOP"
                                  << "\n";
                        debug_out.close();
                    }
#endif
#endif
                    //					nodes_vector.clear();
                    for (size_t k(0); k < msh_nod_vec.size(); k++)
                    {
                        //						std::cout << "\t" << k << "\t"
                        //<< nodes_vector_old[k] << "\t" <<
                        // msh_nod_vec[k]
                        //<<
                        //"\n";
                        nodes_vector.push_back(msh_nod_vec[k]);
                    }
                    size_t nodes_vector_length(nodes_vector.size());

                    if (bc->isConstrainedBC() && nodes_vector_length > 0)
                    {
                        for (std::size_t i = 0;
                             i < bc->getNumberOfConstrainedBCs();
                             i++)
                        {
                            const Constrained& temp(bc->getConstrainedBC(i));
                            if (temp.constrainedVariable ==
                                ConstrainedVariable::VELOCITY)
                            {
                                // calculate normals of triangles
                                sfc->calculateTriangleNormals();
                            }
                        }
                    }

                    if (bc->getProcessDistributionType() ==
                        FiniteElement::LINEAR)
                    {
                        std::vector<CGLPolyline*>::iterator p =
                            m_surface->polyline_of_surface_vector.begin();
                        node_value.resize(nodes_vector_length);
                        p = m_surface->polyline_of_surface_vector.begin();
                        while (p != m_surface->polyline_of_surface_vector.end())
                        {
                            m_polyline = *p;
                            for (size_t i(0); i < bc->getDistribedBC().size();
                                 i++)
                            {
                                for (size_t j = 0;
                                     j < m_polyline->point_vector.size();
                                     j++)
                                    if (bc->getPointsWithDistribedBC()[i] ==
                                        m_polyline->point_vector[j]->id)
                                    {
                                        if (fabs(bc->getDistribedBC()[i]) <
                                            MKleinsteZahl)
                                            bc->getDistribedBC()[i] = 1.0e-20;
                                        m_polyline->point_vector[j]->setPropert(
                                            bc->getDistribedBC()[i]);
                                        break;
                                    }
                            }
                            p++;
                        }
                        // WW
                        node_value.resize(nodes_vector_length);
                        bc->SurfaceInterpolation(
                            pcs, nodes_vector, node_value);  // WW
                    }

                    for (size_t i = 0; i < nodes_vector_length; i++)
                    {
                        m_node_value = new CBoundaryConditionNode();
                        m_node_value->msh_node_number = -1;
                        m_node_value->msh_node_number =
                            nodes_vector[i] + ShiftInNodeVector;  // nodes[i];

                        // nodes[i];
                        m_node_value->geo_node_number = nodes_vector[i];
                        // YD/WW
                        m_node_value->pcs_pv_name = _pcs_pv_name;

                        if (bc->getProcessDistributionType() ==
                            FiniteElement::LINEAR)
                        {
                            m_node_value->node_value = node_value[i];
                        }
                        else
                        {
                            if (bc->getProcessDistributionType() ==
                                FiniteElement::GRADIENT)
                            {  // 6/2012 JOD
                                m_node_value->node_value =
                                    bc->gradient_ref_depth_gradient *
                                        (bc->gradient_ref_depth -
                                         m_msh
                                             ->nod_vector[m_node_value
                                                              ->geo_node_number]
                                             ->getData()[2]) +
                                    bc->gradient_ref_depth_value;
                            }
                            else
                            {
                                // 25.08.2011. WW
                                if (bc->getProcessDistributionType() ==
                                    FiniteElement::FUNCTION)
                                {
                                    //                            a_node =
                                    //                            m_msh->nod_vector[m_node_value->geo_node_number];
                                    //                            m_node_value->node_value
                                    //                            =
                                    //                            bc->dis_linear_f->getValue(a_node->X(),a_node->Y(),a_node->Z());
                                    double const* const coords(
                                        m_msh
                                            ->nod_vector[m_node_value
                                                             ->geo_node_number]
                                            ->getData());
                                    m_node_value->node_value =
                                        bc->dis_linear_f->getValue(
                                            coords[0], coords[1], coords[2]);
                                }
                                else
                                {
                                    if (bc->getProcessDistributionType() ==
                                        FiniteElement::SWITCH)
                                    {
                                        m_node_value->node_value =
                                            bc->getSwitchBC().switchOnValue;
                                    }
                                    else
                                        m_node_value->node_value =
                                            bc->geo_node_value;
                                }
                            }
                        }
                        m_node_value->CurveIndex = bc->getCurveIndex();
                        // OK
                        bc->node_number_vector = nodes_vector;
                        pcs->bc_node.push_back(bc);  // WW
                        // WW
                        pcs->bc_node_value.push_back(m_node_value);
                        // WW group_vector.push_back(m_node_value);
                        // WW bc_group_vector.push_back(bc); //OK

                        if (bc->isConstrainedBC() == true &&
                            nodes_vector_length > 0)
                        {
                            for (std::size_t k = 0;
                                 k < bc->getNumberOfConstrainedBCs();
                                 k++)
                            {
                                const Constrained& temp(bc->getConstrainedBC(
                                    k));  // delete object, else previous
                                          // elements will reside here.
                                if (temp.constrainedVariable ==
                                    ConstrainedVariable::VELOCITY)
                                {
                                    double const* const coords(
                                        m_msh
                                            ->nod_vector[m_node_value
                                                             ->geo_node_number]
                                            ->getData());
                                    // works only for planar surfaces since
                                    // the normal is constant for all triangles
                                    int triangle_id(
                                        sfc->getTriangleIDOfPoint(coords));
                                    if (triangle_id != -1)
                                        m_node_value->SetNormalVector(
                                            sfc->getTriangleNormal(
                                                triangle_id));
                                    else
                                        std::cout
                                            << "Could not find current BC node "
                                            << m_node_value->geo_node_number
                                            << " on given SURFACE "
                                            << m_surface->name << std::endl;
                                }
                            }
                        }
                    }
                    node_value.clear();
                }
            }
            //------------------------------------------------------------------
            // Material domain
            //			if (bc->geo_type_name.find("MATERIAL_DOMAIN") == 0) {
            //				GEOGetNodesInMaterialDomain(m_msh, bc->_geo_type,
            //						nodes_vector, quadratic);
            //				for (i = 0; i < (long) nodes_vector.size(); i++) {
            //					m_node_value = new CBoundaryConditionNode();
            //					m_node_value->msh_node_number = -1;
            //					m_node_value->msh_node_number = nodes_vector[i]
            //							+ ShiftInNodeVector; //nodes[i];
            //					m_node_value->geo_node_number = nodes_vector[i];
            ////nodes[i]; 					m_node_value->node_value =
            /// bc->geo_node_value;
            //					m_node_value->pcs_pv_name = pcs_pv_name; //YD/WW
            //					m_node_value->CurveIndex = bc->getCurveIndex();
            //					pcs->bc_node.push_back(bc); //WW
            //					pcs->bc_node_value.push_back(m_node_value); //WW
            //					//WW group_vector.push_back(m_node_value);
            //					//WW bc_group_vector.push_back(bc); //OK
            //				}
            //			}
            //------------------------------------------------------------------
            // MSH types //OK4105
            if (bc->getMeshTypeName().compare("NODE") == 0)
            {
                m_node_value = new CBoundaryConditionNode;
                m_node_value->msh_node_number = bc->getMeshNodeNumber();
                m_node_value->geo_node_number = bc->getMeshNodeNumber();
                m_node_value->node_value = bc->geo_node_value;
                m_node_value->CurveIndex = bc->getCurveIndex();
                // YD/WW
                m_node_value->pcs_pv_name = _pcs_pv_name;
                pcs->bc_node.push_back(bc);  // WW
                // WW
                pcs->bc_node_value.push_back(m_node_value);
                // WW group_vector.push_back(m_node_value);
                // WW bc_group_vector.push_back(bc); //OK
            }
            //------------------------------------------------------------------
            // FCT types //OK
            if (bc->fct_name.size() > 0)
                // WW
                for (size_t i = 0; i < pcs->bc_node_value.size(); i++)
                {
                    pcs->bc_node_value[i]->fct_name = bc->fct_name;
                    pcs->bc_node_value[i]->msh_node_number_subst =
                        msh_node_number_subst;
                }
            // WW fct_name = bc->fct_name;
            //------------------------------------------------------------------
        }  // PCS
        ++p_bc;
    }  // list

    clock_t end_time(clock());

    ScreenMessage("\t[BC] set BC took %0.3e\n",
                  (end_time - start_time) / (double)(CLOCKS_PER_SEC));

    start_time = clock();
    // SetTransientBCtoNodes  10/2008 WW/CB Implementation
    p_bc = bc_list.begin();
    while (p_bc != bc_list.end())
    {
        CBoundaryCondition* bc(*p_bc);
        if (!bc->time_dep_interpol)  // WW/CB
        {
            ++p_bc;
            continue;
        }
        if (bc->getProcess() == pcs)
            //................................................................
            if (bc->getGeoType() == GEOLIB::POLYLINE)
            {
                // CC
                m_polyline = GEOGetPLYByName(bc->geo_name);

                if (m_polyline)
                {
                    // WW
                    if (bc->getProcessDistributionType() ==
                        FiniteElement::LINEAR)
                    {
                        // TF
                        double msh_min_edge_length = m_msh->getMinEdgeLength();
                        m_msh->setMinEdgeLength(m_polyline->epsilon);
                        std::vector<size_t> my_nodes_vector;
                        GEOLIB::Polyline const* ply(
                            static_cast<GEOLIB::Polyline const*>(
                                bc->getGeoObj()));
                        m_msh->GetNODOnPLY(ply, my_nodes_vector);
                        m_msh->setMinEdgeLength(msh_min_edge_length);

                        nodes_vector.clear();
                        for (size_t k(0); k < my_nodes_vector.size(); k++)
                            nodes_vector.push_back(my_nodes_vector[k]);

                        pcs->bc_transient_index.push_back(
                            (long)pcs->bc_node.size());
                        for (size_t i = 0; i < nodes_vector.size(); i++)
                        {
                            m_node_value = new CBoundaryConditionNode();
                            m_node_value->msh_node_number = -1;
                            m_node_value->msh_node_number =
                                nodes_vector[i] + ShiftInNodeVector;
                            m_node_value->geo_node_number = nodes_vector[i];
                            m_node_value->node_value = 0.0;
                            // YD/WW
                            m_node_value->pcs_pv_name = _pcs_pv_name;
                            m_node_value->CurveIndex = bc->getCurveIndex();
                            // WW
                            pcs->bc_node.push_back(bc);
                            // WW
                            pcs->bc_node_value.push_back(m_node_value);
                        }
                        node_value.clear();
                    }
                    //................................................................
                    // delete(values);
                    Free(nodes);
                }  // if(m_ply)
            }
        //------------------------------------------------------------------
        // PCS
        ++p_bc;
    }  // list
    /* // Make the following as comment by WW
       // Test
       long no_bc = (long)pcs->bc_node_value.size();
       if(no_bc<1)
       cout << "Warning: no boundary conditions specified for " << pcs_type_name
       << "\n";
     */
    if (std::fabs(value_offset) > 0.0)
    {
        for (std::size_t i = previous_size; i < pcs->bc_node_value.size(); i++)
        {
            pcs->bc_node_value[i]->node_value_offset = value_offset;
        }
    }

    end_time = clock();
    ScreenMessage("\t[BC] set transient BC took %0.3e\n",
                  (end_time - start_time) / (double)(CLOCKS_PER_SEC));
}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Get
   Task: set boundary conditions
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CBoundaryConditionsGroup* CBoundaryConditionsGroup::Get(
    const std::string& pcs_name)
{
    CBoundaryConditionsGroup* m_bc_group = NULL;
    std::list<CBoundaryConditionsGroup*>::const_iterator p_bc_group =
        bc_group_list.begin();
    while (p_bc_group != bc_group_list.end())
    {
        m_bc_group = *p_bc_group;
        if (m_bc_group->group_name.compare(pcs_name) == 0)
            return m_bc_group;
        ++p_bc_group;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
   last modification:
**************************************************************************/
CBoundaryConditionsGroup* BCGetGroup(const std::string& pcs_type_name,
                                     const std::string& pcs_pv_name)
{
    std::list<CBoundaryConditionsGroup*>::const_iterator it =
        bc_group_list.begin();
    while (it != bc_group_list.end())
    {
        if (((*it)->getProcessTypeName().compare(pcs_type_name) == 0) &&
            ((*it)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
            return *it;
        ++it;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void BCDelete()
{
    CBoundaryCondition* m_bc = NULL;
    std::list<CBoundaryCondition*>::const_iterator p = bc_list.begin();
    while (p != bc_list.end())
    {
        // bc_list.remove(*p);
        m_bc = *p;
        delete m_bc;
        ++p;
    }
    bc_list.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void BCGroupDelete()
{
    CBoundaryConditionsGroup* m_bc_group = NULL;
    std::list<CBoundaryConditionsGroup*>::const_iterator p =
        bc_group_list.begin();
    while (p != bc_group_list.end())
    {
        m_bc_group = *p;
        delete m_bc_group;
        // bc_group_list.remove(*p);
        ++p;
    }
    bc_group_list.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void BCGroupDelete(const std::string& pcs_type_name,
                   const std::string& pcs_pv_name)
{
    std::list<CBoundaryConditionsGroup*>::iterator p = bc_group_list.begin();
    while (p != bc_group_list.end())
    {
        if (((*p)->getProcessTypeName().compare(pcs_type_name) == 0) &&
            ((*p)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
        {
            delete *p;
            bc_group_list.erase(p);
            return;
        }
        ++p;
    }
}

/**************************************************************************
   FEMLib-Method:
   07/2007 OK Implementation
**************************************************************************/
CBoundaryCondition* BCGet(const std::string& pcs_type_name)
{
    CBoundaryCondition* m_bc = NULL;

    std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
    while (p_bc != bc_list.end())
    {
        m_bc = *p_bc;
        if (m_bc->getProcessType() ==
            FiniteElement::convertProcessType(pcs_type_name))
            return m_bc;
        ++p_bc;
    }
    return NULL;
}

/**************************************************************************
   ROCKFLOW - Funktion:
   Programming:
   11/2007 WW Implementation
**************************************************************************/
void CBoundaryCondition::SurfaceInterpolation(
    CRFProcess* m_pcs,
    std::vector<long>& nodes_on_sfc,
    std::vector<double>& node_value_vector)
{
    long i, j, k, l;

    //----------------------------------------------------------------------
    // Interpolation of polygon values to nodes_on_sfc
    int nPointsPly = 0;
    double Area1, Area2;
    double Tol = 1e-9;
    // NW. Default tolerance is 1e-9 but it can be changed in a BC file.
    if (this->epsilon != -1)
        Tol = this->epsilon;
    bool Passed;
    double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];
    //
    CGLPolyline* m_polyline = NULL;
    Surface* m_surface = NULL;
    m_surface = GEOGetSFCByName(geo_name);  // CC

    // list<CGLPolyline*>::const_iterator p =
    // m_surface->polyline_of_surface_list.begin();
    std::vector<CGLPolyline*>::iterator p =
        m_surface->polyline_of_surface_vector.begin();

    for (j = 0; j < (long)nodes_on_sfc.size(); j++)
    {
        //      pn[0] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->X();
        //      pn[1] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Y();
        //      pn[2] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Z();
        double const* const pn(
            m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->getData());
        node_value_vector[j] = 0.0;
        Passed = false;
        // nodes close to first polyline
        p = m_surface->polyline_of_surface_vector.begin();
        while (p != m_surface->polyline_of_surface_vector.end())
        {
            m_polyline = *p;
            // Gravity center of this polygon
            for (i = 0; i < 3; i++)
                gC[i] = 0.0;
            vn[2] = 0.0;
            nPointsPly = (int)m_polyline->point_vector.size();
            if (m_polyline->point_vector.front() ==
                m_polyline->point_vector.back())
                nPointsPly -= 1;
            for (i = 0; i < nPointsPly; i++)
            {
                gC[0] += m_polyline->point_vector[i]->x;
                gC[1] += m_polyline->point_vector[i]->y;
                gC[2] += m_polyline->point_vector[i]->z;
                vn[2] += m_polyline->point_vector[i]->getPropert();
            }
            for (i = 0; i < 3; i++)
                gC[i] /= (double)nPointsPly;
            // BC value at center is an average of all point values of polygon
            vn[2] /= (double)nPointsPly;
            // Area of this polygon by the gravity center
            for (i = 0; i < nPointsPly; i++)
            {
                p1[0] = m_polyline->point_vector[i]->x;
                p1[1] = m_polyline->point_vector[i]->y;
                p1[2] = m_polyline->point_vector[i]->z;
                k = i + 1;
                if (i == nPointsPly - 1)
                    k = 0;
                p2[0] = m_polyline->point_vector[k]->x;
                p2[1] = m_polyline->point_vector[k]->y;
                p2[2] = m_polyline->point_vector[k]->z;
                vn[0] = m_polyline->point_vector[i]->getPropert();
                vn[1] = m_polyline->point_vector[k]->getPropert();

                Area1 = fabs(ComputeDetTri(p1, gC, p2));

                Area2 = 0.0;
                // Check if pn is in the triangle by points (p1, gC, p2)
                Area2 = fabs(ComputeDetTri(p2, gC, pn));
                unit[0] = fabs(ComputeDetTri(gC, p1, pn));
                unit[1] = fabs(ComputeDetTri(p1, p2, pn));
                Area2 += unit[0] + unit[1];
                if (fabs(Area1 - Area2) < Tol)
                {
                    // Interpolation within a triangle (p1,p2,gC)
                    // Shape function
                    for (l = 0; l < 2; l++)
                        unit[l] /= Area1;
                    ShapeFunctionTri(NTri, unit);
                    for (l = 0; l < 3; l++)
                        node_value_vector[j] += vn[l] * NTri[l];
                    Passed = true;
                    break;
                }
            }
            //
            p++;
            if (Passed)
                break;
        }  // while
    }      // j
}
