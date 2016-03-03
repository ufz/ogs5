/*
 * FEMEnums.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FEMEnums.h"
#include <cstdlib>
#include <iostream>

namespace FiniteElement
{
ProcessType convertProcessType ( const std::string& pcs_type_string)
{
	if (pcs_type_string.compare ("LIQUID_FLOW") == 0)
		return LIQUID_FLOW;
	if (pcs_type_string.compare ("FLUID_FLOW") == 0)
		return FLUID_FLOW;
	if (pcs_type_string.compare ("TWO_PHASE_FLOW") == 0)
		return TWO_PHASE_FLOW;
	if (pcs_type_string.compare ("RICHARDS_FLOW") == 0)
		return RICHARDS_FLOW;
	if (pcs_type_string.compare ("OVERLAND_FLOW") == 0)
		return OVERLAND_FLOW;
	if (pcs_type_string.compare ("GROUNDWATER_FLOW") == 0)
		return GROUNDWATER_FLOW;
	if (pcs_type_string.compare ("HEAT_TRANSPORT") == 0)
		return HEAT_TRANSPORT;
	if (pcs_type_string.compare ("DEFORMATION") == 0)
		return DEFORMATION;
	if (pcs_type_string.compare ("DEFORMATION_FLOW") == 0)
		return DEFORMATION_FLOW;
	if (pcs_type_string.compare ("DEFORMATION_DYNAMIC") == 0)
		return DEFORMATION_DYNAMIC;
	if (pcs_type_string.compare ("MASS_TRANSPORT") == 0)
		return MASS_TRANSPORT;
	if (pcs_type_string.compare ("MULTI_PHASE_FLOW") == 0)
		return MULTI_PHASE_FLOW;
	if (pcs_type_string.compare ("DEFORMATION_H2") == 0)
		return DEFORMATION_H2;
	if (pcs_type_string.compare ("AIR_FLOW") == 0)
		return AIR_FLOW;
	if (pcs_type_string.compare ("FLUID_MOMENTUM") == 0)
		return FLUID_MOMENTUM;
	if (pcs_type_string.compare ("RANDOM_WALK") == 0)
		return RANDOM_WALK;
	if (pcs_type_string.compare ("FLUX") == 0)
		return FLUX;
	if (pcs_type_string.compare ("PS_GLOBAL") == 0)
		return PS_GLOBAL;
	if (pcs_type_string.compare ("NO_PCS") == 0)
		return NO_PCS;
	if (pcs_type_string.compare ("MULTI_COMPONENTIAL_FLOW") == 0)
		return MULTI_COMPONENTIAL_FLOW;
	if (pcs_type_string.compare ("TNEQ") == 0)
		return TNEQ;
	if (pcs_type_string.compare ("TES") == 0)
		return TES;
	//else
		//std::cout << "WARNING in convertProcessType: process type #" << pcs_type_string <<
		//"# unknown" << "\n";
	return INVALID_PROCESS;
}


std::string convertProcessTypeToString ( ProcessType pcs_type )
{
	switch (pcs_type)
	{
	case LIQUID_FLOW:
		return "LIQUID_FLOW";
	case FLUID_FLOW:
		return "FLUID_FLOW";
	case TWO_PHASE_FLOW:
		return "TWO_PHASE_FLOW";
	case RICHARDS_FLOW:
		return "RICHARDS_FLOW";
	case OVERLAND_FLOW:
		return "OVERLAND_FLOW";
	case GROUNDWATER_FLOW:
		return "GROUNDWATER_FLOW";
	case HEAT_TRANSPORT:
		return "HEAT_TRANSPORT";
	case DEFORMATION:
		return "DEFORMATION";
	case DEFORMATION_FLOW:
		return "DEFORMATION_FLOW";
	case DEFORMATION_DYNAMIC:
		return "DEFORMATION_DYNAMIC";
	case MASS_TRANSPORT:
		return "MASS_TRANSPORT";
	case MULTI_PHASE_FLOW:
		return "MULTI_PHASE_FLOW";
	case DEFORMATION_H2:
		return "DEFORMATION_H2";
	case AIR_FLOW:
		return "AIR_FLOW";
	case FLUID_MOMENTUM:
		return "FLUID_MOMENTUM";
	case RANDOM_WALK:
		return "RANDOM_WALK";
	case FLUX:
		return "FLUX";
	case PS_GLOBAL:
		return "PS_GLOBAL";
	case MULTI_COMPONENTIAL_FLOW:
		return "MULTI_COMPONENTIAL_FLOW";
	case TNEQ:
		return "TNEQ";
	case TES:
		return "TES";
	case NO_PCS:
		return "NO_PCS";
	default:
		return "INVALID_PROCESS";
	}
}


bool isFlowProcess (ProcessType pcs_type)
{
	switch (pcs_type)
	{
	case LIQUID_FLOW:
	case FLUID_FLOW:
	case RICHARDS_FLOW:
	case GROUNDWATER_FLOW:
	case PS_GLOBAL:
	case MULTI_PHASE_FLOW:
	case DEFORMATION_FLOW:
	case DEFORMATION_H2:
	case TWO_PHASE_FLOW:
	case OVERLAND_FLOW:
	case AIR_FLOW:
	case MULTI_COMPONENTIAL_FLOW:
	case TNEQ:
	case TES:
		return true;
	default:
		return false;
	}
}

bool isMultiFlowProcess (ProcessType pcs_type)
{
	if (pcs_type == PS_GLOBAL ||
		pcs_type == MULTI_PHASE_FLOW ||
		pcs_type == TWO_PHASE_FLOW ||
		pcs_type == DEFORMATION_H2)
		return true;
	return false;
}

bool isDeformationProcess (ProcessType pcs_type)
{
	if (pcs_type == DEFORMATION || pcs_type == DEFORMATION_H2 ||
	    pcs_type == DEFORMATION_FLOW || pcs_type == DEFORMATION_DYNAMIC)
		return true;
	return false;
}

const std::list<std::string> getAllProcessNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != PROCESS_END)
	{
		enum_names.push_back( convertProcessTypeToString(static_cast<ProcessType>(count++)) );
	}
	return enum_names;
}
PrimaryVariable convertPrimaryVariable ( const std::string& pcs_pv_string )
{
	if (pcs_pv_string.compare ("PRESSURE1") == 0)
		return PRESSURE;
	if (pcs_pv_string.compare ("PRESSURE2") == 0)
		return PRESSURE2;
	if (pcs_pv_string.compare ("PRESSURE_RATE1") == 0)
		return PRESSURE_RATE1;
	if (pcs_pv_string.compare ("SATURATION1") == 0)
		return SATURATION;
	if (pcs_pv_string.compare ("SATURATION2") == 0)
		return SATURATION2;
	if (pcs_pv_string.compare ("TEMPERATURE1") == 0)
		return TEMPERATURE;
	if (pcs_pv_string.compare ("TEMPERATURE2") == 0)
		return TEMPERATURE2;
	if (pcs_pv_string.compare ("DISPLACEMENT_X1") == 0)
		return DISPLACEMENT_X;
	if (pcs_pv_string.compare ("DISPLACEMENT_Y1") == 0)
		return DISPLACEMENT_Y;
	if (pcs_pv_string.compare ("DISPLACEMENT_Z1") == 0)
		return DISPLACEMENT_Z;
	if (pcs_pv_string.compare ("DISPLACEMENT_N") == 0)
		return DISPLACEMENT_N;
	if (pcs_pv_string.compare ("CONCENTRATION1") == 0)
		return CONCENTRATION;
	if (pcs_pv_string.compare ("HEAD") == 0)
		return HEAD;
	if (pcs_pv_string.compare ("VELOCITY_DM_X") == 0)
		return VELOCITY_DM_X;
	if (pcs_pv_string.compare ("VELOCITY_DM_Y") == 0)
		return VELOCITY_DM_Y;
	if (pcs_pv_string.compare ("VELOCITY_DM_Z") == 0)
		return VELOCITY_DM_Z;
	if (pcs_pv_string.compare ("VELOCITY1_X") == 0)
		return VELOCITY1_X;
	if (pcs_pv_string.compare ("VELOCITY1_Y") == 0)
		return VELOCITY1_Y;
	if (pcs_pv_string.compare ("VELOCITY1_Z") == 0)
		return VELOCITY1_Z;
	if (pcs_pv_string.compare ("STRESS_XX") == 0)
		return STRESS_XX;
	if (pcs_pv_string.compare ("STRESS_XY") == 0)
		return STRESS_XY;
	if (pcs_pv_string.compare ("STRESS_XZ") == 0)
		return STRESS_XZ;
	if (pcs_pv_string.compare ("STRESS_YY") == 0)
		return STRESS_YY;
	if (pcs_pv_string.compare ("STRESS_YZ") == 0)
		return STRESS_YZ;
	if (pcs_pv_string.compare ("STRESS_ZZ") == 0)
		return STRESS_ZZ;
	if (pcs_pv_string.compare ("ACCELERATION_X1") == 0)
		return ACCELERATION_X1;
	if (pcs_pv_string.compare ("ACCELERATION_Y1") == 0)
		return ACCELERATION_Y1;
	if (pcs_pv_string.compare ("ACCELERATION_Z1") == 0)
		return ACCELERATION_Z1;
	if (pcs_pv_string.compare ("EXCAVATION") == 0)
		return EXCAVATION;
	if (pcs_pv_string.compare ("STRAIN_XX") == 0)
		return STRAIN_XX;
	if (pcs_pv_string.compare ("STRAIN_XY") == 0)
		return STRAIN_XY;
	if (pcs_pv_string.compare ("STRAIN_XZ") == 0)
		return STRAIN_XZ;
	if (pcs_pv_string.compare ("STRAIN_YY") == 0)
		return STRAIN_YY;
	if (pcs_pv_string.compare ("STRAIN_YZ") == 0)
		return STRAIN_YZ;
	if (pcs_pv_string.compare ("STRAIN_ZZ") == 0)
		return STRAIN_ZZ;
	if (pcs_pv_string.compare ("STRAIN_PLS") == 0)
		return STRAIN_PLS;
	if (pcs_pv_string.compare ("CARBON1") == 0)
		return CARBON1;
	if (pcs_pv_string.compare ("WATER1") == 0)
		return WATER1;
		if (pcs_pv_string.compare ("METHANE1") == 0)
		return METHANE1;
			if (pcs_pv_string.compare ("NITROGEN1") == 0)
		return NITROGEN1;
	//else
	//{
		//std::cout << "convertPrimaryVariable #" << pcs_pv_string << "# not found" << "\n";
		//exit (1);
	//}
	return INVALID_PV;
}

std::string convertPrimaryVariableToString ( PrimaryVariable pcs_pv )
{
	if (pcs_pv == PRESSURE)
		return "PRESSURE1";
	if (pcs_pv == PRESSURE2)
		return "PRESSURE2";
	if (pcs_pv == TEMPERATURE)
		return "TEMPERATURE1";
	if (pcs_pv == TEMPERATURE1)
		return "TEMPERATURE1";
	if (pcs_pv == TEMPERATURE2)
		return "TEMPERATURE2";
	if (pcs_pv == CARBON1)
		return "CARBON1";
	if (pcs_pv == WATER1)
		return "WATER1";
	if (pcs_pv == METHANE1)
		return "METHANE1";
	if (pcs_pv == NITROGEN1)
		return "NITROGEN1";
	if (pcs_pv == PRESSURE_RATE1)
		return "PRESSURE_RATE1";
	if (pcs_pv == SATURATION)
		return "SATURATION1";
	if (pcs_pv == SATURATION2)
		return "SATURATION2";
	if (pcs_pv == DISPLACEMENT_X)
		return "DISPLACEMENT_X1";
	if (pcs_pv == DISPLACEMENT_Y)
		return "DISPLACEMENT_Y1";
	if (pcs_pv == DISPLACEMENT_Z)
		return "DISPLACEMENT_Z1";
	if (pcs_pv == DISPLACEMENT_N)
		return "DISPLACEMENT_N";
	if (pcs_pv == CONCENTRATION)
		return "CONCENTRATION1";
	if (pcs_pv == HEAD)
		return "HEAD";
	if (pcs_pv == VELOCITY_DM_X)
		return "VELOCITY_DM_X";
	if (pcs_pv == VELOCITY_DM_Y)
		return "VELOCITY_DM_Y";
	if (pcs_pv == VELOCITY_DM_Z)
		return "VELOCITY_DM_Z";
	if (pcs_pv == VELOCITY1_X)
		return "VELOCITY1_X";
	if (pcs_pv == VELOCITY1_Y)
		return "VELOCITY1_Y";
	if (pcs_pv == VELOCITY1_Z)
		return "VELOCITY1_Z";
	if (pcs_pv == STRESS_XX)
		return "STRESS_XX";
	if (pcs_pv == STRESS_XY)
		return "STRESS_XY";
	if (pcs_pv == STRESS_XZ)
		return "STRESS_XZ";
	if (pcs_pv == STRESS_YY)
		return "STRESS_YY";
	if (pcs_pv == STRESS_YZ)
		return "STRESS_YZ";
	if (pcs_pv == STRESS_ZZ)
		return "STRESS_ZZ";
	if (pcs_pv == STRAIN_XX) return "STRAIN_XX";
	if (pcs_pv == STRAIN_XY) return "STRAIN_XY";
	if (pcs_pv == STRAIN_XZ) return "STRAIN_XZ";
	if (pcs_pv == STRAIN_YY) return "STRAIN_YY";
	if (pcs_pv == STRAIN_YZ) return "STRAIN_YZ";
	if (pcs_pv == STRAIN_ZZ) return "STRAIN_ZZ";
	if (pcs_pv == STRAIN_PLS) return "STRAIN_PLS";
	if (pcs_pv == ACCELERATION_X1)
		return "ACCELERATION_X1";
	if (pcs_pv == ACCELERATION_Y1)
		return "ACCELERATION_Y1";
	if (pcs_pv == ACCELERATION_Z1)
		return "ACCELERATION_Z1";
	if (pcs_pv == EXCAVATION)
		return "EXCAVATION";
	return "INVALID_PRIMARY_VARIABLE";
}

const std::list<std::string> getAllPrimaryVariableNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != PV_END)
	{
		enum_names.push_back( convertPrimaryVariableToString(static_cast<PrimaryVariable>(count++)) );
	}
	return enum_names;
}

DistributionType convertDisType(const std::string& dis_type_string)
{
	if (dis_type_string.compare("CONSTANT") == 0)
		return CONSTANT;
	if (dis_type_string.compare("ANALYTICAL") == 0)
		return ANALYTICAL;
	if (dis_type_string.compare("AVERAGE") == 0)
		return AVERAGE;
	if (dis_type_string.compare("CONSTANT_GEO") == 0)
		return CONSTANT_GEO;
	if (dis_type_string.compare("GRADIENT") == 0)
		return GRADIENT;
	if (dis_type_string.compare("RESTART") == 0)
		return RESTART;
	if (dis_type_string.compare("LINEAR") == 0)
		return LINEAR;
	if (dis_type_string.compare("POINT") == 0)
		return POINT;
	if (dis_type_string.compare("CONSTANT_NEUMANN") == 0)
		return CONSTANT_NEUMANN;
	if (dis_type_string.compare("LINEAR_NEUMANN") == 0)
		return LINEAR_NEUMANN;
	if (dis_type_string.compare("NORMALDEPTH") == 0)
		return NORMALDEPTH;
	if (dis_type_string.compare("CRITICALDEPTH") == 0)
		return CRITICALDEPTH;
	if (dis_type_string.compare("GREEN_AMPT") == 0)
		return GREEN_AMPT;
	if (dis_type_string.compare("SYSTEM_DEPENDENT") == 0)
		return SYSTEM_DEPENDENT;
	if (dis_type_string.compare("PRECIPITATION") == 0)
		return PRECIPITATION;
	if (dis_type_string.compare("DIRECT") == 0)
		return DIRECT;
	if (dis_type_string.compare("RECHARGE_DIRECT") == 0)
		return RECHARGE_DIRECT;
	if (dis_type_string.compare("DOMAIN") == 0)
		return NODESCONSTANT;
	if (dis_type_string.compare("CLIMATE") == 0)
		return CLIMATE;
	if (dis_type_string.compare("RECHARGE") == 0)	//MW
		return RECHARGE;
	if (dis_type_string.compare("FUNCTION") == 0)
		return FUNCTION;                              //24.08.2011. WW
	if (dis_type_string.compare("TRANSFER_SURROUNDING") == 0)
		return TRANSFER_SURROUNDING;
	else
	{
		std::cout << "convertDisType #" << dis_type_string << "# not found"
		          << "\n";
		exit(1);
	}
	return INVALID_DIS_TYPE;
}

std::string convertDisTypeToString(DistributionType dis_type)
{
	if (dis_type == ANALYTICAL)
		return "ANALYTICAL";
	if (dis_type == AVERAGE)
		return "AVERAGE";
	if (dis_type == CONSTANT)
		return "CONSTANT";
	if (dis_type == CONSTANT_GEO)
		return "CONSTANT_GEO";
	if (dis_type == GRADIENT)
		return "GRADIENT";
	if (dis_type == RESTART)
		return "RESTART";
	if (dis_type == LINEAR)
		return "LINEAR";
	if (dis_type == POINT)
		return "POINT";
	if (dis_type == CONSTANT_NEUMANN)
		return "CONSTANT_NEUMANN";
	if (dis_type == LINEAR_NEUMANN)
		return "LINEAR_NEUMANN";
	if (dis_type == NORMALDEPTH)
		return "NORMALDEPTH";
	if (dis_type == CRITICALDEPTH)
		return "CRITICALDEPTH";
	if (dis_type == GREEN_AMPT)
		return "GREEN_AMPT";
	if (dis_type == SYSTEM_DEPENDENT)
		return "SYSTEM_DEPENDENT";
	if (dis_type == PRECIPITATION)
		return "PRECIPITATION";
	if (dis_type == DIRECT)
		return "DIRECT";
	if (dis_type == RECHARGE_DIRECT)
		return "RECHARGE_DIRECT";
	if (dis_type == NODESCONSTANT)
		return "DOMAIN";
	if (dis_type == CLIMATE)
		return "CLIMATE";
	if (dis_type == RECHARGE)	//MW
			return "RECHARGE";
	if (dis_type == FUNCTION)
		return "FUNCTION";         //24.08.2011. WW
	if (dis_type == TRANSFER_SURROUNDING)
		return "TRANSFER_SURROUNDING";

	return "INVALID_DIS_TYPE";
}

const std::list<std::string> getAllDistributionNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != DIS_END)
	{
		enum_names.push_back( convertDisTypeToString(static_cast<DistributionType>(count++)) );
	}
	return enum_names;
}

ErrorMethod convertErrorMethod(const std::string& error_method_string)
{
	if (error_method_string.compare("LMAX") == 0)
		return LMAX;
	if (error_method_string.compare("ENORM") == 0)
		return ENORM;
	if (error_method_string.compare("EVNORM") == 0)
		return EVNORM;
	if (error_method_string.compare("ERNORM") == 0)
		return ERNORM;
	if (error_method_string.compare("BNORM") == 0)
		return BNORM;
	else
	{
		std::cout << "convertErrorMethod #" << error_method_string << "# not found"<< "\n";
		exit(1);
	}
	return INVALID_ERROR_METHOD;
}

FrictionPhase convertFrictionPhase( const std::string& friction_string)
{
	if (friction_string.compare("SOLID") == 0)
		return SOLID;
	if (friction_string.compare("FLUID") == 0)
		return FLUID;
	if (friction_string.compare("NONE") == 0)
		return NONE;

	std::cout << "Convert error: " << friction_string << " not found. \n";
	return INVALID_FRICTION_TYPE;

}
std::string convertFrictionPhaseToString(FrictionPhase friction_phase)
{
	if (friction_phase == SOLID)
		return "SOLID";
	if (friction_phase == FLUID)
		return "FLUID";
	if (friction_phase == NONE)
		return "NONE";

	std::cout << "Invalid friction_phase type. \n";
	return "INVALID_FRICTION_TYPE";
}

SolidReactiveSystem convertSolidReactiveSystem( const std::string& reactive_string)
{
	if (reactive_string.compare("INERT") == 0)
		return INERT;
	if (reactive_string.compare("SINUSOIDAL") == 0)
		return SINUSOIDAL;
	if (reactive_string.compare("CaOH2") == 0)
		return CaOH2;
	if (reactive_string.compare("Mn3O4") == 0)
		return Mn3O4;
	if (reactive_string.compare("Z13XBF") == 0)
		return Z13XBF;

	std::cout << "Convert error: " << reactive_string << " not found. \n";
	return INVALID_REACTIVE_SYSTEM;
}

std::string convertSolidReactiveSystemToString(SolidReactiveSystem reactive_system)
{
	if (reactive_system == INERT)
		return "INERT";
	if (reactive_system == SINUSOIDAL)
		return "SINUSOIDAL";
	if (reactive_system == CaOH2)
		return "CaOH2";
	if (reactive_system == Mn3O4)
		return "Mn3O4";
	if (reactive_system == Z13XBF)
		return "Z13XBF";

	std::cout << "Invalid reactive system type. \n";
	return "INVALID_REACTIVE_SYSTEM";
}

} // end namespace FiniteElement

TimType::type convertTimType(const std::string& str)
{
    if (str.compare("STEADY") == 0)
        return TimType::STEADY;
    if (str.compare("TRANSIENT") == 0)
        return TimType::TRANSIENT;
    if (str.compare("PURERWPT") == 0)
        return TimType::PURERWPT;
    return TimType::INVALID_TIM_TYPE;
}

std::string convertTimTypeToString(TimType::type type)
{
    if (type == TimType::STEADY)
        return "STEADY";
    if (type == TimType::TRANSIENT)
        return "TRANSIENT";
    if (type == TimType::PURERWPT)
        return "PURERWPT";
    return "INVALID_TIM_TYPE";
}

IterationType::type convertIterationType(const std::string& str)
{
	if (str.find("LINEAR")!=std::string::npos)
		return IterationType::LINEAR;
	else if (str.find("NONLIN")!=std::string::npos)
		return IterationType::NONLINEAR;
	else if (str.find("COUPLED")!=std::string::npos)
		return IterationType::COUPLED;
	else
		return IterationType::INVALID;

}

std::string convertIterationTypeToString(IterationType::type itr_type)
{
    if (itr_type == IterationType::LINEAR)
        return "LINEAR";
    else if (itr_type == IterationType::NONLINEAR)
        return "NONLINEAR";
    else if (itr_type == IterationType::COUPLED)
        return "COUPLED";
    return "INVALID";
}

TimeControlType::type convertTimeControlType(const std::string &str)
{
    if (str == "STEPS")
        return TimeControlType::FIXED_STEPS;
    else if (str == "PI_AUTO_STEP_SIZE")
        return TimeControlType::PI_AUTO_STEP_SIZE;
    else if (str == "DYNAMIC_VARIABLE")
        return TimeControlType::DYNAMIC_VARIABLE;
    else if (str.find("DYNAMIC_COURANT")!=std::string::npos)
        return TimeControlType::DYNAMIC_COURANT;
    else if (str == "DYNAMIC_PRESSURE")
        return TimeControlType::DYNAMIC_PRESSURE;
    else if (str == "STEP_SIZE_RESTRICTION")
        return TimeControlType::STEP_SIZE_RESTRICTION;
    else if (str == "NEUMANN")
        return TimeControlType::NEUMANN;
    else if (str == "ERROR_CONTROL_ADAPTIVE")
        return TimeControlType::ERROR_CONTROL_ADAPTIVE;
    else if (str.find("SELF_ADAPTIVE")!=std::string::npos)
        return TimeControlType::SELF_ADAPTIVE;
    else if (str == "STABLE_ERROR_ADAPTIVE")
        return TimeControlType::STABLE_ERROR_ADAPTIVE;

    return TimeControlType::INVALID;
}

std::string convertTimeControlTypeToString(TimeControlType::type tc_type)
{
    if (tc_type == TimeControlType::FIXED_STEPS)
        return "STEPS";
    else if (tc_type == TimeControlType::PI_AUTO_STEP_SIZE)
        return "PI_AUTO_STEP_SIZE";
    else if (tc_type == TimeControlType::DYNAMIC_VARIABLE)
        return "DYNAMIC_VARIABLE";
    else if (tc_type == TimeControlType::DYNAMIC_COURANT)
        return "DYNAMIC_COURANT";
    else if (tc_type == TimeControlType::DYNAMIC_PRESSURE)
        return "DYNAMIC_PRESSURE";
    else if (tc_type == TimeControlType::STEP_SIZE_RESTRICTION)
        return "STEP_SIZE_RESTRICTION";
    else if (tc_type == TimeControlType::NEUMANN)
        return "NEUMANN";
    else if (tc_type == TimeControlType::ERROR_CONTROL_ADAPTIVE)
        return "ERROR_CONTROL_ADAPTIVE";
    else if (tc_type == TimeControlType::SELF_ADAPTIVE)
        return "SELF_ADAPTIVE";
    else if (tc_type == TimeControlType::STABLE_ERROR_ADAPTIVE)
        return "STABLE_ERROR_ADAPTIVE";
    return "INVALID";

}

ConstrainedType::type convertConstrainedType(const std::string &str)
{
	if (str.compare("SMALLER") == 0)
		return ConstrainedType::SMALLER;
	else if (str.compare("GREATER") == 0)
		return ConstrainedType::GREATER;
	else if (str.compare("POSITIVE") == 0)
		return ConstrainedType::POSITIVE;
	else if (str.compare("NEGATIVE") == 0)
		return ConstrainedType::NEGATIVE;
	return ConstrainedType::INVALID_CONSTRAINED_TYPE;
}

std::string convertConstrainedTypeToString(ConstrainedType::type constrained_type)
{
	if (constrained_type == ConstrainedType::SMALLER)
		return "SMALLER";
	else if (constrained_type == ConstrainedType::GREATER)
		return "GREATER";
	else if (constrained_type == ConstrainedType::POSITIVE)
		return "POSITIVE";
	else if (constrained_type == ConstrainedType::NEGATIVE)
		return "NEGATIVE";
	return "INVALID_CONSTRAINED_TYPE";
}


ConstrainedVariable::type convertConstrainedVariable(const std::string &str)
{
	if (str.compare("VELOCITY") == 0)
		return ConstrainedVariable::VELOCITY;
	return ConstrainedVariable::INVALID_CONSTRAINED_VARIABLE;
}

std::string convertConstrainedVariableToString(ConstrainedVariable::type constrained_variable)
{
	if (constrained_variable == ConstrainedVariable::VELOCITY)
		return "VELOCITY";
	return "INVALID_CONSTRAINED_VARIABLE";
}

