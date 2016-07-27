#include "BHE_Net_ELE_HeatPump.h"
#include "tools.h"

using namespace BHE;

BHE_Net_ELE_HeatPump::BHE_Net_ELE_HeatPump(std::string & name)
    : BHE_Net_ELE_Abstract(name, BHE_NET_ELE::BHE_NET_HEATPUMP, 1, 1)
{
    _heat_pump_BC_type = HEAT_PUMP_BOUND_POWER_FIXED_DT; 
}

double BHE_Net_ELE_HeatPump::set_BC(double T_in, double current_time)
{
	double T_out = 0.0;
	double power_hp = 0.0;
	double power_bhe = 0.0;
	int flag_valid = false;
	double COP = 0.0;
	double power_el = 0.0;
	double delta_T = 0.0;

	switch (_heat_pump_BC_type)
	{
	case HEAT_PUMP_BOUND_POWER_FIXED_DT:
		break;
	case HEAT_PUMP_BOUND_POWER_FIXED_FLOWRATE:
		double rho_cp_u = _fluid_density * _fluid_heat_capacity * _flowrate;

		if (_power_curve_idx <= 0)
			power_hp = _power_val;
		else
			power_hp = GetCurveValue(_power_curve_idx, 0, current_time, &flag_valid);
		COP = GetCurveValue(_cop_curve_idx, 0, T_in, &flag_valid);
		power_bhe = power_hp * (COP - 1.0) / COP;
		// also how much power from electricity
		power_el = power_hp - power_bhe;
		if (fabs(power_hp) < 0.1)
		{
			T_out = T_in;
		}
		else
		{
			delta_T = -power_bhe / rho_cp_u;
			T_out = T_in - delta_T;
		}
		std::cout << "heat pump: " << this->get_ele_name() << ", T_in: " << T_in << ", T_out: " << T_out << std::endl;
		std::cout << "COP: " << COP << ", Q_bhe: " << power_bhe << ", Q_elect: " << power_el << std::endl;
		break;
	}

	return T_out;
}

double BHE_Net_ELE_HeatPump::get_RHS_value()
{
    double rt_RHS_val = 0.0;

    // depending on the boundary condition, 
    // calculate the RHS value
    switch (_heat_pump_BC_type)
    {
    case HEAT_PUMP_BOUND_POWER_FIXED_DT:
        rt_RHS_val = _delta_T_val;
        break;
    case HEAT_PUMP_BOUND_POWER_FIXED_FLOWRATE:
        // TODO
        break;
    }
    return rt_RHS_val;
}