/**
* \file BHE_Net_ELE_HeatPump.h
* 2015/08/28 HS inital implementation
* borehole heat exchanger network element heat pump class
*
*
*/

#ifndef BHE_NET_ELE_HEATPUMP_H
#define BHE_NET_ELE_HEATPUMP_H

#include "BHE_Net_ELE_Abstract.h"

namespace BHE  // namespace of borehole heat exchanger
{
enum HEAT_PUMP_BOUNDARY_TYPE {
	HEAT_PUMP_BOUND_POWER_FIXED_FLOWRATE,
	HEAT_PUMP_BOUND_POWER_FIXED_DT,
};

class BHE_Net_ELE_HeatPump : public BHE_Net_ELE_Abstract 
{
public: 
    BHE_Net_ELE_HeatPump(std::string & name);

	double set_BC(double T_in, double current_time);

	double get_flowrate()
	{
		return _flowrate;
	}
                    
    double get_RHS_value();

    void set_delta_T_val(double new_T_val)
    {
        _delta_T_val = new_T_val; 
    }

    double get_delta_T_val()
    {
        return _delta_T_val;
    }

    void set_heat_pump_BC_type(HEAT_PUMP_BOUNDARY_TYPE type)
    {
        _heat_pump_BC_type = type; 
    }

	void set_power_val(double power_val)
	{
		_power_val = power_val;
	}

	void set_flowrate(double flowrate)
	{
		_flowrate = flowrate;
	}

	void set_power_curve_idx(int power_curve_idx)
	{
		_power_curve_idx = power_curve_idx;
	}

	void set_cop_curve_idx(int cop_curve_idx)
	{
		_cop_curve_idx = cop_curve_idx;
	}

	void set_fluid_density(double density)
	{
		_fluid_density = density;
	}

	void set_fluid_heat_capacity(double cp)
	{
		_fluid_heat_capacity = cp;
	}

private:

    /**
      * T_in - T_out value
      */
    double _delta_T_val; 

	double _power_val;
	double _flowrate;
	int _power_curve_idx;
	int _cop_curve_idx;
	double _fluid_density;
	double _fluid_heat_capacity;

    HEAT_PUMP_BOUNDARY_TYPE _heat_pump_BC_type; 
};

}

#endif