/**
* \file BHE_1U.cpp
* 2014/06/04 HS inital implementation
* class of borehole heat exchanger with 1 U-tube
*
*/

#include "BHE_1U.h"
#include "makros.h"
#include "tools.h"

using namespace BHE;

/**
* return  the thermal resistance for the inlet pipline
* idx is the index, when 2U case,
* 0 - the first u-tube
* 1 - the second u-tube
*/
double BHE_1U::get_thermal_resistance_fig(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

/**
* return the thermal resistance for the outlet pipline
* idx is the index, when 2U case,
* 0 - the first u-tube
* 1 - the second u-tube
*/
double BHE_1U::get_thermal_resistance_fog(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

/**
* return the thermal resistance
*/
double BHE_1U::get_thermal_resistance(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

void BHE_1U::set_T_in_out_global_idx(std::size_t start_idx)
{
    // T_in
    this->set_T_in_global_index(start_idx);
    // T_out
    this->set_T_out_global_index(start_idx + 1);
}

void BHE_1U::set_T_in_out_bottom_global_idx(std::size_t dof_bhe)
{
    std::size_t start_idx;
    std::size_t global_idx_T_in_bottom; 

    // calculating
    start_idx = this->get_T_in_global_index(); 
    //global_idx_T_in_bottom = start_idx + dof_bhe - 3; // 1U BHE, the order is: T_in, T_out, T_g1, T_g2. 
	global_idx_T_in_bottom = start_idx + dof_bhe - 4; // 1U BHE, the order is: T_in, T_out, T_g1, T_g2. 

    // T_in at the bottom
    this->set_T_in_bottom_global_index(global_idx_T_in_bottom);
    // T_out at the bottom
    this->set_T_out_bottom_global_index(global_idx_T_in_bottom + 1);
}

/**
* calculate thermal resistance
*/
void BHE::BHE_1U::calc_thermal_resistances()
{
	// thermal resistance due to thermal conductivity of the pip wall material
	// Eq. 36 in Diersch_2011_CG
	_R_adv_i1 = 1.0 / (_Nu(0) * lambda_r * PI);
	_R_adv_o1 = 1.0 / (_Nu(1) * lambda_r * PI);


	// thermal resistance due to the grout transition
	double chi;
	double d0; // the average outer diameter of the pipes
	double s; // diagonal distances of pipes
    double R_adv, R_con;

	d0 = 2.0 * r_outer;
	s = omega * std::sqrt(2);
    // Eq. 49
    _R_con_a_i1 = _R_con_a_o1 = std::log(r_outer / r_inner) / (2.0 * PI * lambda_p);
	// Eq. 51
	chi = std::log(std::sqrt(D*D + 2 * d0*d0) / 2 / d0) / std::log(D / std::sqrt(2) / d0);
    if (use_ext_therm_resis)
    {
        R_adv = 0.5 * (_R_adv_i1 + _R_adv_o1); 
        R_con = 0.5 * (_R_con_a_i1 + _R_con_a_o1); 
        _R_g = 2 * ext_Rb - R_adv - R_con; 
    }
	else
    {
        // Eq. 52
        _R_g = acosh((D*D + d0*d0 - omega*omega) / (2 * D*d0)) / (2 * PI * lambda_g) * (1.601 - 0.888 * omega / D);
    }
	_R_con_b = chi * _R_g;
	// Eq. 29 and 30
	if (user_defined_therm_resis)
	{
		_R_fig = ext_Rfig;
		_R_fog = ext_Rfog;
	}
	else
	{
		_R_fig = _R_adv_i1 + _R_con_a_i1 + _R_con_b;
		_R_fog = _R_adv_o1 + _R_con_a_o1 + _R_con_b;
	}

	// thermal resistance due to grout-soil exchange
	if (user_defined_therm_resis)
		_R_gs = ext_Rgs;
	else
		_R_gs = (1 - chi)*_R_g;

	// thermal resistance due to inter-grout exchange
	double R_ar;
    if (use_ext_therm_resis)
    {
        R_ar = ext_Ra - 2 * (R_adv + R_con);
    }
    else
    {
        R_ar = acosh((2.0*omega*omega - d0*d0) / d0 / d0) / (2.0 * PI * lambda_g);
    }
    
	if (user_defined_therm_resis)
		_R_gg = ext_Rgg1;
	else
		_R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);

	if (!std::isfinite(_R_gg))
    {
        std::cout << "Error!!! Grout Thermal Resistance is an infinite number! The simulation will be stopped! \n" ;
        exit(1);
    }

	// check if constraints regarding negative thermal resistances are violated
	// apply correction procedure
	// Section (1.5.5) in FEFLOW White Papers Vol V.
	double constraint = 1.0 / ((1.0 / _R_gg) + (1.0 / (2.0 * _R_gs)));
	int count = 0;
	while (constraint < 0.0)
	{
		if (user_defined_therm_resis || use_ext_therm_resis)
		{
			std::cout << "Error!!! Constraints on thermal resistances are violated! Correction procedure can't be applied due to user defined thermal resistances! The simulation will be stopped! \n";
			exit(1);
		}
		if (count == 0)
		{
			chi *= 0.66;
			_R_gs = (1 - chi)*_R_g;
			_R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
		}
		if (count == 1)
		{
			chi *= 0.5;
			_R_gs = (1 - chi)*_R_g;
			_R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
		}
		if (count == 2)
		{
			chi = 0.0;
			_R_gs = (1 - chi)*_R_g;
			_R_gg = 2.0 * _R_gs * (R_ar - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar + 2.0 * chi * _R_g);
			break;
		}
		std::cout << "Warning! Correction procedure was applied due to negative thermal resistance! Correction step #" << count << "\n";
		constraint = 1.0 / ((1.0 / _R_gg) + (1.0 / (2.0 * _R_gs)));
		count++;
	}

	// print R and phi values
	//std::cout << "Rfig =" << _R_fig << " Rfog =" << _R_fog << " Rgg =" << _R_gg << " Rgs =" << _R_gs << "\n";
	double phi_fig = 1.0 / (_R_fig * S_i);
	double phi_fog = 1.0 / (_R_fog * S_o);
	double phi_gg = 1.0 / (_R_gg * S_g1);
	double phi_gs = 1.0 / (_R_gs * S_gs);
	//std::cout << "phi_fig =" << phi_fig << " phi_fog =" << phi_fog << " phi_gg =" << phi_gg << " phi_gs =" << phi_gs << "\n";
}

/**
* Nusselt number calculation
*/
void BHE_1U::calc_Nu()
{
	// see Eq. 32 in Diersch_2011_CG

	double tmp_Nu = 0.0;
	double gamma, xi;
	double d;

	d = 2.0 * r_inner;

	if (_Re < 2300.0)
	{
		tmp_Nu = 4.364;
	}
	else if (_Re >= 2300.0 && _Re < 10000.0)
	{
		gamma = (_Re - 2300) / (10000 - 2300);

		tmp_Nu = (1.0 - gamma) * 4.364;
		tmp_Nu += gamma * ((0.0308 / 8.0 * 1.0e4 * _Pr) / (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d / L, 2.0 / 3.0)));

	}
	else if (_Re > 10000.0)
	{
		xi = pow(1.8 * std::log10(_Re) - 1.5, -2.0);
		tmp_Nu = (xi / 8.0 * _Re * _Pr) / (1.0 + 12.7 * std::sqrt(xi / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d / L, 2.0 / 3.0));
	}

	_Nu(0) = tmp_Nu;
	_Nu(1) = tmp_Nu;
}

/**
* Renolds number calculation
*/
void BHE_1U::calc_Re()
{
    double u_norm, d;
    u_norm = std::abs(_u(0));
	d = 2.0 * r_inner; // inner diameter of the pipeline

	_Re = u_norm * d / (mu_r / rho_r);
}

/**
* Prandtl number calculation
*/
void BHE_1U::calc_Pr()
{
	_Pr = mu_r * heat_cap_r / lambda_r;
}

/**
* calculate heat transfer coefficient
*/
void BHE_1U::calc_heat_transfer_coefficients()
{
    _PHI_fig = 1.0 / _R_fig ;
    _PHI_fog = 1.0 / _R_fog ;
    _PHI_gg = 1.0 / _R_gg ;
    _PHI_gs = 1.0 / _R_gs ;
}

/**
* flow velocity inside the pipeline
*/
void BHE_1U::calc_u()
{
	double tmp_u;

	tmp_u = Q_r / ( PI * r_inner * r_inner);

	_u(0) = tmp_u;
	_u(1) = tmp_u;
}

double BHE_1U::get_mass_coeff(std::size_t idx_unknown)
{
    double mass_coeff = 0.0; 

    switch (idx_unknown)
    {
    case 0:  // i1
        mass_coeff = rho_r * heat_cap_r * CSA_i; 
        break;
    case 1:  // o1
        mass_coeff = rho_r * heat_cap_r * CSA_o;
        break;
    case 2:  // g1
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1;
        break;
    case 3:  // g2
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2;
        break;
    default:
        break;
    }

    return mass_coeff; 
}

void BHE_1U::get_laplace_matrix(std::size_t idx_unknown, Eigen::MatrixXd & mat_laplace)
{
	// Here we calculates the laplace coefficients in the governing 
	// equations of BHE. These governing equations can be found in 
	// 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
	// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22. 
	double laplace_coeff(0.0);
    mat_laplace.setZero(); 

	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 19
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i;
		break;
	case 1:
		// pipe o1, Eq. 20
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o;
		break;
	case 2:
		// pipe g1, Eq. 21
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g1;
		break;
	case 3:
		// pipe g2, Eq. 22
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g2;
		break;
	default:
		std::cout << "Error !!! The index passed to get_laplace_coeff for BHE is not correct. \n";
		exit(1);
		break;
	}

    mat_laplace(0, 0) = laplace_coeff; 
    mat_laplace(1, 1) = laplace_coeff;
    mat_laplace(2, 2) = laplace_coeff;
}

void BHE_1U::get_advection_vector(std::size_t idx_unknown, Eigen::VectorXd & vec_advection)
{
    double advection_coeff(0);
    vec_advection.setZero();

	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 19
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_i;
        // z direction 
        vec_advection(2) = -1.0 * advection_coeff;
		break;
	case 1:
		// pipe o1, Eq. 20
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_o;
        // z direction 
        vec_advection(2) = 1.0 * advection_coeff;
		break;
	case 2:
		// grout g1, Eq. 21
		advection_coeff = 0.0;
		break;
	case 3:
		// grout g2, Eq. 22
		advection_coeff = 0.0;
		break;
	default:
		std::cout << "Error !!! The index passed to get_advection_coeff for BHE is not correct. \n";
		exit(1);
		break;
	}
}

double BHE_1U::get_boundary_heat_exchange_coeff(std::size_t idx_unknown)
{
    // Here we calculates the boundary heat exchange coefficients 
    // in the governing equations of BHE. 
    // These governing equations can be found in 
    // 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97. 

    double exchange_coeff(0);

    switch (idx_unknown)
    {
    case 0:
        // PHI_fig
        exchange_coeff = _PHI_fig;
        break;
    case 1:
        // PHI_fog
        exchange_coeff = _PHI_fog;
        break;
    case 2:
        // PHI_gg
        exchange_coeff = _PHI_gg;
        break;
    case 3:
        // PHI_gs
        exchange_coeff = _PHI_gs;
        break;
    default:
        std::cout << "Error !!! The index passed to get_boundary_heat_exchange_coeff for BHE is not correct. \n";
        exit(1);
        break;
    }
    return exchange_coeff;
}

int BHE_1U::get_loc_shift_by_pv(FiniteElement::PrimaryVariable pv_name)
{
    int idx(0);

    if (pv_name == FiniteElement::TEMPERATURE_IN_1)
        idx = 0;
    else if (pv_name == FiniteElement::TEMPERATURE_OUT_1)
        idx = 1;
    else if (pv_name == FiniteElement::TEMPERATURE_G_1)
        idx = 2;
    else if (pv_name == FiniteElement::TEMPERATURE_G_2)
        idx = 3; 
    
    return idx;
}

double BHE_1U::get_Tin_by_Tout(double T_out, double current_time = -1.0)
{
    double T_in(0.0);
    double power_tmp(0.0); 
    double building_power_tmp(0.0);
    double power_elect_tmp(0.0);
    int flag_valid = false; 
    double Q_r_tmp(0.0);
    double COP_tmp(0.0);
	double fac_dT = 1.0;

    switch (this->get_bound_type())
    {
    case BHE_BOUND_POWER_IN_WATT:
		if (use_flowrate_curve)
		{
			Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time, &flag_valid);
			update_flow_rate(Q_r_tmp);
		}
		else
			Q_r_tmp = Q_r;
        T_in = power_in_watt_val / Q_r_tmp / heat_cap_r / rho_r + T_out;
        break;
    case BHE_BOUND_FIXED_TEMP_DIFF: 
		if (use_flowrate_curve)
		{
			Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time, &flag_valid);
			update_flow_rate(Q_r_tmp);
		}
        T_in = T_out + delta_T_val;
        break; 
    case BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT:
        // get the power value in the curve
        power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid);
		if (power_tmp < 0)
			fac_dT = -1.0;
		else
			fac_dT = 1.0;
		// if power value exceeds threshold, calculate new values
		if (fabs(power_tmp) > threshold)
		{
			// calculate the corresponding flow rate needed
			// using the defined delta_T value
			Q_r_tmp = power_tmp / (fac_dT*delta_T_val) / heat_cap_r / rho_r;
			// update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out + (fac_dT*delta_T_val);
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
		}
		else
		{
			Q_r_tmp = 1.0e-12; // this has to be a small value to avoid division by zero
			// update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out;
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
		}
        break; 
    case BHE::BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_DT: 
        // get the building power value in the curve
        building_power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid); 
		if (building_power_tmp <= 0.0)
		{
			// get COP value based on T_out in the curve
			COP_tmp = GetCurveValue(_heating_cop_curve_idx, 0, T_out, &flag_valid);
			// now calculate how much power needed from BHE
			power_tmp = building_power_tmp * (COP_tmp - 1.0) / COP_tmp;
			// also how much power from electricity
			power_elect_tmp = building_power_tmp - power_tmp;
			// print the amount of power needed
			//std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp << ", Q_elect: " << power_elect_tmp << std::endl;
			fac_dT = -1.0;
		}
		else
		{
			// get COP value based on T_out in the curve
			COP_tmp = GetCurveValue(_cooling_cop_curve_idx, 0, T_out, &flag_valid);
			// now calculate how much power needed from BHE
			power_tmp = building_power_tmp * (COP_tmp + 1.0) / COP_tmp;
			// also how much power from electricity
			power_elect_tmp = -building_power_tmp + power_tmp;
			// print the amount of power needed
			//std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp << ", Q_elect: " << power_elect_tmp << std::endl;
			fac_dT = 1.0;
		}
        // if power value exceeds threshold, calculate new values
        if (fabs(power_tmp) > threshold)
        {
            // calculate the corresponding flow rate needed
            // using the defined delta_T value
            Q_r_tmp = power_tmp / (fac_dT*delta_T_val) / heat_cap_r / rho_r;
            // update all values dependent on the flow rate
            update_flow_rate(Q_r_tmp);
            // calculate the new T_in
            T_in = T_out + (fac_dT*delta_T_val);
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
        }
        else
        {
            Q_r_tmp = 1.0e-12; // this has to be a small value to avoid division by zero
            // update all values dependent on the flow rate
            update_flow_rate(Q_r_tmp);
            // calculate the new T_in
            T_in = T_out;
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
        }
        break;
    case BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE:
        // get the building power value in the curve
        building_power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid);
		if (building_power_tmp <= 0)
		{
			// get COP value based on T_out in the curve
			COP_tmp = GetCurveValue(_heating_cop_curve_idx, 0, T_out, &flag_valid);
			// now calculate how much power needed from BHE
			power_tmp = building_power_tmp * (COP_tmp - 1.0) / COP_tmp;
			// also how much power from electricity
			power_elect_tmp = building_power_tmp - power_tmp;
			// print the amount of power needed
			//std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp << ", Q_elect: " << power_elect_tmp << std::endl;
		}
		else
		{
			// get COP value based on T_out in the curve
			COP_tmp = GetCurveValue(_cooling_cop_curve_idx, 0, T_out, &flag_valid);
			// now calculate how much power needed from BHE
			power_tmp = building_power_tmp * (COP_tmp + 1.0) / COP_tmp;
			// also how much power from electricity
			power_elect_tmp = -building_power_tmp + power_tmp;
			// print the amount of power needed
			//std::cout << "COP: " << COP_tmp << ", Q_bhe: " << power_tmp << ", Q_elect: " << power_elect_tmp << std::endl;
		}
		// Assign Qr whether from curve or fixed value
		if (use_flowrate_curve)
		{
			Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time, &flag_valid);
			update_flow_rate(Q_r_tmp);
		}
		else
			Q_r_tmp = Q_r;
		if(fabs(power_tmp) < threshold)
		{
			Q_r_tmp = 1.0e-12; // this has to be a small value to avoid division by zero
							   // update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out;
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
		}
		else
		{
			Q_r_tmp = Q_r;
			update_flow_rate(Q_r_tmp);
			// calculate the dT value based on fixed flow rate
			delta_T_val = power_tmp / Q_r_tmp / heat_cap_r / rho_r;
			// calcuate the new T_in 
			T_in = T_out + delta_T_val;
		}
        break;
    case BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE: 
        // get the power value in the curve
        power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid);
		// Assign Qr whether from curve or fixed value
		if (use_flowrate_curve)
		{
			Q_r_tmp = GetCurveValue(flowrate_curve_idx, 0, current_time, &flag_valid);
			update_flow_rate(Q_r_tmp);
		}
		else
			Q_r_tmp = Q_r;
        // calculate the dT value based on fixed flow rate
		if (fabs(power_tmp) < threshold)
		{
			Q_r_tmp = 1.0e-12; // this has to be a small value to avoid division by zero
							   // update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out;
			// print out updated flow rate
			//std::cout << "Qr: " << Q_r_tmp << std::endl;
		}
		else
		{
			Q_r_tmp = Q_r;
			update_flow_rate(Q_r_tmp);
			// calculate the dT value based on fixed flow rate
			delta_T_val = power_tmp / Q_r_tmp / heat_cap_r / rho_r;
			// calcuate the new T_in 
			T_in = T_out + delta_T_val;
		}
        break; 
    default:
        T_in = T_out; 
        break;
    }

    return T_in; 
}

