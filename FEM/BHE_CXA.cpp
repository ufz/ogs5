/**
* \file BHE_CXA.cpp
* 2014/06/04 HS inital implementation
* class of borehole heat exchanger with coaxial annular pipeline
*
*/

#include "BHE_CXA.h"
#include "makros.h"
#include "tools.h"

using namespace BHE;
/**
* return the thermal resistance for the inlet pipline
* idx is the index, when 2U case,
* 0 - the first u-tube
* 1 - the second u-tube
*/
double BHE_CXA::get_thermal_resistance_fig(std::size_t idx = 0)
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
double BHE_CXA::get_thermal_resistance_fog(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

/**
* return the thermal resistance
*/
double BHE_CXA::get_thermal_resistance(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

void BHE_CXA::set_T_in_out_global_idx(std::size_t start_idx)
{
    // 
    this->set_T_in_global_index(start_idx);
    this->set_T_out_global_index(start_idx + 1);
}

void BHE::BHE_CXA::set_T_in_out_bottom_global_idx(std::size_t dof_bhe)
{
    std::size_t start_idx;
    std::size_t global_idx_T_in_bottom;

    // calculating
    start_idx = this->get_T_in_global_index();
    global_idx_T_in_bottom = start_idx + dof_bhe - 3; // CXA BHE, the order is: T_in, T_out, T_g. 

    // T_in at the bottom
    this->set_T_in_bottom_global_index(global_idx_T_in_bottom);
    // T_out at the bottom
    this->set_T_out_bottom_global_index(global_idx_T_in_bottom + 1);
}

/**
* calculate thermal resistance
*/
void BHE_CXA::calc_thermal_resistances()
{
	double Nu_in, Nu_out; 
	double d_o1, d_i1, d_h;
	double chi;
	double _R_con_i1, _R_con_o1;

	Nu_in = _Nu(0);
	Nu_out = _Nu(1);
	d_o1 = 2.0 * (r_inner + b_in);
	d_i1 = 2.0 * r_outer;
	d_h = 2.0 * (r_outer - (r_inner + b_in));

	// thermal resistance due to advective flow of refrigerant in the pipes
	// Eq. 58, 59, and 60 in Diersch_2011_CG
	_R_adv_o1 = 1.0 / (Nu_out * lambda_r * PI);
	_R_adv_a_i1 = 1.0 / (Nu_in * lambda_r * PI) * ( d_h / d_o1) ;
	_R_adv_b_i1 = 1.0 / (Nu_in * lambda_r * PI) * ( d_h / d_i1);

	// thermal resistance due to thermal conductivity of the pip wall material
	// Eq. 66 in Diersch_2011_CG
	_R_con_i1 = std::log( (r_outer + b_out) / r_outer) / (2.0 * PI * lambda_p);
	_R_con_o1 = std::log( (r_inner + b_out) / r_inner) / (2.0 * PI * lambda_p);

	// thermal resistance due to the grout transition
	d_i1 = 2.0 * (r_outer + b_out);
	// Eq. 68
	chi = std::log(std::sqrt(D*D + d_i1*d_i1) / std::sqrt(2) / d_i1) / std::log(D / d_i1);
    if (use_ext_therm_resis)
    {
        _R_g = ext_Rb - _R_adv_b_i1 - _R_con_i1; 
    }
    else
    {
        // Eq. 69
        _R_g = std::log(D / d_i1) / (2.0 * PI * lambda_g);
    }
	// Eq. 67
	_R_con_b = chi * _R_g;
    if (use_ext_therm_resis)
    {
        _R_ff = ext_Ra;
    }
	else if (user_defined_therm_resis)
	{
		_R_ff = ext_Rgg1; // Attention! Here ext_Rgg1 is treated as Rff for coaxial type
	}
	else
    {
        // Eq. 56 
        _R_ff = _R_adv_o1 + _R_adv_a_i1 + _R_con_o1;
    }

    // Eq. 57
	if (user_defined_therm_resis)
		_R_fig = ext_Rfig;
	else
		_R_fig = _R_adv_b_i1 + _R_con_i1 + _R_con_b;
	// thermal resistance due to grout-soil exchange
	if (user_defined_therm_resis)
		_R_gs = ext_Rgs;
	else
		_R_gs = (1 - chi)*_R_g;

	if ( !std::isfinite(_R_gs) )
    {
        std::cout << "Error!!! Grout Thermal Resistance is an infinite number! The simulation will be stopped! \n";
        exit(1);
    }

}

/**
* Nusselt number calculation
*/
void BHE_CXA::calc_Nu()
{
	// see Eq. 32 in Diersch_2011_CG

	double Nu_in(0.0), Nu_out(0.0);
	double gamma, xi;
	double d_o1, d_i1, d_h;

	d_o1 = 2.0 * r_inner;
	d_i1 = 2.0 * r_outer;

	// first calculating Nu_out
	if (_Re_o1 < 2300.0)
	{
		Nu_out = 4.364;
	}
	else if (_Re_o1 >= 2300.0 && _Re_o1 < 10000.0)
	{
		gamma = (_Re_o1 - 2300) / (10000 - 2300);

		Nu_out = (1.0 - gamma) * 4.364;
		Nu_out += gamma * ((0.0308 / 8.0 * 1.0e4 * _Pr) / (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d_o1 / L, 2.0 / 3.0)));

	}
	else if (_Re_o1 > 10000.0)
	{
		xi = pow(1.8 * std::log10(_Re_o1) - 1.5, -2.0);
		Nu_out = (xi / 8.0 * _Re_o1 * _Pr) / (1.0 + 12.7 * std::sqrt(xi / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d_o1 / L, 2.0 / 3.0));
	}

	d_o1 = 2.0 * (r_inner + b_in);
	d_h = 2.0 * (r_outer - (r_inner + b_in));
	// then calculating Nu_in
	if (_Re_i1 < 2300.0)
	{
		Nu_in = 3.66;
		Nu_in += (4.0 - 0.102 / (d_o1 / d_i1 + 0.02)) * pow(d_o1 / d_i1, 0.04);
	}
	else if (_Re_i1 >= 2300.0 && _Re_i1 < 10000.0)
	{
		gamma = (_Re_i1 - 2300) / (10000 - 2300);

		Nu_in = (1.0 - gamma) * (3.66 + (4.0 - 0.102 / (d_i1 / d_o1 + 0.02))) * pow(d_i1 / d_o1, 0.04);
		Nu_in += gamma * ((0.0308 / 8.0 * 1.0e4 * _Pr) / (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d_h / L, 2.0 / 3.0)) * ((0.86 * std::pow(d_i1 / d_o1, 0.84) + 1.0 - 0.14*std::pow(d_i1 / d_o1, 0.6)) / (1.0 + d_i1 / d_o1)));

	}
	else if (_Re_i1 > 10000.0)
	{
		xi = pow(1.8 * std::log10(_Re_i1) - 1.5, -2.0);
		Nu_in = (xi / 8.0 * _Re_i1 * _Pr) / (1.0 + 12.7 * std::sqrt(xi / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d_h / L, 2.0 / 3.0)) * ((0.86 * std::pow(d_o1 / d_i1, 0.84) + 1.0 - 0.14*std::pow(d_o1 / d_i1, 0.6)) / (1.0 + d_o1 / d_i1));
	}

	// _Nu(0) is Nu_in, and _Nu(1) is Nu_out
	_Nu(0) = Nu_in;
	_Nu(1) = Nu_out;
}

/**
* Renolds number calculation
*/
void BHE_CXA::calc_Re()
{
	double d_o1, d_h;

	d_o1 = 2.0 * r_inner; // inner diameter of the pipeline
	d_h = 2.0 * (r_outer - (r_inner + b_in));

	// _u(0) is u_in, and _u(1) is u_out
	_Re_o1 = _u(1) * d_o1 / (mu_r / rho_r);
	_Re_i1 = _u(0) * d_h  / (mu_r / rho_r);
}

/**
* Prandtl number calculation
*/
void BHE_CXA::calc_Pr()
{
	_Pr = mu_r * heat_cap_r / lambda_r;
}

/**
* calculate heat transfer coefficient
*/
void BHE_CXA::calc_heat_transfer_coefficients()
{
	_PHI_fig = 1.0 / _R_fig ;
	_PHI_ff  = 1.0 / _R_ff ;
	_PHI_gs = 1.0 / _R_gs ;
}

/**
* flow velocity inside the pipeline
*/
void BHE_CXA::calc_u()
{
	double u_in, u_out;

	u_in = Q_r / (PI * (r_outer * r_outer - (r_inner + b_in) * (r_inner + b_in)));
	u_out = Q_r / (PI * r_inner * r_inner);

	_u(0) = u_in;
	_u(1) = u_out;
}

double BHE_CXA::get_mass_coeff(std::size_t idx_unknown)
{
    double mass_coeff = 0.0;

    switch (idx_unknown)
    {
    case 0:  // i1
        mass_coeff = rho_r * heat_cap_r * CSA_i;
        break;
    case 1:  // i2
        mass_coeff = rho_r * heat_cap_r * CSA_o;
        break;
    case 2:  // o1
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g;
        break;
    default:
        break;
    }

    return mass_coeff;
}

void BHE_CXA::get_laplace_matrix(std::size_t idx_unknown, Eigen::MatrixXd & mat_laplace)
{
	// Here we calculates the laplace coefficients in the governing 
	// equations of BHE. These governing equations can be found in 
	// 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
	// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 23-25. 
	double laplace_coeff(0.0);
    mat_laplace.setZero();

	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 23
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i;
		break;
	case 1:
		// pipe o1, Eq. 24
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o;
		break;
	case 2:
		// pipe g1, Eq. 25
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g;
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

void BHE_CXA::get_advection_vector(std::size_t idx_unknown, Eigen::VectorXd & vec_advection)
{
    double advection_coeff(0);
    vec_advection.setZero();

	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 23
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_i;
        // z direction 
        vec_advection(2) = -1.0 * advection_coeff;
		break;
	case 1:
		// pipe o1, Eq. 24
        advection_coeff = rho_r * heat_cap_r * _u(1) * CSA_o;
        // z direction 
        vec_advection(2) = advection_coeff;
		break;
	case 2:
		// pipe g1, Eq. 25
		advection_coeff = 0.0;
		break;
	default:
		std::cout << "Error !!! The index passed to get_advection_coeff for BHE is not correct. \n";
		exit(1);
		break;
	}
}

double BHE_CXA::get_boundary_heat_exchange_coeff(std::size_t idx_unknown)
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
        // PHI_ff
        exchange_coeff = _PHI_ff;
        break;
    case 2:
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

int BHE_CXA::get_loc_shift_by_pv(FiniteElement::PrimaryVariable pv_name)
{
    int idx(0);

    if (pv_name == FiniteElement::TEMPERATURE_IN_1)
        idx = 0;
    else if (pv_name == FiniteElement::TEMPERATURE_OUT_1)
        idx = 1;
    else if (pv_name == FiniteElement::TEMPERATURE_G_1)
        idx = 2;

    return idx;
}

double BHE_CXA::get_Tin_by_Tout(double T_out, double current_time = -1.0)
{
    double T_in(0.0);
    double power_tmp(0.0);
    int flag_valid = true;
	double Q_r_tmp(0.0);

    switch (this->get_bound_type())
    {
    case BHE_BOUND_POWER_IN_WATT:
        T_in = power_in_watt_val / Q_r / heat_cap_r / rho_r + T_out;
        break;
    case BHE_BOUND_FIXED_TEMP_DIFF:
        T_in = T_out + delta_T_val;
        break;
	case BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT:
		// get the power value in the curve
		power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid);
		// if power value exceeds threshold, calculate new values
		if (fabs(power_tmp) > threshold)
		{
			// calculate the corresponding flow rate needed
			// using the defined delta_T value
			Q_r_tmp = power_tmp / delta_T_val / heat_cap_r / rho_r;
			// update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out + delta_T_val;
		}
		else
		{
			Q_r_tmp = 1.0e-06; // this has to be a small value to avoid division by zero
			// update all values dependent on the flow rate
			update_flow_rate(Q_r_tmp);
			// calculate the new T_in
			T_in = T_out;
		}
		break;
    case BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE:
        // get the power value in the curve
        power_tmp = GetCurveValue(power_in_watt_curve_idx, 0, current_time, &flag_valid);
        // calculate the dT value based on fixed flow rate
        delta_T_val = power_tmp / Q_r / heat_cap_r / rho_r;
        // calcuate the new T_in 
        T_in = T_out + delta_T_val;
        break;
    default:
        T_in = T_out;
        break;
    }

    return T_in;
}

