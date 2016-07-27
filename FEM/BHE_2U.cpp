/**
* \file BHE_2U.cpp
* 2014/06/04 HS inital implementation
* class of borehole heat exchanger with 2 U-tube
*
*/

#include "BHE_2U.h"
#include "makros.h"
#include "tools.h"

using namespace BHE; 

/**
  * return the thermal resistance for the inlet pipline
  * idx is the index, when 2U case,
  * 0 - the first u-tube
  * 1 - the second u-tube
  */
double BHE_2U::get_thermal_resistance_fig(std::size_t idx = 0)
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
double BHE_2U::get_thermal_resistance_fog(std::size_t idx = 0)
{
	// TODO
	return 0.0;
}

/**
  * return the thermal resistance
  */
double BHE_2U::get_thermal_resistance(std::size_t idx = 0)
{
	// TODO
	return 0.0; 
}

void BHE_2U::set_T_in_out_global_idx(std::size_t start_idx)
{
    // T_in_1
    this->set_T_in_global_index(start_idx, 0);
    // T_out_1
    this->set_T_out_global_index(start_idx + 1, 0); // TO CHECK
    // T_in_2
    this->set_T_in_global_index(start_idx + 2, 1); // TO CHECK
    // T_out_2
    this->set_T_out_global_index(start_idx + 3, 1); // TO CHECK
}

void BHE::BHE_2U::set_T_in_out_bottom_global_idx(std::size_t dof_bhe)
{
    std::size_t start_idx;
    std::size_t global_idx_T_in_bottom;

    // calculating
    start_idx = this->get_T_in_global_index();
    global_idx_T_in_bottom = start_idx + dof_bhe - 8; // 2U BHE, the order is: T_in_1, T_out_1, T_in_2, T_out_2, T_g1, T_g2, T_g3, T_g4. 

    // T_in_1 at the bottom
    this->set_T_in_bottom_global_index(global_idx_T_in_bottom, 0);
    // T_out_1 at the bottom
    this->set_T_out_bottom_global_index(global_idx_T_in_bottom + 1, 0);
    // T_in_2 at the bottom
    this->set_T_in_bottom_global_index(global_idx_T_in_bottom + 2, 1);
    // T_out_2 at the bottom
    this->set_T_out_bottom_global_index(global_idx_T_in_bottom + 3, 1);
}

/**
  * calculate thermal resistance
  */
void BHE_2U::calc_thermal_resistances()
{
	// thermal resistance due to advective flow of refrigerant in the pipes
	// Eq. 31 in Diersch_2011_CG
	_R_adv_i1 = 1.0 / (_Nu(0) * lambda_r * PI); 
	_R_adv_o1 = 1.0 / (_Nu(1) * lambda_r * PI);
	_R_adv_i2 = 1.0 / (_Nu(2) * lambda_r * PI);
	_R_adv_o2 = 1.0 / (_Nu(3) * lambda_r * PI);

	// thermal resistance due to thermal conductivity of the pip wall material
	// Eq. 36 in Diersch_2011_CG
	double _R_con_a;
	_R_con_a = std::log(r_outer / r_inner) / ( 2.0 * PI * lambda_p ); 

	// thermal resistance due to the grout transition
	double chi;
	double d0; // the average outer diameter of the pipes
	double s; // diagonal distances of pipes
    double R_adv, R_con; 
	d0 = 2.0 * r_inner; 
	s = omega * std::sqrt(2); 
	chi = std::log(std::sqrt(D*D + 4 * d0*d0) / 2 / std::sqrt(2) / d0) / std::log(D / 2 / d0);
    // Eq. 36
    _R_con_a_i1 = _R_con_a_i2 = _R_con_a_o1 = _R_con_a_o2 = std::log(r_outer / r_inner) / (2.0 * PI * lambda_p);

    if (use_ext_therm_resis)
    {
        R_adv = 0.25 * (_R_adv_i1 + _R_adv_i2 + _R_adv_o1 + _R_adv_o2);
        R_con = 0.25 * (_R_con_a_i1 + _R_con_a_i2 + _R_con_a_o1 + _R_con_a_o2);
        _R_g = 4 * ext_Rb - R_adv - R_con;
    }
    else
    {
        _R_g = acosh((D*D + d0*d0 - s*s) / (2 * D*d0)) / (2 * PI * lambda_g * lambda_g) * (3.098 - 4.432 * s / D + 2.364 * s * s / D / D);
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
		_R_fig = _R_adv_i1 + _R_adv_i2 + _R_con_a_i1 + _R_con_a_i2 + _R_con_b;
		_R_fog = _R_adv_o1 + _R_adv_o2 + _R_con_a_o1 + _R_con_a_o2 + _R_con_b;
	}


	// thermal resistance due to grout-soil exchange
	if (user_defined_therm_resis)
		_R_gs = ext_Rgs;
	else
		_R_gs = (1 - chi)*_R_g;

	// thermal resistance due to inter-grout exchange
	double R_ar_1, R_ar_2; 
    if (use_ext_therm_resis)
    {
        R_ar_1 = (2.0 + std::sqrt(2.0)) * _R_g * (ext_Ra - R_adv - R_con) / (_R_g + ext_Ra - R_adv - R_con);
        R_ar_2 = std::sqrt(2.0) * R_ar_1;
    }
    else
    {
        R_ar_1 = acosh((s*s - d0*d0) / d0 / d0) / (2.0 * PI * lambda_g);
        R_ar_2 = acosh((2.0*s*s - d0*d0) / d0 / d0) / (2.0 * PI * lambda_g);
    }
	if (user_defined_therm_resis)
	{
		_R_gg_1 = ext_Rgg1;
		_R_gg_2 = ext_Rgg2;
	}
	else
	{
		_R_gg_1 = 2.0 * _R_gs * (R_ar_1 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_1 + 2.0 * chi * _R_g);
		_R_gg_2 = 2.0 * _R_gs * (R_ar_2 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_2 + 2.0 * chi * _R_g);
	}
	

	if (!std::isfinite(_R_gg_1) || !std::isfinite(_R_gg_2))
    {
        std::cout << "Error!!! Grout Thermal Resistance is an infinite number! The simulation will be stopped! \n";
        exit(1);
    }

	// check if constraints regarding negative thermal resistances are violated
	// apply correction procedure
	// Section (1.5.5) in FEFLOW White Papers Vol V.
	double constraint1 = 1.0 / ((1.0 / _R_gg_1) + (1.0 / (2.0 * _R_gs)));
	double constraint2 = 1.0 / ((1.0 / _R_gg_2) + (1.0 / (2.0 * _R_gs)));
	int count = 0;
	while (constraint1 < 0.0 || constraint2 < 0.0)
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
			_R_gg_1 = 2.0 * _R_gs * (R_ar_1 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_1 + 2.0 * chi * _R_g);
			_R_gg_2 = 2.0 * _R_gs * (R_ar_2 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_2 + 2.0 * chi * _R_g);
		}
		if (count == 1)
		{
			chi *= 0.5;
			_R_gs = (1 - chi)*_R_g;
			_R_gg_1 = 2.0 * _R_gs * (R_ar_1 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_1 + 2.0 * chi * _R_g);
			_R_gg_2 = 2.0 * _R_gs * (R_ar_2 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_2 + 2.0 * chi * _R_g);
		}
		if (count == 2)
		{
			chi = 0.0;
			_R_gs = (1 - chi)*_R_g;
			_R_gg_1 = 2.0 * _R_gs * (R_ar_1 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_1 + 2.0 * chi * _R_g);
			_R_gg_2 = 2.0 * _R_gs * (R_ar_2 - 2.0 * chi * _R_g) / (2.0 * _R_gs - R_ar_2 + 2.0 * chi * _R_g);
			break;
		}
		std::cout << "Warning! Correction procedure was applied due to negative thermal resistance! Correction step #" << count << "\n";
		constraint1 = 1.0 / ((1.0 / _R_gg_1) + (1.0 / (2.0 * _R_gs)));
		constraint2 = 1.0 / ((1.0 / _R_gg_2) + (1.0 / (2.0 * _R_gs)));
		count++;
	}
}

/**
  * Nusselt number calculation
  */
void BHE_2U::calc_Nu()
{
	// see Eq. 32 in Diersch_2011_CG

	double tmp_Nu = 0.0;
	double gamma, xi; 
	double d; 

	d = 2.0 * r_inner; 

	if ( _Re < 2300.0 )
	{
		tmp_Nu = 4.364; 
	}
	else if ( _Re >= 2300.0 && _Re < 10000.0 )
	{
		gamma = (_Re - 2300) / (10000 - 2300);

		tmp_Nu = (1.0 - gamma) * 4.364; 
		tmp_Nu += gamma * ( ( 0.0308 / 8.0 * 1.0e4 * _Pr ) / ( 1.0 + 12.7 * std::sqrt(0.0308/8.0) * ( std::pow(_Pr, 2.0/3.0) - 1.0 ) ) * ( 1.0 + std::pow( d / L, 2.0/3.0 ) ) ); 

	}
	else if ( _Re > 10000.0 )
	{
		xi = pow(1.8 * std::log10(_Re) - 1.5, -2.0); 
		tmp_Nu = (xi / 8.0 * _Re * _Pr) / (1.0 + 12.7 * std::sqrt(xi / 8.0) * (std::pow(_Pr, 2.0 / 3.0) - 1.0)) * (1.0 + std::pow(d / L, 2.0 / 3.0));
	}

	_Nu(0) = tmp_Nu; 
	_Nu(1) = tmp_Nu;
	_Nu(2) = tmp_Nu;
	_Nu(3) = tmp_Nu;
}

/**
  * Renolds number calculation
  */
void BHE_2U::calc_Re()
{
	double u_norm, d; 
	u_norm = _u.norm();
	d = 2.0 * r_inner; // inner diameter of the pipeline

	_Re = u_norm * d / (mu_r / rho_r) ; 
}

/**
  * Prandtl number calculation
  */
void BHE_2U::calc_Pr()
{
	_Pr = mu_r * heat_cap_r / lambda_r; 
}

/**
  * calculate heat transfer coefficient
  */
void BHE_2U::calc_heat_transfer_coefficients()
{
	_PHI_fig = 1.0 / _R_fig ;
	_PHI_fog = 1.0 / _R_fog ;
	_PHI_gg_1 = 1.0 / _R_gg_1 ;
	_PHI_gg_2 = 1.0 / _R_gg_2 ;
	_PHI_gs = 1.0 / _R_gs ;
}

/**
  * flow velocity inside the pipeline
  */
void BHE_2U::calc_u()
{
	double tmp_u; 

	// which discharge type it is? 
	if (_discharge_type == BHE::BHE_DISCHARGE_TYPE_PARALLEL)
	{
		tmp_u = Q_r / (2.0 * PI * r_inner * r_inner);
	}
	else  // serial discharge type
	{
		tmp_u = Q_r / (PI * r_inner * r_inner);
	}

	_u(0) = tmp_u; 
	_u(1) = tmp_u;
	_u(2) = tmp_u;
	_u(3) = tmp_u;
}

double BHE_2U::get_mass_coeff(std::size_t idx_unknown)
{
    double mass_coeff = 0.0;

    switch (idx_unknown)
    {
    case 0:  // i1
        mass_coeff = rho_r * heat_cap_r * CSA_i;
        break;
    case 1:  // i2
        mass_coeff = rho_r * heat_cap_r * CSA_i;
        break;
    case 2:  // o1
        mass_coeff = rho_r * heat_cap_r * CSA_o;
        break;
    case 3:  // o2
        mass_coeff = rho_r * heat_cap_r * CSA_o;
        break;
    case 4:  // g1
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1;
        break;
    case 5:  // g2
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g1;
        break;
    case 6:  // g3
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2;
        break;
    case 7:  // g4
        mass_coeff = (1.0 - porosity_g) * rho_g * heat_cap_g * CSA_g2;
        break;
    default:
        break;
    }

    return mass_coeff;
}

void BHE_2U::get_laplace_matrix(std::size_t idx_unknown, Eigen::MatrixXd & mat_laplace)
{
	// Here we calculates the laplace coefficients in the governing 
	// equations of BHE. These governing equations can be found in 
	// 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
	// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 10-18. 
    double laplace_coeff(0.0);
    mat_laplace.setZero();

	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 18
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i;
		break;
	case 1:
		// pipe i2, Eq. 18
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_i;
		break;
	case 2:
		// pipe o1,  Eq. 18
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o;
		break;
	case 3:
		// pipe o2,  Eq. 18
        laplace_coeff = (lambda_r + rho_r * heat_cap_r * alpha_L * _u.norm()) * CSA_o;
		break;
	case 4:
		// pipe g1, Eq. 14
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g1;
		break;
	case 5:
		// pipe g2, Eq. 15
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g1;
		break;
	case 6:
		// pipe g3, Eq. 16
        laplace_coeff = (1.0 - porosity_g) * lambda_g * CSA_g2;
		break;
	case 7:
		// pipe g4, Eq. 17
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

void BHE_2U::get_advection_vector(std::size_t idx_unknown, Eigen::VectorXd & vec_advection)
{
    vec_advection.setZero();

    double advection_coeff(0);
	switch (idx_unknown)
	{
	case 0:
		// pipe i1, Eq. 10
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_i;
        // z direction 
        vec_advection(2) = -1.0 * advection_coeff;
		break;
	case 1:
		// pipe i2, Eq. 11
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_i;
        // z direction 
        vec_advection(2) = -1.0 * advection_coeff;
		break;
	case 2:
		// pipe o1, Eq. 12
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_o;
        // z direction 
        vec_advection(2) = advection_coeff;
		break;
	case 3:
		// pipe o2, Eq. 13
        advection_coeff = rho_r * heat_cap_r * _u(0) * CSA_o;
        // z direction 
        vec_advection(2) = advection_coeff;
		break;
	case 4:
		// pipe g1, Eq. 14
		advection_coeff = 0.0;
		break;
	case 5:
		// pipe g2, Eq. 15
		advection_coeff = 0.0;
		break;
	case 6:
		// pipe g3, Eq. 16
		advection_coeff = 0.0;
		break;
	case 7:
		// pipe g4, Eq. 17
		advection_coeff = 0.0;
		break;
	default:
		std::cout << "Error !!! The index passed to get_advection_coeff for BHE is not correct. \n";
		exit(1);
		break;
	}
}

double BHE_2U::get_boundary_heat_exchange_coeff(std::size_t idx_unknown)
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
        // PHI_gg_1
        exchange_coeff = _PHI_gg_1;
        break;
    case 3:
        // PHI_gg_2
        exchange_coeff = _PHI_gg_2;
        break;
    case 4:
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

int BHE_2U::get_loc_shift_by_pv(FiniteElement::PrimaryVariable pv_name)
{
    int idx(0);

    if (pv_name == FiniteElement::TEMPERATURE_IN_1)
        idx = 0;
    else if (pv_name == FiniteElement::TEMPERATURE_IN_2)
        idx = 1;
    else if (pv_name == FiniteElement::TEMPERATURE_OUT_1)
        idx = 2;
    else if (pv_name == FiniteElement::TEMPERATURE_OUT_2)
        idx = 3;
    else if (pv_name == FiniteElement::TEMPERATURE_G_1)
        idx = 4;
    else if (pv_name == FiniteElement::TEMPERATURE_G_2)
        idx = 5;
    else if (pv_name == FiniteElement::TEMPERATURE_G_3)
        idx = 6;
    else if (pv_name == FiniteElement::TEMPERATURE_G_4)
        idx = 7;

    return idx;
}

double BHE_2U::get_Tin_by_Tout(double T_out, double current_time = -1.0)
{
    double T_in(0.0);
    double power_tmp(0.0);
    int flag_valid = true;

    switch (this->get_bound_type())
    {
    case BHE_BOUND_POWER_IN_WATT:
        T_in = power_in_watt_val / Q_r / heat_cap_r / rho_r + T_out;
        break;
    case BHE_BOUND_FIXED_TEMP_DIFF:
        T_in = T_out + delta_T_val;
        break;
    case BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT:
        // TODO
        std::cout << "BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT feature has not been implemented yet. " << std::endl;
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

