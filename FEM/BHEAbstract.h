/**
* \file BHEAbstract.h
* 2014/06/04 HS inital implementation
* borehole heat exchanger abstract class
*
* 1) Diersch_2011_CG
* Two very important references to understand this class implementations are: 
* H.-J.G. Diersch, D. Bauer, W. Heidemann, W. Rühaak, P. Schätzl, 
* Finite element modeling of borehole heat exchanger systems: 
* Part 1. Fundamentals, Computers & Geosciences, 
* Volume 37, Issue 8, August 2011, Pages 1122-1135, ISSN 0098-3004, 
* http://dx.doi.org/10.1016/j.cageo.2010.08.003.
*
* 2) FEFLOW_2014_Springer
* FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport in Porous and Fractured Media
* Diersch, Hans-Joerg, 2014, XXXV, 996 p, Springer. 
* 
*/

#ifndef BHE_ABSTRACT_H
#define BHE_ABSTRACT_H

#include <iostream>
#include "Eigen/Eigen"
#include "../GEO/Polyline.h"
#include "FEMEnums.h"
//#include "tools.h" // HS: needed for the function GetCurveValue() 
#include <math.h>
#include "BHE_Net_ELE_Abstract.h"

namespace BHE  // namespace of borehole heat exchanger
{
	/**
	  * list the types of borehole heat exchanger
	  */
	enum BHE_TYPE {
		BHE_TYPE_2U,   // two u-tube borehole heat exchanger
		BHE_TYPE_1U,   // one u-tube borehole heat exchanger
		BHE_TYPE_CXC,  // coaxial pipe with annualar inlet
		BHE_TYPE_CXA	  // coaxial pipe with centreed inlet
	};

    enum BHE_BOUNDARY_TYPE {
        BHE_BOUND_FIXED_INFLOW_TEMP, 
        BHE_BOUND_FIXED_INFLOW_TEMP_CURVE,
        BHE_BOUND_POWER_IN_WATT,
        BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT,
        BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_DT,
        BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE,
        BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE,
        BHE_BOUND_FIXED_TEMP_DIFF,
    };

	/**
	  * discharge type of the 2U BHE
	  */
	enum BHE_DISCHARGE_TYPE {
		BHE_DISCHARGE_TYPE_PARALLEL,   // parallel discharge
		BHE_DISCHARGE_TYPE_SERIAL	   // serial discharge
	};

    class BHEAbstract : public BHE_Net_ELE_Abstract
	{
	public:
		/**
		  * constructor
		  */
		BHEAbstract(BHE_TYPE my_type, const std::string name, BHE_BOUNDARY_TYPE my_bound_type = BHE_BOUND_FIXED_INFLOW_TEMP, bool if_use_ext_Ra_Rb = false, bool user_defined_R_vals = false, int bhe_heating_cop_curve_idx = -1, int bhe_cooling_cop_curve_idx = -1, bool if_flowrate_curve = false, int n_T_in = 1, int n_T_out = 1)
			: BHE_Net_ELE_Abstract(name, BHE::BHE_NET_ELE::BHE_NET_BOREHOLE, n_T_in, n_T_out), type(my_type), _name(name), bound_type(my_bound_type), use_ext_therm_resis(if_use_ext_Ra_Rb), user_defined_therm_resis(user_defined_R_vals), _heating_cop_curve_idx(bhe_heating_cop_curve_idx), _cooling_cop_curve_idx(bhe_cooling_cop_curve_idx), use_flowrate_curve(if_flowrate_curve)
		{};

		/**
		  * destructor
		  */
		virtual ~BHEAbstract() {}; 

		/**
		  * return the number of unknowns needed for this BHE
		  * abstract function, need to be realized. 
		  */
		virtual std::size_t get_n_unknowns() = 0;

        /**
          * return the number of boundary heat exchange terms for this BHE
          * abstract function, need to be realized.
          */
        virtual std::size_t get_n_heat_exchange_terms() = 0;

        /**
          *
          */
        virtual void set_T_in_out_global_idx(std::size_t start_idx) = 0; 

        /**
          *
          */
        virtual void set_T_in_out_bottom_global_idx(std::size_t dof_BHE) = 0;

		/**
		  * return the type of the BHE
		  */
		BHE_TYPE get_type() { return type; };

        /**
          * return the type of boundary condition on this BHE
          */
        BHE_BOUNDARY_TYPE get_bound_type() { return bound_type; };

        /**
          * return the name of the BHE
          */
        const std::string get_name() { return _name;  };

		/**
		  * return the thermal resistance for the inlet pipline
		  * idx is the index, when 2U case, 
		  * 0 - the first u-tube
		  * 1 - the second u-tube
		  * needs to be overwritten
		  */
		virtual double get_thermal_resistance_fig(std::size_t idx = 0) = 0;

		/**
		  * return the thermal resistance for the outlet pipline
		  * idx is the index, when 2U case,
		  * 0 - the first u-tube
		  * 1 - the second u-tube
		  * needs to be overwritten
		  */
		virtual double get_thermal_resistance_fog(std::size_t idx = 0) = 0;

		/**
		  * return the thermal resistance
		  */
		virtual double get_thermal_resistance(std::size_t idx = 0) = 0;

		/**
		  * initialization calcultion,
		  * need to be overwritten.
		  */
		virtual void initialize()
		{
			calc_u();
			calc_Re();
			calc_Pr();
			calc_Nu();
			calc_thermal_resistances();
			calc_heat_transfer_coefficients();
		};

        /**
          * update all parameters based on the new flow rate
          * not necessarily needs to be overwritten.
          */
        virtual void update_flow_rate(double new_flow_rate)
        {
            Q_r = new_flow_rate; 
            calc_u();
            calc_Re();
            calc_Pr();
            calc_Nu();
            calc_thermal_resistances();
            calc_heat_transfer_coefficients();
        };

		virtual void update_flowrate_from_curve(double current_time) = 0;

		/**
		  * thermal resistance calculation, 
		  * need to be overwritten. 
		  */
		virtual void calc_thermal_resistances() = 0; 

		/**
		  * Nusselt number calculation,
		  * need to be overwritten.
		  */
		virtual void calc_Nu() = 0;

		/**
		  * Renolds number calculation,
		  * need to be overwritten.
		  */
		virtual void calc_Re() = 0;

		/**
		  * Prandtl number calculation,
		  * need to be overwritten.
		  */
		virtual void calc_Pr() = 0;

		/**
		  * flow velocity inside the pipeline
		  * need to be overwritten.
		  */
		virtual void calc_u() = 0;

		/**
		  * heat transfer coefficient,
		  * need to be overwritten.
		  */
		virtual void calc_heat_transfer_coefficients() = 0;

        /**
          * return the coeff of mass matrix, 
          * depending on the index of unknown. 
          */
        virtual double get_mass_coeff(std::size_t idx_unknown) = 0; 

        /**
          * return the coeff of laplace matrix,
          * depending on the index of unknown.
          */
        virtual void get_laplace_matrix(std::size_t idx_unknown, Eigen::MatrixXd & mat_laplace) = 0;

        /**
          * return the coeff of advection matrix,
          * depending on the index of unknown.
          */
        virtual void get_advection_vector(std::size_t idx_unknown, Eigen::VectorXd & vec_advection) = 0;

        /**
          * return the coeff of boundary heat exchange matrix,
          * depending on the index of unknown.
          */
        virtual double get_boundary_heat_exchange_coeff(std::size_t idx_unknown) = 0;

        /**
          * return the shift index based on primary variable value
          */
        virtual int get_loc_shift_by_pv(FiniteElement::PrimaryVariable pv_name) = 0;

        /**
          * return the number of grout zones in this BHE.  
          */
        virtual std::size_t get_n_grout_zones(void) = 0; 

        /**
          * return the inflow temperature based on outflow temperature and fixed power.
          */
        virtual double get_Tin_by_Tout(double T_in, double current_time) = 0;

        /**
          * get the polyline geometry 
          * that is representing this BHE. 
          */
        const GEOLIB::Polyline* get_geo_ply() { return _geo_ply; }

        /**
          * set the polyline geometry
          * that is representing this BHE.
          */
        void set_geo_ply(const GEOLIB::Polyline* ply) { _geo_ply = ply; }

        /**
        * set the polyline geometry
        * that is representing this BHE.
        */
        void set_ply_eps(double eps) { _ply_eps = eps; }

        /**
        * set the polyline geometry
        * that is representing this BHE.
        */
        double get_ply_eps(void) { return _ply_eps; }

		/**
		  * total refrigerant flow discharge of BHE
		  * unit is m^3/sec 
		  */
		double Q_r; 

		/**
		  * radius of the pipline inner side
		  * unit is m
		  */
		double r_inner;

		/**
		  * radius of the pipline outer side
		  * unit is m
		  */
		double r_outer;

		/**
		  * pipe-in wall thickness
		  * unit is m
		  */
		double b_in;

		/**
		  * pipe-out wall thickness
		  * unit is m
		  */
		double b_out;

		/**
		  * dynamics viscosity of the refrigerant
		  * unit is kg m-1 sec-1
		  */
		double mu_r;

		/**
		  * density of the refrigerant
		  * unit is kg m-3
		  */
		double rho_r;

		/**
		  * longitudinal dispersivity of the
		  * referigerant flow in the pipeline
		  */
		double alpha_L;

        /**
          * density of the grout
          * unit is kg m-3
          */
        double rho_g;

		/**
		  * porosity of the grout
		  * unit is [-]
		  */
		double porosity_g;

		/**
		  * specific heat capacity of the refrigerant
		  * unit is m^2 sec^-2 K^-1
		  */
		double heat_cap_r;

        /**
          * specific heat capacity of the grout
          * unit is m^2 sec^-2 K^-1
          */
        double heat_cap_g;

		/**
		  * thermal conductivity of the refrigerant
		  * unit is kg m sec^-3 K^-1
		  */
		double lambda_r;

		/**
		  * thermal conductivity of the pipe wall
		  * unit is kg m sec^-3 K^-1
		  */
		double lambda_p; 

		/**
		  * thermal conductivity of the grout
		  * unit is kg m sec^-3 K^-1
		  */
		double lambda_g;

		/**
		  * length/depth of the BHE
		  * unit is m
		  */
		double L; 

		/**
		  * diameter of the BHE
		  * unit is m
		  */
		double D;

        /**
          * power extracted from or injected into the BHE
          * unit is Watt
          * if value positive, then injecting power
          * if value negative, then extracting power
          */
        double power_in_watt_val; 

        /**
          * index of the power in watt curve
          */
        std::size_t power_in_watt_curve_idx; 

        /**
          * temperature difference between inflow and 
          * outflow pipelines
          */
        double delta_T_val;

		/**
		  * threshold Q value for switching off the BHE
		  * when using the Q_curve_fixed_dT B.C.
		  */
		double threshold;

        /**
          * whether or not using external given borehole thermal resistance values Ra, Rb
          */
        bool use_ext_therm_resis; 

        /**
          * external given borehole internal thermal resistance value
          */
        double ext_Ra; 

        /**
          * external given borehole thermal resistance value
          */
        double ext_Rb; 

		/**
		* whether or not using user defined borehole thermal resistance Rfig, Rfog, Rgg, Rgs
		*/
		bool user_defined_therm_resis;

		/**
		* external given borehole internal thermal resistance value
		*/
		double ext_Rfig;

		/**
		* external given borehole internal thermal resistance value
		*/
		double ext_Rfog;

		/**
		* external given borehole internal thermal resistance value
		*/
		double ext_Rgg1;

		/**
		* external given borehole internal thermal resistance value
		*/
		double ext_Rgg2;

		/**
		* external given borehole internal thermal resistance value
		*/
		double ext_Rgs;

        /**
        * heating COP curve index
        */
        const int _heating_cop_curve_idx;

		/**
		* cooling COP curve index
		*/
		const int _cooling_cop_curve_idx;

		/**
		* use refrigerant flow rate curve
		*/
		bool use_flowrate_curve;

		/**
		* refrigerant flow rate curve
		*/
		int flowrate_curve_idx;

        /**
          * for BHEs, the RHS value is zero 
          */
        double get_RHS_value()
        {
            return 0; 
        }

	private:

		/**
		  * the type of the BHE
		  */
		const BHE_TYPE type;

        /**
          * the type of the boundary condition on this BHE
          */
        const BHE_BOUNDARY_TYPE bound_type; 

        /**
          * the polyline geometry representing the BHE
          */
        const GEOLIB::Polyline* _geo_ply;

        /**
          * name of the borehole heat exchanger
          */
        const std::string _name;

        /**
          * epsilon value of the BHE polyline
          */
        double _ply_eps; 
	};

}  // end of namespace
#endif
