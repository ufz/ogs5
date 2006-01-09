/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: TIM
   Task: class implementation
   Programing:
   06/2004 OK Implementation
   12/2008 WW PI time step control
   last modified:
**************************************************************************/
#ifndef rf_tim_new_INC
#define rf_tim_new_INC
// C++ STL
//#include <fstream>
#include <iostream>
//#include <string>
#include <vector>
#include "makros.h" // JT2012
#include "FEMEnums.h"
// kg44 this is to synchronize time step size  for PETSC
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#include "petscksp.h"
#include "petscmat.h"
typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;
#endif
using namespace std;
//----------------------------------------------------------------
class CRFProcess;                                 //21.08.2008. WW
class CTimeDiscretization
{
private:
	double safty_coe;
	double dt_sum;                        // 17.09.2007 WW
	// For PI time step control. Aug-Nov.2008. by WW
	//Begin of data section for PI Time control ------------------------
public:                                           //OK
	double this_stepsize;
	double relative_error;
	double absolute_error;
	double reject_factor; //for automatic timestep control: BG
private:
	double h_min;
	double h_max;
	double hacc;
	double erracc;
	double dynamic_control_tolerance[DOF_NUMBER_MAX+1];	//JT2012
	std::string dynamic_control_error_method;			//JT2012
	int dynamic_time_buffer;							//JT2012
	double dynamic_minimum_suggestion;				    //JT2012
	double dynamic_failure_threshold;					//JT2012
public:                                           //OK
	int PI_tsize_ctrl_type;
private:
	std::vector<double> critical_time;
	//End of data section for PI Time control ------------------------
	friend bool IsSynCron();              //WW
public:
	std::string file_base_name;
	// TIM
	std::vector<double>time_step_vector;
	std::vector<int> time_adapt_tim_vector;
	std::vector<double>time_adapt_coe_vector;
	//WW vector<double>fixed_point_vector;

	//WW vector<double> time_step_target_vector; // kg44 for adaptive steps..intermediate time target that need to be reached
	double time_start;
	double time_end;
	double time_current;
	double time_control_manipulate;       //CMCD
	double next_active_time;
	double last_active_time;
	double recommended_time_step;
	double dt_failure_reduction_factor;
	int step_current;
	bool repeat;                          //OK/YD
	bool time_active;					//JT2012
	bool time_independence;				//JT2012
	bool last_dt_accepted;				//JT2012
	bool minimum_dt_reached;			//JT2012
	long accepted_step_count;			//JT2012
	long rejected_step_count;			//JT2012
	//
	// PCS
	std::string pcs_type_name;            //OK
	// NUM
	std::string time_type_name;           //OK
	TimeControlType::type time_control_type;
	std::string time_unit;                //WW
	double iter_times;                    //YD
	double multiply_coef;                 //YD
	double max_time_step;                 //YD
	double min_time_step;
	double initial_step_size;
	IterationType::type adapt_itr_type;
	size_t last_rejected_timestep;
	size_t stay_steps_after_rejection;
	double desired_error;
	double max_increase;
	double min_increase;
	double last_time_step_length;
	double dampening;
	double SEA_a, SEA_b, SEA_c;

	//
	//WW double minish; // JOD
	//WW int sub_steps; // JOD 4.7.10
	bool Write_tim_discrete;              //YD
	std::fstream* tim_discrete;           //YD
	double nonlinear_iteration_error;     //OK/YD
	//WW double max_adaptive_factor; // kg44
	//WW double max_adaptive_concentration_change; // kg44
public:
	CTimeDiscretization(void);
	//21.08.2008. WW
	CTimeDiscretization(const CTimeDiscretization& a_tim, std::string pcsname);
	~CTimeDiscretization(void);
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::fstream*);
	double time_step_length_neumann;      //YD
	double time_step_length;              //YD
	double CalcTimeStep(double crt_time = 0.0); // Add argument double crt_time. 25.08.2008. WW
	double FirstTimeStepEstimate();
	double AdaptiveFirstTimeStepEstimate();
	// For PI time step control. Aug-Nov.2008. by WW
	//Begin of function section for PI Time control ------------------------
	int GetPITimeStepCrtlType() const {return PI_tsize_ctrl_type; }
	double GetTimeStep() const {return this_stepsize; }
	double GetEndTime() const { return time_end; }
	void SetTimeStep( double hnew)  {this_stepsize = hnew; }
	double GetRTol() const { return relative_error; }
	double GetATol() const { return absolute_error; }
	double GetMaximumTSizeRestric() const { return h_max; }
	double GetMinimumTSizeRestric() const { return h_min; }
	double GetHacc() const { return hacc; }
	double GetErracc() const { return erracc; }
	void SetHacc(const double hacc_val)  { hacc = hacc_val; }
	void setErracc(const double erracc_val) { erracc = erracc_val; }
	//
	//JT: done differently now. // double CheckTime(double const c_time, const double dt0);
	void FillCriticalTime();              //21.08.2008.
	//Begin of function section for PI Time control ------------------------
	double ErrorControlAdaptiveTimeControl();
	double NeumannTimeControl();
	double SelfAdaptiveTimeControl();
	double StableErrorAdaptive( void );
	bool SEA_parameters_are_bad( void );
	void SEA_calc_parameters( void );
	double SEA_zbrent(const double tol);
	double SEA_func(double const c);
	inline double SEA_SIGN(const double a, const float b);
	double DynamicVariableTimeControl(); //JT2012
	double DynamicTimeSmoothing(double suggested_time_step_change);		//JT2012
	//
	// Dynamic time control methods JT2012
	bool isDynamicTimeFailureSuggested(CRFProcess *this_pcs=NULL);
	//
#ifdef GEM_REACT
	double MaxTimeStep();
#endif
	//
	//WW bool GetTimeStepTargetVector(); // kg44
	double CheckCourant();                //CMCD
};

extern std::vector<CTimeDiscretization*> time_vector;
extern bool TIMRead(std::string);
// extern CTimeDiscretization* TIMGet(const string&);
//kg44 const string made trouble for me
extern CTimeDiscretization* TIMGet(const std::string& pcs_type_name);

extern void TIMWrite(std::string);
extern bool IsSynCron();                          //WW
extern void TIMDelete();
extern void TIMDelete(std::string);
#define TIM_FILE_EXTENSION ".tim"
//ToDo
extern double aktuelle_zeit;
extern size_t aktueller_zeitschritt;
extern double dt;
extern int rwpt_numsplits;                        // JT 2010, for specifying sub time step for random walker in .tim input file
#endif
