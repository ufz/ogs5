/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: NUM
   Task: class implementation
   Programing:
   11/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_num_new_INC
#define rf_num_new_INC

#include "makros.h"  // JT
#include "FEMEnums.h"

#define NUM_FILE_EXTENSION ".num"
// C++ STL
#include "prototyp.h"
#include <fstream>
#include <string>
#include <vector>

//----------------------------------------------------------------
class CNumerics
{
public:
    int getNumIntegrationSamplePoints() const { return ele_gauss_points; }
    // method
    std::string method_name;  // OK
    // PCS
    std::string pcs_type_name;
    // RENUMBER
    int renumber_method;
    int renumber_parameter;
    //
    // LS - Linear Solver
    int ls_method;
    int ls_max_iterations;
    int ls_error_method;
    double ls_error_tolerance;
    double ls_theta;
    int ls_precond;
    int ls_storage_method;
    std::string ls_extra_arg;  // NW
    //
    // NLS - Non-linear Solver
    std::string nls_method_name;
    int nls_method;        // Picard or Newton
    int nls_error_method;  // WW
    int nls_max_iterations;
    double nls_relaxation;
    double
        nls_error_tolerance[DOF_NUMBER_MAX];  // JT2012: array function of dof
    double nls_plasticity_local_tolerance;
    void setNonLinearErrorMethod(FiniteElement::ErrorMethod nls_method)
    {
        _pcs_nls_error_method = nls_method;
    }
    FiniteElement::ErrorMethod getNonLinearErrorMethod() const
    {
        return _pcs_nls_error_method;
    }
    //
    // CPL WW
    std::string cpl_variable;      // MB
    std::string cpl_process;       // JT
    std::string cpl_variable_JOD;  // JT->JOD. This one defaults to FLUX. I'm
                                   // not sure what you want to do with it, but
                                   // cpl_variable must default to "NONE".
    int cpl_max_iterations;
    int cpl_min_iterations;  // JT2012
    double
        cpl_error_tolerance[DOF_NUMBER_MAX];  // JT2012: array function of dof
    bool cpl_error_specified;                 // JT2012
    bool cpl_master_process;
    void setCouplingErrorMethod(FiniteElement::ErrorMethod cpl_method)
    {
        _pcs_cpl_error_method = cpl_method;
    }
    FiniteElement::ErrorMethod getCouplingErrorMethod() const
    {
        return _pcs_cpl_error_method;
    }
    // local_picard1
    double local_picard1_tolerance;
    int local_picard1_max_iterations;
    // velocity update within picard
    int update_velocity_within_nonlinear;
    // ELE
    int ele_gauss_points;  // probably element-type-wise
    int ele_mass_lumping;
    int ele_upwind_method;  // CB
    double ele_upwinding;
    int ele_supg_method;              // NW
    int ele_supg_method_length;       // NW
    int ele_supg_method_diffusivity;  // NW
    // FEM-FCT
    int fct_method;                    // NW
    unsigned int fct_prelimiter_type;  // NW
    double fct_const_alpha;            // NW
    // Deformation
    int GravityProfile;
    // LAGRANGE method //OK
    double lag_quality;
    int lag_max_steps;
    double lag_local_eps;
    int lag_time_weighting;
    double lag_min_weight;
    int lag_use_matrix;
    int lag_vel_method;
    double newton_damping_factor;
    double newton_damping_tolerance;
    double nls_abs_residual_tolerance;
    double nls_abs_unknown_tolerance;
    double nls_rel_unknown_tolerance;
    //
    // Configure
    void NumConfigure(bool overall_coupling_exists);  // JT2012
    //
    // Dynamics
    bool CheckDynamic();
    double GetDynamicDamping_beta1() const { return DynamicDamping[0]; }
    double GetDynamicDamping_beta2() const { return DynamicDamping[1]; }
    double GetDynamicDamping_bbeta() const { return DynamicDamping[2]; }
    //
    /// For GMRES. WW
    long Get_m() const { return m_cols; }
    CNumerics(std::string);
    ~CNumerics(void);
    std::ios::pos_type Read(std::ifstream*);
    void Write(std::fstream*);

#ifdef USE_PETSC
    const char* getLinearSolverName() const { return lsover_name.c_str(); }
    const char* getPreconditionerName() const { return pres_name.c_str(); }
#endif

private:
    // cf. Computational Geomechanics pp.62 WW
    double* DynamicDamping;
    /// For GMRES solver. 30.06.2010. WW
    long m_cols;
    FiniteElement::ErrorMethod _pcs_nls_error_method;
    FiniteElement::ErrorMethod _pcs_cpl_error_method;

#ifdef USE_PETSC
    std::string lsover_name;  // WW
    std::string pres_name;
#endif
};

extern std::vector<CNumerics*> num_vector;
extern bool NUMRead(std::string);
extern void NUMWrite(std::string);
extern void NUMDelete();
extern CNumerics* NUMGet(std::string);

//////////////////////////////////////////////////////////////////////////
// SOLVER
//////////////////////////////////////////////////////////////////////////
struct LINEAR_SOLVER_PROPERTIES
{
    char* name;
    long type;
    long maxiter;
    double eps;
    long nom;
    long precond;
    long store;
    long criterium;
    double theta;
    long repeat;
    double time;
    long kind;
    long level;
};

struct LINEAR_SOLVER
{
    char pcs_type_name[80];
    void* matrix;
    double* b;
    double* x;
    long dim;
    long store;
    LINEAR_SOLVER_PROPERTIES* lsp;
    char* lsp_name;
    double* xx;
    double* r;
    double* memory;
    int memory_number;
    double** new_memory;
    void (*init_function)();
    void (*assemble_function)(double*, double*, double);
    long assemble_index;
    long level;
    IntFuncDXDXL LinearSolver;
    long master_iter;
    long iter_sum;
    char* system_time_name;
    char* system_time_assemble_function_name;
    char* name_group_ls;
    char* name_ls;
    long number_ls;
    long num_of_unknowns_ls;
    // OK UNKNOWN_LINEAR_SOLVER **unknown_ls;
    int unknown_vector_dimension; /* nodal degree of freedom */
    int*
        unknown_vector_indeces; /* pointer of field
                                   unknown_vector_index[unknown_vector_dimension]
                                 */
    long*
        unknown_node_numbers; /* pointer of field
                                 unknown_node_numbers[unknown_vector_dimension]
                               */
    int*
        unknown_update_methods; /* pointer of field
                                   unknown_update_methods[unknown_vector_dimension]
                                 */
};

#ifdef USE_MPI  // WW
extern LINEAR_SOLVER* InitVectorLinearSolver(LINEAR_SOLVER*);
#endif
//#ifndef NEW_EQS                                   //WW 07.11.2008
#if !defined(NEW_EQS) && \
    !defined(            \
        USE_PETSC)  // && !defined(other parallel solver) //WW. 04.10.2012

extern LINEAR_SOLVER* InitLinearSolver(LINEAR_SOLVER*);
//
extern void SetLinearSolverType(LINEAR_SOLVER*, CNumerics*);
extern LINEAR_SOLVER* InitializeLinearSolver(LINEAR_SOLVER*, CNumerics*);
extern LINEAR_SOLVER* InitMemoryLinearSolver(LINEAR_SOLVER*, int);
extern LINEAR_SOLVER* SetMemoryZeroLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER* SetZeroLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER* DestroyLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER* DestroyMemoryLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER* SetLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER* CreateLinearSolver(long store, long dim);
extern LINEAR_SOLVER* CreateLinearSolverDim(long store,
                                            int unknown_vector_dimension,
                                            long dim);
extern void ConfigSolverProperties(void);
//
//
extern int GetUnknownVectorDimensionLinearSolver(LINEAR_SOLVER*);
#endif  // ifndef NEW_EQS //WW 07.11.2008

//////////////////////////////////////////////////////////////////////////
// NUM
//////////////////////////////////////////////////////////////////////////
extern double GetNumericalTimeCollocation(char* name);
extern int GetNumericsGaussPoints(int element_dimension);
extern double NUMCalcIterationError(double* new_iteration,
                                    double* old_iteration, double* reference,
                                    long length, int method);
#endif
