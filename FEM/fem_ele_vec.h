/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   Class element declaration
   class for finite element.
   Designed and programmed by WW, 06/2004
 */
#ifndef fem_dm_INC
#define fem_dm_INC

#include <vector>

#include "fem_ele.h"
#include "FEMEnums.h"
#include "Material/Solid/BGRaCreep.h"

namespace SolidProp
{
class CSolidProperties;
}

namespace MaterialLib
{
namespace solid
{
class BGRaCreep;
}
}
class CRFProcess;
class CFluidProperties;
class CMediumProperties;

namespace process
{
class CRFProcessDeformation;
}

namespace Math_Group
{
class Matrix;
}

namespace MeshLib
{
class CElem;
}
namespace FiniteElement
{
enum ProcessType;

// Vector for storing element values
class ElementValue_DM
{
public:
    ElementValue_DM(CElem* ele, const int NGP, bool HM_Staggered);
    ~ElementValue_DM();
    void ResetStress(bool cpl_loop);
    /// \param last_step The last time step or the end of the program.
    void Write_BIN(std::fstream& os, const bool last_step = false);
    void Read_BIN(std::fstream& is);
    void ReadElementStressASCI(std::fstream& is);
    double MeanStress(const int gp);

private:
    // Friend class
    friend class MaterialLib::solid::BGRaCreep;
    friend class SolidProp::CSolidProperties;
    friend class process::CRFProcessDeformation;
    friend class ::CMediumProperties;
    friend class CFiniteElementVec;
    Math_Group::Matrix* Stress0;  // Initial stress
    Math_Group::Matrix* Stress;
    Math_Group::Matrix* Stress_i;
    Math_Group::Matrix* Stress_j;
    Math_Group::Matrix* pStrain;
    Math_Group::Matrix* y_surface;
    // Preconsolidation pressure
    Math_Group::Matrix* prep0;
    Math_Group::Matrix* e_i;  // Void ratio
    // Variables of single yield surface model
    Math_Group::Matrix* xi;          // Rotational hardening variables
    Math_Group::Matrix* MatP;        // Material parameters
    Math_Group::Matrix* Strain_Kel;  // TN - Burgers model dev. Kelvin strain
    Math_Group::Matrix* Strain_Max;  // TN - Burgers model dev. Maxwell strain
    Math_Group::Matrix* Strain_pl;   // TN - Minkley - deviatoric plastic strain
    Math_Group::Matrix* Strain_t_ip;  // TN - introduced to save previous strain
                                      // in all integration points
    Math_Group::Matrix* e_pl;
    Math_Group::Matrix* ev_loc_nr_res;
    Math_Group::Matrix* lambda_pl;  // TN - Minkley - volumetric plastic strain,
                                    // plastic arc length, plastic multiplier
    Math_Group::Matrix* Strain;     // NB - Strain tensor for reload-feature
    // Discontinuity
    double disp_j;
    double tract_j;
    bool Localized;
    Math_Group::Matrix* NodesOnPath;
    double* orientation;

    Math_Group::Matrix* scalar_aniso_comp;  // WX:11.2011
    Math_Group::Matrix* scalar_aniso_tens;  // WX:11.2011 for aniso. plas.
};

// Derived element for deformation calculation
class CFiniteElementVec : public CElement
{
public:
    CFiniteElementVec(process::CRFProcessDeformation* dm_pcs,
                      const int C_Sys_Flad, const int order = 2);
    ~CFiniteElementVec();

    // Set memory for local matrices
    void SetMemory();

    // Compute the local finite element matrices and vectors
    void LocalAssembly(const int update);
    // Assemble local matrices and vectors to the global system
    bool GlobalAssembly();

    // Compute strains
    void ComputeStrain(const int ip);

    // Set material data
    void SetMaterial();

    // Get strain
    double* GetStrain() const { return dstrain; }

    ProcessType getCoupledFlowProcessType() const { return _flow_type; }

    //----------- Enhanced element -----------------------
    // Geometry related
    bool LocalAssembly_CheckLocalization(CElem* MElement);
    int IntersectionPoint(const int O_edge, const double* NodeA, double* NodeB);
    //----------- End of enhanced element ----------------
private:
    process::CRFProcessDeformation* pcs;
    ::CRFProcess* h_pcs;
    ::CRFProcess* t_pcs;
    // excavation
    bool excavation;  // 12.2009. WW
    //
    const int ns;  // Number of stresses components

    // Flow coupling
    ProcessType _flow_type;

    // Primary value indeces
    // Column index in the node value table
    int idx_P, idx_P0, idx_P1, idx_P1_0, idx_P2;
    int idx_T0, idx_T1;
    int idx_S0, idx_S, idx_Snw;
    int idx_pls;
    int idx_p1_ini;
    int idx_p2_ini;
    // Displacement column indices in the node value table
    int* _nodal_stress_indices;
    int* _nodal_strain_indices;

    // B matrix
    Matrix* B_matrix;
    Matrix* B_matrix_T;
    std::vector<Matrix*> vec_B_matrix;    // NW
    std::vector<Matrix*> vec_B_matrix_T;  // NW

    // Elastic constitutive matrix
    Matrix* De;
    // Consistent tangential matrix
    Matrix* ConsistDep;

    // Local matrices and vectors
    Matrix* AuxMatrix;
    Matrix* AuxMatrix2;  // NW
    Matrix* Stiffness;

    Matrix* PressureC;
    Matrix* PressureC_S;     // Function of S
    Matrix* PressureC_S_dp;  // Function of S and ds_dp
    Vec* RHS;
    // Global RHS. 08.2010. WW
    double* b_rhs;

    //------ Material -------
    SolidProp::CSolidProperties* smat;
    CFluidProperties* m_mfp;  // Fluid coupling
    // Medium property
    CMediumProperties* m_mmp;  // Fluid coupling
    double CalDensity();

    //  Stresses:
    //  s11, s22, s33, s12, s13, s23
    double* dstress;
    //  Strains:
    //  s11, s22, s33, s12, s13, s23
    double* dstrain;
    double* stress_ne;
    double* strain_ne;
    double* stress0;
    // Results, displacements
    //  u_x1, u_x2, u_x3, ..., u_xn,
    //  u_y1, u_y2, u_y3, ..., u_yn,
    //  u_z1, u_z2, u_z3, ..., u_zn
    double* Disp;

    // Temperatures of nodes
    double* nodal_dT;
    double* T1;

    double _wettingS;

    // Element value
    ElementValue_DM* eleV_DM;

    //------ Enhanced element ------
    // Jump flag of element nodes
    bool* NodesInJumpedA;
    // Regular enhanced strain matrix
    Matrix* Ge;
    // Singular enhanced strain matrix
    Matrix* Pe;
    // 2. For enhanced strain approach
    Matrix *BDG, *PDB, *DtD, *PeDe;  // For enhanced strain element
    // Additional node. Normally, the gravity center
    double* X0;
    // Normal to the discontinuity surface
    double* n_jump;
    // principle stresses
    double* pr_stress;
    // Compute principle stresses
    double ComputePrincipleStresses(const double* Stresses);
    // Compute principle stresses
    double ComputeJumpDirectionAngle(const double* Mat);
    //------ End of enhanced element ------

    // Form B matric
    void setB_Matrix(const int LocalIndex);
    // Form the tanspose of B matric
    void setTransB_Matrix(const int LocalIndex);
    //
    void ComputeMatrix_RHS(const double fkt, const Matrix* p_D);

    // Temporarily used variables
    double *Sxx, *Syy, *Szz, *Sxy, *Sxz, *Syz, *pstr;

    /// Extrapolation
    bool RecordGuassStrain(const int gp, const int gp_r, const int gp_s,
                           int gp_t);
    // Effective strain
    double CalcStrain_v();
    void ExtropolateGuassStrain();
    void ExtropolateGuassStress();

    // Compute the local finite element matrices
    void LocalAssembly_continuum(const int update);
    void LocalAssembly_EnhancedStrain(const int update);

    // Assembly local stiffness matrix
    void GlobalAssembly_Stiffness();
    void GlobalAssembly_PressureCoupling(Matrix* pCMatrix, double fct,
                                         const int phase = 0);
    void GlobalAssembly_RHS();

    //----------- Enhanced element ----------------
    void CheckNodesInJumpedDomain();
    // Compute the regular enhanced strain matrix
    void ComputeRESM(const double* tangJump = NULL);
    // Compute the singular enhanced strain matrix
    void ComputeSESM(const double* tangJump = NULL);

    friend class process::CRFProcessDeformation;

    double* _nodal_p1;
    double* _nodal_p2;

    double* _nodal_cp0;  /// capillary pressure.
    double* _nodal_dcp;  /// capillary pressure increment.

    // Auxiliary vector
    double* _nodal_S0;
    double* _nodal_S;
    double* AuxNodal1;

    // Dynamic
    // Damping parameters
    const bool dynamic;
    Matrix* Mass;  // For dynamic analysis
    int* Idx_Vel;
    // Auxiliary vector
    Vec* dAcceleration;
    const double beta2, bbeta1;
    void ComputeMass();
};
}  // namespace FiniteElement

extern std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;
#endif
