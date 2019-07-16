/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: Boundary Conditions
   Task: class implementation
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_bc_new_INC
#define rf_bc_new_INC

//#include <list>
//#include <fstream>
//#include <string>
//#include <vector>

namespace FileIO
{
class BoundaryConditionIO;
}

// new GEOLIB
#include "DistributionInfo.h"  // TF
#include "GEOObjects.h"
#include "GeoInfo.h"             // TF
#include "LinearFunctionData.h"  // TF
#include "ProcessInfo.h"         // KR
#include "Constrained.h"
#include "SwitchBC.h"

// GEOLib
//#include "geo_ply.h"
// MSHLib
//#include "msh_lib.h"
// PCSLib
//#include "rf_pcs.h"
namespace MeshLib
{
class CFEMesh;
}

struct TimeInterval
{
    TimeInterval(const double start_time_, const double end_time_)
        : start_time(start_time_), end_time(end_time_)
    {
    }

    bool isInTimeInterval(const double time) const
    {
        if (time < start_time)
            return false;
        if (time > end_time)
            return false;

        return true;
    }

    const double start_time;
    const double end_time;
};

class BoundaryCondition;

class CBoundaryCondition : public ProcessInfo,
                           public GeoInfo,
                           public DistributionInfo
{
public:
    //	CBoundaryCondition(ProcessInfo const& process_info,
    //			GeoInfo const& geo_info,
    //			DistributionInfo const& distribution_info,
    //			);
    friend class CBoundaryConditionsGroup;
    friend class FileIO::BoundaryConditionIO;
    CBoundaryCondition();
    CBoundaryCondition(const BoundaryCondition* bc);

    ~CBoundaryCondition();
    //      void Write(std::fstream*) const;
    void WriteTecplot(std::fstream*) const;

    /**
     * reads a boundary condition from stream
     * @param in input file stream for reading
     * @param geo_obj pointer to the geometric object manager
     * @param unique_fname the project name
     * @param valid after return the variable valid contains the status of the
     * object, valid is false if there occured an error while reading the data,
     * else true
     * @return the position in the stream after the boundary condition
     */
    // TF
    std::ios::pos_type Read(std::ifstream* in,
                            const GEOLIB::GEOObjects& geo_obj,
                            const std::string& unique_fname, bool& valid);

    /**
     * ToDo remove after transition to new GEOLIB - REMOVE CANDIDATE
     * getGeoName returns a string used as id for geometric entity
     * @return the value of attribute geo_name in case of
     * geo_type_name == POLYLINE or geo_type_name = SURFACE
     * If geo_type_name == POINT the id of the point is returned.
     */
    const std::string& getGeoName() const;  // TF 05/2010

    int getCurveIndex() const  // TF 05/2010
    {
        return _curve_index;
    }

    bool isPeriodic() const  // TF 07/2010
    {
        return _periodic;
    }
    double getPeriodeTimeLength() const  // TF 07/2010
    {
        return _periode_time_length;
    }
    double getPeriodePhaseShift() const  // TF 07/2010
    {
        return _periode_phase_shift;
    }

    bool isInTimeInterval(const double time) const
    {
        // No period defined. That means the time is always in period.
        if (_time_interval == NULL)
            return true;

        return _time_interval->isInTimeInterval(time);
    }

    const std::vector<int>& getPointsWithDistribedBC() const
    {
        return _PointsHaveDistribedBC;
    }
    const std::vector<double>& getDistribedBC() const { return _DistribedBC; }
    std::vector<double>& getDistribedBC() { return _DistribedBC; }
    double getGeoNodeValue() const { return geo_node_value; }
    // KR

    const std::vector<std::string>& getPointsFCTNames() const
    {
        return _PointsFCTNames;
    }
    size_t getMeshNodeNumber() const { return _msh_node_number; }
    const std::string& getMeshTypeName() const { return _msh_type_name; }
    int getExcav() { return bcExcav; }  // WX:12.2010 get bc excav model
    int getExcavMatGr()
    {
        return MatGr;
    }  // WX:12.2010 get excav material group
    int getNoDispIncre() { return NoDispIncre; };  // WX:12.2012
    // give head bc for PRESSURE1 primary variable	//MW
    int getPressureAsHeadModel() const { return _pressure_as_head_model; }
    // return given density
    double getPressureAsHeadDensity() const
    {
        return _pressure_as_head_density;
    };
    // constrain a BC by other process
    bool isConstrainedBC() const { return _isConstrainedBC; }
    Constrained const& getConstrainedBC(std::size_t i) const
    {
        return _constrainedBC[i];
    }
    std::size_t getNumberOfConstrainedBCs() const
    {
        return _constrainedBC.size();
    }
    bool isSeepageBC() const { return _isSeepageBC; }
    bool isSwitchBC() const { return _isSwitchBC; }
    SwitchBC const& getSwitchBC() const { return _switchBC; }
private:
    std::vector<std::string> _PointsFCTNames;
    std::vector<int> _PointsHaveDistribedBC;
    std::vector<double> _DistribedBC;

    // GEO
    /**
     * the id of the geometric object as string REMOVE CANDIDATE
     */
    std::string geo_name;  // TF 05/2010
    std::string geo_type_name;

    std::string fname;  // 27.02.2009. WW
    int _curve_index;   // Time function index

    // DIS
    std::vector<long> node_number_vector;
    std::vector<double> node_value_vector;
    long geo_node_number;
    double geo_node_value;

    double _periode_phase_shift;  // JOD
    double _periode_time_length;  // JOD
    bool _periodic;               // JOD

    double gradient_ref_depth;  // 6/2012 JOD
    double gradient_ref_depth_value;
    double gradient_ref_depth_gradient;

    double node_value_cond;  // OK
    double condition;        // OK
    double epsilon;  // NW. temporally set here for surface interpolation
    bool time_dep_interpol;

    // FCT
    std::string fct_name;
    bool conditional;

    LinearFunctionData* dis_linear_f;  // 24.8.2011. WW

    // WW
    void SurfaceInterpolation(CRFProcess* m_pcs,
                              std::vector<long>& nodes_on_sfc,
                              std::vector<double>& node_value_vector);
    inline void DirectAssign(CRFProcess* m_pcs, long ShiftInNodeVector);
    // 19.03.2009. WW
    inline void PatchAssign(long ShiftInNodeVector);

    // MSH
    long _msh_node_number;
    std::string _msh_type_name;  // OK4105

    // copy values   SB 09.2012
    std::string copy_geom;
    std::string copy_geom_name;

    // Excavation WX:12.2010
    int bcExcav;
    int MatGr;
    // no displacement increment 12.2012
    int NoDispIncre;
    // give head bc for PRESSURE1 primary variable	//MW
    int _pressure_as_head_model;
    // given density for pressure_as_head BC
    double _pressure_as_head_density;
    // constrain a BC by other process
    bool _isConstrainedBC;
    std::vector<Constrained> _constrainedBC;
    bool _isSeepageBC;

    bool _isSwitchBC;
    SwitchBC _switchBC;

    TimeInterval* _time_interval;
};

class CBoundaryConditionNode  // OK raus
{
public:
    long geo_node_number;
    long msh_node_number;
    long msh_node_number_subst;  // WW

    double node_value;
    double node_value_offset;
    double node_value_pre_calc;
    int CurveIndex;           // Time dependent function index
    std::string pcs_pv_name;  // YD/WW
    //
    std::string fct_name;  // WW
    // FCT
    int conditional;  // OK
    std::string bc_node_copy_geom;
    std::string bc_node_copy_geom_name;
    CBoundaryConditionNode();

    void SetNormalVector(double const* const normal_vector);
    double const* GetNormalVector() const;

    // 25.08.2011. WW
    void Read(std::istream& is);
    void Write(std::ostream& os) const;

private:
    double _normal_vector[3];
};

class CBoundaryConditionsGroup
{
public:
    CBoundaryConditionsGroup(void);
    ~CBoundaryConditionsGroup(void);

    void Set(CRFProcess* pcs, int ShiftInNodeVector, const double value_offset,
             const std::string& this_pv_name = "");
    CBoundaryConditionsGroup* Get(const std::string&);

    const std::string& getProcessTypeName() const { return _pcs_type_name; }
    void setProcessTypeName(const std::string& pcs_type_name)
    {
        _pcs_type_name = pcs_type_name;
    }
    const std::string& getProcessPrimaryVariableName() const
    {
        return _pcs_pv_name;
    }
    void setProcessPrimaryVariableName(const std::string& pcs_pv_name)
    {
        if (_pcs_type_name.find("MASS_TRANSPORT") == std::string::npos)
            _pcs_pv_name = pcs_pv_name;
        else
            _pcs_pv_name = "CONCENTRATION1";
    }
    long msh_node_number_subst;  // WW
    std::string fct_name;        // OK

    MeshLib::CFEMesh* m_msh;  // OK
    // WW std::vector<CBoundaryCondition*>bc_group_vector; //OK
    // WW double GetConditionalNODValue(int,CBoundaryCondition*); //OK
    int time_dep_bc;

private:
    std::string group_name;
    std::string _pcs_type_name;  // OK
    std::string _pcs_pv_name;    // OK
};

//========================================================================
#define BC_FILE_EXTENSION ".bc"
extern std::list<CBoundaryConditionsGroup*> bc_group_list;
extern CBoundaryConditionsGroup* BCGetGroup(const std::string& pcs_type_name,
                                            const std::string& pcs_pv_name);
extern std::list<CBoundaryCondition*> bc_list;

/**
 * read boundary conditions from file
 * @param file_base_name the base name of the file (without extension)
 * @param geo_obj the geometric object managing geometric entities
 * @param unique_name the (unique) name of the project
 * @return false, if the file can not opened, else true
 */
bool BCRead(std::string const& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

extern void BCWrite(std::string const&);
extern void BCDelete();
extern void BCGroupDelete(const std::string& pcs_type_name,
                          const std::string& pcs_pv_name);
extern void BCGroupDelete(void);
// OK
extern CBoundaryCondition* BCGet(const std::string&, const std::string&,
                                 const std::string&);
extern CBoundaryCondition* BCGet(std::string);  // OK

// ToDo
extern void ScalingDirichletBoundaryConditions(const double factor);
#endif
