/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: Initial Conditions IC
   Task: class implementation
   Programing:
   08/2004 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_ic_new_INC
#define rf_ic_new_INC

#define IC_FILE_EXTENSION ".ic"

// C++ STL
//#include <fstream>
//#include <string>
#include <vector>

// FEM
#include "DistributionInfo.h"    // TF
#include "GeoInfo.h"             // TF
#include "LinearFunctionData.h"  // TF
#include "ProcessInfo.h"         // KR

//#include "rf_pcs.h"

class CNodeValue;

namespace GEOLIB
{
class GEOObjects;
}
namespace MeshLib
{
class CFEMesh;
}
class InitialCondition;

/**
 * class for handling initial conditions
 */
class CInitialCondition : public ProcessInfo,
                          public GeoInfo,
                          public DistributionInfo
{
private:
    size_t SubNumber;               // WW
    std::vector<int> subdom_index;  // WW
    std::vector<double> subdom_ic;  // WW
    std::string fname;              // 17.11.2009. PCH

    LinearFunctionData* dis_linear_f;  // 24.8.2011. WW

    // REMOVE CANDIDATE
    std::string geo_name;   // TF 05/2010
    double geo_node_value;  // KR
public:
    const std::string& getGeoName() const { return geo_name; }  // KR
    double getGeoNodeValue() const { return geo_node_value; }   // KR
    int GetNumDom() const { return (int)subdom_index.size(); }  // WW
    int GetDomain(const int dom_index) const
    {
        return subdom_index[dom_index];
    }  // WW
    // int mat_type; //MX
    // DIS
    // KR std::vector<CNodeValue*> node_value_vector;
    void SetDomain(int);
    void SetByNodeIndex(int);  // 19.11.2009 PCH
    void SetPolyline(int);
    void SetSurface(int);
    void SetPoint(int);
    void StoreInitialValues();  // JOD 2014-11-10
    bool storeValues;
    // void SetMaterialDomain(int); //MX
    double gradient_ref_depth;
    double gradient_ref_depth_value;
    double gradient_ref_depth_gradient;
    std::string rfr_file_name;  // OK
    CInitialCondition();
    CInitialCondition(const InitialCondition* ic);
    ~CInitialCondition();
    /**
     * read initial condition from stream
     * @param in input stream from file
     * @param geo_obj object of class GEOObjects that manages the geometric
     * entities
     * @param unique_name the name of the project to access the right geometric
     * entities
     * @return the new position in the stream after reading
     */
    std::ios::pos_type Read(std::ifstream* in,
                            const GEOLIB::GEOObjects& geo_obj,
                            const std::string& unique_name);
    void Write(std::fstream*) const;
    void Set(int);
    void SetEle(int);        // MX
    void SetDomainEle(int);  // MX
    LinearFunctionData* getLinearFunction() const { return dis_linear_f; }
    MeshLib::CFEMesh* m_msh;
};

class CInitialConditionGroup
{
public:
    std::string pcs_type_name;  // OK
    std::string pcs_pv_name;    // OK
    std::vector<CNodeValue*> group_vector;
};

extern std::vector<CInitialConditionGroup*> ic_group_vector;
extern std::vector<CInitialCondition*> ic_vector;
/**
 * read file that stores initial conditions
 * @param file_base_name base file name (without extension) containing the
 * initial conditions
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if initial conditions found in file, else false
 */
bool ICRead(const std::string& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);
extern void ICWrite(std::string);
extern void ICDelete();
extern CInitialCondition* ICGet(std::string);  // OK
#endif
