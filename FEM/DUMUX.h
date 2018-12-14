/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// Classes for interface GeoSys - DUMUX
#include <vector>

// class CRFProcess;
//#include "rf_pcs.h"
#include "fem_ele.h"
#include "fem_ele_std.h"

#include "geo_pnt.h"

class CReadTextfiles_DuMux
{
public:
    std::vector<std::string> Data;
    std::vector<std::vector<std::string> > Data_separated;
    long NumberOfRows;
    std::vector<std::string> SplittedString;
    std::vector<std::string> Header;

    CReadTextfiles_DuMux();   // Konstruktor
    ~CReadTextfiles_DuMux();  // Desturktor

    bool Read_Text(std::string Filename);

    void SplitStrings(const std::string str, std::string delimiter);

    bool Read_SeparatedText(std::string Filename, std::string delimiter);
};

class CWriteTextfiles_DuMux
{
public:
    CWriteTextfiles_DuMux();   // Konstruktor
    ~CWriteTextfiles_DuMux();  // Desturktor

    void Write_Text(std::string Filename, std::vector<std::string> Text);
};

class PointDuMux : public GEOLIB::Point
{
public:
    PointDuMux(double x, double y, double z, double temperature,
               double CO2_in_liquid, double NaCl_in_liquid)
        : GEOLIB::Point(x, y, z),
          _temperature(temperature),
          _CO2_in_liquid(CO2_in_liquid),
          _NaCl_in_Liquid(NaCl_in_liquid)
    {
    }

    PointDuMux(double const* const coords, double temperature,
               double CO2_in_liquid, double NaCl_in_liquid)
        : GEOLIB::Point(coords),
          _temperature(temperature),
          _CO2_in_liquid(CO2_in_liquid),
          _NaCl_in_Liquid(NaCl_in_liquid)
    {
    }

    double getTemperature() const { return _temperature; }
    double getCO2InLiquid() const { return _CO2_in_liquid; }
    double getNaClInLiquid() const { return _NaCl_in_Liquid; }
    void setTemperature(double temperature) { _temperature = temperature; }
    void setCO2InLiquid(double CO2_in_liquid)
    {
        _CO2_in_liquid = CO2_in_liquid;
    }
    void setNaClInLiquid(double NaCl_in_Liquid)
    {
        _NaCl_in_Liquid = NaCl_in_Liquid;
    }
    std::vector<double>& getPhasePressure() { return _phase_pressure; }
    std::vector<double> const& getPhasePressure() const
    {
        return _phase_pressure;
    }
    std::vector<double>& getPhaseSaturation() { return _phase_saturation; }
    std::vector<double> const& getPhaseSaturation() const
    {
        return _phase_saturation;
    }
    std::vector<double>& getPhaseDensity() { return _phase_density; }
    std::vector<double> const& getPhaseDensity() const
    {
        return _phase_density;
    }
    std::vector<std::vector<double> >& getQ() { return _q; }
    std::vector<std::vector<double> > const& getQ() const { return _q; }

private:
    double _temperature;
    double _CO2_in_liquid;
    double _NaCl_in_Liquid;
    std::vector<double> _phase_pressure;
    std::vector<double> _phase_saturation;
    std::vector<double> _phase_density;
    std::vector<std::vector<double> > _q;
};

// class CPointData_DuMux {
// public:
//	double x;
//	double y;
//	double z;
//	double temperature;
//	double CO2inLiquid;
//	double NaClinLiquid;
//	std::vector <double> phase_pressure;
//	std::vector <double> phase_saturation;
//	std::vector <double> phase_density;
//	std::vector <std::vector <double> > q;
//
//	CPointData_DuMux() {x = 0; y = 0; z = 0;}
//	~CPointData_DuMux() {}
//};

class CDUMUXData
{
public:
    CDUMUXData();
    ~CDUMUXData();

    //	std::vector <CPointData_DuMux*> NodeData;
    std::vector<PointDuMux*> NodeData;
    std::vector<std::string> Phases;
    int dim;
    int ProcessIndex_CO2inLiquid;
    int ProcessIndex_NaClinLiquid;
    bool Windows_System;
    bool UsePrecalculatedFiles;
    double Molweight_CO2;  // [g/mol]
    double TotalSimulationTime;
    std::string
        dissolved_co2_pcs_name_DUMUX;  // Keyword DISSOLVED_CO2_PCS_NAME; Name
                                       // of MASS_TRANSPORT Process which is
    // used to store total dissolved CO2 from DUMUX

    // CFiniteElementStd* GetAssembler() {return fem; }

    bool CheckIfFileExists(std::string strFilename);

    std::string AddZero(double Number, int Places, bool before);

    bool MakeNodeVector(void);

    void ExecuteDuMux(CRFProcess* m_pcs, std::string folder);

    int WriteInputForDuMux(CRFProcess* m_pcs, std::string Pathname,
                           long Timestep);

    void ReadDuMuxData(CRFProcess* m_pcs, std::string Pathname, long Timestep);

    void WriteDataToGeoSys(CRFProcess* m_pcs);

    int RunDuMux(long Timestep, CRFProcess* m_pcs);
};
