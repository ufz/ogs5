/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// Classes for interface GeoSys - ECLIPSE
#include <vector>
//#define MAX_ZEILEN 512

class CRFProcess;
#include "msh_faces.h"
#include "rf_pcs.h"

typedef struct
{
    double Number;
    int Exponent;
} typeExponentialNumber;
// KB: neue Struktur zum Einlesen der Brunnendatei, siehe
// CECLIPSEData::ReadWellData()
struct structWell /* deklariert den Strukturtyp well, wird beim Einlesen von
                     .well-File neu aufgerufen */
{
    std::string name;
    std::vector<double> time;
    std::vector<std::string> rate;
    std::string phase;
    std::string open_flag;
    std::string control_mode;
    // std::string dummy_rate;
    // std::string dummy_zeile;
    // std::string rest_zeile;
};

class CECLIPSEBlock
{
public:
    long index;
    long row;
    long column;
    long layer;
    std::vector<double> x_coordinates;
    std::vector<double> y_coordinates;
    std::vector<double> z_coordinates;
    double x_barycentre;
    double y_barycentre;
    double z_barycentre;
    int active;
    double volume;
    std::vector<long> connected_faces;

    std::vector<long> NeighbourElement;
    std::vector<long> ConnectedBoundaryCondition;

    CECLIPSEBlock(long Nodelength, long Facelength);
    ~CECLIPSEBlock();

    // void CalcBarycentre(void);
    void CalculateFaceCentres(void);
};

class CReadTextfiles_ECL
{
public:
    std::vector<std::string> Data;
    std::vector<std::vector<std::string> > Data_separated;
    long NumberOfRows;
    std::vector<std::string> SplittedString;
    std::vector<std::string> Header;

    CReadTextfiles_ECL();   // Konstruktor
    ~CReadTextfiles_ECL();  // Desturktor

    bool Read_Text(std::string Filename);

    void SplitStrings(const std::string str, std::string delimiter);

    bool Read_SeparatedText(std::string Filename, std::string delimiter);

    // vector <string> Header;
    //   vector <string> Data;
    //   long NumberOfRows;

    // CReadTextfiles();		//Konstruktor
    //~CReadTextfiles();		//Desturktor

    // bool Read_Text(std::string Filename);
};

class CWriteTextfiles_ECL
{
public:
    CWriteTextfiles_ECL();   // Konstruktor
    ~CWriteTextfiles_ECL();  // Desturktor

    void Write_Text(std::string Filename, std::vector<std::string> Text);
};

class CPointData_ECL
{
public:
    double x;
    double y;
    double z;
    // double Flow[3];
    double pressure;
    double temperature;
    // double Gas_dissolved;
    double CO2inLiquid;
    double NaClinLiquid;
    double deltaDIC;
    double deltaSat;
    double deltaPress;
    double VaporComponentMassFraction;
    std::vector<double> phase_pressure;
    std::vector<double> phase_saturation;
    std::vector<double> phase_density;
    std::vector<std::vector<double> > q;

    CPointData_ECL()
    {
        x = 0;
        y = 0;
        z = 0;
    }
    CPointData_ECL(double const* const pnt)
    {
        x = pnt[0];
        y = pnt[1];
        z = pnt[2];
    }
    ~CPointData_ECL() {}
};

class CBoundaryConditions
{
public:
    int index;
    int number;
    long connected_element;
    std::string boundary_position;
    double value[4];
    CBoundaryConditions()
    {
        index = -1;
        number = -1;
        connected_element = -1;
        boundary_position = "";
    }
    ~CBoundaryConditions() {}
};

class CECLIPSEData
{
public:
    long elements;
    long rows;
    long columns;
    long layers;
    long times;
    int numberOutputParameters;
    int activeCells;
    bool Radialmodell;
    bool Radial_J;
    bool Radial_I;
    bool RadialModellIpos;
    bool RadialModellJpos;
    double Molweight_CO2;   // [g/mol]
    double Molweight_H2O;   // [g/mol]
    double Molweight_NaCl;  // [g/mol]
    double SurfaceCO2Density;
    bool E100;
    bool phase_shift_flag;
    double sumCO2removed;
    int ProcessIndex_CO2inLiquid;
    int ProcessIndex_CO2inGas;  // KB
    double actual_time;
    bool Windows_System;
    bool existWells;
    bool UsePrecalculatedFiles;
    bool UseSaveEclipseDataFiles;
    bool TempIncludeFile;  // WTP bool flag to select temperature exchange with
                           // ecl and geosys
    std::string
        dissolved_co2_pcs_name_ECL;  // Keyword DISSOLVED_CO2_PCS_NAME, Name of
                                     // MASS_TRANSPORT Process which is
    // used to store total dissolved CO2 from ECLIPSE
    std::string dissolved_co2_ingas_pcs_name_ECL;  // KB
    std::vector<CECLIPSEBlock*> eclgrid;
    std::vector<std::string> SplittedString;
    std::vector<std::string> Variables;
    std::vector<CFaces*> faces;
    std::vector<CPointData_ECL*> NodeData;
    std::vector<CBoundaryConditions*> BC;
    std::vector<structWell*> ecl_well;

    double** Data;                       // array of points, times and variables
    std::vector<std::string> WellRates;  // KB, abspeichern der neuen Raten f�r
                                         // einen bestimmten Zeitschritt

    std::vector<long> output_x;
    std::vector<long> output_y;
    std::vector<long> output_z;
    std::vector<long> output_time;
    std::vector<long> CorrespondingEclipseElement;
    std::vector<long> CorrespondingGeosysElement;
    std::vector<std::string> Phases;
    std::vector<std::string> Components;
    long a[8][2];  // 2D Array um Keywords abzuspeichern
    std::vector<bool> eclipse_ele_active_flag;  // CB
    bool PoroPermIncludeFile;
    CECLIPSEData();
    ~CECLIPSEData();

    int GetVariableIndex(std::string Variablename);

    void SplitStrings(const std::string str, std::string delimiter);

    double Round(double Number, int Decimalplaces);

    typeExponentialNumber RoundEXP(double Number, int Decimalplaces);

    std::string AddZero(double Number, int Places, bool before);

    bool CheckIfFileExists(std::string strFilename);

    bool ReplaceASectionInFile(std::string Filename, std::string Keyword,
                               std::vector<std::string> Data,
                               bool CheckLengthOfSection);
    bool WriteIncludeFile(std::string Filename, std::string Keyword,
                          std::vector<std::string> Data, bool append);  // CB
    bool ReplaceWellRate(std::string Filename, std::string Keyword_well);

    int WriteDataBackToEclipse(CRFProcess* m_pcs, std::string projectname);

    std::string ExecuteEclipse(long Timestep, CRFProcess* m_pcs,
                               std::string folder);

    void ReadEclipseGrid(std::string Filename);

    void DetermineNeighbourElements(std::string Filename);

    bool ReadBoundaryData(int index_boundary, std::vector<std::string> Data);

    int ReadDataFromInputFile(std::string Filename);

    bool ReadPositionBoundaryCondition(std::string Filename);

    bool CorrespondingElements(void);

    bool CompareElementsGeosysEclipse(void);

    // TF commented out method since we have already functions for computing the
    // distance between points double CalculateDistanceBetween2Points(double
    // Point1[3], double Point2[3]);

    bool CreateFaces(void);

    bool ConnectFacesToElements(void);

    // bool MakeNodeVector(CRFProcess *m_pcs, std::string path, int timestep,
    // int phase_index);
    bool MakeNodeVector(void);

    void ReadEclipseData(std::string Pathname, long timestep);

    void CalculateRSfromMassFraction_E300();

    bool GetFlowForFaces(int phase_index);

    bool GetVelForFaces(void);

    bool CalcBlockBudget(int phase_index);

    void InterpolateDataFromFacesToNodes(long ele_nr, double* n_vel_x,
                                         double* n_vel_y, double* n_vel_z,
                                         int phase_index);

    void InterpolateDataFromBlocksToNodes(CRFProcess* m_pcs, std::string path,
                                          int phase_index);

    void InterpolateGeosysVelocitiesToNodes(CRFProcess* m_pcs, double* vel_nod,
                                            long node);

    void WriteDataToGeoSys(CRFProcess* m_pcs, std::string folder);

    void SaveEclipseDataFile(
        long Timestep,
        CRFProcess* m_pcs);  // WTP function to save a copy of the .data file

    bool CleanUpEclipseFiles(std::string folder, std::string projectname);

    int RunEclipse(long Timestep, CRFProcess* m_pcs);

    void ReadWellData(std::string Filename_Wells);

    void WriteOutput_2DSection(long Timestep);
};
