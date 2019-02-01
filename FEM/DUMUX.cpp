/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//#include "stdafx.h"
#include "DUMUX.h"
#include "rf_pcs.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>   // Datei streams
#include <iostream>  // Bildschirmausgabe
#include <sstream>
#include <sstream>  // string streams (in)
#include <string>
#include <sys/types.h>
#include <vector>
//#include "dirent.h"
#include "display.h"
#include <errno.h>
#include <math.h>
//#include "Windows.h"
#include "geo_mathlib.h"
#include "mathlib.h"
#include "rf_mmp_new.h"  // MAT
#include <ctime>
#include <sys/stat.h>  //for check if files exist

CReadTextfiles_DuMux::CReadTextfiles_DuMux(void)
{
    this->NumberOfRows = 0;
}

CReadTextfiles_DuMux::~CReadTextfiles_DuMux(void) {}

using std::cout;
using std::string;
using std::vector;
// using std:"\n";

/*-------------------------------------------------------------------------
   GeoSys - Function: SplitStrings
   Task: Separate a string with a delimiter
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CReadTextfiles_DuMux::SplitStrings(const string str, string delimiter)
{
    this->SplittedString.clear();

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        this->SplittedString.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Read_Text
   Task: Reads a textfile into a vector of strings
   Return: sucess
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CReadTextfiles_DuMux::Read_Text(string Filename)
{
    char Line[MAX_ZEILEN];
    bool error = false;
    bool abort = false;
    string tempstring;

    // .data ... provides the filename as it is necessary for C, ios::in ...
    // reads the file
    std::ifstream in_file(Filename.data(), std::ios::in);
    if (in_file.fail())
    {
        error = true;
        cout << "The file " << Filename.data() << " can not be opened!"
             << "\n";
    }
    else
    {
        error = false;
        while ((abort == false) && (in_file.eof() == false))
        {
            tempstring.clear();
            in_file.getline(Line, MAX_ZEILEN);
            tempstring = Line;
            // You could basically use this later to avoid line lenghts of
            // pre-defined length only getline(in_file,tempstring);
            if (tempstring.length() == MAX_ZEILEN - 1)
            {
                cout << " Error - increase MAX_ZEILEN in order to read ECLIPSE "
                        "data file "
                     << "\n";
                cout << " Or shorten the line in " << Filename.data() << ": "
                     << tempstring << " to " << MAX_ZEILEN << " characters "
                     << "\n";
                exit(0);
            }
            if (tempstring.compare("#STOP") != 0)
            {
                this->Data.push_back(Line);
                this->NumberOfRows = this->NumberOfRows + 1;
            }
            else
                abort = true;
        }
    }
    return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Read_SeparatedText
   Task: Reads a textfile into a vector of strings
   Return: sucess
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CReadTextfiles_DuMux::Read_SeparatedText(string Filename, string delimiter)
{
    char Line[MAX_ZEILEN];
    bool error = false;
    bool abort = false;
    string tempstring;

    // .data ... provides the filename as it is necessary for C, ios::in ...
    // reads the file
    std::ifstream in_file(Filename.data(), std::ios::in);
    if (in_file.fail())
    {
        error = true;
        cout << "The file " << Filename.data() << " can not be opened!"
             << "\n";
        // system("Pause");
    }
    else
    {
        error = false;
        this->NumberOfRows = 0;
        while ((abort == false) && (in_file.eof() == false))
        {
            tempstring.clear();
            in_file.getline(Line, MAX_ZEILEN);
            tempstring = Line;

            if ((this->NumberOfRows == 0) && (this->Header.size() == 0))
            {
                // Reading the header
                this->SplitStrings(tempstring, delimiter);
                for (int i = 0; i < int(this->SplittedString.size()); i++)
                    this->Header.push_back(this->SplittedString[i]);
                this->SplittedString.clear();
            }
            else
            {
                this->SplitStrings(tempstring, delimiter);
                // check if it is the same number of columns
                if (this->SplittedString.size() > 0)
                {
                    if (this->SplittedString.size() != this->Header.size())
                    {
                        error = true;
                        cout << "The number of columns in the textfile is not "
                                "constant!"
                             << "\n";
                        // system("Pause");
                    }
                    else
                    {
                        this->NumberOfRows = this->NumberOfRows + 1;
                        this->Data_separated.push_back(this->SplittedString);
                    }
                }
                else
                    abort = true;

                this->SplittedString.clear();
            }
        }
    }
    return error;
}

CWriteTextfiles_DuMux::CWriteTextfiles_DuMux(void) {}

CWriteTextfiles_DuMux::~CWriteTextfiles_DuMux(void) {}

/*-------------------------------------------------------------------------
   GeoSys - Function: Write_Text
   Task: Writes a vector of strings to a text file
   Return: sucess
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CWriteTextfiles_DuMux::Write_Text(string Filename, vector<string> Text)
{
    std::ofstream textfile;
    textfile.open(Filename.data(), std::ios::out);

    for (int i = 0; i < int(Text.size()); i++)
    {
        if (i == int(Text.size()) - 1)
            textfile << Text[i].data();
        else
            textfile << Text[i].data() << "\n";
    }
    textfile.close();
}

/*-------------------------------------------------------------------------
   Constructor and Destructor of the class CECLIPSEData
   -------------------------------------------------------------------------*/
CDUMUXData::CDUMUXData(void)
{
    this->Phases.clear();
    this->NodeData.clear();
    this->ProcessIndex_CO2inLiquid = -1;
    this->ProcessIndex_NaClinLiquid = -1;
    this->TotalSimulationTime = 0;
}

CDUMUXData::~CDUMUXData(void) {}

/*-------------------------------------------------------------------------
   GeoSys - Function: CheckIfFileExists
   Task: Check if the file exists
   Return: true if exist
   Programming: 11/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CDUMUXData::CheckIfFileExists(string strFilename)
{
    // code source: http://www.techbytes.ca/techbyte103.html (6.11.2009) no
    // restriction for use

    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0)
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    else
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;

    return blnReturn;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AddZero
   Task: Adds zero before or after the given number
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
string CDUMUXData::AddZero(double Number, int Places, bool before)
{
    string tempstring;
    string Result;
    std::stringstream NumberString;

    NumberString << Number;
    tempstring = NumberString.str();
    Result = tempstring;

    if (before == true)
        while (int(Result.length()) < Places)
        {
            if (Result.substr(0, 1) == "-")
                Result = "-0" + Result.substr(2, Result.length());
            else
                Result = "0" + Result;
        }
    else
        while (int(Result.length()) < Places)
        {
            if (Result.find(".") != std::string::npos)
                Result = Result + "0";
            else
                Result = Result + ".";
        }
    return Result;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: MakeNodeVector
   Task: Creates the Nodevector for storing the DUMUX data
   Return: true or false
   Programming: 09/2009 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
bool CDUMUXData::MakeNodeVector(void)
{
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    //	CNode* m_node = NULL;
    //	CPointData_DuMux* m_NodeData = NULL;
    //	m_NodeData = new CPointData_DuMux;
    //	PointDuMux* node_data (NULL);
    vector<double> temp_q;

    for (int i = 0; i < dim; i++)
        temp_q.push_back(-1.0E+99);

    if (this->NodeData.size() < 1)

        for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
        {
            // create new instance of CPointData
            //			m_NodeData = new CPointData_DuMux;
            PointDuMux* node_data(new PointDuMux(
                m_msh->nod_vector[i]->getData(),
                -1.0E+99,  // temperature
                -1.0E+99,  // CO2 in liquid
                -1.0E+99   // NaCl in liquid
                ));
            // Get the node
            //			m_node = m_msh->nod_vector[i];
            //			m_NodeData->x = m_node->X();
            //			m_NodeData->y = m_node->Y();
            //			m_NodeData->z = m_node->Z();
            //			m_NodeData->phase_pressure.resize(this->Phases.size());
            //			m_NodeData->phase_saturation.resize(this->Phases.size());
            //			m_NodeData->phase_density.resize(this->Phases.size());
            for (size_t j = 0; j < this->Phases.size(); j++)
                //				m_NodeData->q.push_back(temp_q);
                node_data->getQ().push_back(temp_q);
            // Set variable to zero
            //			m_NodeData->temperature = -1.0E+99;
            //			m_NodeData->CO2inLiquid = -1.0E+99;
            //			m_NodeData->NaClinLiquid = -1.0E+99;
            for (size_t k = 0; k < this->Phases.size(); k++)
            {
                //				m_NodeData->phase_pressure[k] = -1.0E+99;
                //				m_NodeData->phase_saturation[k] = -1.0E+99;
                node_data->getPhasePressure()[k] = -1.0E+99;
                node_data->getPhaseSaturation()[k] = -1.0E+99;
            }
            // transfer Data to node
            //			this->NodeData.push_back(m_NodeData);
            this->NodeData.push_back(node_data);
        }
    return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: WriteInputForDuMux
   Task: Writes the necessary input for the DuMux run of one timestep
   Return: nothing
   Programming: 11/2010 BG
   Modification:
   -------------------------------------------------------------------------*/
int CDUMUXData::WriteInputForDuMux(CRFProcess* m_pcs,
                                   string Folder,
                                   long Timestep)
{
    CWriteTextfiles_DuMux* TextFile;
    vector<string> vec_string;
    string tempstring;
    string DOScommand;
    std::ostringstream temp;
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    double value;
    double timestep_length;
    // CRFProcess *n_pcs = NULL;
    int indexConcentration_DIC = 0;
    // int indexConcentration_NaCl_dissolved;

    // delete old input files
    if (this->Windows_System == true)
        DOScommand = "del " + Folder + "\\dataForDumux.dat";
    else
        DOScommand = "rm " + Folder + "/dataForDumux.dat";

    if (system(DOScommand.c_str()))
        Display::DisplayMsgLn("Could not delete input file! ");

    // Read length of current timestep and recalculate it to days
    if (m_pcs->Tim->time_unit == "DAY")
        timestep_length = m_pcs->Tim->time_step_length * 60 * 60 * 24;
    else
    {
        if (m_pcs->Tim->time_unit == "MINUTE")
            timestep_length = m_pcs->Tim->time_step_length * 60;
        else
        {
            if (m_pcs->Tim->time_unit == "SECOND")
                timestep_length = m_pcs->Tim->time_step_length;
            else
            {
                cout << "This time unit was not considered yet"
                     << "\n";
                // system("Pause");
                exit(0);
            }
        }
    }

    // write timestep and length of timestep
    // header
    vec_string.clear();
    // tempstring = "Timestep Stepsize_[s] Restartfile";
    // vec_string.push_back(tempstring);
    // timestep
    temp.str("");
    temp.clear();
    temp << m_pcs->Tim->step_current - 1;
    tempstring = temp.str();
    tempstring += " ";
    temp.str("");
    temp.clear();
    temp << timestep_length;
    tempstring += temp.str();
    // define Restart name
    if (m_pcs->Tim->step_current > 1)
        // tempstring +=  " Restart";
        tempstring += " " + m_pcs->simulator_path + "_time=" +
                      this->AddZero(this->TotalSimulationTime, 3, false) +
                      "_rank=00000.drs";
    // tempstring += AddZero(m_pcs->Tim->step_current - 2,5,true) + ".txt";
    vec_string.push_back(tempstring);
    this->TotalSimulationTime += timestep_length;

    if (m_pcs->Tim->step_current > 1)
    {
        // write dissolved CO2 and NaCl concentration of each node
        // header
        tempstring = "Node X Y Z X_CO2 X_NaCl rho_liquid rho_gas";
        vec_string.push_back(tempstring);
        // data

        // get index of species concentration in nodevaluevector of this process
        indexConcentration_DIC =
            pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
                pcs_vector[this->ProcessIndex_CO2inLiquid]
                    ->pcs_primary_function_name[0]) +
            1;  // +1: new timelevel
        // indexConcentration_NaCl_dissolved =
        // pcs_vector[indexProcess_NaCl_dissolved]->GetNodeValueIndex(pcs_vector[indexProcess_NaCl_dissolved]->pcs_primary_function_name[0])
        // + 1; // +1: new timelevel
        for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
        {
            double const* const pnt(m_msh->nod_vector[i]->getData());
            temp.precision(12);
            temp.str("");
            temp.clear();
            temp << i;
            tempstring = temp.str();
            temp.str("");
            temp.clear();
            temp << pnt[0];
            tempstring += " " + temp.str();
            temp.str("");
            temp.clear();
            temp << pnt[1];
            tempstring += " " + temp.str();
            temp.str("");
            temp.clear();
            temp << pnt[2];
            tempstring += " " + temp.str();

            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value =
            // pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
            // indexConcentration_DIC) * (this->Molweight_CO2 / 1000) /
            // this->NodeData[i]->phase_density[0];
            value = pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(
                        i, indexConcentration_DIC) *
                    (this->Molweight_CO2 / 1000) /
                    this->NodeData[i]->getPhaseDensity()[0];
            // cout << i << " X_CO2 " << value << " C_CO2 " <<
            // pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
            // indexConcentration_DIC) << " Dichte " <<
            // this->NodeData[i]->phase_density[0] << " X_CO2_alt " <<
            // this->NodeData[i]->CO2inLiquid << "\n";
            temp.str("");
            temp.clear();
            temp << value;
            tempstring += " " + temp.str();

            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->NaClinLiquid;
            value = this->NodeData[i]->getNaClInLiquid();
            temp.str("");
            temp.clear();
            temp << value;
            tempstring += " " + temp.str();
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->phase_density[0];
            value = this->NodeData[i]->getPhaseDensity()[0];
            temp.str("");
            temp.clear();
            temp << value;
            tempstring += " " + temp.str();
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->phase_density[1];
            value = this->NodeData[i]->getPhaseDensity()[1];
            temp.str("");
            temp.clear();
            temp << value;
            tempstring += " " + temp.str();

            vec_string.push_back(tempstring);
        }
    }
    TextFile = new CWriteTextfiles_DuMux;
    // int position = int(path.find_last_of("\\"));
    // string path_new;
    // path_new = path.substr(0,position);
    // position = int(path_new.find_last_of("\\"));
    // path_new = path_new.substr(0,position);
    if (this->Windows_System == true)
        TextFile->Write_Text(Folder + "\\dataForDumux.dat", vec_string);
    else
        TextFile->Write_Text(Folder + "/dataForDumux.dat", vec_string);

    //--------------------------------------------------------------------------------------------
    // Test output
    MeshLib::CElem* m_ele = NULL;
    MeshLib::CNode* m_node = NULL;
    CMediumProperties* m_mat_mp = NULL;
    double node_volume;
    // int position;
    double porosity;
    double concentration_CO2_water;
    int group;
    double mass_CO2_gas, mass_CO2_water;

    vec_string.clear();

    tempstring =
        "Knoten; X; Y; Z; P1; P2; S1; S2; rho1; rho2; xCO2_liq; mCO2_gas; "
        "mCO2_liq";

    vec_string.push_back(tempstring);
    // Loop over all nodes
    for (size_t i = 0; i < this->NodeData.size(); i++)
    {
        m_node = m_msh->nod_vector[i];  // get element
        node_volume = 0;
        concentration_CO2_water = 0;
        for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
        {
            m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];
            // get the phase volume of current element elem
            group = m_ele->GetPatchIndex();
            m_mat_mp = mmp_vector[group];
            porosity = m_mat_mp->Porosity(
                m_ele->GetIndex(),
                1);  // CB Now provides also heterogeneous porosity, model 11
            node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
        }

        temp.str("");
        temp.clear();
        temp << i;
        tempstring = temp.str();
        // TF commented out since we want to use the improved PointDuMux class
        //		temp.str(""); temp.clear(); temp << this->NodeData[i]->x;
        // tempstring += "; " + temp.str(); 		temp.str(""); temp.clear();
        // temp
        // << this->NodeData[i]->y; tempstring += "; " + temp.str();
        // temp.str(""); temp.clear(); temp << this->NodeData[i]->z; tempstring
        // += "; " + temp.str();
        for (size_t k(0); k < 3; k++)
        {
            temp.str("");
            temp.clear();
            temp << (*(this->NodeData[i]))[k];
            tempstring += "; " + temp.str();
        }

        // TF commented out since we want to use the improved PointDuMux class
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_pressure[0]; tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_pressure[1]; tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_saturation[0];
        // tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_saturation[1];
        // tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_density[0]; tempstring
        //+=
        //"; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_density[1]; tempstring
        //+=
        //"; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->CO2inLiquid; tempstring += ";
        //"
        //+ temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhasePressure()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhasePressure()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseSaturation()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseSaturation()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseDensity()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseDensity()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getCO2InLiquid();
        tempstring += "; " + temp.str();

        concentration_CO2_water =
            pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(
                i, indexConcentration_DIC);
        // TF commented out since we want to use the improved PointDuMux class
        //		mass_CO2_gas = node_volume *
        // this->NodeData[i]->phase_saturation[1] *
        // this->NodeData[i]->phase_density[1];
        //		mass_CO2_water = node_volume *
        // this->NodeData[i]->phase_saturation[0] * concentration_CO2_water *
        // this->Molweight_CO2 * 0.001;
        mass_CO2_gas = node_volume *
                       this->NodeData[i]->getPhaseSaturation()[1] *
                       this->NodeData[i]->getPhaseDensity()[1];
        mass_CO2_water = node_volume *
                         this->NodeData[i]->getPhaseSaturation()[0] *
                         concentration_CO2_water * this->Molweight_CO2 * 0.001;

        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << mass_CO2_gas;
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << mass_CO2_water;
        tempstring += "; " + temp.str();

        vec_string.push_back(tempstring);
    }  // end node loop

    // Test Output
    int position;
    string aus_file;
    if (m_pcs->DuMuxData->Windows_System == true)
        position = Folder.find_last_of("\\");
    else
        position = Folder.find_last_of("/");
    string path = Folder.substr(0, position);
    if (m_pcs->DuMuxData->Windows_System == true)
        position = path.find_last_of("\\");
    else
        position = path.find_last_of("/");
    path = path.substr(0, position);
    // temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
    if (m_pcs->DuMuxData->Windows_System == true)
    {
        aus_file = path;
        aus_file += "\\CheckDataWroteIn_";
        temp.str("");
        temp.clear();
        temp << this->AddZero(Timestep, 4, true);
        aus_file += temp.str();
        aus_file += ".csv";
    }
    else
    {
        aus_file = path;
        aus_file += "/CheckDataWroteIn_";
        temp.str("");
        temp.clear();
        temp << this->AddZero(Timestep, 4, true);
        aus_file += temp.str();
        aus_file += ".csv";
    }
    std::ofstream aus(aus_file.data(), std::ios::out);
    for (unsigned int i = 0; i < vec_string.size(); i++)
        aus << vec_string[i] << "\n";
    aus.close();

    return 1;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReadDuMuxData
   Task: Read the data of the current time step of the DuMux model from output
   file Return: nothing Programming: 09/2010 BG Modification:
   -------------------------------------------------------------------------*/
void CDUMUXData::ReadDuMuxData(CRFProcess* m_pcs,
                               string Filename,
                               long Timestep)
{
    vector<string> files = vector<string>();
    CReadTextfiles_DuMux* TextFile;
    string tempstring;
    // WW bool saturation_water,  saturation_gas,  bool saturation_oil;
    clock_t start, finish;
    double time;
    int column;
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    MeshLib::CElem* m_ele = NULL;

    start = clock();

    // WW saturation_water =  saturation_gas =    saturation_oil = false;

    cout << "        ReadDuMuxData() ";

    // read dimension of the model
    m_ele = m_msh->ele_vector[0];  // get element
    dim = m_ele->GetDimension();

    // create correct filename with the timestep
    Filename = Filename.replace(
        Filename.length() - 8, 4, this->AddZero(Timestep - 1, 4, true));

    // Reads the text file
    bool Error;
    TextFile = new CReadTextfiles_DuMux;
    Error = TextFile->Read_SeparatedText(Filename, " ");
    if (Error == true)
    {
        cout << "The program is canceled"
             << "\n";
        // system("Pause");
        exit(0);
    }

    // index    x  y  z  p1   p2    S1  S2  qx1  qy1  qz1  qx2  qy2  qz2
    // X_CO2inBrine X_NaClInBrine  T rho_liquid rho_gas check the header of the
    // text file
    if (this->Phases.size() == 1)
        if ((TextFile->Header[0] != "index") || (TextFile->Header[1] != "x") ||
            (TextFile->Header[2] != "y") || (TextFile->Header[3] != "z") ||
            (TextFile->Header[4] != "p") || (TextFile->Header[5] != "qx") ||
            (TextFile->Header[6] != "qy") || (TextFile->Header[7] != "qz"))
        {
            cout << "The header of the DUMUX result file does not fit to the "
                    "definition!"
                 << "\n";
            // system("Pause");
            exit(0);
        }
    if (this->Phases.size() == 2)
        if ((TextFile->Header[0] != "index") || (TextFile->Header[1] != "x") ||
            (TextFile->Header[2] != "y") || (TextFile->Header[3] != "z") ||
            (TextFile->Header[4] != "p1") || (TextFile->Header[5] != "p2") ||
            (TextFile->Header[6] != "S1") || (TextFile->Header[7] != "S2") ||
            (TextFile->Header[8] != "qx1") || (TextFile->Header[9] != "qy1") ||
            (TextFile->Header[10] != "qz1") ||
            (TextFile->Header[11] != "qx2") ||
            (TextFile->Header[12] != "qy2") ||
            (TextFile->Header[13] != "qz2") ||
            (TextFile->Header[14] != "X_CO2inBrine") ||
            (TextFile->Header[15] != "X_NaClInBrine") ||
            (TextFile->Header[16] != "T") ||
            (TextFile->Header[17] != "rho_liquid") ||
            (TextFile->Header[18] != "rho_gas"))
        {
            cout << "The header of the DUMUX result file does not fit to the "
                    "definition!"
                 << "\n";
            // system("Pause");
            exit(0);
        }

    // Create the nodeData structure and check the coordinates of the nodes
    // between OGS and DUMUX
    if (Timestep == 1)
    {
        this->MakeNodeVector();

        // check the number of nodes
        if (int(this->NodeData.size()) != TextFile->NumberOfRows)
        {
            cout << "The number of nodes is not equal between OGS and DUMUX! "
                 << this->NodeData.size() << ", " << TextFile->NumberOfRows
                 << "\n";
            // system("Pause");
            exit(0);
        }

        for (long i = 0; i < TextFile->NumberOfRows; i++)
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			if ((atof(TextFile->Data_separated[i][1].data()) !=
            // this->NodeData[i]->x) ||
            //(atof(TextFile->Data_separated[i][2].data()) !=
            // this->NodeData[i]->y) ||
            //(atof(TextFile->Data_separated[i][3].data()) !=
            // this->NodeData[i]->z)) {
            if ((atof(TextFile->Data_separated[i][1].data()) !=
                 (*(this->NodeData[i]))[0]) ||
                (atof(TextFile->Data_separated[i][2].data()) !=
                 (*(this->NodeData[i]))[1]) ||
                (atof(TextFile->Data_separated[i][3].data()) !=
                 (*(this->NodeData[i]))[2]))
            {
                cout << "The node coordinates are not equal between OGS and "
                        "DUMUX!"
                     << "\n";
                // system("Pause");
                exit(0);
            }
    }

    // Read the DUMUX data
    for (long i = 0; i < TextFile->NumberOfRows; i++)
    {
        column = 3;
        double Multiplier = 1;
        for (int j = 0; j < int(this->Phases.size()); j++)
        {
            column += 1;
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			this->NodeData[i]->phase_pressure[j] =
            // atof(TextFile->Data_separated[i][column].data()) *
            // Multiplier;
            this->NodeData[i]->getPhasePressure()[j] =
                atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        }
        if (this->Phases.size() > 1)
            for (int j = 0; j < int(this->Phases.size()); j++)
            {
                column += 1;
                // TF commented out since we want to use the improved PointDuMux
                // class
                //				this->NodeData[i]->phase_saturation[j] =
                // atof(TextFile->Data_separated[i][column].data())
                //*
                // Multiplier;
                this->NodeData[i]->getPhaseSaturation()[j] =
                    atof(TextFile->Data_separated[i][column].data()) *
                    Multiplier;
            }
        for (int j = 0; j < int(this->Phases.size()); j++)
            for (int k = 0; k < 3; k++)
            {
                column += 1;
                // TF commented out since we want to use the improved PointDuMux
                // class
                //				this->NodeData[i]->q[j][k] =
                // atof(TextFile->Data_separated[i][column].data()) *
                // Multiplier;
                this->NodeData[i]->getQ()[j][k] =
                    atof(TextFile->Data_separated[i][column].data()) *
                    Multiplier;
            }
        column += 1;
        // TF commented out since we want to use the improved PointDuMux class
        //		this->NodeData[i]->CO2inLiquid =
        // atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        this->NodeData[i]->setCO2InLiquid(
            atof(TextFile->Data_separated[i][column].data()) * Multiplier);
        // cout << i << " " << this->NodeData[i]->CO2inLiquid << "\n";
        column += 1;
        // TF commented out since we want to use the improved PointDuMux class
        //		this->NodeData[i]->NaClinLiquid =
        // atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        this->NodeData[i]->setNaClInLiquid(
            atof(TextFile->Data_separated[i][column].data()) * Multiplier);
        column += 1;
        // TF commented out since we want to use the improved PointDuMux class
        //		this->NodeData[i]->temperature =
        // atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        this->NodeData[i]->setTemperature(
            atof(TextFile->Data_separated[i][column].data()) * Multiplier);
        column += 1;
        // TF commented out since we want to use the improved PointDuMux class
        //		this->NodeData[i]->phase_density[0] =
        // atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        this->NodeData[i]->getPhaseDensity()[0] =
            atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        column += 1;
        // TF commented out since we want to use the improved PointDuMux class
        //		this->NodeData[i]->phase_density[1] =
        // atof(TextFile->Data_separated[i][column].data()) * Multiplier;
        this->NodeData[i]->getPhaseDensity()[1] =
            atof(TextFile->Data_separated[i][column].data()) * Multiplier;
    }
    // Release Textfile object
    delete (TextFile);

    finish = clock();
    time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
    cout << "                    Time: " << time << " seconds."
         << "\n";

    //--------------------------------------------------------------------------------------------
    // Test output
    MeshLib::CNode* m_node = NULL;
    CMediumProperties* m_mat_mp = NULL;
    std::ostringstream temp;
    double mass_CO2_gas, mass_CO2_water;  // unused:, mass_CO2;
    double node_volume;
    vector<string> vec_string;
    // int position;
    double porosity;
    double concentration_CO2_water;
    // CRFProcess *n_pcs = NULL;
    int group;

    tempstring =
        "Knoten; X; Y; Z; P1; P2; S1; S2; rho1; rho2; xCO2_liq; mCO2_gas; "
        "mCO2_liq";

    vec_string.push_back(tempstring);
    // Loop over all nodes
    for (size_t i = 0; i < this->NodeData.size(); i++)
    {
        m_node = m_msh->nod_vector[i];  // get element
        node_volume = 0;
        concentration_CO2_water = 0;
        for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
        {
            m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];
            // get the phase volume of current element elem
            group = m_ele->GetPatchIndex();
            m_mat_mp = mmp_vector[group];
            porosity = m_mat_mp->Porosity(
                m_ele->GetIndex(),
                1);  // CB Now provides also heterogeneous porosity, model 11
            node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
        }

        temp.str("");
        temp.clear();
        temp << i;
        tempstring = temp.str();
        // TF commented out since we want to use the improved PointDuMux class
        //		temp.str(""); temp.clear(); temp << this->NodeData[i]->x;
        // tempstring += "; " + temp.str(); 		temp.str(""); temp.clear();
        // temp
        // << this->NodeData[i]->y; tempstring += "; " + temp.str();
        // temp.str(""); temp.clear(); temp << this->NodeData[i]->z; tempstring
        // += "; " + temp.str();
        for (size_t k(0); k < 3; k++)
        {
            temp.str("");
            temp.clear();
            temp << (*(this->NodeData[i]))[k];
            tempstring += "; " + temp.str();
        }

        // TF commented out since we want to use the improved PointDuMux class
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_pressure[0]; tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_pressure[1]; tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_saturation[0];
        // tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_saturation[1];
        // tempstring
        //+= "; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_density[0]; tempstring
        //+=
        //"; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->phase_density[1]; tempstring
        //+=
        //"; " + temp.str();
        //		temp.str(""); temp.clear(); temp.precision(12); temp <<
        // this->NodeData[i]->CO2inLiquid; tempstring += ";
        //"
        //+ temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhasePressure()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhasePressure()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseSaturation()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseSaturation()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseDensity()[0];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getPhaseDensity()[1];
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << this->NodeData[i]->getCO2InLiquid();
        tempstring += "; " + temp.str();

        // TF commented out since we want to use the improved PointDuMux class
        //		concentration_CO2_water = this->NodeData[i]->CO2inLiquid *
        // this->NodeData[i]->phase_density[0] / (this->Molweight_CO2 * 1e-3);
        //		mass_CO2_gas = node_volume *
        // this->NodeData[i]->phase_saturation[1] *
        // this->NodeData[i]->phase_density[1];
        //		mass_CO2_water = node_volume *
        // this->NodeData[i]->phase_saturation[0] * concentration_CO2_water *
        // this->Molweight_CO2 * 0.001;
        concentration_CO2_water = this->NodeData[i]->getCO2InLiquid() *
                                  this->NodeData[i]->getPhaseDensity()[0] /
                                  (this->Molweight_CO2 * 1e-3);
        mass_CO2_gas = node_volume *
                       this->NodeData[i]->getPhaseSaturation()[1] *
                       this->NodeData[i]->getPhaseDensity()[1];
        mass_CO2_water = node_volume *
                         this->NodeData[i]->getPhaseSaturation()[0] *
                         concentration_CO2_water * this->Molweight_CO2 * 0.001;

        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << mass_CO2_gas;
        tempstring += "; " + temp.str();
        temp.str("");
        temp.clear();
        temp.precision(12);
        temp << mass_CO2_water;
        tempstring += "; " + temp.str();

        vec_string.push_back(tempstring);
    }  // end node loop

    // Test Output
    int position;
    string aus_file;
    if (m_pcs->DuMuxData->Windows_System == true)
        position = Filename.find_last_of("\\");
    else
        position = Filename.find_last_of("/");
    string path = Filename.substr(0, position);
    if (m_pcs->DuMuxData->Windows_System == true)
        position = path.find_last_of("\\");
    else
        position = path.find_last_of("/");
    path = path.substr(0, position);
    // temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
    if (m_pcs->DuMuxData->Windows_System == true)
    {
        aus_file = path;
        aus_file += "\\CheckDataRedIn_";
        temp.str("");
        temp.clear();
        temp << this->AddZero(Timestep, 4, true);
        aus_file += temp.str();
        aus_file += ".csv";
    }
    else
    {
        aus_file = path;
        aus_file += "/CheckDataRedIn_";
        temp.str("");
        temp.clear();
        temp << this->AddZero(Timestep, 4, true);
        aus_file += temp.str();
        aus_file += ".csv";
    }
    std::ofstream aus(aus_file.data(), std::ios::out);
    for (unsigned int i = 0; i < vec_string.size(); i++)
        aus << vec_string[i] << "\n";
    aus.close();
}

/*-------------------------------------------------------------------------
   GeoSys - Function: WriteDataToGeoSys
   Task: Writes pressure, saturation and velocities to OGS nodes
   Return: true or false
   Programming: 08/2010 SB
   Modification:
   -------------------------------------------------------------------------*/
void CDUMUXData::WriteDataToGeoSys(CRFProcess* m_pcs)
{
    CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
    MeshLib::CElem* m_ele = NULL;
    CFiniteElementStd* fem;
    long index;
    double value = 0;
    // double n_vel_x[8], n_vel_y[8], n_vel_z[8];
    Math_Group::vec<long> nod_index(8);
    // ---- Gauss integral
    // int gp, gp_r=0, gp_s=0, gp_t=0;
    // double coef = 0.0, fkt=0.0;
    // int numberGaussPoints;
    // static double temp_vel[3]={0.0,0.0,0.0};
    // long Index;
    int i_ind[3];

    fem = m_pcs->fem;

    for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
    {
        // TF abbreviation
        PointDuMux const* const pnt_dumux(this->NodeData[i]);
        // pressure of wetting phase

        index =
            m_pcs->GetNodeValueIndex("PRESSURE1") + 1;  //+1... new time level
        // TF commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->getPhasePressure()[1] -
        // this->NodeData[i]->getPhasePressure()[0];
        value =
            pnt_dumux->getPhasePressure()[1] - pnt_dumux->getPhasePressure()[0];
        m_pcs->SetNodeValue(i, index, value);
        // transfer of velocities to OGS nodes
        index = m_pcs->GetNodeValueIndex("VELOCITY_X1");  //+1... new time level
        // TF commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->q[0][0];
        value = pnt_dumux->getQ()[0][0];
        m_pcs->SetNodeValue(i, index, value);
        // transfer of velocities to OGS nodes
        index = m_pcs->GetNodeValueIndex("VELOCITY_Y1");  //+1... new time level
        // TF commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->q[0][1];
        value = pnt_dumux->getQ()[0][1];
        m_pcs->SetNodeValue(i, index, value);
        // transfer of velocities to OGS nodes
        index = m_pcs->GetNodeValueIndex("VELOCITY_Z1");  //+1... new time level
        // TF commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->q[0][2];
        value = pnt_dumux->getQ()[0][2];
        m_pcs->SetNodeValue(i, index, value);
        index = m_pcs->GetNodeValueIndex("DENSITY1");  //+1... new time level
        // TF commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->phase_density[0];
        value = pnt_dumux->getPhaseDensity()[0];
        m_pcs->SetNodeValue(i, index, value);

        if (this->Phases.size() == 2)
        {
            // pressure of nonwetting phase
            index = m_pcs->GetNodeValueIndex("PRESSURE2") +
                    1;  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->phase_pressure[1];
            value = pnt_dumux->getPhasePressure()[1];
            m_pcs->SetNodeValue(i, index, value);
            // saturation of wetting phase
            index = m_pcs->GetNodeValueIndex("SATURATION1") +
                    1;  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->phase_saturation[0];
            value = pnt_dumux->getPhaseSaturation()[0];
            m_pcs->SetNodeValue(i, index, value);
            // saturation of nonwetting phase
            // index = m_pcs->GetNodeValueIndex("SATURATION2") + 1; //+1... new
            // time level value = this->NodeData[i]->phase_saturation[1];
            // m_pcs->SetNodeValue(i,index,value);
            // transfer of velocities to OGS nodes
            index =
                m_pcs->GetNodeValueIndex("VELOCITY_X2");  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->q[1][0];
            value = pnt_dumux->getQ()[1][0];
            m_pcs->SetNodeValue(i, index, value);
            // transfer of velocities to OGS nodes
            index =
                m_pcs->GetNodeValueIndex("VELOCITY_Y2");  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->q[1][1];
            value = pnt_dumux->getQ()[1][1];
            m_pcs->SetNodeValue(i, index, value);
            // transfer of velocities to OGS nodes
            index =
                m_pcs->GetNodeValueIndex("VELOCITY_Z2");  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->q[1][2];
            value = pnt_dumux->getQ()[1][2];
            m_pcs->SetNodeValue(i, index, value);
            index =
                m_pcs->GetNodeValueIndex("DENSITY2");  //+1... new time level
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			value = this->NodeData[i]->phase_density[1];
            value = pnt_dumux->getPhaseDensity()[1];
            m_pcs->SetNodeValue(i, index, value);
        }
    }

    // assign dissolved gas to GeoSys nodes
    CRFProcess* n_pcs = NULL;
    int indexConcentration;
    if (this->ProcessIndex_CO2inLiquid == -1)
        for (int i = 0; i < int(pcs_vector.size()); i++)
        {
            n_pcs = pcs_vector[i];
            // identify your process and store idx of pcs-vector
            // if ((n_pcs->nod_val_name_vector[0] == "C(4)") ||
            //    (n_pcs->nod_val_name_vector[0] == "CO2_w")) //
            //    "C(4)"...Phreeqc, "CO2_w"...Chemapp
            if (n_pcs->nod_val_name_vector[0] ==
                this->dissolved_co2_pcs_name_DUMUX)
                this->ProcessIndex_CO2inLiquid = i;
            if ((n_pcs->nod_val_name_vector[0] == "NaClinLiquid"))
                this->ProcessIndex_NaClinLiquid = i;
        }
    if (this->ProcessIndex_CO2inLiquid == -1)
    {
        cout << "In the model exists dissolved gas but there is no dissolved C "
                "in water defined. Please ad the mass "
                "transport for dissolved CO2!"
             << "\n";
        // system("Pause");
        exit(0);
    }

    // get index of species concentration in nodevaluevector of this process
    indexConcentration =
        pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
            pcs_vector[this->ProcessIndex_CO2inLiquid]
                ->pcs_primary_function_name[0]) +
        1;  // +1: new timelevel
    for (unsigned long i = 0;
         i < pcs_vector[this->ProcessIndex_CO2inLiquid]->nod_val_vector.size();
         i++)
    {
        // TF abbreviation
        PointDuMux const* const pnt_dumux(this->NodeData[i]);
        // recalculate dissolve gas: c_CO2 [mol/m�] = w_CO2 [kg_CO2 / kg_liquid]
        // * density_liq[kg/m�] / (Molweight_CO2 [g/mol] * 1e-3 [kg/g]) TF
        // commented out since we want to use the improved PointDuMux class
        //		value = this->NodeData[i]->CO2inLiquid *
        // this->NodeData[i]->phase_density[0] / (this->Molweight_CO2 *
        // 1e-3);
        value = pnt_dumux->getCO2InLiquid() * pnt_dumux->getPhaseDensity()[0] /
                (this->Molweight_CO2 * 1e-3);
        // cout << " c_CO2 " << value << " Dichte " <<
        // this->NodeData[i]->phase_density[0] << "\n";
        pcs_vector[this->ProcessIndex_CO2inLiquid]->SetNodeValue(
            i, indexConcentration, value);
        // cout << "Node: " << i << " Druck: " << m_pcs->GetNodeValue(i,
        // index_pressure2) << " RS: " << this->NodeData[i]->Gas_dissolved << "
        // Dichte: " << this->SurfaceCO2Density << " C(4): " << value << "\n";
        if (value < 0)
        {
            // TF commented out since we want to use the improved PointDuMux
            // class
            //			cout << "Node: " << i << " Druck: " <<
            // this->NodeData[i]->phase_pressure[1] << " X_CO2: " <<
            // this->NodeData[i]->CO2inLiquid << " Dichte: " <<
            // this->NodeData[i]->phase_density[0] << " C(4): " << value
            //<< "\n";
            cout << "Node: " << i
                 << " Druck: " << pnt_dumux->getPhasePressure()[1]
                 << " X_CO2: " << this->NodeData[i]->getCO2InLiquid()
                 << " Dichte: " << pnt_dumux->getPhaseDensity()[0]
                 << " dissolved CO2: " << value << "\n";
            cout << "  Fehler in Berechnung von DIC: " << value << "\n";
        }
    }
    // Loop over all elements to calculate gauss point velocities from node
    // values

    //????????????????????????????????
    // warum muss Koordinatenausrichtung �berpr�ft werden
    // warum ist in Cal_GP_Velocity_FM der Prozess notwendig?????
    //????????????????????????????

    for (int k = 0; k < int(this->Phases.size()); k++)
    {
        if (k == 0)
        {
            i_ind[0] = m_pcs->GetNodeValueIndex(
                "VELOCITY_X1");  // get index of velocity
            i_ind[1] = m_pcs->GetNodeValueIndex(
                "VELOCITY_Y1");  // get index of velocity
            i_ind[2] = m_pcs->GetNodeValueIndex(
                "VELOCITY_Z1");  // get index of velocity

            ////  check all possibilities for grid orientation (ccord_flag)
            // int coordinateflag = this->m_msh->GetCoordinateFlag();
            // i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY_X1")+1; // get
            // index of velocity if(coordinateflag == 11) i_ind[0] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Y1")+1; // get index of
            // velocity if(coordinateflag == 12) i_ind[0] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Z1")+1; // get index of
            // velocity i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY_Y1")+1;
            // // get index of velocity if(coordinateflag == 22) i_ind[1] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Z1")+1; // get index of
            // velocity i_ind[2] = m_pcs_fm->GetNodeValueIndex("VELOCITY_Z1")+1;
            // // get index of velocity
        }
        else
        {
            i_ind[0] = m_pcs->GetNodeValueIndex(
                "VELOCITY_X2");  // get index of velocity
            i_ind[1] = m_pcs->GetNodeValueIndex(
                "VELOCITY_Y2");  // get index of velocity
            i_ind[2] = m_pcs->GetNodeValueIndex(
                "VELOCITY_Z2");  // get index of velocity

            //  check all possibilities for grid orientation (ccord_flag)
            // int coordinateflag = this->m_msh->GetCoordinateFlag();
            // i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY_X2")+1; // get
            // index of velocity if(coordinateflag == 11) i_ind[0] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Y2")+1; // get index of
            // velocity if(coordinateflag == 12) i_ind[0] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Z2")+1; // get index of
            // velocity i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY_Y2")+1;
            // // get index of velocity if(coordinateflag == 22) i_ind[1] =
            // m_pcs_fm->GetNodeValueIndex("VELOCITY_Z2")+1; // get index of
            // velocity i_ind[2] = m_pcs_fm->GetNodeValueIndex("VELOCITY_Z2")+1;
            // // get index of velocity
        }
        if ((i_ind[0] < 0) || (i_ind[1] < 0) || (i_ind[2] < 0))
            cout << " Error - wrong index in Cal_GP_Velocity_FM "
                 << "\n";

        // Test Output
        // vector <string> vec_string;
        // string tempstring;
        // ostringstream temp;
        // vec_string.push_back("Element; X; Y; Z; v_Geosys_x; v_Geosys_y;
        // v_Geosys_z; v_DuMux_x; v_DuMux_y; v_DuMux_z");

        // Loop over all elements
        for (long i = 0; i < (long)m_msh->ele_vector.size(); i++)
        {
            m_ele = m_msh->ele_vector[i];  // get element
            if (m_ele->GetMark())          // Marked for use
            {  // Configure Element for interpolation of node velocities to GP
               // velocities
                fem->ConfigElement(m_ele);

                std::string tempstring;
                tempstring = "";
                // Test Output
                // temp.str(""); temp.clear(); temp << i; tempstring =
                // temp.str(); double* gc = m_ele->GetGravityCenter();
                // temp.str(""); temp.clear(); temp << gc[0]; tempstring += "; "
                // + temp.str(); temp.str(""); temp.clear(); temp << gc[1];
                // tempstring += "; " + temp.str(); temp.str(""); temp.clear();
                // temp << gc[2]; tempstring += "; " + temp.str();

                tempstring =
                    tempstring + fem->Cal_GP_Velocity_DuMux(i_ind, m_pcs, k);

                // vec_string.push_back(tempstring);
            }
        }  // end element loop

        // Test Output
        // int timestep = m_pcs->Tim->step_current;
        // if (timestep == 1 || timestep % 10 == 0 ) {
        //	string path;
        //	path = m_pcs->simulator_model_path;
        //	int position = path.find_last_of("\\");
        //	path = path.substr(0,position);
        //	position = path.find_last_of("\\");
        //	path = path.substr(0,position);
        //	temp.str(""); temp.clear(); temp << timestep; tempstring =
        // temp.str(); 	string aus_file = path + "\\CheckVelocity_" + Phases[k]
        // +
        //"_" + tempstring + ".csv"; 	ofstream aus;
        //	aus.open(aus_file.data(),ios::out);
        //	for (int i = 0; i < vec_string.size(); i++) {
        //		aus << vec_string[i] << "\n";
        //	}
        //	aus.close();
        //}
    }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ExecuteDuMux
   Task: starts DuMux
   Return: nothing
   Programming: 09/2010 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
void CDUMUXData::ExecuteDuMux(CRFProcess* m_pcs, string folder)
{
    string DuMuxExe;
    std::ostringstream temp;
    string tempstring;
    vector<string> vec_string;

    DuMuxExe = m_pcs->simulator_path;

    if (CheckIfFileExists(DuMuxExe) == false)
    {
        cout << "The DuMux executable could not be found! (" << DuMuxExe << ")"
             << "\n";
        // system("Pause");
        exit(0);
    }

    // make string for external program call to DuMux
    if (UsePrecalculatedFiles == false)
        if (this->Windows_System == false)
        {
            tempstring = DuMuxExe + " " + folder + "DuMux.dgf" + " " + folder +
                         "dataForDumux.dat";
            cout << tempstring << "\n";
            if (system(tempstring.c_str()))
            {
                Display::DisplayMsgLn("Warnung: DuMux doesn't run properly!!! ");
                exit(0);
            }
        }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RunEclipse
   Task: Preprocessing, running Eclipse, Postprocessing
   Return: nothing
   Programming: 09/2010 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
int CDUMUXData::RunDuMux(long Timestep, CRFProcess* m_pcs)
{
    string tempstring;
    string projectname;
    string Filename;
    string Executable_Filename;
    std::ostringstream temp;
    string DOScommand;
    string Executable_Folder;
    string folder;
    int position;
    clock_t start, finish;
    clock_t start_execute, finish_execute;
    double time;
    bool ReadPrecalculatedFiles;
    string root_folder, geosys_folder;

    this->Molweight_CO2 = 44.009;  // [g/mol]

    start = clock();

    ReadPrecalculatedFiles = false;

    Executable_Filename =
        m_pcs->simulator_model_path;  // path to eclipse input data
    if (Timestep == 1)
    {
        position = int(Filename.find_last_of("\\"));
        if (position >= 0)
            this->Windows_System = true;
        else
            this->Windows_System = false;
    }

    if (this->Windows_System == true)
        position = int(Executable_Filename.find_last_of("\\"));
    else
        position = int(Executable_Filename.find_last_of("/"));

    Executable_Folder = Executable_Filename.substr(0, position + 1);

    // Output of dissolved CO2, NaCl and timestep
    if (WriteInputForDuMux(m_pcs, Executable_Folder, Timestep) == 0)
    {
        cout << "There are problems with writing input data for DuMux! The run "
                "is terminated."
             << "\n";
        // system("Pause");
        exit(0);
    }

    Filename = m_pcs->simulator_model_path;  // path to eclipse input data
    if (this->Windows_System == true)
        position = int(Filename.find_last_of("\\"));
    else
        position = int(Filename.find_last_of("/"));
    folder = Filename.substr(0, position + 1);

    cout << "      Delete old result files "
         << "\n";

    if (ReadPrecalculatedFiles == false)
    {
        // Delete Resultfiles
        if (this->Windows_System == true)
            DOScommand = "del " + folder + "\\dataForGeoSys*.dat";
        else
            DOScommand = "rm " + folder + "dataForGeoSys*.dat";

        if (system(DOScommand.c_str()))
            Display::DisplayMsgLn("Could not delete input file! ");
    }

    // check if dumux folder is subfolder of the geosys folder
    if (Timestep == 1)
    {
        if (this->Windows_System == true)
            position = int(Filename.find_last_of("\\"));
        else
            position = int(Filename.find_last_of("/"));
        root_folder = Filename.substr(0, position);
        if (this->Windows_System == true)
            position = int(root_folder.find_last_of("\\"));
        else
            position = int(root_folder.find_last_of("/"));
        root_folder = Filename.substr(0, position);

        if (this->Windows_System == true)
            position = int(m_pcs->file_name_base.find_last_of("\\"));
        else
            position = int(m_pcs->file_name_base.find_last_of("/"));
        geosys_folder = m_pcs->file_name_base.substr(0, position);
        if (root_folder != geosys_folder)
        {
            cout << "Warning: The DuMux simulator model path is not part of "
                    "the GeoSys model path!!!"
                 << "\n";
            cout << root_folder << "\n";
            cout << geosys_folder << "\n";
            // system("Pause");
        }
    }

    cout << "      RunDuMux() called "
         << "\n";

    // Execute DuMux, try several times in case of problems finding the license
    // server
    start_execute = clock();

    // define number of phases
    if (Timestep == 1)
    {
        if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
            this->Phases.push_back("Water");
        if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
        {
            this->Phases.push_back("WATER");
            this->Phases.push_back("GAS");
        }
    }

    int number_loops = 0;
    int maximum_loops = 10;
    if (ReadPrecalculatedFiles == false)
    {
        do
        {
            cout << number_loops + 1 << ". trial"
                 << "\n";
            ExecuteDuMux(m_pcs, folder);
            // check if DuMux has run properly
            number_loops += 1;
            Filename = folder + "dataForGeoSys" +
                       AddZero(m_pcs->Tim->step_current - 1, 4, true) + ".dat";
        } while ((CheckIfFileExists(Filename) == false) &&
                 (number_loops <= maximum_loops));
    }

    finish_execute = clock();
    time = (double(finish_execute) - double(start_execute)) / CLOCKS_PER_SEC;
    cout << "\n";
    cout << "  Timestep: " << Timestep << "\n";
    cout << "        RunDUMUX() called                   Time: " << time
         << " seconds."
         << "\n";

    if (number_loops > maximum_loops)
    {
        cout << "The DuMux execution does not work after " << number_loops
             << " trials!"
             << "\n";
        // system("Pause");
        exit(0);
    }

    // Read Data for the current timestep from the output file of DuMux
    this->ReadDuMuxData(m_pcs, Filename, Timestep);

    WriteDataToGeoSys(m_pcs);

    if (ReadPrecalculatedFiles == false)
    {
        // CleanUpEclipseFiles(folder, projectname);
    }

    finish = clock();
    time = (double(finish) - double(start)) / CLOCKS_PER_SEC;

    cout << "        Time for this timestep: " << time << " seconds."
         << "\n";

    return 1;
}

// ToDO:
//	- Konstruktor nur einmal aufrufen
//	- in st-file: Masse einlesen und mit aktueller Dichte in Volumenstrom
// umrechnen
