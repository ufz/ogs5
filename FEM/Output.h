/**
 * \file FEM/Output.h
 * 05/04/2011 LB Refactoring: Moved from rf_out_new.h
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include "makros.h"
#include "DistributionInfo.h"
#include "GeoInfo.h"
#include "ProcessInfo.h"

#include <iostream>
#include <vector>

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
#include "mpi.h"
#endif

namespace MeshLib
{class CFEMesh;
}
namespace GEOLIB
{class GEOObjects;
}
class CVTK;

class COutput : public GeoInfo, public ProcessInfo, public DistributionInfo
{
	friend class LegacyVtkInterface;

public:
	COutput();
	COutput (size_t id);
	/**
	 * method initializes process and mesh attributes
	 */
	void init ();
        void CreateVTKInstance(void); //WW
	~COutput(void);

	/**
	 * scaling factor for values
	 * @param amplifier - a double value for scaling data
	 */
	void setAmplifier(double amplifier)
	{
		out_amplifier = amplifier;
	}

	CRFProcess* GetPCS(const std::string&); //OK
	CRFProcess* GetPCS();                 // 09/2010 TF
	CRFProcess* GetPCS_ELE(const std::string&); //OK

	/**
	 * checking the consistency of the output data as specified in the input file
	 * This means up to now, that data for missing processes is not written.
	 */
	void checkConsistency();              // CB (refactored by TF)
    void setInternalVarialbeNames(MeshLib::CFEMesh *msh);

	void GetNodeIndexVector(std::vector<int>&); //OK
	void SetNODFluxAtPLY();               //OK

	// ELE values
	const std::vector<std::string>& getElementValueVector () const { return _ele_value_vector; }
	//OK
	void GetELEValuesIndexVector(std::vector<int>&);

	/**
	 *
	 * @return
	 */
	const std::vector<std::string>& getRandomWalkParticleTracingValueVector() const
	{
		return
		        _rwpt_value_vector;
	}

	/**
	 * ToDo remove after transition to new GEOLIB - REMOVE CANDIDATE
	 * getGeoName returns a string used as id for geometric entity
	 * @return the value of attribute geo_name in case of
	 * geo_type_name == POLYLINE or geo_type_name = SURFACE
	 * If geo_type_name == POINT the id of the point is returned.
	 */
	const std::string& getGeoName() const; // TF 05/2010

	MeshLib::CFEMesh* getMesh ()                   // TF
	{
		return m_msh;
	}

	/**
	 * read from file stream
	 * @param in input file stream
	 * @param geo_obj object of class GEOObjects that manages the geometric entities
	 * @param unique_name the name of the project to access the right geometric entities
	 * @return the new position in the stream after reading
	 */
	std::ios::pos_type Read(std::ifstream& in, const GEOLIB::GEOObjects& geo_obj,
	                        const std::string& unique_name);

	void Write(std::fstream*);

	// TF not used (at the moment?) REMOVE CANDIDATE
	//    int GetPointClose(CGLPoint);
	void WriteTimeCurveData(std::fstream &);
	void WriteTimeCurveHeader(std::fstream &);
	void NODWriteDOMDataTEC();
	void WriteTECHeader(std::fstream&, int, std::string);
	void WriteTECNodeData(std::fstream&);
	void WriteTECElementData(std::fstream&, int);
	double NODWritePLYDataTEC(int);
	void NODWritePNTDataTEC(double, int);
	void ELEWriteDOMDataTEC();
	void WriteELEValuesTECHeader(std::fstream&);
	void WriteELEValuesTECData(std::fstream&);
	void NODWriteSFCDataTEC(int);
	void NODWriteSFCAverageDataTEC(double, int); //OK
	void WriteRFO();                      //OK
	void WriteRFOHeader(std::fstream&);   //OK
	void WriteRFONodes(std::fstream&);    //OK
	void WriteRFOElements(std::fstream&); //OK
	void WriteRFOValues(std::fstream&);   //OK
	void NODWriteLAYDataTEC(int);         //OK
	void ELEWriteSFC_TEC();               //OK
	void ELEWriteSFC_TECHeader(std::fstream&); //OK
	void ELEWriteSFC_TECData(std::fstream&); //OK
	void CalcELEFluxes();
	void ELEWritePLY_TEC();               //OK
	void ELEWritePLY_TECHeader(std::fstream&); //OK
	void ELEWritePLY_TECData(std::fstream&); //OK
	void TIMValue_TEC(double);            //OK
	void TIMValues_TEC(double tim_value[5], std::string *header, int dimension);   //BG 04/2011 added for more than 1 value per time
	double NODFlux(long);                 //OK
	void PCONWriteDOMDataTEC();           //MX
	void WriteTECNodePCONData(std::fstream &); //MX


	void NODWriteTotalFlux(double, int);	// JOD 2014-11-10
	void NODWritePointsCombined(double);	// 6/2012 JOD
	void NODWritePrimaryVariableList(double);	// JOD 2014-11-10
	void CalculateTotalFlux(MeshLib::CFEMesh*, std::vector<long>&, std::vector<double>&, std::vector<double>&); // JOD 2014-11-10
	void SetTotalFluxNodes(std::vector<long>& nodes_vector); //JOD 2014-11-10
	void SetTotalFluxNodesPLY(std::vector<long>& nodes_vector); // JOD 2014-11-10
	void SetTotalFluxNodesSURF(std::vector<long>& nodes_vector); // JOD 2014-11-10
	void SetTotalFluxNodesDOM(std::vector<long>& nodes_vector); // JOD 2014-11-10
    //------------------------------------------------------
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
	void setMPI_Info(const int rank, const int size, std::string rank_str);
	/// Head for binary output for parallel computing. 01.2014. WW
	void NODDomainWriteBinary_Header();
	/// Binary output for parallel computing. 01.2014. WW
	void NODDomainWriteBinary();
#endif
	void WritePetrelElementData(int timestepnumber);

	void setTime (double time) { _time = time; }
	/**
	 * get time returns the value of attribute time
	 * @return
	 */
	double getTime () const { return _time; }

	const std::vector<double>& getTimeVector () const { return time_vector; }
    const std::string& getFileBaseName () const { return file_base_name; }

    /**
     * @brief sets file_base_name to the full path corresponding to the given base name.
     *
     * The function internally uses the defaultOutputPath as set as a commandline argument.
     */
    void setFileBaseName(const std::string& fn);

	size_t getNSteps () const { return nSteps; }
	/**
	 * constructs/adds the output file name using geo_name,
	 * process type, mesh type
	 * @param fname a reference to the constructed file name
	 * @param geo switch on/off geo info in file name (default = on)
	 * @param process switch on/off process info in file name (default = on)
	 * @param mesh switch on/off mesh info in file name (default = on)
	 */
	void addInfoToFileName(std::string& fname, bool geo = true, bool process =
	                               true, bool mesh = true, const std::string &ext=TEC_FILE_EXTENSION) const; // 09/2010 TF

	std::vector<std::string> _nod_value_vector;
    std::vector<std::string> _alias_nod_value_vector;
	// MAT values
	std::vector<std::string> mmp_value_vector; //OK
	std::vector<std::string> mfp_value_vector; //OK

	CRFProcess* m_pcs;                    //OK

	//	std::vector<double>& getRWPTTimeVector () { return rwpt_time_vector; }
	std::vector<double>& getRWPTTimeVector () { return time_vector; }
    bool VARIABLESHARING;						 // Coordinates of each node as well as connection list is stored only for the first time step; BG: 05/2011

private:
	friend void OUTData(double, int step, bool force_output);

	//	std::vector<double> rwpt_time_vector; //JT, needed because outputs are treated differently in RWPT

	// MSH
	std::string msh_type_name;            //OK

	// TIM
	std::string tim_type_name;            // STEPS or TIMES ?
	std::vector<double> time_vector;
	double _time;

	/**
	 * the position in the global vector out_vector, used only in NODWritePLYDataTEC
	 */
	size_t _id;

	std::string file_base_name;
	double out_amplifier;                 //WW to amplify output
	                                      //WW/OK

	MeshLib::CFEMesh* m_msh;
	int nSteps;                           // After each nSteps, make output

	CVTK* vtk;
	// GEO
	/**
	 * the id of the geometric object as string REMOVE CANDIDATE
	 */
	std::string geo_name;                 // TF 05/2010

	// File status
	bool _new_file_opened;                //WW

	// DAT
	/**
	 * this attribute stores the output format
	 */
	std::string dat_type_name;

	// ELE value
	std::vector<std::string> _ele_value_vector;

	// RWPT values
	std::vector<std::string> _rwpt_value_vector;

	// PCON values
	std::vector<std::string> _pcon_value_vector;

	/// Tecplot share zone
	bool tecplot_zone_share; // 10.2012. WW

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
	int mrank;
	int msize;
	std::string mrank_str;

	int int_disp;
	MPI_Offset offset;

	unsigned domain_output_counter; // WW 04.2014

	void setDataArrayDisp();
#endif
};
#endif // OUTPUT_H
