/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib-Object:
   Task: Random Walk - an alternative for FDM or FEM of transport equation
   Programing:
   07/2005 PCH Implementation
**************************************************************************/
#ifndef break_RWPT                                //OK
#define break_RWPT                                //OK

#include "rf_mmp_new.h"
#include "rfmat_cp.h"

#define ALLOW_PARTICLES_GO_OUTSIDE

class Particle
{
public:
	// Position Vector
	double x;
	double y;
	double z;
	int on_boundary; //YS

	// Velocity Vector
	double Vx;
	double Vy;
	double Vz;

	// Konductivity
	double K;

	// Dispersion coefficient tensor
	double D[9];

	// Derivative of velcities
	double dVxdx;
	double dVydy;
	double dVzdz;

	// Time
	double t;
	double StartingTime;                  // JT 2010, added for continuous particle boundary, segemented by StartingTime

	// the element it belongs to
	int elementIndex;
	int edgeIndex;
	// particle identity
	int identity;

	// Constructor
	Particle(void);

	// Some operator overloading
	Particle& operator=(const Particle& B)
	{
		x = B.x;
		y = B.y;
		z = B.z;
		Vx = B.Vx;
		Vy = B.Vy;
		Vz = B.Vz;
		K = B.K;
		t = B.t;
		elementIndex = B.elementIndex;
		edgeIndex = B.edgeIndex;
		identity = B.identity;
		StartingTime = B.StartingTime;
		dVxdx = B.dVxdx;
		dVydy = B.dVydy;
		dVzdz = B.dVzdz;
        on_boundary = B.on_boundary;

		for(int i = 0; i < 9; ++i)
			D[i] = B.D[i];

		return *this;
	}
};

class Trace
{
public:
	// I am only bookkeeping the current and future positions
	// to save memory.
	Particle Past;
	Particle Now;

	// I am going to use the system default constructor and destructor.
};

class RandomWalk
{
public:
	RandomWalk(int srand_seed);
	~RandomWalk(void);

	// This will be determined from user's input
	// This may be a global variable somewhere in application-level later.
	// For now, I just created here.
	int numOfParticles;
	int UniformOrNormal;
	int leavingParticles;
	int srand_seed;
	int RWPTMode;                         // 0: Advection and dispersion for homogeneous media
	// 1: Advection and dispersion for heterogeneous media
	// 2: Advection only for homogeneous media
	// 3: Advection only for heterogeneous media
	// 4: Dispersion only for homogeneous media
	// 5: Dispersion only for heterogeneous media
	int PURERWPT;                         // 0: Defualt - Velocity solved by GeoSys
	// 1: Velocity fields on nodes are given in a separate file.
	// 2: Velocity solved as in FDM approach
	int FDMIndexSwitch;                   // 0: Build
	// 1: No need to build

	int GridOption;
	double CurrentTime;
	double* ChanceOfIrreversed;

	class FDMIndex
	{
public:
		int i;
		int j;
		int k;
		int eleIndex;

		FDMIndex(void)
		{
			i = j = k = eleIndex = -10;
		}
	};

	class Position
	{
public:
		double p[3];
	};

	class Pathline
	{
public:
		std::vector<Position> path;
	};

	std::vector<Trace> X;
	std::vector<Pathline> pathline;

	double Marsaglia(void);               // N(0,1) sample generator
	int IsTheParticleInThisElement(Particle* A);

	void InterpolateVelocityOfTheParticleByInverseDistance(Particle* A);
	void InterpolateVelocityOfTheParticleByBilinear(int option, Particle* A);
	double* InterpolateLocationOfTheParticleByBilinear(Particle* A, double dt);
	void InterpolateVelocity(Particle* A);
	void TracePathlineInThisElement(Particle* A);
	int IndexOfTheElementThatThisParticleBelong(int option, Particle* A);

	double randomMinusOneToOne(void);     // create uniform random number between -1 and 1
	double randomZeroToOne(void);         // create uniform random number between 0 and 1

	void AdvanceToNextTimeStep(double dt,double ctime);
	void AdvanceBySplitTime(double dt, int numOfSplit);
	void TraceStreamline(void);
	void GetDisplacement(Particle* B,
	                     double* Z,
	                     double* V,
	                     double* dD,
	                     double time,
	                     double* dsp);

	void RandomlyDriftAway(Particle* A, double dt, double* delta, int type);
	int RandomWalkDrift(double* Z, int type);
	void SolveDispersionCoefficient(Particle* A);
	void RandomWalkOutput(double,int);    //JT 2010

	int SolveForNextPosition(Particle* A, Particle* B);

	int SolveForTwoIntersectionsInTheElement(Particle* A, double* p1, double* p2, int axis);
	int SolveForDerivativeOfVelocity(Particle* A);
	int SolveForDisplacementByDerivativeOfDispersion(Particle* A, double* dD);
	double SolveDistanceBetweenTwoPoints(double* p1, double* p2);

	int GetTheElementOfTheParticleFromNeighbor(Particle* A);
	int GetTheElementOfTheParticle(Particle* Pold, Particle* Pnew);
	int SearchElementFromNeighbor(Particle* A);
	int CheckElementIndex(Particle* A);

	int IsParticleOnTheEdge(Particle* A);

	// Fracture Network
	void DoJointEffectOfElementInitially(void);
	void MakeRoulette(double* fit, double* roulette, int numOfCases);
	int Select(double* roulette, int numOfCases);
	int RouletteWheelSelection(double* chances, int numOfCases);

	// Transform coordinates
	void ToTheXYPlane(MeshLib::CElem* E, double* X);
	void ToTheXYPlane(int idx, double* X);
	void ToTheRealPlane(MeshLib::CElem* E, double* X);
	void ToTheRealPlane(int idx, double* X);
	void SolveAnglesOfTheElment(MeshLib::CElem* E);
	void IsoparametricMappingQuadfromPtoR(int index, double* R);
	void IsoparametricMappingQuadfromRtoP(int index, double* P);
	double TA(double Gx,double uT, double uA, double xTA);
	double TB(double Gx,double uT, double uA, double uB, double xTB);
	double Tmin(double* a, int* idx);

	// Read in velocity fields from an separate file
	int ReadInVelocityFieldOnNodes(std::string file_base_name);
	void buildFDMIndex(void);
	void RecordPath(int no, Particle* P);

	// Rate-limited reaction - sorption and desorption
	double Two_rateModel(double A, double k1, double k2, double t);
	CRFProcess *getFlowPCS() const {return flow_pcs;}

protected:
	FiniteElement::CFiniteElementStd* fem;

private:
	CRFProcess* m_pcs;
    CRFProcess* flow_pcs; //FM_TEST
	CFEMesh* m_msh;

	std::vector<FDMIndex> indexFDM;
	std::vector<std::string>rwpt_out_strings; //JT
	int nx;
	int ny;
	int nz;
	int neFDM;                            // JT
	double dx;
	double dy;
	double dz;
	double xrw_range;                     //JT
	double yrw_range;
	double zrw_range;

	double ComputeVolume(Particle* A, MeshLib::CElem* m_ele);
	double ComputeVolume(Particle* A, Particle* element, MeshLib::CElem* m_ele);
	void CopyParticleCoordToArray(Particle* A,
	                              double* x1buff,
	                              double* x2buff,
	                              double* x3buff,
	                              double* x4buff);

	void GetNodeOfMiniFEMforTheEdge(MeshLib::CNode* theNode,
	                                MeshLib::CEdge* theEdge,
	                                Particle* A);
	  void CheckBoundary2D(Particle* A, Particle* B);
	  void CheckBoundary3D(Particle* A, Particle* B);

	int G_intersect_line_segments (
	        double ax1,double ay1, double ax2,double ay2,
	        double bx1,double by1, double bx2,double by2,
	        double* ra,double* rb,
	        double* x,double* y);
	int G_intersect_line_segments_3D(
	        double* pl1,double* pl2, double* pp1,double* pp2,
	        double* pp3,double* pi);
	void ConcPTFile(const char* file_name);

	/**
	 * Select the mesh whose process name has the mesh for Fluid_Momentum
	 * @return
	 */
	CFEMesh* selectMeshForFluidMomentumProcess ();
};

extern void PCTRead(std::string);
extern void DATWriteParticleFile(int);
#endif                                            //OK
