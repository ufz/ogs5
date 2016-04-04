/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "math.h"
#include "msh_faces.h"
#include "stdlib.h"

// MathLib
#include "MathTools.h"

/*-------------------------------------------------------------------------
   Constructor and Destructor of the class CECLIPSEData
   -------------------------------------------------------------------------*/
CPlaneEquation::CPlaneEquation(void)
{
}

CPlaneEquation::~CPlaneEquation(void)
{
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CalculatePlaneEquationFrom3Points
   Task: Calculates the plane equations as well as the normal vector of the face
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CPlaneEquation::CalculatePlaneEquationFrom3Points(const double Point1[3],
                                                       const double Point2[3],
                                                       const double Point3[3])
{
	// Set the initial point
	this->Point[0] = Point1[0];
	this->Point[1] = Point1[1];
	this->Point[2] = Point1[2];

	// Calculate vectors
	this->vector1[0] = Point2[0] - Point1[0];
	this->vector1[1] = Point2[1] - Point1[1];
	this->vector1[2] = Point2[2] - Point1[2];

	this->vector2[0] = Point3[0] - Point1[0];
	this->vector2[1] = Point3[1] - Point1[1];
	this->vector2[2] = Point3[2] - Point1[2];

	// Calculate Normal vector
	this->normal_vector[0] = this->vector1[1] * this->vector2[2] - this->vector1[2] * this->vector2[1];
	this->normal_vector[1] = this->vector1[2] * this->vector2[0] - this->vector1[0] * this->vector2[2];
	this->normal_vector[2] = this->vector1[0] * this->vector2[1] - this->vector1[1] * this->vector2[0];

	// norm the normal vector to 1
	double norm = sqrt(MathLib::scpr(this->normal_vector, this->normal_vector, 3));
	this->normal_vector[0] = 1. / norm * this->normal_vector[0];
	this->normal_vector[1] = 1. / norm * this->normal_vector[1];
	this->normal_vector[2] = 1. / norm * this->normal_vector[2];

	// cout << Point1[0] << ", " << Point1[1] << ", " << Point1[2] << " " << Point2[0] << ", " << Point2[1] << ", " <<
	// Point2[2] << " " << Point3[0] << ", " << Point3[1] << ", " << Point3[2];
	// Calculate Lambda for Normal equation of the plane (= scalar product of normal vector and one point of the plane)
	this->Lambda_NormalEquation = this->Point[0] * this->normal_vector[0] + this->Point[1] * this->normal_vector[1]
	                              + this->Point[2] * this->normal_vector[2];
}
/*-------------------------------------------------------------------------
   GeoSys - Function: CheckIfPointInPlane
   Task: Check if a point lies within the plane
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CPlaneEquation::CheckIfPointInPlane(const double Point[3])
{
	(void)Point;
	// the point is included in the normal equation of the plane, it is part of the plane if the equation is fullfilled
	if (this->Lambda_NormalEquation
	    == this->Point[0] * this->normal_vector[0] + this->Point[1] * this->normal_vector[1]
	           + this->Point[2] * this->normal_vector[2])
		return 1;
	else
		return 0;
}

CFaces::CFaces(int number_phases)
{
	CFlowData* m_FlowData;

	this->gravity_centre[0] = 0.;
	this->gravity_centre[1] = 0.;
	this->gravity_centre[2] = 0.;
	this->nnodes = 0;
	// this->phases.resize(number_phases);
	for (int i = 0; i < number_phases; i++)
	{
		m_FlowData = new CFlowData;
		this->phases.push_back(m_FlowData);
		// this->phases[i]->q_norm = 0.0;
		// this->phases[i]->q[0] = 0.0;
		// this->phases[i]->q[1] = 0.0;
		// this->phases[i]->q[2] = 0.0;
	}
}

CFaces::~CFaces(void)
{
}
/*-------------------------------------------------------------------------
   GeoSys - Function: FaceCentre
   Task: Calculates the coordinates of the centre of each face of a EclipseBlock; Assumes at the moment regular
   hexahedron
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CFaces::Calculate_FaceGravityCentre(const double Point1[3],
                                         const double Point2[3],
                                         const double Point3[3],
                                         const double Point4[3])
{
	// Order of the faces: 0..left(x); 1..right(x); 2..front(y); 3..back(y); 4..bottom(z); 5..top(z)
	if (this->connected_nodes.size() > 4)
	{
		std::cout << "Error: The face has more than 4 corner points!"
		          << "\n";
		exit(1);
	}

	// FaceGravityCentre is calculated for quadrilateral element by calculating the gravitycentre for each triangle
	// the gravity-centre of the quadrilateral is the area weighted average of both centres of the triangles
	// https://lists.cs.columbia.edu/pipermail/acis-alliance/2003-September/000171.html
	// http://www.matheboard.de/archive/51412/thread.html
	// Centroid of the first triangle
	double xc1 = (Point1[0] + Point2[0] + Point3[0]) / 3.;
	double yc1 = (Point1[1] + Point2[1] + Point3[1]) / 3.;
	double zc1 = (Point1[2] + Point2[2] + Point3[2]) / 3.;
	//	vector point1-point2
	double vec_a[3];
	vec_a[0] = Point2[0] - Point1[0];
	vec_a[1] = Point2[1] - Point1[1];
	vec_a[2] = Point2[2] - Point1[2];
	//	vector point1-point3
	double vec_b[3];
	vec_b[0] = Point3[0] - Point1[0];
	vec_b[1] = Point3[1] - Point1[1];
	vec_b[2] = Point3[2] - Point1[2];

	// Area of the first triangle = half of the cross product of the 2 vectors
	double a1 = 0.5 * sqrt((vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]) * (vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1])
	                       + (vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2]) * (vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2])
	                       + (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]) * (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]));

	// Centroid of the second triangle

	///�nderung Bastian Point1 -> Point2 ge�ndert
	double xc2 = (Point1[0] + Point4[0] + Point3[0]) / 3.;
	double yc2 = (Point1[1] + Point4[1] + Point3[1]) / 3.;
	double zc2 = (Point1[2] + Point4[2] + Point3[2]) / 3.;
	//	vector point1-point4
	double vec_c[3];
	vec_c[0] = Point4[0] - Point1[0];
	vec_c[1] = Point4[1] - Point1[1];
	vec_c[2] = Point4[2] - Point1[2];

	// Area of the second triangle
	double a2 = 0.5 * sqrt((vec_c[1] * vec_b[2] - vec_c[2] * vec_b[1]) * (vec_c[1] * vec_b[2] - vec_c[2] * vec_b[1])
	                       + (vec_c[2] * vec_b[0] - vec_c[0] * vec_b[2]) * (vec_c[2] * vec_b[0] - vec_c[0] * vec_b[2])
	                       + (vec_c[0] * vec_b[1] - vec_c[1] * vec_b[0]) * (vec_c[0] * vec_b[1] - vec_c[1] * vec_b[0]));

	// Centroid of the 4-sided polygon
	if ((a1 + a2) > 0)
	{
		this->gravity_centre[0] = (xc1 * a1 + xc2 * a2) / (a1 + a2);
		this->gravity_centre[1] = (yc1 * a1 + yc2 * a2) / (a1 + a2);
		this->gravity_centre[2] = (zc1 * a1 + zc2 * a2) / (a1 + a2);
		this->face_area = a1 + a2;
		return true;
	}
	else
		return false;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: SetNodes
   Task: Sets the nodes which form the face
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CFaces::SetNodes(MeshLib::CNode* Point1, MeshLib::CNode* Point2, MeshLib::CNode* Point3, MeshLib::CNode* Point4)
{
	this->connected_nodes.push_back(Point1);
	this->connected_nodes.push_back(Point2);
	this->connected_nodes.push_back(Point3);
	this->connected_nodes.push_back(Point4);
	this->nnodes = 4;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: SetElements
   Task: Sets the elements connceted to the face
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CFaces::SetElements(std::vector<long> element_indices)
{
	this->connected_elements = element_indices;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: Calculate_q_along_axis
   Task: Calculates the x-, y and z-components of a norm of a given vector with the normal vector
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CFaces::Calculate_components_of_a_vector(int flag, int phase_index, bool Radialmodell)
{
	double* normal_vector;
	normal_vector = this->PlaneEquation->GetNormalVector();
	// contoll that normalvector is positiv in all directions, is only valid if it is no radial model
	if (Radialmodell == false)
	{
		if ((this->model_axis == "I+") || (this->model_axis == "I-"))
			if (normal_vector[0] < 0.0)
			{
				std::cout << "Error at the normal vector. The i-component is not positiv!"
				          << "\n";
				// system("Pause");
				exit(1);
			}
		if ((this->model_axis == "J+") || (this->model_axis == "J-"))
			if (normal_vector[1] < 0.0)
			{
				std::cout << "Error at the normal vector. The j-component is not positiv!"
				          << "\n";
				// system("Pause");
				exit(1);
			}
	}
	if ((this->model_axis == "K+") || (this->model_axis == "K-"))
		if (normal_vector[2] < 0.0)
		{
			std::cout << "Error at the normal vector. The k-component is not positiv!"
			          << "\n";
			// system("Pause");
			exit(1);
		}

	for (int i = 0; i < 3; i++)
	{
		// cout << this->q[i] << " " << normal_vector[i] << " ";
		if (flag == 0)
			this->phases[phase_index]->q[i] = this->phases[phase_index]->q_norm * normal_vector[i];
		if (flag == 1)
			this->vel[i] = this->v_norm * fabs(normal_vector[i]);
		// cout << this->q[i] << "\n";
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFace
   Task: Creates the face and its equations, necessary for coupling with Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CFaces::CreateFace(MeshLib::CNode* Point1, MeshLib::CNode* Point2, MeshLib::CNode* Point3, MeshLib::CNode* Point4)
{
	// Set nodes of the face
	this->SetNodes(Point1, Point2, Point3, Point4);

	// Calculate the plane equation
	double const* const coord_Point1(this->connected_nodes[0]->getData());
	double const* const coord_Point2(this->connected_nodes[1]->getData());
	double const* const coord_Point3(this->connected_nodes[2]->getData());
	double const* const coord_Point4(this->connected_nodes[3]->getData());

	this->PlaneEquation = new CPlaneEquation();
	this->PlaneEquation->CalculatePlaneEquationFrom3Points(coord_Point1, coord_Point2, coord_Point3);

	// check if 4. point lies within the plane
	if (this->PlaneEquation->CheckIfPointInPlane(coord_Point4) == false)
		return false;

	// Calculate gravity centre of the face
	// the order of the point was changed to calculate the correct gravity centre
	if (this->Calculate_FaceGravityCentre(coord_Point1, coord_Point2, coord_Point4, coord_Point3) == false)
		return false;

	return true;
}
