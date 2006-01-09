/*
 * InverseDistanceInterpolation.h
 *
 * 2012/08/15 KR Initial implementation
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef INVERSEDISTANCEINTERPOLATION_H_
#define INVERSEDISTANCEINTERPOLATION_H_

#include <limits>
#include <vector>

#include "Point.h"
#include "Matrix.h"

namespace MathLib
{

/**
 * Inverse distance interpolation, calculating distances from all domain points to
 * all supporting points.
 * The template classes T1 and T2 are expecting GEOLIB::Point Types.
 */
template<class T1, class T2>
class InverseDistanceInterpolation
{
public:
	InverseDistanceInterpolation(const std::vector<T1> domain_points,
								 const std::vector<T2> supporting_points);

	virtual ~InverseDistanceInterpolation();

	/// Returns the distance between point pnt_idx and supporting point support_idx.
	double getDistance(size_t pnt_idx, size_t support_idx) const;

	/// Returns to sum of distances of point pnt_idx to all supporting points.
	double getSumOfDistances(size_t pnt_idx) const;

	const T1 getDomainPoint(size_t pnt_idx) const { return _domain_points[pnt_idx]; };

	const T2 getSupportingPoint(size_t support_idx) const { return _supporting_points[support_idx]; };

	size_t getNDomainPoints() const { return _domain_points.size(); };

	size_t getNSupportingPoints() const { return _supporting_points.size(); };

	/// Allows the change the interpolation exponent (default value is 2).
	void setInterpolationExponent(double exponent) { _interpolation_exponent = exponent; };

private:
	/// Calculates the distances from all domain points to all supporting points.
	void calculateDistancesToSupportingPoints();

	double _interpolation_exponent;
	const std::vector<T1> _domain_points;
	const std::vector<T2> _supporting_points;
	std::vector<double> _sum_of_distances;
	Matrix<double> _distances;

};

template<class T1, class T2> InverseDistanceInterpolation<T1, T2>::InverseDistanceInterpolation(const std::vector<T1> domain_points,
	                                                       const std::vector<T2> supporting_points)
	: _interpolation_exponent(2.0),
	  _domain_points (domain_points),
	  _supporting_points (supporting_points),
	  _sum_of_distances ( std::vector<double>(_domain_points.size(), 0.0) ),
	  _distances ( Matrix<double>(_domain_points.size(), _supporting_points.size()) )
{
	calculateDistancesToSupportingPoints();
}

template<class T1, class T2> InverseDistanceInterpolation<T1, T2>::~InverseDistanceInterpolation()
{}

template<class T1, class T2> double InverseDistanceInterpolation<T1, T2>::getDistance(size_t pnt_idx, size_t support_idx) const
{
	if ( pnt_idx < _domain_points.size() && support_idx < _supporting_points.size() )
		return _distances(pnt_idx, support_idx);
	return -1;
}

template<class T1, class T2> double InverseDistanceInterpolation<T1, T2>::getSumOfDistances(size_t pnt_idx) const
{
	if ( pnt_idx < _domain_points.size() )
		return _sum_of_distances[pnt_idx];
	return -1;
}

template<class T1, class T2> void InverseDistanceInterpolation<T1, T2>::calculateDistancesToSupportingPoints()
{
	const size_t nPoints (_domain_points.size());
	for (size_t n=0; n<nPoints;n++)
	{
		const double *coords (_domain_points[n]->getData());

		double sum (0.0);
		std::vector<double> dummy;
		//DistanceToWeatherStation foo;

		double x (coords[0]);
		double y (coords[1]);

		const size_t nSupportingPoints (_supporting_points.size());
		for (size_t i=0; i<nSupportingPoints; i++) // Distance to each weather station
		{
			// ignore z-coordinate for the moment, assuming elevation is ignorable compared to x-y-distance
			double dist = std::max(DBL_MIN, sqrt( ((*_supporting_points[i])[0]-x) * ((*_supporting_points[i])[0]-x)
					  	                        + ((*_supporting_points[i])[1]-y) * ((*_supporting_points[i])[1]-y) ) ); // good old pythagoras
			dist = 1.0 / pow(dist, _interpolation_exponent);
			_distances(n,i) = dist;
			sum += dist;
		}

		this->_sum_of_distances[n] = sum;
	}

}

} // end namespace MathLib

#endif /* INVERSEDISTANCEINTERPOLATION_H_ */
