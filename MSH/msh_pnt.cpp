/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
#include <string>
#include <vector>
using namespace std;

/**************************************************************************
   MshLib-Method:  tricircumcenter3d()   Find the circumcenter of a triangle in 3D
   Task: The result is returned both in terms of xyz coordinates
      coordinates, relative to the triangle's point `a' (that is, `a' is
      the origin of both coordinate systems).  Hence, the xyz coordinates
      returned are NOT absolute; one must add the coordinates of `a' to
      find the absolute coordinates of the circumcircle.  However, this means
      that the result is frequently more accurate than would be possible if
      absolute coordinates were returned, due to limited floating-point
      precision.  In general, the circumradius can be computed much more
      accurately.
   Programing:
   05/2004 TK Implementation
**************************************************************************/
void TriCircumcenter3d(double a[3],double b[3],double c[3],double circumcenter[3])
{
	double xba, yba, zba, xca, yca, zca;
	double balength, calength;
	double xcrossbc, ycrossbc, zcrossbc;
	double denominator;
	double xcirca, ycirca, zcirca;

	/* Use coordinates relative to point `a' of the triangle. */
	xba = b[0] - a[0];
	yba = b[1] - a[1];
	zba = b[2] - a[2];
	xca = c[0] - a[0];
	yca = c[1] - a[1];
	zca = c[2] - a[2];
	/* Squares of lengths of the edges incident to `a'. */
	balength = xba * xba + yba * yba + zba * zba;
	calength = xca * xca + yca * yca + zca * zca;

	/* Cross product of these edges. */
#ifdef EXACT
	/* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
	/*   to ensure a correctly signed (and reasonably accurate) result, */
	/*   avoiding any possibility of division by zero.                  */
	xcrossbc = orient2d(b[1], b[2], c[1], c[2], a[1], a[2]);
	ycrossbc = orient2d(b[2], b[0], c[2], c[0], a[2], a[0]);
	zcrossbc = orient2d(b[0], b[1], c[0], c[1], a[0], a[1]);
#else
	/* Take your chances with floating-point roundoff. */
	xcrossbc = yba * zca - yca * zba;
	ycrossbc = zba * xca - zca * xba;
	zcrossbc = xba * yca - xca * yba;
#endif

	/* Calculate the denominator of the formulae. */
	denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
	                     zcrossbc * zcrossbc);

	/* Calculate offset (from `a') of circumcenter. */
	xcirca = ((balength * yca - calength * yba) * zcrossbc -
	          (balength * zca - calength * zba) * ycrossbc) * denominator;
	ycirca = ((balength * zca - calength * zba) * xcrossbc -
	          (balength * xca - calength * xba) * zcrossbc) * denominator;
	zcirca = ((balength * xca - calength * xba) * ycrossbc -
	          (balength * yca - calength * yba) * xcrossbc) * denominator;
	circumcenter[0] = xcirca + a[0];
	circumcenter[1] = ycirca + a[1];
	circumcenter[2] = zcirca + a[2];
}

