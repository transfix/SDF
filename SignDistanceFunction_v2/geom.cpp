/*
  Copyright (c): Xiaoyu Zhang (xiaoyu@csusm.edu)

  This file is part of sdf (signed distance function).

  sdf is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  sdf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with sdf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "DistanceTransform.h"
#include <math.h>

static const double TOLERANCE = 1e-6;
 
static bool isEqual (double one, double two)
{
	if ( (-1*TOLERANCE <= (one-two)) && ((one-two) <= TOLERANCE) )
		return true;
	return false;
}

static bool isZero(double num)
{
	if ( (-1*TOLERANCE <= num) && (num <= TOLERANCE) )
		return true;
	return false;
}

static bool isNegative(double num)
{
	if (num <0)	return true;
	return false;
}

static bool isBetween(double one, double two, double num)
{
	if ( ((one<=num) && (num<=two)) || ((isEqual(num, one)) || (isEqual(num, two))) )
		return true;
	return false;
}

static bool isZero(const Point3f& one)
{
	double val = sqrt(one[0]*one[0] + one[1]*one[1] + one[2]*one[2]);
	
	if (isZero(val))
		return true;
	return false;
}

static bool isSame(const Point3f& one, const Point3f& two)
{
	if (isZero(one - two))
		return true;
	return false;
}

//static bool edgeIntersectSquare(const Point3f& p0, const Point3f& p1, 
static bool edgeIntersectSquare(double y0, double z0, double y1, double z1,
								double ymin, double ymax, double zmin, double zmax, int flag)
{
//	double y0, y1, z0, z1;
//	switch(flag) {
//	case 1:			// Y-Z plane
//		y0 = p0[1]; y1 = p1[1];
//		z0 = p0[2]; z1 = p1[2];
//		break;
//	case 2:			// X-Z plane
//		y0 = p0[0]; y1 = p1[0];
//		z0 = p0[2]; z1 = p1[2];
//		break;
//	case 3:			// X-Y plane
//		y0 = p0[0]; y1 = p1[0];
//		z0 = p0[1]; z1 = p1[1];
//	default:
//		assert(0);
//	}

	double t;
	if (!isZero(y1 - y0)) {
		t = (ymin - y0) / (y1 - y0);
		double z = z0 + t * (z1 - z0);
		if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, z)))	return true;

		t = (ymax - y0) / (y1 - y0);
		z = z0 + t * (z1 - z0);
		if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, z)))	return true;
	}

	if (!isZero(z1 - z0)) {
		t = (zmin - z0) / (z1 - z0);
		double y = y0 + t*(y1 - y0);
		if ((isBetween(0.0, 1.0, t)) && (isBetween(ymin, ymax, y)))	return true;

		t = (ymax - y0) / (y1 - y0);
		y = y0 + t*(y1 - y0);
		if ((isBetween(0.0, 1.0, t)) && (isBetween(ymin, ymax, y)))	return true;
	}

	return false;
}

static bool edgeIntersectFace(double x0, double y0, double z0, double x1, double y1, double z1,
							  double xmin, double ymin, double ymax, double zmin, double zmax)
{
	double t;
	if (isEqual(x0, xmin) && isEqual(x1, xmin)) { //edge (p0, p1) is in the x = xmin plane
		return edgeIntersectSquare(y0, z0, y1, z1, ymin, ymax, zmin, zmax, 1);
	} else if (!isEqual(x0, x1)) {
		t = (xmin - x0) / (x1 - x0);
		if(isBetween(0, 1, t)) {
			double y = y0 + t*(y1 - y0);
			double z = z0 + t*(z1 - z0);
			if ((isBetween(ymin, ymax, y)) && isBetween(zmin, zmax, z)) {
				return true;
			}
		}
	}	
	return false;
}

static bool edgeIntersectCube(const Point3f& p0, const Point3f& p1, const Point3f& lower,
							  const Point3f& upper)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = lower[0]; xmax = upper[0];
	ymin = lower[1]; ymax = upper[1];
	zmin = lower[2]; zmax = upper[2];
	
	//1) face 1 x = xmin...
	if(edgeIntersectFace(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2],xmin,
		ymin, ymax, zmin, zmax)) {
		return true;
	}

	//1) b> face x = xmax...
	if(edgeIntersectFace(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], xmax,
		ymin, ymax, zmin, zmax)) {
		return true;
	}

	//1) face y = ymin...
	if (edgeIntersectFace(p0[1], p0[0], p0[2], p1[1], p1[0], p1[2], ymin,
		xmin, xmax, zmin, zmax)) {
		return true;
	}

	//1) face y = ymax ...
	if (edgeIntersectFace(p0[1], p0[0], p0[2], p1[1], p1[0], p1[2], ymax,
		xmin, xmax, zmin, zmax)) {
		return true;
	}

	//1) face z = zmin...
	if (edgeIntersectFace(p0[2], p0[0], p0[1], p1[2], p1[0], p1[1], zmin,
		xmin, xmax, ymin, ymax)) {
		return true;
	}

	//1) face z = zmax...
	if (edgeIntersectFace(p0[2], p0[0], p0[1], p1[2], p1[0], p1[1], zmax,
		xmin, xmax, ymin, ymax)) {
		return true;
	}

	return false;
}

bool DistanceTransform::TriangleCubeIntersection(int nt, const Point3f& lower, const Point3f& upper)
{
	Point3f vert0, vert1, vert2;
	p_Surf->getTriVerts(nt, vert0, vert1, vert2);
	double x0, y0, z0, x1, y1, z1, x2, y2, z2;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = lower[0]; xmax = upper[0];
	ymin = lower[1]; ymax = upper[1];
	zmin = lower[2]; zmax = upper[2];

	// find if completely on some side !
	x0 = vert0[0]; 
	y0 = vert0[1]; 
	z0 = vert0[2]; 

	x1 = vert1[0]; 
	y1 = vert1[1];
	z1 = vert1[2];

	x2 = vert2[0]; 
	y2 = vert2[1];
	z2 = vert2[2];

	// If all 3 vertices on same side, then return false;
	if( (x0<xmin) && (x1<xmin) && (x2<xmin) ) return false;
	if( (x0>xmax) && (x1>xmax) && (x2>xmax) ) return false;

	if( (y0<ymin) && (y1<ymin) && (y2<ymin) ) return false;
	if( (y0>ymax) && (y1>ymax) && (y2>ymax) ) return false;

	if( (z0<zmin) && (z1<zmin) && (z2<zmin) ) return false;
	if( (z0>zmax) && (z1>zmax) && (z2>zmax) ) return false;

	//if either of the 3 points are inside the cube, then the Triangle
	// MUST intersect the cube....
	if ( (isBetween(xmin, xmax, x0)) && (isBetween(ymin, ymax, y0)) 
		&& (isBetween(zmin, zmax, z0)) ) return true;
	if ( (isBetween(xmin, xmax, x1)) && (isBetween(ymin, ymax, y1)) 
		&& (isBetween(zmin, zmax, z1)) ) return true;
	if ( (isBetween(xmin, xmax, x2)) && (isBetween(ymin, ymax, y2)) 
		&& (isBetween(zmin, zmax, z2)) ) return true;

	// then the dammed cases when the Triangle intersects any face (edge) 
	// of the cube or vice versa...
	//A) Triangle with the cube...
	// Triangle edge 1
	if (edgeIntersectCube(vert0, vert1, lower, upper)) {
		return true;
	}

	// Triangle edge 2
	if (edgeIntersectCube(vert0, vert2, lower, upper)) {
		return true;
	}

	// Triangle edge 3
	if (edgeIntersectCube(vert1, vert2, lower, upper)) {
		return true;
	}

	///////////////////////////
	//Then the case where the 12 Cube edges intersect the Triangle...
	Point3f p0, p1;

	// x direction
	p0[0] = xmin;
	p1[0] = xmax;

	p0[1] = p1[1] = ymin;
	p0[2] = p1[2] = zmin;
	
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}
	p0[1] = p1[1] = ymax;
	p0[2] = p1[2] = zmin;
	
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}

	p0[1] = p1[1] = ymin;
	p0[2] = p1[2] = zmax;
	
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}
	
	p0[1] = p1[1] = ymax;
	p0[2] = p1[2] = zmax;
	
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}
	
	// y direction
	p0[1] = ymin;
	p1[1] = ymax;
	
	p0[0] = p1[0] = xmin;
	p0[2] = p1[2] = zmin;
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}

	p0[0] = p1[0] = xmax;
	p0[2] = p1[2] = zmin;
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}

	p0[0] = p1[0] = xmin;
	p0[2] = p1[2] = zmax;
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}

	p0[0] = p1[0] = xmax;
	p0[2] = p1[2] = zmax;
	if (rayTriangleIntersection(nt, p0, p1) <= 1) {
		return true;
	}

	// z direction
	p0[2] = zmin;
	p1[2] = zmax;

	p0[0] = p1[0] = xmin;
	p0[1] = p1[1] = ymin;
	if (rayTriangleIntersection(nt, p0, p1)<= 1) {
		return true;
	}

	p0[0] = p1[0] = xmax;
	p0[1] = p1[1] = ymin;
	if (rayTriangleIntersection(nt, p0, p1)<= 1) {
		return true;
	}

	p0[0] = p1[0] = xmin;
	p0[1] = p1[1] = ymax;
	if (rayTriangleIntersection(nt, p0, p1)<= 1) {
		return true;
	}

	p0[0] = p1[0] = xmax;
	p0[1] = p1[1] = ymax;
	if (rayTriangleIntersection(nt, p0, p1)<= 1) {
		return true;
	}

	return false;
}