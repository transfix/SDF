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
#pragma once

#include "FaceVertSet3D.h"
#include "reg3data.h"
#include "RawivParser.h"

class Cell {
public:
	Cell() {
		
	}
	~Cell() {
		
	}
	
	bool hasTriangle() const {
		return triList.length() > 0;
	}

	dynamic_array<int> triList;
};

class DistanceTransform
{
	Reg3Data<float> *p_Data;
	Cell*	p_Cells;
	FaceVertSet3D* p_Surf;
public:
	// sx, sy, sz: scale factor of the bounding volume
	DistanceTransform(FaceVertSet3D& fvs, int dim[3], float Distance = 20, float sx=2, float sy=2, float sz=2);

	DistanceTransform(FaceVertSet3D& fvs, const Reg3Data<float>& reg3);

	~DistanceTransform(void);

	void transform();

	void writeRawiv(const char* fname) {
		RawivParser parser;
		parser.write(*p_Data, fname);
	}

	//////////////////////////////////////////////////////////////////////////
	void transform1D(int len, float f[], float d[], int parent[], float span);
//private:
	const static float MAX_FLOAT;
	///
	bool intersectCell(const Point3f& v0, const Vector3f& norm, int i, int j, int k);

	//////////////////////////////////////////////////////////////////////////
	bool nearSurface(int i, int j, int k);

	//////////////////////////////////////////////////////////////////////////
	float computeNearDistance(int ix, int iy, int iz, Point3f& pnt);

	//////////////////////////////////////////////////////////////////////////
	double distance2Triangle(const Point3f& pnt, int nt, Point3f& ne);

	//////////////////////////////////////////////////////////////////////////
	double distance2Edge(const Point3f& pnt, const Point3f& v1, const Point3f& v2, Point3f& ne);

	//////////////////////////////////////////////////////////////////////////
	bool pointInTriangle(const Point3f& res, const Point3f& v0, const Point3f& v1, 
						 const Point3f& v2, const Vector3f& norm);

	//////////////////////////////////////////////////////////////////////////
	int nearestPlane(const dynamic_array<int>& triList, const Point3f& vert);

	///
	double rayTriangleIntersection(int nt, const Point3f& begin, const Point3f& end);

	/**
	 * @return -1, if point is inside, return 1 if point is outside
	 */
	int inOrOut(const dynamic_array<int>& triList, const Point3f& vert, const Point3f& nearPnt);

	/************************************************************************/
	/* @return 1 if they intersect                                          */
	/************************************************************************/
	bool TriangleCubeIntersection(int nt, const Point3f& lower, const Point3f& upper);

	const Reg3Data<float> & getReg3Data() {
		return *p_Data;
	}
private:
	void init();
};
