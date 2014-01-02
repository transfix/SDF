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
/******
 *
 * geom.h : Commonly used geometric primitives
 *
 * (c) 2000 Xiaoyu Zhang
 *
 * v. 1.0
 ******/

#ifndef ZVIS_GEOMETRY_H
#define	ZVIS_GEOMETRY_H

#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define    EPS  0.00001

struct Point3f {
	float pnt[3];
	
	Point3f(float x = 0, float y = 0, float z = 0) {
		pnt[0] = x;
		pnt[1] = y;
		pnt[2] = z;
	}

	Point3f(const float* coord) {
		for(int i = 0; i < 3; i++) {
			pnt[i] = coord[i];
		}
	}
	
	const float operator []  (int i) const {
		return pnt[i];
	}

	float& operator [] (int i) {
		return pnt[i];
	}

	Point3f& operator = (const Point3f& p) {
		if(this == &p) return *this;
		for(int i = 0; i < 3; i++) {
			pnt[i] = p.pnt[i];
		}
		return *this;
	}

	friend const Point3f operator - (const Point3f& p) {
		return Point3f(-p.pnt[0], -p.pnt[1], -p.pnt[2]);
	}

	friend bool operator == (const Point3f& p1, const Point3f& p2) {
		return (p1.pnt[0] == p2.pnt[0] && p1.pnt[1] == p2.pnt[1] && p1.pnt[2] == p2.pnt[2]);
	}

	friend const Point3f operator - (const Point3f p1, const Point3f p2) {
		return Point3f(p1.pnt[0]-p2.pnt[0], p1.pnt[1]-p2.pnt[1], p1.pnt[2]-p2.pnt[2]);
	}

	friend float pointDistance(const Point3f& p1, const Point3f& p2) {
		return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + 
					(p1[1]-p2[1])*(p1[1]-p2[1]) + 
					(p1[2]-p2[2])*(p1[2]-p2[2]));
	}
};

struct Point4f {
	Point4f(float x = 0, float y = 0, float z = 0, float w = 0) {
		pnt[0] = x;
		pnt[1] = y;
		pnt[2] = z;
		pnt[3] = w;
	}

	Point4f(const float* coord) {
		for(int i = 0; i < 4; i++) {
			pnt[i] = coord[i];
		}
	}

	float pnt[4];

	const float operator []  (int i) const {
		return pnt[i];
	}

	float& operator [] (int i) {
		return pnt[i];
	}

	Point4f& operator = (const Point4f& p) {
		if(this == &p) return *this;
		for(int i = 0; i < 4; i++) {
			pnt[i] = p.pnt[i];
		}
		return *this;
	}

	bool operator == (const Point4f& p) const {
		return (pnt[0] == p.pnt[0] && pnt[1] == p.pnt[1] && pnt[2] == p.pnt[2] && pnt[3] == p.pnt[3]);
	}
};

struct TriId3i {
	unsigned int id[3];

	TriId3i(int v1, int v2, int v3) {
		id[0] = v1;
		id[1] = v2;
		id[2] = v3;
	}

	const unsigned int operator []  (int i) const {
		return id[i];
	}

	unsigned int& operator [] (int i) {
		return id[i];
	}
};

typedef Point3f   Vector3f;
typedef Point4f   Vector4f;

inline void NormalizeVector3f(Vector3f& v)
{
    float  len;
    len = (float)sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (len == 0) {
        return;
    }
    v[0] = v[0]/len;
    v[1] = v[1]/len;
    v[2] = v[2]/len;
}

static float DotProduct(const Vector3f& v1, const Vector3f& v2) {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


typedef struct _Ray {
	Point3f  orig;				/* starting point of ray */
	Vector3f dir;				/* ray direction */
} Ray;

struct BoundingBox {
	Point3f lower;
	Point3f upper;
	
	BoundingBox(float lx = -1, float ly = -1, float lz = -1,
				float ux = 1, float uy = 1, float uz = 1) 
				: lower(lx, ly, lz), upper(ux, uy, uz) {
	}

	BoundingBox& operator = (const BoundingBox& box) {
		if(this == &box) return *this;
		lower = box.lower;
		upper = box.upper;
		return *this;
	}
};

struct BoundingSphere {
	Point3f center;
	float	radius;

	BoundingSphere() : center(0, 0, 0), radius(1) {
	}

	/// assignment operator
	BoundingSphere& operator = (const BoundingSphere& sphere) {
		if(this == &sphere) return *this;
		center = sphere.center;
		radius = sphere.radius;
		return *this;
	}
};



struct VertIdPair {
  Point3f cord;
  int idx;
};

struct LTVert{
  bool operator()(const VertIdPair& p1, const VertIdPair& p2) const
    {
       return ((p1.cord[2] < p2.cord[2]) ||
	       ((p1.cord[2] == p2.cord[2]) && (p1.cord[1] < p2.cord[1])) ||
	       ((p1.cord[2] == p2.cord[2]) && (p1.cord[1] == p2.cord[1]) && (p1.cord[0] < p2.cord[0])));
    }
  
  bool operator()(const Point3f& p1, const Point3f& p2) const
    {
       return ((p1[2] < p2[2]) ||
	       ((p1[2] == p2[2]) && (p1[1] < p2[1])) ||
	       ((p1[2] == p2[2]) && (p1[1] == p2[1]) && (p1[0] < p2[0])));
    }
};

#endif  	// ZVIS_GEOMETRY_H
