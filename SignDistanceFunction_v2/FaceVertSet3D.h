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

#ifndef FACEVERTSET3D_H
#define FACEVERTSET3D_H

#include <set>

#include "geom.h"
#include "dynarray.h"

/**
 * A 3D triangular mesh represented as face vertex set. 
 */
class FaceVertSet3D {
	friend class Geom3DParser;
	friend class DistanceTransform;
public:    
	/**
	 * Constructor 
	 * @note 
	 */
	FaceVertSet3D(int nv = 0, int nt = 0, Point3f *verts = 0, TriId3i *tids = 0, Vector3f *norms = 0);

	/**
	 * Destructor 
	 */
	virtual ~FaceVertSet3D();

	/**
	 * Get the number of vertices	
	 */
	int vertCount() const {
		return nvert;
	}

	/**
	 *	Get the number of triangles
	 */
	int triCount() const {
		return ntri;
	}

	/** 
	 * Add a vertex with the given position and normal
	 * @return The id of the newly added vertex
	 */ 
	int AddVert(const Point3f& pos, const Vector3f& norm);

	int addVert(float x, float y, float z);

	/** 
	 *	Add a vertex to the mesh if it doesn't exist now. Otherwise return the id of the vertex
	 * @return The id of the vertex
	 */
	int AddVertUnique(const Point3f& pos, const Vector3f& norm);

	/** 
	 * Add a triangle with given vert indices
	 * @return The id of the new triangle
	 */
	int AddTri(const TriId3i& id);

	/**
	 * Remove unessential member fields to reduce memory usage.	
	 */
	void compact();

	/**
	 * Set up bounding box of the mesh
	 */ 
	void buildBBox();

	/**
	 *	Get the bounding box of the mesh
	 */
	virtual BoundingBox getExtent() const {
		return bbox;
	}

	void getTriVerts(int nt, Point3f& v0, Point3f& v1, Point3f& v2) const {
		v0 = (*pVerts)[(*pTris)[nt][0]];
		v1 = (*pVerts)[(*pTris)[nt][1]];
		v2 = (*pVerts)[(*pTris)[nt][2]];
	}

	void getTriNormal(int nt, Vector3f& norm) {
		assert(nt < ntri && nt >= 0);
		norm = (*pTriNorms)[nt];
	}

	TriId3i getTriId(int nt) {
		return (*pTris)[nt];
	}
	
	/**
	 *	Set the bounding box of the mesh
	 */
	virtual void setExtent(const BoundingBox& _box) {
		bbox = _box;
	}

	// from Node3D
	virtual void render();

	/**
	 * Compute vertex normals of the mesh
	 * @param force: If force is true, force to recompute normals
	 */ 
	virtual void ComputeNormal(bool force = false);

	void computeTriNormals();

	void flipTriNormals() {
		for (int i = 0; i < ntri; i++) {
			(*pTriNorms)[i] = -(*pTriNorms)[i];
		}
	}

	void invertNormal() {
		for(int i = 0; i < nvert; i++) {
			(*pNorms)[i] = -((*pNorms)[i]);
		}
	}
private:
	/**
	 * The number of vertices and triangles	
	 */ 
	int nvert, ntri;

	// Bounding box of the mesh
	BoundingBox     bbox;

	// Vertex array, normal arry, and triangle index array
	dynamic_array<Point3f>      *pVerts;			// vertex array
	dynamic_array<Vector3f>     *pNorms;			// normal array 
	dynamic_array<TriId3i>      *pTris;				// triangle index array
	dynamic_array<Vector3f>		*pTriNorms;			// triangle normal
	std::set<VertIdPair, LTVert>	*pVSet;
};
#endif //FACEVERTSET3D_H
