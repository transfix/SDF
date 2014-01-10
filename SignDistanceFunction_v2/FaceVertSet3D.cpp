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
#ifdef WIN32
#include <windows.h>
#endif

#if 0
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "mtxlib.h"
#include "FaceVertSet3D.h"

FaceVertSet3D::FaceVertSet3D(int nv, int nt, Point3f *verts, TriId3i *tids, Vector3f *norms)
{
	nvert = nv;
	ntri = nt;
	pVerts = new dynamic_array<Point3f>(nvert, verts);
	pTris = new dynamic_array<TriId3i>(ntri, tids);
	pTriNorms = new dynamic_array<Vector3f>(ntri);

	if(norms) {
        pNorms = new dynamic_array<Vector3f>(nv, norms);
    } else {
        pNorms = new dynamic_array<Vector3f>;
    }
	pVSet = NULL;
	if (nv > 2 ) {
		buildBBox();
	}
}

FaceVertSet3D::~FaceVertSet3D()
{
	delete pVerts;
	delete pTris;
	delete pNorms;
	delete pVSet;
	delete pTriNorms;
}

int FaceVertSet3D::AddVert(const Point3f& pos, const Vector3f& norm)
{
	int n = nvert++;
	pVerts->insert(pos);
	pNorms->insert(norm);

	return n;
}

int FaceVertSet3D::addVert(float x, float y, float z)
{
	int n = nvert++;
	Point3f p(x, y, z);
	Vector3f norm;
	pVerts->insert(p);
	pNorms->insert(norm);
	
	return n;
}

int FaceVertSet3D::AddVertUnique(const Point3f& pos, const Vector3f& norm)
{
	if (!pVSet) {
		pVSet = new std::set<VertIdPair, LTVert>();
	}
	VertIdPair vip;
	vip.idx = nvert;
	vip.cord = pos;
	std::set<VertIdPair, LTVert>::iterator it = pVSet->find(vip);
	if(it != pVSet->end()) {
		return it->idx;
	}
	pVSet->insert(vip);
	return AddVert(pos, norm);
}

int FaceVertSet3D::AddTri(const TriId3i& id)
{
	int n = ntri ++;
	pTris->insert(id);

	return n;
}

void FaceVertSet3D::compact()
{
	//pVSet->erase(pVSet->begin(), pVSet->end());
	delete pVSet;
	pVSet = 0;
}

void FaceVertSet3D::buildBBox()
{
	int i, j;
	Point3f min, max;
	for(i = 0; i < 3; i++) {
		min[i] = (*pVerts)[0][i];
		max[i] = (*pVerts)[0][i];
	}

	for(i = 1; i < nvert; i++) {
		for(j = 0; j < 3; j++) {
			if((*pVerts)[i][j] < min[j]) {
				min[j] = (*pVerts)[i][j];
				//printf("i = %d, j = %d, c = %f\n", i, j, (*pVerts)[i][j]);
			}
			else if((*pVerts)[i][j] > max[j]) 
				max[j] = (*pVerts)[i][j];
		}
	}
	bbox.lower = min;
	bbox.upper = max;
}

#if 0
void FaceVertSet3D::render()
{
	glPushMatrix();

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	
	glVertexPointer(3, GL_FLOAT, 0, pVerts->data());
	glNormalPointer(GL_FLOAT, 0, pNorms->data());
	glDrawElements(GL_TRIANGLES, 3*triCount(), GL_UNSIGNED_INT, pTris->data());
	
//	glBegin(GL_TRIANGLES);
//	for (int i = 0 ; i < ntri; i++ ) {
//		for (int j= 0 ; j < 3; j++) {
//			int vid = (*pTris)[i][j];
//			glNormal3f((*pNorms)[vid][0], (*pNorms)[vid][1], (*pNorms)[vid][2]);
//          glVertex3f((*pVerts)[vid][0], (*pVerts)[vid][1], (*pVerts)[vid][2]);
//		}
//	}
// 	glEnd();
	glPopMatrix();
}
#endif

void FaceVertSet3D::ComputeNormal(bool force) 
{
	// If normals exist and not forced to recompute, just return
	if (!force && pNorms) {
		return;
	}

	// compute vertex normals by averaging triangle normals
	if (pNorms) {
		delete pNorms;
		pNorms = 0;
	}
	Vector3f *vnorms = new Vector3f[nvert];		// vertex normals
	Vector3f *fnorms = new Vector3f[ntri];		// triangle normals
	int i;
	//compute face normals
	for (i = 0; i < ntri; i++) {
		vector3 v1, v2, v3;
		v1.x = (*pVerts)[(*pTris)[i][1]][0] - (*pVerts)[(*pTris)[i][0]][0];
		v2.x = (*pVerts)[(*pTris)[i][2]][0] - (*pVerts)[(*pTris)[i][0]][0];
		v1.y = (*pVerts)[(*pTris)[i][1]][1] - (*pVerts)[(*pTris)[i][0]][1];
		v2.y = (*pVerts)[(*pTris)[i][2]][1] - (*pVerts)[(*pTris)[i][0]][1];
		v1.z = (*pVerts)[(*pTris)[i][1]][2] - (*pVerts)[(*pTris)[i][0]][2];
		v2.z = (*pVerts)[(*pTris)[i][2]][2] - (*pVerts)[(*pTris)[i][0]][2];
		v3 = CrossProduct(v1, v2);
		v3.normalize();
		fnorms[i][0] = v3.x;
		fnorms[i][1] = v3.y;
		fnorms[i][2] = v3.z;
	}
	//sum for vertex nomals
	for (i = 0; i < ntri; i++) {
		for (int j = 0; j < 3; j++) {
			vnorms[(*pTris)[i][0]][j] += fnorms[i][j];
			vnorms[(*pTris)[i][1]][j] += fnorms[i][j];
			vnorms[(*pTris)[i][2]][j] += fnorms[i][j];
		}
	}
	delete[] fnorms;
	for (i = 0; i < nvert; i++) {
		NormalizeVector3f(vnorms[i]);
	}
	pNorms = new dynamic_array<Vector3f>(nvert, vnorms);
	delete[] vnorms;
}

void FaceVertSet3D::computeTriNormals()
{
	int i;
	//compute face normals
	for (i = 0; i < ntri; i++) {
		vector3 v1, v2, v3;
		v1.x = (*pVerts)[(*pTris)[i][1]][0] - (*pVerts)[(*pTris)[i][0]][0];
		v2.x = (*pVerts)[(*pTris)[i][2]][0] - (*pVerts)[(*pTris)[i][0]][0];
		v1.y = (*pVerts)[(*pTris)[i][1]][1] - (*pVerts)[(*pTris)[i][0]][1];
		v2.y = (*pVerts)[(*pTris)[i][2]][1] - (*pVerts)[(*pTris)[i][0]][1];
		v1.z = (*pVerts)[(*pTris)[i][1]][2] - (*pVerts)[(*pTris)[i][0]][2];
		v2.z = (*pVerts)[(*pTris)[i][2]][2] - (*pVerts)[(*pTris)[i][0]][2];
		v3 = CrossProduct(v1, v2);
		v3.normalize();
		Vector3f vec(v3.x, v3.y, v3.z);
		(*pTriNorms)[i] = vec;
	}
}
