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
#include "mtxlib.h"

#include <set>

const float DistanceTransform::MAX_FLOAT = (float)1.0e12;

void DistanceTransform::init() {
	int i, j, k;
	float orig[3], span[3];
	int dim[3];
	p_Data->getOrig(orig);
	p_Data->getSpan(span);
	p_Data->getDim(dim);

	p_Cells = new Cell[p_Data->getNCells()];
	// build up the triangle list for cells
	for(int nt = 0; nt < p_Surf->triCount(); nt++) {
		Point3f v0, v1, v2;
		p_Surf->getTriVerts(nt,v0,v1,v2);
		float lower[3], upper[3];
		for(j = 0; j < 3; j++) {
			lower[j] = v0[j];
			upper[j] = v0[j];
		}
		for(j = 0; j < 3; j++) {
			if(v1[j] < lower[j]) lower[j] = v1[j];
			else if(v1[j] > upper[j]) upper[j] = v1[j];
			if(v2[j] < lower[j]) lower[j] = v2[j];
			else if(v2[j] > upper[j]) upper[j] = v2[j];
		}
	
		int minId[3], maxId[3];
		for(j = 0; j < 3; j++) {
			minId[j] = (int)floor((lower[j]-orig[j]) / span[j]);
			maxId[j] = (int)floor((upper[j]-orig[j]) / span[j]);
		}
		for(k = minId[2]; k <= maxId[2]; k++) {
			for(j = minId[1]; j <= maxId[1]; j++) {
				for (i = minId[0]; i <= maxId[0]; i++) {
					int nc = p_Data->index2cell(i, j, k);
					assert(nc < p_Data->getNCells() && nc >= 0);
					Vector3f norm;
					p_Surf->getTriNormal(nt, norm);
					if(intersectCell(v0, norm, i, j, k)) {
						Point3f lower, upper;
						lower[0] = orig[0] + i*span[0];
						lower[1] = orig[1] + j*span[1];
						lower[2] = orig[2] + k*span[2];
						upper[0] = lower[0] + span[0];
						upper[1] = lower[1] + span[1];
						upper[2] = lower[2] + span[2];
						if (TriangleCubeIntersection(nt, lower, upper)) {
							p_Cells[nc].triList.insert(nt);
						}			
					}
				}
			}
		}	
	}

	// Initialize the distance function of points near the surface to 0
	for (k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {
			for(i = 0; i < dim[0]; i++) {
				int nv = p_Data->index2vert(i, j, k);
				if(nearSurface(i, j, k)) {
					p_Data->setValue(nv, 0);
				} else {
					p_Data->setValue(nv, MAX_FLOAT);
				}
			}
		}
	}
	p_Data->init();
	printf("distance min: %f, max = %f\n", p_Data->getFuncMin(), p_Data->getFuncMax());
}

DistanceTransform::DistanceTransform(FaceVertSet3D& fvs, int dim[3], float dist,
									 float sx, float sy, float sz)
{
	float factor[3];
	factor[0] = sx; factor[1] = sy; factor[2] = sz;
	p_Data = new Reg3Data<float> (dim, NULL);
	p_Surf = &fvs;
	
	BoundingBox bbox = p_Surf->getExtent();
	float center[3], ext[3], orig[3], span[3];
	int i, j, k;
	float max_fac = 0, max_ext = 0;
	for(i = 0; i < 3; i++) {
		center[i] = (bbox.lower[i] + bbox.upper[i]) / 2;
		ext[i] = bbox.upper[i] - bbox.lower[i];
		if(ext[i] > max_ext) max_ext = ext[i];
//		float fac = ext[i] / (dim[i] - 1);
//		if(fac > max_fac) max_fac = fac;
//		orig[i] = center[i] - factor[i]*ext[i] / 2;
//		span[i] = factor[i]*ext[i] / (dim[i] - 1);
	}	
	const float DISTANCE = dist*1.1;
	for(i = 0; i < 3; i++) {
		if (max_ext *(factor[i]-1) < 2.2 * DISTANCE) {
			factor[i] = 1 + (2.2*DISTANCE) / max_ext;
		}
	}
	printf("scaling factors: %f %f %f\n", factor[0], factor[1], factor[2]);

	for(i = 0; i < 3; i++) {
		orig[i] = center[i] - factor[i]*max_ext / 2;
		span[i] = factor[i]*max_ext / (dim[i] - 1);
	}
	//float scale = 0.5f / max_fac;
	//float trans[3];
	//for(i = 0; i < 3; i++) {
	//	trans[i] = (dim[i]-1) / 2 - (center[i]-bbox.lower[i])*scale;
	//}
	//for(i = 0; i < p_Surf->vertCount(); i++) {
	//	for(j = 0; j < 3; j++) {
	//		(*(p_Surf->pVerts))[i][j] = ((*(p_Surf->pVerts))[i][j] - bbox.lower[j])*scale + trans[j];
	//	}
	//}
	//orig[0] = orig[1] = orig[2] = 0;
	//span[0] = span[1] = span[2] = 1;

	p_Data->setOrig(orig);
	p_Data->setSpan(span);
#ifdef _DEBUG
	printf("lower: %f %f %f\n", bbox.lower[0], bbox.lower[1], bbox.lower[2]);
	printf("upper: %f %f %f\n", bbox.upper[0], bbox.upper[1], bbox.upper[2]);
	printf("center: %f %f %f\n", center[0], center[1], center[2]);
	printf("orig: %f %f %f\n", orig[0], orig[1], orig[2]);
	printf("span: %f %f %f\n", span[0], span[1], span[2]);
#endif
	init();
}

DistanceTransform::DistanceTransform(FaceVertSet3D &fvs, const Reg3Data<float>& reg3)
{
	p_Data = new Reg3Data<float>(reg3);
	p_Surf = &fvs;

	init();
}

DistanceTransform::~DistanceTransform(void)
{
	if (p_Data)
	{
		delete p_Data;
	}
	if(p_Cells) delete[] p_Cells;
}

void DistanceTransform::transform()
{
	int nverts = p_Data->getNVerts();
	int *parent = new int[nverts];
	int i, j, k;
	for(i = 0; i < nverts; i++) {
		parent[i] = i;
	}

	int dim[3];
	float orig[3], span[3];
	p_Data->getDim(dim);
	p_Data->getSpan(span);
	p_Data->getOrig(orig);

	// 1D distance transform along x

	float *f = new float[dim[0]];
	float *d = new float[dim[0]];
	int *p = new int[dim[0]];
	for (k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {	
			int index = p_Data->index2vert(0, j, k);
			for(i = 0; i < dim[0]; i++) {
				f[i] = p_Data->getValue(index+i);
				p[i] = parent[index+i];
			}
			transform1D(dim[0], f, d, p, span[0]);
			for(i = 0; i < dim[0]; i++) {
				p_Data->setValue(index+i, d[i]);
				parent[index+i] = p[i];
			}
		}
	}

	delete[] f;
	delete[] d;
	delete[] p;

	f = new float[dim[1]];
	d = new float[dim[1]];
	p = new int[dim[1]];
	// 1D distance transform along y
	for(k = 0; k < dim[2]; k++) {
		for(i = 0; i < dim[0]; i++) {
			int index = p_Data->index2vert(i, 0, k);
			for(j = 0; j < dim[1]; j++) {
				f[j] = p_Data->getValue(index+dim[0]*j);
				p[j] = parent[index+dim[0]*j];
			}
			transform1D(dim[1], f, d, p, span[1]);
			for(j = 0; j < dim[1]; j++){
				p_Data->setValue(index+dim[0]*j, d[j]);
				parent[index+dim[0]*j] = p[j];
			}
		}
	}

	delete[] f;
	delete[] d;
	delete[] p;
	
	f = new float[dim[2]];
	d = new float[dim[2]];
	p = new int[dim[2]];
	
	for(j = 0; j < dim[1]; j++) {
		for (i = 0; i < dim[0]; i++) {
			int index = p_Data->index2vert(i, j, 0);
			for(k = 0; k < dim[2]; k++) {
				f[k] = p_Data->getValue(index+k*dim[0]*dim[1]);
				p[k] = parent[index+k*dim[0]*dim[1]];
			}
			transform1D(dim[2], f, d, p, span[2]);
			for(k = 0; k < dim[2]; k++) {
				p_Data->setValue(index+k*dim[0]*dim[1], d[k]);
				parent[index+k*dim[0]*dim[1]] = p[k];
			}
		}
	}
	delete[] f;
	delete[] d;
	delete[] p;

//	for(i = 0; i < nverts; i++) {
//		if (i%1000 == 1) {
//			
//			int ix, iy, iz;
//			p_Data->vert2index(parent[i], ix, iy, iz);
//			int px, py, pz;
//			p_Data->vert2index(i, px, py, pz);
//			float d = (px-ix)*(px-ix)*span[0]*span[0] +
//					  (py-iy)*(py-iy)*span[1]*span[1] +
//					  (pz-iz)*(pz-iz)*span[2]*span[2];
//			printf("vert : %d, val = %f, dis = %f\n", i, p_Data->getValue(i), d);
//			printf("ix: %d, iy: %d, iz: %d\n", ix, iy, iz);
//			assert(nearSurface(ix, iy, iz));
//		}
// 	}

	Point3f *closest = new Point3f[nverts];
	// Compute the distance functions for points near the surface
	for (k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {
			for (i = 0; i < dim[0]; i++) {
				int nv = p_Data->index2vert(i, j, k);
				//if (nearSurface(i, j, k)) {
				if(parent[nv] == nv) {	
					float d = computeNearDistance(i, j, k, closest[nv]);
					p_Data->setValue(nv, d);
				}
			}
		}
	}

	// update the distances for other points
	for (k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {
			for (i = 0; i < dim[0]; i++) {
				int nv = p_Data->index2vert(i, j, k);
				if (parent[nv] != nv) {
					Point3f p(orig[0]+i*span[0], orig[1]+j*span[1], orig[2]+k*span[2]);
					float d = pointDistance(p, closest[parent[nv]]);
					assert(parent[parent[nv]] == parent[nv]);
					if(p_Data->getValue(parent[nv]) < 0) {
						d = -d;
					}
//#ifdef _DEBUG
//					if(p_Data->isBoundary(nv) && d < 0) {
//						int pid = parent[nv];
//						int ix, iy, iz;
//						float distn = p_Data->getValue(pid);
//						printf("WARNING: wrongly classified point %d, close to (%f %f %f), dist=%f\n", pid, 
//								closest[pid][0], closest[pid][1], closest[pid][2], distn);
//						p_Data->vert2index(pid, ix, iy, iz);
//						float d = computeNearDistance(ix, iy, iz, closest[pid]);
//					}
//#endif
					p_Data->setValue(nv, d);
				}
			}
		}
	}
	delete[] closest;

	p_Data->init();
	printf("min dist = %f , max dist = %f\n", p_Data->getFuncMin(), p_Data->getFuncMax());
	delete[] parent;
}

float DistanceTransform::computeNearDistance(int ix, int iy, int iz, Point3f& nearPnt)
{
	std::set<int> triSet;
	int dim[3];
	p_Data->getDim(dim);
	for(int kk = -1; kk <= 0; kk++) {
		for(int jj = -1; jj <= 0; jj++) {
			for(int ii = -1; ii <= 0; ii++) {
				int ci, cj, ck;
				ci = ix + ii;
				cj = iy + jj;
				ck = iz + kk;
				if(ci < 0 || ci >= dim[0]-1 || cj < 0 || cj >= dim[1]-1 || ck < 0 || ck >= dim[2]-1) {
					continue;
				}
				int nc = p_Data->index2cell(ci, cj, ck);
				for(int i = 0; i < p_Cells[nc].triList.length(); i++) {
					triSet.insert(p_Cells[nc].triList[i]);	
				}
			}
		}
	}
	std::set<int>::iterator it = triSet.begin();
//	if (ix % 10 ==1 && iy % 10 == 1 && iz % 10 == 1) {
//		printf("vert: %d %d %d has %d triangles\n", ix, iy, iz, triSet.size());
//		while(it != triSet.end()) {
//			printf("\t tri %d\n", (*it));
//			++it;
//		}
//	}
	int nclose = 0;
	dynamic_array<int> nearTri;
	double dis = MAX_FLOAT;
	double TOLERANCE = 1.0e-6;
	Point3f vert;
	float orig[3], span[3];
	p_Data->getOrig(orig);
	p_Data->getSpan(span);
	vert[0] = orig[0] + ix*span[0];
	vert[1] = orig[1] + iy*span[1];
	vert[2] = orig[2] + iz*span[2];
	while(it != triSet.end()) {
		Point3f myNear;
		int tid = *it;
		assert(tid >= 0);
		double d = distance2Triangle(vert, tid, myNear);
		if(d <= dis - TOLERANCE) {	// new closest intersection
			dis = d;
			nearPnt = myNear;
			nearTri.clear();
			nclose = 0;
			nearTri[nclose++] = tid;
		} else if(fabs(d - dis) < TOLERANCE) {
			if(pointDistance(nearPnt, myNear) < TOLERANCE) { // same intersection for another triangle
				nearTri[nclose++] = tid;
			}
		}
		++it;
	}
	// set the sign of the distance according to (p - p')*n to nearest plane
//	printf("vert: (%f %f %f) has distance %f and nearest point (%f %f %f) on %d triangles\n", 
//			vert[0], vert[1], vert[2], dis, nearPnt[0], nearPnt[1], nearPnt[2], nearTri.length());
	 
	int id;
	if(nearTri.length() == 1) {
		id = 0;
		Vector3f trinorm;
		p_Surf->getTriNormal(nearTri[id], trinorm);
		Vector3f diff = vert - nearPnt;
		if (DotProduct(diff, trinorm) < 0) {
			dis = - dis;
		}
	} else {
		int sgn = inOrOut(nearTri, vert, nearPnt); 
		dis *= sgn;
	}
	
	return dis;
}

//!! this function doesn't work properly
int DistanceTransform::nearestPlane(const dynamic_array<int>& triList, const Point3f& vert)
{
	int i, j, len = triList.length();
	assert(len > 1);

	for (i = 0; i < len; i++) {
		Vector3f tnorm;
		Point3f v0, v1, v2;
		p_Surf->getTriNormal(triList[i], tnorm);
		p_Surf->getTriVerts(triList[i], v0, v1, v2);
		Vector3f diff = vert - v0;
		float dot = DotProduct(tnorm, diff);
		// it may be unstable if the vert is very close to the nearest corner
		const double TOLERANCE = 1.0e-6;
		if(fabs(dot) < TOLERANCE) {
			//printf("WARNING: ");
			continue;
		}
		bool nearest = true;
		for(j = 0; j < len; j++) {
			if (j == i) {
				continue;
			}
			Point3f anv0, anv1, anv2;
			p_Surf->getTriVerts(triList[j], anv0, anv1, anv2);
			float d0 = DotProduct(tnorm, anv0 - v0);
			float d1 = DotProduct(tnorm, anv1 - v0);
			float d2 = DotProduct(tnorm, anv2 - v0);
			float d;
			if(fabs(d0) > fabs(d1)) {
				if(fabs(d0) > fabs(d2)) {
					d = d0;
				} else {
					d = d2;
				}
			} else {
				if (fabs(d1) > fabs(d2)) {
					d = d1;
				} else {
					d = d2;
				}
			}
			if( d * dot > 0) {	// triangle i is not the nearest plane
				nearest = false;
				break;
			}
		}
		if (nearest) {
			return i;
		}

	}
	// should not come here
	printf("WARNING: trouble in finding the nearest plane\n");
/*
	for(i = 0; i < len; i++) {
		TriId3i tid3 = p_Surf->getTriId(triList[i]);
		printf("vert (%f %f %f) 's %dth triangle has vertices %d %d %d\n", vert[0], vert[1],
				vert[2], i, tid3[0], tid3[1], tid3[2]);
	}
*/	

	return 0;
}

int DistanceTransform::inOrOut(const dynamic_array<int> &triList, const Point3f& vert, 
							   const Point3f& nearPnt)
{
	int i, j, len = triList.length();
	assert(len > 1);

	Point3f center, cpnt;
	Point3f v0, v1, v2;
	p_Surf->getTriVerts(triList[0], v0, v1, v2);
	double tnear = 1;
	int nearTri = 0;
	for(i = 0; i < 3; i++) {
		center[i] = (v0[i] + v1[i] + v2[i])/ 3;
		cpnt[i] = 0.9 * nearPnt[i] + 0.1 * center[i];
	}
	for(i = 1; i < len; i++) {
		double t = rayTriangleIntersection(triList[i], vert, cpnt);
		if(t < tnear) {
			tnear = t;
			nearTri = i;
		}
	}

	Vector3f trinorm;
	p_Surf->getTriNormal(triList[nearTri], trinorm);
	p_Surf->getTriVerts(triList[nearTri], v0, v1, v2);
	Vector3f diff = vert - v0;
	if (DotProduct(diff, trinorm) < 0) {
		return -1;
	}
	return 1;
}

	
static double signedDist2Plane(const Point3f& pnt, const Point3f& v0, const Vector3f& norm)
{
	Vector3f vdiff;
	vdiff[0] = pnt[0] - v0[0];
	vdiff[1] = pnt[1] - v0[1];
	vdiff[2] = pnt[2] - v0[2];

	return DotProduct(vdiff, norm);
}

double DistanceTransform::rayTriangleIntersection(int nt, const Point3f& begin, const Point3f& end)
{
	const double TOLERANCE = 1.0e-6;
	Vector3f norm;
	Point3f v0, v1, v2;
	p_Surf->getTriNormal(nt, norm);
	p_Surf->getTriVerts(nt, v0, v1, v2);

	double dbeg = signedDist2Plane(begin, v0, norm);
	double dend = signedDist2Plane(end, v0, norm);
	if(fabs(dbeg) < TOLERANCE) {
		if(pointInTriangle(begin, v0, v1, v2, norm)) return 0;
	}
	if(fabs(dend) < TOLERANCE) {
		if(pointInTriangle(end, v0, v1, v2, norm)) return 1;
	}
	if(dbeg * dend > 0) return MAX_FLOAT;
	if(fabs(dbeg) < TOLERANCE && fabs(dend) < TOLERANCE) return 0.5;  // a hack that needs to fixed later
	double denom = DotProduct(norm, end - begin);
	double t = - dbeg / denom;
	assert (t >= 0); 
	if(t > 1) t = 1;
	Point3f intersect;
	intersect[0] = begin[0] + t*(end[0] - begin[0]);
	intersect[1] = begin[1] + t*(end[1] - begin[1]);
	intersect[2] = begin[2] + t*(end[2] - begin[2]);
	if(pointInTriangle(intersect, v0, v1, v2, norm)) return t;
	return MAX_FLOAT;
}

double DistanceTransform::distance2Triangle(const Point3f& pnt, int nt, Point3f& nearPnt)
{
	double  dist, temp[3];
	Point3f v0, v1, v2;
	p_Surf->getTriVerts(nt, v0, v1, v2);
	Vector3f norm;
	p_Surf->getTriNormal(nt, norm);
	vector3 vnorm(norm[0], norm[1], norm[2]), vdiff;
	
	Point3f edgePnt;

	//1) First compute the shortest distance from the pnt to the plane of the Triangle.
	vdiff[0] = pnt[0] - v0[0];
	vdiff[1] = pnt[1] - v0[1];
	vdiff[2] = pnt[2] - v0[2];
	dist = DotProduct(vdiff, vnorm);
	nearPnt[0] = pnt[0] - dist * norm[0];
	nearPnt[1] = pnt[1] - dist * norm[1];
	nearPnt[2] = pnt[2] - dist * norm[2];

	//2) Then check if the projected point is within the triangle or not. 
	//	  if yes, then the above is the correct shortest distance.
	if(pointInTriangle(nearPnt, v0, v1, v2, norm)) {
		return fabs(dist);
	}

    //3) now, the closest point must on the edge.
	// find the nearest point on each edge.
	dist = distance2Edge(pnt, v0, v1, nearPnt);
	double edgeDist = distance2Edge(pnt, v1, v2, edgePnt);
	if(edgeDist < dist) {
		dist = edgeDist;
		nearPnt = edgePnt;
	}
	edgeDist = distance2Edge(pnt, v2, v0, edgePnt);
	if(edgeDist < dist) {
		dist = edgeDist;
		nearPnt = edgePnt;
	}
	
	assert(dist >= 0);
	return dist;
}

double DistanceTransform::distance2Edge(const Point3f& pnt, const Point3f& v1, const Point3f& v2, Point3f& ne)
{
	vector3 edge(v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]);
	float len = edge.length();
	edge.normalize();

	vector3 vec(pnt[0]-v1[0], pnt[1]-v1[1], pnt[2]-v1[2]);
	float dot = DotProduct(vec, edge);

	if(dot < 0) {
		ne = v1;
		return pointDistance(pnt, v1);
	} else if(dot > len) {
		ne = v2;
		return pointDistance(pnt, v2);
	}
	ne[0] = v1[0] + dot * edge[0]; 
	ne[1] = v1[1] + dot * edge[1];
	ne[2] = v1[2] + dot * edge[2];
	return pointDistance(pnt, ne);
}

bool DistanceTransform::pointInTriangle(const Point3f& res, const Point3f& vert0, 
										const Point3f& vert1, const Point3f& vert2, const Vector3f& norm)
{
	double alpha, beta;
	double u0, u1, u2, v0, v1, v2;
	int index;

	if(fabs(norm[0]) > fabs(norm[1])) {
		if(fabs(norm[0]) > fabs(norm[2])) {
			index = 0;
		} else {
			index = 2;
		}
	} else {
		if(fabs(norm[1]) > fabs(norm[2])) {
			index = 1;
		} else {
			index = 2;
		}
	}
	
	if( index == 0 ) {		// project to y-z plane
		u0 = res[1] - vert0[1]; 
		u1 = vert1[1] - vert0[1]; 
		u2 = vert2[1] - vert0[1]; 

		v0 = res[2] - vert0[2]; 
		v1 = vert1[2] - vert0[2]; 
		v2 = vert2[2] - vert0[2]; 
	}
	else if ( index == 1 ) {	// project to x-z plane
		u0 = res[2] - vert0[2]; 
		u1 = vert1[2] - vert0[2]; 
		u2 = vert2[2] - vert0[2];

		v0 = res[0] - vert0[0];
		v1 = vert1[0] - vert0[0];
		v2 = vert2[0] - vert0[0];
	} else 	{	// project to x-y plane
		u0 = res[0] - vert0[0];
		u1 = vert1[0] - vert0[0];
		u2 = vert2[0] - vert0[0];

		v0 = res[1] - vert0[1]; 
		v1 = vert1[1] - vert0[1]; 
		v2 = vert2[1] - vert0[1]; 
	}

	alpha = ( u0*v2 - v0*u2 ) / ( u1*v2 - v1*u2 );
	if( alpha < 0 ) return false;

	beta = ( u1*v0 - v1*u0 ) / ( u1*v2 - v1*u2 );
	if( beta < 0 ) return false;

	if( alpha+beta > 1) return false;

	return true;
}

bool DistanceTransform::intersectCell(const Point3f& v0, const Vector3f& norm, int i, int j, int k)
{
	vector3 n(norm[0], norm[1], norm[2]);

	float xyz[3], orig[3], span[3];
	p_Data->getOrig(orig);
	p_Data->getSpan(span);
	xyz[0] = orig[0] + i*span[0];
	xyz[1] = orig[1] + j*span[1];
	xyz[2] = orig[2] + k*span[2];
	vector3 p(xyz[0]-v0[0], xyz[1]-v0[1], xyz[2]-v0[2]);
	float dot = DotProduct(n, p);

	for(int kk = 0; kk <= 1; kk++) {
		for(int jj = 0; jj <= 1; jj++) {
			for(int ii = 0; ii <= 1; ii++) {
				xyz[0] = orig[0] + (i+ii)*span[0];
				xyz[1] = orig[1] + (j+jj)*span[1];
				xyz[2] = orig[2] + (k+kk)*span[2];
				vector3 p(xyz[0]-v0[0], xyz[1]-v0[1], xyz[2]-v0[2]);
				float dot2 = DotProduct(n, p);
				if(dot*dot2 <= 0) {		// two vertices on the opposite sides of the triangle plane
					return true;
				}
			}
		}
	}
	return false;
}

bool DistanceTransform::nearSurface(int i, int j, int k)
{
	int dim[3];
	p_Data->getDim(dim);
	for(int kk = -1; kk <= 0; kk++) {
		for(int jj = -1; jj <= 0; jj++) {
			for(int ii = -1; ii <= 0; ii++) {
				int ci, cj, ck;
				ci = i + ii;
				cj = j + jj;
				ck = k + kk;
				if(ci < 0 || ci >= dim[0]-1 || cj < 0 || cj >= dim[1]-1 || ck < 0 || ck >= dim[2]-1) {
					continue;
				}
				int nc = p_Data->index2cell(ci, cj, ck);
				if(p_Cells[nc].triList.length() > 0) return true;
			}
		}
	}
	return false;
}

void DistanceTransform::transform1D(int n, float f[], float d[], int parent[], float span)
{
	int i, k = 0;		// Index of rightmost parabola in lower envelope
	int *v = new int[n];
	int *temp = new int[n];
	float *z = new float[n+1];
	float s2 = span*span;
 
	float min = f[0];
	d[0] = f[0];
	for (i = 1; i < n; i++) {
		if(f[i] < min) min = f[i];
		d[i] = f[i];
	}
	if(min >= MAX_FLOAT) return;		// every number is max, no need to transform.

	for(i = 0; i < n; i++) {
		if(f[i] < MAX_FLOAT) f[i] = f[i] / s2 ;
		temp[i] = parent[i];
	}

	v[0] = 0;
	z[0] = -MAX_FLOAT;
	z[1] = MAX_FLOAT;

	for(i = 1; i < n; i++) {
		//intersection of ith parabola and the current right most parabola 
		float s; 
		if(f[i] < MAX_FLOAT) {
			s = ((f[i] + i*i) - (f[v[k]] + v[k]*v[k])) / (2*i-2*v[k]);
		} else {
			s = MAX_FLOAT;
		}
		
		while (s <= z[k]) {
			k --;		// delete one parabola from the lower envelop
			s = ((f[i] + i*i) - (f[v[k]] + v[k]*v[k])) / (2*i-2*v[k]);
		}
		k ++;
		v[k] = i;
		z[k] = s;
		z[k+1] = MAX_FLOAT;
	}

	k = 0;
	for(i = 0; i < n; i++) {
		while (z[k+1] < i) {
			k++;
		}
		if(f[v[k]] < MAX_FLOAT) {
			d[i] = (i - v[k])*(i-v[k]) + f[v[k]];
			d[i] *= s2;
			parent[i] = temp[v[k]];
		} else {
			d[i] = MAX_FLOAT;
		}
		
	}

	delete[] temp;
	delete[] v;
	delete[] z;
}
