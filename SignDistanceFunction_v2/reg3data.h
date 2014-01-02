#ifndef _REG3_DATA_H
#define _REG3_DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "cubes.h"

#define error(x) {fprintf(stderr, "%s\n", x); exit(1);}

/**
 * Reg3Data: A class representing a scalar field on a regular 3D grid.
 *
 * @note This class is based on many previous versions of regular 3D data class.
 */

template <class T>
class Reg3Data {
public:
	friend class Reg3Parser;
	friend class RawivParser;
	friend class RawVParser;
	friend class SimplexGraph;
	friend class BSplineInterpolation;

	// Construct from a data array
	Reg3Data(int dim[3], T* array) {
		for(int i = 0; i < 3; i++) {
			m_dim[i] = dim[i];
			m_orig[i] = 0;
			m_span[i] = 1;
		}
		p_data = new T[dim[0]*dim[1]*dim[2]];
		p_grad = NULL;
		if(array) memcpy(p_data, array, sizeof(T)*dim[0]*dim[1]*dim[2]);
		init();
	}

	Reg3Data(const Reg3Data<T>& data) {
		int i;
		for(i = 0; i < 3; i++) {
			m_dim[i] = data.m_dim[i];
			m_orig[i] = data.m_orig[i];
			m_span[i] = data.m_span[i];
		}
		m_fmin = data.m_fmin;
		m_fmax = data.m_fmax;
		p_data = NULL;
		p_grad = NULL;

		if(data.p_data) {
			p_data = new T[getNVerts()];
			memcpy(p_data, data.p_data, sizeof(T)*getNVerts());
		}
		m_dirty = data.m_dirty;
	}

	// Default Constructor
	Reg3Data() {
		m_dim[0] = m_dim[1] = m_dim[2] = 0;
		m_orig[0] = m_orig[1] = m_orig[2] = 0;
		m_span[0] = m_span[1] = m_span[2] = 1;
		p_data = 0;
		p_grad = 0;
		m_dirty = true;
	}

	// Destructor
	virtual ~Reg3Data() {
		if (p_data) {
			delete[] p_data;
		}
		if (p_grad) {
			free(p_grad);
		}
	}

	void getDim(int dim[3]) const {
		for (int i = 0; i < 3; i++)	dim[i] = m_dim[i];
	}

	int getNVerts() const {
		return m_dim[0]*m_dim[1]*m_dim[2];
	}

	int getNCells() const {
		return(m_dim[0] - 1)*(m_dim[1] - 1)*(m_dim[2] - 1);
	}

	void getOrig(float orig[3]) const {
		for (int i = 0; i < 3; i++)	orig[i] = m_orig[i];
	}

	void setOrig(const float orig[3]) {
		for (int i = 0; i < 3; i++)	m_orig[i] = orig[i];
	}

	void getSpan(float span[3]) const {
		for (int i = 0; i < 3; i++)	span[i] = m_span[i];
	}

	void setSpan(const float span[3]) {
		for (int i = 0; i < 3; i++)	m_span[i] = span[i];
	}

	float getFuncMin() {
		if(m_dirty) init();
		return m_fmin;
	}

	float getFuncMax() {
		if(m_dirty) init();
		return m_fmax;
	}

	void getFuncMinMax(float& min, float& max) const {
		min = m_fmin; max = m_fmax;
	}

	void getBoundingBox(float min[3], float max[3]) const {
		for (int i = 0; i < 3; i++) {
			min[i] = m_orig[i];
			max[i] = m_orig[i] + (m_dim[i]-1)*m_span[i];
		}
	}

	void getVertPos(int vid, float& x, float& y, float& z)
	{
		int xx,yy,zz;
		vert2index(vid, xx, yy, zz);

		x = m_orig[0]+xx*m_span[0];
		y = m_orig[1]+yy*m_span[1];
		z = m_orig[2]+zz*m_span[2];
	}

	float getValue(int nv) const {
		assert(nv < m_dim[0]*m_dim[1]*m_dim[2]);
		return p_data[nv];
	}

	void setValue(int nv, T x) {
		assert(nv < m_dim[0]*m_dim[1]*m_dim[2]);
		p_data[nv] = x;
		m_dirty = true;
	}

	float operator [] (int nv) const {
		return getValue(nv);
	}

	float getValue(int i, int j, int k) const {
		return getValue(index2vert(i, j, k));
	}

	void getVertGrad(int i, int j, int k, float grad[3]) const; 

	void getVertGrad(int v, float grad[3]) const;

	void setVertGrad(int v, float grad[3]) {
		assert(hasGradient());
		p_grad[v][0] = grad[0];
		p_grad[v][1] = grad[1];
		p_grad[v][2] = grad[2];
	}

	void getCellValues(int i, int j, int k, float vals[8]) const;

	void getCellValues(int c, float vals[8]) {
		int i, j, k;
		cell2index(c, i, j, k);
		getCellValues(i, j, k, vals);
	}

	void getCellGrads(int i, int j, int k, float grads[8][3]);

	void getCellGrads(int c, float grads[8][3]) {
		int i, j, k;
		cell2index(c, i, j, k);
		getCellGrads(i, j, k, grads);
	}

	/// Euclidean distance between vertex u and v
	float pointDistance(int u, int v) {
		int idxu[3], idxv[3];
		vert2index(u, idxu[0], idxu[1], idxu[2]);
		vert2index(v, idxv[0], idxv[1], idxv[2]);

		float d2 = 0;
		for (int i = 0; i < 3; i++) {
			float delta = (idxu[i] - idxv[i])*m_span[i];
			d2 += delta*delta;
		}
		return (float)sqrt(d2);
	}

	bool isValidVtx(int vtx_id, float thredshold = 0) {
		if (p_data[vtx_id] >= thredshold) return true;
		return false;
	}

	bool isBoundary(int v) 
	{
		int x,y,z;
		vert2index(v, x, y, z);
		if ( (x<=0 || x>= m_dim[0]-1) ||
			 (y<=0 || y>= m_dim[1]-1) ||
			 (z<=0 || z>= m_dim[2]-1) )	return true;
		return false;
	}

	int getNeighborVtx(int vtx , int* neighbor_vtx);


	/**
	 * Compute the gradient vector at all grid points.
	 * case 0: simple finite difference
	 * case 1: cubic B-spline approximation
	 */
	void calcGradient(int mode = -1);

	int index2vert(int i, int j, int k) const {
		if ( ((i < 0) || ( i >= m_dim[0])) ||
			 ((j < 0) || ( j >= m_dim[1])) ||
			 ((k < 0) || ( k >= m_dim[2])) ) return -1;
		return(k*m_dim[1]*m_dim[0]+j*m_dim[0]+i);
	}

	void vert2index(int v, int &i, int &j, int &k) const {
		i = v % m_dim[0];
		j = (v / m_dim[0]) % m_dim[1];
		k = v / (m_dim[0]*m_dim[1]);
	}

	int index2cell(int i, int j, int k) const {
		 if ( ((i < 0) || ( i >= m_dim[0]-1)) ||
             ((j < 0) || ( j >= m_dim[1]-1)) ||
             ((k < 0) || ( k >= m_dim[2]-1)) ) return -1;
		return( i + j * (m_dim[0]-1) + k *(m_dim[0]-1)*(m_dim[1]-1));       
	}

	void cell2index(int c, int &i, int &j, int &k) const {
		i = c % (m_dim[0]-1);
		j = (c / (m_dim[0]-1)) % (m_dim[1]-1);
		k = c / ((m_dim[0]-1) * (m_dim[1]-1));
	}
	void init();

protected:
	// min max values of the scalar function
	float       m_fmin, m_fmax; 
	int         m_dim[3];
	float       m_orig[3], m_span[3];
	T*          p_data;
	float       (*p_grad)[3];
private:
	bool		m_dirty;
	bool hasGradient() const {
		return(p_grad != NULL);
	}

};

template <class T>
void Reg3Data<T>::init()
{
	m_fmin = m_fmax = p_data[0];
	int i, n = m_dim[0]*m_dim[1]*m_dim[2];

	for (i = 1; i < n; i++) {
		if (p_data[i] < m_fmin)	m_fmin = p_data[i];
		else if (p_data[i] > m_fmax) m_fmax = p_data[i];
	}
	m_dirty = false;
}

template <class T>
void Reg3Data<T>::getVertGrad(int i, int j, int k, float grad[3]) const
{
	int v = index2vert(i, j, k); 
	return getVertGrad(v, grad);
}

template <class T>
void Reg3Data<T>::getVertGrad(int v, float grad[3]) const
{
	assert(hasGradient());
	// read from saved gradients
	grad[0] = p_grad[v][0];
	grad[1] = p_grad[v][1];
	grad[2] = p_grad[v][2];
} 

template<class T>
int Reg3Data<T>::getNeighborVtx(int vtx , int* neighbor_vtx)
{
	int i , j , k , idx;
	int nbr_idx = 0;
	int nidx=0;
	int nbr_vtx[100];

	vert2index(vtx , i , j ,k);
	nbr_vtx[nbr_idx++] = index2vert(i-1 , j   , k  );
	nbr_vtx[nbr_idx++] = index2vert(i+1 , j   , k  );
	nbr_vtx[nbr_idx++] = index2vert(i   , j-1 , k  );
	nbr_vtx[nbr_idx++] = index2vert(i   , j+1 , k  );
	nbr_vtx[nbr_idx++] = index2vert(i   , j   , k-1);
	nbr_vtx[nbr_idx++] = index2vert(i   , j   , k+1);

	nbr_vtx[nbr_idx++] = index2vert(i+1 , j+1 , k  );
	nbr_vtx[nbr_idx++] = index2vert(i-1 , j-1 , k  );
	nbr_vtx[nbr_idx++] = index2vert(i   , j+1 , k-1);
	nbr_vtx[nbr_idx++] = index2vert(i   , j-1 , k+1);
	nbr_vtx[nbr_idx++] = index2vert(i+1 , j   , k-1);
	nbr_vtx[nbr_idx++] = index2vert(i-1 , j   , k+1);
	nbr_vtx[nbr_idx++] = index2vert(i+1 , j+1 , k-1);
	nbr_vtx[nbr_idx++] = index2vert(i-1 , j-1 , k+1);

	for (idx = 0 ; idx < nbr_idx ; idx++) {
		if (nbr_vtx[idx] != -1)	neighbor_vtx[nidx++] = nbr_vtx[idx];
	}
	return nidx;
}

/************************************************************************
* mode = 0: simple finite differences
* mode = 1: smoothed finite differences (Bspline?)
* Default:  all gradients are (0, 0, 1)                                                                 
************************************************************************/
template <class T>
void Reg3Data<T>::calcGradient(int mode)
{
	if (hasGradient()) return;

	// Use B-spline approximation to calculate gradients
	float v[27];
	int ix[3], iy[3], iz[3];
	int i, j, k, l, m, n, t;

	p_grad = (float (*)[3]) malloc(sizeof(float[3]) * getNVerts()); 
	memset(p_grad, 0, sizeof(float)*3*getNVerts());
	switch (mode) {
	case 0:
		for (k = 0; k < m_dim[2]; k++) {
			for (j = 0; j < m_dim[1]; j++) {
				for (i = 0; i < m_dim[0]; i++) {
					int vid = index2vert(i, j, k);
					if (i==0) {
						// use right difference
						p_grad[vid][0] = getValue(vid+1) - getValue(vid);
					} else if (i == m_dim[0]-1) {
						// use left difference
						p_grad[vid][0] = getValue(vid) - getValue(vid-1);
					} else {
						// use central difference
						p_grad[vid][0] = (getValue(vid+1) - getValue(vid-1)) / 2;
					}

					if (j == 0) {
						p_grad[vid][1] = getValue(vid+m_dim[0]) - getValue(vid);
					} else if (j  ==  m_dim[1]-1) {
						p_grad[vid][1] = getValue(vid) - getValue(vid-m_dim[0]);
					} else {
						p_grad[vid][1] = (getValue(vid+m_dim[0]) - getValue(vid-m_dim[0])) / 2;
					}

					if (k == 0) {
						p_grad[vid][2] = getValue(vid+m_dim[0]*m_dim[1]) - getValue(vid);
					} else if (k == m_dim[2]-1) {
						p_grad[vid][2] = getValue(vid) - getValue(vid-m_dim[1]*m_dim[0]);
					} else {
						p_grad[vid][2] = (getValue(vid+m_dim[1]*m_dim[0]) - getValue(vid-m_dim[1]*m_dim[0])) / 2;
					}

					p_grad[vid][0] /= m_span[0];
					p_grad[vid][1] /= m_span[1];
					p_grad[vid][2] /= m_span[2];
				}
			}
		}

		break;
	case 1:
		for (k = 0; k < m_dim[2]; k++)
			for (j = 0; j < m_dim[1]; j++)
				for (i = 0; i < m_dim[0]; i++) {
					ix[0] = (i-1 >= 0)? i-1:0;
					ix[1] = i;
					ix[2] = (i+1 < m_dim[0])? i+1:i;
					iy[0] = (j-1 >= 0)? j-1:0;
					iy[1] = j;
					iy[2] = (j+1 < m_dim[1])? j+1:j;
					iz[0] = (k-1 >= 0)? k-1:0;
					iz[1] = k;
					iz[2] = (k+1 < m_dim[2])? k+1:k;

					t = 0;
					for (n = 0; n < 3; n++) {
						for (m = 0; m < 3; m++) {
							for (l = 0; l < 3; l++) {
								v[t] = getValue(ix[l], iy[m], iz[n]);
								t++;
							}
						}
					}

					int vid = index2vert(i, j, k);
					for (l = 0; l < 27; l++) {
						p_grad[vid][0] += x_grad_mask[l]*v[l];
						p_grad[vid][1] += y_grad_mask[l]*v[l];
						p_grad[vid][2] += z_grad_mask[l]*v[l];
					}
					p_grad[vid][0] /= m_span[0];
					p_grad[vid][1] /= m_span[1];
					p_grad[vid][2] /= m_span[2];
				}
		break;
	default:
		int vid = 0;
		for (k = 0; k < m_dim[2]; k++) {
			for (j = 0; j < m_dim[1]; j++) {
				for (i = 0; i < m_dim[0]; i++) {
					p_grad[vid][0] = 0;
					p_grad[vid][1] = 0;
					p_grad[vid][2] = 1;
					vid++;
				}
			}
		}
	}
}
#endif


