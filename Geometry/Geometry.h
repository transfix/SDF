/*
  Copyright 2002-2003 The University of Texas at Austin
  
    Authors: Anthony Thane <thanea@ices.utexas.edu>
             Vinay Siddavanahalli <skvinay@cs.utexas.edu>
    Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// Geometry.h: interface for the Geometry class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GEOMETRY_H__A99E4F7E_DEED_48DB_9F49_0DBB21EC4211__INCLUDED_)
#define AFX_GEOMETRY_H__A99E4F7E_DEED_48DB_9F49_0DBB21EC4211__INCLUDED_

#include <vector>
using std::vector;

#include <cstdio>

class Geometry  
{
public:
	Geometry();
	Geometry(const Geometry& copy);
	Geometry& operator=(const Geometry& copy);
	virtual ~Geometry();

	void ZeroMembers();
	void InitializeColors();

	void CopyPoints(const Geometry& copy);
	void CopyPointColors(const Geometry& copy);
	void CopyLines(const Geometry& copy);
	void CopyLineColors(const Geometry& copy);
	void CopyQuads(const Geometry& copy);
	void CopyQuadTexCoords(const Geometry& copy);
	void CopyTris(const Geometry& copy);
	void CopyTriTexCoords(const Geometry& copy);
	void CopyTriVertColors(const Geometry& copy);
	void CopyColors(const Geometry& copy);

	void ClearGeometry();
	void ClearPoints();
	void ClearPointColors();
	void ClearLines();
	void ClearLineColors();
	void ClearQuads();
	void ClearQuadTexCoords();
	void ClearTris();
	void ClearTriTexCoords();
	void ClearTriTexCoords2D();
	void ClearTriVertColors();

	void AllocatePoints(unsigned int NumPoints);
	void AllocatePointColors();
	void AllocateLines(unsigned int NumLineVerts, unsigned int NumLines);
	void AllocateLineColors();
	void AllocateQuads(unsigned int NumQuadVerts, unsigned int NumQuads);
	void AllocateQuadTexCoords();
	void AllocateTris(unsigned int NumTriVerts, unsigned int NumTris);
	void AllocateTriTexCoords();
	void AllocateTriTexCoords2D();
	void AllocateTriVertColors();
	void allocateTriangleDerivates();

	void CalculateQuadSmoothNormals();
	void CalculateTriSmoothNormals();
	void CalculateQuadFlatNormals();
	void CalculateTriFlatNormals();
	void SetTriNormalsReady();
	void SetQuadNormalsReady();
	void computeDerivatives();
	void computeTriangleMeshDerivatives();
	bool printTriangleDerivatives( std::FILE* fp );

	void setTriangleColorsByNormals();
	void setTriangleColorsByMeanCurv();
	void setTriangleColorsByGaussianCurv();
	void setTriangleColors( int renderingMode );

	bool printDerivatives( const char* filename );
	void setColors( int renderingMode );

	void SetDiffusedColor( float r, float g, float b );
	void SetSpecularColor( float r, float g, float b );
	void SetAmbientColor( float r, float g, float b );
	void SetShininess( float s );

	void SetWireframeWidth( float wireframeWidth );
	void setWireframeColor( float r, float g, float b );

	void setPointSize( float pointSize );
	void setPointColor( float r, float g, float b );
	bool useUniquePointColors();

	void SetLineWidth( float lineWidth );
	void setLineColor( float r, float g, float b );
	bool useUniqueLineColors();

	bool useUniqueTriangleColors();
	bool useUniqueWireframeColors();

	void GetReadyToDrawWire();
	void GetReadyToDrawSmooth();
	void GetReadyToDrawFlat();

	void getExtents( double* minx, double* miny, double* minz, double* maxx, double* maxy, double* maxz);

	bool separateComponents( vector<Geometry *>* components );
	Geometry* merge( Geometry* geometry );

	void print();
	static bool getPointInTriangle( float x1, float y1, float z1, 
									float x2, float y2, float z2, 
									float x3, float y3, float z3, 
									float *xinterp, float *yinterp, float *zinterp );
	float m_DiffuseColor[4];
	float m_SpecularColor[4];
	float m_AmbientColor[4];
	float m_Shininess;


	float* m_Points;
	unsigned int m_NumPoints;
	float m_UniquePointColors[3];
	float* m_PointColorsArray;
	float m_PointSize;
	bool m_RenderPoints;
	bool m_UsePointColors;

	float* m_LineVerts;
	float m_UniqueLineColors[3];
	float* m_LineColors;
	unsigned int m_NumLineVerts;
	unsigned int* m_Lines;
	unsigned int m_NumLines;
	float m_LineWidth;
	bool m_RenderLines;
	bool m_UseLineColors;

	float* m_QuadVerts;
	float* m_QuadVertNormals;
	float* m_QuadVertTexCoords;
	unsigned int m_NumQuadVerts;
	unsigned int* m_Quads;
	unsigned int m_NumQuads;
	float* m_QuadFlatVerts;
	float* m_QuadFlatNormals;
	float* m_QuadFlatTexCoords;
	bool m_bQuadFlatNormalsReady;
	bool m_bQuadSmoothNormalsReady;
	bool m_RenderQuads;

	float* m_TriVerts;
	float* m_TriVertNormals;
	float* m_TriVertColors;
	float* m_TriVertTexCoords;
	float* m_TriVertTexCoords2D;
	unsigned int m_NumTriVerts;
	unsigned int* m_Tris;
	unsigned int m_NumTris;
	float* m_TriFlatVerts;
	float* m_TriFlatNormals;
	float* m_TriFlatTexCoords;
	bool m_bTriFlatNormalsReady;
	bool m_bTriSmoothNormalsReady;
	bool m_RenderTriangles;
	bool m_UseTriangleColors;
	bool m_UseWireframeColors;
	float* m_TriMeanCurv;
	float* m_TriGaussianCurv;
	float* m_K1;
	float* m_K2;
	bool m_TriangleDerivativesValid;

	float m_WireframeWidth;
	float m_UniqueWireframeColors[3];
	bool m_RenderWireframe;

	float m_Center[3];
	float m_Min[3];
	float m_Max[3];
	bool m_bExtentsReady;
	void CalculateExtents();

protected:
	void mergeComponents( int* componentIndex, int v1, int v2, int n );
	bool getNextComponentSize( int* numVerts, int* numTris, int* nextComponentIndex, int* componentIndex );
	bool addComponent( Geometry* component, int nextComponentIndex, int* componentIndex );
	void resetUsedComponentIndex( int nextComponentIndex, int* componentIndex );

private:
	void CalculateTriangleNormal(float* norm, unsigned int v0, unsigned int v1, unsigned int v2);
	void CalculateQuadNormal(float* norm, unsigned int v0, unsigned int v1, unsigned int v2);
	
};

#endif // !defined(AFX_GEOMETRY_H__A99E4F7E_DEED_48DB_9F49_0DBB21EC4211__INCLUDED_)
