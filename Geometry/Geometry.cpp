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

// Geometry.cpp: implementation of the Geometry class.
//
//////////////////////////////////////////////////////////////////////

#include "Geometry.h"
#include <math.h>
#include "MeshDerivatives.h"

#include <cstdlib>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Geometry::Geometry()
{
	ZeroMembers();
	InitializeColors();
}

Geometry::~Geometry()
{
	ClearGeometry();
}

void Geometry::ZeroMembers()
{
	m_Points = 0;
	m_NumPoints = 0;
	m_UniquePointColors[0] = 0;
	m_UniquePointColors[1] = 0;
	m_UniquePointColors[2] = 0;
	m_PointColorsArray = 0;
	m_PointSize = 1.0;
	m_RenderPoints = false;
	m_UsePointColors = false;

	m_LineVerts = 0;
	m_UniqueLineColors[0] = 0.5;
	m_UniqueLineColors[1] = 0.5;
	m_UniqueLineColors[2] = 0.5;
	m_LineColors = 0;
	m_NumLineVerts = 0;
	m_Lines = 0;
	m_NumLines = 0;
	m_RenderLines = false;
	m_LineWidth = 1;
	m_UseLineColors = false;

	m_QuadVerts = 0;
	m_QuadVertNormals = 0;
	m_QuadVertTexCoords = 0;
	m_NumQuadVerts = 0;
	m_Quads = 0;
	m_NumQuads = 0;
	m_QuadFlatVerts = 0;
	m_QuadFlatNormals = 0;
	m_QuadFlatTexCoords = 0;
	m_bQuadFlatNormalsReady = 0;
	m_bQuadSmoothNormalsReady = 0;
	m_RenderQuads = false;

	m_TriVerts = 0;
	m_TriVertNormals = 0;
	m_TriVertTexCoords = 0;
	m_TriVertTexCoords2D = 0;
	m_TriVertColors = 0;
	m_NumTriVerts = 0;
	m_Tris = 0;
	m_NumTris = 0;
	m_TriFlatVerts = 0;
	m_TriFlatNormals = 0;
	m_TriFlatTexCoords = 0;
	m_bTriFlatNormalsReady = 0;
	m_bTriSmoothNormalsReady = 0;
	m_RenderTriangles = false;
	m_UseTriangleColors = false;
	m_UseWireframeColors = false;
	m_TriMeanCurv = 0;
	m_TriGaussianCurv = 0;
	m_K1 = 0;
	m_K2 = 0;
	m_TriangleDerivativesValid = false;

	m_UniqueWireframeColors[0] = 0.5;
	m_UniqueWireframeColors[1] = 0.5;
	m_UniqueWireframeColors[2] = 0.5;
	m_RenderWireframe = false;

	m_bExtentsReady = false;
}

void Geometry::InitializeColors()
{
	m_WireframeWidth = 1.0;

	m_DiffuseColor[0] = 0.5;
	m_DiffuseColor[1] = 0.5;
	m_DiffuseColor[2] = 0.5;
	m_DiffuseColor[3] = 1.0;

	m_SpecularColor[0] = 1.0;
	m_SpecularColor[1] = 1.0;
	m_SpecularColor[2] = 1.0;
	m_SpecularColor[3] = 1.0;

	m_AmbientColor[0] = 0.0;
	m_AmbientColor[1] = 0.0;
	m_AmbientColor[2] = 0.0;
	m_AmbientColor[3] = 0.0;

	m_Shininess = 30.0;
}

Geometry::Geometry(const Geometry& copy)
{
	ZeroMembers();
	CopyPoints(copy);
	CopyPointColors(copy);
	CopyLines(copy);
	CopyLineColors(copy);
	CopyQuads(copy);
	CopyQuadTexCoords(copy);
	CopyTris(copy);
	CopyTriTexCoords(copy);
	CopyTriVertColors(copy);
	CopyColors(copy);
}

Geometry& Geometry::operator=(const Geometry& copy)
{
	if (this!=&copy) {
		ClearGeometry();
		CopyPoints(copy);
		CopyPointColors(copy);
		CopyLines(copy);
		CopyQuads(copy);
		CopyQuadTexCoords(copy);
		CopyTris(copy);
		CopyTriTexCoords(copy);
		CopyTriVertColors(copy);
		CopyColors(copy);
	}
	return *this;
}

void Geometry::CopyPoints(const Geometry& copy)
{
	AllocatePoints(copy.m_NumPoints);
	unsigned int c;
	for (c=0; c<m_NumPoints*3; c++) { // for every coordinate of every point
		m_Points[c] = copy.m_Points[c];
	}
	m_UniquePointColors[0] = copy.m_UniquePointColors[0];
	m_UniquePointColors[1] = copy.m_UniquePointColors[1];
	m_UniquePointColors[2] = copy.m_UniquePointColors[2];
	m_RenderPoints = copy.m_RenderPoints;
	m_UsePointColors = copy.m_UsePointColors;
}

void Geometry::CopyPointColors(const Geometry& copy)
{
	unsigned int c;
	if (copy.m_PointColorsArray) {
		AllocatePointColors();
		for (c=0; c<m_NumPoints*3; c++) {
			m_PointColorsArray[c] = copy.m_PointColorsArray[c];
		}
	}
	else {
		m_PointColorsArray = 0;
	}
}

void Geometry::CopyLines(const Geometry& copy)
{
	AllocateLines(copy.m_NumLineVerts, copy.m_NumLines);	
	unsigned int c;
	for (c=0; c<m_NumLineVerts*3; c++) { // for every coordinate of every LineVert
		m_LineVerts[c] = copy.m_LineVerts[c];
	}
	for (c=0; c<m_NumLines*2; c++) { // for every Vert of every Line
		m_Lines[c] = copy.m_Lines[c];
	}
	m_UniqueLineColors[0] = copy.m_UniqueLineColors[0];
	m_UniqueLineColors[1] = copy.m_UniqueLineColors[1];
	m_UniqueLineColors[2] = copy.m_UniqueLineColors[2];
	m_RenderLines = copy.m_RenderLines;
	m_LineWidth = copy.m_LineWidth;
	m_UseLineColors = copy.m_UseLineColors;
}

void Geometry::CopyLineColors(const Geometry& copy)
{
	unsigned int c;
	if (copy.m_LineColors) {
		AllocateLineColors();
		for (c=0; c<m_NumLineVerts*3; c++) {
			m_LineColors[c] = copy.m_LineColors[c];
		}
	}
	else {
		m_LineColors = 0;
	}
}

void Geometry::CopyQuads(const Geometry& copy)
{
	AllocateQuads(copy.m_NumQuadVerts, copy.m_NumQuads);	
	unsigned int c;
	for (c=0; c<m_NumQuadVerts*3; c++) { // for every coordinate of every QuadVert
		m_QuadVerts[c] = copy.m_QuadVerts[c];
	}
	for (c=0; c<m_NumQuadVerts*3; c++) { // for every coordinate of every QuadVertNormal
		m_QuadVertNormals[c] = copy.m_QuadVertNormals[c];
	}
	for (c=0; c<m_NumQuads*4; c++) { // for every Vert of every Quad
		m_Quads[c] = copy.m_Quads[c];
	}
	for (c=0; c<m_NumQuads*4*3; c++) { // for every coordinate of every vertex of every flat quad
		m_QuadFlatVerts[c] = copy.m_QuadFlatVerts[c];
		m_QuadFlatNormals[c] = copy.m_QuadFlatNormals[c];
	}
	m_bQuadFlatNormalsReady = copy.m_bQuadFlatNormalsReady;
	m_bQuadSmoothNormalsReady = copy.m_bQuadSmoothNormalsReady;
	m_RenderQuads = copy.m_RenderQuads;
}

void Geometry::CopyQuadTexCoords(const Geometry& copy)
{
	unsigned int c;
	if (copy.m_QuadVertTexCoords) {
		AllocateQuadTexCoords();
		for (c=0; c<m_NumQuadVerts*3; c++) {
			m_QuadVertTexCoords[c] = copy.m_QuadVertTexCoords[c];
		}
		for (c=0; c<m_NumQuads*4*3; c++) {
			m_QuadFlatTexCoords[c] = copy.m_QuadFlatTexCoords[c];
		}
	}
	else {
		m_QuadVertTexCoords = 0;
		m_QuadFlatTexCoords = 0;
	}
}

void Geometry::CopyTris(const Geometry& copy)
{
	AllocateTris(copy.m_NumTriVerts, copy.m_NumTris);
	unsigned int c;
	for (c=0; c<m_NumTriVerts*3; c++) { // for every coordinate of every TriVert
		m_TriVerts[c] = copy.m_TriVerts[c];
	}
	for (c=0; c<m_NumTriVerts*3; c++) { // for every coordinate of every TriVertNormal
		m_TriVertNormals[c] = copy.m_TriVertNormals[c];
	}
	for (c=0; c<m_NumTris*3; c++) { // for every Vert of every Tri
		m_Tris[c] = copy.m_Tris[c];
	}
	for (c=0; c<m_NumTris*3*3; c++) { // for every coordinate of every vertex of every flat Tri
		m_TriFlatVerts[c] = copy.m_TriFlatVerts[c];
		m_TriFlatNormals[c] = copy.m_TriFlatNormals[c];
	}
	m_bTriFlatNormalsReady = copy.m_bTriFlatNormalsReady;
	m_bTriSmoothNormalsReady = copy.m_bTriSmoothNormalsReady;

	m_UniqueWireframeColors[0] = copy.m_UniqueWireframeColors[0];
	m_UniqueWireframeColors[1] = copy.m_UniqueWireframeColors[1];
	m_UniqueWireframeColors[2] = copy.m_UniqueWireframeColors[2];
	m_UseWireframeColors = copy.m_UseWireframeColors;
	m_RenderTriangles = copy.m_RenderTriangles;
	m_RenderWireframe = copy.m_RenderWireframe;

	m_UseTriangleColors = copy.m_UseTriangleColors;

	if( copy.m_TriMeanCurv )
	{
		if( m_TriMeanCurv ) delete []m_TriMeanCurv; m_TriMeanCurv = 0;
		m_TriMeanCurv = new float[m_NumTriVerts];
		for (c=0; c<m_NumTriVerts; c++) { 
			m_TriMeanCurv[c] = copy.m_TriMeanCurv[c];
		}
	}
	if( copy.m_TriGaussianCurv )
	{
		if( m_TriGaussianCurv ) delete []m_TriGaussianCurv; m_TriGaussianCurv = 0;
		m_TriGaussianCurv = new float[m_NumTriVerts];
		for (c=0; c<m_NumTriVerts; c++) { 
			m_TriGaussianCurv[c] = copy.m_TriGaussianCurv[c];
		}
	}
	if( copy.m_K1 )
	{
		if( m_K1 ) delete []m_K1; m_K1 = 0;
		m_K1 = new float[m_NumTriVerts*3];
		for (c=0; c<m_NumTriVerts*3; c++) { // for every coordinate of every TriVertNormal
			m_K1[c] = copy.m_K1[c];
		}
	}
	if( copy.m_K2 )
	{
		if( m_K2 ) delete []m_K2; m_K2 = 0;
		m_K2 = new float[m_NumTriVerts*3];
		for (c=0; c<m_NumTriVerts*3; c++) { // for every coordinate of every TriVertNormal
			m_K2[c] = copy.m_K2[c];
		}
	}

	m_TriangleDerivativesValid = copy.m_TriangleDerivativesValid;
}

void Geometry::CopyTriTexCoords(const Geometry& copy)
{
	unsigned int c;
	if (copy.m_TriVertTexCoords) {
		AllocateTriTexCoords();
		for (c=0; c<m_NumTriVerts*3; c++) {
			m_TriVertTexCoords[c] = copy.m_TriVertTexCoords[c];
		}
		for (c=0; c<m_NumTris*3*3; c++) {
			m_TriFlatTexCoords[c] = copy.m_TriFlatTexCoords[c];
		}
	}
	else {
		m_TriVertTexCoords = 0;
		m_TriFlatTexCoords = 0;
	}
	if (copy.m_TriVertTexCoords2D) {
		AllocateTriTexCoords2D();
		for (c=0; c<m_NumTriVerts*2; c++) {
			m_TriVertTexCoords[c] = copy.m_TriVertTexCoords[c];
		}
	}
	else {
		m_TriVertTexCoords2D = 0;
	}
}

void Geometry::CopyTriVertColors(const Geometry& copy)
{
	unsigned int c;
	if (copy.m_TriVertColors) {
		AllocateTriVertColors();
		for (c=0; c<m_NumTriVerts*3; c++) {
			m_TriVertColors[c] = copy.m_TriVertColors[c];
		}
	}
	else {
		m_TriVertColors = 0;
	}
}

void Geometry::CopyColors(const Geometry& copy)
{
	m_WireframeWidth = copy.m_WireframeWidth;

	m_DiffuseColor[0] = copy.m_DiffuseColor[0];
	m_DiffuseColor[1] = copy.m_DiffuseColor[1];
	m_DiffuseColor[2] = copy.m_DiffuseColor[2];
	m_DiffuseColor[3] = copy.m_DiffuseColor[3];

	m_SpecularColor[0] = copy.m_SpecularColor[0];
	m_SpecularColor[1] = copy.m_SpecularColor[1];
	m_SpecularColor[2] = copy.m_SpecularColor[2];
	m_SpecularColor[3] = copy.m_SpecularColor[3];

	m_AmbientColor[0] = copy.m_AmbientColor[0];
	m_AmbientColor[1] = copy.m_AmbientColor[1];
	m_AmbientColor[2] = copy.m_AmbientColor[2];
	m_AmbientColor[3] = copy.m_AmbientColor[3];

	m_Shininess = copy.m_Shininess;
}

void Geometry::ClearGeometry()
{
	ClearPoints();
	ClearPointColors();
	ClearLines();
	ClearLineColors();
	ClearQuads();
	ClearQuadTexCoords();
	ClearTris();
	ClearTriTexCoords();
	ClearTriTexCoords2D();
	ClearTriVertColors();
	m_bExtentsReady = false;
}

void Geometry::ClearPoints()
{
	delete [] m_Points;
	m_Points = 0;
	m_NumPoints = 0;
}

void Geometry::ClearPointColors()
{
	if( m_PointColorsArray ) { delete [] m_PointColorsArray; m_PointColorsArray = 0; }
}

void Geometry::ClearLines()
{
	delete [] m_LineVerts;
	m_LineVerts = 0;
	m_NumLineVerts = 0;
	delete [] m_Lines;
	m_Lines = 0;
	m_NumLines = 0;
}

void Geometry::ClearLineColors()
{
	delete [] m_LineColors;
	m_LineColors = 0;
}

void Geometry::ClearQuads()
{
	delete [] m_QuadVerts;
	m_QuadVerts = 0;
	delete [] m_QuadVertNormals;
	m_QuadVertNormals = 0;
	m_NumQuadVerts = 0;
	delete [] m_Quads;
	m_Quads = 0;
	m_NumQuads = 0;
	delete [] m_QuadFlatVerts;
	m_QuadFlatVerts = 0;
	delete [] m_QuadFlatNormals;
	m_QuadFlatNormals = 0;
	m_bQuadFlatNormalsReady = false;
	m_bQuadSmoothNormalsReady = false;
}

void Geometry::ClearQuadTexCoords()
{
	delete [] m_QuadVertTexCoords;
	m_QuadVertTexCoords = 0;
	delete [] m_QuadFlatTexCoords;
	m_QuadFlatTexCoords = 0;
}

void Geometry::ClearTris()
{
	delete [] m_TriVerts;
	m_TriVerts = 0;
	delete [] m_TriVertNormals;
	m_TriVertNormals = 0;
	m_NumTriVerts = 0;
	delete [] m_Tris;
	m_Tris = 0;
	m_NumTris = 0;
	delete [] m_TriFlatVerts;
	m_TriFlatVerts = 0;
	delete [] m_TriFlatNormals;
	m_TriFlatNormals = 0;
	m_bTriFlatNormalsReady = false;
	m_bTriSmoothNormalsReady = false;

	delete [] m_TriMeanCurv;
	m_TriMeanCurv = 0;
	delete [] m_TriGaussianCurv;
	m_TriGaussianCurv = 0;
	delete [] m_TriMeanCurv;
	m_TriMeanCurv = 0;
	delete [] m_K1;
	m_K1 = 0;
	delete [] m_K2;
	m_K2 = 0;
}

void Geometry::ClearTriTexCoords()
{
	delete [] m_TriVertTexCoords;
	m_TriVertTexCoords = 0;
	delete [] m_TriFlatTexCoords;
	m_TriFlatTexCoords = 0;
}

void Geometry::ClearTriTexCoords2D()
{
	delete [] m_TriVertTexCoords2D;
	m_TriVertTexCoords2D = 0;
}

void Geometry::ClearTriVertColors()
{
	delete [] m_TriVertColors;
	m_TriVertColors = 0;
}

void Geometry::AllocatePoints(unsigned int NumPoints)
{
	ClearPoints();
	m_Points = new float[NumPoints*3];
	m_NumPoints = NumPoints;
	m_bExtentsReady = false;
}

void Geometry::AllocatePointColors()
{
	ClearPointColors();
	m_PointColorsArray = new float[m_NumPoints*3];
	unsigned int c;
	for (c=0; c<m_NumPoints; c++) {
		m_PointColorsArray[c*3+0] = m_DiffuseColor[0];
		m_PointColorsArray[c*3+1] = m_DiffuseColor[1];
		m_PointColorsArray[c*3+2] = m_DiffuseColor[2];
	}
}

void Geometry::AllocateLines(unsigned int NumLineVerts, unsigned int NumLines)
{
	ClearLines();
	m_LineVerts = new float[NumLineVerts*3];
	m_NumLineVerts = NumLineVerts;
	m_Lines = new unsigned int[NumLines*2];
	m_NumLines = NumLines;
	m_bExtentsReady = false;
}

void Geometry::AllocateLineColors()
{
	ClearLineColors();
	m_LineColors = new float[m_NumLineVerts*3];
	unsigned int c;
	for (c=0; c<m_NumLineVerts; c++) {
		m_LineColors[c*3+0] = m_DiffuseColor[0];
		m_LineColors[c*3+1] = m_DiffuseColor[1];
		m_LineColors[c*3+2] = m_DiffuseColor[2];
	}
}

void Geometry::AllocateQuads(unsigned int NumQuadVerts, unsigned int NumQuads)
{
	ClearQuads();
	m_QuadVerts = new float[NumQuadVerts*3];
	m_QuadVertNormals = new float[NumQuadVerts*3];
	m_NumQuadVerts = NumQuadVerts;
	m_Quads = new unsigned int[NumQuads*4];
	m_NumQuads = NumQuads;
	m_QuadFlatVerts = new float[NumQuads*4*3];
	m_QuadFlatNormals = new float[NumQuads*4*3];
	unsigned int c;
	for (c=0; c<NumQuads*4*3; c++) {
		m_QuadFlatVerts[c] = 0.0;
		m_QuadFlatNormals[c] = 0.0;
	}
	m_bExtentsReady = false;
}

void Geometry::AllocateQuadTexCoords()
{
	ClearQuadTexCoords();
	m_QuadVertTexCoords = new float[m_NumQuadVerts*3];
	unsigned int c;
	for (c=0; c<m_NumQuadVerts*3; c++) {
		m_QuadVertTexCoords[c] = 0.0f;
	}
	m_QuadFlatTexCoords = new float[m_NumQuads*4*3];
	for (c=0; c<m_NumQuads*4*3; c++) {
		m_QuadFlatTexCoords[c] = 0.0f;
	}
}

void Geometry::AllocateTris(unsigned int NumTriVerts, unsigned int NumTris)
{
	ClearTris();
	m_TriVerts = new float[NumTriVerts*3];
	m_TriVertNormals = new float[NumTriVerts*3];
	m_NumTriVerts = NumTriVerts;
	m_Tris = new unsigned int[NumTris*3];
	m_NumTris = NumTris;
	m_TriFlatVerts = new float[NumTris*3*3];
	m_TriFlatNormals = new float[NumTris*3*3];
	unsigned int c;
	for (c=0; c<NumTris*3*3; c++) {
		m_TriFlatVerts[c] = 0.0;
		m_TriFlatNormals[c] = 0.0;
	}
	m_bExtentsReady = false;
}

void Geometry::AllocateTriTexCoords()
{
	ClearTriTexCoords();
	m_TriVertTexCoords = new float[m_NumTriVerts*3];
	unsigned int c;
	for (c=0; c<m_NumTriVerts*3; c++) {
		m_TriVertTexCoords[c] = 0.0f;
	}
	m_TriFlatTexCoords = new float[m_NumTris*3*3];
	for (c=0; c<m_NumTris*3*3; c++) {
		m_TriFlatTexCoords[c] = 0.0f;
	}
}

void Geometry::AllocateTriTexCoords2D()
{
	ClearTriTexCoords2D();
	m_TriVertTexCoords2D = new float[m_NumTriVerts*2];

	if( !m_bExtentsReady ) CalculateExtents();

	unsigned int c;
	for (c=0; c<m_NumTriVerts; c++) 
	{
		m_TriVertTexCoords2D[c*2+0] = (m_TriVerts[c*3+0] - m_Min[0]) / (m_Max[0] - m_Min[0]);
		m_TriVertTexCoords2D[c*2+1] = (m_TriVerts[c*3+1] - m_Min[1]) / (m_Max[1] - m_Min[1]);
	}
}

void Geometry::AllocateTriVertColors()
{
	ClearTriVertColors();
	m_TriVertColors = new float[m_NumTriVerts*3];
	unsigned int c;
	for (c=0; c<m_NumTriVerts*3; c++) {
		m_TriVertColors[c] = 0.0f;
	}
}

static void cross(float* dest, const float* v1, const float* v2)
{

	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

static void normalize(float* v)
{
	float len = (float)sqrt(
		v[0] * v[0] +
		v[1] * v[1] +
		v[2] * v[2]);
	if (len!=0.0) {
		v[0]/=len; //*biggestDim;
		v[1]/=len; //*biggestDim;
		v[2]/=len; //*biggestDim;
	}
	else {
		v[0] = 1.0;
	}
}

static void add(float* dest, const float* v)
{
	dest[0] += v[0];
	dest[1] += v[1];
	dest[2] += v[2];
}

static void set(float* dest, const float* v)
{
	dest[0] = v[0];
	dest[1] = v[1];
	dest[2] = v[2];
}

void Geometry::CalculateQuadSmoothNormals()
{
	unsigned int c, v0, v1, v2, v3;
	float normal[3];

	float zero[3] = {0.0f, 0.0f, 0.0f};

	for (c=0; c<m_NumQuadVerts; c++) {
		set(m_QuadVertNormals+c*3, zero);
	}

		
	// for each Quadangle
	for (c=0; c<m_NumQuads; c++) {
		v0 = m_Quads[c*4+0];
		v1 = m_Quads[c*4+1];
		v2 = m_Quads[c*4+2];
		v3 = m_Quads[c*4+3];
		CalculateQuadNormal(normal, v0, v1, v3);
		add(m_QuadVertNormals+v0*3, normal);
		add(m_QuadVertNormals+v1*3, normal);
		add(m_QuadVertNormals+v3*3, normal);
		CalculateQuadNormal(normal, v2, v3, v1);
		add(m_QuadVertNormals+v2*3, normal);
		add(m_QuadVertNormals+v3*3, normal);
		add(m_QuadVertNormals+v1*3, normal);
	}
	
	// normalize the vectors
	for (c=0; c<m_NumQuadVerts; c++) {
		normalize(m_QuadVertNormals+c*3);
	}
	m_bQuadSmoothNormalsReady = true;
}

void Geometry::CalculateTriSmoothNormals()
{
	unsigned int c, v0, v1, v2;
	float normal[3];
		
	float zero[3] = {0.0f, 0.0f, 0.0f};

	for (c=0; c<m_NumTriVerts; c++) {
		set(m_TriVertNormals+c*3, zero);
	}

	// for each triangle
	for (c=0; c<m_NumTris; c++) {
		v0 = m_Tris[c*3+0];
		v1 = m_Tris[c*3+1];
		v2 = m_Tris[c*3+2];
		CalculateTriangleNormal(normal, v0, v1, v2);
		add(m_TriVertNormals+v0*3, normal);
		add(m_TriVertNormals+v1*3, normal);
		add(m_TriVertNormals+v2*3, normal);
	}
	
	// normalize the vectors
	for (c=0; c<m_NumTriVerts; c++) {
		normalize(m_TriVertNormals+c*3);
	}
	m_bTriSmoothNormalsReady = true;
}

void Geometry::CalculateQuadFlatNormals()
{
	unsigned int c, v0, v1, v2, v3;
	float normal1[3], normal2[3];
	for (c=0; c<m_NumQuads; c++) { // for every quad
		v0 = m_Quads[c*4+0];
		v1 = m_Quads[c*4+1];
		v2 = m_Quads[c*4+2];
		v3 = m_Quads[c*4+3];
		CalculateQuadNormal(normal1, v0, v1, v3);
		CalculateQuadNormal(normal2, v2, v3, v1);
		add(normal1, normal2);
		normalize(normal1);
		set(m_QuadFlatVerts+c*12+0*3, m_QuadVerts+v0*3); // set the verts
		set(m_QuadFlatVerts+c*12+1*3, m_QuadVerts+v1*3);
		set(m_QuadFlatVerts+c*12+2*3, m_QuadVerts+v2*3);
		set(m_QuadFlatVerts+c*12+3*3, m_QuadVerts+v3*3);
		if (m_QuadVertTexCoords) {
			set(m_QuadFlatTexCoords+c*12+0*3, m_QuadVertTexCoords+v0*3); // set the tex coords
			set(m_QuadFlatTexCoords+c*12+1*3, m_QuadVertTexCoords+v1*3);
			set(m_QuadFlatTexCoords+c*12+2*3, m_QuadVertTexCoords+v2*3);
			set(m_QuadFlatTexCoords+c*12+3*3, m_QuadVertTexCoords+v3*3);
		}
		set(m_QuadFlatNormals+c*12+0*3, normal1); // set the normals
		set(m_QuadFlatNormals+c*12+1*3, normal1);
		set(m_QuadFlatNormals+c*12+2*3, normal1);
		set(m_QuadFlatNormals+c*12+3*3, normal1);
	}
	m_bQuadFlatNormalsReady = true;
}

void Geometry::CalculateTriFlatNormals()
{
	unsigned int c, v0, v1, v2;
	float normal[3];
	for (c=0; c<m_NumTris; c++) { // for every triangle
		v0 = m_Tris[c*3+0];
		v1 = m_Tris[c*3+1];
		v2 = m_Tris[c*3+2];
		CalculateTriangleNormal(normal, v0, v1, v2);
		normalize(normal);
		set(m_TriFlatVerts+c*9+0*3, m_TriVerts+v0*3); // set the verts
		set(m_TriFlatVerts+c*9+1*3, m_TriVerts+v1*3);
		set(m_TriFlatVerts+c*9+2*3, m_TriVerts+v2*3);
		if (m_TriVertTexCoords) {
			set(m_TriFlatTexCoords+c*9+0*3, m_TriVertTexCoords+v0*3); // set the tex coords
			set(m_TriFlatTexCoords+c*9+1*3, m_TriVertTexCoords+v1*3);
			set(m_TriFlatTexCoords+c*9+2*3, m_TriVertTexCoords+v2*3);
		}
		set(m_TriFlatNormals+c*9+0*3, normal); // set the normals
		set(m_TriFlatNormals+c*9+1*3, normal);
		set(m_TriFlatNormals+c*9+2*3, normal);
	}
	m_bTriFlatNormalsReady = true;
}

void Geometry::SetTriNormalsReady()
{
	m_bTriSmoothNormalsReady = true;
}

void Geometry::SetQuadNormalsReady()
{
	m_bQuadSmoothNormalsReady = true;
}

void Geometry::SetDiffusedColor( float r, float g, float b )
{
	m_DiffuseColor[0] = r;
	m_DiffuseColor[1] = g;
	m_DiffuseColor[2] = b;
}

void Geometry::SetSpecularColor( float r, float g, float b )
{
	m_SpecularColor[0] = r;
	m_SpecularColor[1] = g;
	m_SpecularColor[2] = b;
}

void Geometry::SetAmbientColor( float r, float g, float b )
{
	m_AmbientColor[0] = r;
	m_AmbientColor[1] = g;
	m_AmbientColor[2] = b;
}

void Geometry::SetShininess( float s )
{
	m_Shininess = s;
}

void Geometry::SetWireframeWidth( float wireframeWidth )
{
	m_WireframeWidth = wireframeWidth;
}

void Geometry::GetReadyToDrawWire()
{
	if (!m_bExtentsReady) {
		CalculateExtents();
	}
}

void Geometry::GetReadyToDrawSmooth()
{
	if (!m_bExtentsReady) {
		CalculateExtents();
	}
	if ((m_NumQuads > 0) && (!m_bQuadSmoothNormalsReady)) {
		CalculateQuadSmoothNormals();
	}
	if ((m_NumTris > 0) && (!m_bTriSmoothNormalsReady)) {
		CalculateTriSmoothNormals();
	}
	GetReadyToDrawFlat();
}

void Geometry::GetReadyToDrawFlat()
{
	if (!m_bExtentsReady) {
		CalculateExtents();
	}
	if ((m_NumQuads > 0) && (!m_bQuadFlatNormalsReady)) {
		CalculateQuadFlatNormals();
	}
	if ((m_NumTris > 0) && (!m_bTriFlatNormalsReady)) {
		CalculateTriFlatNormals();
	}
}

static void CheckMin(bool& initialized, float* CurrentMin, float* Check)
{
	if (initialized) {
		CurrentMin[0] = (CurrentMin[0]<Check[0]?CurrentMin[0]:Check[0]);
		CurrentMin[1] = (CurrentMin[1]<Check[1]?CurrentMin[1]:Check[1]);
		CurrentMin[2] = (CurrentMin[2]<Check[2]?CurrentMin[2]:Check[2]);
	}
	else {
		CurrentMin[0] = Check[0];
		CurrentMin[1] = Check[1];
		CurrentMin[2] = Check[2];
		initialized = true;
	}
}

static void CheckMax(bool& initialized, float* CurrentMax, float* Check)
{
	if (initialized) {
		CurrentMax[0] = (CurrentMax[0]>Check[0]?CurrentMax[0]:Check[0]);
		CurrentMax[1] = (CurrentMax[1]>Check[1]?CurrentMax[1]:Check[1]);
		CurrentMax[2] = (CurrentMax[2]>Check[2]?CurrentMax[2]:Check[2]);
	}
	else {
		CurrentMax[0] = Check[0];
		CurrentMax[1] = Check[1];
		CurrentMax[2] = Check[2];
		initialized = true;
	}
}

void Geometry::CalculateExtents()
{
	unsigned int c;
	bool initialized = false;

	for(c=0;c<m_NumPoints;c++)
	{
		CheckMin(initialized, m_Min, m_Points+c*3);
		CheckMax(initialized, m_Max, m_Points+c*3);
	}
	for(c=0;c<m_NumLineVerts;c++)
	{
		CheckMin(initialized, m_Min, m_LineVerts+c*3);
		CheckMax(initialized, m_Max, m_LineVerts+c*3);
	}
	for(c=0;c<m_NumTriVerts;c++)
	{
		CheckMin(initialized, m_Min, m_TriVerts+c*3);
		CheckMax(initialized, m_Max, m_TriVerts+c*3);
	}
	for(c=0;c<m_NumQuadVerts;c++)
	{
		CheckMin(initialized, m_Min, m_QuadVerts+c*3);
		CheckMax(initialized, m_Max, m_QuadVerts+c*3);
	}

	m_Center[0]=( m_Max[0] + m_Min[0] ) / 2.0f;
	m_Center[1]=( m_Max[1] + m_Min[1] ) / 2.0f;
	m_Center[2]=( m_Max[2] + m_Min[2] ) / 2.0f;
	m_bExtentsReady = true;
}

void Geometry::CalculateTriangleNormal(float* norm, unsigned int v0, unsigned int v1, unsigned int v2)
{
	float vec1[3], vec2[3];
	vec1[0] = vec2[0] = -m_TriVerts[v0*3+0];
	vec1[1] = vec2[1] = -m_TriVerts[v0*3+1];
	vec1[2] = vec2[2] = -m_TriVerts[v0*3+2];
	
	vec1[0] += m_TriVerts[v1*3+0];
	vec1[1] += m_TriVerts[v1*3+1];
	vec1[2] += m_TriVerts[v1*3+2];
	
	vec2[0] += m_TriVerts[v2*3+0];
	vec2[1] += m_TriVerts[v2*3+1];
	vec2[2] += m_TriVerts[v2*3+2];
	
	cross(norm, vec1, vec2);
}

void Geometry::CalculateQuadNormal(float* norm, unsigned int v0, unsigned int v1, unsigned int v2)
{
	float vec1[3], vec2[3];
	vec1[0] = vec2[0] = -m_QuadVerts[v0*3+0];
	vec1[1] = vec2[1] = -m_QuadVerts[v0*3+1];
	vec1[2] = vec2[2] = -m_QuadVerts[v0*3+2];
	
	vec1[0] += m_QuadVerts[v1*3+0];
	vec1[1] += m_QuadVerts[v1*3+1];
	vec1[2] += m_QuadVerts[v1*3+2];
	
	vec2[0] += m_QuadVerts[v2*3+0];
	vec2[1] += m_QuadVerts[v2*3+1];
	vec2[2] += m_QuadVerts[v2*3+2];


	
	cross(norm, vec1, vec2);
}

/////////////////////// separate components ////////////////

void Geometry::mergeComponents( int* componentIndex, int v1, int v2, int n )
{
	if( !componentIndex || v1 >= n || v2 >= n ) return;

	if( componentIndex[v2] == -1 )
		componentIndex[v2] = componentIndex[v1];
	else if( componentIndex[v2] != componentIndex[v1] )
	{
		int j;
		for( j=0; j<n; j++ )
		{
			if( componentIndex[j] == componentIndex[v2] ) componentIndex[j] = componentIndex[v1];
		}
	}
}

bool Geometry::getNextComponentSize( int* numVerts, int* numTris, int* nextComponentIndex, int* componentIndex )
{
	// find the first index != -1. Find the number of vertices and triangles in it. Return the index also
	if( !numVerts || !numTris || !nextComponentIndex || !componentIndex ) return false;

	{
		int i;
		*nextComponentIndex = -1;
		for( i=0; i<m_NumTriVerts; i++ )
		{
			if( componentIndex[i] != -1 ) 
			{
				*nextComponentIndex = componentIndex[i];
				break;
			}
		}
		if( (*nextComponentIndex) == -1 ) return false;
	}

	*numVerts = 0;
	*numTris = 0;
	{
		// mark every triangle belonging to this component
		int i;
		for( i=0; i<m_NumTris; i++ )
		{
			if( componentIndex[m_Tris[3*i+0]] == (*nextComponentIndex) )
				(*numTris)++;
		}

		for( i=0; i<m_NumTriVerts; i++ )
		{
			if( componentIndex[i] == (*nextComponentIndex) ) (*numVerts)++;
		}
	}

	return true;
}

bool Geometry::addComponent( Geometry* component, int nextComponentIndex, int* componentIndex )
{
	if( !component ) return false;

	int* vertexIndex = new int[m_NumTriVerts];
	{
		// for every vertex in the new component, assign it an index
		int nV = m_NumTriVerts;
		int i;
		int nextVertexIndex = 0;

		for( i=0; i<nV; i++ ) 
		{
			vertexIndex[i] = -1;

			if( componentIndex[i] == nextComponentIndex )
			{
				vertexIndex[i] = nextVertexIndex;
				nextVertexIndex++;
			}
		}
	}
	{
		// assign triangles
		int i;
		int curTri = 0;
		for( i=0; i<m_NumTris; i++ )
		{
			if( componentIndex[m_Tris[3*i+0]] == nextComponentIndex )
			{
				int v1 = m_Tris[3*i+0];
				int v2 = m_Tris[3*i+1];
				int v3 = m_Tris[3*i+2];

				component->m_Tris[3*curTri+0] = vertexIndex[v1];
				component->m_Tris[3*curTri+1] = vertexIndex[v2];
				component->m_Tris[3*curTri+2] = vertexIndex[v3];

				curTri++;
			}
		}
	}
	{
		// assign vertices
		int nV = m_NumTriVerts;
		int i;
		for( i=0; i<nV; i++ )
		{
			if( vertexIndex[i] == -1 ) continue;

			component->m_TriVerts[3*vertexIndex[i] + 0] = m_TriVerts[3*i+0];
			component->m_TriVerts[3*vertexIndex[i] + 1] = m_TriVerts[3*i+1];
			component->m_TriVerts[3*vertexIndex[i] + 2] = m_TriVerts[3*i+2];
		}

	}
	return true;
}

void Geometry::resetUsedComponentIndex( int nextComponentIndex, int* componentIndex )
{
	// for all indices == nextComponentIndex, set it to -1;
	int i;
	for( i=0; i<m_NumTriVerts; i++ )
	{
		if( componentIndex[i] == nextComponentIndex )
			componentIndex[i] = -1;
	}
}

/*!
	Each component is written out as a file.
*/
bool Geometry::separateComponents( 	vector<Geometry *>* components )
{
	if(  !components ) return false;

	int* componentIndex = new int[m_NumTriVerts];
	int i;

	for(i=0; i<m_NumTriVerts; i++ ) componentIndex[i] = -1;
	
	int tempComponentIndex = 0;
	for(i=0; i<m_NumTris; i++ )
	{
		int v1 = m_Tris[3*i+0];
		int v2 = m_Tris[3*i+1];
		int v3 = m_Tris[3*i+2];

		if( componentIndex[v1] != -1 )
		{
			mergeComponents( componentIndex, v1, v2, m_NumTriVerts );
			mergeComponents( componentIndex, v1, v3, m_NumTriVerts );
		}
		if( componentIndex[v2] != -1 )
		{
			mergeComponents( componentIndex, v2, v1, m_NumTriVerts );
			mergeComponents( componentIndex, v2, v3, m_NumTriVerts );
		}
		if( componentIndex[v3] != -1 )
		{
			mergeComponents( componentIndex, v3, v1, m_NumTriVerts );
			mergeComponents( componentIndex, v3, v2, m_NumTriVerts );
		}
		if( componentIndex[v1] == -1 && componentIndex[v2] == -1 && componentIndex[v3] == -1 )
		{
			componentIndex[v1] = componentIndex[v2] = componentIndex[v3] = tempComponentIndex;
			tempComponentIndex++;
		}
	}

	int numberOfComponents = 0;
	while ( true )
	{
		int numVerts = 0;
		int numTris = 0;
		int nextComponentIndex = -1;

		if( ! getNextComponentSize( &numVerts, &numTris, &nextComponentIndex, componentIndex ) ) break;
		if( numVerts == 0 || numTris == 0 || nextComponentIndex == -1 ) break;

		Geometry* component = new Geometry();
		component->AllocateTris( numVerts, numTris );
		if( m_TriVertColors )
			component->AllocateTriVertColors();

		addComponent( component, nextComponentIndex, componentIndex );

		components->push_back( component );

		resetUsedComponentIndex( nextComponentIndex, componentIndex );
		numberOfComponents++;
	}

	delete []componentIndex;

	if( numberOfComponents == 0 ) return false;
	return true;
}

Geometry* Geometry::merge( Geometry* geometry )
{
	if( !geometry ) return 0;

	Geometry* mergedGeometry = new Geometry;

	mergedGeometry->AllocateTris( m_NumTriVerts + geometry->m_NumTriVerts, m_NumTris + geometry->m_NumTris );
	mergedGeometry->AllocateTriVertColors();

	int i;

	if( m_TriVerts && geometry->m_TriVerts )
	{
		for( i=0; i<m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriVerts[3*i+0] = m_TriVerts[3*i+0];
			mergedGeometry->m_TriVerts[3*i+1] = m_TriVerts[3*i+1];
			mergedGeometry->m_TriVerts[3*i+2] = m_TriVerts[3*i+2];
		}
		for( i=m_NumTriVerts; i<m_NumTriVerts+geometry->m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriVerts[3*i+0] = geometry->m_TriVerts[3*(i-m_NumTriVerts)+0];
			mergedGeometry->m_TriVerts[3*i+1] = geometry->m_TriVerts[3*(i-m_NumTriVerts)+1];
			mergedGeometry->m_TriVerts[3*i+2] = geometry->m_TriVerts[3*(i-m_NumTriVerts)+2];
		}
	}

	if( m_Tris && geometry->m_Tris )
	{
		for( i=0; i<m_NumTris; i++ )
		{
			mergedGeometry->m_Tris[3*i+0] = m_Tris[3*i+0];
			mergedGeometry->m_Tris[3*i+1] = m_Tris[3*i+1];
			mergedGeometry->m_Tris[3*i+2] = m_Tris[3*i+2];
		}
		for( i=m_NumTris; i<m_NumTris+geometry->m_NumTris; i++ )
		{
			mergedGeometry->m_Tris[3*i+0] = m_NumTriVerts+geometry->m_Tris[3*(i-m_NumTris)+0];
			mergedGeometry->m_Tris[3*i+1] = m_NumTriVerts+geometry->m_Tris[3*(i-m_NumTris)+1];
			mergedGeometry->m_Tris[3*i+2] = m_NumTriVerts+geometry->m_Tris[3*(i-m_NumTris)+2];
		}
	}

	if( m_TriVertNormals && geometry->m_TriVertNormals )
	{
		for( i=0; i<m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriVertNormals[3*i+0] = m_TriVertNormals[3*i+0];
			mergedGeometry->m_TriVertNormals[3*i+1] = m_TriVertNormals[3*i+1];
			mergedGeometry->m_TriVertNormals[3*i+2] = m_TriVertNormals[3*i+2];
		}
		for( i=m_NumTriVerts; i<m_NumTriVerts+geometry->m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriVertNormals[3*i+0] = geometry->m_TriVertNormals[3*(i-m_NumTriVerts)+0];
			mergedGeometry->m_TriVertNormals[3*i+1] = geometry->m_TriVertNormals[3*(i-m_NumTriVerts)+1];
			mergedGeometry->m_TriVertNormals[3*i+2] = geometry->m_TriVertNormals[3*(i-m_NumTriVerts)+2];
		}
	}

	if( m_TriMeanCurv && geometry->m_TriMeanCurv )
	{
		for( i=0; i<m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriMeanCurv[i] = m_TriMeanCurv[i];
		}
		for( i=m_NumTriVerts; i<m_NumTriVerts+geometry->m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriMeanCurv[i] = geometry->m_TriMeanCurv[(i-m_NumTriVerts)];
		}
	}
	if( m_TriGaussianCurv && geometry->m_TriGaussianCurv )
	{
		for( i=0; i<m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriGaussianCurv[i] = m_TriGaussianCurv[i];
		}
		for( i=m_NumTriVerts; i<m_NumTriVerts+geometry->m_NumTriVerts; i++ )
		{
			mergedGeometry->m_TriGaussianCurv[i] = geometry->m_TriGaussianCurv[(i-m_NumTriVerts)];
		}
	}
	if( m_K1 && geometry->m_K1 )
	{
		for( i=0; i<m_NumTriVerts*3; i++ )
		{
			mergedGeometry->m_K1[i] = m_K1[i];
		}
		for( i=m_NumTriVerts*3; i<(m_NumTriVerts+geometry->m_NumTriVerts)*3; i++ )
		{
			mergedGeometry->m_K1[i] = geometry->m_K1[(i-m_NumTriVerts*3)];
		}
	}
	if( m_K2 && geometry->m_K2 )
	{
		for( i=0; i<m_NumTriVerts*3; i++ )
		{
			mergedGeometry->m_K2[i] = m_K2[i];
		}
		for( i=m_NumTriVerts*3; i<(m_NumTriVerts+geometry->m_NumTriVerts)*3; i++ )
		{
			mergedGeometry->m_K2[i] = geometry->m_K2[(i-m_NumTriVerts*3)];
		}
	}

	if( !m_TriangleDerivativesValid || !geometry->m_TriangleDerivativesValid )
		m_TriangleDerivativesValid = false;

	return mergedGeometry;
}

void Geometry::getExtents( double* minx, double* miny, double* minz, double* maxx, double* maxy, double* maxz)
{
	if( !minx || !miny || !minz || !maxx || !maxy || !maxz ) return;

	int i;
	if( m_NumPoints )
	{
		*minx = *maxx = m_Points[0];
		*miny = *maxy = m_Points[1];
		*minz = *maxz = m_Points[2];

		for( i=1; i<m_NumPoints; i++ )
		{
			if( m_Points[3*i+0] < *minx ) *minx = m_Points[3*i+0];
			if( m_Points[3*i+1] < *miny ) *miny = m_Points[3*i+1];
			if( m_Points[3*i+2] < *minz ) *minz = m_Points[3*i+2];

			if( m_Points[3*i+0] > *maxx ) *maxx = m_Points[3*i+0];
			if( m_Points[3*i+1] > *maxy ) *maxy = m_Points[3*i+1];
			if( m_Points[3*i+2] > *maxz ) *maxz = m_Points[3*i+2];

		}
	}
	if( m_NumLineVerts )
	{
		*minx = *maxx = m_LineVerts[0];
		*miny = *maxy = m_LineVerts[1];
		*minz = *maxz = m_LineVerts[2];

		for( i=1; i<m_NumLineVerts; i++ )
		{
			if( m_LineVerts[3*i+0] < *minx ) *minx = m_LineVerts[3*i+0];
			if( m_LineVerts[3*i+1] < *miny ) *miny = m_LineVerts[3*i+1];
			if( m_LineVerts[3*i+2] < *minz ) *minz = m_LineVerts[3*i+2];

			if( m_LineVerts[3*i+0] > *maxx ) *maxx = m_LineVerts[3*i+0];
			if( m_LineVerts[3*i+1] > *maxy ) *maxy = m_LineVerts[3*i+1];
			if( m_LineVerts[3*i+2] > *maxz ) *maxz = m_LineVerts[3*i+2];

		}
	}
	if( m_NumTriVerts )
	{
		*minx = *maxx = m_TriVerts[0];
		*miny = *maxy = m_TriVerts[1];
		*minz = *maxz = m_TriVerts[2];

		for( i=1; i<m_NumTriVerts; i++ )
		{
			if( m_TriVerts[3*i+0] < *minx ) *minx = m_TriVerts[3*i+0];
			if( m_TriVerts[3*i+1] < *miny ) *miny = m_TriVerts[3*i+1];
			if( m_TriVerts[3*i+2] < *minz ) *minz = m_TriVerts[3*i+2];

			if( m_TriVerts[3*i+0] > *maxx ) *maxx = m_TriVerts[3*i+0];
			if( m_TriVerts[3*i+1] > *maxy ) *maxy = m_TriVerts[3*i+1];
			if( m_TriVerts[3*i+2] > *maxz ) *maxz = m_TriVerts[3*i+2];

		}
	}
}

void Geometry::print()
{
	int i;
	printf("Number of points: %u\n", m_NumPoints );

	for( i=0; i<m_NumPoints; i++ )
	{
		printf("%u [%f %f %f]\n", i, m_Points[3*i+0], m_Points[3*i+1], m_Points[3*i+2] );
	}
	
	printf("Number of lines: %u\n", m_NumLines );

	for( i=0; i<m_NumLines; i++ )
	{
		unsigned int idx1 = m_Lines[2*i+0];
		unsigned int idx2 = m_Lines[2*i+1];

		printf("%u [%u %u]: [%f %f %f] - [%f %f %f]\n", i, 
			idx1, idx2,
			m_LineVerts[3*idx1+0], m_LineVerts[3*idx1+1], m_LineVerts[3*idx1+2],
			m_LineVerts[3*idx2+0], m_LineVerts[3*idx2+1], m_LineVerts[3*idx2+2] 
			);
	}

	printf("Number of triangles: %u\n", m_NumTris );

	for( i=0; i<m_NumTris; i++ )
	{
		unsigned int idx1 = m_Tris[3*i+0];
		unsigned int idx2 = m_Tris[3*i+1];
		unsigned int idx3 = m_Tris[3*i+2];

		printf("%u [%u %u %u]: [%f %f %f] - [%f %f %f] - [%f %f %f]\n", i, 
			idx1, idx2, idx3,
			m_TriVerts[3*idx1+0], m_TriVerts[3*idx1+1], m_TriVerts[3*idx1+2],
			m_TriVerts[3*idx2+0], m_TriVerts[3*idx2+1], m_TriVerts[3*idx2+2],
			m_TriVerts[3*idx3+0], m_TriVerts[3*idx3+1], m_TriVerts[3*idx3+2] 
			);
	}
}

void Geometry::setPointSize( float pointSize )
{
	if( pointSize > 0 ) // do I need this condition ?
		m_PointSize = pointSize;
}

void Geometry::setPointColor( float r, float g, float b )
{
	m_UniquePointColors[0] = r;
	m_UniquePointColors[1] = g;
	m_UniquePointColors[2] = b;
}

void Geometry::setLineColor( float r, float g, float b )
{
	m_UniqueLineColors[0] = r;
	m_UniqueLineColors[1] = g;
	m_UniqueLineColors[2] = b;
}

void Geometry::SetLineWidth( float lineWidth )
{
	m_LineWidth = lineWidth;
}


void Geometry::setWireframeColor( float r, float g, float b )
{
	m_UniqueWireframeColors[0] = r;
	m_UniqueWireframeColors[1] = g;
	m_UniqueWireframeColors[2] = b;
}

bool Geometry::useUniquePointColors()
{
	return 	m_UsePointColors;
}

bool Geometry::useUniqueLineColors()
{
	return 	m_UseLineColors;
}


bool Geometry::useUniqueTriangleColors()
{
	return 	m_UseTriangleColors;
}

bool Geometry::useUniqueWireframeColors()
{
	return 	m_UseWireframeColors;
}

void Geometry::allocateTriangleDerivates()
{
	if( m_NumTriVerts <= 2 ) return;
	if( m_NumTris <= 0 ) return;

	if( m_TriMeanCurv ) { delete []m_TriMeanCurv; m_TriMeanCurv = 0; }
	m_TriMeanCurv = new float[ m_NumTriVerts ];
	if( m_TriGaussianCurv ) { delete []m_TriGaussianCurv; m_TriGaussianCurv = 0; }
	m_TriGaussianCurv = new float[ m_NumTriVerts ];
	if( m_K1 ) { delete []m_K1; m_K1 = 0; }
	m_K1 = new float[ m_NumTriVerts*3 ];
	if( m_K2 ) { delete []m_K2; m_K2 = 0; }
	m_K2 = new float[ m_NumTriVerts*3 ];
}

void Geometry::computeTriangleMeshDerivatives()
{
	if( m_TriangleDerivativesValid ) return;
	allocateTriangleDerivates();
	MeshDerivatives* meshDerivatives = new MeshDerivatives( m_NumTris, m_NumTriVerts, m_Tris, m_TriVerts, 
		m_TriVertNormals, m_TriMeanCurv, m_TriGaussianCurv, m_K1, m_K2);
	meshDerivatives->computeDerivatives();
	m_TriangleDerivativesValid = true;
	delete meshDerivatives;
}

void Geometry::computeDerivatives()
{
	computeTriangleMeshDerivatives();
}

bool Geometry::printTriangleDerivatives( FILE* fp )
{
	if( !fp ) return false;
	if( !m_TriangleDerivativesValid ) return false; // we wont bother trying to compute them here

	if( !m_TriVertNormals || !m_TriMeanCurv || !m_TriGaussianCurv || !m_K1 || !m_K2 ) return false;

	int i;
	// print normals
	fprintf(fp, "Normals:\n");
	for( i=0; i<m_NumTriVerts; i++ )
	{
		fprintf(fp, "%8.3f %8.3f %8.3f\n", m_TriVertNormals[3*i+0], m_TriVertNormals[3*i+1], m_TriVertNormals[3*i+2] );
	}

	// print mean curv and gauss curv
	fprintf(fp, "Mean, Gaussian curvatures:\n");
	for( i=0; i<m_NumTriVerts; i++ )
	{
		fprintf(fp, "%10.7f %10.7f\n", m_TriMeanCurv[i], m_TriGaussianCurv[i] );
	}

	// K1 and K2 are not currently computed SKVINAY
	return true;
}

bool Geometry::printDerivatives( const char* filename )
{
	if( !filename ) return false;
	FILE* fp = fopen( filename, "w" );
	if( !fp ) return false;

	bool ret = printTriangleDerivatives( fp );

	fclose( fp );
	return ret;
}

void Geometry::setTriangleColorsByNormals()
{
	computeDerivatives();
	if( !m_TriangleDerivativesValid ) return; // we wont bother trying to compute them here

	AllocateTriVertColors();
	int i;
	for( i=0; i<m_NumTriVerts*3; i++ )
	{
		m_TriVertColors[i] = m_TriVertNormals[i]/2.0 + 0.5;
	}
}

void Geometry::setTriangleColorsByMeanCurv()
{
	computeDerivatives();
	if( !m_TriangleDerivativesValid ) return; // we wont bother trying to compute them here

	AllocateTriVertColors();
	int i;
	for( i=0; i<m_NumTriVerts; i++ )
	{
		double h = m_TriMeanCurv[i];
		double r = 1.0, g = 1.0, b = 1.0;
		if( h < 0 ) 
		{
			r = 1.0 + h*5;
			b = 1.0 + h*5;
			if( h < -0.2 )
			{
				r = 0.0; b = 0.0; 
			}
		}
		if( h > 0 ) 
		{
			g = 1.0 - h*5;
			b = 1.0 - h*5;
			if( h > 0.2 )
			{
				g = 0.0; b = 0.0; 
			}
		}

		m_TriVertColors[3*i+0] = r;
		m_TriVertColors[3*i+1] = g;
		m_TriVertColors[3*i+2] = b;
	}
}

void Geometry::setTriangleColorsByGaussianCurv()
{
	computeDerivatives();
	if( !m_TriangleDerivativesValid ) return; // we wont bother trying to compute them here

	AllocateTriVertColors();
	int i;
	for( i=0; i<m_NumTriVerts; i++ )
	{
		double k = m_TriGaussianCurv[i];
		double r = 1.0, g = 1.0, b = 1.0;
		if( k < 0 ) 
		{
			double k2 = -1*sqrt(-k);
			r = 1.0 + k2*5;
			b = 1.0 + k2*5;
			if( k2 < -0.2 )
			{
				r = 0.0; b = 0.0; 
			}
		}
		if( k > 0 ) 
		{
			double k2 = sqrt(k);
			g = 1.0 - k2*5;
			b = 1.0 - k2*5;
			if( k2 > 0.2 )
			{
				g = 0.0; b = 0.0; 
			}
		}

		m_TriVertColors[3*i+0] = r;
		m_TriVertColors[3*i+1] = g;
		m_TriVertColors[3*i+2] = b;
	}
}

void Geometry::setTriangleColors( int renderingMode )
{
	switch( renderingMode )
	{
	case 0:
		setTriangleColorsByNormals();
		return;
	case 1:
		setTriangleColorsByMeanCurv();
		return;
	case 2:
		setTriangleColorsByGaussianCurv();
		return;
	}
}

void Geometry::setColors( int renderingMode )
{
	setTriangleColors( renderingMode );
}

bool Geometry::getPointInTriangle( float x1, float y1, float z1, 
									float x2, float y2, float z2, 
									float x3, float y3, float z3, 
									float *xinterp, float *yinterp, float *zinterp )
{
	float v1x = x2-x1;
	float v1y = y2-y1;
	float v1z = z2-z1;

	float v2x = x3-x1;
	float v2y = y3-y1;
	float v2z = z3-z1;

	float r1 = (float)rand()/RAND_MAX;
	float r2 = (float)rand()/RAND_MAX;

	// this is not uniformly random. SKVINAY
	(*xinterp) = r1*v1x + (1.0-r1)*r2*v2x + x1;
	(*yinterp) = r1*v1y + (1.0-r1)*r2*v2y + y1;
	(*zinterp) = r1*v1z + (1.0-r1)*r2*v2z + z1;

	return true;
}




























