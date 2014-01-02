/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
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

// RawFile.cpp: implementation of the RawFile class.
//
//////////////////////////////////////////////////////////////////////

#include "RawFile.h"
#include <Geometry.h>
#include <stdio.h>
//#include <qfileinfo.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>
#include <cstring>
#include <cstdlib>

RawFile RawFile::ms_RawFileRepresentative;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RawFile::RawFile()
{

}

RawFile::~RawFile()
{

}

Geometry* RawFile::loadFile(const string& fileName)
{
	bool zeroSeen = false, maxSeen = false;
	//QFileInfo fileInfo(fileName);

	//if (!fileInfo.exists()) {
	//	qDebug("File does not exist");
	//	return 0;
	//}
	struct stat fst;
	if (stat(fileName.c_str(), &fst) == -1) {
		printf("File does not exist\n");
		return 0;
	}

	//QString absFilePath = fileInfo.absFilePath();

	//FILE* fp = fopen(absFilePath, "r");
	FILE* fp = fopen(fileName.c_str(), "r");
	if (!fp) {
		//qDebug("Error opening file");
		printf("Error opening file\n");
		return 0;
	}

	// get the number of verts and polygons
	int numverts, numpolygons;
	if (2!=fscanf(fp, "%d %d", &numverts, &numpolygons)) {
		printf("Error reading in number of verts and polygons\n");
		fclose(fp);
		return 0;
	}
	// make sure the number of verts & polygons are positive
	if (numverts<0 || numpolygons<0) {
		printf("Negative number of verts or polygons\n");
		fclose(fp);
		return 0;
	}

	std::vector<double> vertPositions;
	int c;
	// read in the verts
	for (c=0; c<numverts; c++) 
	{
		float x, y, z;
		// read in a single vert
		if (3!=fscanf(fp, "%f %f %f\n", &x, &y, &z) )
		{
			printf("Error reading in vert # %d\n", c);
			fclose(fp);
			return 0;
		}
		vertPositions.push_back( x );
		vertPositions.push_back( y );
		vertPositions.push_back( z );
	}

	std::vector<unsigned int> pointIndices;
	std::vector<unsigned int> lineIndices;
	std::vector<unsigned int> triangleIndices;

	// read in the polygons
	for (c=0; c<numpolygons; c++) 
	{
		char polyIndices[1024];
		fgets( polyIndices, 1023, fp );
		std::vector<unsigned int> indexList;
		{
			char * pch;
			pch = strtok (polyIndices," \t\n");
			while (pch != NULL)
			{
				unsigned int idx = atoi(pch);
				if( idx == 0 ) zeroSeen = true;
				if( idx == numverts ) maxSeen = true;
				if( zeroSeen && maxSeen )
				{
					printf("Indices were both 0 and maxNumVerts!\n");
					fclose( fp );
					return 0;
				}
				if( idx > numverts )
				{
					printf("Index out of bounds\n");
					fclose( fp );
					return 0;
				}
				indexList.push_back( idx );
				pch = strtok (NULL, " \t\n");
			}
		}
		if( indexList.size() == 0 )
		{
		}
		else if( indexList.size() == 1 )
		{
			pointIndices.push_back( indexList[0] );
		}
		else if( indexList.size() == 2 )
		{
			lineIndices.push_back( indexList[0] );
			lineIndices.push_back( indexList[1] );
		}
		else
		{
			int n = indexList.size();
			int i;
			unsigned int firstIdx = indexList[0];
			unsigned int prevIdx = indexList[1];
			for( i=2; i<n; i++ )
			{
				triangleIndices.push_back( firstIdx );
				triangleIndices.push_back( prevIdx );
				triangleIndices.push_back( indexList[i] );
				prevIdx = indexList[i];
			}
		}
		
	}

	if( maxSeen )
	{
		int i;
		for( i=0; i<pointIndices.size();    i++ )    pointIndices[i]--;
		for( i=0; i<lineIndices.size();     i++ )     lineIndices[i]--;
		for( i=0; i<triangleIndices.size(); i++ ) triangleIndices[i]--;
	}

	if( (vertPositions.size()==0) || ((pointIndices.size() == 0) && (lineIndices.size() == 0) && (triangleIndices.size() == 0)) )
	{
		return 0;
	}

	// initialize the geometry
	Geometry* geometry = new Geometry;
	if( pointIndices.size() )    geometry->AllocatePoints(pointIndices.size());
	if( lineIndices.size() )     geometry->AllocateLines(numverts, lineIndices.size() / 2);
	if( triangleIndices.size() ) geometry->AllocateTris(numverts, triangleIndices.size() / 3);

	// copy vertices, unfortunately, need to duplicate
	if( pointIndices.size() )
	{
		int i;
		for( i=0; i<pointIndices.size(); i++ )
		{
			unsigned int pIdx = pointIndices[i];
			geometry->m_Points[3*i+0] = vertPositions[3*pIdx+0];
			geometry->m_Points[3*i+1] = vertPositions[3*pIdx+1];
			geometry->m_Points[3*i+2] = vertPositions[3*pIdx+2];
		}
	}
	
	if( lineIndices.size() )
	{
		int i;
		for( i=0; i<vertPositions.size(); i++ ) geometry->m_LineVerts[i] = vertPositions[i];
		for( i=0; i<lineIndices.size(); i++ ) geometry->m_Lines[i] = lineIndices[i];
	}

	if( triangleIndices.size() )
	{
		int i;
		for( i=0; i<vertPositions.size(); i++ ) geometry->m_TriVerts[i] = vertPositions[i];
		for( i=0; i<triangleIndices.size(); i++ ) geometry->m_Tris[i] = triangleIndices[i];
	}
	
	fclose(fp);

	//geometry->print();
	return geometry;
}

bool RawFile::checkType(const string& fileName)
{
	return false;
}

bool RawFile::saveFile(const Geometry* geometry, const string& fileName)
{

	FILE* fp;
	// open the file
	fp = fopen(fileName.c_str(), "w");
	if (!fp) {
		printf("Error opening the output file\n");
		return false;
	}

	// write the number of verts & tris
	if (0>=fprintf(fp, "%d %d\n", geometry->m_NumTriVerts, geometry->m_NumTris)) {
		//qDebug("Error writing the number of verts and tris");
		printf("Error writing the number of verts and tris\n");
		return false;
	}

	unsigned int c;
	// write out the verts
	for (c=0; c<geometry->m_NumTriVerts; c++) {
		if (0>=fprintf(fp, "%f %f %f\n", 
			(geometry->m_TriVerts[c*3+0]),
			(geometry->m_TriVerts[c*3+1]),
			(geometry->m_TriVerts[c*3+2]))) {
			//qDebug("Error writing out vert # %d", c);
			printf("Error writing out vert # %d\n", c);
			fclose(fp);
			return false;
		}
	}
	// write out the tris
	for (c=0; c<geometry->m_NumTris; c++) {
		if (0>=fprintf(fp, "%d %d %d\n", 
			(geometry->m_Tris[c*3+0]),
			(geometry->m_Tris[c*3+1]),
			(geometry->m_Tris[c*3+2]))) {
			//qDebug("Error writing out tri # %d", c);
			printf("Error writing out tri # %d\n", c);
			fclose(fp);
			return false;
		}
	}
	fclose(fp);
	return true;
}

GeometryFileType* RawFile::getRepresentative()
{
	return &ms_RawFileRepresentative;
}

