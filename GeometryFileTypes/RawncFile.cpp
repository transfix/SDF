// RawncFile.cpp: implementation of the RawncFile class.
//
//////////////////////////////////////////////////////////////////////

#include "RawncFile.h"
#include <Geometry.h>
#include <stdio.h>
//#include <qfileinfo.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>
#include <cstring>
#include <cstdlib>

RawncFile RawncFile::ms_RawncFileRepresentative;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RawncFile::RawncFile()
{

}

RawncFile::~RawncFile()
{

}

Geometry* RawncFile::loadFile(const string& fileName)
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
	std::vector<double> vertColors;
	std::vector<double> vertNormals;

	int c;
	// read in the verts
	for (c=0; c<numverts; c++) 
	{
		float x, y, z, nx, ny, nz, r, g, b;
		// read in a single vert
		if (9!=fscanf(fp, "%f %f %f %f %f %f %f %f %f\n", &x, &y, &z, &nx, &ny, &nz, &r, &g, &b) )
		{
			printf("Error reading in vert # %d\n", c);
			fclose(fp);
			return 0;
		}
		vertPositions.push_back( x );
		vertPositions.push_back( y );
		vertPositions.push_back( z );

		vertNormals.push_back( nx );
		vertNormals.push_back( ny );
		vertNormals.push_back( nz );

		vertColors.push_back( r );
		vertColors.push_back( g );
		vertColors.push_back( b );
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

	fclose(fp);
	
	// initialize the geometry
	Geometry* geometry = new Geometry;

	if( pointIndices.size() )
	{
		geometry->AllocatePoints(pointIndices.size());
		geometry->AllocatePointColors();
	}
	if( lineIndices.size() )
	{
		geometry->AllocateLines(numverts, lineIndices.size() / 2);
		geometry->AllocateLineColors();
	}
	if( triangleIndices.size() ) 
	{ 
		geometry->AllocateTris(numverts, triangleIndices.size() / 3);
		geometry->AllocateTriVertColors();
	}

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

			geometry->m_PointColorsArray[3*i+0] = vertColors[3*pIdx+0];
			geometry->m_PointColorsArray[3*i+1] = vertColors[3*pIdx+1];
			geometry->m_PointColorsArray[3*i+2] = vertColors[3*pIdx+2];
		}
	}
	
	if( lineIndices.size() )
	{
		int i;
		for( i=0; i<vertPositions.size(); i++ ) geometry->m_LineVerts[i] = vertPositions[i];
		for( i=0; i<vertColors.size(); i++ ) geometry->m_LineColors[i] = vertColors[i];
		for( i=0; i<lineIndices.size(); i++ ) geometry->m_Lines[i] = lineIndices[i];
	}

	if( triangleIndices.size() )
	{
		int i;
		for( i=0; i<vertPositions.size(); i++ ) geometry->m_TriVerts[i] = vertPositions[i];
		for( i=0; i<vertColors.size(); i++ ) geometry->m_TriVertColors[i] = vertColors[i];
		for( i=0; i<vertNormals.size(); i++ ) geometry->m_TriVertNormals[i] = vertNormals[i];
		for( i=0; i<triangleIndices.size(); i++ ) geometry->m_Tris[i] = triangleIndices[i];
	}

	geometry->SetTriNormalsReady();

	return geometry;
}

bool RawncFile::checkType(const string& fileName)
{
	return false;
}

bool RawncFile::saveMolCenters(const Geometry* geometry)
{
	FILE* fp1;
	FILE* fp2;
	// open the file
	fp1 = fopen("2bg9_trans.txt", "w");
	fp2 = fopen("1c2b_trans.txt", "w");
	if (!fp1 || !fp2) {
		printf("Error opening the output file\n");
		return false;
	}

	unsigned int c;
	// write out the tris
	for (c=0; c<geometry->m_NumTris; c++) 
	{
		unsigned int v1 = geometry->m_Tris[c*3+0];
		unsigned int v2 = geometry->m_Tris[c*3+1];
		unsigned int v3 = geometry->m_Tris[c*3+2];

		float x1 = geometry->m_TriVerts[v1*3+0];
		float y1 = geometry->m_TriVerts[v1*3+1];
		float z1 = geometry->m_TriVerts[v1*3+2];

		float x2 = geometry->m_TriVerts[v2*3+0];
		float y2 = geometry->m_TriVerts[v2*3+1];
		float z2 = geometry->m_TriVerts[v2*3+2];

		float x3 = geometry->m_TriVerts[v3*3+0];
		float y3 = geometry->m_TriVerts[v3*3+1];
		float z3 = geometry->m_TriVerts[v3*3+2];
		
		float r_val = (float)rand()/(RAND_MAX);
		if( r_val < 0.6 )
		{
			float xinterp = 0, yinterp = 0, zinterp = 0;
			Geometry::getPointInTriangle( x1, y1, z1, x2, y2, z2, x3, y3, z3, &xinterp, &yinterp, &zinterp );
			float nx, ny, nz;
			nx = geometry->m_TriFlatNormals[c*3+0];
			ny = geometry->m_TriFlatNormals[c*3+1];
			nz = geometry->m_TriFlatNormals[c*3+2];

			if( r_val < 0.4 )
			{
				fprintf( fp1, "%f %f %f %f %f %f\n", 
					xinterp, yinterp, zinterp, 
					nx, ny, nz );
			}
			else
			{
				float height = 400;
				xinterp += height*nx;
				yinterp += height*ny;
				zinterp += height*nz;

				fprintf( fp2, "%f %f %f %f %f %f\n", 
					xinterp, yinterp, zinterp, 
					nx, ny, nz );
			}
		}
	}
	fclose(fp1);
	fclose(fp2);

	return true;
}

bool RawncFile::saveFile(const Geometry* geometry, const string& fileName)
{
	saveMolCenters(geometry);

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
		if (geometry->m_TriVertColors) {
			if (0>=fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", 
				(geometry->m_TriVerts[c*3+0]),
				(geometry->m_TriVerts[c*3+1]),
				(geometry->m_TriVerts[c*3+2]),
				(geometry->m_TriVertNormals[c*3+0]),
				(geometry->m_TriVertNormals[c*3+1]),
				(geometry->m_TriVertNormals[c*3+2]),
				(geometry->m_TriVertColors[c*3+0]),
				(geometry->m_TriVertColors[c*3+1]),
				(geometry->m_TriVertColors[c*3+2]))) {
				//qDebug("Error writing out vert # %d", c);
				printf("Error writing out vert # %d\n", c);
				fclose(fp);
				return false;
			}
		}
		else {
			if (0>=fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", 
				(geometry->m_TriVerts[c*3+0]),
				(geometry->m_TriVerts[c*3+1]),
				(geometry->m_TriVerts[c*3+2]),
				(geometry->m_TriVertNormals[c*3+0]),
				(geometry->m_TriVertNormals[c*3+1]),
				(geometry->m_TriVertNormals[c*3+2]),
				(1.0),
				(1.0),
				(1.0))) {
				//qDebug("Error writing out vert # %d", c);
				printf("Error writing out vert # %d\n", c);
				fclose(fp);
				return false;
			}
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


GeometryFileType* RawncFile::getRepresentative()
{
	return &ms_RawncFileRepresentative;
}

