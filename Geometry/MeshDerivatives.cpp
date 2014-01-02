// MeshDerivatives.cpp: implementation of the MeshDerivatives class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshDerivatives.h"
#include <math.h>

#include <cstdio>

#define PI 3.1415926535

FILE* fpTheta;
FILE* fpAMixed;

MeshDerivatives::MeshDerivatives(
	unsigned int numTris, 
	unsigned int numTriVerts, 
	unsigned int* tris, 
	float* triVerts, 
	float* triVertNormals, 
	float* triMeanCurv, 
	float* triGaussianCurv, 
	float* k1, 
	float* k2)
{
	m_NumTris = numTris;
	m_NumTriVerts = numTriVerts;
	m_Tris = tris;
	m_TriVerts = triVerts;
	m_TriVertNormals = triVertNormals;
	m_TriMeanCurv = triMeanCurv;
	m_TriGaussianCurv = triGaussianCurv; 
	m_K1 = k1;
	m_K2 = k2;

	if( m_NumTriVerts > 2  && m_NumTris > 0)
	{
		m_FaceSet = new vector<int>[m_NumTriVerts];
		m_Angles = new double[m_NumTris*3];
	}
	else
		m_FaceSet = 0;
}

MeshDerivatives::~MeshDerivatives()
{
	delete []m_FaceSet; m_FaceSet = 0;
	delete []m_Angles; m_Angles = 0;
}

bool MeshDerivatives::computeFaceSet()
{
	if( (m_NumTris<=0) || 
		!m_Tris ||
		!m_FaceSet )
		return false;

	int i;
	for( i=0; i<m_NumTris; i++ )
	{
		unsigned int idx1 = m_Tris[3*i+0];
		unsigned int idx2 = m_Tris[3*i+1];
		unsigned int idx3 = m_Tris[3*i+2];

		m_FaceSet[idx1].push_back(i);
		m_FaceSet[idx2].push_back(i);
		m_FaceSet[idx3].push_back(i);
	}

	return true;
}

// be consistant with ordering! 
bool MeshDerivatives::findVerts( unsigned int triIndex, unsigned int vertIndex, unsigned int* vleft, unsigned int* vright )
{
	if( !vleft || !vright ) return false;
	if( triIndex > m_NumTris ) return false;

	if( m_Tris[3*triIndex+0] == vertIndex )
	{
		*vleft  = m_Tris[3*triIndex+2];
		*vright = m_Tris[3*triIndex+1];
		return true;
	}
	if( m_Tris[3*triIndex+1] == vertIndex )
	{
		*vleft  = m_Tris[3*triIndex+0];
		*vright = m_Tris[3*triIndex+2];
		return true;
	}
	if( m_Tris[3*triIndex+2] == vertIndex )
	{
		*vleft  = m_Tris[3*triIndex+1];
		*vright = m_Tris[3*triIndex+0];
		return true;
	}

	return true;
}

double MeshDerivatives::getEdgeLengthSQ( unsigned int v1, unsigned int v2)
{
	double x1 = m_TriVerts[3*v1+0];
	double y1 = m_TriVerts[3*v1+1];
	double z1 = m_TriVerts[3*v1+2];

	double x2 = m_TriVerts[3*v2+0];
	double y2 = m_TriVerts[3*v2+1];
	double z2 = m_TriVerts[3*v2+2];

	return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

double MeshDerivatives::voronoiArea(unsigned int vertIndex, unsigned int vleft, unsigned int vright, 
									double angleVLeft, double angleVRight )
{
	double PR_sq = getEdgeLengthSQ(vertIndex, vleft);
	double PQ_sq = getEdgeLengthSQ(vertIndex, vright);
	double cotQ = 1.0/tan(angleVLeft);
	double cotR = 1.0/tan(angleVRight);

	return 1.0/8.0*(PR_sq*cotQ + PQ_sq*cotR);
}

bool MeshDerivatives::isObtuseTriangle( double angle1, double angle2, double angle3 )
{
	if( angle1 < 0 || angle1 > PI / 2.0 ) return true;
	if( angle2 < 0 || angle2 > PI / 2.0 ) return true;
	if( angle3 < 0 || angle3 > PI / 2.0 ) return true;

	return false;
}

double MeshDerivatives::getTriangleArea( unsigned int v1, unsigned int v2, unsigned int v3 )
{
	double x1 = m_TriVerts[3*v1+0];
	double y1 = m_TriVerts[3*v1+1];
	double z1 = m_TriVerts[3*v1+2];

	double x2 = m_TriVerts[3*v2+0];
	double y2 = m_TriVerts[3*v2+1];
	double z2 = m_TriVerts[3*v2+2];

	double x3 = m_TriVerts[3*v3+0];
	double y3 = m_TriVerts[3*v3+1];
	double z3 = m_TriVerts[3*v3+2];

	double a = getEdgeLength( v1, v2 );
	double b = getEdgeLength( v1, v3 );
	double c = getEdgeLength( v2, v3 );

	double s = (a+b+c)/2.0;

	return sqrt(s*(s-a)*(s-b)*(s-c));
}

double MeshDerivatives::areaOfInfinitesimalTriangle( int triangleIdx, unsigned int vertIndex, unsigned int vleft, unsigned int vright )
{
	double angleVertIndex = 0;
	double angleVLeft = 0;
	double angleVRight = 0;

	if( m_Tris[3*triangleIdx+0] == vertIndex )
		angleVertIndex = m_Angles[triangleIdx*3+0];
	if( m_Tris[3*triangleIdx+0] == vleft )
		angleVLeft = m_Angles[triangleIdx*3+0];
	if( m_Tris[3*triangleIdx+0] == vright )
		angleVRight = m_Angles[triangleIdx*3+0];

	if( m_Tris[3*triangleIdx+1] == vertIndex )
		angleVertIndex = m_Angles[triangleIdx*3+1];
	if( m_Tris[3*triangleIdx+1] == vleft )
		angleVLeft = m_Angles[triangleIdx*3+1];
	if( m_Tris[3*triangleIdx+1] == vright )
		angleVRight = m_Angles[triangleIdx*3+1];

	if( m_Tris[3*triangleIdx+2] == vertIndex )
		angleVertIndex = m_Angles[triangleIdx*3+2];
	if( m_Tris[3*triangleIdx+2] == vleft )
		angleVLeft = m_Angles[triangleIdx*3+2];
	if( m_Tris[3*triangleIdx+2] == vright )
		angleVRight = m_Angles[triangleIdx*3+2];

	double f = 180.0/PI;
	if( isObtuseTriangle( angleVertIndex, angleVLeft, angleVRight ) )
	{
		double triArea = getTriangleArea( vertIndex, vleft, vright );
		if( angleVertIndex >= 0 && angleVertIndex <= PI / 2.0 )
		{
			fprintf( fpAMixed, "%10.7lf T/4: %10.7lf <%10.7lf, %10.7lf, %10.7lf = %10.7lf>\n", triArea, triArea / 4.0,
				angleVertIndex*f, angleVLeft*f, angleVRight*f, angleVertIndex*f + angleVLeft*f + angleVRight*f );
			return triArea / 4.0;
		}
		else
		{
			fprintf( fpAMixed, "%10.7lf T/2: %10.7lf <%10.7lf, %10.7lf, %10.7lf = %10.7lf>\n", triArea, triArea / 2.0,
				angleVertIndex*f, angleVLeft*f, angleVRight*f, angleVertIndex*f + angleVLeft*f + angleVRight*f );
			return triArea / 2.0;
		}
	}
	double triArea = getTriangleArea( vertIndex, vleft, vright );
	double vArea = voronoiArea(vertIndex, vleft, vright, angleVLeft, angleVRight);
	fprintf( fpAMixed, "%10.7lf v : %10.7lf <%10.7lf, %10.7lf, %10.7lf = %10.7lf>\n", triArea, vArea ,
				angleVertIndex*f, angleVLeft*f, angleVRight*f, angleVertIndex*f + angleVLeft*f + angleVRight*f );

	return vArea;
}

double MeshDerivatives::computeArea( unsigned int vertIndex )
{
	int numtris = m_FaceSet[vertIndex].size();
	if( numtris < 1 ) return 0; //wierd!

	double area = 0;

	int i;
	for( i=0; i<numtris; i++ )
	{
		unsigned int vleft, vright;
		if( !findVerts( m_FaceSet[vertIndex][i], vertIndex, &vleft, &vright ) ) 
		{
			continue;
		}
		area += areaOfInfinitesimalTriangle( m_FaceSet[vertIndex][i], vertIndex, vleft, vright );
	}

	return area;
}

double MeshDerivatives::computeThetaSum( unsigned int vertIndex )
{
	if( m_FaceSet[vertIndex].size() < 1 ) return 0; //wierd!

	//loop through all faces. Find which of 3 vertices this is. Pick up that angle and add to sum

	double thetaSum = 0;
	int i;
	for( i=0; i<m_FaceSet[vertIndex].size(); i++ )
	{
		unsigned int triangleIdx = m_FaceSet[vertIndex][i];
		if( m_Tris[3*triangleIdx+0] == vertIndex )
			thetaSum += m_Angles[triangleIdx*3+0];
		if( m_Tris[3*triangleIdx+1] == vertIndex )
			thetaSum += m_Angles[triangleIdx*3+1];
		if( m_Tris[3*triangleIdx+2] == vertIndex )
			thetaSum += m_Angles[triangleIdx*3+2];
	}
	return thetaSum;
}

void MeshDerivatives::unionVertex( unsigned int v, vector<unsigned int> *neighboringVertices )
{
	int j;
	for( j=0; j<neighboringVertices->size(); j++ )
	{
		if( neighboringVertices->at(j) == v ) return;
	}
	neighboringVertices->push_back( v );
}

void MeshDerivatives::getNeighbors( unsigned int vertIndex, vector<unsigned int> *neighboringVertices )
{
	int i;
	for( i=0; i<m_FaceSet[vertIndex].size(); i++ )
	{
		unsigned int v1 = m_Tris[m_FaceSet[vertIndex][i]*3+0];
		unsigned int v2 = m_Tris[m_FaceSet[vertIndex][i]*3+1];
		unsigned int v3 = m_Tris[m_FaceSet[vertIndex][i]*3+2];
		if( vertIndex != v1 )
			unionVertex( v1, neighboringVertices );
		if( vertIndex != v2 )
			unionVertex( v2, neighboringVertices );
		if( vertIndex != v3 )
			unionVertex( v3, neighboringVertices );
	}
}

void MeshDerivatives::getNeighborsAngles( unsigned int x_i, unsigned int x_j, double* alpha, double* beta, int* count )
{
	int i;
	for( i=0; i<m_FaceSet[x_i].size(); i++ )
	{
		unsigned int triIdx = m_FaceSet[x_i][i];

		unsigned int v1 = m_Tris[triIdx*3+0];
		unsigned int v2 = m_Tris[triIdx*3+1];
		unsigned int v3 = m_Tris[triIdx*3+2];

		if( (v1 == x_i || v1 == x_j) && (v2 == x_i || v2 == x_j) )
		{
			if( *count == 0 ) 
			{
				*alpha = m_Angles[triIdx*3+2];
				(*count)++;
			}
			else
			{
				*beta = m_Angles[triIdx*3+2];
				(*count)++;
				return;
			}
		}

		if( (v1 == x_i || v1 == x_j) && (v3 == x_i || v3 == x_j) )
		{
			if( *count == 0 ) 
			{
				*alpha = m_Angles[triIdx*3+1];
				(*count)++;
			}
			else
			{
				*beta = m_Angles[triIdx*3+1];
				(*count)++;
				return;
			}
		}

		if( (v2 == x_i || v2 == x_j) && (v3 == x_i || v3 == x_j) )
		{
			if( *count == 0 ) 
			{
				*alpha = m_Angles[triIdx*3+0];
				(*count)++;
			}
			else
			{
				*beta = m_Angles[triIdx*3+0];
				(*count)++;
				return;
			}
		}
	}
}

bool MeshDerivatives::getLaplaceBeltramiOperator(double* Kx, double* Ky, double* Kz, unsigned int vertIndex, double A_mixed)
{
	if( !Kx || !Ky || !Kz ) return false;

	vector<unsigned int> neighboringVertices; 

	getNeighbors( vertIndex, &neighboringVertices );

	if( neighboringVertices.size() < 2 ) return false;

	int i;
	for( i=0; i<neighboringVertices.size(); i++ )
	{
		// find the other 2 vertices and the respective angles. Add to K
		// if we can only get one other neighbor, do some funky approximation SKVINAY
		double alpha, beta;
		int count = 0;
		getNeighborsAngles( vertIndex, neighboringVertices[i], &alpha, &beta, &count );
		if( count == 0 ) 
		{
			continue;
		}
		double x1 = m_TriVerts[vertIndex*3+0];
		double y1 = m_TriVerts[vertIndex*3+1];
		double z1 = m_TriVerts[vertIndex*3+2];

		double x2 = m_TriVerts[neighboringVertices[i]*3+0];
		double y2 = m_TriVerts[neighboringVertices[i]*3+1];
		double z2 = m_TriVerts[neighboringVertices[i]*3+2];

		if( count == 1 ) 
		{
			(*Kx) += (1.0/tan(alpha))*(x2-x1);
			(*Ky) += (1.0/tan(alpha))*(y2-y1);
			(*Kz) += (1.0/tan(alpha))*(z2-z1);
		}
		else // better be 2!
		{
			(*Kx) += (1.0/tan(alpha) + 1.0/tan(beta))*(x2-x1);
			(*Ky) += (1.0/tan(alpha) + 1.0/tan(beta))*(y2-y1);
			(*Kz) += (1.0/tan(alpha) + 1.0/tan(beta))*(z2-z1);
		}
	}

	*Kx = (*Kx) / (2.0*A_mixed);
	*Ky = (*Ky) / (2.0*A_mixed);
	*Kz = (*Kz) / (2.0*A_mixed);
	return true;
}

bool MeshDerivatives::computeDerivatives( unsigned int vertIndex )
{
	if( m_FaceSet[vertIndex].size() < 1 ) return false; //wierd!

	double A_mixed = computeArea( vertIndex );
	if( A_mixed == 0 ) return false;

	double theta_sum = computeThetaSum( vertIndex );

	fprintf( fpTheta, "%8.3lf\n", theta_sum*180.0/PI );
	fprintf( fpAMixed, "%10.7lf\n", A_mixed );
	if( theta_sum < 2*PI ) 
	{
		printf("This is messed if input is a sphere\n");
	}
	// compute Gaussian curvature
	m_TriGaussianCurv[vertIndex] = (2.0*PI - theta_sum) / A_mixed;

	// for each vertex neighbor, compute alphas, betas, xj-xi, and update K(x) - laplace beltrami operator
	double Kx = 0, Ky = 0, Kz = 0;
	if( !getLaplaceBeltramiOperator(&Kx, &Ky, &Kz, vertIndex, A_mixed) ) return false;

	double norm = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
	m_TriMeanCurv[vertIndex] = 0.5*(norm);

	if( norm )
	{
		m_TriVertNormals[vertIndex*3 + 0] = Kx / norm;
		m_TriVertNormals[vertIndex*3 + 1] = Ky / norm;
		m_TriVertNormals[vertIndex*3 + 2] = Kz / norm;
	}
	else
	{
		m_TriVertNormals[vertIndex*3 + 0] = 1;
		m_TriVertNormals[vertIndex*3 + 1] = 0;
		m_TriVertNormals[vertIndex*3 + 2] = 0;
	}
	// compute normal vector, mean curvature

	return true;
}

double MeshDerivatives::getEdgeLength( unsigned int v1, unsigned int v2 )
{
	return sqrt(getEdgeLengthSQ(v1, v2));
}

bool MeshDerivatives::computeAngles()
{
	int i;
	for( i=0; i<m_NumTris; i++ )
	{
		unsigned int A = m_Tris[3*i+0];
		unsigned int B = m_Tris[3*i+1];
		unsigned int C = m_Tris[3*i+2];

		double a = getEdgeLength( B, C );
		double b = getEdgeLength( A, C );
		double c = getEdgeLength( A, B );
		m_Angles[i*3+0] = acos((b*b + c*c - a*a) / (2.0*b*c));
		m_Angles[i*3+1] = acos((a*a + c*c - b*b) / (2.0*a*c));
		m_Angles[i*3+2] = acos((a*a + b*b - c*c) / (2.0*a*b));
	}

/*	{
		FILE* fp = fopen("angles.txt", "w" );
		int i;
		for( i=0; i<m_NumTris; i++ )
		{
			fprintf( fp, "%8.3lf %8.3lf% 8.3lf\n", m_Angles[i*3+0], m_Angles[i*3+1], m_Angles[i*3+2] );
		}
		fclose(fp);
	}*/
	return true;
}

bool MeshDerivatives::computeDerivatives()
{
	fpTheta = fopen( "thetasum.txt","w");
	fpAMixed = fopen("amixed.txt","w");

	if( (m_NumTris<=0) || 
		(m_NumTriVerts<=2) || 
		!m_Tris ||
		!m_TriVerts ||
		!m_TriVertNormals ||
		!m_TriMeanCurv ||
		!m_TriGaussianCurv ||
		!m_K1 ||
		!m_K2 ||
		!m_FaceSet ||
		!m_Angles )
		return false;

	// compute face set for each vertex
	if( !computeFaceSet() ) return false;
	
	if( !computeAngles() ) return false;

	// compute n, H, G, k1, k2 for each vertex. (right now im not implementing the minimization for k1, k2 - SKVINAY)
	int i;
	for( i=0; i<m_NumTriVerts; i++ )
	{
		computeDerivatives( i );
	}

	fclose( fpTheta );
	fclose( fpAMixed );
	return true;
}
