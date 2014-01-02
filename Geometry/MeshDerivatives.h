// MeshDerivatives.h: interface for the MeshDerivatives class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CCV_TEXMOL_GEOMETRY_MESH_DERIVATIVES_H
#define CCV_TEXMOL_GEOMETRY_MESH_DERIVATIVES_H

#include <vector>

using namespace std;
// this class is based on the theory in 
// "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds"
// Mark Meyer, Mathieu Desbrun, Peter Schroder, Alan H Barr.
class MeshDerivatives  
{
public:
	MeshDerivatives( 
		unsigned int numTris, 
		unsigned int numTriVerts, 
		unsigned int* tris, 
		float* triVerts, 
		float* m_TriVertNormals, 
		float* triMeanCurv, 
		float* triGaussianCurv, 
		float* k1, 
		float* k2);
	virtual ~MeshDerivatives();

	bool computeDerivatives();

protected:
	bool computeFaceSet();
	bool computeDerivatives( unsigned int vertIndex );
	double computeArea( unsigned int vertIndex );
	double computeThetaSum( unsigned int vertIndex );
	bool findVerts( unsigned int triIndex, unsigned int vertIndex, unsigned int* vleft, unsigned int* vright );
	double areaOfInfinitesimalTriangle( int triangleIdx, unsigned int vertIndex, unsigned int vleft, unsigned int vright );
	double getEdgeLength( unsigned int v1, unsigned int v2 );
	double voronoiArea(	unsigned int vertIndex, unsigned int vleft, unsigned int vright, 
						double angleVLeft, double angleVRight );
	double getEdgeLengthSQ( unsigned int v1, unsigned int v2);
	bool isObtuseTriangle( double angle1, double angle2, double angle3 );
	double getTriangleArea( unsigned int vertIndex, unsigned int vleft, unsigned int vright );
	bool computeAngles();
	void getNeighbors( unsigned int vertIndex, vector<unsigned int> *neighboringVertices );
	bool getLaplaceBeltramiOperator(double* Kx, double* Ky, double* Kz, unsigned int vertIndex, double A_mixed);
	void unionVertex( unsigned int v, vector<unsigned int> *neighboringVertices );
	void getNeighborsAngles( unsigned int x_i, unsigned int x_j, double* alpha, double* beta, int* count );

	unsigned int m_NumTris;
	unsigned int m_NumTriVerts;
	unsigned int* m_Tris;
	float* m_TriVerts;
	float* m_TriVertNormals;
	float* m_TriMeanCurv;
	float* m_TriGaussianCurv;
	float* m_K1;
	float* m_K2;

	vector<int>* m_FaceSet;
	double* m_Angles;
};

#endif 