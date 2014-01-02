/*****************************************************************************/
/*                             ______________________                        */
/*                            / _ _ _ _ _ _ _ _ _ _ _)                       */
/*            ____  ____  _  / /__  __  _____  __                            */
/*           (_  _)( ___)( \/ /(  \/  )(  _  )(  )                           */
/*             )(   )__)  )  (  )    (  )(_)(  )(__                          */
/*            (__) (____)/ /\_)(_/\/\_)(_____)(____)                         */
/*            _ _ _ _ __/ /                                                  */
/*           (___________/                     ___  ___                      */
/*                                      \  )| |   ) _ _|\   )                */
/*                                 ---   \/ | |  / |___| \_/                 */
/*                                                       _/                  */
/*                                                                           */
/*   Copyright (C) The University of Texas at Austin                         */
/*                                                                           */
/*     Authors:    Vinay Siddavanahalli <skvinay@cs.utexas.edu>   2004-2005  */
/*     Authors:     Anthony Thane        <thanea@ices.utexas.edu> 2003-2003  */
/*                                                                           */
/*     Principal Investigator: Chandrajit Bajaj <bajaj@ices.utexas.edu>      */
/*                                                                           */
/*         Professor of Computer Sciences,                                   */
/*         Computational and Applied Mathematics Chair in Visualization,     */
/*         Director, Computational Visualization Center (CVC),               */
/*         Institute of Computational Engineering and Sciences (ICES)        */
/*         The University of Texas at Austin,                                */
/*         201 East 24th Street, ACES 2.324A,                                */
/*         1 University Station, C0200                                       */
/*         Austin, TX 78712-0027                                             */
/*         http://www.cs.utexas.edu/~bajaj                                   */
/*                                                                           */
/*         http://www.ices.utexas.edu/CVC                                    */
/*  This software comes with a license. Using this code implies that you     */
/*  read, understood and agreed to all the terms and conditions in that      */
/*  license.                                                                 */
/*                                                                           */
/*  We request that you agree to acknowledge the use of the software that    */
/*  results in any published work, including scientific papers, films and    */
/*  videotapes by citing the reference listed below                          */
/*                                                                           */
/*    C. Bajaj, P. Djeu, V. Siddavanahalli, A. Thane,                        */
/*    Interactive Visual Exploration of Large Flexible Multi-component       */
/*    Molecular Complexes,                                                   */
/*    Proc. of the Annual IEEE Visualization Conference, October 2004,       */
/*    Austin, Texas, IEEE Computer Society Press, pp. 243-250.               */
/*                                                                           */
/*****************************************************************************/

// Ray.cpp: implementation of the Ray class.
//
//////////////////////////////////////////////////////////////////////

#include "Ray.h"

#include <math.h>

using CCVOpenGLMath::Vector; 
using CCVOpenGLMath::Ray; 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Ray::Ray() : m_Origin(0.0f, 0.0f, 0.0f, 1.0f), m_Dir(0.0f, 0.0f, 1.0f, 0.0f)
{
}

Ray::Ray(const Vector& origin, const Vector& dir) : m_Origin(origin), m_Dir(dir)
{
	
}

Ray::~Ray()
{

}

Vector Ray::getPointOnRay(float t) const
{
	return m_Origin+m_Dir*t;
}

float Ray::nearestTOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToXAxis(Origin);
	float t = -(ray.m_Origin[1]*ray.m_Dir[1] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[2]*ray.m_Dir[2]) );
	return t;
}

float Ray::nearestTOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToYAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[0]*ray.m_Dir[0]+ray.m_Dir[2]*ray.m_Dir[2]));
	return t;
}

float Ray::nearestTOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToZAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[1]*ray.m_Dir[1])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[0]*ray.m_Dir[0]));
	return t;
}

Vector Ray::nearestPointOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnXAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[1] = Origin[1];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnYAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnZAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[1] = Origin[1];
	return result;
}

float Ray::distanceToXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

float Ray::distanceToYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[0]-ray.m_Origin[0]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[0]*ray.m_Dir[0] )
		);
}

float Ray::distanceToZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[0]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[0] ) /
		(float)sqrt( ray.m_Dir[0]*ray.m_Dir[0] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

/*********************************************************/
/*                                                       */
/*  Returns false if there is no intersection.           */
/*  Else, it returns both the points and the values of   */
/*  the parameter where the intersections took place.    */
/*                                                       */
/*********************************************************/
bool Ray::intersectSphere( Vector center, float radius, Vector *point1, Vector* point2, float *distance1, float* distance2 )
{
	if( !point1 || !point2 ) return false;
	if( radius <= 0 ) return false;


	/// solve quadratic equation /////
	////	A = Xd^2 + Yd^2 + Zd^2
	////	B = 2 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc))
	////	C = (X0 - Xc)^2 + (Y0 - Yc)^2 + (Z0 - Zc)^2 - Sr^2
	///////////////////////////////////


	float A =	m_Dir[0]*m_Dir[0] + 
				m_Dir[1]*m_Dir[1] + 
				m_Dir[2]*m_Dir[2];
	  
	float B = 2* (	m_Dir[0] * (m_Origin[0] - center[0]) +
					m_Dir[1] * (m_Origin[1] - center[1]) +
					m_Dir[2] * (m_Origin[2] - center[2]) );

	float C = (m_Origin[0] - center[0])*(m_Origin[0] - center[0]) +
			  (m_Origin[1] - center[1])*(m_Origin[1] - center[1]) +
			  (m_Origin[2] - center[2])*(m_Origin[2] - center[2]) -
			  radius*radius;

	float discriminant = B*B - 4*A*C;
	if( discriminant < 0 ) return false;

	*distance1 = (float)(( -B - sqrt( discriminant ) ) / ( 4.0 * A * C ));
	*distance2 = (float)(( -B + sqrt( discriminant ) ) / ( 4.0 * A * C ));

	*point1 = m_Origin + m_Dir * (*distance1);
	*point2 = m_Origin + m_Dir * (*distance2);

	return true;
}