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

// Quaternion.cpp: implementation of the Quaternion class.
//
//////////////////////////////////////////////////////////////////////

#include "Quaternion.h"
#include "Vector.h"
#include "Ray.h"
#include "Matrix.h"

#include <math.h>

using CCVOpenGLMath::Tuple;
using CCVOpenGLMath::Vector;
using CCVOpenGLMath::Ray;
using CCVOpenGLMath::Quaternion;
using CCVOpenGLMath::Matrix;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Quaternion::Quaternion(float w, float x, float y, float z) : Tuple(w,x,y,z)
{
}

Quaternion::Quaternion() : Tuple(1.0, 0.0, 0.0, 0.0)
{

}

Quaternion::Quaternion(const Vector& vec)
{
	p[0] = 0.0f;
	p[1] = vec[0];
	p[2] = vec[1];
	p[3] = vec[2];
}

Quaternion::~Quaternion()
{

}

Quaternion::Quaternion(const Quaternion& copy): Tuple(copy)
{

}

Quaternion& Quaternion::operator=(const Quaternion& copy)
{
	if (this!=&copy) {
		set(copy);
	}
	return *this;
}

Quaternion& Quaternion::set(float w, float x, float y, float z)
{
	Tuple::set(w,x,y,z);
	return *this;
}

Quaternion& Quaternion::set(float* array)
{
	Tuple::set(array);
	return *this;
}

Quaternion& Quaternion::set(const Quaternion& copy)
{
	Tuple::set(copy);
	return *this;
}

Quaternion Quaternion::operator*(const Quaternion& quat) const
{
	return Quaternion(
		p[0]*quat[0] - p[1]*quat[1] - p[2]*quat[2] - p[3]*quat[3],
		p[0]*quat[1] + p[1]*quat[0] + p[2]*quat[3] - p[3]*quat[2],
		p[0]*quat[2] - p[1]*quat[3] + p[2]*quat[0] + p[3]*quat[1],
		p[0]*quat[3] + p[1]*quat[2] - p[2]*quat[1] + p[3]*quat[0]
		);
}

Quaternion& Quaternion::preMultiply(const Quaternion& quat)
{
	return set(quat*(*this));
}

Quaternion& Quaternion::postMultiply(const Quaternion& quat)
{
	return set((*this)*quat);
}

Quaternion& Quaternion::rotate(float angle, float x, float y, float z)
{
	return preMultiply(rotation(angle, x, y, z));
}

Quaternion Quaternion::operator*(float scalar) const
{
	return Quaternion(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]*scalar);
}

Quaternion& Quaternion::operator*=(float scalar)
{
	return set(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]*scalar);
}

Quaternion Quaternion::operator/(float scalar) const
{
	return Quaternion(p[0]/scalar, p[1]/scalar, p[2]/scalar, p[3]/scalar);
}

Quaternion& Quaternion::operator/=(float scalar)
{
	return set(p[0]/scalar, p[1]/scalar, p[2]/scalar, p[3]/scalar);
}

Quaternion& Quaternion::normalize()
{
	float n = norm();
	return (*this)/=n;
}

Quaternion Quaternion::conjugate() const
{
	return Quaternion(p[0], -p[1], -p[2], -p[3]);
}

Quaternion Quaternion::inverse() const
{
	return conjugate()/norm();
}

float Quaternion::norm() const
{
	return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
}

Vector Quaternion::applyRotation(const Vector& vec) const
{
	Quaternion result = (*this) * Quaternion(vec) * (conjugate());
	return Vector(result[1], result[2], result[3], vec[3]);
}

Ray Quaternion::applyRotation(const Ray& ray) const
{
	Quaternion origin = (*this) * Quaternion(ray.m_Origin) * (conjugate());
	Quaternion dir = (*this) * Quaternion(ray.m_Dir) * (conjugate());
	return Ray(Vector(origin[1], origin[2], origin[3], ray.m_Origin[3]),
		Vector(dir[1], dir[2], dir[3], ray.m_Dir[3]));
}

Matrix Quaternion::buildMatrix() const
{
	float w = p[0];
	float x = p[1];
	float y = p[2];
	float z = p[3];
	return Matrix(
		1.0f-2.0f*y*y-2.0f*z*z,	2.0f*x*y-2.0f*w*z,		2.0f*x*z + 2.0f*w*y, 0.0f,
		2.0f*x*y + 2.0f*w*z,	1.0f - 2.0f*x*x - 2.0f*z*z,	2.0f*y*z - 2.0f*w*x, 0.0f,
        2.0f*x*z - 2.0f*w*y,	2.0f*y*z + 2.0f*w*x,		1.0f - 2.0f*x*x - 2.0f*y*y, 0.0f,
		0.0f,0.0f,0.0f,1.0f
		);
}

Quaternion Quaternion::rotation(float angle, float x, float y, float z)
{
	float len = (float)sqrt(x*x+y*y+z*z);
	if (len!=0.0) {
		len = (float)(sin(angle/2.0f)/len);
		return Quaternion((float) cos(angle/2.0f), x*len, y*len, z*len);
	}
	else {
		return Quaternion();
	}
}

Quaternion Quaternion::rotation(float angle, const Vector& axis)
{
	float len = (float)sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
	if (len!=0.0) {
		len = (float)(sin(angle/2.0f)/len);
		return Quaternion((float) cos(angle/2.0f), axis[0]*len, axis[1]*len, axis[2]*len);
	}
	else {
		return Quaternion();
	}
}

Quaternion Quaternion::power( double scalar )
{
	float Dest[4];
	double theta;
	if (p[0]>=0.9999f) {
		theta = 0;
	}
	else if (p[0]<=-0.9999f) {
		theta = 2.0*3.1415926535897932384626433832795;
	}
	else {
		theta = acos(p[0]);
	}
	double u[3];
	double scale = p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	scale = sqrt(scale);
	if (p[1]==0.0f && p[2]==0.0f && p[3]==0.0f) {
		u[0] = 0.0;
		u[1] = 0.0;
		u[2] = 0.0;
	}
	else {
		u[0] = p[1]/scale;
		u[1] = p[2]/scale;
		u[2] = p[3]/scale;
	}
	Dest[0] = (float)cos(scalar*theta);
	Dest[1] = (float)(u[0] * sin(scalar*theta));
	Dest[2] = (float)(u[1] * sin(scalar*theta));
	Dest[3] = (float)(u[2] * sin(scalar*theta));
	return Quaternion( Dest[0], Dest[1], Dest[2], Dest[3] );
}

