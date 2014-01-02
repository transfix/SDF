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

// Vector.cpp: implementation of the Vector class.
//
//////////////////////////////////////////////////////////////////////

#include "Vector.h"
#include <math.h>

using CCVOpenGLMath::Tuple;
using CCVOpenGLMath::Vector;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const float EPS = 0.00001f;

Vector::Vector(float x, float y, float z, float w) : Tuple(x,y,z,w)
{
}

Vector::Vector() : Tuple()
{
}

Vector::Vector(float* array)
{
	set(array);
}

Vector::~Vector()
{

}

Vector::Vector(const Vector& copy): Tuple(copy)
{
}

Vector& Vector::operator=(const Vector& copy)
{
	if (this!=&copy) {
		set(copy);
	}
	return *this;
}


Vector& Vector::set(float x, float y, float z, float w)
{
	Tuple::set(x,y,z,w);
	return *this;
}

Vector& Vector::set(float* array)
{
	Tuple::set(array);
	return *this;
}

Vector& Vector::set(const Vector& copy)
{
	Tuple::set(copy);
	return *this;
}

Vector Vector::cross(const Vector& vec) const
{
	return Vector(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

Vector& Vector::crossEquals(const Vector& vec)
{
	return set(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

float Vector::dot(const Vector& vec) const
{
	return p[0]*vec[0] + p[1]*vec[1] + p[2]*vec[2] + p[3]*vec[3]; 
}


Vector Vector::operator+(const Vector vec) const
{
	return Vector(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector& Vector::operator+=(const Vector vec)
{
	return set(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector Vector::operator-(const Vector vec) const
{
	return Vector(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector& Vector::operator-=(const Vector vec)
{
	return set(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector Vector::operator*(float scalar) const
{
	return Vector(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector& Vector::operator*=(float scalar)
{
	return set(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector Vector::operator-() const
{
	return Vector(-p[0], -p[1], -p[2], p[3]);
}

Vector& Vector::normalize()
{
	if ((float)fabs(p[3])<=EPS) {
		float length = (float)sqrt((double)(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
		return set(p[0]/length,p[1]/length,p[2]/length,0.0f);
	}
	else {
		return set(p[0]/p[3], p[1]/p[3], p[2]/p[3], 1.0f);
	}
}

float Vector::norm()
{
	return (float)sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

bool Vector::isBad()
{
	return (p[0]==0.0f && p[1]==0.0f && p[2]==0.0f && p[3]==0.0f);
}

Vector Vector::badVector()
{
	return Vector(0.0f, 0.0f, 0.0f, 0.0f);
}

Vector* Vector::clone() const
{
	return new Vector(*this);
}

bool Vector::getCorners(double* min, double* max, CCVOpenGLMath::Vector* vCorner)
{
	if( !min || !max || !vCorner ) return false;

	vCorner[0].set(min[0], min[1], min[2], 1);
	vCorner[1].set(max[0], min[1], min[2], 1);
	vCorner[2].set(min[0], max[1], min[2], 1);
	vCorner[3].set(max[0], max[1], min[2], 1);
	vCorner[4].set(min[0], min[1], max[2], 1);
	vCorner[5].set(max[0], min[1], max[2], 1);
	vCorner[6].set(min[0], max[1], max[2], 1);
	vCorner[7].set(max[0], max[1], max[2], 1);

	return true;
}
