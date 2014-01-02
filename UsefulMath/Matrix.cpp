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
// Matrix.cpp: implementation of the Matrix class.
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"
#include "Vector.h"
#include "Ray.h"
#include <math.h>

using CCVOpenGLMath::Vector;
using CCVOpenGLMath::Ray;
//using CCVOpenGLMath::Quaternion;
using CCVOpenGLMath::Matrix;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Matrix::Matrix()
{
	set(1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix::Matrix(
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33
	)
{
	m[0]=m00; m[1]=m10; m[2]=m20; m[3]=m30;
	m[4]=m01; m[5]=m11; m[6]=m21; m[7]=m31;
	m[8]=m02; m[9]=m12; m[10]=m22; m[11]=m32;
	m[12]=m03; m[13]=m13; m[14]=m23; m[15]=m33;
}

Matrix::~Matrix()
{

}

Matrix::Matrix(const Matrix& copy)
{
	set(copy);
}

Matrix& Matrix::operator=(const Matrix& copy)
{
	return set(copy);
}


Matrix& Matrix::set (
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33
	)
{
	m[0]=m00; m[1]=m10; m[2]=m20; m[3]=m30;
	m[4]=m01; m[5]=m11; m[6]=m21; m[7]=m31;
	m[8]=m02; m[9]=m12; m[10]=m22; m[11]=m32;
	m[12]=m03; m[13]=m13; m[14]=m23; m[15]=m33;
	return *this;
}

Matrix& Matrix::set(const Matrix& copy)
{
	if (this!=&copy) {
		set(
			copy.m[0], copy.m[4], copy.m[8], copy.m[12],
			copy.m[1], copy.m[5], copy.m[9], copy.m[13],
			copy.m[2], copy.m[6], copy.m[10], copy.m[14],
			copy.m[3], copy.m[7], copy.m[11], copy.m[15]
			);
	}
	return *this;
}

Matrix& Matrix::reset()
{
	set(
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0);
	return *this;
}

const float* Matrix::getMatrix() const
{
	return m;
}

Vector Matrix::operator*(const Vector& vec) const
{
	return Vector(
		m[0]*vec[0]+m[4]*vec[1]+m[8]*vec[2]+m[12]*vec[3],
		m[1]*vec[0]+m[5]*vec[1]+m[9]*vec[2]+m[13]*vec[3],
		m[2]*vec[0]+m[6]*vec[1]+m[10]*vec[2]+m[14]*vec[3],
		m[3]*vec[0]+m[7]*vec[1]+m[11]*vec[2]+m[15]*vec[3]
		);
}

Ray Matrix::operator*(const Ray& ray) const
{
	return Ray((*this)*ray.m_Origin, (*this)*ray.m_Dir);
}

Matrix Matrix::operator*(const Matrix& mat) const
{
	return Matrix(
		m[0]*mat.m[0]+m[4]*mat.m[1]+m[8]*mat.m[2]+m[12]*mat.m[3],
		m[0]*mat.m[4]+m[4]*mat.m[5]+m[8]*mat.m[6]+m[12]*mat.m[7],
		m[0]*mat.m[8]+m[4]*mat.m[9]+m[8]*mat.m[10]+m[12]*mat.m[11],
		m[0]*mat.m[12]+m[4]*mat.m[13]+m[8]*mat.m[14]+m[12]*mat.m[15],

		m[1]*mat.m[0]+m[5]*mat.m[1]+m[9]*mat.m[2]+m[13]*mat.m[3],
		m[1]*mat.m[4]+m[5]*mat.m[5]+m[9]*mat.m[6]+m[13]*mat.m[7],
		m[1]*mat.m[8]+m[5]*mat.m[9]+m[9]*mat.m[10]+m[13]*mat.m[11],
		m[1]*mat.m[12]+m[5]*mat.m[13]+m[9]*mat.m[14]+m[13]*mat.m[15],

		m[2]*mat.m[0]+m[6]*mat.m[1]+m[10]*mat.m[2]+m[14]*mat.m[3],
		m[2]*mat.m[4]+m[6]*mat.m[5]+m[10]*mat.m[6]+m[14]*mat.m[7],
		m[2]*mat.m[8]+m[6]*mat.m[9]+m[10]*mat.m[10]+m[14]*mat.m[11],
		m[2]*mat.m[12]+m[6]*mat.m[13]+m[10]*mat.m[14]+m[14]*mat.m[15],

		m[3]*mat.m[0]+m[7]*mat.m[1]+m[11]*mat.m[2]+m[15]*mat.m[3],
		m[3]*mat.m[4]+m[7]*mat.m[5]+m[11]*mat.m[6]+m[15]*mat.m[7],
		m[3]*mat.m[8]+m[7]*mat.m[9]+m[11]*mat.m[10]+m[15]*mat.m[11],
		m[3]*mat.m[12]+m[7]*mat.m[13]+m[11]*mat.m[14]+m[15]*mat.m[15]
		);
}

Matrix& Matrix::preMultiplication(const Matrix& mat)
{
	return set(mat*(*this));
}

Matrix& Matrix::postMultiplication(const Matrix& mat)
{
	return set((*this)*mat);
}

Matrix Matrix::inverse() const
{
	Matrix ret;
	float det = determinant();
	if (det!=0.0) {
		ret.set(
			(get(1, 2)*get(2, 3)*get(3, 1) - get(1, 3)*get(2, 2)*get(3, 1) + get(1, 3)*get(2, 1)*get(3, 2) - get(1, 1)*get(2, 3)*get(3, 2) - get(1, 2)*get(2, 1)*get(3, 3) + get(1, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(2, 2)*get(3, 1) - get(0, 2)*get(2, 3)*get(3, 1) - get(0, 3)*get(2, 1)*get(3, 2) + get(0, 1)*get(2, 3)*get(3, 2) + get(0, 2)*get(2, 1)*get(3, 3) - get(0, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(1, 3)*get(3, 1) - get(0, 3)*get(1, 2)*get(3, 1) + get(0, 3)*get(1, 1)*get(3, 2) - get(0, 1)*get(1, 3)*get(3, 2) - get(0, 2)*get(1, 1)*get(3, 3) + get(0, 1)*get(1, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(2, 1) - get(0, 2)*get(1, 3)*get(2, 1) - get(0, 3)*get(1, 1)*get(2, 2) + get(0, 1)*get(1, 3)*get(2, 2) + get(0, 2)*get(1, 1)*get(2, 3) - get(0, 1)*get(1, 2)*get(2, 3))/det,

			(get(1, 3)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 0)*get(3, 2) + get(1, 0)*get(2, 3)*get(3, 2) + get(1, 2)*get(2, 0)*get(3, 3) - get(1, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 2)*get(3, 0) + get(0, 3)*get(2, 0)*get(3, 2) - get(0, 0)*get(2, 3)*get(3, 2) - get(0, 2)*get(2, 0)*get(3, 3) + get(0, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 0)*get(3, 2) + get(0, 0)*get(1, 3)*get(3, 2) + get(0, 2)*get(1, 0)*get(3, 3) - get(0, 0)*get(1, 2)*get(3, 3))/det,
			(get(0, 2)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 2)*get(2, 0) + get(0, 3)*get(1, 0)*get(2, 2) - get(0, 0)*get(1, 3)*get(2, 2) - get(0, 2)*get(1, 0)*get(2, 3) + get(0, 0)*get(1, 2)*get(2, 3))/det,

			(get(1, 1)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 1)*get(3, 0) + get(1, 3)*get(2, 0)*get(3, 1) - get(1, 0)*get(2, 3)*get(3, 1) - get(1, 1)*get(2, 0)*get(3, 3) + get(1, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 3)*get(2, 1)*get(3, 0) - get(0, 1)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 0)*get(3, 1) + get(0, 0)*get(2, 3)*get(3, 1) + get(0, 1)*get(2, 0)*get(3, 3) - get(0, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 1)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 1)*get(3, 0) + get(0, 3)*get(1, 0)*get(3, 1) - get(0, 0)*get(1, 3)*get(3, 1) - get(0, 1)*get(1, 0)*get(3, 3) + get(0, 0)*get(1, 1)*get(3, 3))/det,
			(get(0, 3)*get(1, 1)*get(2, 0) - get(0, 1)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 0)*get(2, 1) + get(0, 0)*get(1, 3)*get(2, 1) + get(0, 1)*get(1, 0)*get(2, 3) - get(0, 0)*get(1, 1)*get(2, 3))/det,

			(get(1, 2)*get(2, 1)*get(3, 0) - get(1, 1)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 0)*get(3, 1) + get(1, 0)*get(2, 2)*get(3, 1) + get(1, 1)*get(2, 0)*get(3, 2) - get(1, 0)*get(2, 1)*get(3, 2))/det,
			(get(0, 1)*get(2, 2)*get(3, 0) - get(0, 2)*get(2, 1)*get(3, 0) + get(0, 2)*get(2, 0)*get(3, 1) - get(0, 0)*get(2, 2)*get(3, 1) - get(0, 1)*get(2, 0)*get(3, 2) + get(0, 0)*get(2, 1)*get(3, 2))/det,
			(get(0, 2)*get(1, 1)*get(3, 0) - get(0, 1)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 0)*get(3, 1) + get(0, 0)*get(1, 2)*get(3, 1) + get(0, 1)*get(1, 0)*get(3, 2) - get(0, 0)*get(1, 1)*get(3, 2))/det,
			(get(0, 1)*get(1, 2)*get(2, 0) - get(0, 2)*get(1, 1)*get(2, 0) + get(0, 2)*get(1, 0)*get(2, 1) - get(0, 0)*get(1, 2)*get(2, 1) - get(0, 1)*get(1, 0)*get(2, 2) + get(0, 0)*get(1, 1)*get(2, 2))/det
		);
	}
	// if det==0.0, ret will be identity
	return ret;
}

Matrix Matrix::inverseTranspose() const
{
	Matrix ret;
	float det = determinant();
	if (det!=0.0) {
		ret.set(
			(get(1, 2)*get(2, 3)*get(3, 1) - get(1, 3)*get(2, 2)*get(3, 1) + get(1, 3)*get(2, 1)*get(3, 2) - get(1, 1)*get(2, 3)*get(3, 2) - get(1, 2)*get(2, 1)*get(3, 3) + get(1, 1)*get(2, 2)*get(3, 3))/det,
			(get(1, 3)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 0)*get(3, 2) + get(1, 0)*get(2, 3)*get(3, 2) + get(1, 2)*get(2, 0)*get(3, 3) - get(1, 0)*get(2, 2)*get(3, 3))/det,
			(get(1, 1)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 1)*get(3, 0) + get(1, 3)*get(2, 0)*get(3, 1) - get(1, 0)*get(2, 3)*get(3, 1) - get(1, 1)*get(2, 0)*get(3, 3) + get(1, 0)*get(2, 1)*get(3, 3))/det,
			(get(1, 2)*get(2, 1)*get(3, 0) - get(1, 1)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 0)*get(3, 1) + get(1, 0)*get(2, 2)*get(3, 1) + get(1, 1)*get(2, 0)*get(3, 2) - get(1, 0)*get(2, 1)*get(3, 2))/det,

			(get(0, 3)*get(2, 2)*get(3, 1) - get(0, 2)*get(2, 3)*get(3, 1) - get(0, 3)*get(2, 1)*get(3, 2) + get(0, 1)*get(2, 3)*get(3, 2) + get(0, 2)*get(2, 1)*get(3, 3) - get(0, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 2)*get(3, 0) + get(0, 3)*get(2, 0)*get(3, 2) - get(0, 0)*get(2, 3)*get(3, 2) - get(0, 2)*get(2, 0)*get(3, 3) + get(0, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(2, 1)*get(3, 0) - get(0, 1)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 0)*get(3, 1) + get(0, 0)*get(2, 3)*get(3, 1) + get(0, 1)*get(2, 0)*get(3, 3) - get(0, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 1)*get(2, 2)*get(3, 0) - get(0, 2)*get(2, 1)*get(3, 0) + get(0, 2)*get(2, 0)*get(3, 1) - get(0, 0)*get(2, 2)*get(3, 1) - get(0, 1)*get(2, 0)*get(3, 2) + get(0, 0)*get(2, 1)*get(3, 2))/det,

			(get(0, 2)*get(1, 3)*get(3, 1) - get(0, 3)*get(1, 2)*get(3, 1) + get(0, 3)*get(1, 1)*get(3, 2) - get(0, 1)*get(1, 3)*get(3, 2) - get(0, 2)*get(1, 1)*get(3, 3) + get(0, 1)*get(1, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 0)*get(3, 2) + get(0, 0)*get(1, 3)*get(3, 2) + get(0, 2)*get(1, 0)*get(3, 3) - get(0, 0)*get(1, 2)*get(3, 3))/det,
			(get(0, 1)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 1)*get(3, 0) + get(0, 3)*get(1, 0)*get(3, 1) - get(0, 0)*get(1, 3)*get(3, 1) - get(0, 1)*get(1, 0)*get(3, 3) + get(0, 0)*get(1, 1)*get(3, 3))/det,
			(get(0, 2)*get(1, 1)*get(3, 0) - get(0, 1)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 0)*get(3, 1) + get(0, 0)*get(1, 2)*get(3, 1) + get(0, 1)*get(1, 0)*get(3, 2) - get(0, 0)*get(1, 1)*get(3, 2))/det,
			
			(get(0, 3)*get(1, 2)*get(2, 1) - get(0, 2)*get(1, 3)*get(2, 1) - get(0, 3)*get(1, 1)*get(2, 2) + get(0, 1)*get(1, 3)*get(2, 2) + get(0, 2)*get(1, 1)*get(2, 3) - get(0, 1)*get(1, 2)*get(2, 3))/det,
			(get(0, 2)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 2)*get(2, 0) + get(0, 3)*get(1, 0)*get(2, 2) - get(0, 0)*get(1, 3)*get(2, 2) - get(0, 2)*get(1, 0)*get(2, 3) + get(0, 0)*get(1, 2)*get(2, 3))/det,
			(get(0, 3)*get(1, 1)*get(2, 0) - get(0, 1)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 0)*get(2, 1) + get(0, 0)*get(1, 3)*get(2, 1) + get(0, 1)*get(1, 0)*get(2, 3) - get(0, 0)*get(1, 1)*get(2, 3))/det,
			(get(0, 1)*get(1, 2)*get(2, 0) - get(0, 2)*get(1, 1)*get(2, 0) + get(0, 2)*get(1, 0)*get(2, 1) - get(0, 0)*get(1, 2)*get(2, 1) - get(0, 1)*get(1, 0)*get(2, 2) + get(0, 0)*get(1, 1)*get(2, 2))/det
		);
	}
	// if det==0.0, ret will be identity
	return ret;
}

float Matrix::determinant() const
{
	double ret;
	ret = 
		get(0, 3) * get(1, 2) * get(2, 1) * get(3, 0)-get(0, 2) * get(1, 3) * get(2, 1) * get(3, 0)-get(0, 3) * get(1, 1) * get(2, 2) * get(3, 0)+get(0, 1) * get(1, 3) * get(2, 2) * get(3, 0)+
		get(0, 2) * get(1, 1) * get(2, 3) * get(3, 0)-get(0, 1) * get(1, 2) * get(2, 3) * get(3, 0)-get(0, 3) * get(1, 2) * get(2, 0) * get(3, 1)+get(0, 2) * get(1, 3) * get(2, 0) * get(3, 1)+
		get(0, 3) * get(1, 0) * get(2, 2) * get(3, 1)-get(0, 0) * get(1, 3) * get(2, 2) * get(3, 1)-get(0, 2) * get(1, 0) * get(2, 3) * get(3, 1)+get(0, 0) * get(1, 2) * get(2, 3) * get(3, 1)+
		get(0, 3) * get(1, 1) * get(2, 0) * get(3, 2)-get(0, 1) * get(1, 3) * get(2, 0) * get(3, 2)-get(0, 3) * get(1, 0) * get(2, 1) * get(3, 2)+get(0, 0) * get(1, 3) * get(2, 1) * get(3, 2)+
		get(0, 1) * get(1, 0) * get(2, 3) * get(3, 2)-get(0, 0) * get(1, 1) * get(2, 3) * get(3, 2)-get(0, 2) * get(1, 1) * get(2, 0) * get(3, 3)+get(0, 1) * get(1, 2) * get(2, 0) * get(3, 3)+
		get(0, 2) * get(1, 0) * get(2, 1) * get(3, 3)-get(0, 0) * get(1, 2) * get(2, 1) * get(3, 3)-get(0, 1) * get(1, 0) * get(2, 2) * get(3, 3)+get(0, 0) * get(1, 1) * get(2, 2) * get(3, 3);
	return (float)ret;
}

Matrix Matrix::rotationX(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, ca, sa, 0.0f,
		0.0f, -sa, ca, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::rotationY(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(ca, 0.0f, -sa, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		sa, 0.0f, ca, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::rotationZ(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(ca, sa, 0.0f, 0.0f,
		-sa, ca, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::translation(float x, float y, float z)
{
	return Matrix(1.0f, 0.0f, 0.0f, x,
		0.0f, 1.0f, 0.0f, y,
		0.0f, 0.0f, 1.0f, z,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::translation(const Vector& vec)
{
	return Matrix(1.0f, 0.0f, 0.0f, vec[0],
		0.0f, 1.0f, 0.0f, vec[1],
		0.0f, 0.0f, 1.0f, vec[2],
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::scale(float x, float y, float z)
{
	return Matrix(x, 0.0f, 0.0f, 0.0f,
		0.0f, y, 0.0f, 0.0f,
		0.0f, 0.0f, z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}





