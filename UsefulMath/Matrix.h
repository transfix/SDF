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

// Matrix.h: interface for the Matrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
#define AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_

namespace CCVOpenGLMath {

class Quaternion;
class Vector;
class Ray;

class Matrix  
{
public:
	Matrix();
	Matrix(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix(const Quaternion& quat);
	virtual ~Matrix();
	Matrix(const Matrix& copy);
	Matrix& operator=(const Matrix& copy);


	Matrix& set (
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix& set(const Matrix& copy);
	Matrix& reset();
	inline float get(int row, int column) const;
	const float* getMatrix() const;

	Vector operator*(const Vector& vec) const;
	Ray operator*(const Ray& ray) const;
	Matrix operator*(const Matrix& mat) const;
	Matrix& preMultiplication(const Matrix& mat);
	Matrix& postMultiplication(const Matrix& mat);

	Matrix inverse() const;
	Matrix inverseTranspose() const;

	float determinant() const;

	static Matrix rotationX(float angle);
	static Matrix rotationY(float angle);
	static Matrix rotationZ(float angle);
	static Matrix translation(float x, float y, float z);
	static Matrix translation(const Vector& vec);
	static Matrix scale(float x, float y, float z);

protected:
	float m[16];

};

inline float Matrix::get(int row, int column) const
{
	return m[row + column*4];
}

};

#endif // !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
