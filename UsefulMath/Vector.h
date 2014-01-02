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

// Vector.h: interface for the Vector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
#define AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_

#include "Tuple.h"

namespace CCVOpenGLMath {

class Vector : public Tuple  
{
public:
	Vector(float x, float y, float z, float w);
	Vector(float* array);
	Vector();
	virtual ~Vector();
	Vector(const Vector& copy);
	Vector& operator=(const Vector& copy);

	Vector& set(float x, float y, float z, float w);
	Vector& set(float* array);
	Vector& set(const Vector& copy);

	Vector cross(const Vector& vec) const;
	Vector& crossEquals(const Vector& vec);
	float dot(const Vector& vec) const;

	Vector operator+(const Vector vec) const;
	Vector& operator+=(const Vector vec);
	Vector operator-(const Vector vec) const;
	Vector& operator-=(const Vector vec);

	Vector operator*(float scalar) const;
	Vector& operator*=(float scalar);

	Vector operator-() const;

	Vector& normalize();
	float norm();

	bool isBad();

	static Vector badVector();
	static bool getCorners(double* min, double* max, CCVOpenGLMath::Vector* vCorner);
	virtual Vector* clone() const;
};

};

#endif // !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
