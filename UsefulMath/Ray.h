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
// Ray.h: interface for the Ray class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAY_H__8760CFFE_A656_460F_A893_4C81FF806264__INCLUDED_)
#define AFX_RAY_H__8760CFFE_A656_460F_A893_4C81FF806264__INCLUDED_

#include "Vector.h"

namespace CCVOpenGLMath {

	class Ray  
	{
	public:
		Ray();
		Ray(const Vector& origin, const Vector& dir);
		virtual ~Ray();

		Vector getPointOnRay(float t) const;

		float nearestTOnXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		float nearestTOnYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		float nearestTOnZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;

		Vector nearestPointOnXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		Vector nearestPointOnYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		Vector nearestPointOnZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;

		float distanceToXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		float distanceToYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
		float distanceToZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;

		bool intersectSphere( Vector center, float radius, Vector *point1, Vector* point2, float *distance1, float* distance2 );

		Vector m_Origin;
		Vector m_Dir;

	};

};

#endif // !defined(AFX_RAY_H__8760CFFE_A656_460F_A893_4C81FF806264__INCLUDED_)
