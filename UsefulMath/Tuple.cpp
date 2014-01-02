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

// Tuple.cpp: implementation of the Tuple class.
//
//////////////////////////////////////////////////////////////////////

#include "Tuple.h"

using CCVOpenGLMath::Tuple;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Tuple::Tuple(float x, float y, float z, float w)
{
	set(x,y,z,w);
}

Tuple::Tuple()
{
	set(0.0, 0.0, 0.0, 0.0);
}

Tuple::~Tuple()
{

}

Tuple::Tuple(const Tuple& copy)
{
	set(copy);
}

Tuple& Tuple::operator=(const Tuple& copy)
{
	return set(copy);
}

Tuple& Tuple::set(float x, float y, float z, float w)
{
	p[0] = x;
	p[1] = y;
	p[2] = z;
	p[3] = w;
	return *this;
}

Tuple& Tuple::set(float* array)
{
	p[0] = array[0];
	p[1] = array[1];
	p[2] = array[2];
	p[3] = array[3];
	return *this;
}


Tuple& Tuple::set(const Tuple& copy)
{
	if (this!=&copy) {
		p[0] = copy.p[0];
		p[1] = copy.p[1];
		p[2] = copy.p[2];
		p[3] = copy.p[3];
	}
	return *this;
}



float& Tuple::operator[](unsigned int i)
{
	return p[i];
}


const float& Tuple::operator[](unsigned int i) const
{
	return p[i];	
}


