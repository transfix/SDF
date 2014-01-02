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
/*     Author:     Lalit Karlapalem <ckl@ices.utexas.edu>         2004-2005  */
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


#ifndef CCV_SDF_HEAD_H
#define CCV_SDF_HEAD_H

#include <vector>

namespace SDFLibrary {

	#define MAX_TRIS_PER_VERT 100
	
	typedef struct _Pt_
	{
		double x;
		double y;
		double z;
		char isNull;
		
	}  myPoint;

	typedef struct _Vt_
	{
		double x;
		double y;
		double z;
		char isNull;

	  //int tris[MAX_TRIS_PER_VERT]; //not more than MAX_TRIS_PER_VERT triangles can share a vertex.
	        std::vector<int> tris;
		int trisUsed; //elements used in the above array.

	}  myVert;


	typedef struct _tri_
	{
		int v1;
		int v2;
		int v3;

		int type; // default = -1; done =1.	wrong =3;
	} triangle;

	typedef struct listnodedef
	{
		int index;	//index of the triangle
		struct listnodedef* next;
	} listnode;

	typedef struct nodedef 
	{
		char useful;	//  0 - no triangles in it	; 1 - there are triangles in it
		char type;		//	0 - interior node		; 1 - leaf node, containing triangles
		long int no;
		listnode* tindex;		
	} cell;

	typedef struct 
	{
		double ox;
		double oy;
		double oz;

		double dx; 
		double dy;
		double dz;
	} ray; 

	typedef struct _voxel_
	{
		float value;
		signed char signe;  //-1 = inside,		1 = outside	
		bool processed;		// 1 = propagated distance FROM here. 0 = not
		int closestV;		//the closest triangle on the surface
	}voxel;

	
}; //namespace SDFLibrary

#endif
