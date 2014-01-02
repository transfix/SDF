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


#ifndef CCV_SDF_COMMON_H
#define CCV_SDF_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <map>

#include "head.h"

namespace SDFLibrary {


bool initSDF();

void readGeom(int nverts, float* verts, int ntris, int* tris);

void adjustData();

void compute();



int isEqual (double one, double two);

int isZero(double num);

int isNegative(double num);

int isBetween(double one, double two, double num);

int isZero(SDFLibrary::myPoint one);

int isSame(SDFLibrary::myPoint one, SDFLibrary::myPoint two);

void init_all_vars();



void propagate_left(int i, int j, int k);

void propagate_bottom(int i, int j, int k);

void propagate_inside(int i, int j, int k);

void propagate_right(int i, int j, int k);

void propagate_top(int i, int j, int k);

void propagate_outside(int i, int j, int k);

void apply_distance_transform(int vi, int vj, int vk);

void insert_bound_vert(int vert);




int index2vert(int i, int j, int k);

void _vert2index(int c, int &i, int &j, int &k);

int index2cell(int i, int j, int k);

void _cell2index(int c, int &i, int &j, int &k);

double xCoord(int i);

double yCoord(int i);

double zCoord(int i);

void object2octree(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int &ci, int &cj, int &ck);

double getTime();

	extern double MAX_DIST;
	extern int size;

	extern triangle* surface;
	extern myVert* vertices;
	extern myPoint* normals;
	extern double* distances;
	extern cell*** sdf;
	extern voxel* values;
	extern int total_points, total_triangles, all_verts_touched;
	extern double minx, miny, minz, maxx, maxy, maxz;

	extern double TOLERANCE;

	extern int octree_depth;
	extern int flipNormals;

	extern bool *bverts;
	extern int *queues;

	extern double minext[3];
	extern double maxext[3];
	extern double span[3];
	
}; //namespace SDFLibrary

#endif