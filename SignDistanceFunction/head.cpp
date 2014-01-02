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

#include <stdio.h>
#include <stdlib.h>

#include "common.h"

namespace SDFLibrary{

	double MAX_DIST;
	int size;
	int total_points, total_triangles, all_verts_touched;
	double minx, miny, minz, maxx, maxy, maxz;
	double TOLERANCE;
	int octree_depth;
	int flipNormals;

	SDFLibrary::triangle* surface;
	SDFLibrary::myVert* vertices;
	SDFLibrary::myPoint* normals;
	SDFLibrary::cell*** sdf;
	SDFLibrary::voxel* values;
	double* distances;
	
	char *ifname;

	bool* bverts;
	int* queues;

	double minext[3];
	double maxext[3];
	double span[3];
};

void SDFLibrary::init_all_vars()
{
	SDFLibrary::TOLERANCE = 1e-5;
	SDFLibrary::size =64;
	SDFLibrary::flipNormals =0;
	
	SDFLibrary::ifname= NULL;
	SDFLibrary::surface = NULL;
	SDFLibrary::vertices = NULL;
	SDFLibrary::normals = NULL;
	SDFLibrary::distances = NULL;
	SDFLibrary::sdf = NULL;
	SDFLibrary::values = NULL;
	SDFLibrary::bverts = NULL;
	SDFLibrary::queues = NULL;

	SDFLibrary::minext[0] = SDFLibrary::minext[1] = SDFLibrary::minext[2] = 10000.00;
	SDFLibrary::maxext[0] = SDFLibrary::maxext[1] = SDFLibrary::maxext[2] = -10000.00;
	SDFLibrary::span[0] = SDFLibrary::span[1] = SDFLibrary::span[2] = 1.0;

	SDFLibrary::total_points = SDFLibrary::total_triangles = SDFLibrary::all_verts_touched= 0;
}