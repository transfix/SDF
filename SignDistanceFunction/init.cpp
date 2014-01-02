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
#include <math.h>

#include "common.h"

using namespace SDFLibrary;

extern void update_bounding_box(long int current_triangle, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int cur_level);
extern void start_fireworks();
extern void write_octree();

int maxInd;

void build_octree()
{
	double t1, t2;

	t1 = getTime();

	for (int i =0; i<total_triangles; i++)
	{
		update_bounding_box((long)i, minext[0], maxext[0], minext[1], maxext[1],minext[2], maxext[2], 0);
		//update_bounding_box((long)i, 0, size, 0, size, 0, size, 0);
		if (i%1000 == 0)
			printf("%d processed in octree\n", i);
	}

	t2 = getTime();
	printf("Octree constructed for the data in %f seconds\n", (t2-t1));

	//write_octree();
}

void process_triangle(int i)
{
	double p1x, p1y, p1z, p2x, p2y, p2z;
	double nx, ny, nz;
	double denom;
	int v1,v2,v3;

	v1 = surface[i].v1;	v2 = surface[i].v2;	v3 = surface[i].v3;

	//assume that the triangles are consistently oriented V1-V2-V3.
	p1x = vertices[v3].x - vertices[v2].x;
	p1y = vertices[v3].y - vertices[v2].y;
	p1z = vertices[v3].z - vertices[v2].z;
	p2x = vertices[v1].x - vertices[v2].x;
	p2y = vertices[v1].y - vertices[v2].y;
	p2z = vertices[v1].z - vertices[v2].z;

	nx = ( ( p1y * p2z ) - ( p1z * p2y ) );
	ny = ( ( p1z * p2x ) - ( p1x * p2z ) );
	nz = ( ( p1x * p2y ) - ( p1y * p2x ) );

	denom = (double)sqrt( nx*nx + ny*ny + nz*nz );

	nx /= denom;
	ny /= denom;
	nz /= denom;
	normals[i].x = (double)nx;	normals[i].y = (double)ny;	normals[i].z = (double)nz;

	//calculate the Distance of the current Triangle from the Origin
	distances[i] =(double) (-1* ( ( nx * vertices[v1].x ) + ( ny * vertices[v1].y ) + ( nz * vertices[v1].z ) ));
	surface[i].type = -1;
}

void reverse_ptrs()
{
	int i, flag=0;

	for (i=0; i<total_triangles; i++)
	{
		process_triangle(i);
		
		//vertices[ surface[i].v1 ].tris [ vertices[ surface[i].v1 ].trisUsed++ ] = i;
		//vertices[ surface[i].v2 ].tris [ vertices[ surface[i].v2 ].trisUsed++ ] = i;
		//vertices[ surface[i].v3 ].tris [ vertices[ surface[i].v3 ].trisUsed++ ] = i;

		vertices[ surface[i].v1 ].tris.push_back(i); vertices[ surface[i].v1 ].trisUsed++;
		vertices[ surface[i].v2 ].tris.push_back(i); vertices[ surface[i].v2 ].trisUsed++;
		vertices[ surface[i].v3 ].tris.push_back(i); vertices[ surface[i].v3 ].trisUsed++;

#if 0
		if (vertices[ surface[i].v1 ].trisUsed >= MAX_TRIS_PER_VERT) 
		{
			printf("more than %d triangles share this vertex... %d for vert=%d\n", MAX_TRIS_PER_VERT, vertices[ surface[i].v1 ].trisUsed, surface[i].v1);
			flag =1;
		}
		if (vertices[ surface[i].v2 ].trisUsed >= MAX_TRIS_PER_VERT) 
		{
			printf("more than %d triangles share this vertex... %d for vert=%d\n", MAX_TRIS_PER_VERT, vertices[ surface[i].v2 ].trisUsed, surface[i].v2);
			flag =1;
		}
		if (vertices[ surface[i].v3 ].trisUsed >= MAX_TRIS_PER_VERT) 
		{
			printf("more than %d triangles share this vertex... %d for vert=%d\n", MAX_TRIS_PER_VERT, vertices[ surface[i].v3 ].trisUsed, surface[i].v3);
			flag =1;
		}

		//If any of these above statements are printed, then please increase the MAX_TRIS_PER_VERT definition in head.h file and try.
		if (flag ==1)
		{
			printf("Please try changing the MAX_TRIS_PER_VERT variable in <head.h> file and rerun\n");
			exit(0);
		}
#endif
	}
}

void SDFLibrary::adjustData()
{
#if 0
	if (minx < minext[0])
	{
		printf("MinX is changed from %f to %f to accomodate the data\n", minext[0], minx);
		minext[0] = minx;
	}

	if (miny < minext[1])
	{
		printf("miny is changed from %f to %f to accomodate the data\n", minext[1], miny);
		minext[1] = miny;
	}

	if (minz < minext[2])
	{
		printf("minz is changed from %f to %f to accomodate the data\n", minext[2], minz);
		minext[2] = minz;
	}

	if (maxz > maxext[0])
	{
		printf("maxz is changed from %f to %f to accomodate the data\n", maxext[0], maxz);
		maxext[0] = maxz;
	}

	if (maxy > maxext[1])
	{
		printf("maxy is changed from %f to %f to accomodate the data\n", maxext[1], maxy);
		maxext[1] = maxy;
	}

	if (maxz > maxext[2])
	{
		printf("maxz is changed from %f to %f to accomodate the data\n", maxext[2], maxz);
		maxext[2] = maxz;
	}
#endif

	//not sure why span is re-calculated here with the denominator != size-1
	//if i change it to size-1, i get errors in the output... -Joe R.
	span[0] = (maxext[0]-minext[0])/(size);
	span[1] = (maxext[1]-minext[1])/(size);
	span[2] = (maxext[2]-minext[2])/(size);

	printf("\n\nSurface Bounding box is: %f %f %f to %f %f %f \n", minx, miny, minz, maxx, maxy, maxz);
	printf("\nVolume Bounding box is %f %f %f to %f %f %f \n", minext[0], minext[1], minext[2], maxext[0], maxext[1], maxext[2]);


	//Then calculate the normals and back-pointers of the triangles. 
	reverse_ptrs();

	//This wud align them in a consistent manner. ie: all out or all in. :-)
	if (flipNormals)
		start_fireworks();

	//Then build the Octree.
	build_octree();
}

bool setOctree_depth()
{
	switch (size) {
	case (16):
		octree_depth = 4;
	break;
	case (32):
		octree_depth = 5;
	break;
	case (64):
		octree_depth = 6;
	break;
	case (128):
		octree_depth = 7;
	break;
	case (256):
		octree_depth = 8;
	break;
	case (512):
		octree_depth = 9;
	break;
	case (1024):
		octree_depth = 10;
	break;

	default:
		printf("This version can only deal with Volumes of sizes 16, 32, 64, 128, 256, 512 or 1024\n");
		return false;
	}

	return true;
}

bool SDFLibrary::initSDF()
{
    int i, j, k;
	
	MAX_DIST =(double) (size * sqrt(3.0));
	minx = miny = minz = 10000.0;
	maxx = maxy = maxz = -10000.0;
	
	maxInd =-1;
		
	if( !setOctree_depth() ) 
	{
		return false;
	}
	sdf = (cell***) malloc( sizeof(cell**) * (size));
    for (i = 0; i < size; i++)
	{
		sdf[i] = (cell**) malloc( sizeof(cell*) * (size));
		for (j = 0; j < size; j++)
		{
			sdf[i][j] = (cell*) malloc( sizeof(cell) * (size));
			for (k = 0; k < size; k++)
			{
				sdf[i][j][k].useful = 0;
				sdf[i][j][k].type = 1;
				sdf[i][j][k].no = 0;
				sdf[i][j][k].tindex = NULL;					
			}
		}
	}

	k = (size+1)*(size+1)*(size+1);
	values = (voxel*)(malloc(sizeof(voxel) * k));
	bverts = (bool*)(malloc(sizeof(bool)*k));
	queues = (int*)(malloc(sizeof(int)*k));

	for (i=0; i<k; i++)
	{
		values[i].value = (float)MAX_DIST;
		values[i].signe = 0;
		values[i].processed =0;
		values[i].closestV = 0;
		bverts[i] = 0;
	}

	return true;
}

void check_bounds(int i)
{
	if (vertices[i].x < minx) minx = (double) vertices[i].x;
	if (vertices[i].y < miny) miny = (double) vertices[i].y;
	if (vertices[i].z < minz) minz = (double) vertices[i].z;

	if (vertices[i].x > maxx) maxx = (double) vertices[i].x;
	if (vertices[i].y > maxy) maxy = (double) vertices[i].y;
	if (vertices[i].z > maxz) maxz = (double) vertices[i].z;
}

void SDFLibrary::readGeom(int nverts, float* verts, int ntris, int* tris)
{
    int	i;

	total_points = nverts;
	total_triangles = ntris;

	printf("vert= %d and tri = %d \n", total_points,total_triangles);

	vertices = (myVert*) malloc (sizeof (myVert) * total_points);
	surface = (triangle*) malloc (sizeof (triangle) * total_triangles);
	normals = (myPoint*) malloc (sizeof (myPoint) * total_triangles);
	distances = (double*) malloc (sizeof (double) * total_triangles);

	for (i=0; i<total_points; i++)
	{
		vertices[i].x = verts[3*i+0];	vertices[i].y = verts[3*i+1];	vertices[i].z = verts[3*i+2];
		check_bounds(i);
		vertices[i].isNull = 0;		vertices[i].trisUsed =0;
		
		if (!(i% 5000))
			printf("still working on points !!!! %d \n",i);
	}
	
	printf("Finished reading the Vertices.. Now reading the Triangles\n");

	for (i=0; i<total_triangles; i++)
	{
		surface[i].v1 = tris[3*i+0];	surface[i].v2 = tris[3*i+1];	surface[i].v3 = tris[3*i+2];

		if (maxInd < surface[i].v1) maxInd = surface[i].v1;
		if (maxInd < surface[i].v2) maxInd = surface[i].v2;
		if (maxInd < surface[i].v3) maxInd = surface[i].v3;

		if (!(i% 5000))
			printf("still working on Triangles !!!! %d \n",i);
	}
	
	printf("Bounding box is: %f %f %f to %f %f %f \n", minx, miny, minz, maxx, maxy, maxz);
}
