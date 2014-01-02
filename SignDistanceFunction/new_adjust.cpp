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

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <math.h>
#include <stdlib.h>

#include "common.h"

using namespace SDFLibrary;

int* neighbors;
int usedNeighs, prevUsed, total_done;

std::map<int,int> myMap;
std::map<int,int>::const_iterator iter;


void re_orient_all()
{
	int i, closestTri;
	int inside_point = -1;
	double err, dist;
	//double lamda = -0.5; //The dist that the pt is going to be moved inside the surface, along the normal.
	
	err = size*size*size;

	for (i=0; i<total_points; i++)
	{
		dist =0.0;

		dist += (vertices[i].x - minx) * (vertices[i].x - minx);
		dist += (vertices[i].y - miny) * (vertices[i].y - miny);
		dist += (vertices[i].z - minz) * (vertices[i].z - minz);

		if (fabs(dist) < err)
		{
			err = fabs(dist);
			inside_point =i;
		}
	}

	printf("min is %f %f %f and closest pt is %f %f %f\n", minx, miny, minz,
			vertices[inside_point].x, vertices[inside_point].y, vertices[inside_point].z);

	//Then chq if Origin (0,0,0) lies on the outisde of the triangle at the Vertex inside_point or not.
	for (i=0; i<total_triangles; i++)
	{
		if ( (surface[i].v1 == inside_point) || (surface[i].v2 == inside_point) ||  (surface[i].v3 == inside_point) )
			break;
	}
	closestTri =i;

	//The closest distance 'tween the surface[closesetTri] and the Origin is only the distance of the plane from O
	if (distances[closestTri] > 0)
	{
		printf("Normals are correctly oriented\n");
	}
	else
	{
		for (i=0; i<total_triangles; i++)
		{
			normals[i].x *= -1;	normals[i].y *= -1;	normals[i].z *= -1; distances[i] *= -1; 
		}
		printf("Normals were flipped again to be correctly oriented\n");
	}
}

int isAligned (int ver1, int ver2)
{
	if (ver1 == 1)
		if (ver2 == 2)	return 1;
		else return 0;

	if (ver1 == 2)
		if (ver2 == 3)	return 1;
		else return 0;

	if (ver1 == 3)
		if (ver2 == 1)	return 1;
		else return 0;

	return -1; //its an error, but, just to make the compiler happy. :-(
}

void exchangeVerts(int tri, int ver1, int ver2)
{
	if (surface[tri].v1 == ver1)
	{
		surface[tri].v1 = ver2;
		if (surface[tri].v2 == ver2)
			surface[tri].v2 = ver1;
		else	surface[tri].v3 = ver1;
	}
	else if (surface[tri].v2 == ver1)
	{
		surface[tri].v2 = ver2;
		if (surface[tri].v1 == ver2)
			surface[tri].v1 = ver1;
		else	surface[tri].v3 = ver1;
	}
	else if (surface[tri].v3 == ver1)
	{
		surface[tri].v3 = ver2;
		if (surface[tri].v1 == ver2)
			surface[tri].v1 = ver1;
		else	surface[tri].v2 = ver1;
	}
}

//also, change the dammed order of the vertices for the triangle.
int triangle_angles(int one, int two, int ver1, int ver2)
{
	int v1, v2, c1, c2;

	v1= v2 = c1 = c2= -1;

	if (surface[one].v1 == ver1)	v1 =1;
	if (surface[one].v1 == ver2)	v2 =1;
	if (surface[one].v2 == ver1)	v1 =2;
	if (surface[one].v2 == ver2)	v2 =2;
	if (surface[one].v3 == ver1)	v1 =3;
	if (surface[one].v3 == ver2)	v2 =3;

	if (surface[two].v1 == ver1)	c1 =1;
	if (surface[two].v1 == ver2)	c2 =1;
	if (surface[two].v2 == ver1)	c1 =2;
	if (surface[two].v2 == ver2)	c2 =2;
	if (surface[two].v3 == ver1)	c1 =3;
	if (surface[two].v3 == ver2)	c2 =3;

	if ( (v1 == -1) || (v2 == -1) || (c1 == -1) || (c2 == -1) )
	{
		printf("some err in <triangle_angles> : %d %d %d %d\n", one, two, ver1, ver2);
		return 1; //wot 2 return... :-(
	}
	
	if (isAligned(v1, v2))
	{
		if (isAligned(c1, c2))
		{
			//problemo.
			exchangeVerts(two, ver1, ver2);
			return 0;
		}
		else
		{
			//no problemo
			return 1;
		}
	}
	else
	{
		if (isAligned(c1, c2))
		{
			//no problemo.
			return 1;
		}
		else
		{
			//problemo
			exchangeVerts(two, ver1, ver2);
			return 0;
		}
	}
}


void insert_tri(int tri)
{
	if (surface[tri].type == -1) return;

	iter = myMap.find(tri);
	if(iter == myMap.end()) //ie not found
	{
		myMap[tri] = tri;
		neighbors[usedNeighs++] = tri;
		total_done++;
	}
}

void align_us(int with, int what, int vert)
{
	int i, j, flag=-1;
	int v1[3], v2[3];
	
	if (surface[what].type != -1) return;

	v1[0] = surface[with].v1;	v1[1] = surface[with].v2;	v1[2] = surface[with].v3;
	v2[0] = surface[what].v1;	v2[1] = surface[what].v2;	v2[2] = surface[what].v3;
	
	for (i=0; i<3; i++)
	{
		if (v1[i] == vert) continue;

		for (j=0; j<3; j++)
		{
			if (v2[j] == vert) continue;

			if (v1[i] == v2[j])
				flag = v1[i];
		}
	}

	if (flag == -1)
		return;

	//then compare the two triangles.
	if (triangle_angles(with, what, vert, flag))
		surface[what].type = surface[with].type;
	else
	{
		normals[what].x *= -1;	normals[what].y *= -1;	normals[what].z *= -1;
		distances[what] *= -1; //need to re-calculate the distances also.
		surface[what].type = !(surface[with].type);
	}

	//Then insert this triangle into the NEIGHBORS array.
	insert_tri(what);
}

void orient_vert(int tri, int vert)
{
	int i;

	for (i=0; i<vertices[vert].trisUsed; i++)
	{
		if (tri != vertices[vert].tris[i])
			align_us(tri, vertices[vert].tris[i], vert);
	}
}

void correct_tri(int tri)
{
	orient_vert(tri, surface[tri].v1);
	orient_vert(tri, surface[tri].v2);
	orient_vert(tri, surface[tri].v3);			
}

//This one is called for each unconnected component. Assume that the first triangle is pointing outwards...
void getNextComponent()
{
	int i;

	for (i=0; i<total_triangles; i++)
	{
		if (surface[i].type == -1)
			break;
	}
	
	surface[i].type =1;
	insert_tri(i);
	prevUsed =usedNeighs;
}

void start_fireworks()
{
	int i, j, lastone;
	int* tarray;

	neighbors = (int*) (malloc(sizeof(int) * total_triangles));
	tarray = (int*) (malloc(sizeof(int) * total_triangles));

	printf("\n<start_fireworks> started...\n");

	myMap.clear();
	lastone = usedNeighs= total_done =0;

	while (1)
	{
		prevUsed = usedNeighs;
		printf("still processing with %d Triangles\n", prevUsed);
		
		if (lastone == prevUsed)
			getNextComponent();
		else
			lastone = prevUsed;

		for (i=0; i<prevUsed; i++)
			correct_tri(neighbors[i]);

		if (total_done == total_triangles)
		{
			printf("The reqd normal flipping is done.\n");
			break;
		}

		j=0;
		for (iter=myMap.begin(); iter!=myMap.end(); ++iter)
		{
			neighbors[j++] = (*iter).first;
		}
		usedNeighs =j;
	}

	free(neighbors);
	free(tarray);
	myMap.clear();

	re_orient_all();

	printf("<start_fireworks> over...\n");
}

