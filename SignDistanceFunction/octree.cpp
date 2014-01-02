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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"

using namespace SDFLibrary;

long int INF = 9999999;
#define PI 3.14159

double n_dotv( double x, double y, double z, ray r , double temp123);
myPoint inbox ( ray r, myPoint p, double dist, double* t);
int ray_polygon_intersection (ray r, int tri);
myPoint normalize(double x, double y, double z) ;
int max_3( double x, double y, double z );
int inside_cube(ray r, double xmin, double xmax, double ymin, double ymax, int flag);
int point_in_polygon(myPoint result, int tri); 
void update_boundary_vertices(int cx, int cy, int cz);

int within( int tri, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax )
{
	int v1, v2, v3;
	double t;
	double x0, y0, z0, x1, y1, z1, x2, y2, z2;
	myPoint p,one,two,three;
	ray r;
	
	v1 = surface[tri].v1;	v2 = surface[tri].v2;	v3 = surface[tri].v3;

	// find if completely on some side !
	x0 = vertices[v1].x;
	y0 = vertices[v1].y;
	z0 = vertices[v1].z;

	x1 = vertices[v2].x;
	y1 = vertices[v2].y;
	z1 = vertices[v2].z;

	x2 = vertices[v3].x;
	y2 = vertices[v3].y;
	z2 = vertices[v3].z;

	// If all 3 vertices on same side, then return 0
	if( (x0<xmin) && (x1<xmin) && (x2<xmin) ) return 0;
	if( (x0>xmax) && (x1>xmax) && (x2>xmax) ) return 0;

	if( (y0<ymin) && (y1<ymin) && (y2<ymin) ) return 0;
	if( (y0>ymax) && (y1>ymax) && (y2>ymax) ) return 0;

	if( (z0<zmin) && (z1<zmin) && (z2<zmin) ) return 0;
	if( (z0>zmax) && (z1>zmax) && (z2>zmax) ) return 0;

	//if either of the 3 points are inside the cube, then the Triangle MUST intersect the cube....
	if ( (isBetween(xmin, xmax, x0)) && (isBetween(ymin, ymax, y0)) && (isBetween(zmin, zmax, z0)) ) return 1;
	if ( (isBetween(xmin, xmax, x1)) && (isBetween(ymin, ymax, y1)) && (isBetween(zmin, zmax, z1)) ) return 1;
	if ( (isBetween(xmin, xmax, x2)) && (isBetween(ymin, ymax, y2)) && (isBetween(zmin, zmax, z2)) ) return 1;

	//then the dammed cases when the Triangle intersects any face (edge) of the cube or vice versa...
	//A) Triangle with the cube...
	r.ox = (double)x0;
	r.oy = (double)y0;
	r.oz = (double)z0;
	r.dx = (double)(x1-x0);
	r.dy = (double)(y1-y0);
	r.dz = (double)(z1-z0);

	//1) a> edge1 with face 1 ...
	one = normalize(xmax- xmin,0,0);
	p = inbox(r,one,-xmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1; 
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//1) b> edge1 with face 4...
	one = normalize(xmax- xmin,0,0);
	p = inbox(r,one,-xmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1; 
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//1) c> edge1 with face 2...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1; 
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//1) d> edge1 with face 5...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//1) e> edge1 with face 3...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	//1) f> edge1 with face 6...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	/////////////////////////////////
	r.ox = (double)x1;
	r.oy = (double)y1;
	r.oz = (double)z1;
	r.dx = (double)(x2-x1);
	r.dy = (double)(y2-y1);
	r.dz = (double)(z2-z1);

	//2) a> edge2 with face 1 ...
	one = normalize(xmax-xmin,0,0);
	p = inbox(r,one,-xmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//2) b> edge2 with face 4...klc
	one = normalize(xmax- xmin,0,0);
	p = inbox(r,one,-xmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//2) c> edge2 with face 2...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//2) d> edge2 with face 5...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//2) e> edge2 with face 3...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	//2) f> edge2 with face 6...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	/////////////////////////////////
	r.ox = (double)x2;
	r.oy = (double)y2;
	r.oz = (double)z2;
	r.dx = (double)(x0-x2);
	r.dy = (double)(y0-y2);
	r.dz = (double)(z0-z2);

	//3)  a> edge3 with face 1 ...
	one = normalize(xmax- xmin,0,0);
	p = inbox(r,one,-xmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//3)  b> edge3 with face 4...
	one = normalize(xmax- xmin,0,0);
	p = inbox(r,one,-xmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, zmin, zmax, 1)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(ymin, ymax, p.y)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//3)  c> edge3 with face 2...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//3)  d> edge3 with face 5...
	one = normalize(0,ymax-ymin,0);
	p = inbox(r,one,-ymax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, xmin, xmax, zmin, zmax, 2)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(zmin, zmax, p.z))
				return 1;
	}

	//3)  e> edge3 with face 3...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmin, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	//3)  f> edge3 with face 6...
	one = normalize(0,0,zmax-zmin);
	p = inbox(r,one,-zmax, &t);
	if((p.isNull))
	{
		if ( (isEqual(1.0, t)) && (inside_cube(r, ymin, ymax, xmin, xmax, 3)) )		return 1;  
	}
	else
	{
		if(isBetween(0.0, 1.0, t))
			if ((isBetween(xmin, xmax, p.x)) && isBetween(ymin, ymax, p.y))
				return 1;
	}

	///////////////////////////
	//Then the case where the Cube intersects the Triangle...
	one.x = (double)x0;
	one.y = (double)y0;
	one.z = (double)z0;
	two.x = (double)x1;
	two.y = (double)y1;
	two.z = (double)z1;
	three.x = (double)x2;
	three.y = (double)y2;
	three.z = (double)z2;

	//2)1 a>.
	r.ox = (double)xmin; 
	r.oy = (double)ymin;
	r.oz = (double)zmin;
	r.dx = (double)0;
	r.dy = (double)0;
	r.dz = (double)(zmax-zmin);
 	if (ray_polygon_intersection(r, tri)) return 1;

	//2)1 b>.
	r.dx = (double)(xmax-xmin);
	r.dy = (double)0;
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)1 c>.
	r.dx = (double)0;
	r.dy = (double)(ymax-ymin);
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;


	//2)2 a>.
	r.ox = (double)xmax;
	r.oy = (double)ymax;
	r.oz = (double)zmax;
	r.dx = (double)0;
	r.dy = (double)0;
	r.dz = (double)(-(zmax-zmin));
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)2 b>.
	r.dx = (double)(-(xmax-xmin));
	r.dy = (double)0;
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)2 c>.
	r.dx = (double)0;
	r.dy = (double)(-(ymax-ymin));
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;


	//2)3 a>.
	r.ox = (double)xmax;
	r.oy = (double)ymax;
	r.oz = (double)zmin;
	r.dx = (double)0;
	r.dy = (double)(-(ymax-ymin));
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)3 c>.
	r.dx = (double)(-(xmax-xmin));
	r.dy = (double)0;
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;


	//2)4 a>.
	r.ox = (double)xmax;
	r.oy = (double)ymin;
	r.oz = (double)zmax;
	r.dx = (double)0;
	r.dy = (double)0;
	r.dz = (double)(-(zmax-zmin));
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)4 b>.
	r.dx = (double)(-(xmax-xmin));
	r.dy = (double)0;
	r.dz = (double)0;
	p.isNull = 0;
	if (ray_polygon_intersection(r, tri)) return 1;

	
	//2)5 a>.
	r.ox = (double)xmin;
	r.oy = (double)ymax;
	r.oz = (double)zmax;
	r.dx = (double)0;
	r.dy = (double)0;
	r.dz = (double)(-(zmax-zmin));
	if (ray_polygon_intersection(r, tri)) return 1;

	//2)5 b>.
	r.dx = (double)0;
	r.dy = (double)(-(ymax-ymin));
	r.dz = (double)0;
	if (ray_polygon_intersection(r, tri)) return 1;

	return 0;
}

void update_bounding_box(long int current_triangle, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int cur_level)
{
	int intersects = 0;
	int i, j, k;

	listnode* l;
	listnode* temp;

	intersects = 0;

	if( within(current_triangle, xmin, xmax, ymin, ymax, zmin, zmax ) ) intersects = 1;

	if( intersects )
	{
		if( cur_level < octree_depth )
		{
			update_bounding_box( current_triangle, xmin,					( xmax + xmin ) /2.0, 	( ymax + ymin ) /2.0, ymax,				zmin,				( zmax+ zmin ) /2.0,  cur_level+1 );
			update_bounding_box( current_triangle, ( xmax + xmin ) /2.0, 	xmax,					( ymax + ymin ) /2.0, ymax,				zmin,				( zmax+ zmin ) /2.0,  cur_level+1 );
			update_bounding_box( current_triangle, ( xmax + xmin ) /2.0, 	xmax,					( ymax + ymin ) /2.0, ymax,				( zmax + zmin ) /2.0,  zmax,				cur_level+1 );
			update_bounding_box( current_triangle, xmin,					( xmax + xmin ) /2.0, 	( ymax + ymin ) /2.0, ymax,				( zmax + zmin ) /2.0,  zmax,				cur_level+1 );
			update_bounding_box( current_triangle, xmin,					( xmax + xmin ) /2.0, 	ymin,				( ymax + ymin ) /2.0, zmin,				( zmax+ zmin ) /2.0,  cur_level+1 );
			update_bounding_box( current_triangle, ( xmax + xmin ) /2.0, 	xmax,					ymin,				( ymax + ymin ) /2.0, zmin,				( zmax+ zmin ) /2.0,  cur_level+1 );
			update_bounding_box( current_triangle, ( xmax + xmin ) /2.0, 	xmax,					ymin,				( ymax + ymin ) /2.0, ( zmax + zmin ) /2.0,  zmax,				cur_level+1 );
			update_bounding_box( current_triangle, xmin,					( xmax + xmin ) /2.0, 	ymin,				( ymax + ymin ) /2.0, ( zmax + zmin ) /2.0,  zmax,				cur_level+1 );
		}
		else
		{
			object2octree(xmin, ymin, zmin, xmax, ymax, zmax, i, j, k);
			
			l = (listnode*) malloc( sizeof( listnode ) );
			l->index = current_triangle;
			l->next = NULL;

			if( sdf[i][j][k].tindex == NULL )
			{
				sdf[i][j][k].useful = 1;
				sdf[i][j][k].tindex = l;
				sdf[i][j][k].no = 1;
				sdf[i][j][k].type =4;
			}
			else
			{
				temp = sdf[i][j][k].tindex;
				l->next = temp;
				sdf[i][j][k].tindex = l;
				sdf[i][j][k].no++;
			}

			update_boundary_vertices(i, j, k);
		}
	}
}

void write_octree()
{
	int i, j, k, sgn;

	FILE* fp = fopen("octree.txt","w");

	for(i=0; i<size; i++)
	{
		fprintf(fp, "%d\n", i);
		for(j=0; j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				sgn = sdf[i][j][k].useful;
				if (sgn== 0)
					fprintf(fp, "0");
				else if (sgn ==1)
					fprintf(fp, "1");
				else 
					fprintf(fp, "%d", sdf[i][j][k].useful);
			}
			fprintf(fp, "\n");	
		}
	}
	fflush(fp);
	fclose(fp);
	printf("octree.txt written \n");
}

int inside_cube(ray r, double ymin, double ymax, double zmin, double zmax, int flag)
{
	double k=1, t=1;
	double xmax, xmin;
	
	//intersect each edge of the tri with ray R n c if the pt of intersection is on the edge segment...
	switch (flag)
	{
	case 1: //ZY
		//a) ymin
		if (! isZero(r.dy))
		{
			t = (double)( (ymin - r.oy)/r.dy );
			k = r.oz + t*r.dz;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, k)))	return 1;
		}

		//b) ymax
		if (! isZero(r.dy))
		{
			t = (double)( (ymax - r.oy)/r.dy );
			k = r.oz + t*r.dz;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, k)))	return 1;
		}

		//c) zmin
		if (! isZero(r.dz))
		{
			t = (double)( (zmin - r.oz)/r.dz );
			k = r.oy + t*r.dy;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(ymin, ymax, k)))	return 1;
		}

		//d) zmax
		if (! isZero(r.dz))
		{
			t = (double)( (zmax- r.oz)/r.dz );
			k = r.oy + t*r.dy;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, k)))	return 1;
		}

		//if ((isBetween(zmin, zmax, r.oz)) && (ymin, ymax, r.oy)) return 1;
		break;

	case 2: //XZ
		xmin = ymin;	xmax = ymax;
		//a) xmin
		if (! isZero(r.dx))
		{
			t = (double)( (xmin - r.ox)/r.dx );
			k = r.oz + t*r.dz;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, k)))	return 1;
		}

		//b) xmax
		if (! isZero(r.dx))
		{
			t = (double)( (xmax - r.ox)/r.dx );
			k = r.oz + t*r.dz;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(zmin, zmax, k)))	return 1;
		}

		//c) zmin
		if (! isZero(r.dz))
		{
			t = (double)( (zmin - r.oz)/r.dz );
			k = r.ox + t*r.dx;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(xmin, xmax, k)))	return 1;
		}

		//d) zmax
		if (! isZero(r.dz))
		{
			t = (double)( (zmax- r.oz)/r.dz );
			k = r.ox + t*r.dx;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(xmin, xmax, k)))	return 1;		
		}

		//if ((isBetween(zmin, zmax, r.oz)) && (xmin, xmax, r.ox)) return 1;
		break;

	case 3: //YX
		xmin = zmin;	xmax = zmax;
		//a) ymin
		if (! isZero(r.dy))
		{
			t = (double)( (ymin - r.oy)/r.dy );
			k = r.ox + t*r.dx;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(xmin, xmax, k)))	return 1;
		}

		//b) ymax
		if (! isZero(r.dy))
		{
			t = (double)( (ymax - r.oy)/r.dy );
			k = r.ox + t*r.dx;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(xmin, xmax, k)))	return 1;
		}

		//c) xmin
		if (! isZero(r.dx))
		{
			t = (double)( (xmin - r.ox)/r.dx );
			k = r.oy + t*r.dy;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(ymin, ymax, k)))	return 1;
		}

		//d) xmax
		if (! isZero(r.dx))
		{
			t = (double)( (xmax- r.ox)/r.dx );
			k = r.oy + t*r.dy;
			if ((isBetween(0.0, 1.0, t)) && (isBetween(ymin, ymax, k)))	return 1;		
		}

		//if ((isBetween(ymin, ymax, r.oy)) && (xmin, xmax, r.ox)) return 1;
		break;

	default:
		printf("unknown case in inside_cube: %d \n", flag);
		return 1;
		break;
	}

	return 0;
}

myPoint inbox (ray r, myPoint p, double dist, double* t)
{
	myPoint result;
	double myt=0.0;

	myt = n_dotv(p.x, p.y, p.z, r, dist);
	result.x = result.y =result.z=0;

	if (myt==INF) //ie the denoms 0 and so the ray is || to the plane...
	{
		//c if the Origin of the ray satisfies the plane eqn or not.
		if ( isZero( (p.x * r.ox) + (p.y * r.oy) + (p.z * r.oz) + (double)dist) )
			*t=1;	//need to chq later if the ray ACTUALLY intersects the cube or not....
		else
			*t=0;

		result.isNull = 1;
		return result;
	}

	result.x = (double)(r.ox + (myt)*r.dx);
	result.y = (double)(r.oy + (myt)*r.dy);
	result.z = (double)(r.oz + (myt)*r.dz);
	result.isNull = 0;
	*t = myt;
	return result;
}

int chqOrientedCorrectly(myPoint* start, myPoint* finish, int tri, ray r)
{
	double dist[2];

	dist[0] = ( ((( start->x)* normals[tri].x +	(start->y)* normals[tri].y +
		(start->z)* normals[tri].z)) + distances[tri]);

	dist[1] = ( ((( finish->x)* normals[tri].x +	(finish->y)* normals[tri].y +
		(finish->z)* normals[tri].z)) + distances[tri]);

	if ( (isZero(dist[0])) || (isZero(dist[1])) )	return 1;
	if (dist[0]*dist[1] <0)	return 1;	
	return 0;
}

int sign3DTest(myPoint d, myPoint a, myPoint b, myPoint c)
{
	double m11, m12, m13, m21, m22, m23, m31, m32, m33;
	double determ;

	m11 = a.x-d.x;		m12 = a.y-d.y;		m13 = a.z-d.z;
	m21 = b.x-d.x;		m22 = b.y-d.y;		m23 = b.z-d.z;
	m31 = c.x-d.x;		m32 = c.y-d.y;		m33 = c.z-d.z;

	determ =  m11*(m22*m33-m23*m32);
	determ -= m12*(m21*m33-m23*m31);
	determ += m13*(m21*m32-m22*m31);
	determ /=6.0;

	if (isZero(determ))	return 0;
	if (isNegative(determ))	return -1;
	return 1;
}

int ray_polygon_intersection (ray r, int tri)
{
	myPoint start, finish, triA, triB, triC;
	int flag, i, j, k;

	start.x = r.ox;							start.y = r.oy;							start.z = r.oz;
	triA.x = vertices[surface[tri].v1].x;	triA.y = vertices[surface[tri].v1].y;	triA.z = vertices[surface[tri].v1].z;
	triB.x = vertices[surface[tri].v2].x;	triB.y = vertices[surface[tri].v2].y;	triB.z = vertices[surface[tri].v2].z;
	triC.x = vertices[surface[tri].v3].x;	triC.y = vertices[surface[tri].v3].y;	triC.z = vertices[surface[tri].v3].z;
	
	//Start from the Origin and shoot to the bounding box.
	if (r.dx >0)							finish.x = (double)(size+1);
	else if (r.dx ==0)						finish.x = (double)(start.x);
	else									finish.x = (double)(0.0);

	if (r.dy >0)							finish.y = (double)(size+1);
	else if (r.dy ==0)						finish.y = (double)(start.y);
	else									finish.y = (double)(0.0);

	if (r.dz >0)							finish.z = (double)(size+1);
	else if (r.dz ==0)						finish.z = (double)(start.z);
	else									finish.z = (double)(0.0);

	if (chqOrientedCorrectly(&start, &finish, tri, r) ==0)
	{
		//Both the points are on the same side of the trianngle. 
		//So, the ray will NEVER intersect the triangle ?
		return 0;
	}

	if (sign3DTest(start, triA, triB, triC) >= 0)
	{
		flag =1;
		i = sign3DTest(finish, triA, triB, start);
		j = sign3DTest(finish, triB, triC, start);
		k = sign3DTest(finish, triC, triA, start);
	}
	else if (sign3DTest(start, triA, triC, triB) >= 0)
	{
		flag =-1;
		i = sign3DTest(finish, triA, triC, start);
		j = sign3DTest(finish, triC, triB, start);
		k = sign3DTest(finish, triB, triA, start);
	}
	else
		printf("wot now?\n");

	//First, intersection at a vertex of the triangle
	/*if ( ( (i==0) && (j==0) )  || 
		 ( (i==0) && (k==0) )  ||
		 ( (j==0) && (k==0) )  )
	{
		//intersects in the corresponding vertex.
		return 3;
	}*/
	if ( (j==0) && (k==0) )				return (flag*31);
	else if ( (i==0) && (k==0) )		return (flag*32);
	else if ( (j==0) && (i==0) )		return (flag*33);

	//Second, interesection on an edge of the triangle
	/*if ( (i==0) && (j==k) ||
		 (j==0) && (i==k) ||
		 (k==0) && (i==j) )
	{
		//intersects in the corresponding edge
		return 2;
	}*/
	if ( (i==0) && (j==k) )				return (flag*21);
	else if ( (j==0) && (i==k) )		return (flag*22);
	else if ( (k==0) && (j==i) )		return (flag*23);

	//Finally, intersection inside the triangle
	if ( (i==j) && (j==k) )
		return 1;

	return 0;
}

double n_dotv( double x, double y, double z, ray r, double temp123)
{
	double ndotv1, ndotv2, dote;
	ndotv1=ndotv2=dote=0.0;

	ndotv1 += x * r.dx;
	ndotv1 += y * r.dy;
	ndotv1 += z * r.dz;
	if (isZero(ndotv1))	 return INF;

	ndotv2 += x * r.ox;
	ndotv2 += y * r.oy;
	ndotv2 += z * r.oz;

	dote = -(ndotv2+ temp123)/(ndotv1);
	return dote;
}

myPoint normalize(double x, double y, double z)
{
	myPoint result;
	double n;

	n = sqrt(x*x  + y*y + z*z);
	result.x = (double)(x/n);
	result.y = (double)(y/n);
	result.z = (double)(z/n);

	return result;
}

int max_3( double x, double y, double z )
{
	if( x < 0 ) x *= -1;
	if( y < 0 ) y *= -1;
	if( z < 0 ) z *= -1;

	if( x > y )
	{
		if( x > z ) return 0;
		return 2;
	}
	if( y > z ) return 1;
	return 2;
}

////////////////////////////////////.

//take a point on the plane and c if it lies in the triangle or not.
int point_in_polygon(myPoint result, int tri)
{
	double alpha, beta;
	double u0, u1, u2, v0, v1, v2;
	int index;
	double p1, p2;
	int i, j;

	//added this
	if ( !isZero((result.x*normals[tri].x) + (result.y*normals[tri].y) + 
			(result.z*normals[tri].z) + (distances[tri])) )
		return 0;

	//now do the point in Triangle test...
	index = max_3( normals[tri].x, normals[tri].y, normals[tri].z );
	if( index == 0 )
	{
		p1 = result.y; p2 = result.z;
		i = 1; j = 2;
		u0 = p1 - vertices[ surface[tri].v1 ].y;
		u1 = vertices[ surface[tri].v2 ].y - vertices[ surface[tri].v1 ].y;
		u2 = vertices[ surface[tri].v3 ].y - vertices[ surface[tri].v1 ].y;

		v0 = p2 - vertices[ surface[tri].v1 ].z;
		v1 = vertices[ surface[tri].v2 ].z - vertices[ surface[tri].v1 ].z;
		v2 = vertices[ surface[tri].v3 ].z - vertices[ surface[tri].v1 ].z;
	}
	else if ( index == 1 )
	{
		p1 = result.z; p2 = result.x;
		i = 2; j = 0;
		u0 = p1 - vertices[ surface[tri].v1 ].z;
		u1 = vertices[ surface[tri].v2 ].z - vertices[ surface[tri].v1 ].z;
		u2 = vertices[ surface[tri].v3 ].z - vertices[ surface[tri].v1 ].z;

		v0 = p2 - vertices[ surface[tri].v1 ].x;
		v1 = vertices[ surface[tri].v2 ].x - vertices[ surface[tri].v1 ].x;
		v2 = vertices[ surface[tri].v3 ].x - vertices[ surface[tri].v1 ].x;
	}
	else
	{
		p1 = result.x; p2 = result.y;
		i = 0; j = 1;
		u0 = p1 - vertices[ surface[tri].v1 ].x;
		u1 = vertices[ surface[tri].v2 ].x - vertices[ surface[tri].v1 ].x;
		u2 = vertices[ surface[tri].v3 ].x - vertices[ surface[tri].v1 ].x;

		v0 = p2 - vertices[ surface[tri].v1 ].y;
		v1 = vertices[ surface[tri].v2 ].y - vertices[ surface[tri].v1 ].y;
		v2 = vertices[ surface[tri].v3 ].y - vertices[ surface[tri].v1 ].y;
	}

	alpha = ( u0*v2 - v0*u2 ) / ( u1*v2 - v1*u2 );
	if( isNegative(alpha) ) return 0;

	beta = ( u1*v0 - v1*u0 ) / ( u1*v2 - v1*u2 );
	if( isNegative(beta)) return 0;

	if( isBetween(0.0, 1.0, alpha+beta))
		return 1;

	return 0;
}


//Add all the vertices of the cell into the boundary_vertices array.
void update_boundary_vertices(int cx, int cy, int cz)
{
	insert_bound_vert(index2vert(cx+0, cy+0, cz+0));
	insert_bound_vert(index2vert(cx+1, cy+0, cz+0));
	insert_bound_vert(index2vert(cx+1, cy+1, cz+0));
	insert_bound_vert(index2vert(cx+0, cy+1, cz+0));

	insert_bound_vert(index2vert(cx+0, cy+0, cz+1));
	insert_bound_vert(index2vert(cx+1, cy+0, cz+1));
	insert_bound_vert(index2vert(cx+1, cy+1, cz+1));
	insert_bound_vert(index2vert(cx+0, cy+1, cz+1));
}