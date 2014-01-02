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

extern int point_in_polygon(myPoint result, int tri);
extern int ray_polygon_intersection (ray r, int tri);
extern int confirm_SDF(int flag);

double point_2_plane(int tri, double i, double j, double k, myPoint* inter);
double getClipPoint(int one, int two, double i, double j, double k, myPoint* inter);

double getEuclideanDist(int on, double vi, double vj, double vk)
{
	myPoint one;

	one.x = vertices[on].x;		one.y = vertices[on].y;		one.z = vertices[on].z;
	return (sqrt( (one.x-vi)*(one.x-vi) + (one.y-vj)*(one.y-vj) + (one.z-vk)*(one.z-vk) ));
}

/*In case you intersect an edge (-21, -22, -23 or 21, 22, 23), then, find the distance of the grid_vertex to that edge
  In case you intersect a vertex, (-31, -32, -33 or 31, 32 33), then, find the distance between the two points*/
int chq_inters(int pts[100], int inters, double  vi, double vj, double vk)
{
	int i, j, counter, ret;
	myPoint interPoint;
	double distns[500];

	for (counter=0,i=1; i<2*inters;)
	{
		switch(pts[i])
		{
		case 1: break;
		case 21: case -23:
			distns[counter] = getClipPoint(surface[pts[i-1]].v1, surface[pts[i-1]].v2, vi, vj, vk, &interPoint);
			counter++;
			break;
		case 22: case -22:
			distns[counter] = getClipPoint(surface[pts[i-1]].v2, surface[pts[i-1]].v3, vi, vj, vk, &interPoint);
			counter++;
			break;
		case 23: case -21:
			distns[counter] = getClipPoint(surface[pts[i-1]].v3, surface[pts[i-1]].v1, vi, vj, vk, &interPoint);
			counter++;
			break;
		case 31: case -33:
			distns[counter] = getEuclideanDist(surface[pts[i-1]].v3, vi, vj, vk);
			counter++;
			break;
		case 32: case -32:
			distns[counter] = getEuclideanDist(surface[pts[i-1]].v1, vi, vj, vk);
			counter++;
			break;
		case 33: case -31:
			distns[counter] = getEuclideanDist(surface[pts[i-1]].v2, vi, vj, vk);
			counter++;
			break;
		}
		i +=2;
	}
	
	for (ret=0,i=0; i<counter; i++)
	{
		for (j=i+1; j<counter; j++)
			if (isEqual(distns[i], distns[j]) )
				ret++;
	}

	if (counter ==1)	ret=0;

	return (inters+ret);
}

int x_assign(int vi, int vj, int vk)
{
	int i, j, k, inters, temp, flag, v, ret;
	listnode* currNode;
	ray r;
	int pts[1000];	//the given ray cant intersect the surface more than 1000 times...

	//pick a direction n shoot rays to the outside BB.
	i=(int)vi;				j= (int)vj;				k= (int)vk;
	r.ox = xCoord(vi);		r.oy = yCoord(vj);		r.oz = zCoord(vk);
	r.dx =1.0f;				r.dy =0.0f;				r.dz =0.0f;
 	inters=0;				v = i;

	for (i=v; i<size; i++)
	{
		if (sdf[i][j][k].type==4)
		{
			currNode = sdf[i][j][k].tindex;

			while(currNode != NULL)
			{
				ret = ray_polygon_intersection(r, currNode->index);

				if (ret !=0)
				{
					for (flag =0,temp=0; temp<2*inters;)
					{
						if (pts[temp] == currNode->index)	flag=1;
						temp +=2;
					}

					if (flag ==0)
					{
						pts[2*inters+0] = currNode->index;
						pts[2*inters+1] = ret;
						inters++;
					}
				}
				currNode = currNode->next;
			}
		}
	}

	return (chq_inters(pts, inters, r.ox, r.oy, r.oz));
}

int y_assign(int vi, int vj, int vk)
{
	int i, j, k, inters, temp, flag, v, ret;
	listnode* currNode;
	ray r;
	int pts[1000];	//the given ray cant intersect the surface more than 1000 times...

	//pick a direction n shoot rays to the outside BB.
	i=(int)vi;				j= (int)vj;				k= (int)vk;
	r.ox = xCoord(vi);		r.oy = yCoord(vj);		r.oz = zCoord(vk);
	r.dx =0.0f;				r.dy =1.0f;				r.dz =0.0f;
	inters=0;				v = j;

	for (j=v; j<size; j++)
	{
		if (sdf[i][j][k].type==4)
		{
			currNode = sdf[i][j][k].tindex;

			while(currNode != NULL)
			{
				ret = ray_polygon_intersection(r, currNode->index);

				if (ret !=0)
				{
					for (flag =0,temp=0; temp<2*inters;)
					{
						if (pts[temp] == currNode->index)	flag=1;
						temp +=2;
					}

					if (flag ==0)
					{
						pts[2*inters+0] = currNode->index;
						pts[2*inters+1] = ret;
						inters++;
					}
				}
				currNode = currNode->next;
			}
		}
	}

	return (chq_inters(pts, inters, r.ox, r.oy, r.oz));
}

int z_assign(int vi, int vj, int vk)
{
	int i, j, k, inters, temp, flag, v, ret;
	listnode* currNode;
	ray r;
	int pts[1000];	//the given ray cant intersect the surface more than 1000 times...

	//pick a direction n shoot rays to the outside BB.
	i=(int)vi;				j= (int)vj;				k= (int)vk;
	r.ox = xCoord(vi);		r.oy = yCoord(vj);		r.oz = zCoord(vk);
	r.dx =0.0f;				r.dy =0.0f;				r.dz =1.0f;
	inters=0;				v = k;

	for (k=v; k<size; k++)
	{
		if (sdf[i][j][k].type==4)
		{
			currNode = sdf[i][j][k].tindex;

			while(currNode != NULL)
			{
				ret = ray_polygon_intersection(r, currNode->index);

				if (ret !=0)
				{
					for (flag =0,temp=0; temp<2*inters;)
					{
						if (pts[temp] == currNode->index)	flag=1;
						temp +=2;
					}

					if (flag ==0)
					{
						pts[2*inters+0] = currNode->index;
						pts[2*inters+1] = ret;
						inters++;
					}
				}
				currNode = currNode->next;
			}
		}
	}

	return (chq_inters(pts, inters, r.ox, r.oy, r.oz));
}

int klc_assign(int vi, int vj, int vk)
{
	int inters[3];

	if ((vi<=0) || (vj<=0) || (vk<=0) || (vi>=size) || (vj>=size) || (vk>=size))
	{
		return 1;
	}
	else
	{
		inters[0] = x_assign(vi, vj, vk);
		inters[1] = y_assign(vi, vj, vk);
		inters[2] = z_assign(vi, vj, vk);
	}

	if ( (inters[0]%2 ==0) && (inters[1]%2 ==0) && (inters[2]%2 ==0) )
		return 1;
	else if ( (inters[0]%2 ==1) && (inters[1]%2 ==1) && (inters[2]%2 ==1) )
		return -1;

	//Else, u have run into some doubleing pt. error and need to count the #of intersections.
	if ( ((inters[0]%2) + (inters[1]%2) + (inters[2]%2)) %2 ==1)
		return 1;
	else
		return -1;
}

double sort_3_distances(double vals[3], myPoint closest[3], myPoint* inter)
{
	double dist;

	if (vals[0] <= vals[1])
	{
		if (vals[0] <=vals[2])
		{
			dist = vals[0];
			inter->x = closest[0].x;		inter->y = closest[0].y;		inter->z = closest[0].z;
		}
		else
		{
			dist = vals[2];
			inter->x = closest[2].x;		inter->y = closest[2].y;		inter->z = closest[2].z;
		}
	}
	else
	{
		if (vals[1] <=vals[2])
		{
			dist = vals[1];
			inter->x = closest[1].x;		inter->y = closest[1].y;		inter->z = closest[1].z;
		}
		else
		{
			dist = vals[2];
			inter->x = closest[2].x;		inter->y = closest[2].y;		inter->z = closest[2].z;
		}
	}
	
	return dist;
}

//project the point proj onto the line i-j and find the nearest point.
double getClipPoint(int one, int two, double i, double j, double k, myPoint* inter)
{
	double denom, theta, tempLen, t;
	double d1, d2, d3, e1, e2, e3;

	//1) Normalize the triangle edge and the line from a vertex to the point.
	d1 = vertices[one].x - vertices[two].x;
	d2 = vertices[one].y - vertices[two].y;
	d3 = vertices[one].z - vertices[two].z;
	denom = d1*d1 + d2*d2 + d3*d3;
    denom = sqrt(denom);
    tempLen = denom; //len of the tri edge
	d1 /=denom;	d2 /=denom;	d3 /=denom;

	e1 = i - vertices[two].x;
	e2 = j - vertices[two].y;
	e3 = k - vertices[two].z;
	denom = e1*e1 + e2*e2 + e3*e3;
	if (isZero(denom))
	{
		//Then, the pt is CLOSE to the vertex TWO.
		inter->x = vertices[two].x;		inter->y = vertices[two].y;		inter->z = vertices[two].z;
		return fabs(denom);
	}
    denom = sqrt(denom); //len of the line from point to tri vertex.
	e1 /=denom;	e2 /=denom;	e3 /=denom;

	//2) Find the angle between these lines.
	theta = e1*d1 + e2*d2 + e3*d3;

	if (isZero(theta))
	{
		d1 = i - vertices[one].x;
		d2 = j - vertices[one].y;
		d3 = k - vertices[one].z;
		theta = sqrt(d1*d1 + d2*d2 + d3*d3);
		
		if (theta <=denom)
		{
			inter->x = vertices[one].x;		inter->y = vertices[one].y;		inter->z = vertices[one].z;	return fabs(theta);
		}
		else
		{
			inter->x = vertices[two].x;		inter->y = vertices[two].y;		inter->z = vertices[two].z;	return fabs(denom);
		}		
	}

    if (theta <0)
	{
		inter->x = vertices[two].x;		inter->y = vertices[two].y;		inter->z = vertices[two].z;
		return (denom);
	}
    else if ((denom * theta) > tempLen)
    {
        d1 = i - vertices[one].x;
		d2 = j - vertices[one].y;
		d3 = k - vertices[one].z;
		theta = d1*d1 + d2*d2 + d3*d3;
		inter->x = vertices[one].x;		inter->y = vertices[one].y;		inter->z = vertices[one].z;	
        return sqrt(theta);
    }
    else
    {
        t = denom*theta;				theta = acos(theta);		
		inter->x = vertices[two].x + t*(vertices[one].x - vertices[two].x);
		inter->y = vertices[two].y + t*(vertices[one].y - vertices[two].y);
		inter->z = vertices[two].z + t*(vertices[one].z - vertices[two].z);
        return fabs((sin(theta) *denom));
    }
}


//compute the least distance of the triangle TRI to the vertex (i,j,k).
double point_2_plane(int tri, double i, double j, double k, myPoint* inter)
{
	double  dist, temp[3];
	int fact;
	myPoint res, val[3];

	//1) First compute the shortest signed distance between the Vertex and the Plane of the Triangle.
	dist = ( ((i* normals[tri].x +	j* normals[tri].y + k* normals[tri].z)) + distances[tri]);

	if (isZero(dist))
	{
		res.x =i;	res.y =j;	res.z =k;
		if(point_in_polygon(res, tri))
		{
			(*inter).x = res.x;		(*inter).y = res.y;			(*inter).z = res.z;
			return fabs(dist);
		}
	}
	if (dist <0) fact = -1;
	else fact = 1;

	//2) Then chq if the projected point is within the triangle or not. if yes, then the above is the correct shortest distance.
	res.x = (double)(i -normals[tri].x*dist);
	res.y = (double)(j -normals[tri].y*dist);
	res.z = (double)(k -normals[tri].z*dist);
	if (point_in_polygon(res, tri))
	{
		(*inter).x = res.x;		(*inter).y = res.y;			(*inter).z = res.z;
		return fabs(dist);
	}

    //3) now, project the planeProj onto the edge "i <-> j" and get the nearest point on the line segment to this point.
	temp[0] = getClipPoint(surface[tri].v1, surface[tri].v2, i, j, k, &val[0]);
    temp[1] = getClipPoint(surface[tri].v3, surface[tri].v2, i, j, k, &val[1]);
    temp[2] = getClipPoint(surface[tri].v1, surface[tri].v3, i, j, k, &val[2]);

    dist = sort_3_distances(temp, val, inter);

	//if (dist >= (MAX_DIST) || dist <= (-1*MAX_DIST))
	//	printf("ERR: vert= <%lf %lf %lf> tri= <%d> dist is %lf, max allowed = %lf\n", i, j, k, tri, dist, MAX_DIST);
	return dist;
}


//compute the distance of each Triangle in the cell (ci,cj,ck) from the Vertex vert.
int each_cell(int ci, int cj, int ck, int vi, int vj, int vk)
{
	int vert, cellnum, currrentTri;
	listnode* currNode;
	double val;
	myPoint temp;
	int res =0;

	vert = index2vert(vi, vj, vk);
	currNode = sdf[ci][cj][ck].tindex;
	cellnum = index2cell(ci, cj, ck);

	while (currNode != NULL)
	{
		currrentTri = currNode->index;
		val = (double)point_2_plane(currrentTri, xCoord(vi), yCoord(vj), zCoord(vk), &temp);
		if (val < values[vert].value )	
		{
			values[vert].value = (float)val;
			values[vert].closestV = currrentTri;
		}

		currNode = currNode->next;
		res =1;
	}

	//values[vert].processed =1;

	if (values[vert].value >= (MAX_DIST) || values[vert].value <= (-1*MAX_DIST))
		printf("err vert= %d %d %d \n", vi, vj, vk);

	return res;
}

//Just compute the Distance Function for the given vert.
void compute_SDF(int i, int j, int k)
{
	int level, ci, cj, ck;

	//Get the corresponding Octree cell.
	if (i == size) ci = i-1;	else	ci=i;
	if (j == size) cj = j-1;	else	cj=j;
	if (k == size) ck = k-1;	else	ck=k;

	level =2;

	for (ci=i-level; ci<=i+level; ci++)
	{
		for (cj=j-level; cj<=j+level; cj++)
		{
			for (ck=k-level; ck<=k+level; ck++)
			{
				if ((ci < 0) || (ci >= size))  continue;
				if ((cj < 0) || (cj >= size))  continue;
				if ((ck < 0) || (ck >= size))  continue;

				if (sdf[ci][cj][ck].useful >0)
					each_cell(ci, cj, ck, i, j, k);
			}
		}
	}
}


void compute_boundarySDF()
{
	int ind, vi, vj, vk;

	for (ind =0; ind<all_verts_touched; ind++)
	{
		_vert2index(queues[ind], vi, vj, vk);
		compute_SDF(vi, vj, vk);
		values[queues[ind]].processed =1;
		if (ind%5000 ==0) printf("computing the boundary SDF %3.4f %% over\n", (100.0*ind/all_verts_touched));
	}
}

void compute_signs()
{
	int i, j, k, sgn, ind;

	printf("\nnow going to compute.\n");
//	FILE* fp = fopen("chq.txt","w");

	for(i=0; i<=size; i++)
	{
//		fprintf(fp, "%d\n", i);
		for(j=0; j<=size; j++)
		{
			for(k=0; k<=size; k++)
			{
				sgn= klc_assign(i, j, k);
				ind = index2vert(i, j, k);
				values[ind].signe = sgn;
/*
				if (sgn== -1)
					fprintf(fp, "0");
				else
					fprintf(fp, "1");
*/
			}
//			fprintf(fp, "\n");	
			//printf("SDF computations: %d %d %d over\n", i, j, k);
		}

		fprintf(stderr,"SDF computations: %5.2f %%\r",(float(i)/float(size))*100.0);
	}
//	fclose(fp);
}


void SDFLibrary::compute()
{
	int grid_pts;
	int i, j, k, m, ind, prevMin, prevVerts;
	double t1, t2 ,t;

	//First compute the SIGN for all the volume grid points.
	t1 = getTime();
	compute_signs();
	t2 = getTime();		t = t2-t1;
	printf("Sign computations done in %f seconds\n", (t2-t1));
	
	//Next, compute the SDF for the boundary voxels.
	t1 = getTime();
	compute_boundarySDF();	
	t2 = getTime();		t += t2-t1;
	printf("Function evaluated at the %d boundary vertices in %f seconds\n",all_verts_touched,  (t2-t1));

	//Finally, propagate the distances to the interior voxels.
	grid_pts = (size+1)*(size+1)*(size+1);
	printf("total grid points: %d and starting with %d points\n", grid_pts, all_verts_touched);

	m=0;	prevMin =0;		prevVerts = all_verts_touched;		t1 = getTime();

	do {
		//FOR EACH VOXEL IN THE BOUNDARY_VERTS
		for (ind=prevMin; ind<prevVerts; ind++)
		{
			_vert2index(queues[ind], i, j, k);

			if ( (prevMin !=0) && (values[queues[ind]].processed == 1)	)	continue;

			apply_distance_transform(i, j, k);

			values[queues[ind]].processed =1;
			if (ind%10000 ==0)	printf("iter#%d: %d processed\n", m, ind);
		}	
		
		prevMin = prevVerts;
		prevVerts = all_verts_touched;
		m++;
		printf("in Iteration# %d, with %d vertices in the queue\n", m, prevVerts);

		if (prevMin == prevVerts)
		{
			printf("SDF propagation saturated. Now, checking for untouched voxels... \n");
			confirm_SDF(0);
			break;
		}

	} while(all_verts_touched != grid_pts);

	t2 = getTime();		t += t2-t1;
	printf("Distance Propagation for %d grid points done in %f seconds\n", all_verts_touched, (t2-t1));
	printf("All of the SDF computations are done in %f seconds!!! \n", t);

	confirm_SDF(1);
}


//The brute force method of computing the SDF. Compute the distance of a grid point from each triangle
/*
void compute()
{
	int i, j, k, t, trid;
	double dist, val;
	myPoint temp;
	double t1, t2, tt;

  	t1 = getTime();
	compute_signs();
	t2 = getTime();		tt = t2-t1;
	printf("Sign computations done in %f seconds\n", (t2-t1));

    t1 = getTime();
	for(i=0; i<=size; i++)
	{
		for(j=0; j<=size; j++)
		{
			for(k=0; k<=size; k++)
			{
				for (dist=1000,t=0; t<total_triangles; t++)
				{
					val = point_2_plane(t, i, j, k, &temp);
					if (val<dist)	
					{
						dist = val;
						trid = t;
					}
				}

				values[index2vert(i, j, k)].value = (float)dist;
				values[index2vert(i, j, k)].closestV = trid;
			}
			printf("SDF computations: %d %d %d over\n", i, j, k);
		}
	}

	t2 = getTime();		tt += t2-t1;
	printf("Distance Field computations done in %f seconds\n", tt);
}
*/
