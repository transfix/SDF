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


#include <math.h>

#include "common.h"

using namespace SDFLibrary;

extern double point_2_plane(int tri, double i, double j, double k, myPoint* inter);
extern void compute_SDF(int vi, int vj, int vk, int ci, int cj, int ck);

double dist_grid_3Dpts(int one, int two)
{
	int vi, vj, vk;
	myPoint temp;

	_vert2index(one, vi, vj, vk);
	return (point_2_plane(values[two].closestV, xCoord(vi), yCoord(vj), zCoord(vk), &temp));
}

void update_distance_2_vertex(int ind, int vi, int vj, int vk)
{
	double val;
	int upd;

	upd = index2vert(vi, vj, vk);
	
	if ( (vi>=0) && (vi<=size) )
	{
		if ( (vj>=0) && (vj<=size) )
		{
			if ( (vk>=0) && (vk<=size) )
			{
				if (values[upd].processed == 1)	return;
		
				val = dist_grid_3Dpts(upd, ind);

				if (val < values[upd].value )
				{
					values[upd].value = (float)val;
					values[upd].closestV = values[ind].closestV;
				}
				insert_bound_vert(upd);
			}	
		}	
	}
}


//Current implementation only does the COMPLETE 3X3 Distance Matrix
void SDFLibrary::apply_distance_transform(int vi, int vj, int vk)
{
	int ind;  //Current vertex

	ind = index2vert(vi, vj, vk);

	//Front Y slice
	update_distance_2_vertex(ind, vi-1, vj-1, vk-1);
	update_distance_2_vertex(ind, vi,   vj-1, vk-1);
	update_distance_2_vertex(ind, vi+1, vj-1, vk-1);
	
	update_distance_2_vertex(ind, vi-1, vj-1, vk);
	update_distance_2_vertex(ind, vi,   vj-1, vk);
	update_distance_2_vertex(ind, vi+1, vj-1, vk);

	update_distance_2_vertex(ind, vi-1, vj-1, vk+1);
	update_distance_2_vertex(ind, vi,   vj-1, vk+1);
	update_distance_2_vertex(ind, vi+1, vj-1, vk+1);

	//Middle Y slice
	update_distance_2_vertex(ind, vi-1, vj,	 vk-1);
	update_distance_2_vertex(ind, vi,   vj,	 vk-1);
	update_distance_2_vertex(ind, vi+1, vj,	 vk-1);
	
	update_distance_2_vertex(ind, vi-1, vj,	 vk);
  //update_distance_2_vertex(ind, vi,   vj,	 vk); //Current vertex
	update_distance_2_vertex(ind, vi+1, vj,	 vk);

	update_distance_2_vertex(ind, vi-1, vj,	 vk+1);
	update_distance_2_vertex(ind, vi,   vj,	 vk+1);
	update_distance_2_vertex(ind, vi+1, vj,	 vk+1);

	//Back Y slice
	update_distance_2_vertex(ind, vi-1, vj+1, vk-1);
	update_distance_2_vertex(ind, vi,   vj+1, vk-1);
	update_distance_2_vertex(ind, vi+1, vj+1, vk-1);
	
	update_distance_2_vertex(ind, vi-1, vj+1, vk);
	update_distance_2_vertex(ind, vi,   vj+1, vk);
	update_distance_2_vertex(ind, vi+1, vj+1, vk);

	update_distance_2_vertex(ind, vi-1, vj+1, vk+1);
	update_distance_2_vertex(ind, vi,   vj+1, vk+1);
	update_distance_2_vertex(ind, vi+1, vj+1, vk+1);
}

void SDFLibrary::insert_bound_vert(int vert)
{
	if(bverts[vert] == 0) //ie not found
	{
		bverts[vert] =1;
		queues[all_verts_touched++] = vert;
	}
}

void propagate_from_here(int vert)
{
	int i, j, k, ci, cj, ck, level, ind, test;
	int MAX_LEVELS =10;

	_vert2index(vert, i, j, k);

	if (i == size) ci = i-1;	else	ci=i;
	if (j == size) cj = j-1;	else	cj=j;
	if (k == size) ck = k-1;	else	ck=k;

	for (level=1; level<MAX_LEVELS; level++)
	{
		for (ci=i-level; ci<=i+level; ci++)
		{
			for (cj=j-level; cj<=j+level; cj++)
			{
				for (ck=k-level; ck<=k+level; ck++)
				{
					if ((ci < 0) || (ci >= size))  continue;
					if ((cj < 0) || (cj >= size))  continue;
					if ((ck < 0) || (ck >= size))  continue;

					ind = index2vert(ci, cj, ck);
					test = values[ind].processed;
					test = values[ind].signe;
					test = (int)values[ind].value;
					if ( (values[ind].processed ==1) && (values[ind].value != MAX_DIST) )
						update_distance_2_vertex(ind, i, j, k);
				}
			}
		}
	}
}

int confirm_SDF(int flag)
{
	int i, grid_pts;

	grid_pts = (size+1)*(size+1)*(size+1);

	for (i=0; i<grid_pts; i++)
	{
		if (isEqual(values[i].value, MAX_DIST))
		{
			if (flag) 
				printf("some error in computing the SDF for vertex %d\n", i);
			propagate_from_here(i);
		}
	}

	return 0;
}