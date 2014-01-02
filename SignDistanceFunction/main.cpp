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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "sdfLib.h"

using namespace SDFLibrary;

void free_memory()
{
	int i, j, k;
	SDFLibrary::listnode* temp;
	SDFLibrary::listnode* currNode;

	printf("starting memory de-allocation\n");

	//1. Octree
	for (i = 0; i < SDFLibrary::size; i++)
	{
		for (j = 0; j < SDFLibrary::size; j++)
		{
			for (k = 0; k < SDFLibrary::size; k++)
			{
				currNode = SDFLibrary::sdf[i][j][k].tindex;

				while(currNode != NULL)
				{
					temp = currNode;
					currNode = currNode->next;
					free(temp);
				}
			}
			free(SDFLibrary::sdf[i][j]);
		}
		free(SDFLibrary::sdf[i]);
	}	
	free(SDFLibrary::sdf);

	free(SDFLibrary::values);

	if (SDFLibrary::vertices != NULL)
		free(SDFLibrary::vertices);

	if (SDFLibrary::surface != NULL)
		free(SDFLibrary::surface);

	if (SDFLibrary::normals != NULL)
		free(SDFLibrary::normals);

	if (SDFLibrary::distances != NULL)
		free(SDFLibrary::distances);

	if (SDFLibrary::queues != NULL)
		free(SDFLibrary::queues);

	if (SDFLibrary::bverts != NULL)
		free(SDFLibrary::bverts);

	printf("Memory de-allocated successfully! \n");
}

void SDFLibrary::setParameters(int Size, int isNormalFlip, float* mins, float* maxs)
{
	//First the default values.
	SDFLibrary::init_all_vars();
	
	//Then, assign the actual input values.
	SDFLibrary::size = Size;
	SDFLibrary::flipNormals = isNormalFlip;

	SDFLibrary::minext[0] = mins[0];	SDFLibrary::minext[1] = mins[1];	SDFLibrary::minext[2] = mins[2];
	SDFLibrary::maxext[0] = maxs[0];	SDFLibrary::maxext[1] = maxs[1];	SDFLibrary::maxext[2] = maxs[2];
	SDFLibrary::span[0] = (maxs[0]-mins[0])/(SDFLibrary::size);
	SDFLibrary::span[1] = (maxs[1]-mins[1])/(SDFLibrary::size);
	SDFLibrary::span[2] = (maxs[2]-mins[2])/(SDFLibrary::size);

	if ((Size!=16) && (Size!=32) &&(Size!=64) && (Size!=128) && (Size!=256) &&(Size!=512) &&(Size!=1024))
	{
		printf("size is incorrect\n");
		exit(1);
	}
}

float* SDFLibrary::computeSDF(int nverts, float* verts, int ntris, int* tris)
{
	int i, numb;
	float* sdfValues =NULL;
	float isoval;

	//Set up the volume grid
	if( !initSDF() ) return 0;

	//Read in the Geometry
	readGeom(nverts, verts, ntris, tris);

	//Setup the Octree
	adjustData();

	//Compute the SDF
    compute();

	//Return the SDF
	numb = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
	sdfValues = (float*)(malloc(sizeof(float)*(numb)));
	isoval = 100.0f;

	for (i=0; i<numb; i++)
		sdfValues[i] = SDFLibrary::values[i].value * SDFLibrary::values[i].signe;

	free_memory();

	return (sdfValues);
}	

RAWIV_header* SDFLibrary::getVolumeInfo()
{
	int i;

	RAWIV_header* ret = (RAWIV_header*)(malloc(sizeof(RAWIV_header)*1));

	for (i=0; i<3; i++)
	{
		ret->minext[i] = SDFLibrary::minext[i];
		ret->maxext[i] = SDFLibrary::maxext[i];
		ret->span[i] = SDFLibrary::span[i];
		ret->origin[i] = 0.0f;
		ret->dim[i] = SDFLibrary::size+1;
	}

	ret->ngridpts = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
	ret->ncells = (SDFLibrary::size)*(SDFLibrary::size)*(SDFLibrary::size);
	ret->size = SDFLibrary::size;

	return ret;
}
