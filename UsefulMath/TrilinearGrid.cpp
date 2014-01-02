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
// TrilinearGrid.cpp: implementation of the TrilinearGrid class.
//
//////////////////////////////////////////////////////////////////////

#include "TrilinearGrid.h"

using CCVOpenGLMath::TrilinearGrid;

TrilinearGrid::TrilinearGrid()
{

}

TrilinearGrid::~TrilinearGrid()
{

}

void TrilinearGrid::cell2xyz(int cell_index, int& x, int& y, int& z, unsigned int* gdim)
{
	x = cell_index % (gdim[0]-1);
	y = (cell_index / (gdim[0]-1)) % (gdim[1]-1);
	z = cell_index / ( (gdim[0]-1) * (gdim[1]-1) );
}

int TrilinearGrid::xyz2cell(int x, int y, int z, unsigned int* gdim)
{
	return x+y*(gdim[0]-1)+z*(gdim[0]-1)*(gdim[1]-1);
}

void TrilinearGrid::vtx2xyz(int vtx_idx,int& x, int& y, int& z, unsigned int* gdim)
{
	x = vtx_idx % (gdim[0]);
	y = (vtx_idx / (gdim[0])) % (gdim[1]);
	z = vtx_idx / ( (gdim[0]) * (gdim[1]) );
}

void TrilinearGrid::getCellVertices(int cell_index,int* cellVertexArray,unsigned int* gdim)
{
	int x,y,z;
	cell2xyz(cell_index,x,y,z,gdim);
	cellVertexArray[0]=TrilinearGrid::xyz2vtx(x  ,y  ,z  ,gdim);
	cellVertexArray[1]=TrilinearGrid::xyz2vtx(x+1,y  ,z  ,gdim);
	cellVertexArray[2]=TrilinearGrid::xyz2vtx(x  ,y+1,z  ,gdim);
	cellVertexArray[3]=TrilinearGrid::xyz2vtx(x+1,y+1,z  ,gdim);
	cellVertexArray[4]=TrilinearGrid::xyz2vtx(x  ,y  ,z+1,gdim);
	cellVertexArray[5]=TrilinearGrid::xyz2vtx(x+1,y  ,z+1,gdim);
	cellVertexArray[6]=TrilinearGrid::xyz2vtx(x  ,y+1,z+1,gdim);
	cellVertexArray[7]=TrilinearGrid::xyz2vtx(x+1,y+1,z+1,gdim);
}

int TrilinearGrid::xyz2vtx(int x, int y, int z, unsigned int* dim)
{
	return x+y*dim[0]+z*dim[0]*dim[1];
}

int TrilinearGrid::getNeighbor(int i, int j, int k, int* neighborArray, unsigned int* dim)
{
	int neighborIndex=0;

	if (i<=0) {
		neighborArray[neighborIndex++]=xyz2vtx(i+1,j,k,dim);
	} else if (i>=(int)(dim[0])-1) {
		neighborArray[neighborIndex++]=xyz2vtx(i-1,j,k,dim);
	} else {
		neighborArray[neighborIndex++]=xyz2vtx(i-1,j,k,dim);
		neighborArray[neighborIndex++]=xyz2vtx(i+1,j,k,dim);
	}

	if (j<=0) {
		neighborArray[neighborIndex++]=xyz2vtx(i,j+1,k,dim);

	} else if (j>=(int)(dim[1])-1) {
		neighborArray[neighborIndex++]=xyz2vtx(i,j-1,k,dim);
	} else {
		neighborArray[neighborIndex++]=xyz2vtx(i,j-1,k,dim);
		neighborArray[neighborIndex++]=xyz2vtx(i,j+1,k,dim);
	}

	if (k<=0) {
		neighborArray[neighborIndex++]=xyz2vtx(i,j,k+1,dim);
	} else if (k>=(int)(dim[2])-1) {
		neighborArray[neighborIndex++]=xyz2vtx(i,j,k-1,dim);
	} else {
		neighborArray[neighborIndex++]=xyz2vtx(i,j,k-1,dim);
		neighborArray[neighborIndex++]=xyz2vtx(i,j,k+1,dim);
	}
	return neighborIndex;
}
